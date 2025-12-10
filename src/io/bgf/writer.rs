use crate::io::error::Error;
use crate::model::metadata::{ResidueCategory, StandardResidue};
use crate::model::topology::ForgedSystem;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;

const DEFAULT_HEADERS: [&str; 2] = ["BIOGRF  332", "FORCEFIELD DREIDING"];
const FORMAT_ATOM: &str =
    "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)";
const FORMAT_CONECT: &str = "FORMAT CONECT (a6,12i6)";

pub fn write<W: Write>(mut writer: W, forged: &ForgedSystem) -> Result<(), Error> {
    let system = &forged.system;
    let metadata = system
        .bio_metadata
        .as_ref()
        .ok_or(Error::MissingMetadata("BGF"))?;

    if metadata.atom_info.len() != system.atoms.len() {
        return Err(Error::Conversion(
            "bio metadata atom_info length must match atom count".into(),
        ));
    }
    if forged.atom_properties.len() != system.atoms.len() {
        return Err(Error::Conversion(
            "atom_properties length must match atom count".into(),
        ));
    }

    let mut indices: Vec<usize> = (0..system.atoms.len()).collect();
    indices.sort_by(|&a, &b| {
        let ia = &metadata.atom_info[a];
        let ib = &metadata.atom_info[b];
        ia.chain_id
            .cmp(&ib.chain_id)
            .then_with(|| ia.residue_id.cmp(&ib.residue_id))
            .then_with(|| ia.insertion_code.cmp(&ib.insertion_code))
            .then_with(|| ia.residue_name.cmp(&ib.residue_name))
            .then_with(|| ia.atom_name.cmp(&ib.atom_name))
    });

    let mut resid_map: HashMap<(char, i32, char), i32> = HashMap::new();
    let mut last_resid_per_chain: HashMap<char, i32> = HashMap::new();
    for &idx in &indices {
        let info = &metadata.atom_info[idx];
        let key = (info.chain_id, info.residue_id, info.insertion_code);
        resid_map.entry(key).or_insert_with(|| {
            let last = last_resid_per_chain.get(&info.chain_id).copied();
            let mut assigned = info.residue_id;
            if let Some(prev) = last.filter(|p| assigned <= *p) {
                assigned = prev + 1;
            }
            last_resid_per_chain.insert(info.chain_id, assigned);
            assigned
        });
    }

    let mut idx_to_serial: HashMap<usize, usize> = HashMap::new();
    for (serial, &idx) in indices.iter().enumerate() {
        idx_to_serial.insert(idx, serial + 1);
    }

    let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
    for bond in &system.bonds {
        let s1 = *idx_to_serial
            .get(&bond.i)
            .ok_or_else(|| Error::Conversion("bond references atom beyond range".into()))?;
        let s2 = *idx_to_serial
            .get(&bond.j)
            .ok_or_else(|| Error::Conversion("bond references atom beyond range".into()))?;
        adjacency.entry(s1).or_default().insert(s2);
        adjacency.entry(s2).or_default().insert(s1);
    }

    for line in DEFAULT_HEADERS {
        writeln!(writer, "{}", line)?;
    }
    writeln!(writer, "{}", FORMAT_ATOM)?;

    for (serial, &idx) in indices.iter().enumerate() {
        let atom = &system.atoms[idx];
        let info = &metadata.atom_info[idx];
        let props = &forged.atom_properties[idx];
        let ff_type = forged
            .atom_types
            .get(props.type_idx)
            .ok_or_else(|| Error::Conversion("atom type index out of range".into()))?;

        let record = match info.category {
            ResidueCategory::Standard if info.standard_name == Some(StandardResidue::HOH) => {
                "HETATM"
            }
            ResidueCategory::Standard => "ATOM  ",
            ResidueCategory::Hetero | ResidueCategory::Ion => "HETATM",
        };

        let resid_out = *resid_map
            .get(&(info.chain_id, info.residue_id, info.insertion_code))
            .ok_or_else(|| Error::Conversion("residue id map missing for atom".into()))?;

        let atoms_connected: usize = adjacency.get(&(serial + 1)).map(|s| s.len()).unwrap_or(0);
        let lone_pairs: usize = 0;

        writeln!(
            writer,
            "{:<6} {:>5} {:<5} {:<3} {:1} {:>5}{:>10.5}{:>10.5}{:>10.5} {:<5}{:>3}{:>2} {:>8.5}{:>10.5}",
            fit_left(record, 6),
            serial + 1,
            fit_left(&info.atom_name, 5),
            fit_left(&info.residue_name, 3),
            info.chain_id,
            resid_out,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            fit_left(ff_type, 5),
            atoms_connected,
            lone_pairs,
            props.charge,
            atom.element.atomic_mass(),
        )?;
    }

    writeln!(writer, "{}", FORMAT_CONECT)?;

    let mut ordered_adj: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (base, neighbors) in adjacency {
        let mut list: Vec<usize> = neighbors.into_iter().collect();
        list.sort_unstable();
        ordered_adj.insert(base, list);
    }

    const MAX_NEIGHBORS: usize = 12;
    for (base, neighbors) in ordered_adj {
        for chunk in neighbors.chunks(MAX_NEIGHBORS) {
            let mut line = format!("CONECT{:>6}", base);
            for &n in chunk {
                line.push_str(&format!("{:>6}", n));
            }
            writeln!(writer, "{}", line)?;
        }
    }

    writeln!(writer, "END")?;
    Ok(())
}

fn fit_left(text: &str, width: usize) -> String {
    let mut s = text.trim().to_string();
    if s.len() > width {
        s.truncate(width);
    }
    format!("{:<width$}", s, width = width)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        metadata::{AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition},
        system::{Bond, System},
        topology::{AtomParam, ForgedSystem, Potentials},
        types::{BondOrder, Element},
    };

    fn sample_forged_system() -> ForgedSystem {
        let atoms = vec![
            Atom::new(Element::N, [0.0, 0.0, 0.0]),
            Atom::new(Element::C, [1.0, 0.0, 0.0]),
            Atom::new(Element::O, [2.0, 0.0, 0.0]),
            Atom::new(Element::C, [3.0, 0.0, 0.0]),
        ];
        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
            Bond::new(2, 3, BondOrder::Single),
        ];
        let meta = BioMetadata {
            atom_info: vec![
                AtomResidueInfo::builder("N", "GLY", 1, 'A')
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("CA", "GLY", 1, 'A')
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("O1", "LIG", 2, 'B')
                    .category(ResidueCategory::Hetero)
                    .position(ResiduePosition::None)
                    .build(),
                AtomResidueInfo::builder("C1", "LIG", 2, 'B')
                    .category(ResidueCategory::Hetero)
                    .position(ResiduePosition::None)
                    .build(),
            ],
        };

        let atom_properties = vec![
            AtomParam {
                charge: -0.3,
                mass: 14.0,
                type_idx: 0,
            },
            AtomParam {
                charge: 0.1,
                mass: 12.0,
                type_idx: 1,
            },
            AtomParam {
                charge: -0.2,
                mass: 16.0,
                type_idx: 2,
            },
            AtomParam {
                charge: 0.0,
                mass: 12.0,
                type_idx: 1,
            },
        ];

        ForgedSystem {
            system: System {
                atoms,
                bonds,
                box_vectors: None,
                bio_metadata: Some(meta),
            },
            atom_types: vec!["N_R".into(), "C_3".into(), "O_2".into()],
            atom_properties,
            potentials: Potentials::default(),
        }
    }

    #[test]
    fn writes_bgf_with_sorted_atoms_and_conect() {
        let forged = sample_forged_system();
        let mut buf = Vec::new();
        write(&mut buf, &forged).expect("write bgf");
        let out = String::from_utf8(buf).expect("utf8");
        let lines: Vec<_> = out.lines().collect();

        assert_eq!(lines[0], "BIOGRF  332");
        assert!(lines.contains(&FORMAT_ATOM));
        assert!(lines.contains(&FORMAT_CONECT));
        assert!(lines.last().unwrap().trim() == "END");

        let atom_lines: Vec<_> = lines
            .iter()
            .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .collect();
        assert_eq!(atom_lines.len(), forged.system.atom_count());

        let fields: Vec<Vec<_>> = atom_lines
            .iter()
            .map(|l| l.split_whitespace().collect())
            .collect();
        assert_eq!(fields[0][0], "ATOM");
        assert_eq!(fields[0][2], "CA");
        assert_eq!(fields[0][3], "GLY");
        assert_eq!(fields[0][4], "A");

        assert_eq!(fields[1][2], "N");
        assert_eq!(fields[1][3], "GLY");
        assert_eq!(fields[1][4], "A");

        assert_eq!(fields[2][0], "HETATM");
        assert_eq!(fields[2][2], "C1");
        assert_eq!(fields[2][3], "LIG");
        assert_eq!(fields[2][4], "B");

        for (i, line) in atom_lines.iter().enumerate() {
            let serial: usize = line[6..12].trim().parse().unwrap();
            assert_eq!(serial, i + 1);
        }

        let conect_lines: Vec<_> = lines.iter().filter(|l| l.starts_with("CONECT")).collect();
        assert!(!conect_lines.is_empty());
        let conect_body = conect_lines
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        assert!(conect_body.contains(" 1") && conect_body.contains(" 2"));
    }

    #[test]
    fn renumbers_insertions_per_chain() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0]),
            Atom::new(Element::C, [1.0, 0.0, 0.0]),
            Atom::new(Element::C, [2.0, 0.0, 0.0]),
        ];
        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Single),
        ];
        let meta = BioMetadata {
            atom_info: vec![
                AtomResidueInfo::builder("C1", "GLY", 10, 'A')
                    .insertion_code_opt(Some(' '))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("C2", "GLY", 10, 'A')
                    .insertion_code_opt(Some('A'))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("C3", "GLY", 11, 'A')
                    .insertion_code_opt(Some(' '))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
            ],
        };

        let atom_properties = vec![
            AtomParam {
                charge: 0.0,
                mass: 12.0,
                type_idx: 0,
            },
            AtomParam {
                charge: 0.0,
                mass: 12.0,
                type_idx: 0,
            },
            AtomParam {
                charge: 0.0,
                mass: 12.0,
                type_idx: 0,
            },
        ];

        let forged = ForgedSystem {
            system: System {
                atoms,
                bonds,
                box_vectors: None,
                bio_metadata: Some(meta),
            },
            atom_types: vec!["C_3".into()],
            atom_properties,
            potentials: Potentials::default(),
        };

        let mut buf = Vec::new();
        write(&mut buf, &forged).expect("write bgf");
        let out = String::from_utf8(buf).expect("utf8");
        let atoms_out: Vec<_> = out
            .lines()
            .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .collect();
        let res_ids: Vec<i32> = atoms_out
            .iter()
            .map(|l| l.split_whitespace().nth(5).unwrap().parse::<i32>().unwrap())
            .collect();
        assert_eq!(res_ids, vec![10, 11, 12]);
    }

    #[test]
    fn treats_water_as_hetatm() {
        let atoms = vec![Atom::new(Element::O, [0.0, 0.0, 0.0])];
        let bonds = vec![];
        let meta = BioMetadata {
            atom_info: vec![
                AtomResidueInfo::builder("O", "HOH", 1, 'A')
                    .standard_name(Some(StandardResidue::HOH))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::None)
                    .build(),
            ],
        };
        let atom_properties = vec![AtomParam {
            charge: 0.0,
            mass: 16.0,
            type_idx: 0,
        }];
        let forged = ForgedSystem {
            system: System {
                atoms,
                bonds,
                box_vectors: None,
                bio_metadata: Some(meta),
            },
            atom_types: vec!["O_2".into()],
            atom_properties,
            potentials: Potentials::default(),
        };

        let mut buf = Vec::new();
        write(&mut buf, &forged).expect("write bgf");
        let out = String::from_utf8(buf).expect("utf8");
        let first = out
            .lines()
            .find(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .unwrap();
        assert!(
            first.starts_with("HETATM"),
            "water should be HETATM: {first}"
        );
    }

    #[test]
    fn errors_when_metadata_missing() {
        let mut forged = sample_forged_system();
        forged.system.bio_metadata = None;
        let err = write(Vec::new(), &forged).expect_err("should fail");
        assert!(matches!(err, Error::MissingMetadata("BGF")));
    }

    #[test]
    fn errors_on_type_index_out_of_range() {
        let mut forged = sample_forged_system();
        forged.atom_types = vec!["N_R".into()];
        let err = write(Vec::new(), &forged).expect_err("type idx fail");
        assert!(matches!(err, Error::Conversion(_)));
    }
}
