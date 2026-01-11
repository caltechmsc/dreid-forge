use crate::io::{BioReader, error::Error, util};
use crate::model::system::System;
use bio_forge as bf;
use std::io::BufRead;

pub fn read<R: BufRead>(builder: BioReader<R>) -> Result<System, Error> {
    let bio_context = bf::io::IoContext::new_default();
    let mut bio_struct = bf::io::read_pdb_structure(builder.reader, &bio_context)?;

    if let Some(clean_config) = builder.clean_config {
        let bf_clean_config = util::to_bf_clean_config(clean_config);
        bf::ops::clean_structure(&mut bio_struct, &bf_clean_config)?;
    }

    bf::ops::repair_structure(&mut bio_struct)?;

    let bf_hydro_config = util::to_bf_hydro_config(builder.protonation_config.clone());
    bf::ops::add_hydrogens(&mut bio_struct, &bf_hydro_config)?;

    if let Some(solvate_config) = builder.solvate_config {
        let bf_solvate_config = util::to_bf_solvate_config(solvate_config);
        bf::ops::solvate_structure(&mut bio_struct, &bf_solvate_config)?;
    }

    let bf_topo = util::build_topology(bio_struct, &builder.topology_config)?;
    let mut system = util::from_bio_topology(bf_topo)?;

    if let Some(metadata) = system.bio_metadata.as_mut() {
        metadata.target_ph = builder.protonation_config.target_ph;
    }

    Ok(system)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{Format, Template, TopologyConfig};
    use crate::model::{
        atom::Atom,
        metadata::{AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition},
        system::{Bond, System},
        types::{BondOrder, Element},
    };
    use std::collections::HashSet;
    use std::io::Cursor;

    fn sample_system() -> (System, Template) {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0]),
            Atom::new(Element::O, [1.0, 0.0, 0.0]),
        ];
        let metadata = BioMetadata {
            atom_info: vec![
                AtomResidueInfo::builder("C1", "LIG", 1, 'A')
                    .category(ResidueCategory::Hetero)
                    .position(ResiduePosition::None)
                    .build(),
                AtomResidueInfo::builder("O1", "LIG", 1, 'A')
                    .category(ResidueCategory::Hetero)
                    .position(ResiduePosition::None)
                    .build(),
            ],
            target_ph: None,
        };
        let bonds = vec![Bond::new(0, 1, BondOrder::Double)];
        let system = System {
            atoms,
            bonds,
            box_vectors: Some([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            bio_metadata: Some(metadata),
        };

        let template = Template::new(
            "LIG",
            vec!["C1".into(), "O1".into()],
            vec![("C1".into(), "O1".into(), bf::BondOrder::Double)],
        );

        (system, template)
    }

    fn assert_systems_equivalent(expected: &System, actual: &System) {
        match (expected.box_vectors, actual.box_vectors) {
            (Some(a), Some(b)) => {
                for i in 0..3 {
                    for j in 0..3 {
                        assert!(
                            (a[i][j] - b[i][j]).abs() < 1e-9,
                            "box vectors differ at [{i}][{j}]: {} vs {}",
                            a[i][j],
                            b[i][j]
                        );
                    }
                }
            }
            (None, None) => {}
            _ => panic!("box vector presence mismatch"),
        }

        let meta_a = expected.bio_metadata.as_ref().expect("expected metadata");
        let meta_b = actual.bio_metadata.as_ref().expect("actual metadata");
        assert_eq!(
            meta_a.atom_info.len(),
            meta_b.atom_info.len(),
            "atom_info len"
        );
        assert_eq!(expected.atom_count(), actual.atom_count(), "atom count");
        assert_eq!(expected.bond_count(), actual.bond_count(), "bond count");

        let key = |info: &AtomResidueInfo| {
            (
                info.chain_id.clone(),
                info.residue_id,
                info.insertion_code,
                info.atom_name.clone(),
            )
        };
        let mut map = std::collections::HashMap::new();
        for (i, info) in meta_a.atom_info.iter().enumerate() {
            map.insert(key(info), i);
        }

        for (j, info_b) in meta_b.atom_info.iter().enumerate() {
            let i = *map.get(&key(info_b)).expect("atom key missing in expected");
            assert_eq!(
                expected.atoms[i].element, actual.atoms[j].element,
                "element mismatch"
            );
            for k in 0..3 {
                assert!((expected.atoms[i].position[k] - actual.atoms[j].position[k]).abs() < 1e-6);
            }
            assert_eq!(meta_a.atom_info[i], *info_b, "metadata mismatch");
        }

        let bonds_a: HashSet<_> = expected.bonds.iter().cloned().collect();
        let bonds_b: HashSet<_> = actual
            .bonds
            .iter()
            .map(|b| {
                let info_i = &meta_b.atom_info[b.i];
                let info_j = &meta_b.atom_info[b.j];
                let i = *map
                    .get(&(
                        info_i.chain_id.clone(),
                        info_i.residue_id,
                        info_i.insertion_code,
                        info_i.atom_name.clone(),
                    ))
                    .unwrap();
                let j = *map
                    .get(&(
                        info_j.chain_id.clone(),
                        info_j.residue_id,
                        info_j.insertion_code,
                        info_j.atom_name.clone(),
                    ))
                    .unwrap();
                Bond::new(i, j, b.order)
            })
            .collect();
        assert_eq!(bonds_a, bonds_b, "bond sets differ");
    }

    #[test]
    fn reads_pdb_topology_roundtrip() {
        let (system, template) = sample_system();
        let bf_topo = util::to_bio_topology(&system).expect("to bio topo");

        let mut buf = Vec::new();
        bf::io::write_pdb_topology(&mut buf, &bf_topo).expect("write bf pdb");

        let roundtrip = BioReader::new(Cursor::new(buf), Format::Pdb)
            .topology(TopologyConfig {
                hetero_templates: vec![template],
                disulfide_bond_cutoff: 2.2,
            })
            .read()
            .expect("read pdb");

        assert_systems_equivalent(&system, &roundtrip);
    }

    #[test]
    fn handles_invalid_input_gracefully() {
        let bogus = b"ATOM".as_slice();
        let system = BioReader::new(Cursor::new(bogus), Format::Pdb)
            .read()
            .expect("parsing should not panic");
        assert_eq!(
            system.atom_count(),
            0,
            "expected empty system for invalid input"
        );
    }
}
