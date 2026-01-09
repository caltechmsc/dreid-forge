use crate::io::{error::Error, util};
use crate::model::system::System;
use bio_forge as bf;
use std::io::Write;

pub fn write<W: Write>(writer: W, system: &System) -> Result<(), Error> {
    if system.bio_metadata.is_none() {
        return Err(Error::MissingMetadata("PDB"));
    }

    let topo = util::to_bio_topology(system)?;
    bf::io::write_pdb_topology(writer, &topo).map_err(Error::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{BioReader, BioWriter, Format, Template, TopologyConfig};
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
                info.chain_id,
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
                        info_i.chain_id,
                        info_i.residue_id,
                        info_i.insertion_code,
                        info_i.atom_name.clone(),
                    ))
                    .unwrap();
                let j = *map
                    .get(&(
                        info_j.chain_id,
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
    fn writes_and_reads_back_consistently() {
        let (system, template) = sample_system();
        let mut buf = Vec::new();
        BioWriter::new(&mut buf, Format::Pdb)
            .write(&system)
            .expect("write pdb");

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
    fn errors_without_metadata() {
        let system = System {
            atoms: vec![Atom::new(Element::C, [0.0, 0.0, 0.0])],
            bonds: vec![],
            box_vectors: None,
            bio_metadata: None,
        };

        let err = BioWriter::new(Vec::new(), Format::Pdb)
            .write(&system)
            .expect_err("missing metadata should fail");
        match err {
            Error::MissingMetadata(fmt) => assert_eq!(fmt, "PDB"),
            other => panic!("unexpected error: {:?}", other),
        }
    }
}
