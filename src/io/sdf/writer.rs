use crate::io::{error::Error, util};
use crate::model::system::System;
use std::io::Write;

pub fn write<W: Write>(mut writer: W, system: &System) -> Result<(), Error> {
    let atom_count = system.atom_count();
    let bond_count = system.bond_count();

    writeln!(writer, "SDF Export")?;
    writeln!(writer, "dreid-forge")?;
    writeln!(writer, "")?;
    writeln!(
        writer,
        "{:>3}{:>3}  0  0  0  0  0  0  0  0  0999 V2000",
        atom_count, bond_count
    )?;

    for atom in &system.atoms {
        writeln!(
            writer,
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
            atom.position[0],
            atom.position[1],
            atom.position[2],
            atom.element.symbol()
        )?;
    }

    for bond in &system.bonds {
        writeln!(
            writer,
            "{:>3}{:>3}{:>3}  0  0  0  0",
            bond.i + 1,
            bond.j + 1,
            util::bond_order_to_ctfile(bond.order)
        )?;
    }

    writeln!(writer, "M  END")?;
    writeln!(writer, "$$$$")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::sdf::reader;
    use crate::model::{
        atom::Atom,
        system::Bond,
        types::{BondOrder, Element},
    };
    use std::io::Cursor;

    #[test]
    fn writes_and_reads_roundtrip() {
        let system = System {
            atoms: vec![
                Atom::new(Element::C, [0.0, 0.0, 0.0]),
                Atom::new(Element::O, [1.0, 0.0, 0.0]),
                Atom::new(Element::H, [0.0, 1.0, 0.0]),
            ],
            bonds: vec![
                Bond::new(0, 1, BondOrder::Double),
                Bond::new(0, 2, BondOrder::Single),
            ],
            box_vectors: None,
            bio_metadata: None,
        };

        let mut buf = Vec::new();
        write(&mut buf, &system).expect("write sdf");
        let parsed = reader::read(Cursor::new(buf)).expect("read sdf");

        assert_eq!(parsed.atom_count(), system.atom_count());
        assert_eq!(parsed.bond_count(), system.bond_count());
        for (a, b) in system.atoms.iter().zip(parsed.atoms.iter()) {
            assert_eq!(a.element, b.element);
            for k in 0..3 {
                assert!((a.position[k] - b.position[k]).abs() < 1e-4);
            }
        }
        assert_eq!(parsed.bonds, system.bonds);
    }
}
