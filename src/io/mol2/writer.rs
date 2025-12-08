use crate::io::{error::Error, util};
use crate::model::system::System;
use std::io::Write;

pub fn write<W: Write>(mut writer: W, system: &System) -> Result<(), Error> {
    let atom_count = system.atom_count();
    let bond_count = system.bond_count();

    writeln!(writer, "@<TRIPOS>MOLECULE")?;
    writeln!(writer, "DREID-FORGE")?;
    writeln!(writer, "{:>5} {:>5} 0 0 0", atom_count, bond_count)?;
    writeln!(writer, "SMALL")?;
    writeln!(writer, "NO_CHARGES")?;
    writeln!(writer, "****")?;
    writeln!(writer)?;

    writeln!(writer, "@<TRIPOS>ATOM")?;
    for (i, atom) in system.atoms.iter().enumerate() {
        let name = system
            .bio_metadata
            .as_ref()
            .and_then(|m| m.atom_info.get(i))
            .map(|info| info.atom_name.as_str())
            .unwrap_or_else(|| atom.element.symbol());
        let atom_type = atom.element.symbol();
        writeln!(
            writer,
            "{:>7} {:<8} {:>10.4} {:>10.4} {:>10.4} {:<6} {:>3} {:<8} {:>8.4}",
            i + 1,
            name,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            atom_type,
            1,
            "SYS",
            0.0
        )?;
    }

    writeln!(writer, "@<TRIPOS>BOND")?;
    for (i, bond) in system.bonds.iter().enumerate() {
        writeln!(
            writer,
            "{:>7} {:>4} {:>4} {}",
            i + 1,
            bond.i + 1,
            bond.j + 1,
            util::bond_order_to_mol2(bond.order)
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::mol2::reader;
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
                Atom::new(Element::O, [1.2, 0.0, 0.0]),
                Atom::new(Element::N, [0.0, 1.1, 0.0]),
            ],
            bonds: vec![
                Bond::new(0, 1, BondOrder::Double),
                Bond::new(0, 2, BondOrder::Single),
            ],
            box_vectors: None,
            bio_metadata: None,
        };

        let mut buf = Vec::new();
        write(&mut buf, &system).expect("write mol2");

        let parsed = reader::read(Cursor::new(buf)).expect("read mol2");
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
