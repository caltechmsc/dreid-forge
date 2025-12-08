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
