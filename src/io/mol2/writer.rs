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
