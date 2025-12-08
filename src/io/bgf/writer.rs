use crate::io::error::Error;
use crate::model::metadata::{ResidueCategory, StandardResidue};
use crate::model::topology::ForgedSystem;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;

const DEFAULT_HEADERS: [&str; 2] = ["BIOGRF  332", "FORCEFIELD DREIDING"];
const FORMAT_ATOM: &str =
    "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)";
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
            if let Some(prev) = last {
                if assigned <= prev {
                    assigned = prev + 1;
                }
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
            .expect("residue id map must exist");

        let atoms_connected: usize = adjacency.get(&(serial + 1)).map(|s| s.len()).unwrap_or(0);
        let lone_pairs: usize = 0;

        writeln!(
            writer,
            "{:<6}{:>5} {:<5} {:>3} {:1} {:>5}{:>10.5}{:>10.5}{:>10.5} {:<5}{:>3}{:>2}{:>8.5}",
            record,
            serial + 1,
            fit_left(&info.atom_name, 5),
            fit_right(&info.residue_name, 3),
            info.chain_id,
            resid_out,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            fit_left(ff_type, 5),
            atoms_connected,
            lone_pairs,
            props.charge,
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

fn fit_right(text: &str, width: usize) -> String {
    let mut s = text.trim().to_string();
    if s.len() > width {
        s.truncate(width);
    }
    format!("{:>width$}", s, width = width)
}
