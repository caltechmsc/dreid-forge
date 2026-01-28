use std::path::Path;

use dreid_forge::io::Format;

pub fn bio_input(path: &Path) -> Option<Format> {
    let ext = path.extension()?.to_str()?.to_lowercase();
    match ext.as_str() {
        "pdb" | "ent" => Some(Format::Pdb),
        "cif" | "mmcif" => Some(Format::Mmcif),
        _ => None,
    }
}

pub fn chem_input(path: &Path) -> Option<Format> {
    let ext = path.extension()?.to_str()?.to_lowercase();
    match ext.as_str() {
        "mol2" => Some(Format::Mol2),
        "sdf" | "mol" => Some(Format::Sdf),
        _ => None,
    }
}

pub fn output(path: &Path) -> Option<Format> {
    let ext = path.extension()?.to_str()?.to_lowercase();
    match ext.as_str() {
        // Bio formats
        "bgf" => Some(Format::Bgf),
        "pdb" | "ent" => Some(Format::Pdb),
        "cif" | "mmcif" => Some(Format::Mmcif),
        // Chem formats
        "mol2" => Some(Format::Mol2),
        "sdf" | "mol" => Some(Format::Sdf),
        _ => None,
    }
}
