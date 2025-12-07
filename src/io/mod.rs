use std::fmt;

pub mod error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Pdb,
    Mmcif,
    Sdf,
    Mol2,
    Bgf,
    Lammps,
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Format::Pdb => write!(f, "PDB"),
            Format::Mmcif => write!(f, "mmCIF"),
            Format::Sdf => write!(f, "SDF"),
            Format::Mol2 => write!(f, "MOL2"),
            Format::Bgf => write!(f, "BGF"),
            Format::Lammps => write!(f, "LAMMPS"),
        }
    }
}
