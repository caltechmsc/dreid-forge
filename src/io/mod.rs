use std::collections::HashSet;
use std::fmt;

pub mod error;
pub mod util;

pub mod pdb;

pub use bio_forge::Template;

#[derive(Debug, Clone, Default)]
pub struct CleanConfig {
    pub remove_water: bool,
    pub remove_ions: bool,
    pub remove_hydrogens: bool,
    pub remove_hetero: bool,
    pub remove_residue_names: HashSet<String>,
    pub keep_residue_names: HashSet<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HisStrategy {
    DirectHID,
    DirectHIE,
    Random,
    #[default]
    HbNetwork,
}

#[derive(Debug, Clone)]
pub struct ProtonationConfig {
    pub target_ph: Option<f64>,
    pub remove_existing_h: bool,
    pub his_strategy: HisStrategy,
}

impl Default for ProtonationConfig {
    fn default() -> Self {
        Self {
            target_ph: None,
            remove_existing_h: true,
            his_strategy: HisStrategy::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Cation {
    Na,
    K,
    Mg,
    Ca,
    Li,
    Zn,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Anion {
    Cl,
    Br,
    I,
    F,
}

#[derive(Debug, Clone)]
pub struct SolvateConfig {
    pub margin: f64,
    pub water_spacing: f64,
    pub vdw_cutoff: f64,
    pub remove_existing: bool,
    pub cations: Vec<Cation>,
    pub anions: Vec<Anion>,
    pub target_charge: i32,
    pub rng_seed: Option<u64>,
}

impl Default for SolvateConfig {
    fn default() -> Self {
        Self {
            margin: 10.0,
            water_spacing: 3.1,
            vdw_cutoff: 2.4,
            remove_existing: true,
            cations: vec![Cation::Na],
            anions: vec![Anion::Cl],
            target_charge: 0,
            rng_seed: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TopologyConfig {
    pub hetero_templates: Vec<Template>,
    pub disulfide_bond_cutoff: f64,
}

impl Default for TopologyConfig {
    fn default() -> Self {
        Self {
            hetero_templates: Vec::new(),
            disulfide_bond_cutoff: 2.2,
        }
    }
}

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
