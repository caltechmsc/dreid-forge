//! Molecular structure file I/O operations.
//!
//! This module provides readers and writers for various molecular structure
//! file formats, supporting both small-molecule chemistry and biomolecular
//! systems with full biological metadata preservation.
//!
//! # Supported Formats
//!
//! | Format | Read | Write | Use Case |
//! |--------|------|-------|----------|
//! | PDB | ✓ | ✓ | Biomolecules |
//! | mmCIF | ✓ | ✓ | Biomolecules (modern PDB) |
//! | MOL2 | ✓ | ✓ | Small molecules, atom types |
//! | SDF | ✓ | ✓ | Small molecules (MDL) |
//! | BGF | — | ✓ | Force field output (DREIDING) |
//! | LAMMPS | — | ✓ | MD simulation input |
//!
//! # Reader Types
//!
//! - [`ChemReader`] — For small molecules (MOL2, SDF). Direct read without
//!   biological processing.
//! - [`BioReader`] — For biomolecules (PDB, mmCIF). Supports structure
//!   cleaning, protonation, solvation, and topology generation.
//!
//! # Writer Types
//!
//! - [`ChemWriter`] — For small molecules (MOL2, SDF).
//! - [`BioWriter`] — For biomolecules (PDB, mmCIF).
//! - [`write_bgf`] — BGF format for force field visualization.
//! - [`write_lammps_package`] — Complete LAMMPS simulation input files.
//!
//! # Configuration Structs
//!
//! - [`CleanConfig`] — Control structure cleaning (remove water, ions, etc.)
//! - [`ProtonationConfig`] — Control hydrogen addition and histidine tautomers
//! - [`SolvateConfig`] — Control water box and ion placement
//! - [`TopologyConfig`] — Control bond perception and templates

use crate::model::system::System;
use std::collections::HashSet;
use std::fmt;
use std::io::{BufRead, Write};

mod error;
mod util;

mod bgf;
mod lammps;
mod mmcif;
mod mol2;
mod pdb;
mod sdf;

pub use error::Error;

pub use bgf::writer::write as write_bgf;

pub use lammps::writer::{
    LammpsConfig, SystemType, write_data_file as write_lammps_data,
    write_package as write_lammps_package, write_settings_file as write_lammps_settings,
};

pub use bio_forge::Template;

/// Reads a MOL2 template file for heterogeneous residue topology.
///
/// Template files define bond connectivity and atom types for non-standard
/// residues (ligands, cofactors, etc.) that are not in the standard residue
/// library. These templates are used during topology building to correctly
/// perceive bonds in hetero groups.
///
/// # Arguments
///
/// * `reader` — Buffered reader containing MOL2 template data
///
/// # Returns
///
/// A [`Template`] that can be passed to [`TopologyConfig`] for use during
/// biomolecular structure reading.
///
/// # Errors
///
/// Returns [`Error`] if the MOL2 data cannot be parsed or contains
/// invalid structure definitions.
pub fn read_mol2_template<R: BufRead>(reader: R) -> Result<Template, Error> {
    bio_forge::io::read_mol2_template(reader).map_err(Error::from)
}

/// Configuration for structure cleaning operations.
///
/// Controls which atoms, residues, or molecule types are removed
/// during biomolecular structure preparation with [`BioReader`].
///
/// # Examples
///
/// ```
/// use dreid_forge::io::CleanConfig;
/// use std::collections::HashSet;
///
/// // Remove water and common ions
/// let config = CleanConfig {
///     remove_water: true,
///     remove_ions: true,
///     ..Default::default()
/// };
///
/// // Keep only specific residues
/// let mut keep = HashSet::new();
/// keep.insert("ALA".to_string());
/// keep.insert("GLY".to_string());
/// let selective = CleanConfig {
///     keep_residue_names: keep,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone, Default)]
pub struct CleanConfig {
    /// Remove water molecules (residue name HOH, WAT, etc.).
    pub remove_water: bool,
    /// Remove monoatomic ions (Na+, Cl-, etc.).
    pub remove_ions: bool,
    /// Remove all hydrogen atoms.
    pub remove_hydrogens: bool,
    /// Remove all hetero groups (non-standard residues).
    pub remove_hetero: bool,
    /// Specific residue names to remove.
    pub remove_residue_names: HashSet<String>,
    /// If non-empty, keep only these residue names (whitelist mode).
    pub keep_residue_names: HashSet<String>,
}

/// Strategy for histidine tautomer assignment.
///
/// Histidine can exist in three protonation states: HID (Nδ protonated),
/// HIE (Nε protonated), or HIP (doubly protonated, +1 charge). This enum
/// controls how the correct tautomer is selected during protonation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HisStrategy {
    /// Always use HID (Nδ protonated) form.
    DirectHID,
    /// Always use HIE (Nε protonated) form.
    DirectHIE,
    /// Randomly assign HID or HIE with equal probability.
    Random,
    /// Analyze hydrogen bonding network to choose optimal tautomer.
    ///
    /// This is the default and recommended strategy for most applications.
    #[default]
    HbNetwork,
}

/// Configuration for hydrogen addition (protonation).
///
/// Controls how hydrogens are added to a biomolecular structure,
/// including pH-dependent protonation states and histidine handling.
///
/// # Examples
///
/// ```
/// use dreid_forge::io::{ProtonationConfig, HisStrategy};
///
/// // Standard protonation at physiological pH
/// let config = ProtonationConfig {
///     target_ph: Some(7.4),
///     remove_existing_h: true,
///     his_strategy: HisStrategy::HbNetwork,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct ProtonationConfig {
    /// Target pH for protonation state calculation.
    ///
    /// If `None`, uses default protonation states in the input structure.
    pub target_ph: Option<f64>,
    /// Remove existing hydrogens before adding new ones.
    pub remove_existing_h: bool,
    /// Strategy for selecting histidine tautomers.
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

/// Cation species for solvation and charge neutralization.
///
/// Represents common cations that can be added to neutralize
/// system charge or achieve a target ionic strength.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Cation {
    /// Sodium ion (Na⁺).
    Na,
    /// Potassium ion (K⁺).
    K,
    /// Magnesium ion (Mg²⁺).
    Mg,
    /// Calcium ion (Ca²⁺).
    Ca,
    /// Lithium ion (Li⁺).
    Li,
    /// Zinc ion (Zn²⁺).
    Zn,
}

/// Anion species for solvation and charge neutralization.
///
/// Represents common anions that can be added to neutralize
/// system charge or achieve a target ionic strength.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Anion {
    /// Chloride ion (Cl⁻).
    Cl,
    /// Bromide ion (Br⁻).
    Br,
    /// Iodide ion (I⁻).
    I,
    /// Fluoride ion (F⁻).
    F,
}

/// Configuration for solvation and ion placement.
///
/// Controls how a water box is added around the solute and how
/// counterions are placed to neutralize the system charge.
///
/// # Examples
///
/// ```
/// use dreid_forge::io::{SolvateConfig, Cation, Anion};
///
/// // Create a 12 Å water box with NaCl ions
/// let config = SolvateConfig {
///     margin: 12.0,
///     cations: vec![Cation::Na],
///     anions: vec![Anion::Cl],
///     target_charge: 0, // neutralize
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct SolvateConfig {
    /// Minimum distance from solute to box edge (Å).
    pub margin: f64,
    /// Grid spacing for water molecule placement (Å).
    pub water_spacing: f64,
    /// Van der Waals cutoff for water clash detection (Å).
    pub vdw_cutoff: f64,
    /// Remove existing water/ions before solvating.
    pub remove_existing: bool,
    /// Cation species to add for neutralization.
    pub cations: Vec<Cation>,
    /// Anion species to add for neutralization.
    pub anions: Vec<Anion>,
    /// Target net charge after ion addition.
    pub target_charge: i32,
    /// Random seed for reproducible ion placement.
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

/// Configuration for molecular topology generation.
///
/// Controls bond perception and residue template matching during
/// biomolecular structure reading.
///
/// # Examples
///
/// ```no_run
/// use dreid_forge::io::{TopologyConfig, read_mol2_template};
/// use std::fs::File;
/// use std::io::BufReader;
///
/// // Load a ligand template for topology building
/// let file = File::open("ligand.mol2")?;
/// let template = read_mol2_template(BufReader::new(file))?;
///
/// let config = TopologyConfig {
///     hetero_templates: vec![template],
///     disulfide_bond_cutoff: 2.2,
/// };
/// # Ok::<(), dreid_forge::io::Error>(())
/// ```
#[derive(Debug, Clone)]
pub struct TopologyConfig {
    /// Templates for non-standard residues (ligands, cofactors).
    ///
    /// Each template defines the expected atom names, bond connectivity,
    /// and atom types for a heterogeneous residue.
    pub hetero_templates: Vec<Template>,
    /// Maximum S-S distance (Å) for disulfide bond detection.
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

/// Molecular structure file format enumeration.
///
/// Identifies the format of input or output files for reader/writer
/// construction and error reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Protein Data Bank format (`.pdb`).
    Pdb,
    /// Macromolecular Crystallographic Information File (`.cif`).
    Mmcif,
    /// MDL Structure Data File (`.sdf`, `.mol`).
    Sdf,
    /// Tripos MOL2 format (`.mol2`).
    Mol2,
    /// Biograf format for force field output (`.bgf`).
    Bgf,
    /// LAMMPS data file format (`.data`, `.in`).
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

/// Reader for small-molecule structure files.
///
/// Provides a simple interface for reading molecular structures from
/// chemistry-focused formats like MOL2 and SDF. These formats do not
/// require biological processing and are read directly.
///
/// For biomolecular formats (PDB, mmCIF) that require structure
/// preparation, use [`BioReader`] instead.
///
/// # Supported Formats
///
/// * [`Format::Mol2`] — Tripos MOL2 with atom types and bonds
/// * [`Format::Sdf`] — MDL SDF/MOL with connection table
///
/// # Examples
///
/// ```no_run
/// use dreid_forge::io::{ChemReader, Format};
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("molecule.mol2")?;
/// let system = ChemReader::new(BufReader::new(file), Format::Mol2).read()?;
///
/// println!("Loaded {} atoms", system.atom_count());
/// # Ok::<(), dreid_forge::io::Error>(())
/// ```
pub struct ChemReader<R: BufRead> {
    reader: R,
    format: Format,
}

impl<R: BufRead> ChemReader<R> {
    /// Creates a new chemistry file reader.
    ///
    /// # Arguments
    ///
    /// * `reader` — Buffered reader containing the file data
    /// * `format` — File format (must be MOL2 or SDF)
    pub fn new(reader: R, format: Format) -> Self {
        Self { reader, format }
    }

    /// Reads and parses the molecular structure.
    ///
    /// # Returns
    ///
    /// A [`System`] containing atoms and bonds from the file.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] if:
    /// * The format is not supported for reading (use [`BioReader`] for PDB/mmCIF)
    /// * The file contains invalid or malformed data
    pub fn read(self) -> Result<System, Error> {
        match self.format {
            Format::Mol2 => mol2::reader::read(self.reader),
            Format::Sdf => sdf::reader::read(self.reader),
            _ => Err(Error::UnsupportedReadFormat(self.format)),
        }
    }
}

/// Writer for small-molecule structure files.
///
/// Provides a simple interface for writing molecular structures to
/// chemistry-focused formats like MOL2 and SDF.
///
/// For biomolecular formats (PDB, mmCIF) that require biological
/// metadata, use [`BioWriter`] instead.
///
/// # Supported Formats
///
/// * [`Format::Mol2`] — Tripos MOL2 with atom types and bonds
/// * [`Format::Sdf`] — MDL SDF/MOL with connection table
///
/// # Examples
///
/// ```no_run
/// use dreid_forge::io::{ChemWriter, Format};
/// use dreid_forge::System;
/// use std::fs::File;
///
/// let system = System::new(); // your molecular system
/// let file = File::create("output.mol2")?;
/// ChemWriter::new(file, Format::Mol2).write(&system)?;
/// # Ok::<(), dreid_forge::io::Error>(())
/// ```
pub struct ChemWriter<W: Write> {
    writer: W,
    format: Format,
}

impl<W: Write> ChemWriter<W> {
    /// Creates a new chemistry file writer.
    ///
    /// # Arguments
    ///
    /// * `writer` — Output destination implementing [`Write`]
    /// * `format` — Output file format (must be MOL2 or SDF)
    pub fn new(writer: W, format: Format) -> Self {
        Self { writer, format }
    }

    /// Writes the molecular structure to the output.
    ///
    /// # Arguments
    ///
    /// * `system` — The molecular system to write
    ///
    /// # Errors
    ///
    /// Returns [`Error`] if:
    /// * The format is not supported for writing
    /// * An I/O error occurs during writing
    pub fn write(self, system: &System) -> Result<(), Error> {
        match self.format {
            Format::Mol2 => mol2::writer::write(self.writer, system),
            Format::Sdf => sdf::writer::write(self.writer, system),
            _ => Err(Error::UnsupportedWriteFormat(self.format)),
        }
    }
}

/// Reader for biomolecular structure files with preparation pipeline.
///
/// Provides a builder-style interface for reading biomolecular structures
/// from PDB or mmCIF format with optional structure preparation steps:
/// cleaning, protonation, solvation, and topology generation.
///
/// # Preparation Pipeline
///
/// 1. **Clean** — Remove unwanted atoms (water, ions, hetero groups)
/// 2. **Repair** — Fix missing atoms/residues
/// 3. **Protonate** — Add hydrogens with appropriate protonation states
/// 4. **Solvate** — Add water box and counterions (optional)
/// 5. **Topology** — Build bond connectivity from residue templates
///
/// # Supported Formats
///
/// * [`Format::Pdb`] — Protein Data Bank format
/// * [`Format::Mmcif`] — Macromolecular CIF format
///
/// # Examples
///
/// ```no_run
/// use dreid_forge::io::{BioReader, Format, CleanConfig, ProtonationConfig};
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("protein.pdb")?;
/// let system = BioReader::new(BufReader::new(file), Format::Pdb)
///     .clean(CleanConfig {
///         remove_water: true,
///         remove_ions: true,
///         ..Default::default()
///     })
///     .protonate(ProtonationConfig::default())
///     .read()?;
///
/// println!("Prepared system with {} atoms", system.atom_count());
/// # Ok::<(), dreid_forge::io::Error>(())
/// ```
pub struct BioReader<R: BufRead> {
    reader: R,
    format: Format,
    clean_config: Option<CleanConfig>,
    protonation_config: ProtonationConfig,
    solvate_config: Option<SolvateConfig>,
    topology_config: TopologyConfig,
}

impl<R: BufRead> BioReader<R> {
    /// Creates a new biomolecular file reader.
    ///
    /// # Arguments
    ///
    /// * `reader` — Buffered reader containing PDB or mmCIF data
    /// * `format` — File format (must be PDB or mmCIF)
    pub fn new(reader: R, format: Format) -> Self {
        Self {
            reader,
            format,
            clean_config: None,
            protonation_config: ProtonationConfig::default(),
            solvate_config: None,
            topology_config: TopologyConfig::default(),
        }
    }

    /// Enables structure cleaning with the given configuration.
    ///
    /// # Arguments
    ///
    /// * `config` — Cleaning options (remove water, ions, etc.)
    pub fn clean(mut self, config: CleanConfig) -> Self {
        self.clean_config = Some(config);
        self
    }

    /// Sets the protonation configuration.
    ///
    /// # Arguments
    ///
    /// * `config` — Protonation options (pH, histidine handling)
    pub fn protonate(mut self, config: ProtonationConfig) -> Self {
        self.protonation_config = config;
        self
    }

    /// Enables solvation with the given configuration.
    ///
    /// # Arguments
    ///
    /// * `config` — Solvation options (water box, ions)
    pub fn solvate(mut self, config: SolvateConfig) -> Self {
        self.solvate_config = Some(config);
        self
    }

    /// Sets the topology building configuration.
    ///
    /// # Arguments
    ///
    /// * `config` — Topology options (templates, disulfide detection)
    pub fn topology(mut self, config: TopologyConfig) -> Self {
        self.topology_config = config;
        self
    }

    /// Reads and prepares the biomolecular structure.
    ///
    /// Applies all configured preparation steps in order:
    /// clean → protonate → solvate → topology.
    ///
    /// # Returns
    ///
    /// A fully prepared [`System`] with biological metadata and bonds.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] if:
    /// * The format is not supported for reading
    /// * Parsing fails due to invalid file content
    /// * Any preparation step encounters an error
    pub fn read(self) -> Result<System, Error> {
        match self.format {
            Format::Pdb => pdb::reader::read(self),
            Format::Mmcif => mmcif::reader::read(self),
            _ => Err(Error::UnsupportedReadFormat(self.format)),
        }
    }
}

/// Writer for biomolecular structure files.
///
/// Writes molecular structures with biological metadata to PDB or
/// mmCIF format. The system must have [`BioMetadata`](crate::BioMetadata)
/// attached for proper residue and chain information output.
///
/// # Supported Formats
///
/// * [`Format::Pdb`] — Protein Data Bank format
/// * [`Format::Mmcif`] — Macromolecular CIF format
///
/// # Examples
///
/// ```no_run
/// use dreid_forge::io::{BioWriter, Format};
/// use dreid_forge::System;
/// use std::fs::File;
///
/// let system = System::new(); // your prepared biomolecular system
/// let file = File::create("output.pdb")?;
/// BioWriter::new(file, Format::Pdb).write(&system)?;
/// # Ok::<(), dreid_forge::io::Error>(())
/// ```
pub struct BioWriter<W: Write> {
    writer: W,
    format: Format,
}

impl<W: Write> BioWriter<W> {
    /// Creates a new biomolecular file writer.
    ///
    /// # Arguments
    ///
    /// * `writer` — Output destination implementing [`Write`]
    /// * `format` — Output file format (must be PDB or mmCIF)
    pub fn new(writer: W, format: Format) -> Self {
        Self { writer, format }
    }

    /// Writes the biomolecular structure to the output.
    ///
    /// # Arguments
    ///
    /// * `system` — The molecular system to write (must have bio_metadata)
    ///
    /// # Errors
    ///
    /// Returns [`Error`] if:
    /// * The format is not supported for writing
    /// * The system lacks required biological metadata
    /// * An I/O error occurs during writing
    pub fn write(self, system: &System) -> Result<(), Error> {
        match self.format {
            Format::Pdb => pdb::writer::write(self.writer, system),
            Format::Mmcif => mmcif::writer::write(self.writer, system),
            _ => Err(Error::UnsupportedWriteFormat(self.format)),
        }
    }
}
