use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(Parser)]
#[command(
    name = "dforge",
    about = "DREIDING force field parameterization",
    version,
    author,
    before_help = crate::display::banner_for_help(),
    propagate_version = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Parameterize biological macromolecules (PDB/mmCIF)
    #[command(visible_alias = "b")]
    Bio(BioArgs),

    /// Parameterize small molecules (MOL2/SDF)
    #[command(visible_alias = "c")]
    Chem(ChemArgs),
}

/// I/O options shared by all commands.
#[derive(Args)]
pub struct IoOptions {
    /// Input file (stdin if omitted, requires --infmt)
    #[arg(short, long, value_name = "FILE")]
    pub input: Option<PathBuf>,

    /// Output file(s), repeatable for multi-format output
    #[arg(short, long, value_name = "FILE", action = clap::ArgAction::Append)]
    pub output: Vec<PathBuf>,

    /// Suppress progress output (for scripting)
    #[arg(short, long)]
    pub quiet: bool,
}

/// Force field options shared by bio and chem commands.
#[derive(Args)]
#[command(next_help_heading = "Force Field Options")]
pub struct ForgeOptions {
    /// Charge calculation method
    #[arg(long, value_name = "METHOD", default_value = "qeq")]
    pub charge: ChargeMethod,

    /// Total system charge (for QEq constraint)
    #[arg(
        long,
        value_name = "Q",
        default_value = "0.0",
        allow_hyphen_values = true
    )]
    pub total_charge: f64,

    /// Bond potential functional form
    #[arg(long, value_name = "TYPE", default_value = "harmonic")]
    pub bond_potential: BondPotential,

    /// Angle potential functional form
    #[arg(long, value_name = "TYPE", default_value = "theta-harmonic")]
    pub angle_potential: AnglePotential,

    /// Van der Waals potential functional form
    #[arg(long, value_name = "TYPE", default_value = "lj")]
    pub vdw_potential: VdwPotential,

    /// Custom typing rules (TOML file)
    #[arg(long, value_name = "FILE")]
    pub rules: Option<PathBuf>,

    /// Custom force field parameters (TOML file)
    #[arg(long, value_name = "FILE")]
    pub params: Option<PathBuf>,
}

/// LAMMPS output options shared by bio and chem commands.
#[derive(Args)]
#[command(next_help_heading = "LAMMPS Output Options")]
pub struct LammpsOptions {
    /// Non-bonded cutoff distance (Å)
    #[arg(long = "lmp-cutoff", value_name = "Å", default_value = "12.0")]
    pub nonbonded_cutoff: f64,

    /// Hydrogen bond inner cutoff (Å)
    #[arg(long = "lmp-hbond-inner", value_name = "Å", default_value = "10.0")]
    pub hbond_inner: f64,

    /// Hydrogen bond outer cutoff (Å)
    #[arg(long = "lmp-hbond-outer", value_name = "Å", default_value = "12.0")]
    pub hbond_outer: f64,

    /// Hydrogen bond angle cutoff (degrees)
    #[arg(long = "lmp-hbond-angle", value_name = "DEG", default_value = "90.0")]
    pub hbond_angle: f64,

    /// Box margin for non-periodic systems (Å)
    #[arg(long = "lmp-margin", value_name = "Å", default_value = "5.0")]
    pub aabb_margin: f64,

    /// Boundary condition type
    #[arg(long = "lmp-system", value_name = "TYPE", default_value = "periodic")]
    pub system_type: SystemBoundary,
}

#[derive(Args)]
pub struct BioArgs {
    #[command(flatten)]
    pub io: IoOptions,

    /// Input format (inferred from extension if not specified)
    #[arg(long = "infmt", value_name = "FORMAT")]
    pub input_format: Option<BioInputFormat>,

    /// Output format for first/only output
    #[arg(long = "outfmt", value_name = "FORMAT")]
    pub output_format: Option<BioOutputFormat>,

    #[command(flatten)]
    pub clean: CleanOptions,

    #[command(flatten)]
    pub protonation: ProtonationOptions,

    #[command(flatten)]
    pub solvation: SolvationOptions,

    #[command(flatten)]
    pub topology: TopologyOptions,

    #[command(flatten)]
    pub forge: ForgeOptions,

    #[command(flatten)]
    pub lammps: LammpsOptions,
}

#[derive(Args)]
#[command(next_help_heading = "Structure Cleaning")]
pub struct CleanOptions {
    /// Remove water molecules
    #[arg(long)]
    pub no_water: bool,

    /// Remove ions
    #[arg(long)]
    pub no_ions: bool,

    /// Remove existing hydrogens
    #[arg(long)]
    pub no_hydrogens: bool,

    /// Remove hetero atoms (HETATM)
    #[arg(long)]
    pub no_hetero: bool,

    /// Remove specific residues (comma-separated names)
    #[arg(long, value_name = "RES", value_delimiter = ',')]
    pub remove: Vec<String>,

    /// Keep only these residues (comma-separated names)
    #[arg(long, value_name = "RES", value_delimiter = ',')]
    pub keep: Vec<String>,
}

#[derive(Args)]
#[command(next_help_heading = "Protonation")]
pub struct ProtonationOptions {
    /// Target pH for protonation state assignment
    #[arg(long, value_name = "PH", default_value = "7.0")]
    pub ph: f64,

    /// Histidine tautomer selection strategy
    #[arg(long, value_name = "STRATEGY", default_value = "network")]
    pub his: HisStrategy,
}

#[derive(Args)]
#[command(next_help_heading = "Solvation")]
pub struct SolvationOptions {
    /// Enable solvation (add water box)
    #[arg(long)]
    pub solvate: bool,

    /// Water box margin around solute (Å)
    #[arg(long = "solv-margin", value_name = "Å", default_value = "10.0")]
    pub box_margin: f64,

    /// Water molecule spacing (Å)
    #[arg(long = "solv-spacing", value_name = "Å", default_value = "3.1")]
    pub spacing: f64,

    /// Minimum solute-water distance (Å)
    #[arg(long = "solv-cutoff", value_name = "Å", default_value = "2.4")]
    pub vdw_cutoff: f64,

    /// Cation type for neutralization
    #[arg(long = "solv-cation", value_name = "ION", default_value = "na")]
    pub cation: Cation,

    /// Anion type for neutralization
    #[arg(long = "solv-anion", value_name = "ION", default_value = "cl")]
    pub anion: Anion,

    /// Target net charge after ion addition
    #[arg(
        long = "solv-charge",
        value_name = "Q",
        default_value = "0",
        allow_hyphen_values = true
    )]
    pub target_charge: i32,

    /// Random seed for reproducible ion placement
    #[arg(long = "solv-seed", value_name = "SEED")]
    pub seed: Option<u64>,
}

#[derive(Args)]
#[command(next_help_heading = "Topology")]
pub struct TopologyOptions {
    /// Disulfide bond detection cutoff (Å)
    #[arg(long = "ss-cutoff", value_name = "Å", default_value = "2.2")]
    pub ss_cutoff: f64,

    /// Hetero residue template files (MOL2)
    #[arg(long = "template", value_name = "FILE", action = clap::ArgAction::Append)]
    pub templates: Vec<PathBuf>,
}

#[derive(Args)]
pub struct ChemArgs {
    #[command(flatten)]
    pub io: IoOptions,

    /// Input format (inferred from extension if not specified)
    #[arg(long = "infmt", value_name = "FORMAT")]
    pub input_format: Option<ChemInputFormat>,

    /// Output format for first/only output
    #[arg(long = "outfmt", value_name = "FORMAT")]
    pub output_format: Option<ChemOutputFormat>,

    #[command(flatten)]
    pub forge: ForgeOptions,

    #[command(flatten)]
    pub lammps: LammpsOptions,
}

#[derive(Clone, Copy, ValueEnum)]
pub enum BioInputFormat {
    Pdb,
    Mmcif,
}

#[derive(Clone, Copy, ValueEnum)]
pub enum ChemInputFormat {
    Mol2,
    Sdf,
}

#[derive(Clone, Copy, ValueEnum)]
pub enum BioOutputFormat {
    /// LAMMPS data file (topology + coordinates)
    #[value(name = "lammps-data", alias = "data")]
    LammpsData,
    /// LAMMPS settings file (force field parameters)
    #[value(name = "lammps-settings", alias = "settings")]
    LammpsSettings,
    /// BGF format (DREIDING native)
    Bgf,
    /// PDB format
    Pdb,
    /// mmCIF format
    Mmcif,
    /// MOL2 format
    Mol2,
    /// SDF format
    Sdf,
}

#[derive(Clone, Copy, ValueEnum)]
pub enum ChemOutputFormat {
    /// LAMMPS data file
    #[value(name = "lammps-data", alias = "data")]
    LammpsData,
    /// LAMMPS settings file
    #[value(name = "lammps-settings", alias = "settings")]
    LammpsSettings,
    /// MOL2 format
    Mol2,
    /// SDF format
    Sdf,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum ChargeMethod {
    /// QEq charge equilibration
    #[default]
    Qeq,
    /// No charge calculation (all zeros + H-bond compensation)
    None,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum BondPotential {
    /// Harmonic bond stretching
    #[default]
    Harmonic,
    /// Morse anharmonic potential
    Morse,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum AnglePotential {
    /// Cosine-harmonic (DREIDING original)
    Cosine,
    /// Theta-harmonic (simple harmonic in angle)
    #[default]
    #[value(name = "theta-harmonic", alias = "theta")]
    ThetaHarmonic,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum VdwPotential {
    /// Lennard-Jones 12-6
    #[default]
    #[value(alias = "lennard-jones")]
    Lj,
    /// Exponential-6 (Buckingham)
    #[value(alias = "buckingham")]
    Exp6,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum SystemBoundary {
    /// Periodic boundary conditions (PPPM)
    #[default]
    Periodic,
    /// Non-periodic (shrink-wrapped)
    #[value(name = "non-periodic", alias = "shrink")]
    NonPeriodic,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum HisStrategy {
    /// Always HID (Nδ protonated)
    Hid,
    /// Always HIE (Nε protonated)
    Hie,
    /// Random selection
    Random,
    /// H-bond network analysis
    #[default]
    Network,
}

#[derive(Clone, Copy, Debug, ValueEnum, Default)]
pub enum Cation {
    #[default]
    Na,
    K,
    Mg,
    Ca,
    Li,
    Zn,
}

#[derive(Clone, Copy, Debug, ValueEnum, Default)]
pub enum Anion {
    #[default]
    Cl,
    Br,
    I,
    F,
}

pub fn parse() -> Cli {
    Cli::parse()
}
