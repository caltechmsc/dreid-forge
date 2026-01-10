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

/// Charge calculation options shared by bio and chem commands.
#[derive(Args)]
#[command(next_help_heading = "Charge Calculation")]
pub struct ChargeOptions {
    /// Charge calculation method
    #[arg(long = "charge", value_name = "METHOD", default_value = "none")]
    pub method: ChargeMethod,

    /// Total system charge (for QEq constraint)
    #[arg(
        long = "total-charge",
        value_name = "Q",
        default_value = "0.0",
        allow_hyphen_values = true
    )]
    pub total_charge: f64,
}

/// Hybrid charge calculation options (bio command only).
#[derive(Args)]
#[command(next_help_heading = "Hybrid Charge Options")]
pub struct HybridChargeOptions {
    /// Protein force field charge scheme
    #[arg(
        long = "protein-scheme",
        value_name = "SCHEME",
        default_value = "amber-ffsb"
    )]
    pub protein_scheme: ProteinScheme,

    /// Nucleic acid force field charge scheme
    #[arg(
        long = "nucleic-scheme",
        value_name = "SCHEME",
        default_value = "amber"
    )]
    pub nucleic_scheme: NucleicScheme,

    /// Water model charge scheme
    #[arg(long = "water-scheme", value_name = "SCHEME", default_value = "tip3p")]
    pub water_scheme: WaterScheme,

    /// Ligand configuration (CHAIN:RESID[:ICODE][:METHOD[:CUTOFF]]), repeatable
    ///
    /// METHOD: vacuum | embedded (default: use --default-ligand-method)
    /// CUTOFF: environment radius in Å for embedded method
    #[arg(long = "ligand", value_name = "CONFIG", action = clap::ArgAction::Append)]
    pub ligands: Vec<String>,

    /// Default ligand QEq method for unlisted ligands
    #[arg(
        long = "default-ligand-method",
        value_name = "METHOD",
        default_value = "embedded"
    )]
    pub default_ligand_method: LigandQeqMethod,

    /// Default embedded QEq cutoff radius for unlisted ligands (Å)
    #[arg(
        long = "default-ligand-cutoff",
        value_name = "Å",
        default_value = "10.0"
    )]
    pub default_ligand_cutoff: f64,
}

/// QEq solver options (advanced tuning).
#[derive(Args)]
#[command(next_help_heading = "QEq Solver Options")]
pub struct QeqSolverOptions {
    /// Convergence tolerance for charge equilibration
    #[arg(long = "qeq-tolerance", value_name = "TOL", default_value = "1e-6")]
    pub tolerance: f64,

    /// Maximum iterations for QEq solver
    #[arg(long = "qeq-max-iter", value_name = "N", default_value = "100")]
    pub max_iterations: u32,

    /// Orbital screening parameter λ (Rappe–Goddard)
    #[arg(long = "qeq-lambda", value_name = "λ", default_value = "0.5")]
    pub lambda_scale: f64,

    /// Enable hydrogen SCF (nonlinear hardness update)
    #[arg(long = "qeq-hydrogen-scf", value_name = "BOOL", default_value = "true")]
    pub hydrogen_scf: bool,

    /// Basis function type for Coulomb integrals
    #[arg(long = "qeq-basis", value_name = "TYPE", default_value = "sto")]
    pub basis_type: BasisType,

    /// SCF damping strategy (none, fixed:<f>, auto, auto:<f>)
    #[arg(long = "qeq-damping", value_name = "STRATEGY", default_value = "auto")]
    pub damping: DampingStrategy,
}

/// Potential function options shared by bio and chem commands.
#[derive(Args)]
#[command(next_help_heading = "Potential Functions")]
pub struct PotentialOptions {
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
    pub charge: ChargeOptions,

    #[command(flatten)]
    pub hybrid: HybridChargeOptions,

    #[command(flatten)]
    pub qeq: QeqSolverOptions,

    #[command(flatten)]
    pub potential: PotentialOptions,

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
    pub charge: ChargeOptions,

    #[command(flatten)]
    pub qeq: QeqSolverOptions,

    #[command(flatten)]
    pub potential: PotentialOptions,

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
    /// No charge calculation (all zeros + H-bond compensation)
    #[default]
    None,
    /// QEq charge equilibration for all atoms
    Qeq,
    /// Hybrid: force field for biomolecules, QEq for ligands
    Hybrid,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum ProteinScheme {
    /// AMBER ff99SB/ff14SB/ff19SB
    #[default]
    #[value(name = "amber-ffsb", alias = "ffsb")]
    AmberFfsb,
    /// AMBER ff03
    #[value(name = "amber-ff03", alias = "ff03")]
    AmberFf03,
    /// CHARMM22/27/36/36m
    Charmm,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum NucleicScheme {
    /// AMBER OL15/OL21/OL24/bsc1/OL3
    #[default]
    Amber,
    /// CHARMM C27/C36
    Charmm,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum WaterScheme {
    /// TIP3P
    #[default]
    Tip3p,
    /// TIP3P-FB
    #[value(name = "tip3p-fb")]
    Tip3pFb,
    /// SPC
    Spc,
    /// SPC/E
    #[value(name = "spc-e")]
    SpcE,
    /// OPC3
    Opc3,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum LigandQeqMethod {
    /// Vacuum QEq (isolated ligand)
    Vacuum,
    /// Embedded QEq (polarized by environment)
    #[default]
    Embedded,
}

#[derive(Clone, Copy, ValueEnum, Default)]
pub enum BasisType {
    /// Gaussian-type orbitals (faster, approximate)
    Gto,
    /// Slater-type orbitals (exact)
    #[default]
    Sto,
}

#[derive(Clone, Debug)]
pub enum DampingStrategy {
    /// No damping (fastest, may not converge)
    None,
    /// Fixed damping factor (0 < d ≤ 1)
    Fixed(f64),
    /// Automatic adaptive damping
    Auto {
        /// Initial damping factor
        initial: f64,
    },
}

impl Default for DampingStrategy {
    fn default() -> Self {
        Self::Auto { initial: 0.5 }
    }
}

impl std::str::FromStr for DampingStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_lowercase();
        if s == "none" {
            Ok(Self::None)
        } else if s == "auto" {
            Ok(Self::Auto { initial: 0.5 })
        } else if let Some(val) = s.strip_prefix("auto:") {
            let initial = val
                .parse::<f64>()
                .map_err(|_| format!("invalid auto damping value: {}", val))?;
            if initial <= 0.0 || initial > 1.0 {
                return Err("auto damping initial must be in (0, 1]".into());
            }
            Ok(Self::Auto { initial })
        } else if let Some(val) = s.strip_prefix("fixed:") {
            let factor = val
                .parse::<f64>()
                .map_err(|_| format!("invalid fixed damping value: {}", val))?;
            if factor <= 0.0 || factor > 1.0 {
                return Err("fixed damping must be in (0, 1]".into());
            }
            Ok(Self::Fixed(factor))
        } else {
            Err(format!(
                "unknown damping strategy: '{}' (use none, auto, auto:<f>, or fixed:<f>)",
                s
            ))
        }
    }
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
