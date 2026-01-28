use dreid_forge::io::{
    Anion as LibAnion, Cation as LibCation, Format, HisStrategy as LibHisStrategy,
};
use dreid_forge::{
    AnglePotentialType, BondPotentialType, ChargeMethod as LibChargeMethod, EmbeddedQeqConfig,
    HybridConfig, LigandChargeConfig, LigandQeqMethod as LibLigandQeqMethod, QeqConfig,
    ResidueSelector, SolverOptions, VdwPotentialType,
};
use dreid_forge::{BasisType as LibBasisType, DampingStrategy as LibDampingStrategy};
use dreid_forge::{
    NucleicScheme as LibNucleicScheme, ProteinScheme as LibProteinScheme,
    WaterScheme as LibWaterScheme,
};

use crate::cli;

pub fn build_chem_charge_method(
    charge: &cli::ChargeOptions,
    qeq: &cli::QeqSolverOptions,
) -> LibChargeMethod {
    match charge.method {
        cli::ChargeMethod::None => LibChargeMethod::None,
        cli::ChargeMethod::Qeq | cli::ChargeMethod::Hybrid => LibChargeMethod::Qeq(QeqConfig {
            total_charge: charge.total_charge,
            solver_options: build_solver_options(qeq),
        }),
    }
}

pub fn build_bio_charge_method(
    charge: &cli::ChargeOptions,
    hybrid: &cli::HybridChargeOptions,
    qeq: &cli::QeqSolverOptions,
) -> LibChargeMethod {
    match charge.method {
        cli::ChargeMethod::None => LibChargeMethod::None,
        cli::ChargeMethod::Qeq => LibChargeMethod::Qeq(QeqConfig {
            total_charge: charge.total_charge,
            solver_options: build_solver_options(qeq),
        }),
        cli::ChargeMethod::Hybrid => {
            LibChargeMethod::Hybrid(build_hybrid_config(charge.total_charge, hybrid, qeq))
        }
    }
}

fn build_hybrid_config(
    total_charge: f64,
    hybrid: &cli::HybridChargeOptions,
    qeq: &cli::QeqSolverOptions,
) -> HybridConfig {
    let solver_options = build_solver_options(qeq);

    let ligand_configs: Vec<LigandChargeConfig> = hybrid
        .ligands
        .iter()
        .filter_map(|s| parse_ligand_config(s, &solver_options, total_charge, hybrid))
        .collect();

    let default_ligand_method = build_default_ligand_method(hybrid, &solver_options, total_charge);

    HybridConfig {
        protein_scheme: hybrid.protein_scheme.into(),
        nucleic_scheme: hybrid.nucleic_scheme.into(),
        water_scheme: hybrid.water_scheme.into(),
        ligand_configs,
        default_ligand_method,
    }
}

fn build_default_ligand_method(
    hybrid: &cli::HybridChargeOptions,
    solver_options: &SolverOptions,
    total_charge: f64,
) -> LibLigandQeqMethod {
    let qeq = QeqConfig {
        total_charge,
        solver_options: *solver_options,
    };

    match hybrid.default_ligand_method {
        cli::LigandQeqMethod::Vacuum => LibLigandQeqMethod::Vacuum(qeq),
        cli::LigandQeqMethod::Embedded => LibLigandQeqMethod::Embedded(EmbeddedQeqConfig {
            cutoff_radius: hybrid.default_ligand_cutoff,
            qeq,
        }),
    }
}

fn parse_ligand_config(
    s: &str,
    solver_options: &SolverOptions,
    total_charge: f64,
    hybrid: &cli::HybridChargeOptions,
) -> Option<LigandChargeConfig> {
    let parts: Vec<&str> = s.split(':').collect();

    if parts.len() < 2 {
        return None;
    }

    let chain = parts[0];
    let res_num: i32 = parts[1].parse().ok()?;

    let (icode, method) = match parts.len() {
        // CHAIN:RESID
        2 => (None, None),
        // CHAIN:RESID:X - X could be icode or method
        3 => {
            let x = parts[2];
            if is_method_keyword(x) {
                (
                    None,
                    Some(parse_ligand_method(
                        x,
                        None,
                        solver_options,
                        total_charge,
                        hybrid,
                    )?),
                )
            } else {
                (Some(x.chars().next()?), None)
            }
        }
        // CHAIN:RESID:X:Y - X could be icode, Y could be method or cutoff
        4 => {
            let x = parts[2];
            let y = parts[3];
            if is_method_keyword(x) {
                // X is method, Y is cutoff
                let cutoff: f64 = y.parse().ok()?;
                (
                    None,
                    Some(parse_ligand_method(
                        x,
                        Some(cutoff),
                        solver_options,
                        total_charge,
                        hybrid,
                    )?),
                )
            } else {
                // X is icode, Y is method
                let icode = x.chars().next()?;
                (
                    Some(icode),
                    Some(parse_ligand_method(
                        y,
                        None,
                        solver_options,
                        total_charge,
                        hybrid,
                    )?),
                )
            }
        }
        // CHAIN:RESID:ICODE:METHOD:CUTOFF (5 parts)
        5 => {
            let icode = parts[2].chars().next()?;
            let method_str = parts[3];
            let cutoff: f64 = parts[4].parse().ok()?;
            (
                Some(icode),
                Some(parse_ligand_method(
                    method_str,
                    Some(cutoff),
                    solver_options,
                    total_charge,
                    hybrid,
                )?),
            )
        }
        _ => return None,
    };

    let selector = ResidueSelector::new(chain, res_num, icode);

    let method =
        method.unwrap_or_else(|| build_default_ligand_method(hybrid, solver_options, total_charge));

    Some(LigandChargeConfig { selector, method })
}

fn is_method_keyword(s: &str) -> bool {
    let lower = s.to_lowercase();
    lower == "vacuum" || lower == "embedded"
}

fn parse_ligand_method(
    method_str: &str,
    cutoff: Option<f64>,
    solver_options: &SolverOptions,
    total_charge: f64,
    hybrid: &cli::HybridChargeOptions,
) -> Option<LibLigandQeqMethod> {
    let qeq = QeqConfig {
        total_charge,
        solver_options: *solver_options,
    };

    match method_str.to_lowercase().as_str() {
        "vacuum" => Some(LibLigandQeqMethod::Vacuum(qeq)),
        "embedded" => {
            let cutoff_radius = cutoff.unwrap_or(hybrid.default_ligand_cutoff);
            Some(LibLigandQeqMethod::Embedded(EmbeddedQeqConfig {
                cutoff_radius,
                qeq,
            }))
        }
        _ => None,
    }
}

fn build_solver_options(qeq: &cli::QeqSolverOptions) -> SolverOptions {
    SolverOptions {
        tolerance: qeq.tolerance,
        max_iterations: qeq.max_iterations,
        lambda_scale: qeq.lambda_scale,
        hydrogen_scf: qeq.hydrogen_scf,
        basis_type: qeq.basis_type.into(),
        damping: (&qeq.damping).into(),
    }
}

impl From<cli::ProteinScheme> for LibProteinScheme {
    fn from(s: cli::ProteinScheme) -> Self {
        match s {
            cli::ProteinScheme::AmberFfsb => Self::AmberFFSB,
            cli::ProteinScheme::AmberFf03 => Self::AmberFF03,
            cli::ProteinScheme::Charmm => Self::Charmm,
        }
    }
}

impl From<cli::NucleicScheme> for LibNucleicScheme {
    fn from(s: cli::NucleicScheme) -> Self {
        match s {
            cli::NucleicScheme::Amber => Self::Amber,
            cli::NucleicScheme::Charmm => Self::Charmm,
        }
    }
}

impl From<cli::WaterScheme> for LibWaterScheme {
    fn from(s: cli::WaterScheme) -> Self {
        match s {
            cli::WaterScheme::Tip3p => Self::Tip3p,
            cli::WaterScheme::Tip3pFb => Self::Tip3pFb,
            cli::WaterScheme::Spc => Self::Spc,
            cli::WaterScheme::SpcE => Self::SpcE,
            cli::WaterScheme::Opc3 => Self::Opc3,
        }
    }
}

impl From<cli::BasisType> for LibBasisType {
    fn from(b: cli::BasisType) -> Self {
        match b {
            cli::BasisType::Gto => Self::Gto,
            cli::BasisType::Sto => Self::Sto,
        }
    }
}

impl From<&cli::DampingStrategy> for LibDampingStrategy {
    fn from(d: &cli::DampingStrategy) -> Self {
        match d {
            cli::DampingStrategy::None => Self::None,
            cli::DampingStrategy::Fixed(f) => Self::Fixed(*f),
            cli::DampingStrategy::Auto { initial } => Self::Auto { initial: *initial },
        }
    }
}

impl From<cli::BondPotential> for BondPotentialType {
    fn from(p: cli::BondPotential) -> Self {
        match p {
            cli::BondPotential::Harmonic => Self::Harmonic,
            cli::BondPotential::Morse => Self::Morse,
        }
    }
}

impl From<cli::AnglePotential> for AnglePotentialType {
    fn from(p: cli::AnglePotential) -> Self {
        match p {
            cli::AnglePotential::Cosine => Self::CosineHarmonic,
            cli::AnglePotential::ThetaHarmonic => Self::ThetaHarmonic,
        }
    }
}

impl From<cli::VdwPotential> for VdwPotentialType {
    fn from(p: cli::VdwPotential) -> Self {
        match p {
            cli::VdwPotential::Lj => Self::LennardJones,
            cli::VdwPotential::Exp6 => Self::Exponential6,
        }
    }
}

impl From<cli::BioInputFormat> for Format {
    fn from(f: cli::BioInputFormat) -> Self {
        match f {
            cli::BioInputFormat::Pdb => Self::Pdb,
            cli::BioInputFormat::Mmcif => Self::Mmcif,
        }
    }
}

impl From<cli::ChemInputFormat> for Format {
    fn from(f: cli::ChemInputFormat) -> Self {
        match f {
            cli::ChemInputFormat::Mol2 => Self::Mol2,
            cli::ChemInputFormat::Sdf => Self::Sdf,
        }
    }
}

impl From<cli::BioOutputFormat> for Format {
    fn from(f: cli::BioOutputFormat) -> Self {
        match f {
            cli::BioOutputFormat::Bgf => Self::Bgf,
            cli::BioOutputFormat::Pdb => Self::Pdb,
            cli::BioOutputFormat::Mmcif => Self::Mmcif,
            cli::BioOutputFormat::Mol2 => Self::Mol2,
            cli::BioOutputFormat::Sdf => Self::Sdf,
        }
    }
}

impl From<cli::ChemOutputFormat> for Format {
    fn from(f: cli::ChemOutputFormat) -> Self {
        match f {
            cli::ChemOutputFormat::Mol2 => Self::Mol2,
            cli::ChemOutputFormat::Sdf => Self::Sdf,
        }
    }
}

impl From<cli::HisStrategy> for LibHisStrategy {
    fn from(s: cli::HisStrategy) -> Self {
        match s {
            cli::HisStrategy::Hid => Self::DirectHID,
            cli::HisStrategy::Hie => Self::DirectHIE,
            cli::HisStrategy::Random => Self::Random,
            cli::HisStrategy::Network => Self::HbNetwork,
        }
    }
}

impl From<cli::Cation> for LibCation {
    fn from(c: cli::Cation) -> Self {
        match c {
            cli::Cation::Na => Self::Na,
            cli::Cation::K => Self::K,
            cli::Cation::Mg => Self::Mg,
            cli::Cation::Ca => Self::Ca,
            cli::Cation::Li => Self::Li,
            cli::Cation::Zn => Self::Zn,
        }
    }
}

impl From<cli::Anion> for LibAnion {
    fn from(a: cli::Anion) -> Self {
        match a {
            cli::Anion::Cl => Self::Cl,
            cli::Anion::Br => Self::Br,
            cli::Anion::I => Self::I,
            cli::Anion::F => Self::F,
        }
    }
}

pub fn potential_display_names(
    potential: &cli::PotentialOptions,
) -> (&'static str, &'static str, &'static str) {
    let bond = match potential.bond_potential {
        cli::BondPotential::Harmonic => "Harmonic",
        cli::BondPotential::Morse => "Morse",
    };
    let angle = match potential.angle_potential {
        cli::AnglePotential::Cosine => "Cosine",
        cli::AnglePotential::ThetaHarmonic => "θ-Harm",
    };
    let vdw = match potential.vdw_potential {
        cli::VdwPotential::Lj => "LJ 12-6",
        cli::VdwPotential::Exp6 => "Exp-6",
    };
    (bond, angle, vdw)
}

pub fn protein_scheme_display_name(scheme: cli::ProteinScheme) -> &'static str {
    match scheme {
        cli::ProteinScheme::AmberFfsb => "AMBER-FFSB",
        cli::ProteinScheme::AmberFf03 => "AMBER-FF03",
        cli::ProteinScheme::Charmm => "CHARMM",
    }
}

pub fn nucleic_scheme_display_name(scheme: cli::NucleicScheme) -> &'static str {
    match scheme {
        cli::NucleicScheme::Amber => "AMBER",
        cli::NucleicScheme::Charmm => "CHARMM",
    }
}

pub fn water_scheme_display_name(scheme: cli::WaterScheme) -> &'static str {
    match scheme {
        cli::WaterScheme::Tip3p => "TIP3P",
        cli::WaterScheme::Tip3pFb => "TIP3P-FB",
        cli::WaterScheme::Spc => "SPC",
        cli::WaterScheme::SpcE => "SPC/E",
        cli::WaterScheme::Opc3 => "OPC3",
    }
}

pub fn ligand_method_display_name(method: cli::LigandQeqMethod) -> &'static str {
    match method {
        cli::LigandQeqMethod::Vacuum => "Vacuum",
        cli::LigandQeqMethod::Embedded => "Embedded",
    }
}
