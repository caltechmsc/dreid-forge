use dreid_forge::io::{
    Anion as LibAnion, Cation as LibCation, Format, HisStrategy as LibHisStrategy, SystemType,
};
use dreid_forge::{
    AnglePotentialType, BondPotentialType, ChargeMethod as LibChargeMethod, QeqConfig,
    VdwPotentialType,
};

use crate::cli;

pub fn build_charge_method(method: cli::ChargeMethod, total: f64) -> LibChargeMethod {
    match method {
        cli::ChargeMethod::Qeq => LibChargeMethod::Qeq(QeqConfig {
            total_charge: total,
            ..Default::default()
        }),
        cli::ChargeMethod::None => LibChargeMethod::None,
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
            cli::BioOutputFormat::LammpsData => Self::LammpsData,
            cli::BioOutputFormat::LammpsSettings => Self::LammpsSettings,
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
            cli::ChemOutputFormat::LammpsData => Self::LammpsData,
            cli::ChemOutputFormat::LammpsSettings => Self::LammpsSettings,
            cli::ChemOutputFormat::Mol2 => Self::Mol2,
            cli::ChemOutputFormat::Sdf => Self::Sdf,
        }
    }
}

impl From<cli::SystemBoundary> for SystemType {
    fn from(s: cli::SystemBoundary) -> Self {
        match s {
            cli::SystemBoundary::Periodic => Self::Periodic,
            cli::SystemBoundary::NonPeriodic => Self::NonPeriodic,
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
    forge: &cli::ForgeOptions,
) -> (&'static str, &'static str, &'static str) {
    let bond = match forge.bond_potential {
        cli::BondPotential::Harmonic => "Harmonic",
        cli::BondPotential::Morse => "Morse",
    };
    let angle = match forge.angle_potential {
        cli::AnglePotential::Cosine => "Cosine",
        cli::AnglePotential::ThetaHarmonic => "θ-Harm",
    };
    let vdw = match forge.vdw_potential {
        cli::VdwPotential::Lj => "LJ 12-6",
        cli::VdwPotential::Exp6 => "Exp-6",
    };
    (bond, angle, vdw)
}
