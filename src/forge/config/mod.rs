//! Configuration types for DREIDING parameterization.
//!
//! This module defines all configuration structures used to control the
//! behavior of the [`forge`](super::forge) function. Settings include
//! potential function types, charge calculation methods, and custom
//! parameter file paths.
//!
//! # Overview
//!
//! - [`ForgeConfig`] — Main configuration struct
//! - [`ChargeMethod`] — Charge calculation method selection
//! - [`QeqConfig`] — QEq charge equilibration settings
//! - [`HybridConfig`] — Hybrid biological/QEq charge assignment
//! - [`BondPotentialType`] — Bond stretching potential selection
//! - [`AnglePotentialType`] — Angle bending potential selection
//! - [`VdwPotentialType`] — Van der Waals potential selection

mod charge;
mod potential;

pub use charge::{
    BasisType, ChargeMethod, DampingStrategy, EmbeddedQeqConfig, HybridConfig, LigandChargeConfig,
    LigandQeqMethod, NucleicScheme, ProteinScheme, QeqConfig, ResidueSelector, SolverOptions,
    WaterScheme,
};
pub use potential::{AnglePotentialType, BondPotentialType, VdwPotentialType};

/// Main configuration for DREIDING force field parameterization.
///
/// Controls all aspects of the parameterization process, including
/// potential function types, charge calculation method, and optional
/// custom parameter files.
///
/// # Examples
///
/// ```
/// use dreid_forge::{ForgeConfig, ChargeMethod, QeqConfig, BondPotentialType};
///
/// // Default configuration (no charges)
/// let default = ForgeConfig::default();
///
/// // Custom configuration with QEq charges and Morse bonds
/// let custom = ForgeConfig {
///     charge_method: ChargeMethod::Qeq(QeqConfig::default()),
///     bond_potential: BondPotentialType::Morse,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct ForgeConfig {
    /// Custom atom typing rules in TOML format.
    ///
    /// If `None`, uses the built-in DREIDING typing rules from `dreid-typer`.
    pub rules: Option<String>,

    /// Custom force field parameters in TOML format.
    ///
    /// If `None`, uses the embedded `default.params.toml` with standard
    /// DREIDING parameters.
    pub params: Option<String>,

    /// Method for calculating partial atomic charges.
    pub charge_method: ChargeMethod,

    /// Type of bond stretching potential to generate.
    pub bond_potential: BondPotentialType,

    /// Type of angle bending potential to generate.
    pub angle_potential: AnglePotentialType,

    /// Type of van der Waals non-bonded potential to generate.
    pub vdw_potential: VdwPotentialType,
}

impl Default for ForgeConfig {
    fn default() -> Self {
        Self {
            rules: None,
            params: None,
            charge_method: ChargeMethod::None,
            bond_potential: BondPotentialType::Harmonic,
            angle_potential: AnglePotentialType::ThetaHarmonic,
            vdw_potential: VdwPotentialType::LennardJones,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_values() {
        let config = ForgeConfig::default();
        assert!(config.rules.is_none());
        assert!(config.params.is_none());
        assert!(matches!(config.charge_method, ChargeMethod::None));
        assert_eq!(config.bond_potential, BondPotentialType::Harmonic);
        assert_eq!(config.angle_potential, AnglePotentialType::ThetaHarmonic);
        assert_eq!(config.vdw_potential, VdwPotentialType::LennardJones);
    }
}
