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
//! - [`BondPotentialType`] — Bond stretching potential selection
//! - [`AnglePotentialType`] — Angle bending potential selection
//! - [`VdwPotentialType`] — Van der Waals potential selection

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
/// // Default configuration
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

/// Method for calculating partial atomic charges.
///
/// Determines how partial charges are assigned to atoms during
/// parameterization. The choice affects both the accuracy of
/// electrostatic interactions and hydrogen bond parameters.
#[derive(Debug, Clone, Default)]
pub enum ChargeMethod {
    /// No charge calculation; all charges remain zero.
    ///
    /// Use this for gas-phase calculations where electrostatics
    /// are not critical, or when charges will be assigned externally.
    #[default]
    None,

    /// Charge equilibration (QEq) method.
    ///
    /// Calculates electronegativity-equalized charges based on
    /// atomic positions and electronegativity parameters.
    Qeq(QeqConfig),
}

/// Configuration for QEq charge equilibration.
///
/// Controls the behavior of the charge equilibration solver,
/// including the target total charge and numerical solver options.
///
/// # Examples
///
/// ```
/// use dreid_forge::QeqConfig;
///
/// // Neutral molecule (default)
/// let neutral = QeqConfig::default();
/// assert_eq!(neutral.total_charge, 0.0);
///
/// // Negatively charged system
/// let anion = QeqConfig {
///     total_charge: -1.0,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct QeqConfig {
    /// Target total charge of the system in elementary charge units.
    ///
    /// The QEq solver will constrain the sum of all partial charges
    /// to equal this value. Default is `0.0` (neutral).
    pub total_charge: f64,

    /// QEq solver options.
    ///
    /// Controls convergence tolerance, iteration limits, and
    /// hydrogen SCF treatment.
    pub solver_options: SolverOptions,
}

impl Default for QeqConfig {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            solver_options: SolverOptions {
                hydrogen_scf: false,
                ..SolverOptions::default()
            },
        }
    }
}

pub use cheq::SolverOptions;

/// Bond stretching potential function type.
///
/// Determines the functional form used for bond stretching terms
/// in the force field. Both options use the same equilibrium bond
/// length but differ in behavior far from equilibrium.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondPotentialType {
    /// Harmonic bond potential.
    #[default]
    Harmonic,

    /// Morse anharmonic potential.
    Morse,
}

/// Angle bending potential function type.
///
/// Determines the functional form used for angle bending terms.
/// The choice affects behavior especially for linear or near-linear
/// angles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AnglePotentialType {
    /// Cosine-harmonic angle potential (DREIDING original).
    CosineHarmonic,

    /// Theta-harmonic angle potential.
    #[default]
    ThetaHarmonic,
}

/// Van der Waals non-bonded potential function type.
///
/// Determines the functional form used for van der Waals (dispersion
/// and repulsion) interactions between non-bonded atom pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VdwPotentialType {
    /// Lennard-Jones 12-6 potential.
    #[default]
    LennardJones,

    /// Exponential-6 (Buckingham) potential.
    Exponential6,
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

    #[test]
    fn qeq_config_default() {
        let qeq = QeqConfig::default();
        assert_eq!(qeq.total_charge, 0.0);
        assert_eq!(qeq.solver_options.tolerance, 1.0e-6);
    }

    #[test]
    fn charge_method_with_qeq() {
        let config = ForgeConfig {
            charge_method: ChargeMethod::Qeq(QeqConfig {
                total_charge: -1.0,
                ..Default::default()
            }),
            ..Default::default()
        };

        if let ChargeMethod::Qeq(qeq) = config.charge_method {
            assert_eq!(qeq.total_charge, -1.0);
        } else {
            panic!("Expected Qeq variant");
        }
    }
}
