//! DREIDING force field parameter definitions and loading.
//!
//! This module handles loading and accessing DREIDING force field parameters.
//! Parameters are defined in TOML format and include atom type properties
//! (bond radii, angles, vdW parameters), global force constants, and
//! hydrogen bond settings.
//!
//! The default parameters are embedded in the library from
//! `resources/default.params.toml`. Custom parameters can be provided
//! via [`ForgeConfig::params`](super::ForgeConfig::params).

use super::error::Error;
use super::intermediate::Hybridization;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::OnceLock;

const DEFAULT_PARAMS_TOML: &str = include_str!("../../resources/default.params.toml");

static DEFAULT_PARAMS: OnceLock<ForceFieldParams> = OnceLock::new();

/// Complete set of DREIDING force field parameters.
///
/// Contains global force constants, per-atom-type parameters, and
/// hydrogen bond settings. Deserialized from TOML format.
#[derive(Debug, Clone, Deserialize)]
pub struct ForceFieldParams {
    /// Global force constants and tolerances.
    pub global: GlobalParams,
    /// Per-atom-type parameters indexed by type name (e.g., "C_3").
    #[serde(default)]
    pub atoms: HashMap<String, AtomTypeParams>,
    /// Hydrogen bond potential parameters.
    #[serde(default)]
    pub hydrogen_bond: HydrogenBondParams,
}

/// Global force field parameters.
///
/// These parameters apply across all atom types and control the
/// base force constants for bonded interactions.
#[derive(Debug, Clone, Deserialize)]
pub struct GlobalParams {
    /// Bond stretching force constant base (kcal/mol/Å²).
    #[serde(default = "default_bond_k")]
    pub bond_k: f64,
    /// Bond dissociation energy base for Morse potential (kcal/mol).
    #[serde(default = "default_bond_d")]
    pub bond_d: f64,
    /// Angle bending force constant (kcal/mol/rad²).
    #[serde(default = "default_angle_k")]
    pub angle_k: f64,
    /// Inversion (improper) force constant (kcal/mol/rad²).
    #[serde(default = "default_inversion_k")]
    pub inversion_k: f64,
    /// Bond length correction for covalent radii sum (Å).
    #[serde(default = "default_bond_delta")]
    pub bond_delta: f64,
}

fn default_bond_k() -> f64 {
    700.0
}
fn default_bond_d() -> f64 {
    70.0
}
fn default_angle_k() -> f64 {
    100.0
}
fn default_inversion_k() -> f64 {
    40.0
}
fn default_bond_delta() -> f64 {
    0.01
}

impl Default for GlobalParams {
    fn default() -> Self {
        Self {
            bond_k: default_bond_k(),
            bond_d: default_bond_d(),
            angle_k: default_angle_k(),
            inversion_k: default_inversion_k(),
            bond_delta: default_bond_delta(),
        }
    }
}

/// Per-atom-type force field parameters.
///
/// Contains geometric and van der Waals parameters specific to
/// each DREIDING atom type.
#[derive(Debug, Clone, Deserialize)]
pub struct AtomTypeParams {
    /// Covalent bond radius (Å).
    pub bond_radius: f64,
    /// Equilibrium bond angle (degrees).
    pub bond_angle: f64,
    /// Van der Waals equilibrium distance (Å).
    pub vdw_r0: f64,
    /// Van der Waals well depth (kcal/mol).
    pub vdw_d0: f64,
    /// Exponential-6 zeta parameter (dimensionless).
    #[serde(default = "default_vdw_zeta")]
    pub vdw_zeta: f64,
}

fn default_vdw_zeta() -> f64 {
    12.0
}

/// Hydrogen bond potential parameters.
///
/// Controls the directional hydrogen bond terms in DREIDING.
#[derive(Debug, Clone, Deserialize)]
pub struct HydrogenBondParams {
    /// Equilibrium H···A distance (Å).
    #[serde(default = "default_hbond_r0")]
    pub r0: f64,
    /// H-bond energy when no explicit charges (kcal/mol).
    #[serde(default = "default_hbond_d0_no_charge")]
    pub d0_no_charge: f64,
    /// H-bond energy with explicit charges (kcal/mol).
    #[serde(default = "default_hbond_d0_explicit")]
    pub d0_explicit: f64,
}

fn default_hbond_r0() -> f64 {
    2.75
}
fn default_hbond_d0_no_charge() -> f64 {
    9.0
}
fn default_hbond_d0_explicit() -> f64 {
    4.0
}

impl Default for HydrogenBondParams {
    fn default() -> Self {
        Self {
            r0: default_hbond_r0(),
            d0_no_charge: default_hbond_d0_no_charge(),
            d0_explicit: default_hbond_d0_explicit(),
        }
    }
}

/// Loads force field parameters from TOML.
///
/// If `custom_toml` is provided, parses it as a complete parameter set.
/// Otherwise, returns a clone of the embedded default parameters.
///
/// # Arguments
///
/// * `custom_toml` — Optional TOML string containing custom parameters
///
/// # Errors
///
/// Returns [`Error::ParameterParse`] if the TOML is malformed.
pub fn load_parameters(custom_toml: Option<&str>) -> Result<ForceFieldParams, Error> {
    match custom_toml {
        Some(toml) => {
            let params: ForceFieldParams = toml::from_str(toml)?;
            Ok(params)
        }
        None => Ok(get_default_parameters().clone()),
    }
}

/// Returns a reference to the embedded default parameters.
///
/// The parameters are parsed once and cached for subsequent calls.
pub fn get_default_parameters() -> &'static ForceFieldParams {
    DEFAULT_PARAMS.get_or_init(|| {
        toml::from_str(DEFAULT_PARAMS_TOML)
            .expect("Failed to parse embedded default parameters. This is a library bug.")
    })
}

/// Torsion potential parameters for a dihedral angle.
#[derive(Debug, Clone, Copy)]
pub struct TorsionParams {
    /// Barrier height (kcal/mol).
    pub v_barrier: f64,
    /// Periodicity (number of minima per 360°).
    pub periodicity: i32,
    /// Phase offset (degrees).
    pub phase_offset: f64,
}

/// Determines torsion parameters based on central bond hybridizations.
///
/// Implements the DREIDING torsion parameter rules based on the
/// hybridization states of the two central atoms (j and k) in the
/// i-j-k-l dihedral.
///
/// # Arguments
///
/// * `j_hyb` — Hybridization of atom j
/// * `k_hyb` — Hybridization of atom k
/// * `j_is_oxygen_column` — Whether atom j is O, S, Se, or Te
/// * `k_is_oxygen_column` — Whether atom k is O, S, Se, or Te
/// * `i_hyb` — Hybridization of atom i (for special cases)
///
/// # Returns
///
/// `Some(TorsionParams)` if a torsion term should be generated,
/// `None` if no torsion is defined for this combination.
pub fn get_torsion_params(
    j_hyb: Hybridization,
    k_hyb: Hybridization,
    j_is_oxygen_column: bool,
    k_is_oxygen_column: bool,
    i_hyb: Hybridization,
) -> Option<TorsionParams> {
    match (j_hyb, k_hyb) {
        (Hybridization::SP3, Hybridization::SP3) => {
            if j_is_oxygen_column && k_is_oxygen_column {
                Some(TorsionParams {
                    v_barrier: 2.0,
                    periodicity: 2,
                    phase_offset: 90.0,
                })
            } else {
                Some(TorsionParams {
                    v_barrier: 2.0,
                    periodicity: 3,
                    phase_offset: 180.0,
                })
            }
        }
        (Hybridization::SP2, Hybridization::SP2) => Some(TorsionParams {
            v_barrier: 45.0,
            periodicity: 2,
            phase_offset: 180.0,
        }),
        (Hybridization::Resonant, Hybridization::Resonant) => Some(TorsionParams {
            v_barrier: 25.0,
            periodicity: 2,
            phase_offset: 180.0,
        }),
        (Hybridization::SP2, Hybridization::Resonant)
        | (Hybridization::Resonant, Hybridization::SP2) => Some(TorsionParams {
            v_barrier: 5.0,
            periodicity: 2,
            phase_offset: 180.0,
        }),
        (Hybridization::SP2 | Hybridization::Resonant, Hybridization::SP3)
        | (Hybridization::SP3, Hybridization::SP2 | Hybridization::Resonant) => {
            let sp3_is_oxygen = if j_hyb == Hybridization::SP3 {
                j_is_oxygen_column
            } else {
                k_is_oxygen_column
            };

            if sp3_is_oxygen {
                Some(TorsionParams {
                    v_barrier: 2.0,
                    periodicity: 2,
                    phase_offset: 180.0,
                })
            } else if !matches!(i_hyb, Hybridization::SP2 | Hybridization::Resonant) {
                Some(TorsionParams {
                    v_barrier: 2.0,
                    periodicity: 3,
                    phase_offset: 180.0,
                })
            } else {
                Some(TorsionParams {
                    v_barrier: 1.0,
                    periodicity: 6,
                    phase_offset: 0.0,
                })
            }
        }
        _ => None,
    }
}

/// Checks if an element belongs to the oxygen column (group 16).
///
/// Returns `true` for O, S, Se, Te which have special torsion rules
/// in DREIDING.
pub fn is_oxygen_column(element: crate::model::types::Element) -> bool {
    use crate::model::types::Element;
    matches!(element, Element::O | Element::S | Element::Se | Element::Te)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::Element;

    #[test]
    fn default_parameters_load_common_types() {
        let params = get_default_parameters();
        assert!(params.atoms.contains_key("C_3"));
        assert!(params.atoms.contains_key("C_2"));
        assert!(params.atoms.contains_key("C_R"));
        assert!(params.atoms.contains_key("O_3"));
        assert!(params.atoms.contains_key("N_3"));
        assert!(params.atoms.contains_key("H_"));
        assert!(params.atoms.contains_key("H_HB"));
    }

    #[test]
    fn custom_parameters_parse_valid_toml() {
        let custom = r#"
            [global]
            bond_k = 800.0

            [atoms.C_3]
            bond_radius = 0.77
            bond_angle = 109.47
            vdw_r0 = 3.90
            vdw_d0 = 0.095
        "#;
        let params = load_parameters(Some(custom)).unwrap();
        assert_eq!(params.global.bond_k, 800.0);
        assert!(params.atoms.contains_key("C_3"));
    }

    #[test]
    fn errors_on_invalid_custom_toml() {
        let invalid = "not valid [[[toml";
        let result = load_parameters(Some(invalid));
        assert!(result.is_err());
    }

    #[test]
    fn global_params_default_values() {
        let global = GlobalParams::default();
        assert_eq!(global.bond_k, 700.0);
        assert_eq!(global.bond_d, 70.0);
        assert_eq!(global.angle_k, 100.0);
        assert_eq!(global.inversion_k, 40.0);
        assert_eq!(global.bond_delta, 0.01);
    }

    #[test]
    fn hydrogen_bond_params_default_values() {
        let hb = HydrogenBondParams::default();
        assert_eq!(hb.r0, 2.75);
        assert_eq!(hb.d0_no_charge, 9.0);
        assert_eq!(hb.d0_explicit, 4.0);
    }

    #[test]
    fn torsion_sp3_sp3_standard() {
        let params = get_torsion_params(
            Hybridization::SP3,
            Hybridization::SP3,
            false,
            false,
            Hybridization::SP3,
        );
        let p = params.expect("should have torsion params");
        assert_eq!(p.periodicity, 3);
        assert_eq!(p.v_barrier, 2.0);
        assert_eq!(p.phase_offset, 180.0);
    }

    #[test]
    fn torsion_sp2_sp2_double_bond() {
        let params = get_torsion_params(
            Hybridization::SP2,
            Hybridization::SP2,
            false,
            false,
            Hybridization::SP3,
        );
        let p = params.expect("should have torsion params");
        assert_eq!(p.periodicity, 2);
        assert_eq!(p.v_barrier, 45.0);
        assert_eq!(p.phase_offset, 180.0);
    }

    #[test]
    fn torsion_resonant_resonant() {
        let params = get_torsion_params(
            Hybridization::Resonant,
            Hybridization::Resonant,
            false,
            false,
            Hybridization::Resonant,
        );
        let p = params.expect("should have torsion params");
        assert_eq!(p.periodicity, 2);
        assert_eq!(p.v_barrier, 25.0);
    }

    #[test]
    fn torsion_sp2_resonant_mixed() {
        let params = get_torsion_params(
            Hybridization::SP2,
            Hybridization::Resonant,
            false,
            false,
            Hybridization::SP3,
        );
        let p = params.expect("should have torsion params");
        assert_eq!(p.v_barrier, 5.0);
        assert_eq!(p.periodicity, 2);
    }

    #[test]
    fn torsion_oxygen_column_sp3_sp3() {
        let params = get_torsion_params(
            Hybridization::SP3,
            Hybridization::SP3,
            true,
            true,
            Hybridization::SP3,
        );
        let p = params.expect("should have torsion params");
        assert_eq!(p.periodicity, 2);
        assert_eq!(p.phase_offset, 90.0);
    }

    #[test]
    fn torsion_none_for_sp1() {
        assert!(
            get_torsion_params(
                Hybridization::SP,
                Hybridization::SP3,
                false,
                false,
                Hybridization::SP3
            )
            .is_none()
        );
    }

    #[test]
    fn oxygen_column_detection() {
        assert!(is_oxygen_column(Element::O));
        assert!(is_oxygen_column(Element::S));
        assert!(is_oxygen_column(Element::Se));
        assert!(is_oxygen_column(Element::Te));
        assert!(!is_oxygen_column(Element::C));
        assert!(!is_oxygen_column(Element::N));
    }
}
