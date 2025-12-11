use super::error::Error;
use super::intermediate::Hybridization;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::OnceLock;

const DEFAULT_PARAMS_TOML: &str = include_str!("../../resources/default.params.toml");

static DEFAULT_PARAMS: OnceLock<ForceFieldParams> = OnceLock::new();

#[derive(Debug, Clone, Deserialize)]
pub struct ForceFieldParams {
    pub global: GlobalParams,
    #[serde(default)]
    pub atoms: HashMap<String, AtomTypeParams>,
    #[serde(default)]
    pub hydrogen_bond: HydrogenBondParams,
}

#[derive(Debug, Clone, Deserialize)]
pub struct GlobalParams {
    #[serde(default = "default_bond_k")]
    pub bond_k: f64,
    #[serde(default = "default_bond_d")]
    pub bond_d: f64,
    #[serde(default = "default_angle_k")]
    pub angle_k: f64,
    #[serde(default = "default_inversion_k")]
    pub inversion_k: f64,
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

#[derive(Debug, Clone, Deserialize)]
pub struct AtomTypeParams {
    pub bond_radius: f64,
    pub bond_angle: f64,
    pub vdw_r0: f64,
    pub vdw_d0: f64,
    #[serde(default = "default_vdw_zeta")]
    pub vdw_zeta: f64,
}

fn default_vdw_zeta() -> f64 {
    12.0
}

#[derive(Debug, Clone, Deserialize)]
pub struct HydrogenBondParams {
    #[serde(default = "default_hbond_r0")]
    pub r0: f64,
    #[serde(default = "default_hbond_d0_no_charge")]
    pub d0_no_charge: f64,
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

pub fn load_parameters(custom_toml: Option<&str>) -> Result<&'static ForceFieldParams, Error> {
    match custom_toml {
        Some(toml) => {
            let params: ForceFieldParams = toml::from_str(toml)?;
            Ok(Box::leak(Box::new(params)))
        }
        None => Ok(get_default_parameters()),
    }
}

pub fn get_default_parameters() -> &'static ForceFieldParams {
    DEFAULT_PARAMS.get_or_init(|| {
        toml::from_str(DEFAULT_PARAMS_TOML)
            .expect("Failed to parse embedded default parameters. This is a library bug.")
    })
}

#[derive(Debug, Clone, Copy)]
pub struct TorsionParams {
    pub v_barrier: f64,
    pub periodicity: i32,
    pub phase_offset: f64,
}

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

pub fn is_oxygen_column(element: crate::model::types::Element) -> bool {
    use crate::model::types::Element;
    matches!(element, Element::O | Element::S | Element::Se | Element::Te)
}
