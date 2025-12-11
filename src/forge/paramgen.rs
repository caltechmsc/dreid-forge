use super::config::{
    AnglePotentialType, BondPotentialType, ChargeMethod, ForgeConfig, VdwPotentialType,
};
use super::error::Error;
use super::intermediate::{IntermediateSystem, PhysicalBondOrderExt};
use super::params::{ForceFieldParams, get_torsion_params, is_oxygen_column};
use crate::model::system::System;
use crate::model::topology::{
    AnglePotential, AtomParam, BondPotential, DihedralPotential, ForgedSystem, HBondPotential,
    ImproperPotential, Potentials, VdwPairPotential,
};
use crate::model::types::Element;
use std::collections::{HashMap, HashSet};

const H_BOND_DONOR_HYDROGEN_TYPE: &str = "H_HB";

fn collect_atom_types(intermediate: &IntermediateSystem) -> (Vec<String>, HashMap<String, usize>) {
    let mut types_set = HashSet::new();
    for atom in &intermediate.atoms {
        types_set.insert(atom.atom_type.clone());
    }

    let mut atom_types: Vec<String> = types_set.into_iter().collect();
    atom_types.sort();

    let type_indices: HashMap<String, usize> = atom_types
        .iter()
        .enumerate()
        .map(|(idx, t)| (t.clone(), idx))
        .collect();

    (atom_types, type_indices)
}

fn generate_atom_properties(
    intermediate: &IntermediateSystem,
    type_indices: &HashMap<String, usize>,
) -> Result<Vec<AtomParam>, Error> {
    intermediate
        .atoms
        .iter()
        .map(|atom| {
            let type_idx = *type_indices.get(&atom.atom_type).ok_or_else(|| {
                Error::Conversion(format!("missing type index for '{}'", atom.atom_type))
            })?;

            Ok(AtomParam {
                charge: atom.charge,
                mass: atom.element.atomic_mass(),
                type_idx,
            })
        })
        .collect()
}

fn generate_bond_potentials(
    intermediate: &IntermediateSystem,
    params: &ForceFieldParams,
    config: &ForgeConfig,
) -> Result<Vec<BondPotential>, Error> {
    intermediate
        .bonds
        .iter()
        .map(|bond| {
            let type_i = &intermediate.atoms[bond.i].atom_type;
            let type_j = &intermediate.atoms[bond.j].atom_type;

            let params_i = params
                .atoms
                .get(type_i)
                .ok_or_else(|| Error::missing_parameter(type_i, "bond radius"))?;
            let params_j = params
                .atoms
                .get(type_j)
                .ok_or_else(|| Error::missing_parameter(type_j, "bond radius"))?;

            let r0 = params_i.bond_radius + params_j.bond_radius - params.global.bond_delta;

            let physical_order = bond.physical_order.ok_or_else(|| {
                Error::Conversion(format!(
                    "bond {}-{} has no physical order assigned",
                    bond.i, bond.j
                ))
            })?;

            let order_mult = physical_order.multiplier();
            let k_force = params.global.bond_k * order_mult;
            let d0 = params.global.bond_d * order_mult;

            Ok(match config.bond_potential {
                BondPotentialType::Harmonic => BondPotential::Harmonic {
                    i: bond.i,
                    j: bond.j,
                    k_force,
                    r0,
                },
                BondPotentialType::Morse => {
                    let alpha = (k_force / (2.0 * d0)).sqrt();
                    BondPotential::Morse {
                        i: bond.i,
                        j: bond.j,
                        r0,
                        d0,
                        alpha,
                    }
                }
            })
        })
        .collect()
}

fn generate_angle_potentials(
    intermediate: &IntermediateSystem,
    params: &ForceFieldParams,
    config: &ForgeConfig,
) -> Result<Vec<AnglePotential>, Error> {
    intermediate
        .angles
        .iter()
        .map(|angle| {
            let type_j = &intermediate.atoms[angle.j].atom_type;
            let params_j = params
                .atoms
                .get(type_j)
                .ok_or_else(|| Error::missing_parameter(type_j, "bond angle"))?;

            let theta0_deg = params_j.bond_angle;
            let theta0_rad = theta0_deg.to_radians();
            let k_force = params.global.angle_k;

            Ok(match config.angle_potential {
                AnglePotentialType::ThetaHarmonic => AnglePotential::ThetaHarmonic {
                    i: angle.i,
                    j: angle.j,
                    k: angle.k,
                    k_force,
                    theta0: theta0_rad,
                },
                AnglePotentialType::CosineHarmonic => {
                    if (theta0_deg - 180.0).abs() < 1e-6 {
                        AnglePotential::CosineHarmonic {
                            i: angle.i,
                            j: angle.j,
                            k: angle.k,
                            k_force,
                            theta0: theta0_rad,
                        }
                    } else {
                        let sin_theta0 = theta0_rad.sin();
                        let c_force = k_force / (sin_theta0 * sin_theta0);
                        AnglePotential::CosineHarmonic {
                            i: angle.i,
                            j: angle.j,
                            k: angle.k,
                            k_force: c_force,
                            theta0: theta0_rad,
                        }
                    }
                }
            })
        })
        .collect()
}

fn generate_dihedral_potentials(
    intermediate: &IntermediateSystem,
    _params: &ForceFieldParams,
) -> Result<Vec<DihedralPotential>, Error> {
    let mut dihedrals = Vec::new();

    let mut bond_torsion_counts: HashMap<(usize, usize), usize> = HashMap::new();
    for dih in &intermediate.dihedrals {
        let key = (dih.j.min(dih.k), dih.j.max(dih.k));
        *bond_torsion_counts.entry(key).or_insert(0) += 1;
    }

    for dih in &intermediate.dihedrals {
        let atom_i = &intermediate.atoms[dih.i];
        let atom_j = &intermediate.atoms[dih.j];
        let atom_k = &intermediate.atoms[dih.k];

        let j_is_o_column = is_oxygen_column(atom_j.element);
        let k_is_o_column = is_oxygen_column(atom_k.element);

        if let Some(torsion) = get_torsion_params(
            atom_j.hybridization,
            atom_k.hybridization,
            j_is_o_column,
            k_is_o_column,
            atom_i.hybridization,
        ) {
            let key = (dih.j.min(dih.k), dih.j.max(dih.k));
            let count = *bond_torsion_counts.get(&key).unwrap_or(&1) as f64;
            let v_normalized = torsion.v_barrier / count;

            dihedrals.push(DihedralPotential {
                i: dih.i,
                j: dih.j,
                k: dih.k,
                l: dih.l,
                v_barrier: v_normalized,
                periodicity: torsion.periodicity,
                phase_offset: torsion.phase_offset,
            });
        }
    }

    Ok(dihedrals)
}

fn generate_improper_potentials(
    intermediate: &IntermediateSystem,
    params: &ForceFieldParams,
) -> Result<Vec<ImproperPotential>, Error> {
    let impropers = intermediate
        .impropers
        .iter()
        .map(|imp| ImproperPotential::Planar {
            i: imp.p1,
            j: imp.center,
            k: imp.p2,
            l: imp.p3,
            k_force: params.global.inversion_k,
            chi0: 0.0,
        })
        .collect();

    Ok(impropers)
}
