//! Force field parameter generation.
//!
//! This module generates all potential energy function parameters from
//! the typed and charged intermediate system. It produces bond, angle,
//! torsion, inversion, van der Waals, and hydrogen bond potentials
//! according to DREIDING force field rules.

use super::config::{
    AnglePotentialType, BondPotentialType, ChargeMethod, ForgeConfig, VdwPotentialType,
};
use super::error::Error;
use super::intermediate::{IntermediateSystem, PhysicalBondOrderExt};
use super::params::{ForceFieldParams, get_torsion_params, is_oxygen_column};
use crate::model::system::System;
use crate::model::topology::{
    AnglePotential, AtomParam, BondPotential, ForgedSystem, HBondPotential, InversionPotential,
    Potentials, TorsionPotential, VdwPairPotential,
};
use crate::model::types::Element;
use dreid_kernel::potentials::bonded::{
    CosineHarmonic, Harmonic, PlanarInversion, ThetaHarmonic, Torsion,
};
use dreid_kernel::potentials::nonbonded::{Buckingham, HydrogenBond, LennardJones};
use std::collections::{HashMap, HashSet};

/// DREIDING atom type for hydrogen bond donor hydrogens.
const H_BOND_DONOR_HYDROGEN_TYPE: &str = "H_HB";

/// Generates all force field parameters from the intermediate system.
///
/// This is the final stage of the parameterization pipeline. Takes the
/// typed and charged intermediate system and produces all potential
/// energy function parameters needed for simulation.
///
/// # Arguments
///
/// * `system` — Original molecular system (cloned into output)
/// * `intermediate` — Typed and charged intermediate representation
/// * `params` — Force field parameter set
/// * `config` — Configuration controlling potential types
///
/// # Returns
///
/// A fully parameterized [`ForgedSystem`].
///
/// # Errors
///
/// Returns [`Error::MissingParameter`] if an atom type lacks required
/// parameters, or [`Error::Conversion`] for internal consistency failures.
pub fn generate_parameters(
    system: &System,
    intermediate: &IntermediateSystem,
    params: &ForceFieldParams,
    config: &ForgeConfig,
) -> Result<ForgedSystem, Error> {
    let (atom_types, type_indices) = collect_atom_types(intermediate);

    let atom_properties = generate_atom_properties(intermediate, &type_indices)?;

    let bonds = generate_bond_potentials(intermediate, params, config)?;

    let angles = generate_angle_potentials(intermediate, params, config)?;

    let torsions = generate_torsion_potentials(intermediate, params)?;

    let inversions = generate_inversion_potentials(intermediate, params)?;

    let vdw_pairs = generate_vdw_potentials(&atom_types, params, config)?;

    let h_bonds = generate_hbond_potentials(intermediate, &type_indices, params, config)?;

    let potentials = Potentials {
        bonds,
        angles,
        torsions,
        inversions,
        vdw_pairs,
        h_bonds,
    };

    Ok(ForgedSystem {
        system: system.clone(),
        atom_types,
        atom_properties,
        potentials,
    })
}

/// Collects unique atom types and builds a type-to-index mapping.
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

/// Generates per-atom properties (charge, mass, type index).
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

/// Generates bond stretching potentials (Harmonic or Morse).
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
            let k = params.global.bond_k * order_mult;
            let de = params.global.bond_d * order_mult;
            let atoms = (bond.i, bond.j);

            Ok(match config.bond_potential {
                BondPotentialType::Harmonic => {
                    let (k_half, r0) = Harmonic::precompute(k, r0);
                    BondPotential::Harmonic { atoms, k_half, r0 }
                }
                BondPotentialType::Morse => {
                    let alpha = (k / (2.0 * de)).sqrt();
                    BondPotential::Morse {
                        atoms,
                        de,
                        r0,
                        alpha,
                    }
                }
            })
        })
        .collect()
}

/// Generates angle bending potentials (Cosine-harmonic/linear or Theta-harmonic).
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
            let k = params.global.angle_k;
            let atoms = (angle.i, angle.j, angle.k);

            Ok(match config.angle_potential {
                AnglePotentialType::Theta => {
                    let (k_half, theta0) = ThetaHarmonic::precompute(k, theta0_deg);
                    AnglePotential::ThetaHarmonic {
                        atoms,
                        k_half,
                        theta0,
                    }
                }
                AnglePotentialType::Cosine => {
                    if (theta0_deg - 180.0).abs() < 1e-6 {
                        AnglePotential::CosineLinear { atoms, k }
                    } else {
                        let sin_theta0 = theta0_deg.to_radians().sin();
                        let c = k / (sin_theta0 * sin_theta0);
                        let (k_half, cos0) = CosineHarmonic::precompute(c, theta0_deg);
                        AnglePotential::CosineHarmonic {
                            atoms,
                            k_half,
                            cos0,
                        }
                    }
                }
            })
        })
        .collect()
}

/// Generates torsion potentials based on hybridization.
fn generate_torsion_potentials(
    intermediate: &IntermediateSystem,
    _params: &ForceFieldParams,
) -> Result<Vec<TorsionPotential>, Error> {
    let mut torsions = Vec::new();

    let mut bond_torsion_counts: HashMap<(usize, usize), usize> = HashMap::new();
    for tor in &intermediate.torsions {
        let key = (tor.j.min(tor.k), tor.j.max(tor.k));
        *bond_torsion_counts.entry(key).or_insert(0) += 1;
    }

    for tor in &intermediate.torsions {
        let atom_i = &intermediate.atoms[tor.i];
        let atom_j = &intermediate.atoms[tor.j];
        let atom_k = &intermediate.atoms[tor.k];

        let j_is_o_column = is_oxygen_column(atom_j.element);
        let k_is_o_column = is_oxygen_column(atom_k.element);

        if let Some(torsion_params) = get_torsion_params(
            atom_j.hybridization,
            atom_k.hybridization,
            j_is_o_column,
            k_is_o_column,
            atom_i.hybridization,
        ) {
            let key = (tor.j.min(tor.k), tor.j.max(tor.k));
            let count = *bond_torsion_counts.get(&key).unwrap_or(&1) as f64;
            let v_barrier = torsion_params.v_barrier / count;

            let (v_half, n, cos_n_phi0, sin_n_phi0) = Torsion::precompute(
                v_barrier,
                torsion_params.periodicity,
                torsion_params.phase_offset,
            );

            torsions.push(TorsionPotential {
                atoms: (tor.i, tor.j, tor.k, tor.l),
                v_half,
                n,
                cos_n_phi0,
                sin_n_phi0,
            });
        }
    }

    Ok(torsions)
}

/// Generates inversion (out-of-plane) potentials for sp² and resonant centers.
fn generate_inversion_potentials(
    intermediate: &IntermediateSystem,
    params: &ForceFieldParams,
) -> Result<Vec<InversionPotential>, Error> {
    let inversions = intermediate
        .inversions
        .iter()
        .map(|inv| {
            let atoms = (inv.center, inv.axis, inv.plane1, inv.plane2);
            // Scale by 1/3: each sp² center generates 3 terms
            let k_scaled = params.global.inversion_k / 3.0;
            let c_half = PlanarInversion::precompute(k_scaled);
            InversionPotential::Planar { atoms, c_half }
        })
        .collect();

    Ok(inversions)
}

/// Generates van der Waals pair potentials (Lennard-Jones or Buckingham).
fn generate_vdw_potentials(
    atom_types: &[String],
    params: &ForceFieldParams,
    config: &ForgeConfig,
) -> Result<Vec<VdwPairPotential>, Error> {
    let mut vdw_pairs = Vec::new();

    for (idx1, type1) in atom_types.iter().enumerate() {
        for (idx2, type2) in atom_types.iter().enumerate().skip(idx1) {
            let params1 = params
                .atoms
                .get(type1)
                .ok_or_else(|| Error::missing_parameter(type1, "vdW parameters"))?;
            let params2 = params
                .atoms
                .get(type2)
                .ok_or_else(|| Error::missing_parameter(type2, "vdW parameters"))?;

            let r0_combined = 0.5 * (params1.vdw_r0 + params2.vdw_r0);
            let d0_combined = (params1.vdw_d0 * params2.vdw_d0).sqrt();

            vdw_pairs.push(match config.vdw_potential {
                VdwPotentialType::LennardJones => {
                    let (d0, r0_sq) = LennardJones::precompute(d0_combined, r0_combined);
                    VdwPairPotential::LennardJones {
                        type1_idx: idx1,
                        type2_idx: idx2,
                        d0,
                        r0_sq,
                    }
                }
                VdwPotentialType::Buckingham => {
                    let zeta = 0.5 * (params1.vdw_zeta + params2.vdw_zeta);
                    let (a, b, c, r_max_sq, two_e_max) =
                        Buckingham::precompute(d0_combined, r0_combined, zeta);
                    VdwPairPotential::Buckingham {
                        type1_idx: idx1,
                        type2_idx: idx2,
                        a,
                        b,
                        c,
                        r_max_sq,
                        two_e_max,
                    }
                }
            });
        }
    }

    Ok(vdw_pairs)
}

/// Generates directional hydrogen bond potentials.
///
/// Creates H-bond terms for all unique donor-hydrogen-acceptor triplets
/// where the hydrogen is of type `H_HB` and the acceptor is N, O, or F.
fn generate_hbond_potentials(
    intermediate: &IntermediateSystem,
    type_indices: &HashMap<String, usize>,
    params: &ForceFieldParams,
    config: &ForgeConfig,
) -> Result<Vec<HBondPotential>, Error> {
    let mut h_bonds = Vec::new();
    let mut seen_triplets = HashSet::new();

    let d_hb_raw = match &config.charge_method {
        ChargeMethod::None => params.hydrogen_bond.d0_no_charge,
        ChargeMethod::Qeq(_) | ChargeMethod::Hybrid(_) => params.hydrogen_bond.d0_explicit,
    };
    let r_hb_raw = params.hydrogen_bond.r0;

    let (d_hb, r_hb_sq) = HydrogenBond::<4>::precompute(d_hb_raw, r_hb_raw);

    let hydrogen_type_idx = match type_indices.get(H_BOND_DONOR_HYDROGEN_TYPE) {
        Some(&idx) => idx,
        None => return Ok(h_bonds),
    };

    for h_atom in &intermediate.atoms {
        if h_atom.atom_type != H_BOND_DONOR_HYDROGEN_TYPE {
            continue;
        }

        if h_atom.neighbors.is_empty() {
            continue;
        }
        let donor_atom = &intermediate.atoms[h_atom.neighbors[0]];
        let donor_type_idx = *type_indices.get(&donor_atom.atom_type).ok_or_else(|| {
            Error::Conversion(format!("missing type index for '{}'", donor_atom.atom_type))
        })?;

        for acc_atom in &intermediate.atoms {
            if is_hbond_acceptor(acc_atom.element) {
                let acceptor_type_idx =
                    *type_indices.get(&acc_atom.atom_type).ok_or_else(|| {
                        Error::Conversion(format!(
                            "missing type index for '{}'",
                            acc_atom.atom_type
                        ))
                    })?;

                let triplet = (donor_type_idx, hydrogen_type_idx, acceptor_type_idx);
                if seen_triplets.insert(triplet) {
                    h_bonds.push(HBondPotential {
                        donor_type_idx,
                        hydrogen_type_idx,
                        acceptor_type_idx,
                        d_hb,
                        r_hb_sq,
                    });
                }
            }
        }
    }

    Ok(h_bonds)
}

/// Checks if an element can act as a hydrogen bond acceptor (N, O, F).
fn is_hbond_acceptor(element: Element) -> bool {
    matches!(element, Element::N | Element::O | Element::F)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forge::intermediate::{
        IntermediateAngle, IntermediateDihedral, IntermediateImproper, IntermediateSystem,
        PhysicalBondOrder,
    };
    use crate::model::atom::Atom;
    use crate::model::system::{Bond, System as ModelSystem};
    use crate::model::types::{BondOrder, Element};

    fn make_water() -> ModelSystem {
        let mut sys = ModelSystem::new();
        sys.atoms.push(Atom::new(Element::O, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [0.96, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.24, 0.93, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys
    }

    fn make_typed_water() -> IntermediateSystem {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();
        int.atoms[0].atom_type = "O_3".to_string();
        int.atoms[1].atom_type = H_BOND_DONOR_HYDROGEN_TYPE.to_string();
        int.atoms[2].atom_type = H_BOND_DONOR_HYDROGEN_TYPE.to_string();
        for bond in &mut int.bonds {
            bond.physical_order = Some(PhysicalBondOrder::Single);
        }
        int.angles.push(IntermediateAngle { i: 1, j: 0, k: 2 });
        int
    }

    fn make_typed_ethane() -> IntermediateSystem {
        use super::super::intermediate::Hybridization;

        let mut sys = ModelSystem::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::C, [1.54, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [1.90, 1.03, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 3, BondOrder::Single));

        let mut int = IntermediateSystem::from_system(&sys).unwrap();
        int.atoms[0].atom_type = "C_3".to_string();
        int.atoms[0].hybridization = Hybridization::SP3;
        int.atoms[1].atom_type = "C_3".to_string();
        int.atoms[1].hybridization = Hybridization::SP3;
        int.atoms[2].atom_type = "H_".to_string();
        int.atoms[2].hybridization = Hybridization::None;
        int.atoms[3].atom_type = "H_".to_string();
        int.atoms[3].hybridization = Hybridization::None;
        for bond in &mut int.bonds {
            bond.physical_order = Some(PhysicalBondOrder::Single);
        }
        int.angles.push(IntermediateAngle { i: 2, j: 0, k: 1 });
        int.angles.push(IntermediateAngle { i: 0, j: 1, k: 3 });
        int.dihedrals.push(IntermediateDihedral {
            i: 2,
            j: 0,
            k: 1,
            l: 3,
        });
        int
    }

    fn make_typed_formaldehyde() -> IntermediateSystem {
        let mut sys = ModelSystem::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::O, [1.2, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.5, 0.9, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.5, -0.9, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Double));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 3, BondOrder::Single));

        let mut int = IntermediateSystem::from_system(&sys).unwrap();
        int.atoms[0].atom_type = "C_2".to_string();
        int.atoms[1].atom_type = "O_2".to_string();
        int.atoms[2].atom_type = "H_".to_string();
        int.atoms[3].atom_type = "H_".to_string();
        int.bonds[0].physical_order = Some(PhysicalBondOrder::Double);
        int.bonds[1].physical_order = Some(PhysicalBondOrder::Single);
        int.bonds[2].physical_order = Some(PhysicalBondOrder::Single);
        int.impropers.push(IntermediateImproper {
            center: 0,
            p1: 1,
            p2: 2,
            p3: 3,
        });
        int
    }

    #[test]
    fn collects_unique_atom_types_sorted() {
        let int = make_typed_water();
        let (types, indices) = collect_atom_types(&int);

        assert_eq!(types.len(), 2);
        assert!(types.contains(&"O_3".to_string()));
        assert!(types.contains(&H_BOND_DONOR_HYDROGEN_TYPE.to_string()));
        assert!(indices.contains_key("O_3"));
        assert!(indices.contains_key(H_BOND_DONOR_HYDROGEN_TYPE));
        assert!(types[0] < types[1] || types[1] < types[0]);
    }

    #[test]
    fn generates_atom_properties_with_correct_mass() {
        let int = make_typed_water();
        let (_, indices) = collect_atom_types(&int);
        let props = generate_atom_properties(&int, &indices).unwrap();

        assert_eq!(props.len(), 3);
        assert!((props[0].mass - 15.999).abs() < 0.01);
        assert!((props[1].mass - 1.008).abs() < 0.01);
        assert!((props[2].mass - 1.008).abs() < 0.01);
    }

    #[test]
    fn generates_harmonic_bond_potentials() {
        let int = make_typed_water();
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig::default();
        let bonds = generate_bond_potentials(&int, params, &config).unwrap();

        assert_eq!(bonds.len(), 2);
        for bond in &bonds {
            match bond {
                BondPotential::Harmonic { k_force, r0, .. } => {
                    assert!(*k_force > 0.0);
                    assert!(*r0 > 0.0);
                }
                _ => panic!("Expected Harmonic potential"),
            }
        }
    }

    #[test]
    fn generates_morse_bond_potentials() {
        let int = make_typed_water();
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig {
            bond_potential: BondPotentialType::Morse,
            ..Default::default()
        };
        let bonds = generate_bond_potentials(&int, params, &config).unwrap();

        assert_eq!(bonds.len(), 2);
        for bond in &bonds {
            match bond {
                BondPotential::Morse { r0, d0, alpha, .. } => {
                    assert!(*r0 > 0.0);
                    assert!(*d0 > 0.0);
                    assert!(*alpha > 0.0);
                }
                _ => panic!("Expected Morse potential"),
            }
        }
    }

    #[test]
    fn generates_theta_harmonic_angle_potentials() {
        let int = make_typed_water();
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig::default();
        let angles = generate_angle_potentials(&int, params, &config).unwrap();

        assert_eq!(angles.len(), 1);
        match &angles[0] {
            AnglePotential::ThetaHarmonic {
                theta0, k_force, ..
            } => {
                let expected = 104.51_f64.to_radians();
                assert!((theta0 - expected).abs() < 0.01);
                assert!(*k_force > 0.0);
            }
            _ => panic!("Expected ThetaHarmonic potential"),
        }
    }

    #[test]
    fn generates_cosine_harmonic_angle_potentials() {
        let int = make_typed_water();
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig {
            angle_potential: AnglePotentialType::CosineHarmonic,
            ..Default::default()
        };
        let angles = generate_angle_potentials(&int, params, &config).unwrap();

        assert_eq!(angles.len(), 1);
        match &angles[0] {
            AnglePotential::CosineHarmonic {
                theta0, k_force, ..
            } => {
                assert!(*theta0 > 0.0);
                assert!(*k_force > 0.0);
            }
            _ => panic!("Expected CosineHarmonic potential"),
        }
    }

    #[test]
    fn generates_dihedral_potentials() {
        let int = make_typed_ethane();
        let params = super::super::params::get_default_parameters();
        let dihedrals = generate_dihedral_potentials(&int, params).unwrap();

        assert_eq!(dihedrals.len(), 1);
        assert!(dihedrals[0].v_barrier > 0.0);
        assert!(dihedrals[0].periodicity > 0);
    }

    #[test]
    fn generates_improper_potentials_for_planar_center() {
        let int = make_typed_formaldehyde();
        let params = super::super::params::get_default_parameters();
        let impropers = generate_improper_potentials(&int, params).unwrap();

        assert_eq!(impropers.len(), 1);
        match &impropers[0] {
            ImproperPotential::Planar { k_force, chi0, .. } => {
                assert_eq!(*k_force, params.global.inversion_k);
                assert_eq!(*chi0, 0.0);
            }
            _ => panic!("Expected Planar improper potential"),
        }
    }

    #[test]
    fn generates_lj_vdw_potentials() {
        let int = make_typed_water();
        let (atom_types, _) = collect_atom_types(&int);
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig::default();
        let vdw = generate_vdw_potentials(&atom_types, params, &config).unwrap();

        assert_eq!(vdw.len(), 3);
        for pair in &vdw {
            match pair {
                VdwPairPotential::LennardJones { sigma, epsilon, .. } => {
                    assert!(*sigma > 0.0);
                    assert!(*epsilon > 0.0);
                }
                _ => panic!("Expected LennardJones potential"),
            }
        }
    }

    #[test]
    fn generates_exp6_vdw_potentials() {
        let int = make_typed_water();
        let (atom_types, _) = collect_atom_types(&int);
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig {
            vdw_potential: VdwPotentialType::Exponential6,
            ..Default::default()
        };
        let vdw = generate_vdw_potentials(&atom_types, params, &config).unwrap();

        assert_eq!(vdw.len(), 3);
        for pair in &vdw {
            match pair {
                VdwPairPotential::Exponential6 { a, b, c, .. } => {
                    assert!(*a > 0.0);
                    assert!(*b > 0.0);
                    assert!(*c > 0.0);
                }
                _ => panic!("Expected Exponential6 potential"),
            }
        }
    }

    #[test]
    fn generates_hbond_potentials_for_h_hb() {
        let int = make_typed_water();
        let (_, type_indices) = collect_atom_types(&int);
        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig::default();
        let hbonds = generate_hbond_potentials(&int, &type_indices, params, &config).unwrap();

        assert!(!hbonds.is_empty());
        assert_eq!(hbonds.len(), 1);
        for hb in &hbonds {
            assert!(hb.d0 > 0.0);
            assert!(hb.r0 > 0.0);
        }
    }

    #[test]
    fn hbond_d0_depends_on_charge_method() {
        let int = make_typed_water();
        let (_, type_indices) = collect_atom_types(&int);
        let params = super::super::params::get_default_parameters();

        let config_no_charge = ForgeConfig::default();
        let hbonds_no =
            generate_hbond_potentials(&int, &type_indices, params, &config_no_charge).unwrap();

        let config_qeq = ForgeConfig {
            charge_method: ChargeMethod::Qeq(Default::default()),
            ..Default::default()
        };
        let hbonds_qeq =
            generate_hbond_potentials(&int, &type_indices, params, &config_qeq).unwrap();

        let config_hybrid = ForgeConfig {
            charge_method: ChargeMethod::Hybrid(Default::default()),
            ..Default::default()
        };
        let hbonds_hybrid =
            generate_hbond_potentials(&int, &type_indices, params, &config_hybrid).unwrap();

        assert_ne!(hbonds_no[0].d0, hbonds_qeq[0].d0);
        assert_eq!(hbonds_no[0].d0, params.hydrogen_bond.d0_no_charge);
        assert_eq!(hbonds_qeq[0].d0, params.hydrogen_bond.d0_explicit);
        assert_eq!(hbonds_hybrid[0].d0, params.hydrogen_bond.d0_explicit);
    }

    #[test]
    fn hbond_acceptor_detection_positive() {
        assert!(is_hbond_acceptor(Element::O));
        assert!(is_hbond_acceptor(Element::N));
        assert!(is_hbond_acceptor(Element::F));
    }

    #[test]
    fn hbond_acceptor_detection_negative() {
        assert!(!is_hbond_acceptor(Element::C));
        assert!(!is_hbond_acceptor(Element::H));
        assert!(!is_hbond_acceptor(Element::S));
        assert!(!is_hbond_acceptor(Element::Cl));
    }

    #[test]
    fn errors_on_missing_atom_type_parameter() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();
        int.atoms[0].atom_type = "Xx_UNKNOWN".to_string();
        int.atoms[1].atom_type = "H_".to_string();
        int.atoms[2].atom_type = "H_".to_string();

        let params = super::super::params::get_default_parameters();
        let config = ForgeConfig::default();
        let result = generate_bond_potentials(&int, params, &config);

        assert!(matches!(result, Err(Error::MissingParameter { .. })));
    }
}
