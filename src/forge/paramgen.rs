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
