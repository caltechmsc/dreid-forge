//! Hybrid biological/QEq charge assignment.
//!
//! This module implements the hybrid charge assignment strategy for
//! biological systems, combining classical force field charges for proteins,
//! nucleic acids, water, and ions with QEq charge equilibration for ligands.

use super::spatial::SpatialGrid;
use crate::forge::config::{
    EmbeddedQeqConfig, HybridConfig, LigandChargeConfig, LigandQeqMethod, QeqConfig,
};
use crate::forge::error::Error;
use crate::forge::intermediate::{IntermediateAtom, IntermediateSystem};
use crate::model::metadata::{AtomResidueInfo, BioMetadata, ResidueCategory, StandardResidue};
use cheq::{ExternalPotential, PointCharge, QEqSolver, get_default_parameters};
use ffcharge::{IonScheme, Position as FfPosition};
use std::collections::HashMap;

/// pH threshold for N-terminal deprotonation (NH₃⁺ → NH₂).
const N_TERMINAL_PKA: f64 = 8.0;
/// pH threshold for C-terminal protonation (COO⁻ → COOH).
const C_TERMINAL_PKA: f64 = 3.1;

/// Assigns charges using the hybrid biological/QEq method.
///
/// # Arguments
///
/// * `system` — Mutable reference to the intermediate system
/// * `config` — Hybrid charge configuration
///
/// # Errors
///
/// Returns [`Error::MissingBioMetadata`] if biological metadata is absent.
/// Returns [`Error::HybridChargeAssignment`] if classical charge lookup fails.
/// Returns [`Error::ChargeCalculation`] if QEq solver fails to converge.
pub fn assign_hybrid_charges(
    system: &mut IntermediateSystem,
    config: &HybridConfig,
) -> Result<(), Error> {
    if !system.has_bio_metadata() {
        return Err(Error::MissingBioMetadata);
    }
    let ph = system.effective_ph();

    let metadata = system.bio_metadata.as_ref().unwrap().clone();

    let classification = classify_atoms(&metadata);
    assign_fixed_charges(system, &metadata, config, ph, &classification)?;

    let ligand_groups = identify_ligand_groups(&metadata, &classification);
    if !ligand_groups.is_empty() {
        assign_ligand_charges(system, config, &classification, &ligand_groups)?;
    }

    Ok(())
}

/// Atom classification for charge assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AtomClass {
    /// Standard amino acid residue.
    Protein,
    /// Standard nucleic acid residue.
    NucleicAcid,
    /// Water molecule (HOH).
    Water,
    /// Monoatomic or polyatomic ion.
    Ion,
    /// Ligand or hetero group (for QEq).
    Ligand,
}

/// Classifies each atom by molecule type.
fn classify_atoms(metadata: &BioMetadata) -> Vec<AtomClass> {
    metadata.atom_info.iter().map(classify_atom).collect()
}

/// Classifies a single atom based on its residue info.
fn classify_atom(info: &AtomResidueInfo) -> AtomClass {
    match info.category {
        ResidueCategory::Ion => AtomClass::Ion,
        ResidueCategory::Standard => {
            if let Some(std_res) = info.standard_name {
                match std_res {
                    StandardResidue::HOH => AtomClass::Water,
                    // Amino acids
                    StandardResidue::ALA
                    | StandardResidue::ARG
                    | StandardResidue::ASN
                    | StandardResidue::ASP
                    | StandardResidue::CYS
                    | StandardResidue::GLN
                    | StandardResidue::GLU
                    | StandardResidue::GLY
                    | StandardResidue::HIS
                    | StandardResidue::ILE
                    | StandardResidue::LEU
                    | StandardResidue::LYS
                    | StandardResidue::MET
                    | StandardResidue::PHE
                    | StandardResidue::PRO
                    | StandardResidue::SER
                    | StandardResidue::THR
                    | StandardResidue::TRP
                    | StandardResidue::TYR
                    | StandardResidue::VAL => AtomClass::Protein,
                    // Nucleotides
                    StandardResidue::A
                    | StandardResidue::C
                    | StandardResidue::G
                    | StandardResidue::U
                    | StandardResidue::I
                    | StandardResidue::DA
                    | StandardResidue::DC
                    | StandardResidue::DG
                    | StandardResidue::DT
                    | StandardResidue::DI => AtomClass::NucleicAcid,
                }
            } else {
                AtomClass::Ligand
            }
        }
        ResidueCategory::Hetero => AtomClass::Ligand,
    }
}

/// Assigns fixed charges to proteins, nucleic acids, water, and ions.
fn assign_fixed_charges(
    system: &mut IntermediateSystem,
    metadata: &BioMetadata,
    config: &HybridConfig,
    ph: f64,
    classification: &[AtomClass],
) -> Result<(), Error> {
    for (idx, (&class, info)) in classification.iter().zip(&metadata.atom_info).enumerate() {
        let charge = match class {
            AtomClass::Protein => lookup_protein_charge(config, info, ph)?,
            AtomClass::NucleicAcid => lookup_nucleic_charge(config, info)?,
            AtomClass::Water => lookup_water_charge(config, info)?,
            AtomClass::Ion => lookup_ion_charge(info)?,
            AtomClass::Ligand => continue,
        };
        system.atoms[idx].charge = charge;
    }
    Ok(())
}

/// Maps our ResiduePosition to ffcharge::Position for proteins.
fn map_residue_position(info: &AtomResidueInfo, ph: f64) -> FfPosition {
    use crate::model::metadata::ResiduePosition;

    match info.position {
        ResiduePosition::NTerminal => {
            if ph < N_TERMINAL_PKA {
                FfPosition::NTerminal // Protonated NH3+
            } else {
                FfPosition::NTerminalDeprotonated // Neutral NH2
            }
        }
        ResiduePosition::CTerminal => {
            if ph < C_TERMINAL_PKA {
                FfPosition::CTerminalProtonated // Protonated COOH
            } else {
                FfPosition::CTerminal // Deprotonated COO-
            }
        }
        ResiduePosition::FivePrime => FfPosition::FivePrime,
        ResiduePosition::ThreePrime => FfPosition::ThreePrime,
        ResiduePosition::Internal | ResiduePosition::None => FfPosition::Middle,
    }
}

/// Looks up protein charge from ffcharge.
fn lookup_protein_charge(
    config: &HybridConfig,
    info: &AtomResidueInfo,
    ph: f64,
) -> Result<f64, Error> {
    let position = map_residue_position(info, ph);

    config
        .protein_scheme
        .charge(position, &info.residue_name, &info.atom_name)
        .map(|c| c as f64)
        .ok_or_else(|| {
            Error::hybrid_charge_assignment(
                info.chain_id.clone(),
                info.residue_id,
                &info.residue_name,
                format!(
                    "protein charge not found for atom '{}' at position {:?}",
                    info.atom_name, position
                ),
            )
        })
}

/// Looks up nucleic acid charge from ffcharge.
fn lookup_nucleic_charge(config: &HybridConfig, info: &AtomResidueInfo) -> Result<f64, Error> {
    let position = map_residue_position(info, 7.0); // Assuming pH 7.0 for nucleic acids (no effect)

    config
        .nucleic_scheme
        .charge(position, &info.residue_name, &info.atom_name)
        .map(|c| c as f64)
        .ok_or_else(|| {
            Error::hybrid_charge_assignment(
                info.chain_id.clone(),
                info.residue_id,
                &info.residue_name,
                format!(
                    "nucleic acid charge not found for atom '{}' at position {:?}",
                    info.atom_name, position
                ),
            )
        })
}

/// Looks up water charge from ffcharge.
fn lookup_water_charge(config: &HybridConfig, info: &AtomResidueInfo) -> Result<f64, Error> {
    config
        .water_scheme
        .charges()
        .ok_or_else(|| {
            Error::hybrid_charge_assignment(
                info.chain_id.clone(),
                info.residue_id,
                &info.residue_name,
                "water charge parameters not found".to_string(),
            )
        })
        .and_then(|charges| {
            let charge = if info.atom_name == "O" {
                charges.o
            } else if info.atom_name == "H1" {
                charges.h1
            } else if info.atom_name == "H2" {
                charges.h2
            } else {
                return Err(Error::hybrid_charge_assignment(
                    info.chain_id.clone(),
                    info.residue_id,
                    &info.residue_name,
                    format!("unknown water atom name: '{}'", info.atom_name),
                ));
            };
            Ok(charge as f64)
        })
}

/// Looks up ion charge from ffcharge.
fn lookup_ion_charge(info: &AtomResidueInfo) -> Result<f64, Error> {
    IonScheme::Classic
        .charge(&info.residue_name)
        .map(|c| c as f64)
        .ok_or_else(|| {
            Error::hybrid_charge_assignment(
                info.chain_id.clone(),
                info.residue_id,
                &info.residue_name,
                format!("ion charge not found for residue '{}'", info.residue_name),
            )
        })
}

/// Represents a group of atoms belonging to the same ligand residue.
#[derive(Debug)]
struct LigandGroup {
    /// Unique residue key (chain_id, residue_id, insertion_code).
    key: (String, i32, Option<char>),
    /// Atom indices in the system.
    atom_indices: Vec<usize>,
}

/// Identifies ligand groups from metadata.
fn identify_ligand_groups(
    metadata: &BioMetadata,
    classification: &[AtomClass],
) -> Vec<LigandGroup> {
    let mut groups: HashMap<(String, i32, Option<char>), Vec<usize>> = HashMap::new();

    for (idx, (&class, info)) in classification.iter().zip(&metadata.atom_info).enumerate() {
        if class == AtomClass::Ligand {
            let key = (info.chain_id.clone(), info.residue_id, info.insertion_code);
            groups.entry(key).or_default().push(idx);
        }
    }

    groups
        .into_iter()
        .map(|(key, atom_indices)| LigandGroup { key, atom_indices })
        .collect()
}

/// Assigns charges to ligand groups using QEq.
fn assign_ligand_charges(
    system: &mut IntermediateSystem,
    config: &HybridConfig,
    classification: &[AtomClass],
    ligand_groups: &[LigandGroup],
) -> Result<(), Error> {
    let positions: Vec<[f64; 3]> = system.atoms.iter().map(|a| a.position).collect();
    let fixed_charge_indices: Vec<usize> = classification
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != AtomClass::Ligand)
        .map(|(i, _)| i)
        .collect();

    let custom_configs: HashMap<(String, i32, Option<char>), &LigandChargeConfig> = config
        .ligand_configs
        .iter()
        .map(|lc| {
            let key = (
                lc.selector.chain_id.clone(),
                lc.selector.residue_id,
                lc.selector.insertion_code,
            );
            (key, lc)
        })
        .collect();

    for group in ligand_groups {
        let (ref chain_id, residue_id, insertion_code) = group.key;

        let default_method = LigandQeqMethod::Vacuum(config.default_ligand_qeq.clone());
        let method = find_ligand_method(&custom_configs, chain_id, residue_id, insertion_code)
            .unwrap_or(&default_method);

        match method {
            LigandQeqMethod::Vacuum(qeq_config) => {
                assign_vacuum_qeq(system, &group.atom_indices, qeq_config)?;
            }
            LigandQeqMethod::Embedded(embedded_config) => {
                assign_embedded_qeq(
                    system,
                    &group.atom_indices,
                    embedded_config,
                    &positions,
                    &fixed_charge_indices,
                )?;
            }
        }
    }

    Ok(())
}

/// Finds the QEq method for a specific ligand.
fn find_ligand_method<'a>(
    custom_configs: &HashMap<(String, i32, Option<char>), &'a LigandChargeConfig>,
    chain_id: &str,
    residue_id: i32,
    insertion_code: Option<char>,
) -> Option<&'a LigandQeqMethod> {
    if let Some(lc) = custom_configs.get(&(chain_id.to_string(), residue_id, insertion_code)) {
        return Some(&lc.method);
    }
    None
}

/// Assigns charges to a ligand using vacuum QEq.
fn assign_vacuum_qeq(
    system: &mut IntermediateSystem,
    atom_indices: &[usize],
    config: &QeqConfig,
) -> Result<(), Error> {
    let params = get_default_parameters();
    let solver = QEqSolver::new(params).with_options(config.solver_options);

    let ligand_atoms: Vec<IntermediateAtom> = atom_indices
        .iter()
        .map(|&i| system.atoms[i].clone())
        .collect();

    let result = solver.solve(&ligand_atoms, config.total_charge)?;

    for (local_idx, &global_idx) in atom_indices.iter().enumerate() {
        system.atoms[global_idx].charge = result.charges[local_idx];
    }

    Ok(())
}

/// Assigns charges to a ligand using embedded QEq.
fn assign_embedded_qeq(
    system: &mut IntermediateSystem,
    atom_indices: &[usize],
    config: &EmbeddedQeqConfig,
    all_positions: &[[f64; 3]],
    fixed_charge_indices: &[usize],
) -> Result<(), Error> {
    let params = get_default_parameters();
    let solver = QEqSolver::new(params).with_options(config.qeq.solver_options);

    let ligand_atoms: Vec<IntermediateAtom> = atom_indices
        .iter()
        .map(|&i| system.atoms[i].clone())
        .collect();
    let ligand_positions: Vec<[f64; 3]> = atom_indices.iter().map(|&i| all_positions[i]).collect();

    let fixed_positions: Vec<[f64; 3]> = fixed_charge_indices
        .iter()
        .map(|&i| all_positions[i])
        .collect();
    let grid = SpatialGrid::from_positions(&fixed_positions, config.cutoff_radius);

    let env_local_indices =
        grid.query_radius_multi(&ligand_positions, &fixed_positions, config.cutoff_radius);

    let point_charges: Vec<PointCharge> = env_local_indices
        .iter()
        .map(|&local_idx| {
            let global_idx = fixed_charge_indices[local_idx];
            let atom = &system.atoms[global_idx];
            PointCharge::new(atom.element.atomic_number(), atom.position, atom.charge)
        })
        .collect();

    let external = ExternalPotential::from_point_charges(point_charges);

    let result = solver.solve_in_field(&ligand_atoms, config.qeq.total_charge, &external)?;

    for (local_idx, &global_idx) in atom_indices.iter().enumerate() {
        system.atoms[global_idx].charge = result.charges[local_idx];
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::metadata::{AtomResidueInfo, ResiduePosition};

    #[test]
    fn classify_protein_atom() {
        let info = AtomResidueInfo::builder("CA", "ALA", 1, "A")
            .standard_name(Some(StandardResidue::ALA))
            .category(ResidueCategory::Standard)
            .build();
        assert_eq!(classify_atom(&info), AtomClass::Protein);
    }

    #[test]
    fn classify_nucleic_acid_atom() {
        let info = AtomResidueInfo::builder("C1'", "DA", 1, "B")
            .standard_name(Some(StandardResidue::DA))
            .category(ResidueCategory::Standard)
            .build();
        assert_eq!(classify_atom(&info), AtomClass::NucleicAcid);
    }

    #[test]
    fn classify_water_atom() {
        let info = AtomResidueInfo::builder("O", "HOH", 1, "W")
            .standard_name(Some(StandardResidue::HOH))
            .category(ResidueCategory::Standard)
            .build();
        assert_eq!(classify_atom(&info), AtomClass::Water);
    }

    #[test]
    fn classify_ion_atom() {
        let info = AtomResidueInfo::builder("NA", "NA", 1, "I")
            .category(ResidueCategory::Ion)
            .build();
        assert_eq!(classify_atom(&info), AtomClass::Ion);
    }

    #[test]
    fn classify_ligand_atom() {
        let info = AtomResidueInfo::builder("C1", "LIG", 1, "L")
            .category(ResidueCategory::Hetero)
            .build();
        assert_eq!(classify_atom(&info), AtomClass::Ligand);
    }

    #[test]
    fn map_protein_position_nterminal_low_ph() {
        let info = AtomResidueInfo::builder("N", "ALA", 1, "A")
            .position(ResiduePosition::NTerminal)
            .build();
        let pos = map_residue_position(&info, 7.0);
        assert_eq!(pos, FfPosition::NTerminal);
    }

    #[test]
    fn lookup_water_charge_oxygen() {
        let config = HybridConfig::default();
        let charge = lookup_water_charge(
            &config,
            &AtomResidueInfo::builder("O", "HOH", 1, "W").build(),
        )
        .unwrap();
        assert!((charge - (-0.834)).abs() < 1e-6);
    }

    #[test]
    fn lookup_water_charge_hydrogen() {
        let config = HybridConfig::default();
        let charge = lookup_water_charge(
            &config,
            &AtomResidueInfo::builder("H1", "HOH", 1, "W").build(),
        )
        .unwrap();
        assert!((charge - 0.417).abs() < 1e-6);
    }

    #[test]
    fn lookup_ion_charge_sodium() {
        let charge =
            lookup_ion_charge(&AtomResidueInfo::builder("NA", "NA", 1, "I").build()).unwrap();
        assert!((charge - 1.0).abs() < 1e-6);
    }

    #[test]
    fn lookup_ion_charge_chloride() {
        let charge =
            lookup_ion_charge(&AtomResidueInfo::builder("CL", "CL", 1, "I").build()).unwrap();
        assert!((charge - (-1.0)).abs() < 1e-6);
    }
}
