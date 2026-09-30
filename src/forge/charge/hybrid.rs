//! Hybrid biological/QEq charge assignment.
//!
//! This module implements the hybrid charge assignment strategy for
//! biological systems, combining classical force field charges for proteins,
//! nucleic acids, water, and ions with QEq charge equilibration for ligands.

use super::spatial::SpatialGrid;
use crate::forge::config::{
    EmbeddedQeqConfig, HybridConfig, LigandChargeConfig, LigandChargeMethod, QeqConfig,
};
use crate::forge::error::Error;
use crate::forge::intermediate::{IntermediateAtom, IntermediateSystem};
use crate::model::metadata::{AtomResidueInfo, BioMetadata, ResidueCategory, StandardResidue};
use crate::model::types::Element;
use cheq::{ExternalPotential, PointCharge, QEqSolver, get_default_parameters};
use ffcharge::{IonScheme, Position as FfPosition};
use std::collections::{HashMap, HashSet};

/// pH threshold for N-terminal deprotonation (NH₃⁺ → NH₂).
const N_TERMINAL_PKA: f64 = 8.0;
/// pH threshold for C-terminal protonation (COO⁻ → COOH).
const C_TERMINAL_PKA: f64 = 3.1;

/// Assigns charges using the hybrid biological/ligand method.
///
/// # Arguments
///
/// * `system` — Mutable reference to the intermediate system
/// * `config` — Hybrid charge configuration
/// * `neutral_termini` — When `true` (MPSim mode), protein N/C termini use the
///   neutral charge sets (`NH₂` / `COOH`) regardless of pH
///
/// # Errors
///
/// Returns [`Error::MissingBioMetadata`] if biological metadata is absent.
/// Returns [`Error::HybridChargeAssignment`] if classical charge lookup fails.
/// Returns [`Error::ChargeCalculation`] if QEq solver fails to converge.
pub fn assign_hybrid_charges(
    system: &mut IntermediateSystem,
    config: &HybridConfig,
    neutral_termini: bool,
    preassigned: &HashMap<usize, f64>,
) -> Result<(), Error> {
    if !system.has_bio_metadata() {
        return Err(Error::MissingBioMetadata);
    }
    let ph = system.effective_ph();

    let metadata = system.bio_metadata.as_ref().unwrap().clone();

    let classification = classify_atoms(&metadata);
    assign_fixed_charges(
        system,
        &metadata,
        config,
        ph,
        &classification,
        neutral_termini,
        preassigned,
    )?;

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
    neutral_termini: bool,
    preassigned: &HashMap<usize, f64>,
) -> Result<(), Error> {
    for (idx, (&class, info)) in classification.iter().zip(&metadata.atom_info).enumerate() {
        // Settled outside the residue libraries, which either have no entry for
        // this atom or hold the modern parameterization rather than the
        // convention being reproduced.
        if let Some(&charge) = preassigned.get(&idx) {
            system.atoms[idx].charge = charge;
            continue;
        }
        let charge = match class {
            AtomClass::Protein => lookup_protein_charge(config, info, ph, neutral_termini)?,
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
///
/// When `neutral_termini` is `true` (MPSim mode), protein N/C termini map to
/// their neutral charge sets (`NH₂` / `COOH`) regardless of pH. This must be
/// paired with the matching neutral terminal topology (see
/// [`crate::forge::mpsim`]) so every terminal residue stays net-neutral.
fn map_residue_position(info: &AtomResidueInfo, ph: f64, neutral_termini: bool) -> FfPosition {
    use crate::model::metadata::ResiduePosition;

    match info.position {
        ResiduePosition::NTerminal => {
            if neutral_termini || ph >= N_TERMINAL_PKA {
                FfPosition::NTerminalDeprotonated // Neutral NH2
            } else {
                FfPosition::NTerminal // Protonated NH3+
            }
        }
        ResiduePosition::CTerminal => {
            if neutral_termini || ph < C_TERMINAL_PKA {
                FfPosition::CTerminalProtonated // Protonated COOH (neutral)
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
    neutral_termini: bool,
) -> Result<f64, Error> {
    let position = map_residue_position(info, ph, neutral_termini);

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
    // Nucleic acids only use 5'/3' positions, which neutral_termini does not affect.
    let position = map_residue_position(info, 7.0, false);

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

        let method = find_ligand_method(&custom_configs, chain_id, residue_id, insertion_code)
            .unwrap_or(&config.default_ligand_method);

        match method {
            LigandChargeMethod::Vacuum(qeq_config) => {
                assign_vacuum_qeq(system, &group.atom_indices, qeq_config)?;
            }
            LigandChargeMethod::Embedded(embedded_config) => {
                assign_embedded_qeq(
                    system,
                    &group.atom_indices,
                    embedded_config,
                    &positions,
                    &fixed_charge_indices,
                )?;
            }
            LigandChargeMethod::Zero => {
                for &idx in &group.atom_indices {
                    system.atoms[idx].charge = 0.0;
                }
            }
            LigandChargeMethod::Formal => {
                assign_formal_charges(system, &group.atom_indices);
            }
        }
    }

    Ok(())
}

/// Finds the charge method for a specific ligand.
fn find_ligand_method<'a>(
    custom_configs: &HashMap<(String, i32, Option<char>), &'a LigandChargeConfig>,
    chain_id: &str,
    residue_id: i32,
    insertion_code: Option<char>,
) -> Option<&'a LigandChargeMethod> {
    if let Some(lc) = custom_configs.get(&(chain_id.to_string(), residue_id, insertion_code)) {
        return Some(&lc.method);
    }
    None
}

/// Assigns integer formal charges to a ligand from its own bond orders.
///
/// The formal charge of an atom is the bond-order sum it carries minus the
/// bond-order sum it would carry while neutral, so a carboxylate oxygen with a
/// single bond takes `-1` and its doubly bonded partner takes `0`. Elements
/// with no well-defined neutral covalent valence, and atoms whose bond orders
/// do not sum to an integer (an aromatic bond contributes 1.5), have no
/// unambiguous formal charge and are left at zero.
///
/// Valences are read from the structure as given, so every hydrogen must be
/// explicit — an atom missing one is indistinguishable from an anion.
fn assign_formal_charges(system: &mut IntermediateSystem, atom_indices: &[usize]) {
    let member: HashSet<usize> = atom_indices.iter().copied().collect();

    let mut order_sum: HashMap<usize, f64> = atom_indices.iter().map(|&i| (i, 0.0)).collect();
    for bond in &system.bonds {
        // A bond leaving the ligand still fills a valence, so it counts whether
        // or not the partner belongs to this group.
        for near in [bond.i, bond.j] {
            if member.contains(&near) {
                *order_sum.entry(near).or_insert(0.0) += bond.order.value();
            }
        }
    }

    for &idx in atom_indices {
        let sum = order_sum.get(&idx).copied().unwrap_or(0.0);
        system.atoms[idx].charge = match neutral_valence(system.atoms[idx].element) {
            Some(valence) if (sum - sum.round()).abs() < BOND_ORDER_EPSILON => {
                sum.round() - valence
            }
            _ => 0.0,
        };
    }
}

/// Bond-order sums within this tolerance of an integer are treated as integral.
const BOND_ORDER_EPSILON: f64 = 1e-6;

/// Bond-order sum an atom of this element carries when it bears no formal
/// charge, or `None` when the element has no single such value.
///
/// Only the main-group elements that appear in drug-like ligands are listed;
/// transition metals and hypervalent centres are deliberately absent, because
/// a formal charge derived from an octet argument would be wrong for them.
fn neutral_valence(element: Element) -> Option<f64> {
    use Element::*;
    Some(match element {
        H | F | Cl | Br | I => 1.0,
        O | S | Se => 2.0,
        N | P | B => 3.0,
        C | Si => 4.0,
        _ => return None,
    })
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

    /// Builds acetate with every hydrogen explicit: `C0` carries a doubly
    /// bonded `O1`, a singly bonded `O2` and the methyl `C3`, which in turn
    /// carries `H4`, `H5` and `H6`. This is the carboxylate motif the formal
    /// charge rule has to get right.
    fn acetate() -> IntermediateSystem {
        use crate::model::types::BondOrder;
        use Element::{C, H, O};

        let mut system = IntermediateSystem {
            atoms: [C, O, O, C, H, H, H]
                .into_iter()
                .map(|e| IntermediateAtom::new(e, [0.0, 0.0, 0.0]))
                .collect(),
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            inversions: Vec::new(),
            bio_metadata: None,
        };
        for (i, j, order) in [
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (3, 5, BondOrder::Single),
            (3, 6, BondOrder::Single),
        ] {
            system
                .bonds
                .push(crate::forge::intermediate::IntermediateBond {
                    i,
                    j,
                    order,
                    physical_order: None,
                });
        }
        system
    }

    #[test]
    fn formal_charge_puts_minus_one_on_the_singly_bonded_oxygen() {
        let mut system = acetate();
        let all: Vec<usize> = (0..system.atoms.len()).collect();
        assign_formal_charges(&mut system, &all);

        let charges: Vec<f64> = system.atoms.iter().map(|a| a.charge).collect();
        // Every atom fills its neutral valence except the carboxylate oxygen,
        // which carries one bond order where oxygen wants two.
        assert_eq!(charges, vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(charges.iter().sum::<f64>(), -1.0);
    }

    #[test]
    fn formal_charge_leaves_aromatic_bond_orders_alone() {
        use crate::model::types::BondOrder;

        let mut system = acetate();
        // A delocalized carboxylate gives both oxygens a bond-order sum of 1.5,
        // which no integer formal charge describes without a Kekule choice.
        system.bonds[0].order = BondOrder::Aromatic;
        system.bonds[1].order = BondOrder::Aromatic;
        let all: Vec<usize> = (0..system.atoms.len()).collect();
        assign_formal_charges(&mut system, &all);

        assert_eq!(system.atoms[1].charge, 0.0);
        assert_eq!(system.atoms[2].charge, 0.0);
    }

    #[test]
    fn formal_charge_counts_bonds_that_leave_the_group() {
        let mut system = acetate();
        // Scoring only the carboxylate head: the bond to the methyl carbon is
        // outside the selection but still fills a valence on C.
        assign_formal_charges(&mut system, &[0, 1, 2]);

        assert_eq!(system.atoms[0].charge, 0.0);
        assert_eq!(system.atoms[2].charge, -1.0);
    }

    #[test]
    fn formal_charge_reports_a_missing_hydrogen_rather_than_hiding_it() {
        let mut system = acetate();
        // Drop one methyl hydrogen without replacing it. The rule has no way to
        // tell an incomplete structure from a genuine carbanion, so the methyl
        // carbon reads as -1 instead of 0. Charges are only meaningful on a
        // fully protonated structure, which is what the pipeline supplies.
        system.bonds.pop();
        assign_formal_charges(&mut system, &[0, 1, 2, 3]);

        assert_eq!(system.atoms[3].charge, -1.0);
    }

    #[test]
    fn neutral_valence_is_absent_for_elements_without_an_octet_argument() {
        assert_eq!(neutral_valence(Element::C), Some(4.0));
        assert_eq!(neutral_valence(Element::N), Some(3.0));
        assert_eq!(neutral_valence(Element::O), Some(2.0));
        assert_eq!(neutral_valence(Element::H), Some(1.0));
        assert_eq!(neutral_valence(Element::Fe), None);
    }

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
    fn map_protein_position_nterminal_normal_ph() {
        let info = AtomResidueInfo::builder("N", "ALA", 1, "A")
            .position(ResiduePosition::NTerminal)
            .build();
        let pos = map_residue_position(&info, 7.0, false);
        assert_eq!(pos, FfPosition::NTerminal);
    }

    #[test]
    fn map_protein_position_nterminal_high_ph() {
        let info = AtomResidueInfo::builder("N", "ALA", 1, "A")
            .position(ResiduePosition::NTerminal)
            .build();
        let pos = map_residue_position(&info, 9.0, false);
        assert_eq!(pos, FfPosition::NTerminalDeprotonated);
    }

    #[test]
    fn map_protein_position_cterminal_normal_ph() {
        let info = AtomResidueInfo::builder("C", "ALA", 1, "A")
            .position(ResiduePosition::CTerminal)
            .build();
        let pos = map_residue_position(&info, 7.0, false);
        assert_eq!(pos, FfPosition::CTerminal);
    }

    #[test]
    fn map_protein_position_cterminal_low_ph() {
        let info = AtomResidueInfo::builder("C", "ALA", 1, "A")
            .position(ResiduePosition::CTerminal)
            .build();
        let pos = map_residue_position(&info, 2.0, false);
        assert_eq!(pos, FfPosition::CTerminalProtonated);
    }

    #[test]
    fn map_protein_termini_neutral_regardless_of_ph_in_mpsim_mode() {
        let n_info = AtomResidueInfo::builder("N", "ALA", 1, "A")
            .position(ResiduePosition::NTerminal)
            .build();
        let c_info = AtomResidueInfo::builder("C", "ALA", 1, "A")
            .position(ResiduePosition::CTerminal)
            .build();

        // At physiological pH, neutral_termini forces the neutral charge sets.
        assert_eq!(
            map_residue_position(&n_info, 7.0, true),
            FfPosition::NTerminalDeprotonated
        );
        assert_eq!(
            map_residue_position(&c_info, 7.0, true),
            FfPosition::CTerminalProtonated
        );
        // Even at extreme pH values the neutral state is kept.
        assert_eq!(
            map_residue_position(&n_info, 1.0, true),
            FfPosition::NTerminalDeprotonated
        );
        assert_eq!(
            map_residue_position(&c_info, 13.0, true),
            FfPosition::CTerminalProtonated
        );
    }

    #[test]
    fn map_nucleic_position_five_prime() {
        let info = AtomResidueInfo::builder("P", "DA", 1, "B")
            .position(ResiduePosition::FivePrime)
            .build();
        let pos = map_residue_position(&info, 7.0, false);
        assert_eq!(pos, FfPosition::FivePrime);
    }

    #[test]
    fn map_nucleic_position_three_prime() {
        let info = AtomResidueInfo::builder("O3'", "DA", 1, "B")
            .position(ResiduePosition::ThreePrime)
            .build();
        let pos = map_residue_position(&info, 7.0, false);
        assert_eq!(pos, FfPosition::ThreePrime);
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
