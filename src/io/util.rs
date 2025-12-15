//! Internal utilities for data model conversion.
//!
//! This module provides conversion functions between dreid-forge's
//! internal data model and bio-forge's data structures. These are
//! used internally by readers and writers and are not part of the
//! public API.

use super::{
    Anion, Cation, CleanConfig, HisStrategy, ProtonationConfig, SolvateConfig, TopologyConfig,
};
use crate::model::{
    atom::Atom,
    metadata::{AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition, StandardResidue},
    system::{Bond, System},
    types::{BondOrder, Element},
};
use bio_forge as bf;
use std::collections::HashSet;
use std::str::FromStr;

/// Errors that occur during data model conversion.
///
/// These errors indicate incompatibilities between the dreid-forge
/// and bio-forge data models, such as unsupported elements or
/// missing required metadata.
#[derive(Debug, thiserror::Error)]
pub enum ConversionError {
    /// Encountered bio-forge's `Unknown` element variant.
    #[error("encountered an 'Unknown' element from bio-forge which is not supported")]
    UnsupportedElement,
    /// Encountered a bond order not representable in this model.
    #[error("encountered an unsupported bond order '{0}' during model conversion")]
    UnsupportedBondOrder(String),
    /// System lacks required biological metadata.
    #[error("inconsistent data: system is missing required BioMetadata for this conversion")]
    MissingBioMetadata,
    /// Residue name not recognized as a standard residue.
    #[error("unknown standard residue name")]
    UnknownStandardResidue,
    /// Residue category not recognized.
    #[error("unknown residue category")]
    UnknownResidueCategory,
    /// Residue position not recognized.
    #[error("unknown residue position")]
    UnknownResiduePosition,
}

/// Converts a dreid-forge [`CleanConfig`] to bio-forge's format.
pub fn to_bf_clean_config(config: CleanConfig) -> bf::ops::CleanConfig {
    bf::ops::CleanConfig {
        remove_water: config.remove_water,
        remove_ions: config.remove_ions,
        remove_hydrogens: config.remove_hydrogens,
        remove_hetero: config.remove_hetero,
        remove_residue_names: config.remove_residue_names,
        keep_residue_names: config.keep_residue_names,
    }
}

/// Converts a dreid-forge [`ProtonationConfig`] to bio-forge's format.
pub fn to_bf_hydro_config(config: ProtonationConfig) -> bf::ops::HydroConfig {
    bf::ops::HydroConfig {
        target_ph: config.target_ph,
        remove_existing_h: config.remove_existing_h,
        his_strategy: match config.his_strategy {
            HisStrategy::DirectHID => bf::ops::HisStrategy::DirectHID,
            HisStrategy::DirectHIE => bf::ops::HisStrategy::DirectHIE,
            HisStrategy::Random => bf::ops::HisStrategy::Random,
            HisStrategy::HbNetwork => bf::ops::HisStrategy::HbNetwork,
        },
    }
}

/// Converts a dreid-forge [`SolvateConfig`] to bio-forge's format.
pub fn to_bf_solvate_config(config: SolvateConfig) -> bf::ops::SolvateConfig {
    bf::ops::SolvateConfig {
        margin: config.margin,
        water_spacing: config.water_spacing,
        vdw_cutoff: config.vdw_cutoff,
        remove_existing: config.remove_existing,
        cations: config.cations.into_iter().map(to_bf_cation).collect(),
        anions: config.anions.into_iter().map(to_bf_anion).collect(),
        target_charge: config.target_charge,
        rng_seed: config.rng_seed,
    }
}

/// Converts a dreid-forge [`Cation`] to bio-forge's format.
pub fn to_bf_cation(cation: Cation) -> bf::ops::Cation {
    match cation {
        Cation::Na => bf::ops::Cation::Na,
        Cation::K => bf::ops::Cation::K,
        Cation::Mg => bf::ops::Cation::Mg,
        Cation::Ca => bf::ops::Cation::Ca,
        Cation::Li => bf::ops::Cation::Li,
        Cation::Zn => bf::ops::Cation::Zn,
    }
}

/// Converts a dreid-forge [`Anion`] to bio-forge's format.
pub fn to_bf_anion(anion: Anion) -> bf::ops::Anion {
    match anion {
        Anion::Cl => bf::ops::Anion::Cl,
        Anion::Br => bf::ops::Anion::Br,
        Anion::I => bf::ops::Anion::I,
        Anion::F => bf::ops::Anion::F,
    }
}

/// Builds a bio-forge [`Topology`](bf::Topology) from a
/// [`Structure`](bf::Structure) and a dreid-forge [`TopologyConfig`].
pub fn build_topology(
    structure: bf::Structure,
    config: &TopologyConfig,
) -> Result<bf::Topology, bf::ops::Error> {
    let mut topology_builder =
        bf::ops::TopologyBuilder::new().disulfide_cutoff(config.disulfide_bond_cutoff);

    for tpl in &config.hetero_templates {
        topology_builder = topology_builder.add_hetero_template(tpl.clone());
    }

    topology_builder.build(structure)
}

/// Converts a bio-forge [`Topology`](bf::Topology) to a dreid-forge [`System`].
pub fn from_bio_topology(bio_topo: bf::Topology) -> Result<System, ConversionError> {
    let bio_struct = bio_topo.structure();
    let atom_count = bio_struct.atom_count();

    let mut atoms = Vec::with_capacity(atom_count);
    let mut metadata = BioMetadata::with_capacity(atom_count);

    for (chain, residue, bio_atom) in bio_struct.iter_atoms_with_context() {
        atoms.push(Atom {
            element: convert_element_from_bf(bio_atom.element)?,
            position: [bio_atom.pos.x, bio_atom.pos.y, bio_atom.pos.z],
        });

        let info = AtomResidueInfo::builder(
            bio_atom.name.clone(),
            residue.name.clone(),
            residue.id,
            chain.id.chars().next().unwrap_or(' '),
        )
        .insertion_code_opt(residue.insertion_code)
        .standard_name(convert_std_res_from_bf(residue.standard_name)?)
        .category(convert_res_cat_from_bf(residue.category)?)
        .position(convert_res_pos_from_bf(residue.position)?)
        .build();

        metadata.atom_info.push(info);
    }

    let bonds = bio_topo
        .bonds()
        .iter()
        .map(|bio_bond| {
            Ok(Bond {
                i: bio_bond.a1_idx,
                j: bio_bond.a2_idx,
                order: convert_bond_order_from_bf(bio_bond.order)?,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(System {
        atoms,
        bonds,
        box_vectors: bio_struct.box_vectors,
        bio_metadata: Some(metadata),
    })
}

/// Converts a dreid-forge [`System`] to a bio-forge [`Topology`](bf::Topology).
pub fn to_bio_topology(system: &System) -> Result<bf::Topology, ConversionError> {
    use indexmap::IndexMap;
    use indexmap::map::Entry;

    let metadata = system
        .bio_metadata
        .as_ref()
        .ok_or(ConversionError::MissingBioMetadata)?;

    type ResidueKey = (char, i32, char);
    let mut chains: IndexMap<char, IndexMap<ResidueKey, bf::Residue>> = IndexMap::new();

    for (atom, info) in system.atoms.iter().zip(metadata.atom_info.iter()) {
        let bio_atom = bf::Atom::new(
            &info.atom_name,
            convert_element_to_bf(atom.element)?,
            bf::Point::new(atom.position[0], atom.position[1], atom.position[2]),
        );

        let residue_key = (info.chain_id, info.residue_id, info.insertion_code);
        let residues = chains.entry(info.chain_id).or_default();

        let residue = match residues.entry(residue_key) {
            Entry::Occupied(e) => e.into_mut(),
            Entry::Vacant(e) => {
                let mut res = bf::Residue::new(
                    info.residue_id,
                    Some(info.insertion_code).filter(|&c| c != ' '),
                    &info.residue_name,
                    convert_std_res_to_bf(info.standard_name)?,
                    convert_res_cat_to_bf(info.category)?,
                );
                res.position = convert_res_pos_to_bf(info.position)?;
                e.insert(res)
            }
        };

        residue.add_atom(bio_atom);
    }

    let mut bio_struct = bf::Structure::new();
    bio_struct.box_vectors = system.box_vectors;

    for (chain_id, residues) in chains {
        let mut chain = bf::Chain::new(&chain_id.to_string());
        for (_, residue) in residues {
            chain.add_residue(residue);
        }
        bio_struct.add_chain(chain);
    }

    let bio_bonds = system
        .bonds
        .iter()
        .map(|bond| {
            Ok(bf::Bond::new(
                bond.i,
                bond.j,
                convert_bond_order_to_bf(bond.order)?,
            ))
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(bf::Topology::new(bio_struct, bio_bonds))
}

/// Converts a bio-forge [`Element`] to a dreid-forge [`Element`].
fn convert_element_from_bf(e: bf::Element) -> Result<Element, ConversionError> {
    if matches!(e, bf::Element::Unknown) {
        return Err(ConversionError::UnsupportedElement);
    }

    Element::from_str(e.symbol()).map_err(|_| ConversionError::UnsupportedElement)
}

/// Converts a dreid-forge [`Element`] to a bio-forge [`Element`].
fn convert_element_to_bf(e: Element) -> Result<bf::Element, ConversionError> {
    bf::Element::from_str(e.symbol()).map_err(|_| ConversionError::UnsupportedElement)
}

/// Converts a bio-forge [`BondOrder`] to a dreid-forge [`BondOrder`].
fn convert_bond_order_from_bf(order: bf::BondOrder) -> Result<BondOrder, ConversionError> {
    match order {
        bf::BondOrder::Single => Ok(BondOrder::Single),
        bf::BondOrder::Double => Ok(BondOrder::Double),
        bf::BondOrder::Triple => Ok(BondOrder::Triple),
        bf::BondOrder::Aromatic => Ok(BondOrder::Aromatic),
    }
}
/// Converts a dreid-forge [`BondOrder`] to a bio-forge [`BondOrder`].
fn convert_bond_order_to_bf(order: BondOrder) -> Result<bf::BondOrder, ConversionError> {
    match order {
        BondOrder::Single => Ok(bf::BondOrder::Single),
        BondOrder::Double => Ok(bf::BondOrder::Double),
        BondOrder::Triple => Ok(bf::BondOrder::Triple),
        BondOrder::Aromatic => Ok(bf::BondOrder::Aromatic),
    }
}

/// Converts an optional bio-forge [`StandardResidue`] to a dreid-forge [`StandardResidue`].
fn convert_std_res_from_bf(
    res: Option<bf::StandardResidue>,
) -> Result<Option<StandardResidue>, ConversionError> {
    Ok(res.map(|v| match v {
        bf::StandardResidue::ALA => StandardResidue::ALA,
        bf::StandardResidue::ARG => StandardResidue::ARG,
        bf::StandardResidue::ASN => StandardResidue::ASN,
        bf::StandardResidue::ASP => StandardResidue::ASP,
        bf::StandardResidue::CYS => StandardResidue::CYS,
        bf::StandardResidue::GLN => StandardResidue::GLN,
        bf::StandardResidue::GLU => StandardResidue::GLU,
        bf::StandardResidue::GLY => StandardResidue::GLY,
        bf::StandardResidue::HIS => StandardResidue::HIS,
        bf::StandardResidue::ILE => StandardResidue::ILE,
        bf::StandardResidue::LEU => StandardResidue::LEU,
        bf::StandardResidue::LYS => StandardResidue::LYS,
        bf::StandardResidue::MET => StandardResidue::MET,
        bf::StandardResidue::PHE => StandardResidue::PHE,
        bf::StandardResidue::PRO => StandardResidue::PRO,
        bf::StandardResidue::SER => StandardResidue::SER,
        bf::StandardResidue::THR => StandardResidue::THR,
        bf::StandardResidue::TRP => StandardResidue::TRP,
        bf::StandardResidue::TYR => StandardResidue::TYR,
        bf::StandardResidue::VAL => StandardResidue::VAL,
        bf::StandardResidue::A => StandardResidue::A,
        bf::StandardResidue::C => StandardResidue::C,
        bf::StandardResidue::G => StandardResidue::G,
        bf::StandardResidue::U => StandardResidue::U,
        bf::StandardResidue::I => StandardResidue::I,
        bf::StandardResidue::DA => StandardResidue::DA,
        bf::StandardResidue::DC => StandardResidue::DC,
        bf::StandardResidue::DG => StandardResidue::DG,
        bf::StandardResidue::DT => StandardResidue::DT,
        bf::StandardResidue::DI => StandardResidue::DI,
        bf::StandardResidue::HOH => StandardResidue::HOH,
    }))
}

/// Converts an optional dreid-forge [`StandardResidue`] to a bio-forge [`StandardResidue`].
fn convert_std_res_to_bf(
    res: Option<StandardResidue>,
) -> Result<Option<bf::StandardResidue>, ConversionError> {
    Ok(res.map(|v| match v {
        StandardResidue::ALA => bf::StandardResidue::ALA,
        StandardResidue::ARG => bf::StandardResidue::ARG,
        StandardResidue::ASN => bf::StandardResidue::ASN,
        StandardResidue::ASP => bf::StandardResidue::ASP,
        StandardResidue::CYS => bf::StandardResidue::CYS,
        StandardResidue::GLN => bf::StandardResidue::GLN,
        StandardResidue::GLU => bf::StandardResidue::GLU,
        StandardResidue::GLY => bf::StandardResidue::GLY,
        StandardResidue::HIS => bf::StandardResidue::HIS,
        StandardResidue::ILE => bf::StandardResidue::ILE,
        StandardResidue::LEU => bf::StandardResidue::LEU,
        StandardResidue::LYS => bf::StandardResidue::LYS,
        StandardResidue::MET => bf::StandardResidue::MET,
        StandardResidue::PHE => bf::StandardResidue::PHE,
        StandardResidue::PRO => bf::StandardResidue::PRO,
        StandardResidue::SER => bf::StandardResidue::SER,
        StandardResidue::THR => bf::StandardResidue::THR,
        StandardResidue::TRP => bf::StandardResidue::TRP,
        StandardResidue::TYR => bf::StandardResidue::TYR,
        StandardResidue::VAL => bf::StandardResidue::VAL,
        StandardResidue::A => bf::StandardResidue::A,
        StandardResidue::C => bf::StandardResidue::C,
        StandardResidue::G => bf::StandardResidue::G,
        StandardResidue::U => bf::StandardResidue::U,
        StandardResidue::I => bf::StandardResidue::I,
        StandardResidue::DA => bf::StandardResidue::DA,
        StandardResidue::DC => bf::StandardResidue::DC,
        StandardResidue::DG => bf::StandardResidue::DG,
        StandardResidue::DT => bf::StandardResidue::DT,
        StandardResidue::DI => bf::StandardResidue::DI,
        StandardResidue::HOH => bf::StandardResidue::HOH,
    }))
}

/// Converts a bio-forge [`ResidueCategory`] to a dreid-forge [`ResidueCategory`].
fn convert_res_cat_from_bf(cat: bf::ResidueCategory) -> Result<ResidueCategory, ConversionError> {
    Ok(match cat {
        bf::ResidueCategory::Standard => ResidueCategory::Standard,
        bf::ResidueCategory::Hetero => ResidueCategory::Hetero,
        bf::ResidueCategory::Ion => ResidueCategory::Ion,
    })
}
/// Converts a dreid-forge [`ResidueCategory`] to a bio-forge [`ResidueCategory`].
fn convert_res_cat_to_bf(cat: ResidueCategory) -> Result<bf::ResidueCategory, ConversionError> {
    Ok(match cat {
        ResidueCategory::Standard => bf::ResidueCategory::Standard,
        ResidueCategory::Hetero => bf::ResidueCategory::Hetero,
        ResidueCategory::Ion => bf::ResidueCategory::Ion,
    })
}

/// Converts a bio-forge [`ResiduePosition`] to a dreid-forge [`ResiduePosition`].
fn convert_res_pos_from_bf(pos: bf::ResiduePosition) -> Result<ResiduePosition, ConversionError> {
    Ok(match pos {
        bf::ResiduePosition::None => ResiduePosition::None,
        bf::ResiduePosition::Internal => ResiduePosition::Internal,
        bf::ResiduePosition::NTerminal => ResiduePosition::NTerminal,
        bf::ResiduePosition::CTerminal => ResiduePosition::CTerminal,
        bf::ResiduePosition::FivePrime => ResiduePosition::FivePrime,
        bf::ResiduePosition::ThreePrime => ResiduePosition::ThreePrime,
    })
}
/// Converts a dreid-forge [`ResiduePosition`] to a bio-forge [`ResiduePosition`].
fn convert_res_pos_to_bf(pos: ResiduePosition) -> Result<bf::ResiduePosition, ConversionError> {
    Ok(match pos {
        ResiduePosition::None => bf::ResiduePosition::None,
        ResiduePosition::Internal => bf::ResiduePosition::Internal,
        ResiduePosition::NTerminal => bf::ResiduePosition::NTerminal,
        ResiduePosition::CTerminal => bf::ResiduePosition::CTerminal,
        ResiduePosition::FivePrime => bf::ResiduePosition::FivePrime,
        ResiduePosition::ThreePrime => bf::ResiduePosition::ThreePrime,
    })
}

/// Attempts to guess the chemical element from a string token.
pub fn guess_element_symbol(token: &str) -> Option<Element> {
    if token.trim().is_empty() {
        return None;
    }

    let mut candidates = Vec::new();
    let mut seen = HashSet::new();
    let push = |s: String, candidates: &mut Vec<String>, seen: &mut HashSet<String>| {
        if seen.insert(s.clone()) {
            candidates.push(s);
        }
    };
    let trimmed = token.trim();
    push(trimmed.to_string(), &mut candidates, &mut seen);

    let alpha_prefix: String = trimmed
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    if !alpha_prefix.is_empty() {
        push(alpha_prefix.clone(), &mut candidates, &mut seen);
        let mut chars = alpha_prefix.chars();
        if let Some(first) = chars.next() {
            let mut canon = String::new();
            canon.push(first.to_ascii_uppercase());
            for c in chars {
                canon.push(c.to_ascii_lowercase());
            }
            push(canon, &mut candidates, &mut seen);
        }
    }

    if let Some(chunk) = trimmed
        .split(['.', '-', '_'])
        .next()
        .filter(|c| !c.is_empty())
    {
        push(chunk.to_string(), &mut candidates, &mut seen);
    }

    let mut chars = trimmed.chars();
    if let Some(first) = chars.next() {
        if let Some(second) = chars.next().filter(|c| c.is_ascii_alphabetic()) {
            push(
                format!(
                    "{}{}",
                    first.to_ascii_uppercase(),
                    second.to_ascii_lowercase()
                ),
                &mut candidates,
                &mut seen,
            );
        }
        push(
            first.to_ascii_uppercase().to_string(),
            &mut candidates,
            &mut seen,
        );
    }

    push(trimmed.to_ascii_uppercase(), &mut candidates, &mut seen);

    for cand in candidates {
        if let Ok(el) = Element::from_str(&cand) {
            return Some(el);
        }
    }
    None
}

/// Converts a CT file bond order integer to a dreid-forge [`BondOrder`].
pub fn bond_order_from_ctfile(value: i32) -> Option<BondOrder> {
    match value {
        1 => Some(BondOrder::Single),
        2 => Some(BondOrder::Double),
        3 => Some(BondOrder::Triple),
        4 => Some(BondOrder::Aromatic),
        _ => None,
    }
}

/// Converts a dreid-forge [`BondOrder`] to a CT file bond order integer.
pub fn bond_order_to_ctfile(order: BondOrder) -> i32 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

/// Converts a MOL2 file bond order token to a dreid-forge [`BondOrder`].
pub fn bond_order_from_mol2(token: &str) -> Option<BondOrder> {
    match token.trim().to_ascii_lowercase().as_str() {
        "1" => Some(BondOrder::Single),
        "2" => Some(BondOrder::Double),
        "3" => Some(BondOrder::Triple),
        "ar" => Some(BondOrder::Aromatic),
        "am" | "du" | "un" | "nc" => Some(BondOrder::Single),
        _ => None,
    }
}

/// Converts a dreid-forge [`BondOrder`] to a MOL2 file bond order token.
pub fn bond_order_to_mol2(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "1",
        BondOrder::Double => "2",
        BondOrder::Triple => "3",
        BondOrder::Aromatic => "ar",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{Anion, Cation, CleanConfig, HisStrategy, ProtonationConfig, SolvateConfig};
    use crate::model::{
        atom::Atom,
        metadata::{
            AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition, StandardResidue,
        },
        system::{Bond, System},
        types::{BondOrder, Element},
    };
    use std::collections::HashSet;

    #[test]
    fn roundtrip_system_to_and_from_bio_preserves_metadata_and_bonds() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0]),
            Atom::new(Element::O, [1.0, 0.0, 0.0]),
            Atom::new(Element::N, [0.0, 1.0, 0.0]),
        ];

        let metadata = BioMetadata {
            atom_info: vec![
                AtomResidueInfo::builder("CA", "ALA", 1, 'A')
                    .standard_name(Some(StandardResidue::ALA))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("CB", "ALA", 1, 'A')
                    .standard_name(Some(StandardResidue::ALA))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build(),
                AtomResidueInfo::builder("N", "GLY", 2, 'B')
                    .insertion_code_opt(Some('B'))
                    .standard_name(Some(StandardResidue::GLY))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::NTerminal)
                    .build(),
            ],
        };

        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ];

        let box_vectors = Some([[10.0_f64, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]);

        let system = System {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            box_vectors,
            bio_metadata: Some(metadata.clone()),
        };

        let bio_topo = to_bio_topology(&system).expect("convert system to bio");
        let roundtrip = from_bio_topology(bio_topo).expect("convert bio back to system");

        assert_eq!(roundtrip.atoms, atoms);
        assert_eq!(roundtrip.bonds, bonds);
        assert_eq!(roundtrip.box_vectors, box_vectors);
        assert_eq!(roundtrip.bio_metadata.as_ref(), Some(&metadata));
    }

    #[test]
    fn converts_clean_config_fields() {
        let mut remove_res = HashSet::new();
        remove_res.insert("HOH".to_string());
        let mut keep_res = HashSet::new();
        keep_res.insert("LIG".to_string());

        let cfg = CleanConfig {
            remove_water: true,
            remove_ions: false,
            remove_hydrogens: true,
            remove_hetero: false,
            remove_residue_names: remove_res.clone(),
            keep_residue_names: keep_res.clone(),
        };

        let bf_cfg = to_bf_clean_config(cfg);
        assert!(bf_cfg.remove_water);
        assert!(!bf_cfg.remove_ions);
        assert!(bf_cfg.remove_hydrogens);
        assert!(!bf_cfg.remove_hetero);
        assert_eq!(bf_cfg.remove_residue_names, remove_res);
        assert_eq!(bf_cfg.keep_residue_names, keep_res);
    }

    #[test]
    fn converts_protonation_config_and_his_strategy() {
        let cfg = ProtonationConfig {
            target_ph: Some(7.4),
            remove_existing_h: false,
            his_strategy: HisStrategy::DirectHIE,
        };

        let bf_cfg = to_bf_hydro_config(cfg);
        assert_eq!(bf_cfg.target_ph, Some(7.4));
        assert!(!bf_cfg.remove_existing_h);
        assert!(matches!(
            bf_cfg.his_strategy,
            bf::ops::HisStrategy::DirectHIE
        ));
    }

    #[test]
    fn converts_solvate_config_and_ions() {
        let cfg = SolvateConfig {
            margin: 5.0,
            water_spacing: 2.8,
            vdw_cutoff: 2.0,
            remove_existing: false,
            cations: vec![Cation::K, Cation::Zn],
            anions: vec![Anion::Br],
            target_charge: -1,
            rng_seed: Some(99),
        };

        let bf_cfg = to_bf_solvate_config(cfg);
        assert_eq!(bf_cfg.margin, 5.0);
        assert_eq!(bf_cfg.water_spacing, 2.8);
        assert_eq!(bf_cfg.vdw_cutoff, 2.0);
        assert!(!bf_cfg.remove_existing);
        assert_eq!(
            bf_cfg.cations,
            vec![bf::ops::Cation::K, bf::ops::Cation::Zn]
        );
        assert_eq!(bf_cfg.anions, vec![bf::ops::Anion::Br]);
        assert_eq!(bf_cfg.target_charge, -1);
        assert_eq!(bf_cfg.rng_seed, Some(99));
    }

    #[test]
    fn rejects_unknown_element_from_bio() {
        let mut residue = bf::Residue::new(1, None, "UNK", None, bf::ResidueCategory::Hetero);
        residue.position = bf::ResiduePosition::None;
        residue.add_atom(bf::Atom::new(
            "X",
            bf::Element::Unknown,
            bf::Point::new(0.0, 0.0, 0.0),
        ));

        let mut chain = bf::Chain::new("A");
        chain.add_residue(residue);

        let mut structure = bf::Structure::new();
        structure.add_chain(chain);

        let topo = bf::Topology::new(structure, Vec::new());

        let err = from_bio_topology(topo).unwrap_err();
        assert!(matches!(err, ConversionError::UnsupportedElement));
    }

    #[test]
    fn guesses_elements_from_common_tokens() {
        assert_eq!(guess_element_symbol("C.3"), Some(Element::C));
        assert_eq!(guess_element_symbol("cl"), Some(Element::Cl));
        assert_eq!(guess_element_symbol("BR"), Some(Element::Br));
        assert_eq!(guess_element_symbol("CA1"), Some(Element::Ca));
        assert_eq!(guess_element_symbol("   n   "), Some(Element::N));
    }

    #[test]
    fn converts_ctfile_bond_orders() {
        assert_eq!(bond_order_from_ctfile(1), Some(BondOrder::Single));
        assert_eq!(bond_order_from_ctfile(4), Some(BondOrder::Aromatic));
        assert_eq!(bond_order_from_ctfile(7), None);
        assert_eq!(bond_order_to_ctfile(BondOrder::Triple), 3);
    }

    #[test]
    fn converts_mol2_bond_orders() {
        assert_eq!(bond_order_from_mol2("1"), Some(BondOrder::Single));
        assert_eq!(bond_order_from_mol2("ar"), Some(BondOrder::Aromatic));
        assert_eq!(bond_order_from_mol2("am"), Some(BondOrder::Single));
        assert_eq!(bond_order_from_mol2("xx"), None);
        assert_eq!(bond_order_to_mol2(BondOrder::Double), "2");
    }
}
