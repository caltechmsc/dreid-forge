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
use std::str::FromStr;

#[derive(Debug, thiserror::Error)]
pub enum ConversionError {
    #[error("encountered an 'Unknown' element from bio-forge which is not supported")]
    UnsupportedElement,
    #[error("encountered an unsupported bond order '{0}' during model conversion")]
    UnsupportedBondOrder(String),
    #[error("inconsistent data: system is missing required BioMetadata for this conversion")]
    MissingBioMetadata,
    #[error("unknown standard residue name")]
    UnknownStandardResidue,
    #[error("unknown residue category")]
    UnknownResidueCategory,
    #[error("unknown residue position")]
    UnknownResiduePosition,
}

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

pub fn to_bf_anion(anion: Anion) -> bf::ops::Anion {
    match anion {
        Anion::Cl => bf::ops::Anion::Cl,
        Anion::Br => bf::ops::Anion::Br,
        Anion::I => bf::ops::Anion::I,
        Anion::F => bf::ops::Anion::F,
    }
}

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

pub fn to_bio_topology(system: &System) -> Result<bf::Topology, ConversionError> {
    use std::collections::BTreeMap;

    let metadata = system
        .bio_metadata
        .as_ref()
        .ok_or(ConversionError::MissingBioMetadata)?;

    let bio_atoms: Vec<bf::Atom> = system
        .atoms
        .iter()
        .zip(metadata.atom_info.iter())
        .map(|(atom, info)| {
            Ok(bf::Atom::new(
                &info.atom_name,
                convert_element_to_bf(atom.element)?,
                bf::Point::new(atom.position[0], atom.position[1], atom.position[2]),
            ))
        })
        .collect::<Result<_, _>>()?;

    type ResidueKey = (char, i32, char);
    let mut residue_atom_indices: BTreeMap<ResidueKey, Vec<usize>> = BTreeMap::new();
    for (i, info) in metadata.atom_info.iter().enumerate() {
        let key = (info.chain_id, info.residue_id, info.insertion_code);
        residue_atom_indices.entry(key).or_default().push(i);
    }

    let mut residues: BTreeMap<ResidueKey, bf::Residue> = BTreeMap::new();
    for (key, mut indices) in residue_atom_indices {
        indices.sort_by(|&a, &b| {
            let ia = &metadata.atom_info[a];
            let ib = &metadata.atom_info[b];
            ia.atom_name.cmp(&ib.atom_name).then_with(|| a.cmp(&b))
        });

        let first_info = &metadata.atom_info[indices[0]];
        let mut residue = bf::Residue::new(
            first_info.residue_id,
            Some(first_info.insertion_code).filter(|&c| c != ' '),
            &first_info.residue_name,
            convert_std_res_to_bf(first_info.standard_name)?,
            convert_res_cat_to_bf(first_info.category)?,
        );
        residue.position = convert_res_pos_to_bf(first_info.position)?;

        for idx in indices {
            residue.add_atom(bio_atoms[idx].clone());
        }
        residues.insert(key, residue);
    }

    let mut chains: BTreeMap<char, bf::Chain> = BTreeMap::new();
    for ((chain_id, _, _), residue) in residues {
        chains
            .entry(chain_id)
            .or_insert_with(|| bf::Chain::new(&chain_id.to_string()))
            .add_residue(residue);
    }

    let mut bio_struct = bf::Structure::new();
    bio_struct.box_vectors = system.box_vectors;
    for (_, chain) in chains {
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

fn convert_element_from_bf(e: bf::Element) -> Result<Element, ConversionError> {
    if matches!(e, bf::Element::Unknown) {
        return Err(ConversionError::UnsupportedElement);
    }

    Element::from_str(e.symbol()).map_err(|_| ConversionError::UnsupportedElement)
}

fn convert_element_to_bf(e: Element) -> Result<bf::Element, ConversionError> {
    bf::Element::from_str(e.symbol()).map_err(|_| ConversionError::UnsupportedElement)
}

fn convert_bond_order_from_bf(order: bf::BondOrder) -> Result<BondOrder, ConversionError> {
    match order {
        bf::BondOrder::Single => Ok(BondOrder::Single),
        bf::BondOrder::Double => Ok(BondOrder::Double),
        bf::BondOrder::Triple => Ok(BondOrder::Triple),
        bf::BondOrder::Aromatic => Ok(BondOrder::Aromatic),
    }
}
fn convert_bond_order_to_bf(order: BondOrder) -> Result<bf::BondOrder, ConversionError> {
    match order {
        BondOrder::Single => Ok(bf::BondOrder::Single),
        BondOrder::Double => Ok(bf::BondOrder::Double),
        BondOrder::Triple => Ok(bf::BondOrder::Triple),
        BondOrder::Aromatic => Ok(bf::BondOrder::Aromatic),
    }
}

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

fn convert_res_cat_from_bf(cat: bf::ResidueCategory) -> Result<ResidueCategory, ConversionError> {
    Ok(match cat {
        bf::ResidueCategory::Standard => ResidueCategory::Standard,
        bf::ResidueCategory::Hetero => ResidueCategory::Hetero,
        bf::ResidueCategory::Ion => ResidueCategory::Ion,
    })
}
fn convert_res_cat_to_bf(cat: ResidueCategory) -> Result<bf::ResidueCategory, ConversionError> {
    Ok(match cat {
        ResidueCategory::Standard => bf::ResidueCategory::Standard,
        ResidueCategory::Hetero => bf::ResidueCategory::Hetero,
        ResidueCategory::Ion => bf::ResidueCategory::Ion,
    })
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        metadata::{
            AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition, StandardResidue,
        },
        system::{Bond, System},
        types::{BondOrder, Element},
    };

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
}
