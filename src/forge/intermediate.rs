use super::error::Error;
use crate::model::system::System;
use crate::model::types::{BondOrder, Element};
use cheq::AtomView;

pub use dreid_typer::{Hybridization, TopologyBondOrder as PhysicalBondOrder};

#[derive(Debug, Clone)]
pub struct IntermediateAtom {
    pub element: Element,
    pub position: [f64; 3],
    pub atom_type: String,
    pub hybridization: Hybridization,
    pub charge: f64,
    pub neighbors: Vec<usize>,
}

impl IntermediateAtom {
    pub fn new(element: Element, position: [f64; 3]) -> Self {
        Self {
            element,
            position,
            atom_type: String::new(),
            hybridization: Hybridization::Unknown,
            charge: 0.0,
            neighbors: Vec::new(),
        }
    }
}

impl AtomView for IntermediateAtom {
    #[inline]
    fn atomic_number(&self) -> u8 {
        self.element.atomic_number()
    }

    #[inline]
    fn position(&self) -> [f64; 3] {
        self.position
    }
}

#[derive(Debug, Clone)]
pub struct IntermediateBond {
    pub i: usize,
    pub j: usize,
    pub order: BondOrder,
    pub physical_order: Option<PhysicalBondOrder>,
}

pub trait PhysicalBondOrderExt {
    fn multiplier(&self) -> f64;
}

impl PhysicalBondOrderExt for PhysicalBondOrder {
    fn multiplier(&self) -> f64 {
        match self {
            PhysicalBondOrder::Single => 1.0,
            PhysicalBondOrder::Double => 2.0,
            PhysicalBondOrder::Triple => 3.0,
            PhysicalBondOrder::Resonant => 1.5,
        }
    }
}

#[derive(Debug, Clone)]
pub struct IntermediateAngle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

#[derive(Debug, Clone)]
pub struct IntermediateDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
}

#[derive(Debug, Clone)]
pub struct IntermediateImproper {
    pub center: usize,
    pub p1: usize,
    pub p2: usize,
    pub p3: usize,
}

#[derive(Debug, Clone)]
pub struct IntermediateSystem {
    pub atoms: Vec<IntermediateAtom>,
    pub bonds: Vec<IntermediateBond>,
    pub angles: Vec<IntermediateAngle>,
    pub dihedrals: Vec<IntermediateDihedral>,
    pub impropers: Vec<IntermediateImproper>,
}

impl IntermediateSystem {
    pub fn from_system(system: &System) -> Result<Self, Error> {
        if system.atoms.is_empty() {
            return Err(Error::EmptySystem);
        }

        let n_atoms = system.atoms.len();

        let mut atoms: Vec<IntermediateAtom> = system
            .atoms
            .iter()
            .map(|a| IntermediateAtom::new(a.element, a.position))
            .collect();

        let mut bonds = Vec::with_capacity(system.bonds.len());
        for bond in &system.bonds {
            if bond.i >= n_atoms || bond.j >= n_atoms {
                return Err(Error::invalid_bond(
                    bond.i,
                    bond.j,
                    format!("atom index out of bounds (n_atoms = {})", n_atoms),
                ));
            }

            let int_bond = IntermediateBond {
                i: bond.i,
                j: bond.j,
                order: bond.order,
                physical_order: None,
            };
            bonds.push(int_bond);

            atoms[bond.i].neighbors.push(bond.j);
            atoms[bond.j].neighbors.push(bond.i);
        }

        Ok(Self {
            atoms,
            bonds,
            angles: Vec::new(),
            dihedrals: Vec::new(),
            impropers: Vec::new(),
        })
    }

    pub fn atoms(&self) -> &[IntermediateAtom] {
        &self.atoms
    }
}
