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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::system::Bond;

    fn make_water() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::O, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [0.96, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.24, 0.93, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys
    }

    fn make_ethane() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::C, [1.54, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, 0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, 0.89]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 3, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 4, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 5, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 6, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 7, BondOrder::Single));
        sys
    }

    #[test]
    fn intermediate_from_water() {
        let water = make_water();
        let int = IntermediateSystem::from_system(&water).unwrap();

        assert_eq!(int.atoms.len(), 3);
        assert_eq!(int.bonds.len(), 2);

        assert!(int.angles.is_empty());
        assert!(int.dihedrals.is_empty());
        assert!(int.impropers.is_empty());

        assert_eq!(int.atoms[0].neighbors.len(), 2);
        assert_eq!(int.atoms[1].neighbors.len(), 1);
        assert_eq!(int.atoms[2].neighbors.len(), 1);
    }

    #[test]
    fn intermediate_from_ethane() {
        let ethane = make_ethane();
        let int = IntermediateSystem::from_system(&ethane).unwrap();

        assert_eq!(int.atoms.len(), 8);
        assert_eq!(int.bonds.len(), 7);

        assert!(int.angles.is_empty());
        assert!(int.dihedrals.is_empty());
        assert!(int.impropers.is_empty());
    }

    #[test]
    fn atom_view_implementation() {
        let atom = IntermediateAtom::new(Element::C, [1.0, 2.0, 3.0]);
        assert_eq!(atom.atomic_number(), 6);
        assert_eq!(atom.position(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn physical_bond_order_multiplier() {
        assert_eq!(PhysicalBondOrder::Single.multiplier(), 1.0);
        assert_eq!(PhysicalBondOrder::Double.multiplier(), 2.0);
        assert_eq!(PhysicalBondOrder::Triple.multiplier(), 3.0);
        assert_eq!(PhysicalBondOrder::Resonant.multiplier(), 1.5);
    }

    #[test]
    fn physical_order_is_none_before_typing() {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::C, [1.5, 0.0, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Aromatic));

        let int = IntermediateSystem::from_system(&sys).unwrap();
        assert!(int.bonds[0].physical_order.is_none());
    }

    #[test]
    fn errors_on_empty_system() {
        let empty = System::new();
        let result = IntermediateSystem::from_system(&empty);
        assert!(matches!(result, Err(Error::EmptySystem)));
    }

    #[test]
    fn errors_on_invalid_bond_index() {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.bonds.push(Bond::new(0, 99, BondOrder::Single));

        let result = IntermediateSystem::from_system(&sys);
        assert!(matches!(
            result,
            Err(Error::InvalidBond { i: 0, j: 99, .. })
        ));
    }

    #[test]
    fn intermediate_atom_new_initializes_defaults() {
        let atom = IntermediateAtom::new(Element::N, [1.0, 2.0, 3.0]);
        assert_eq!(atom.element, Element::N);
        assert_eq!(atom.atomic_number(), 7);
        assert_eq!(atom.position, [1.0, 2.0, 3.0]);
        assert!(atom.atom_type.is_empty());
        assert_eq!(atom.charge, 0.0);
        assert!(atom.neighbors.is_empty());
    }

    #[test]
    fn intermediate_bond_preserves_order() {
        let bond = IntermediateBond {
            i: 0,
            j: 1,
            order: BondOrder::Double,
            physical_order: Some(PhysicalBondOrder::Double),
        };
        assert_eq!(bond.order, BondOrder::Double);
        assert_eq!(bond.physical_order, Some(PhysicalBondOrder::Double));
    }
}
