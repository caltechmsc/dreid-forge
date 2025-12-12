//! Intermediate representations for parameterization pipeline.
//!
//! This module defines the internal data structures used during the
//! DREIDING parameterization process. The [`IntermediateSystem`] holds
//! atoms with assigned types and charges, bonds with physical orders,
//! and enumerated angles, dihedrals, and impropers.
//!
//! These structures are not part of the public API and are used only
//! within the forge module to pass data between pipeline stages.

use super::error::Error;
use crate::model::system::System;
use crate::model::types::{BondOrder, Element};
use cheq::AtomView;

pub use dreid_typer::{Hybridization, TopologyBondOrder as PhysicalBondOrder};

/// Intermediate atom representation with typing and charge information.
///
/// Extends the basic atom data with fields populated during parameterization:
/// atom type string, hybridization state, partial charge, and neighbor list.
#[derive(Debug, Clone)]
pub struct IntermediateAtom {
    /// Chemical element.
    pub element: Element,
    /// Cartesian coordinates in Ångströms.
    pub position: [f64; 3],
    /// DREIDING atom type (e.g., "C_3", "O_2").
    pub atom_type: String,
    /// Hybridization state (SP, SP2, SP3, Resonant, etc.).
    pub hybridization: Hybridization,
    /// Partial atomic charge in elementary charge units.
    pub charge: f64,
    /// Indices of bonded neighbor atoms.
    pub neighbors: Vec<usize>,
}

impl IntermediateAtom {
    /// Creates a new intermediate atom with default (empty) type and zero charge.
    ///
    /// # Arguments
    ///
    /// * `element` — Chemical element
    /// * `position` — Cartesian coordinates [x, y, z] in Ångströms
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

/// Intermediate bond representation with physical bond order.
///
/// Stores both the original [`BondOrder`] from the input system and
/// the physical bond order assigned by the typer (Single, Double,
/// Triple, or Resonant).
#[derive(Debug, Clone)]
pub struct IntermediateBond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Original bond order from input.
    pub order: BondOrder,
    /// Physical bond order assigned by typer.
    pub physical_order: Option<PhysicalBondOrder>,
}

/// Extension trait for physical bond order multiplier.
pub trait PhysicalBondOrderExt {
    /// Returns the bond order multiplier for force constant scaling.
    ///
    /// - Single → 1.0
    /// - Double → 2.0
    /// - Triple → 3.0
    /// - Resonant → 1.5
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

/// Intermediate angle (three-atom bend) representation.
#[derive(Debug, Clone)]
pub struct IntermediateAngle {
    /// First atom index.
    pub i: usize,
    /// Central atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
}

/// Intermediate proper dihedral (four-atom torsion) representation.
#[derive(Debug, Clone)]
pub struct IntermediateDihedral {
    /// First atom index.
    pub i: usize,
    /// Second atom index (part of central bond).
    pub j: usize,
    /// Third atom index (part of central bond).
    pub k: usize,
    /// Fourth atom index.
    pub l: usize,
}

/// Intermediate improper dihedral (out-of-plane) representation.
#[derive(Debug, Clone)]
pub struct IntermediateImproper {
    /// Central atom (typically sp² center).
    pub center: usize,
    /// First peripheral atom.
    pub p1: usize,
    /// Second peripheral atom.
    pub p2: usize,
    /// Third peripheral atom.
    pub p3: usize,
}

/// Complete intermediate system for parameterization pipeline.
///
/// Contains all atoms, bonds, and enumerated internal coordinates
/// (angles, dihedrals, impropers) needed for parameter generation.
#[derive(Debug, Clone)]
pub struct IntermediateSystem {
    /// All atoms with typing and charge information.
    pub atoms: Vec<IntermediateAtom>,
    /// All bonds with physical bond orders.
    pub bonds: Vec<IntermediateBond>,
    /// All angle bend terms.
    pub angles: Vec<IntermediateAngle>,
    /// All proper dihedral (torsion) terms.
    pub dihedrals: Vec<IntermediateDihedral>,
    /// All improper dihedral (out-of-plane) terms.
    pub impropers: Vec<IntermediateImproper>,
}

impl IntermediateSystem {
    /// Creates an intermediate system from a molecular [`System`].
    ///
    /// Converts atoms and bonds from the input system, building the
    /// neighbor lists for each atom. Angles, dihedrals, and impropers
    /// are left empty to be populated by the typer.
    ///
    /// # Arguments
    ///
    /// * `system` — The molecular system to convert
    ///
    /// # Errors
    ///
    /// Returns [`Error::EmptySystem`] if the system has no atoms, or
    /// [`Error::InvalidBond`] if a bond references out-of-bounds indices.
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

    /// Returns a slice of all intermediate atoms.
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
