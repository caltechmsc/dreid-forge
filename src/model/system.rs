//! Molecular system representation with atoms, bonds, and metadata.
//!
//! This module provides the core data structures for representing
//! complete molecular systems:
//!
//! - [`Bond`] — Connection between two atoms with associated bond order
//! - [`System`] — Complete molecular system with atoms, bonds, periodicity,
//!   and optional biological metadata

use super::atom::Atom;
use super::metadata::BioMetadata;
use super::types::BondOrder;

/// A chemical bond between two atoms.
///
/// Represents a covalent bond with canonical ordering (i ≤ j) to ensure
/// consistent hashing and equality comparisons. This allows bonds to be
/// stored in hash-based collections without duplicates regardless of the
/// order atoms are specified.
///
/// # Fields
///
/// * `i` — Index of the first atom (always ≤ j)
/// * `j` — Index of the second atom (always ≥ i)
/// * `order` — Bond multiplicity ([`BondOrder`])
///
/// # Examples
///
/// ```
/// use dreid_forge::{Bond, BondOrder};
///
/// // Order is automatically canonicalized
/// let bond1 = Bond::new(5, 2, BondOrder::Single);
/// let bond2 = Bond::new(2, 5, BondOrder::Single);
///
/// assert_eq!(bond1.i, 2); // smaller index
/// assert_eq!(bond1.j, 5); // larger index
/// assert_eq!(bond1, bond2); // equivalent regardless of input order
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    /// Index of the first atom (canonically the smaller index).
    pub i: usize,
    /// Index of the second atom (canonically the larger or equal index).
    pub j: usize,
    /// Bond order (single, double, triple, or aromatic).
    pub order: BondOrder,
}

impl Bond {
    /// Creates a new bond with canonical atom ordering.
    ///
    /// Automatically orders atom indices so that `i ≤ j`, ensuring
    /// consistent representation regardless of the order arguments
    /// are provided.
    ///
    /// # Arguments
    ///
    /// * `idx1` — Index of the first atom
    /// * `idx2` — Index of the second atom
    /// * `order` — Bond multiplicity
    ///
    /// # Returns
    ///
    /// A new [`Bond`] with canonically ordered indices.
    pub fn new(idx1: usize, idx2: usize, order: BondOrder) -> Self {
        if idx1 <= idx2 {
            Self {
                i: idx1,
                j: idx2,
                order,
            }
        } else {
            Self {
                i: idx2,
                j: idx1,
                order,
            }
        }
    }
}

/// A complete molecular system.
///
/// Contains all information needed to represent a molecular structure:
/// atoms with their positions, connectivity through bonds, optional
/// periodic boundary conditions, and optional biological metadata for
/// biomolecular systems.
///
/// # Fields
///
/// * `atoms` — Vector of atoms with element types and positions
/// * `bonds` — Vector of bonds defining molecular connectivity
/// * `box_vectors` — Optional 3×3 matrix defining periodic cell vectors
/// * `bio_metadata` — Optional biological annotation (residues, chains)
///
/// # Examples
///
/// ```
/// use dreid_forge::{Atom, Bond, BondOrder, Element, System};
///
/// // Build a simple CO molecule
/// let mut system = System::new();
/// system.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
/// system.atoms.push(Atom::new(Element::O, [1.128, 0.0, 0.0]));
/// system.bonds.push(Bond::new(0, 1, BondOrder::Triple));
///
/// assert_eq!(system.atom_count(), 2);
/// assert_eq!(system.bond_count(), 1);
/// assert!(!system.is_periodic());
/// ```
#[derive(Debug, Clone, Default)]
pub struct System {
    /// Atoms in the system with element types and positions.
    pub atoms: Vec<Atom>,
    /// Bonds defining molecular connectivity.
    pub bonds: Vec<Bond>,
    /// Periodic cell vectors as row-major 3×3 matrix, if periodic.
    pub box_vectors: Option<[[f64; 3]; 3]>,
    /// Biological metadata (residue info, chains) for biomolecules.
    pub bio_metadata: Option<BioMetadata>,
}

impl System {
    /// Creates a new empty molecular system.
    ///
    /// Initializes a system with no atoms, no bonds, no periodic
    /// boundaries, and no biological metadata.
    ///
    /// # Returns
    ///
    /// An empty [`System`] instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the number of atoms in the system.
    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Returns the number of bonds in the system.
    #[inline]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Returns `true` if the system has periodic boundary conditions.
    ///
    /// A system is periodic if [`box_vectors`](Self::box_vectors) is `Some`.
    #[inline]
    pub fn is_periodic(&self) -> bool {
        self.box_vectors.is_some()
    }

    /// Returns `true` if the system has biological metadata attached.
    ///
    /// Biological metadata includes residue names, chain IDs, and
    /// other annotations from PDB/mmCIF files.
    #[inline]
    pub fn has_bio_metadata(&self) -> bool {
        self.bio_metadata.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::metadata::BioMetadata;
    use crate::model::types::BondOrder;
    use crate::model::types::Element;
    use std::collections::HashSet;

    #[test]
    fn bond_new_canonical_order() {
        let b1 = Bond::new(5, 2, BondOrder::Single);
        assert_eq!(b1.i, 2);
        assert_eq!(b1.j, 5);

        let b2 = Bond::new(2, 5, BondOrder::Single);
        assert_eq!(b1, b2);
    }

    #[test]
    fn bond_hashset_dedupes_equivalent_bonds() {
        let b1 = Bond::new(3, 1, BondOrder::Double);
        let b2 = Bond::new(1, 3, BondOrder::Double);
        let mut hs = HashSet::new();
        hs.insert(b1);
        hs.insert(b2);
        assert_eq!(hs.len(), 1);
    }

    #[test]
    fn system_new_and_counts() {
        let s = System::new();
        assert_eq!(s.atom_count(), 0);
        assert_eq!(s.bond_count(), 0);
        assert!(!s.is_periodic());
        assert!(!s.has_bio_metadata());
    }

    #[test]
    fn system_add_atoms_and_bonds_and_flags() {
        let mut s = System::new();
        s.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        s.atoms.push(Atom::new(Element::O, [1.2, 0.0, 0.0]));
        assert_eq!(s.atom_count(), 2);

        s.bonds.push(Bond::new(0, 1, BondOrder::Double));
        assert_eq!(s.bond_count(), 1);

        assert!(!s.is_periodic());
        s.box_vectors = Some([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]);
        assert!(s.is_periodic());

        assert!(!s.has_bio_metadata());
        s.bio_metadata = Some(BioMetadata::new());
        assert!(s.has_bio_metadata());
    }

    #[test]
    fn debug_and_clone_system() {
        let mut s = System::new();
        s.atoms.push(Atom::new(Element::H, [0.0, 0.0, 0.0]));
        s.bonds.push(Bond::new(0, 0, BondOrder::Single));
        let dbg = format!("{:?}", s.clone());
        assert!(dbg.contains("System"));
        assert!(dbg.contains("atoms"));
        assert!(dbg.contains("bonds"));
    }
}
