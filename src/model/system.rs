use super::atom::Atom;
use super::metadata::BioMetadata;
use super::types::BondOrder;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub order: BondOrder,
}

impl Bond {
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

#[derive(Debug, Clone, Default)]
pub struct System {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub box_vectors: Option<[[f64; 3]; 3]>,
    pub bio_metadata: Option<BioMetadata>,
}

impl System {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    #[inline]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    #[inline]
    pub fn is_periodic(&self) -> bool {
        self.box_vectors.is_some()
    }

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
