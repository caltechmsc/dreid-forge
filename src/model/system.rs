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
            Self { i: idx1, j: idx2, order }
        } else {
            Self { i: idx2, j: idx1, order }
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
