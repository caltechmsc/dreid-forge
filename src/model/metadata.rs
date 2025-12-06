#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomResidueInfo {
    pub atom_name: String,
    pub residue_name: String,
    pub residue_id: i32,
    pub chain_id: char,
    pub insertion_code: char,
}

impl AtomResidueInfo {
    pub fn new(
        atom_name: impl Into<String>,
        residue_name: impl Into<String>,
        residue_id: i32,
        chain_id: char,
        insertion_code: Option<char>,
    ) -> Self {
        Self {
            atom_name: atom_name.into(),
            residue_name: residue_name.into(),
            residue_id,
            chain_id,
            insertion_code: insertion_code.unwrap_or(' '),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BioMetadata {
    pub atom_info: Vec<AtomResidueInfo>,
}

impl BioMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            atom_info: Vec::with_capacity(capacity),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_residue_info_new_and_fields() {
        let info = AtomResidueInfo::new("CA", "ALA", 42, 'A', Some('x'));
        assert_eq!(info.atom_name, "CA");
        assert_eq!(info.residue_name, "ALA");
        assert_eq!(info.residue_id, 42);
        assert_eq!(info.chain_id, 'A');
        assert_eq!(info.insertion_code, 'x');
    }

    #[test]
    fn atom_residue_info_default_insertion_code() {
        let info = AtomResidueInfo::new("N", "GLY", 1, 'B', None);
        assert_eq!(info.insertion_code, ' ');
    }

    #[test]
    fn atom_residue_info_clone_and_eq() {
        let a = AtomResidueInfo::new("O", "HOH", 5, 'C', None);
        let b = a.clone();
        assert_eq!(a, b);
    }

    #[test]
    fn bio_metadata_new_and_default() {
        let bm = BioMetadata::new();
        assert!(bm.atom_info.is_empty());
        assert_eq!(bm, BioMetadata::default());
    }

    #[test]
    fn bio_metadata_with_capacity_and_push() {
        let mut bm = BioMetadata::with_capacity(4);
        assert!(bm.atom_info.capacity() >= 4);

        let info1 = AtomResidueInfo::new("CA", "ALA", 1, 'A', None);
        let info2 = AtomResidueInfo::new("CB", "ALA", 1, 'A', None);
        bm.atom_info.push(info1.clone());
        bm.atom_info.push(info2.clone());

        assert_eq!(bm.atom_info.len(), 2);
        assert_eq!(bm.atom_info[0], info1);
        assert_eq!(bm.atom_info[1], info2);
    }

    #[test]
    fn debug_contains_expected_fields() {
        let info = AtomResidueInfo::new("C1", "LIG", -1, 'Z', Some('A'));
        let bm = BioMetadata {
            atom_info: vec![info.clone()],
        };
        let s_info = format!("{:?}", info);
        let s_bm = format!("{:?}", bm);
        assert!(s_info.contains("atom_name"));
        assert!(s_info.contains("residue_name"));
        assert!(s_info.contains("residue_id"));
        assert!(s_info.contains("chain_id"));
        assert!(s_info.contains("insertion_code"));
        assert!(s_bm.contains("AtomResidueInfo"));
        assert!(s_bm.contains("LIG"));
    }
}
