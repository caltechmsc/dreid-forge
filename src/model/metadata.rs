#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StandardResidue {
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL,
    A,
    C,
    G,
    U,
    I,
    DA,
    DC,
    DG,
    DT,
    DI,
    HOH,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueCategory {
    Standard,
    Hetero,
    Ion,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResiduePosition {
    None,
    Internal,
    NTerminal,
    CTerminal,
    FivePrime,
    ThreePrime,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomResidueInfo {
    pub atom_name: String,
    pub residue_name: String,
    pub residue_id: i32,
    pub chain_id: char,
    pub insertion_code: char,
    pub standard_name: Option<StandardResidue>,
    pub category: ResidueCategory,
    pub position: ResiduePosition,
}

impl AtomResidueInfo {
    pub fn new(
        atom_name: impl Into<String>,
        residue_name: impl Into<String>,
        residue_id: i32,
        chain_id: char,
        insertion_code: Option<char>,
        standard_name: Option<StandardResidue>,
        category: ResidueCategory,
        position: ResiduePosition,
    ) -> Self {
        Self {
            atom_name: atom_name.into(),
            residue_name: residue_name.into(),
            residue_id,
            chain_id,
            insertion_code: insertion_code.unwrap_or(' '),
            standard_name,
            category,
            position,
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
    fn atom_residue_info_new_and_all_fields() {
        let info = AtomResidueInfo::new(
            "CA",
            "ALA",
            42,
            'A',
            Some('x'),
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
            ResiduePosition::Internal,
        );
        assert_eq!(info.atom_name, "CA");
        assert_eq!(info.residue_name, "ALA");
        assert_eq!(info.residue_id, 42);
        assert_eq!(info.chain_id, 'A');
        assert_eq!(info.insertion_code, 'x');
        assert_eq!(info.standard_name, Some(StandardResidue::ALA));
        assert_eq!(info.category, ResidueCategory::Standard);
        assert_eq!(info.position, ResiduePosition::Internal);
    }

    #[test]
    fn atom_residue_info_defaults_and_clone() {
        let info = AtomResidueInfo::new(
            "N",
            "GLY",
            1,
            'B',
            None,
            None,
            ResidueCategory::Hetero,
            ResiduePosition::None,
        );
        assert_eq!(info.insertion_code, ' ');
        let cloned = info.clone();
        assert_eq!(info, cloned);
    }

    #[test]
    fn bio_metadata_new_and_capacity() {
        let mut bm = BioMetadata::with_capacity(4);
        assert!(bm.atom_info.capacity() >= 4);

        let info1 = AtomResidueInfo::new(
            "CA",
            "ALA",
            1,
            'A',
            None,
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
            ResiduePosition::Internal,
        );
        let info2 = AtomResidueInfo::new(
            "CB",
            "ALA",
            1,
            'A',
            None,
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
            ResiduePosition::Internal,
        );
        bm.atom_info.push(info1.clone());
        bm.atom_info.push(info2.clone());

        assert_eq!(bm.atom_info.len(), 2);
        assert_eq!(bm.atom_info[0], info1);
        assert_eq!(bm.atom_info[1], info2);
    }

    #[test]
    fn debug_contains_expected_fields() {
        let info = AtomResidueInfo::new(
            "C1",
            "LIG",
            -1,
            'Z',
            Some('A'),
            None,
            ResidueCategory::Hetero,
            ResiduePosition::NTerminal,
        );
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
        assert!(s_info.contains("standard_name") || s_info.contains("StandardResidue"));
        assert!(s_info.contains("category") || s_info.contains("ResidueCategory"));
        assert!(s_info.contains("position") || s_info.contains("ResiduePosition"));
        assert!(s_bm.contains("AtomResidueInfo"));
        assert!(s_bm.contains("LIG"));
    }
}
