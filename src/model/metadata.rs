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
        insertion_code: Option<char>
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
