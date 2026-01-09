//! Biological metadata for biomolecular systems.
//!
//! This module provides data structures for storing PDB/mmCIF-style
//! annotations that describe the biological context of atoms within
//! protein, nucleic acid, and other biomolecular structures.
//!
//! # Overview
//!
//! - [`StandardResidue`] — Enumeration of standard amino acids, nucleotides, and water
//! - [`ResidueCategory`] — Classification as standard, hetero, or ion
//! - [`ResiduePosition`] — Terminal position annotations (N/C-terminal, 5'/3'-end)
//! - [`AtomResidueInfo`] — Per-atom biological context (residue name, chain, etc.)
//! - [`BioMetadata`] — Collection of atom-level biological annotations
//!
//! # Builder Pattern
//!
//! Use [`AtomResidueInfo::builder`] to construct instances with optional fields:
//!
//! ```
//! use dreid_forge::{AtomResidueInfo, StandardResidue, ResidueCategory, ResiduePosition};
//!
//! let info = AtomResidueInfo::builder("CA", "ALA", 42, 'A')
//!     .standard_name(Some(StandardResidue::ALA))
//!     .category(ResidueCategory::Standard)
//!     .position(ResiduePosition::NTerminal)
//!     .build();
//!
//! assert_eq!(info.atom_name, "CA");
//! assert_eq!(info.chain_id, 'A');
//! ```

/// Standard residue types from PDB/mmCIF nomenclature.
///
/// Covers the 20 canonical amino acids, common nucleotides (RNA and DNA),
/// and water (HOH).
///
/// # Amino Acids
///
/// Three-letter codes for all 20 standard amino acids:
/// ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE,
/// LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL.
///
/// # Nucleotides
///
/// - RNA: A, C, G, U, I (inosine)
/// - DNA: DA, DC, DG, DT, DI
///
/// # Solvent
///
/// - HOH: Water molecule
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StandardResidue {
    /// Alanine
    ALA,
    /// Arginine
    ARG,
    /// Asparagine
    ASN,
    /// Aspartic acid
    ASP,
    /// Cysteine
    CYS,
    /// Glutamine
    GLN,
    /// Glutamic acid
    GLU,
    /// Glycine
    GLY,
    /// Histidine
    HIS,
    /// Isoleucine
    ILE,
    /// Leucine
    LEU,
    /// Lysine
    LYS,
    /// Methionine
    MET,
    /// Phenylalanine
    PHE,
    /// Proline
    PRO,
    /// Serine
    SER,
    /// Threonine
    THR,
    /// Tryptophan
    TRP,
    /// Tyrosine
    TYR,
    /// Valine
    VAL,
    /// Adenosine (RNA)
    A,
    /// Cytidine (RNA)
    C,
    /// Guanosine (RNA)
    G,
    /// Uridine (RNA)
    U,
    /// Inosine (RNA)
    I,
    /// Deoxyadenosine (DNA)
    DA,
    /// Deoxycytidine (DNA)
    DC,
    /// Deoxyguanosine (DNA)
    DG,
    /// Deoxythymidine (DNA)
    DT,
    /// Deoxyinosine (DNA)
    DI,
    /// Water molecule
    HOH,
}

/// Classification of a residue's chemical nature.
///
/// Used to distinguish between standard biomolecular components
/// and non-standard entities like ligands or ions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueCategory {
    /// Standard amino acid or nucleotide.
    Standard,
    /// Hetero atom group (ligand, modified residue, etc.).
    Hetero,
    /// Monoatomic or polyatomic ion.
    Ion,
}

/// Position of a residue within its polymer chain.
///
/// Indicates whether the residue is at a terminal position,
/// which affects protonation states and hydrogen bonding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResiduePosition {
    /// No special terminal position (or not applicable).
    None,
    /// Internal residue within the chain.
    Internal,
    /// N-terminal residue of a protein chain.
    NTerminal,
    /// C-terminal residue of a protein chain.
    CTerminal,
    /// 5' end of a nucleic acid strand.
    FivePrime,
    /// 3' end of a nucleic acid strand.
    ThreePrime,
}

/// Biological annotation for a single atom.
///
/// Stores PDB/mmCIF-style metadata that describes an atom's
/// context within a biomolecular structure, including residue
/// identity, chain assignment, and terminal position.
///
/// # Fields
///
/// * `atom_name` — PDB atom name (e.g., "CA", "N", "O")
/// * `residue_name` — Three-letter residue code (e.g., "ALA", "HOH")
/// * `residue_id` — Residue sequence number
/// * `chain_id` — Single-character chain identifier
/// * `insertion_code` — PDB insertion code (usually ' ')
/// * `standard_name` — Parsed [`StandardResidue`] if recognized
/// * `category` — Residue classification
/// * `position` — Terminal position within chain
///
/// # Construction
///
/// Use [`AtomResidueInfo::builder`] for convenient construction
/// with default values for optional fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomResidueInfo {
    /// PDB-style atom name (e.g., "CA", "CB", "N").
    pub atom_name: String,
    /// Three-letter residue code (e.g., "ALA", "LIG").
    pub residue_name: String,
    /// Residue sequence number.
    pub residue_id: i32,
    /// Single-character chain identifier.
    pub chain_id: char,
    /// PDB insertion code for residue numbering conflicts.
    pub insertion_code: char,
    /// Parsed standard residue type, if recognized.
    pub standard_name: Option<StandardResidue>,
    /// Classification of the residue (standard, hetero, ion).
    pub category: ResidueCategory,
    /// Terminal position within the polymer chain.
    pub position: ResiduePosition,
}

impl AtomResidueInfo {
    /// Creates a builder for constructing an [`AtomResidueInfo`].
    ///
    /// # Arguments
    ///
    /// * `atom_name` — PDB-style atom name
    /// * `residue_name` — Three-letter residue code
    /// * `residue_id` — Residue sequence number
    /// * `chain_id` — Single-character chain identifier
    ///
    /// # Returns
    ///
    /// An [`AtomResidueBuilder`] with default values for optional fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_forge::{AtomResidueInfo, ResidueCategory, ResiduePosition};
    ///
    /// let info = AtomResidueInfo::builder("CA", "GLY", 1, 'A')
    ///     .category(ResidueCategory::Standard)
    ///     .position(ResiduePosition::Internal)
    ///     .build();
    ///
    /// assert_eq!(info.residue_name, "GLY");
    /// ```
    pub fn builder(
        atom_name: impl Into<String>,
        residue_name: impl Into<String>,
        residue_id: i32,
        chain_id: char,
    ) -> AtomResidueBuilder {
        AtomResidueBuilder {
            atom_name: atom_name.into(),
            residue_name: residue_name.into(),
            residue_id,
            chain_id,
            insertion_code: None,
            standard_name: None,
            category: ResidueCategory::Standard,
            position: ResiduePosition::None,
        }
    }
}

/// Builder for constructing [`AtomResidueInfo`] instances.
///
/// Provides a fluent API for setting optional fields with sensible defaults.
#[derive(Debug, Clone)]
pub struct AtomResidueBuilder {
    atom_name: String,
    residue_name: String,
    residue_id: i32,
    chain_id: char,
    insertion_code: Option<char>,
    standard_name: Option<StandardResidue>,
    category: ResidueCategory,
    position: ResiduePosition,
}

impl AtomResidueBuilder {
    /// Sets the PDB insertion code.
    ///
    /// # Arguments
    ///
    /// * `code` — Optional insertion code character
    pub fn insertion_code_opt(mut self, code: Option<char>) -> Self {
        self.insertion_code = code;
        self
    }

    /// Sets the standard residue type.
    ///
    /// # Arguments
    ///
    /// * `name` — Optional [`StandardResidue`] variant
    pub fn standard_name(mut self, name: Option<StandardResidue>) -> Self {
        self.standard_name = name;
        self
    }

    /// Sets the residue category.
    ///
    /// # Arguments
    ///
    /// * `category` — [`ResidueCategory`] classification
    pub fn category(mut self, category: ResidueCategory) -> Self {
        self.category = category;
        self
    }

    /// Sets the terminal position.
    ///
    /// # Arguments
    ///
    /// * `position` — [`ResiduePosition`] within the chain
    pub fn position(mut self, position: ResiduePosition) -> Self {
        self.position = position;
        self
    }

    /// Builds the final [`AtomResidueInfo`] instance.
    ///
    /// Applies default values for unset optional fields:
    /// * `insertion_code` defaults to `' '`
    ///
    /// # Returns
    ///
    /// A fully constructed [`AtomResidueInfo`].
    pub fn build(self) -> AtomResidueInfo {
        AtomResidueInfo {
            atom_name: self.atom_name,
            residue_name: self.residue_name,
            residue_id: self.residue_id,
            chain_id: self.chain_id,
            insertion_code: self.insertion_code.unwrap_or(' '),
            standard_name: self.standard_name,
            category: self.category,
            position: self.position,
        }
    }
}

/// Collection of biological metadata for all atoms in a system.
///
/// Stores per-atom biological annotations in a vector that parallels
/// the atom indices in a [`System`](crate::System). Each entry provides
/// PDB/mmCIF-style information about the atom's residue context.
///
/// # Examples
///
/// ```
/// use dreid_forge::{AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition};
///
/// let mut metadata = BioMetadata::with_capacity(100);
///
/// let info = AtomResidueInfo::builder("N", "ALA", 1, 'A')
///     .category(ResidueCategory::Standard)
///     .position(ResiduePosition::NTerminal)
///     .build();
///
/// metadata.atom_info.push(info);
/// assert_eq!(metadata.atom_info.len(), 1);
///
/// // Set target pH for charge assignment
/// metadata.target_ph = Some(7.0);
/// ```
#[derive(Debug, Clone, Default)]
pub struct BioMetadata {
    /// Per-atom biological annotations, indexed by atom order.
    pub atom_info: Vec<AtomResidueInfo>,
    /// Target pH for protonation state determination.
    pub target_ph: Option<f64>,
}

impl PartialEq for BioMetadata {
    fn eq(&self, other: &Self) -> bool {
        self.atom_info == other.atom_info && self.target_ph == other.target_ph
    }
}

impl Eq for BioMetadata {}

impl BioMetadata {
    /// Creates an empty [`BioMetadata`] container.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a [`BioMetadata`] with pre-allocated capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` — Number of atoms to pre-allocate space for
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            atom_info: Vec::with_capacity(capacity),
            target_ph: None,
        }
    }

    /// Returns the target pH, defaulting to physiological pH (7.4) if not set.
    pub fn effective_ph(&self) -> f64 {
        self.target_ph.unwrap_or(7.4)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_residue_info_new_and_all_fields() {
        let info = AtomResidueInfo::builder("CA", "ALA", 42, 'A')
            .insertion_code_opt(Some('x'))
            .standard_name(Some(StandardResidue::ALA))
            .category(ResidueCategory::Standard)
            .position(ResiduePosition::Internal)
            .build();
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
        let info = AtomResidueInfo::builder("N", "GLY", 1, 'B')
            .category(ResidueCategory::Hetero)
            .position(ResiduePosition::None)
            .build();
        assert_eq!(info.insertion_code, ' ');
        let cloned = info.clone();
        assert_eq!(info, cloned);
    }

    #[test]
    fn atom_residue_info_accepts_into_inputs() {
        let atom_name = String::from("O1");
        let residue_name = "LIG";
        let info = AtomResidueInfo::builder(atom_name, residue_name, 7, 'Z')
            .insertion_code_opt(Some('1'))
            .category(ResidueCategory::Hetero)
            .position(ResiduePosition::CTerminal)
            .build();
        assert_eq!(info.atom_name, "O1");
        assert_eq!(info.residue_name, "LIG");
        assert_eq!(info.residue_id, 7);
        assert_eq!(info.chain_id, 'Z');
        assert_eq!(info.insertion_code, '1');
        assert_eq!(info.standard_name, None);
        assert_eq!(info.category, ResidueCategory::Hetero);
        assert_eq!(info.position, ResiduePosition::CTerminal);
    }

    #[test]
    fn bio_metadata_new_and_capacity() {
        let mut bm = BioMetadata::with_capacity(4);
        assert!(bm.atom_info.capacity() >= 4);

        let info1 = AtomResidueInfo::builder("CA", "ALA", 1, 'A')
            .standard_name(Some(StandardResidue::ALA))
            .category(ResidueCategory::Standard)
            .position(ResiduePosition::Internal)
            .build();
        let info2 = AtomResidueInfo::builder("CB", "ALA", 1, 'A')
            .standard_name(Some(StandardResidue::ALA))
            .category(ResidueCategory::Standard)
            .position(ResiduePosition::Internal)
            .build();
        bm.atom_info.push(info1.clone());
        bm.atom_info.push(info2.clone());

        assert_eq!(bm.atom_info.len(), 2);
        assert_eq!(bm.atom_info[0], info1);
        assert_eq!(bm.atom_info[1], info2);
    }

    #[test]
    fn bio_metadata_target_ph() {
        let mut bm = BioMetadata::new();
        assert!(bm.target_ph.is_none());
        assert!((bm.effective_ph() - 7.4).abs() < f64::EPSILON);

        bm.target_ph = Some(6.5);
        assert!((bm.effective_ph() - 6.5).abs() < f64::EPSILON);
    }

    #[test]
    fn bio_metadata_equality() {
        let info = AtomResidueInfo::builder("N", "ALA", 1, 'A')
            .category(ResidueCategory::Standard)
            .build();

        let bm1 = BioMetadata {
            atom_info: vec![info.clone()],
            target_ph: Some(7.0),
        };
        let bm2 = BioMetadata {
            atom_info: vec![info.clone()],
            target_ph: Some(7.0),
        };
        let bm3 = BioMetadata {
            atom_info: vec![info],
            target_ph: Some(7.5),
        };

        assert_eq!(bm1, bm2);
        assert_ne!(bm1, bm3);
    }

    #[test]
    fn debug_contains_expected_fields() {
        let info = AtomResidueInfo::builder("C1", "LIG", -1, 'Z')
            .insertion_code_opt(Some('A'))
            .category(ResidueCategory::Hetero)
            .position(ResiduePosition::NTerminal)
            .build();
        let bm = BioMetadata {
            atom_info: vec![info.clone()],
            target_ph: None,
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
