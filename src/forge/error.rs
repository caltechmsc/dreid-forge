//! Error types for DREIDING parameterization.
//!
//! This module defines the error type used throughout the forge module.
//! Errors are categorized by source: parameter parsing, atom typing,
//! charge calculation, and missing force field parameters.

use thiserror::Error;

/// Errors that can occur during DREIDING parameterization.
///
/// This enum covers all failure modes when running the [`forge`](super::forge)
/// function, including configuration parsing, atom typing, charge calculation,
/// and parameter lookup failures.
#[derive(Debug, Error)]
pub enum Error {
    /// Failed to parse force field parameters TOML.
    #[error("failed to parse force field parameters: {0}")]
    ParameterParse(#[from] toml::de::Error),

    /// Failed to parse custom atom typing rules.
    #[error("failed to parse custom typing rules: {0}")]
    RuleParse(String),

    /// Atom typing assignment failed.
    ///
    /// Occurs when the molecular graph cannot be typed due to unsupported
    /// elements, unusual bonding patterns, or rule matching failures.
    #[error("atom typing failed: {0}")]
    AtomTyping(String),

    /// Charge calculation failed.
    ///
    /// Occurs when the QEq solver fails to converge or encounters
    /// invalid input data.
    #[error("charge calculation failed: {0}")]
    ChargeCalculation(String),

    /// Hybrid charge method requires biological metadata.
    ///
    /// Occurs when [`ChargeMethod::Hybrid`](crate::ChargeMethod::Hybrid) is selected
    /// but the input system has no [`BioMetadata`](crate::BioMetadata).
    #[error("hybrid charge method requires biological metadata (bio_metadata is None)")]
    MissingBioMetadata,

    /// Hybrid charge assignment failed for a specific residue.
    #[error(
        "hybrid charge assignment failed for residue {residue_name} at chain {chain_id} residue {residue_id}: {detail}"
    )]
    HybridChargeAssignment {
        /// Chain identifier.
        chain_id: String,
        /// Residue sequence number.
        residue_id: i32,
        /// Residue name.
        residue_name: String,
        /// Description of the problem.
        detail: String,
    },

    /// Required force field parameter not found.
    ///
    /// Occurs when an atom type is assigned but no corresponding
    /// parameters exist in the force field parameter file.
    #[error("missing force field parameter for atom type '{atom_type}': {detail}")]
    MissingParameter {
        /// The atom type that is missing parameters.
        atom_type: String,
        /// Description of which parameter is missing.
        detail: String,
    },

    /// Invalid bond definition in the input system.
    #[error("invalid bond between atoms {i} and {j}: {detail}")]
    InvalidBond {
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
        /// Description of the problem.
        detail: String,
    },

    /// The input system contains no atoms.
    #[error("input system is empty: at least one atom is required")]
    EmptySystem,

    /// Internal data conversion error.
    #[error("internal conversion error: {0}")]
    Conversion(String),
}

impl From<dreid_typer::TyperError> for Error {
    fn from(e: dreid_typer::TyperError) -> Self {
        Error::AtomTyping(e.to_string())
    }
}

impl From<cheq::CheqError> for Error {
    fn from(e: cheq::CheqError) -> Self {
        Error::ChargeCalculation(e.to_string())
    }
}

impl Error {
    /// Creates a [`HybridChargeAssignment`](Error::HybridChargeAssignment) error.
    ///
    /// # Arguments
    ///
    /// * `chain_id` — Chain identifier
    /// * `residue_id` — Residue sequence number
    /// * `residue_name` — Residue name
    /// * `details` — Description of the problem
    ///
    /// # Returns
    ///
    /// A [`HybridChargeAssignment`](Error::HybridChargeAssignment) error variant.
    pub fn hybrid_charge_assignment(
        chain_id: impl Into<String>,
        residue_id: i32,
        residue_name: &str,
        details: impl Into<String>,
    ) -> Self {
        Self::HybridChargeAssignment {
            chain_id: chain_id.into(),
            residue_id,
            residue_name: residue_name.to_string(),
            detail: details.into(),
        }
    }

    /// Creates a [`MissingParameter`](Error::MissingParameter) error.
    ///
    /// # Arguments
    ///
    /// * `atom_type` — The atom type that is missing parameters
    /// * `details` — Description of which parameter is missing
    ///
    /// # Returns
    ///
    /// A [`MissingParameter`](Error::MissingParameter) error variant.
    pub fn missing_parameter(atom_type: &str, details: impl Into<String>) -> Self {
        Self::MissingParameter {
            atom_type: atom_type.to_string(),
            detail: details.into(),
        }
    }

    /// Creates an [`InvalidBond`](Error::InvalidBond) error.
    ///
    /// # Arguments
    ///
    /// * `i` — First atom index
    /// * `j` — Second atom index
    /// * `details` — Description of the bond problem
    ///
    /// # Returns
    ///
    /// An [`InvalidBond`](Error::InvalidBond) error variant.
    pub fn invalid_bond(i: usize, j: usize, details: impl Into<String>) -> Self {
        Self::InvalidBond {
            i,
            j,
            detail: details.into(),
        }
    }
}
