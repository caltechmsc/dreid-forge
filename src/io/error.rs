//! I/O error types for file operations.
//!
//! This module defines the error type used by all readers and writers
//! in the [`io`](super) module. Errors are categorized by source:
//! parse failures, unsupported formats, missing metadata, and underlying
//! I/O or library errors.

use super::Format;
use thiserror::Error;

/// Errors that can occur during molecular structure I/O operations.
///
/// This enum covers all failure modes when reading or writing molecular
/// structure files, including format parsing, compatibility issues,
/// and underlying I/O failures.
///
/// # Error Sources
///
/// * **I/O** — File system or stream errors
/// * **Parse** — Malformed or invalid file content
/// * **Unsupported** — Format not available for the requested operation
/// * **Metadata** — Missing required biological annotations
/// * **Conversion** — Data model translation failures
/// * **BioForge** — Errors from the bio-forge library
#[derive(Debug, Error)]
pub enum Error {
    /// Standard I/O operation failed.
    #[error("I/O operation failed: {source}")]
    Io {
        /// The underlying I/O error.
        #[from]
        source: std::io::Error,
    },

    /// Failed to parse file content.
    ///
    /// Includes format identification, approximate line number,
    /// and a description of the parsing issue.
    #[error("failed to parse {format} data: {details} (at line ~{line})")]
    Parse {
        /// The file format being parsed.
        format: Format,
        /// Approximate line number where the error occurred.
        line: usize,
        /// Description of the parsing problem.
        details: String,
    },

    /// The requested format does not support reading.
    #[error("the '{0}' format is not supported for this read operation")]
    UnsupportedReadFormat(Format),

    /// The requested format does not support writing.
    #[error("the '{0}' format is not supported for this write operation")]
    UnsupportedWriteFormat(Format),

    /// Missing biological metadata required for the output format.
    ///
    /// Formats like PDB, mmCIF, and BGF require residue and chain
    /// information that is stored in [`BioMetadata`](crate::BioMetadata).
    #[error("missing biological metadata required for writing to the {0} format")]
    MissingMetadata(&'static str),

    /// Biomolecular preparation step failed.
    ///
    /// Wraps errors from cleaning, protonation, solvation, or
    /// topology building operations.
    #[error("biomolecular preparation failed: {0}")]
    BioForgePreparation(String),

    /// Error from the bio-forge I/O layer.
    #[error("underlying bio-forge I/O error: {0}")]
    BioForgeIo(String),

    /// Data model conversion failed.
    ///
    /// Occurs when translating between dreid-forge and bio-forge
    /// data structures.
    #[error("failed to convert data model: {0}")]
    Conversion(String),
}

impl From<bio_forge::io::Error> for Error {
    fn from(e: bio_forge::io::Error) -> Self {
        Error::BioForgeIo(e.to_string())
    }
}

impl From<bio_forge::ops::Error> for Error {
    fn from(e: bio_forge::ops::Error) -> Self {
        Error::BioForgePreparation(e.to_string())
    }
}

impl From<super::util::ConversionError> for Error {
    fn from(e: super::util::ConversionError) -> Self {
        Error::Conversion(e.to_string())
    }
}

impl Error {
    /// Creates a new parse error with context.
    ///
    /// # Arguments
    ///
    /// * `format` — The file format being parsed
    /// * `line` — Approximate line number of the error
    /// * `details` — Description of the parsing problem
    ///
    /// # Returns
    ///
    /// A [`Parse`](Error::Parse) error variant.
    pub fn parse(format: Format, line: usize, details: impl Into<String>) -> Self {
        Self::Parse {
            format,
            line,
            details: details.into(),
        }
    }
}
