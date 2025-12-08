use super::Format;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("I/O operation failed: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },

    #[error("failed to parse {format} data: {details} (at line ~{line})")]
    Parse {
        format: Format,
        line: usize,
        details: String,
    },

    #[error("the '{0}' format is not supported for this read operation")]
    UnsupportedReadFormat(Format),

    #[error("the '{0}' format is not supported for this write operation")]
    UnsupportedWriteFormat(Format),

    #[error("missing biological metadata required for writing to the {0} format")]
    MissingMetadata(&'static str),

    #[error("biomolecular preparation failed: {0}")]
    BioForgePreparation(String),

    #[error("underlying bio-forge I/O error: {0}")]
    BioForgeIo(String),

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
    pub fn parse(format: Format, line: usize, details: impl Into<String>) -> Self {
        Self::Parse {
            format,
            line,
            details: details.into(),
        }
    }
}
