use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("failed to parse force field parameters: {0}")]
    ParameterParse(#[from] toml::de::Error),

    #[error("failed to parse custom typing rules: {0}")]
    RuleParse(String),

    #[error("atom typing failed: {0}")]
    AtomTyping(String),

    #[error("charge calculation failed: {0}")]
    ChargeCalculation(String),

    #[error("missing force field parameter for atom type '{atom_type}': {detail}")]
    MissingParameter { atom_type: String, detail: String },

    #[error("invalid bond between atoms {i} and {j}: {detail}")]
    InvalidBond { i: usize, j: usize, detail: String },

    #[error("input system is empty: at least one atom is required")]
    EmptySystem,

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
    pub fn missing_parameter(atom_type: &str, details: impl Into<String>) -> Self {
        Self::MissingParameter {
            atom_type: atom_type.to_string(),
            detail: details.into(),
        }
    }

    pub fn invalid_bond(i: usize, j: usize, details: impl Into<String>) -> Self {
        Self::InvalidBond {
            i,
            j,
            detail: details.into(),
        }
    }
}
