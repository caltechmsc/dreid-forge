use std::path::PathBuf;

use dreid_forge::io::Format;

#[derive(Debug, Clone)]
pub struct OutputSpec {
    /// Path to write to, or `None` for stdout.
    pub path: Option<PathBuf>,
    /// Output format.
    pub format: Format,
}
