use std::path::{Path, PathBuf};

pub fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let stem = path.file_stem().unwrap_or_default();
    path.with_file_name(format!("{}{}", stem.to_string_lossy(), suffix))
}
