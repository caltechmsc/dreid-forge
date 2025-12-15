mod banner;
mod error;
mod progress;
mod tables;

pub use banner::{banner_for_help, print_banner};
pub use error::print_error;
pub use progress::Progress;
pub use tables::{print_atom_types, print_chain_breakdown, print_parameters, print_structure_info};

#[derive(Debug, Clone, Copy)]
pub struct Context {
    pub interactive: bool,
}

impl Context {
    pub fn detect() -> Self {
        Self {
            interactive: crate::io::stderr_is_tty(),
        }
    }

    pub fn with_quiet(self, quiet: bool) -> Self {
        if quiet {
            Self { interactive: false }
        } else {
            self
        }
    }
}
