//! MPSim compatibility configuration.
//!
//! MPSim is a legacy EM/MM engine used within the group, fed by a
//! hand-maintained BGF conversion. Its DREIDING inputs follow three
//! conventions that a modern parameterization does not, and this adapter
//! reproduces them:
//!
//! 1. **Hydrogen-bond hydrogen naming** — MPSim expects the force-field type
//!    `H___A` where modern DREIDING assigns `H_HB` to polar (donor) hydrogens.
//! 2. **Chain ends** — every chain end is capped in its neutral, uncharged
//!    state regardless of the modeled pH, and a chain *break* is capped the
//!    same way as a chain end. The C side becomes an aldehyde (`–CHO`) rather
//!    than the carboxylic acid a real C-terminus carries.
//! 3. **Ligand labelling** — every ligand is relabelled to the single residue
//!    `RES 999`, whatever the input called it.
//!
//! None of the three is the better chemistry; each exists so output can be laid
//! beside those reference files atom for atom. They apply only when
//! [`MpsimConfig`] is attached to [`ForgeConfig`](super::ForgeConfig), and none
//! of them alters the underlying DREIDING parameterization.

/// Residue name every ligand takes under the MPSim convention.
pub const MPSIM_LIGAND_RESIDUE_NAME: &str = "RES";
/// Residue sequence number every ligand takes under the MPSim convention.
pub const MPSIM_LIGAND_RESIDUE_ID: i32 = 999;

/// Compatibility adapter for the legacy MPSim EM/MM engine.
///
/// When attached to [`ForgeConfig`](super::ForgeConfig) via
/// [`mpsim`](super::ForgeConfig::mpsim), the [`forge`](crate::forge) pipeline
/// applies MPSim's DREIDING conventions. All behaviors are enabled by default
/// and can be toggled independently.
///
/// # Examples
///
/// ```
/// use dreid_forge::{ForgeConfig, MpsimConfig};
///
/// // Enable the full MPSim adapter with default behavior.
/// let config = ForgeConfig {
///     mpsim: Some(MpsimConfig::default()),
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct MpsimConfig {
    /// Rename hydrogen-bond hydrogens from DREIDING `H_HB` to MPSim `H___A`.
    ///
    /// The rename is applied only to the emitted force-field type names; it
    /// does not affect parameter lookup, hydrogen-bond term generation, or
    /// atom typing internally.
    pub rename_hb_hydrogen: bool,

    /// Cap every chain end, and every chain break, in its neutral state.
    ///
    /// Independent of pH, and applied identically to a real chain end and to a
    /// break where the model simply stops:
    ///
    /// | side   | DREID-Forge      | MPSim                     |
    /// | ------ | ---------------- | ------------------------- |
    /// | N      | `–NH₃⁺` / `–NH`  | `–NH₂`                    |
    /// | C      | `–COO⁻` / `–CO`  | `–CHO` (cap hydrogen `HC`) |
    ///
    /// The C side is the notable one: DREID-Forge builds the carboxylic acid a
    /// real C-terminus carries, while MPSim's converter caps the backbone
    /// carbon directly, leaving an aldehyde. Matching charge sets are applied so
    /// every capped residue stays net-neutral.
    ///
    /// Nucleic-acid 5′/3′ termini are untouched.
    pub neutral_termini: bool,

    /// Relabel every ligand to [`MPSIM_LIGAND_RESIDUE_NAME`]
    /// [`MPSIM_LIGAND_RESIDUE_ID`].
    ///
    /// DREID-Forge otherwise keeps whatever the input called the ligand, since
    /// a structure may hold several distinct ligands that one shared label
    /// collapses into a single residue.
    pub relabel_ligand: bool,
}

impl Default for MpsimConfig {
    fn default() -> Self {
        Self {
            rename_hb_hydrogen: true,
            neutral_termini: true,
            relabel_ligand: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_enables_all_behaviors() {
        let config = MpsimConfig::default();
        assert!(config.rename_hb_hydrogen);
        assert!(config.neutral_termini);
        assert!(config.relabel_ligand);
    }
}
