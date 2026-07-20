//! MPSim compatibility configuration.
//!
//! MPSim is a legacy EM/MM engine used within the group. Its DREIDING
//! implementation differs from the modern convention in two ways that this
//! adapter reconciles:
//!
//! 1. **Hydrogen-bond hydrogen naming** — MPSim expects the force-field type
//!    `H___A` where modern DREIDING assigns `H_HB` to polar (donor) hydrogens.
//! 2. **Protein chain termini** — MPSim inputs use the neutral, uncharged
//!    protonation state for both the N-terminus (–NH₂) and the C-terminus
//!    (–COOH), regardless of the modeled pH.
//!
//! Enabling [`MpsimConfig`] on [`ForgeConfig`](super::ForgeConfig) applies
//! these conventions as a self-contained post-/pre-processing layer without
//! altering the underlying DREIDING parameterization.

/// Compatibility adapter for the legacy MPSim EM/MM engine.
///
/// When attached to [`ForgeConfig`](super::ForgeConfig) via
/// [`mpsim`](super::ForgeConfig::mpsim), the [`forge`](crate::forge) pipeline
/// applies MPSim's DREIDING conventions. Both behaviors are enabled by
/// default and can be toggled independently.
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

    /// Force neutral, uncharged N/C protein termini.
    ///
    /// When `true`, the C-terminus is normalized to the protonated carboxyl
    /// (`–COOH`, adding `HOXT`) and the N-terminus to the neutral amine
    /// (`–NH₂`, removing the extra `H3`), independent of pH. The hybrid
    /// charge assignment likewise uses the neutral terminal charge sets so
    /// each terminal residue carries no net formal charge.
    pub neutral_termini: bool,
}

impl Default for MpsimConfig {
    fn default() -> Self {
        Self {
            rename_hb_hydrogen: true,
            neutral_termini: true,
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
    }
}
