//! Charge calculation method configurations.
//!
//! This module defines the charge assignment strategies available in the
//! DREIDING parameterization pipeline.

pub use ffcharge::NucleicScheme;
pub use ffcharge::ProteinScheme;
pub use ffcharge::WaterScheme;

pub use cheq::{BasisType, DampingStrategy, SolverOptions};

/// Method for calculating partial atomic charges.
///
/// Determines how partial charges are assigned to atoms during
/// parameterization. The choice affects both the accuracy of
/// electrostatic interactions and hydrogen bond parameters.
#[derive(Debug, Clone, Default)]
pub enum ChargeMethod {
    /// No charge calculation; all charges remain zero.
    ///
    /// Use this for gas-phase calculations where electrostatics
    /// are not critical, or when charges will be assigned externally.
    #[default]
    None,

    /// Charge equilibration (QEq) method for all atoms.
    ///
    /// Calculates electronegativity-equalized charges based on
    /// atomic positions and electronegativity parameters. Best suited
    /// for small molecules or systems without biological metadata.
    Qeq(QeqConfig),

    /// Hybrid charge assignment for biological systems.
    ///
    /// Combines classical force field charges for biomolecules with
    /// QEq for ligands and hetero groups. Requires biological metadata
    /// to be present in the input system.
    Hybrid(HybridConfig),
}

/// Configuration for QEq charge equilibration.
///
/// Controls the behavior of the charge equilibration solver,
/// including the target total charge and numerical solver options.
///
/// # Examples
///
/// ```
/// use dreid_forge::QeqConfig;
///
/// // Neutral molecule (default)
/// let neutral = QeqConfig::default();
/// assert_eq!(neutral.total_charge, 0.0);
///
/// // Negatively charged system
/// let anion = QeqConfig {
///     total_charge: -1.0,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct QeqConfig {
    /// Target total charge of the system in elementary charge units.
    ///
    /// The QEq solver will constrain the sum of all partial charges
    /// to equal this value. Default is `0.0` (neutral).
    pub total_charge: f64,

    /// QEq solver options.
    ///
    /// Controls convergence tolerance, iteration limits, and
    /// hydrogen SCF treatment.
    pub solver_options: SolverOptions,
}

impl Default for QeqConfig {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            solver_options: SolverOptions::default(),
        }
    }
}

/// Configuration for hybrid biological/QEq charge assignment.
///
/// Specifies how charges are assigned to different molecule types
/// in a biological system. Standard residues (proteins, nucleic acids,
/// water, ions) receive classical force field charges, while hetero
/// groups (ligands) can use vacuum or embedded QEq.
///
/// # Molecule Classification
///
/// The hybrid method classifies atoms based on [`ResidueCategory`](crate::ResidueCategory):
///
/// | Category | Charge Source |
/// |----------|---------------|
/// | Standard amino acid | [`protein_scheme`](Self::protein_scheme) |
/// | Standard nucleotide | [`nucleic_scheme`](Self::nucleic_scheme) |
/// | Water (HOH) | [`water_scheme`](Self::water_scheme) |
/// | Ion | Formal charges (classical, integral) |
/// | Hetero | Per-ligand QEq via [`ligand_configs`](Self::ligand_configs) |
#[derive(Debug, Clone, Default)]
pub struct HybridConfig {
    /// Force field scheme for protein residue charges.
    ///
    /// Default is [`ProteinScheme::AmberFFSB`] (AMBER ff99SB/ff14SB/ff19SB).
    pub protein_scheme: ProteinScheme,

    /// Force field scheme for nucleic acid residue charges.
    ///
    /// Default is [`NucleicScheme::Amber`] (AMBER OL15/OL21/OL24/bsc1/OL3).
    pub nucleic_scheme: NucleicScheme,

    /// Water model for solvent charges.
    ///
    /// Default is [`WaterScheme::Tip3p`].
    pub water_scheme: WaterScheme,

    /// Per-ligand charge configuration.
    ///
    /// Each entry specifies a residue selector and the QEq method to use.
    /// Ligands not explicitly listed will use the [`default_ligand_qeq`](Self::default_ligand_qeq).
    pub ligand_configs: Vec<LigandChargeConfig>,

    /// Default QEq configuration for ligands not in [`ligand_configs`](Self::ligand_configs).
    ///
    /// Default is vacuum QEq with neutral total charge.
    pub default_ligand_qeq: QeqConfig,
}

/// Residue selector for identifying specific residues.
///
/// Used to target specific ligands or hetero groups for custom
/// charge assignment in hybrid mode.
///
/// # Examples
///
/// ```
/// use dreid_forge::ResidueSelector;
///
/// // Select residue 500 in chain A
/// let selector = ResidueSelector::new("A", 500, None);
///
/// // Select residue 100 with insertion code 'B' in chain L
/// let with_icode = ResidueSelector::new("L", 100, Some('B'));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueSelector {
    /// Chain identifier.
    pub chain_id: String,
    /// Residue sequence number.
    pub residue_id: i32,
    /// Optional insertion code for distinguishing residues with same ID.
    pub insertion_code: Option<char>,
}

impl ResidueSelector {
    /// Creates a new residue selector.
    ///
    /// # Arguments
    ///
    /// * `chain_id` — Chain identifier
    /// * `residue_id` — Residue sequence number
    /// * `insertion_code` — Optional insertion code
    pub fn new(chain_id: impl Into<String>, residue_id: i32, insertion_code: Option<char>) -> Self {
        Self {
            chain_id: chain_id.into(),
            residue_id,
            insertion_code,
        }
    }

    /// Checks if this selector matches an atom's residue info.
    ///
    /// # Arguments
    ///
    /// * `chain_id` — Chain ID to match
    /// * `residue_id` — Residue ID to match
    /// * `insertion_code` — Insertion code to match
    ///
    /// # Returns
    ///
    /// `true` if the selector matches the given residue info,
    /// `false` otherwise.
    pub fn matches(&self, chain_id: &str, residue_id: i32, insertion_code: Option<char>) -> bool {
        self.chain_id == chain_id
            && self.residue_id == residue_id
            && (self.insertion_code.is_none() || self.insertion_code == insertion_code)
    }
}

/// Charge configuration for a specific ligand residue.
///
/// Combines a residue selector with the QEq method to use for
/// that specific ligand.
#[derive(Debug, Clone)]
pub struct LigandChargeConfig {
    /// Selector identifying the target residue.
    pub selector: ResidueSelector,
    /// QEq method to use for this ligand.
    pub method: LigandQeqMethod,
}

/// QEq method variant for ligand charge assignment.
///
/// Ligands can use either vacuum QEq (isolated) or embedded QEq
/// (polarized by surrounding fixed charges).
#[derive(Debug, Clone)]
pub enum LigandQeqMethod {
    /// Vacuum QEq: ligand treated as isolated molecule.
    ///
    /// Best for ligands in solution or far from biomolecular surfaces.
    Vacuum(QeqConfig),

    /// Embedded QEq: ligand polarized by surrounding fixed charges.
    ///
    /// The ligand's charge distribution is influenced by nearby
    /// protein/nucleic acid atoms within the cutoff radius.
    Embedded(EmbeddedQeqConfig),
}

impl Default for LigandQeqMethod {
    fn default() -> Self {
        Self::Vacuum(QeqConfig::default())
    }
}

/// Configuration for embedded QEq calculations.
///
/// In embedded QEq, atoms from surrounding biomolecules (with fixed charges)
/// contribute an external electrostatic potential that polarizes the
/// ligand's charge distribution.
#[derive(Debug, Clone)]
pub struct EmbeddedQeqConfig {
    /// Cutoff radius for including environment atoms (Ångströms).
    ///
    /// Atoms within this distance from any ligand atom contribute
    /// to the external electrostatic potential. Default is `10.0` Å.
    pub cutoff_radius: f64,

    /// QEq solver configuration for the ligand.
    pub qeq: QeqConfig,
}

impl Default for EmbeddedQeqConfig {
    fn default() -> Self {
        Self {
            cutoff_radius: 10.0,
            qeq: QeqConfig::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn charge_method_default_is_none() {
        assert!(matches!(ChargeMethod::default(), ChargeMethod::None));
    }

    #[test]
    fn qeq_config_default_values() {
        let config = QeqConfig::default();
        assert_eq!(config.total_charge, 0.0);
        assert!(config.solver_options.hydrogen_scf);
    }

    #[test]
    fn hybrid_config_default_schemes() {
        let config = HybridConfig::default();
        assert_eq!(config.protein_scheme, ProteinScheme::AmberFFSB);
        assert_eq!(config.nucleic_scheme, NucleicScheme::Amber);
        assert_eq!(config.water_scheme, WaterScheme::Tip3p);
        assert!(config.ligand_configs.is_empty());
    }

    #[test]
    fn residue_selector_matches() {
        let selector = ResidueSelector::new("A", 100, None);
        assert!(selector.matches("A", 100, None));
        assert!(selector.matches("A", 100, Some('B')));
        assert!(!selector.matches("B", 100, None));
        assert!(!selector.matches("A", 101, None));

        let with_icode = ResidueSelector::new("A", 100, Some('X'));
        assert!(with_icode.matches("A", 100, Some('X')));
        assert!(!with_icode.matches("A", 100, None));
        assert!(!with_icode.matches("A", 100, Some('Y')));
    }

    #[test]
    fn ligand_qeq_method_default_is_vacuum() {
        assert!(matches!(
            LigandQeqMethod::default(),
            LigandQeqMethod::Vacuum(_)
        ));
    }

    #[test]
    fn embedded_qeq_config_default_cutoff() {
        let config = EmbeddedQeqConfig::default();
        assert_eq!(config.cutoff_radius, 10.0);
    }

    #[test]
    fn hybrid_config_with_custom_ligand() {
        let config = HybridConfig {
            ligand_configs: vec![LigandChargeConfig {
                selector: ResidueSelector::new("A", 500, None),
                method: LigandQeqMethod::Embedded(EmbeddedQeqConfig {
                    cutoff_radius: 8.0,
                    qeq: QeqConfig {
                        total_charge: -1.0,
                        ..Default::default()
                    },
                }),
            }],
            ..Default::default()
        };

        assert_eq!(config.ligand_configs.len(), 1);
        assert!(config.ligand_configs[0].selector.matches("A", 500, None));

        if let LigandQeqMethod::Embedded(embedded) = &config.ligand_configs[0].method {
            assert_eq!(embedded.cutoff_radius, 8.0);
            assert_eq!(embedded.qeq.total_charge, -1.0);
        } else {
            panic!("Expected Embedded variant");
        }
    }
}
