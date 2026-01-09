//! DREIDING force field parameterization engine.
//!
//! This module provides the core functionality for assigning DREIDING force field
//! parameters to molecular systems. The parameterization pipeline includes atom
//! typing, partial charge calculation, and generation of all bonded/non-bonded
//! potential energy function parameters.

mod charge;
mod config;
mod error;
mod intermediate;
mod paramgen;
mod params;
mod typer;

pub use config::{
    AnglePotentialType, BasisType, BondPotentialType, ChargeMethod, DampingStrategy,
    EmbeddedQeqConfig, ForgeConfig, HybridConfig, LigandChargeConfig, LigandQeqMethod,
    NucleicScheme, ProteinScheme, QeqConfig, ResidueSelector, SolverOptions, VdwPotentialType,
    WaterScheme,
};
pub use error::Error;

use crate::model::system::System;
use crate::model::topology::ForgedSystem;

/// Parameterizes a molecular system with DREIDING force field parameters.
///
/// This is the main entry point for force field parameterization. It takes
/// a molecular [`System`] with atoms and bonds, and produces a [`ForgedSystem`]
/// containing all parameters needed for molecular dynamics simulation.
///
/// # Arguments
///
/// * `system` — The molecular system to parameterize (must have at least one atom)
/// * `config` — Configuration controlling potential types and charge method
///
/// # Returns
///
/// A [`ForgedSystem`] containing:
/// - Original system structure
/// - Assigned atom type names
/// - Per-atom charges and masses
/// - All potential energy function parameters
///
/// # Errors
///
/// Returns [`Error`] if:
/// - The system is empty (no atoms)
/// - Atom typing fails (unsupported element or bonding pattern)
/// - Force field parameters are missing for an atom type
/// - Charge calculation fails
///
/// # Examples
///
/// ```
/// use dreid_forge::{System, Atom, Bond, Element, BondOrder};
/// use dreid_forge::{forge, ForgeConfig, ForgeError};
///
/// let mut system = System::new();
/// system.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
/// system.atoms.push(Atom::new(Element::C, [1.54, 0.0, 0.0]));
/// system.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let forged = forge(&system, &ForgeConfig::default())?;
/// assert!(!forged.potentials.bonds.is_empty());
/// # Ok::<(), ForgeError>(())
/// ```
pub fn forge(system: &System, config: &ForgeConfig) -> Result<ForgedSystem, Error> {
    let ff_params = params::load_parameters(config.params.as_deref())?;

    let mut intermediate = intermediate::IntermediateSystem::from_system(system)?;

    typer::assign_atom_types(&mut intermediate, config.rules.as_deref())?;

    charge::assign_charges(&mut intermediate, &config.charge_method)?;

    let forged = paramgen::generate_parameters(system, &intermediate, &ff_params, config)?;

    Ok(forged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::system::Bond;
    use crate::model::types::{BondOrder, Element};

    fn make_water() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::O, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [0.9575, 0.0, 0.0]));
        sys.atoms
            .push(Atom::new(Element::H, [-0.2399, 0.9272, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys
    }

    fn make_ethane() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::C, [1.54, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, 0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, 0.89]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 3, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 4, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 5, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 6, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 7, BondOrder::Single));
        sys
    }

    fn make_benzene() -> System {
        let mut sys = System::new();
        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            sys.atoms.push(Atom::new(
                Element::C,
                [1.4 * angle.cos(), 1.4 * angle.sin(), 0.0],
            ));
        }
        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            sys.atoms.push(Atom::new(
                Element::H,
                [2.5 * angle.cos(), 2.5 * angle.sin(), 0.0],
            ));
        }
        for i in 0..6 {
            sys.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
            sys.bonds.push(Bond::new(i, i + 6, BondOrder::Single));
        }
        sys
    }

    #[test]
    fn forges_water_with_default_config() {
        let water = make_water();
        let config = ForgeConfig::default();
        let forged = forge(&water, &config).expect("should forge water");

        assert_eq!(forged.atom_properties.len(), 3);
        assert_eq!(forged.atom_types.len(), 2);
        assert_eq!(forged.potentials.bonds.len(), 2);
        assert_eq!(forged.potentials.angles.len(), 1);
        assert!(forged.potentials.dihedrals.is_empty());
    }

    #[test]
    fn forges_ethane_with_dihedrals() {
        let ethane = make_ethane();
        let config = ForgeConfig::default();
        let forged = forge(&ethane, &config).expect("should forge ethane");

        assert_eq!(forged.atom_properties.len(), 8);
        assert_eq!(forged.potentials.bonds.len(), 7);
        assert!(!forged.potentials.angles.is_empty());
        assert!(!forged.potentials.dihedrals.is_empty());
    }

    #[test]
    fn forges_benzene_with_impropers() {
        let benzene = make_benzene();
        let config = ForgeConfig::default();
        let forged = forge(&benzene, &config).expect("should forge benzene");

        assert_eq!(forged.atom_properties.len(), 12);
        assert!(forged.atom_types.contains(&"C_R".to_string()));
        assert!(!forged.potentials.impropers.is_empty());
    }

    #[test]
    fn forges_with_qeq_charges() {
        let water = make_water();
        let config = ForgeConfig {
            charge_method: ChargeMethod::Qeq(Default::default()),
            ..Default::default()
        };
        let forged = forge(&water, &config).expect("should forge with QEq");

        let o_charge = forged.atom_properties[0].charge;
        let h1_charge = forged.atom_properties[1].charge;
        let h2_charge = forged.atom_properties[2].charge;

        assert!(o_charge < 0.0, "oxygen should be negative");
        assert!(h1_charge > 0.0, "hydrogen should be positive");
        assert!(
            (o_charge + h1_charge + h2_charge).abs() < 1e-9,
            "total charge should be zero"
        );
    }

    #[test]
    fn forges_with_morse_bonds() {
        let water = make_water();
        let config = ForgeConfig {
            bond_potential: BondPotentialType::Morse,
            ..Default::default()
        };
        let forged = forge(&water, &config).expect("should forge with Morse");

        for bond in &forged.potentials.bonds {
            assert!(matches!(
                bond,
                crate::model::topology::BondPotential::Morse { .. }
            ));
        }
    }

    #[test]
    fn forges_with_cosine_harmonic_angles() {
        let water = make_water();
        let config = ForgeConfig {
            angle_potential: AnglePotentialType::CosineHarmonic,
            ..Default::default()
        };
        let forged = forge(&water, &config).expect("should forge with CosineHarmonic");

        for angle in &forged.potentials.angles {
            assert!(matches!(
                angle,
                crate::model::topology::AnglePotential::CosineHarmonic { .. }
            ));
        }
    }

    #[test]
    fn forges_with_exp6_vdw() {
        let water = make_water();
        let config = ForgeConfig {
            vdw_potential: VdwPotentialType::Exponential6,
            ..Default::default()
        };
        let forged = forge(&water, &config).expect("should forge with Exp6");

        for vdw in &forged.potentials.vdw_pairs {
            assert!(matches!(
                vdw,
                crate::model::topology::VdwPairPotential::Exponential6 { .. }
            ));
        }
    }

    #[test]
    fn errors_on_empty_system() {
        let empty = System::new();
        let config = ForgeConfig::default();
        let result = forge(&empty, &config);

        assert!(matches!(result, Err(Error::EmptySystem)));
    }

    #[test]
    fn errors_on_invalid_custom_params() {
        let water = make_water();
        let config = ForgeConfig {
            params: Some("invalid [[[ toml".to_string()),
            ..Default::default()
        };
        let result = forge(&water, &config);

        assert!(matches!(result, Err(Error::ParameterParse(_))));
    }

    #[test]
    fn preserves_original_system_in_output() {
        let water = make_water();
        let config = ForgeConfig::default();
        let forged = forge(&water, &config).expect("should forge");

        assert_eq!(forged.system.atoms.len(), water.atoms.len());
        assert_eq!(forged.system.bonds.len(), water.bonds.len());
    }

    #[test]
    fn generates_vdw_pairs_for_all_type_combinations() {
        let water = make_water();
        let config = ForgeConfig::default();
        let forged = forge(&water, &config).expect("should forge");

        let n_types = forged.atom_types.len();
        let expected_pairs = n_types * (n_types + 1) / 2;
        assert_eq!(forged.potentials.vdw_pairs.len(), expected_pairs);
    }
}
