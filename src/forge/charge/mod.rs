//! Partial charge calculation for molecular systems.
//!
//! This module handles the assignment of partial atomic charges using
//! multiple methods, including global QEq, hybrid (biological + QEq), and
//! no charges.

mod hybrid;
mod qeq;
mod spatial;

use super::config::ChargeMethod;
use super::error::Error;
use super::intermediate::IntermediateSystem;

/// Assigns partial charges to atoms based on the configured method.
///
/// Modifies the `charge` field of each atom in the intermediate system
/// according to the specified charge method.
///
/// # Arguments
///
/// * `system` — Mutable reference to the intermediate system
/// * `method` — Charge calculation method to use
///
/// # Errors
///
/// Returns [`Error`] if:
/// - QEq solver fails to converge ([`Error::ChargeCalculation`])
/// - Hybrid method is used without biological metadata ([`Error::MissingBioMetadata`])
/// - Classical charge lookup fails ([`Error::HybridChargeAssignment`])
pub fn assign_charges(system: &mut IntermediateSystem, method: &ChargeMethod) -> Result<(), Error> {
    match method {
        ChargeMethod::None => Ok(()),
        ChargeMethod::Qeq(config) => qeq::assign_qeq_charges(system, config),
        ChargeMethod::Hybrid(config) => hybrid::assign_hybrid_charges(system, config),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forge::config::QeqConfig;
    use crate::model::atom::Atom;
    use crate::model::system::{Bond, System};
    use crate::model::types::{BondOrder, Element};

    fn make_water() -> System {
        let mut sys = System::new();
        let bond_length = 0.9575;
        let angle = 104.45f64.to_radians();

        sys.atoms.push(Atom::new(Element::O, [0.0, 0.0, 0.0]));
        sys.atoms
            .push(Atom::new(Element::H, [bond_length, 0.0, 0.0]));
        sys.atoms.push(Atom::new(
            Element::H,
            [bond_length * angle.cos(), bond_length * angle.sin(), 0.0],
        ));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys
    }

    #[test]
    fn no_charge_method() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();

        assign_charges(&mut int, &ChargeMethod::None).unwrap();

        for atom in &int.atoms {
            assert_eq!(atom.charge, 0.0);
        }
    }

    #[test]
    fn qeq_charge_assignment() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();

        let qeq_config = QeqConfig::default();
        assign_charges(&mut int, &ChargeMethod::Qeq(qeq_config)).unwrap();

        assert!(int.atoms[0].charge < 0.0);
        assert!(int.atoms[1].charge > 0.0);
        assert!(int.atoms[2].charge > 0.0);

        let total: f64 = int.atoms.iter().map(|a| a.charge).sum();
        assert!((total).abs() < 1e-9);
    }

    #[test]
    fn qeq_with_net_charge() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();

        let qeq_config = QeqConfig {
            total_charge: -1.0,
            ..Default::default()
        };
        assign_charges(&mut int, &ChargeMethod::Qeq(qeq_config)).unwrap();

        let total: f64 = int.atoms.iter().map(|a| a.charge).sum();
        assert!((total + 1.0).abs() < 1e-9);
    }
}
