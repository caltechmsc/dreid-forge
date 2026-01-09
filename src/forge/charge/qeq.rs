//! QEq (Charge Equilibration) charge assignment.
//!
//! This module implements global QEq charge assignment for molecular systems.
//! It uses the `cheq` library to compute electronegativity-equalized partial
//! charges based on atomic positions.

use crate::forge::config::QeqConfig;
use crate::forge::error::Error;
use crate::forge::intermediate::IntermediateSystem;
use cheq::{QEqSolver, get_default_parameters};

/// Assigns charges using the QEq (charge equilibration) method.
///
/// Uses the `cheq` library to compute electronegativity-equalized
/// partial charges based on atomic positions.
///
/// # Arguments
///
/// * `system` — Mutable reference to the intermediate system
/// * `config` — QEq configuration (total charge, solver options)
///
/// # Errors
///
/// Returns [`Error::ChargeCalculation`] if the QEq solver fails to converge.
pub fn assign_qeq_charges(
    system: &mut IntermediateSystem,
    config: &QeqConfig,
) -> Result<(), Error> {
    let params = get_default_parameters();
    let solver = QEqSolver::new(params).with_options(config.solver_options);

    let atoms = system.atoms();
    let result = solver.solve(atoms, config.total_charge)?;

    for (atom, &charge) in system.atoms.iter_mut().zip(result.charges.iter()) {
        atom.charge = charge;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::system::{Bond, System};
    use crate::model::types::{BondOrder, Element};

    fn make_methane() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [1.09, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.363, 1.028, 0.0]));
        sys.atoms
            .push(Atom::new(Element::H, [-0.363, -0.514, 0.890]));
        sys.atoms
            .push(Atom::new(Element::H, [-0.363, -0.514, -0.890]));

        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 3, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 4, BondOrder::Single));
        sys
    }

    #[test]
    fn qeq_methane_neutral() {
        let methane = make_methane();
        let mut int = IntermediateSystem::from_system(&methane).unwrap();

        assign_qeq_charges(&mut int, &QeqConfig::default()).unwrap();

        let c_charge = int.atoms[0].charge;
        let h_charges: Vec<f64> = int.atoms[1..].iter().map(|a| a.charge).collect();

        assert!(c_charge < 0.0, "carbon should be negative in CH4");
        for &h in &h_charges {
            assert!(h > 0.0, "hydrogens should be positive in CH4");
        }

        let total: f64 = int.atoms.iter().map(|a| a.charge).sum();
        assert!(total.abs() < 1e-9, "total charge should be zero");
    }

    #[test]
    fn qeq_charged_system() {
        let methane = make_methane();
        let mut int = IntermediateSystem::from_system(&methane).unwrap();

        let config = QeqConfig {
            total_charge: 1.0,
            ..Default::default()
        };
        assign_qeq_charges(&mut int, &config).unwrap();

        let total: f64 = int.atoms.iter().map(|a| a.charge).sum();
        assert!((total - 1.0).abs() < 1e-9, "total charge should be +1");
    }
}
