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
