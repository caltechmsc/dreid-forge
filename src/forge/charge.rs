use super::config::{ChargeMethod, QeqConfig};
use super::error::Error;
use super::intermediate::IntermediateSystem;
use cheq::{QEqSolver, get_default_parameters};

pub fn assign_charges(system: &mut IntermediateSystem, method: &ChargeMethod) -> Result<(), Error> {
    match method {
        ChargeMethod::None => Ok(()),
        ChargeMethod::Qeq(config) => assign_qeq_charges(system, config),
    }
}

fn assign_qeq_charges(system: &mut IntermediateSystem, config: &QeqConfig) -> Result<(), Error> {
    let params = get_default_parameters();
    let solver = QEqSolver::new(params).with_options(config.solver_options);

    let atoms = system.atoms();

    let result = solver.solve(atoms, config.total_charge)?;

    for (atom, &charge) in system.atoms.iter_mut().zip(result.charges.iter()) {
        atom.charge = charge;
    }

    Ok(())
}
