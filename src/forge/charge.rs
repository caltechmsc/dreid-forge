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

#[cfg(test)]
mod tests {
    use super::*;
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
