#[derive(Debug, Clone)]
pub struct ForgeConfig {
    pub rules: Option<String>,
    pub params: Option<String>,
    pub charge_method: ChargeMethod,
    pub bond_potential: BondPotentialType,
    pub angle_potential: AnglePotentialType,
    pub vdw_potential: VdwPotentialType,
}

impl Default for ForgeConfig {
    fn default() -> Self {
        Self {
            rules: None,
            params: None,
            charge_method: ChargeMethod::None,
            bond_potential: BondPotentialType::Harmonic,
            angle_potential: AnglePotentialType::ThetaHarmonic,
            vdw_potential: VdwPotentialType::LennardJones,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub enum ChargeMethod {
    #[default]
    None,
    Qeq(QeqConfig),
}

#[derive(Debug, Clone)]
pub struct QeqConfig {
    pub total_charge: f64,
    pub solver_options: SolverOptions,
}

impl Default for QeqConfig {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            solver_options: SolverOptions {
                hydrogen_scf: false,
                ..SolverOptions::default()
            },
        }
    }
}

pub use cheq::SolverOptions;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondPotentialType {
    #[default]
    Harmonic,
    Morse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AnglePotentialType {
    CosineHarmonic,
    #[default]
    ThetaHarmonic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VdwPotentialType {
    #[default]
    LennardJones,
    Exponential6,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_values() {
        let config = ForgeConfig::default();
        assert!(config.rules.is_none());
        assert!(config.params.is_none());
        assert!(matches!(config.charge_method, ChargeMethod::None));
        assert_eq!(config.bond_potential, BondPotentialType::Harmonic);
        assert_eq!(config.angle_potential, AnglePotentialType::ThetaHarmonic);
        assert_eq!(config.vdw_potential, VdwPotentialType::LennardJones);
    }

    #[test]
    fn qeq_config_default() {
        let qeq = QeqConfig::default();
        assert_eq!(qeq.total_charge, 0.0);
        assert_eq!(qeq.solver_options.tolerance, 1.0e-6);
    }

    #[test]
    fn charge_method_with_qeq() {
        let config = ForgeConfig {
            charge_method: ChargeMethod::Qeq(QeqConfig {
                total_charge: -1.0,
                ..Default::default()
            }),
            ..Default::default()
        };

        if let ChargeMethod::Qeq(qeq) = config.charge_method {
            assert_eq!(qeq.total_charge, -1.0);
        } else {
            panic!("Expected Qeq variant");
        }
    }
}
