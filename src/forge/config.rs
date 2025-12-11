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
