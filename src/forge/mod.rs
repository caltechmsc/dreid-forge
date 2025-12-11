mod charge;
mod config;
mod error;
mod intermediate;
mod paramgen;
mod params;
mod typer;

pub use config::{
    AnglePotentialType, BondPotentialType, ChargeMethod, ForgeConfig, QeqConfig, SolverOptions,
    VdwPotentialType,
};
pub use error::Error;

use crate::model::system::System;
use crate::model::topology::ForgedSystem;

pub fn forge(system: &System, config: &ForgeConfig) -> Result<ForgedSystem, Error> {
    let ff_params = params::load_parameters(config.params.as_deref())?;

    let mut intermediate = intermediate::IntermediateSystem::from_system(system)?;

    typer::assign_atom_types(&mut intermediate, config.rules.as_deref())?;

    charge::assign_charges(&mut intermediate, &config.charge_method)?;

    let forged = paramgen::generate_parameters(system, &intermediate, ff_params, config)?;

    Ok(forged)
}
