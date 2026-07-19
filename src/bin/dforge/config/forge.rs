use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use dreid_forge::{ForgeConfig, MpsimConfig};

use crate::cli::{
    ChargeOptions, HybridChargeOptions, MpsimOptions, PotentialOptions, QeqSolverOptions,
};
use crate::util::convert::{build_bio_charge_method, build_chem_charge_method};

pub fn build_bio_forge_config(
    charge: &ChargeOptions,
    hybrid: &HybridChargeOptions,
    qeq: &QeqSolverOptions,
    potential: &PotentialOptions,
    mpsim: &MpsimOptions,
) -> Result<ForgeConfig> {
    let rules = load_optional_file(&potential.rules, "typing rules")?;
    let params = load_optional_file(&potential.params, "force field parameters")?;

    Ok(ForgeConfig {
        rules,
        params,
        // MPSim's neutral termini need force-field charges to be meaningful, so
        // default to the hybrid method when the user left --charge unset.
        charge_method: build_bio_charge_method(charge, hybrid, qeq, mpsim.enabled),
        bond_potential: potential.bond_potential.into(),
        angle_potential: potential.angle_potential.into(),
        vdw_potential: potential.vdw_potential.into(),
        mpsim: build_mpsim_config(mpsim),
    })
}

pub fn build_chem_forge_config(
    charge: &ChargeOptions,
    qeq: &QeqSolverOptions,
    potential: &PotentialOptions,
    mpsim: &MpsimOptions,
) -> Result<ForgeConfig> {
    let rules = load_optional_file(&potential.rules, "typing rules")?;
    let params = load_optional_file(&potential.params, "force field parameters")?;

    Ok(ForgeConfig {
        rules,
        params,
        charge_method: build_chem_charge_method(charge, qeq),
        bond_potential: potential.bond_potential.into(),
        angle_potential: potential.angle_potential.into(),
        vdw_potential: potential.vdw_potential.into(),
        mpsim: build_mpsim_config(mpsim),
    })
}

fn build_mpsim_config(mpsim: &MpsimOptions) -> Option<MpsimConfig> {
    mpsim.enabled.then(MpsimConfig::default)
}

fn load_optional_file(
    path: &Option<std::path::PathBuf>,
    description: &str,
) -> Result<Option<String>> {
    path.as_ref()
        .map(|p| read_config_file(p, description))
        .transpose()
}

fn read_config_file(path: &Path, description: &str) -> Result<String> {
    fs::read_to_string(path)
        .with_context(|| format!("Failed to read {} file: {}", description, path.display()))
}
