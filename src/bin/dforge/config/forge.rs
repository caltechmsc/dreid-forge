use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use dreid_forge::ForgeConfig;

use crate::cli::ForgeOptions;
use crate::util::convert::build_charge_method;

pub fn build_forge_config(opts: &ForgeOptions) -> Result<ForgeConfig> {
    let rules = opts
        .rules
        .as_ref()
        .map(|p| read_config_file(p, "typing rules"))
        .transpose()?;

    let params = opts
        .params
        .as_ref()
        .map(|p| read_config_file(p, "force field parameters"))
        .transpose()?;

    Ok(ForgeConfig {
        rules,
        params,
        charge_method: build_charge_method(opts.charge, opts.total_charge),
        bond_potential: opts.bond_potential.into(),
        angle_potential: opts.angle_potential.into(),
        vdw_potential: opts.vdw_potential.into(),
    })
}

fn read_config_file(path: &Path, description: &str) -> Result<String> {
    fs::read_to_string(path)
        .with_context(|| format!("Failed to read {} file: {}", description, path.display()))
}
