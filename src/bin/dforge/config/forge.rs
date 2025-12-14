use dreid_forge::ForgeConfig;

use crate::cli::ForgeOptions;
use crate::util::convert::build_charge_method;

pub fn build_forge_config(opts: &ForgeOptions) -> ForgeConfig {
    ForgeConfig {
        rules: opts
            .rules
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        params: opts
            .params
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        charge_method: build_charge_method(opts.charge, opts.total_charge),
        bond_potential: opts.bond_potential.into(),
        angle_potential: opts.angle_potential.into(),
        vdw_potential: opts.vdw_potential.into(),
    }
}
