use dreid_forge::io::LammpsConfig;

use crate::cli::LammpsOptions;

pub fn build_lammps_config(opts: &LammpsOptions) -> LammpsConfig {
    LammpsConfig {
        nonbonded_cutoff: opts.nonbonded_cutoff,
        hbond_cutoff_inner: opts.hbond_inner,
        hbond_cutoff_outer: opts.hbond_outer,
        hbond_angle_cutoff: opts.hbond_angle,
        aabb_margin: opts.aabb_margin,
        system_type: opts.system_type.into(),
    }
}
