use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result};

use dreid_forge::io::{
    CleanConfig, ProtonationConfig, SolvateConfig, Template, TopologyConfig, read_mol2_template,
};

use crate::cli::{CleanOptions, ProtonationOptions, SolvationOptions, TopologyOptions};

pub fn build_clean_config(opts: &CleanOptions) -> CleanConfig {
    CleanConfig {
        remove_water: opts.no_water,
        remove_ions: opts.no_ions,
        remove_hydrogens: opts.no_hydrogens,
        remove_hetero: opts.no_hetero,
        remove_residue_names: opts.remove.iter().cloned().collect(),
        keep_residue_names: opts.keep.iter().cloned().collect(),
    }
}

pub fn build_protonate_config(opts: &ProtonationOptions) -> ProtonationConfig {
    ProtonationConfig {
        target_ph: opts.ph,
        remove_existing_h: true,
        his_strategy: opts.his.into(),
    }
}

pub fn build_solvate_config(opts: &SolvationOptions) -> Option<SolvateConfig> {
    if !opts.solvate {
        return None;
    }

    Some(SolvateConfig {
        margin: opts.box_margin,
        water_spacing: opts.spacing,
        vdw_cutoff: opts.vdw_cutoff,
        remove_existing: true,
        cations: vec![opts.cation.into()],
        anions: vec![opts.anion.into()],
        target_charge: opts.target_charge,
        rng_seed: opts.seed,
    })
}

pub fn build_topology_config(opts: &TopologyOptions) -> Result<TopologyConfig> {
    let templates = load_templates(&opts.templates)?;

    Ok(TopologyConfig {
        hetero_templates: templates,
        disulfide_bond_cutoff: opts.ss_cutoff,
    })
}

fn load_templates(paths: &[impl AsRef<Path>]) -> Result<Vec<Template>> {
    paths
        .iter()
        .map(|path| {
            let path = path.as_ref();
            let file = File::open(path)
                .with_context(|| format!("Failed to open template: {}", path.display()))?;
            read_mol2_template(BufReader::new(file))
                .with_context(|| format!("Failed to parse template: {}", path.display()))
        })
        .collect()
}
