use crate::io::{BioReader, error::Error, util};
use crate::model::system::System;
use bio_forge as bf;
use std::io::BufRead;

pub fn read<R: BufRead>(builder: BioReader<R>) -> Result<System, Error> {
    let bio_context = bf::io::IoContext::new_default();
    let mut bio_struct = bf::io::read_mmcif_structure(builder.reader, &bio_context)?;

    if let Some(clean_config) = builder.clean_config {
        let bf_clean_config = util::to_bf_clean_config(clean_config);
        bf::ops::clean_structure(&mut bio_struct, &bf_clean_config)?;
    }

    bf::ops::repair_structure(&mut bio_struct)?;

    let bf_hydro_config = util::to_bf_hydro_config(builder.protonation_config);
    bf::ops::add_hydrogens(&mut bio_struct, &bf_hydro_config)?;

    if let Some(solvate_config) = builder.solvate_config {
        let bf_solvate_config = util::to_bf_solvate_config(solvate_config);
        bf::ops::solvate_structure(&mut bio_struct, &bf_solvate_config)?;
    }

    let bf_topo = util::build_topology(bio_struct, &builder.topology_config)?;
    let system = util::from_bio_topology(bf_topo)?;

    Ok(system)
}
