use crate::io::{error::Error, util};
use crate::model::system::System;
use bio_forge as bf;
use std::io::Write;

pub fn write<W: Write>(writer: W, system: &System) -> Result<(), Error> {
    let topo = util::to_bio_topology(system)?;
    bf::io::write_mmcif_topology(writer, &topo).map_err(Error::from)
}
