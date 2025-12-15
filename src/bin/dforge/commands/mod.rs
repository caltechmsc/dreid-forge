mod bio;
mod chem;

use bio::run_bio;
use chem::run_chem;

use anyhow::Result;

use crate::cli::Command;
use crate::display::Context;

pub fn dispatch(command: Command, ctx: Context) -> Result<()> {
    match command {
        Command::Bio(args) => run_bio(args, ctx),
        Command::Chem(args) => run_chem(args, ctx),
    }
}
