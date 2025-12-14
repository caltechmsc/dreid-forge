use std::process::ExitCode;

mod cli;
mod commands;
mod config;
mod display;
mod io;

fn main() -> ExitCode {
    let cli = cli::parse();
    let ctx = display::Context::detect().with_quiet(match &cli.command {
        cli::Command::Bio(args) => args.io.quiet,
        cli::Command::Chem(args) => args.io.quiet,
    });

    if ctx.interactive {
        display::print_banner();
    }

    match commands::dispatch(cli.command, ctx) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            display::print_error(&e);
            ExitCode::FAILURE
        }
    }
}
