use std::io::{self, Write};
use std::sync::LazyLock;

const VERSION: &str = env!("CARGO_PKG_VERSION");

const BANNER_ART: &str = r#"
   ██████╗ ██████╗ ███████╗██╗██████╗
   ██╔══██╗██╔══██╗██╔════╝██║██╔══██╗
   ██║  ██║██████╔╝█████╗  ██║██║  ██║
   ██║  ██║██╔══██╗██╔══╝  ██║██║  ██║
   ██████╔╝██║  ██║███████╗██║██████╔╝
   ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝╚═════╝
                    ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
                    ██╔════╝██╔═══██╗██╔══██╗██╔════╝ ██╔════╝
                    █████╗  ██║   ██║██████╔╝██║  ███╗█████╗
                    ██╔══╝  ██║   ██║██╔══██╗██║   ██║██╔══╝
                    ██║     ╚██████╔╝██║  ██║╚██████╔╝███████╗
                    ╚═╝      ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚══════╝

                ╔═╗╔═╗╦  ╔╦╗╔═╗╔═╗╦ ╦  ╔╦╗╔═╗╔═╗
                ║  ╠═╣║   ║ ║╣ ║  ╠═╣  ║║║╚═╗║
                ╚═╝╩ ╩╩═╝ ╩ ╚═╝╚═╝╩ ╩  ╩ ╩╚═╝╚═╝

   ───────────────────────────────────────────────────────────
         Tony Kan  ·  Ted Yu  ·  William A. Goddard III
   ───────────────────────────────────────────────────────────
   "#;

static BANNER_FOR_HELP: LazyLock<String> = LazyLock::new(|| format!("\n{BANNER_ART}"));

pub fn banner_for_help() -> &'static str {
    &BANNER_FOR_HELP
}

pub fn print_banner() {
    let mut stderr = io::stderr().lock();
    let _ = writeln!(stderr);
    let _ = writeln!(stderr, "{BANNER_ART}");
    let _ = writeln!(stderr);
    let _ = writeln!(
        stderr,
        "   DREIDING Force Field Parameterization              v{VERSION}"
    );
    let _ = writeln!(stderr);
}
