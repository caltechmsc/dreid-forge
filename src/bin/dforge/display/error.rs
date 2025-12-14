use std::io::{self, Write};

use anyhow::Error;

use crate::util::text::wrap;

#[rustfmt::skip]
pub fn print_error(err: &Error) {
    let mut stderr = io::stderr().lock();

    let _ = writeln!(stderr);
    let _ = writeln!(stderr, "   ╔══════════════════════════════════════════════════════════════╗");
    let _ = writeln!(stderr, "   ║  ✗ Error                                                     ║");
    let _ = writeln!(stderr, "   ╟──────────────────────────────────────────────────────────────╢");

    let msg = err.to_string();
    for line in wrap(&msg, 59) {
        let _ = writeln!(stderr, "   ║  {:<59} ║", line);
    }

    let mut source = err.source();
    while let Some(cause) = source {
        let _ = writeln!(stderr, "   ╟──────────────────────────────────────────────────────────────╢");
        let _ = writeln!(stderr, "   ║  Caused by:                                                  ║");
        for line in wrap(&cause.to_string(), 59) {
            let _ = writeln!(stderr, "   ║    {:<57} ║", line);
        }
        source = cause.source();
    }

    if let Some(hints) = get_hints(err) {
        let _ = writeln!(stderr, "   ╟──────────────────────────────────────────────────────────────╢");
        let _ = writeln!(stderr, "   ║  Hints:                                                      ║");
        for hint in hints {
            let wrapped = wrap(&hint, 57);
            if let Some((first, rest)) = wrapped.split_first() {
                let _ = writeln!(stderr, "   ║    • {:<55} ║", first);
                for line in rest {
                    let _ = writeln!(stderr, "   ║      {:<55} ║", line);
                }
            }
        }
    }

    let _ = writeln!(stderr, "   ╚══════════════════════════════════════════════════════════════╝");
    let _ = writeln!(stderr);
}

fn get_hints(err: &Error) -> Option<Vec<String>> {
    use dreid_forge::ForgeError;
    use dreid_forge::io::Error as IoError;

    let msg = err.to_string().to_lowercase();

    let mut hints = Vec::new();

    if let Some(io_err) = err.downcast_ref::<IoError>() {
        match io_err {
            IoError::Io { source: _ } => {
                hints.push("I/O error: check file path and permissions".into());
                hints.push("Verify the file exists and is readable".into());
            }
            IoError::Parse {
                format,
                line,
                details,
            } => {
                hints.push(format!("Parsing error in {} at ~line {}", format, line));
                hints
                    .push("Inspect the file around the reported line for malformed entries".into());
                hints.push("Try specifying the input format explicitly with --infmt".into());
                hints.push(format!("Parser details: {}", details));
            }
            IoError::MissingMetadata(fmt) => {
                hints.push(format!(
                    "Output format '{}' requires biological metadata",
                    fmt
                ));
                hints.push("Try using bio mode or provide required metadata fields".into());
            }
            IoError::BioForgePreparation(_) | IoError::BioForgeIo(_) => {
                hints.push("Biomolecular preparation failed — check residue templates and protonation settings".into());
                hints.push("Ensure hetero residue templates are available or remove/replace unknown residues".into());
            }
            IoError::Conversion(_) => {
                hints.push("Data model conversion failed — check file format compatibility".into());
            }
            _ => {}
        }
    }

    if let Some(f_err) = err.downcast_ref::<ForgeError>() {
        match f_err {
            ForgeError::ParameterParse(_) => {
                hints.push("Failed to parse force field parameters; check TOML syntax".into());
            }
            ForgeError::RuleParse(_) => {
                hints.push("Custom typing rules failed to parse; verify TOML structure".into());
            }
            ForgeError::AtomTyping(_) => {
                hints.push("Atom typing failed — try providing custom typing rules or inspect unusual residues".into());
            }
            ForgeError::MissingParameter { atom_type, detail } => {
                hints.push(format!("Missing parameter for atom type '{}'", atom_type));
                hints.push(format!("Detail: {}", detail));
                hints.push("Provide a params file with the missing entries or use --params to load custom parameters".into());
            }
            ForgeError::ChargeCalculation(_) => {
                hints.push(
                    "Charge calculation failed — try --charge none to skip or check atom types"
                        .into(),
                );
            }
            ForgeError::EmptySystem => {
                hints.push(
                    "Input system appears empty — verify the input file contains atoms".into(),
                );
            }
            _ => {}
        }
    }

    if msg.contains("no such file") || msg.contains("not found") {
        hints.push("Check that the file path is correct".into());
        hints.push("Verify the file exists and is readable".into());
    }

    if msg.contains("permission denied") {
        hints.push("Check file permissions with `ls -la`".into());
        hints.push("Try running with appropriate permissions".into());
    }

    if msg.contains("invalid format") || msg.contains("parse error") {
        hints.push("Verify the input file format is correct".into());
        hints.push("Use --format to explicitly specify the format".into());
    }

    if msg.contains("pdb") {
        hints.push("Ensure the PDB file follows standard format".into());
        hints.push("Check for missing ATOM/HETATM records".into());
    }

    if msg.contains("mmcif") || msg.contains("cif") {
        hints.push("Ensure proper _atom_site loop structure".into());
        hints.push("Check for malformed data blocks".into());
    }

    if msg.contains("mol2") {
        hints.push("Verify @<TRIPOS>ATOM and @<TRIPOS>BOND sections".into());
        hints.push("Check atom type assignments".into());
    }

    if msg.contains("sdf") || msg.contains("mol file") {
        hints.push("Ensure proper V2000/V3000 format".into());
        hints.push("Check the counts line in the header".into());
    }

    if msg.contains("charge") || msg.contains("qeq") {
        hints.push("Check that all atoms have valid types".into());
        hints.push("Try using --charge none to skip".into());
    }

    if msg.contains("topology") || msg.contains("bond") {
        hints.push("Verify connectivity information is present".into());
        hints.push("Check bond distance cutoffs".into());
    }

    if msg.contains("empty") || msg.contains("no atoms") {
        hints.push("The input file may be empty or corrupt".into());
        hints.push("Check that the file contains valid data".into());
    }

    if hints.is_empty() { None } else { Some(hints) }
}
