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

    if let Some(hints) = HintCollector::collect(err) {
        let _ = writeln!(stderr, "   ╟──────────────────────────────────────────────────────────────╢");
        let _ = writeln!(stderr, "   ║  Hints:                                                      ║");
        for hint in hints {
            let wrapped = wrap(&hint, 55);
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

struct HintCollector {
    hints: Vec<String>,
    has_typed_hints: bool,
}

impl HintCollector {
    fn new() -> Self {
        Self {
            hints: Vec::new(),
            has_typed_hints: false,
        }
    }

    fn collect(err: &Error) -> Option<Vec<String>> {
        let mut collector = Self::new();

        collector.collect_io_hints(err);
        collector.collect_forge_hints(err);

        if !collector.has_typed_hints {
            collector.collect_fallback_hints(err);
        }

        if collector.hints.is_empty() {
            None
        } else {
            Some(collector.hints)
        }
    }

    fn add(&mut self, hint: impl Into<String>) {
        self.hints.push(hint.into());
    }

    fn mark_typed(&mut self) {
        self.has_typed_hints = true;
    }

    fn collect_io_hints(&mut self, err: &Error) {
        use dreid_forge::io::Error as IoError;

        let Some(io_err) = err.downcast_ref::<IoError>() else {
            return;
        };

        self.mark_typed();

        match io_err {
            IoError::Io { source } => {
                self.collect_std_io_hints(source);
            }

            IoError::Parse { format, line, .. } => {
                self.add(format!(
                    "Parser encountered an issue near line {} in {} format",
                    line, format
                ));
                self.add("Inspect the file around that line for malformed entries");
                self.add("Try specifying --infmt to ensure correct format detection");
                self.add_format_specific_parse_hints(*format);
            }

            IoError::UnsupportedReadFormat(fmt) => {
                self.add(format!("The '{}' format cannot be used for input", fmt));
                self.add("Supported input formats: pdb, mmcif, mol2, sdf");
            }

            IoError::UnsupportedWriteFormat(fmt) => {
                self.add(format!("The '{}' format cannot be used for output", fmt));
                self.add("Supported output formats: pdb, mmcif, mol2, sdf, bgf");
            }

            IoError::MissingMetadata(fmt) => {
                self.add(format!(
                    "The '{}' format requires biological metadata (residues, chains)",
                    fmt
                ));
                self.add("Use 'bio' mode to read structures with metadata preservation");
                self.add("Or choose a format that doesn't require metadata (mol2, sdf)");
            }

            IoError::BioForgePreparation(msg) => {
                self.add("Biomolecular preparation failed during structure processing");
                self.collect_bioforge_prep_hints(msg);
            }

            IoError::BioForgeIo(msg) => {
                self.add("Underlying bio-forge I/O operation failed");
                self.collect_bioforge_io_hints(msg);
            }

            IoError::Conversion(msg) => {
                self.add("Data model conversion failed between internal formats");
                self.collect_conversion_hints(msg);
            }
        }
    }

    fn collect_std_io_hints(&mut self, source: &std::io::Error) {
        use std::io::ErrorKind;

        match source.kind() {
            ErrorKind::NotFound => {
                self.add("File or directory not found");
                self.add("Check the path spelling and ensure the file exists");
            }

            ErrorKind::PermissionDenied => {
                self.add("Permission denied accessing the file");
                self.add("Check file permissions with `ls -la`");
                self.add("Ensure you have read/write access as needed");
            }

            ErrorKind::AlreadyExists => {
                self.add("File already exists");
                self.add("Use a different output path or remove the existing file");
            }

            ErrorKind::InvalidInput => {
                self.add("Invalid input provided to I/O operation");
                self.add("Check that the path is valid and properly formatted");
            }

            ErrorKind::InvalidData => {
                self.add("File contains invalid or corrupt data");
                self.add("Verify the file is not truncated or corrupted");
            }

            ErrorKind::UnexpectedEof => {
                self.add("Unexpected end of file encountered");
                self.add("The file may be truncated or incomplete");
            }

            ErrorKind::WriteZero => {
                self.add("Failed to write data (disk full?)");
                self.add("Check available disk space");
            }

            ErrorKind::BrokenPipe => {
                self.add("Broken pipe — output consumer terminated");
                self.add("This may occur when piping to commands like `head`");
            }

            _ => {
                self.add("I/O operation failed");
                self.add("Check file path, permissions, and disk space");
            }
        }
    }

    fn add_format_specific_parse_hints(&mut self, format: dreid_forge::io::Format) {
        use dreid_forge::io::Format;

        match format {
            Format::Pdb => {
                self.add("PDB: Check ATOM/HETATM record formatting (columns 1-80)");
                self.add("PDB: Ensure proper spacing in coordinate fields");
            }

            Format::Mmcif => {
                self.add("mmCIF: Verify _atom_site loop structure and data alignment");
                self.add("mmCIF: Check for unquoted strings containing special chars");
            }

            Format::Mol2 => {
                self.add("MOL2: Verify @<TRIPOS>ATOM and @<TRIPOS>BOND sections");
                self.add("MOL2: Check column alignment in atom records");
            }

            Format::Sdf => {
                self.add("SDF: Verify V2000/V3000 format in the counts line");
                self.add("SDF: Check atom block and bond block formatting");
            }

            Format::Bgf => {
                // This is a write-only format, parsing hints not applicable
            }
        }
    }

    fn collect_bioforge_prep_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("residue") || msg_lower.contains("template") {
            self.add("Unknown residue encountered without a template");
            self.add("Provide a template via --template <file.mol2>");
            self.add("Or use --no-hetero to remove all hetero atoms");
        } else if msg_lower.contains("proton") || msg_lower.contains("hydrogen") {
            self.add("Protonation step failed");
            self.add("Use --no-hydrogens to remove existing hydrogens first");
        } else if msg_lower.contains("topology") || msg_lower.contains("bond") {
            self.add("Topology building failed");
            self.add("Check --ss-cutoff for disulfide bond detection");
        } else if msg_lower.contains("solv") || msg_lower.contains("water") {
            self.add("Solvation step failed");
            self.add("Check --solvate and solvation-related options");
        } else {
            self.add("Check residue names and ensure templates are available");
            self.add("Use --template to provide hetero residue templates");
        }
    }

    fn collect_bioforge_io_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("parse") || msg_lower.contains("read") {
            self.add("Failed to parse input file via bio-forge");
            self.add("Verify the input file format is correct");
        } else if msg_lower.contains("write") {
            self.add("Failed to write output via bio-forge");
            self.add("Check output path and permissions");
        } else {
            self.add("Bio-forge encountered an I/O issue");
            self.add("Verify input file format and output path");
        }
    }

    fn collect_conversion_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("coord") || msg_lower.contains("position") {
            self.add("Coordinate data conversion failed");
            self.add("Input may contain invalid or missing coordinates");
        } else if msg_lower.contains("bond") || msg_lower.contains("connect") {
            self.add("Bond/connectivity conversion failed");
            self.add("Input may have inconsistent bonding information");
        } else {
            self.add("Internal data structure conversion failed");
            self.add("This may indicate a bug — please report if reproducible");
        }
    }

    fn collect_forge_hints(&mut self, err: &Error) {
        use dreid_forge::ForgeError;

        let Some(forge_err) = err.downcast_ref::<ForgeError>() else {
            return;
        };

        self.mark_typed();

        match forge_err {
            ForgeError::ParameterParse(_) => {
                self.add("Force field parameter file has invalid TOML syntax");
                self.add("Check for missing quotes, brackets, or invalid values");
                self.add("Validate TOML syntax with an online validator");
            }

            ForgeError::RuleParse(msg) => {
                self.add("Custom atom typing rules failed to parse");
                self.add("Verify the typing rules TOML structure");
                self.collect_rule_parse_hints(msg);
            }

            ForgeError::AtomTyping(msg) => {
                self.add("Atom type assignment failed for one or more atoms");
                self.collect_atom_typing_hints(msg);
            }

            ForgeError::ChargeCalculation(msg) => {
                self.add("QEq charge calculation failed");
                self.collect_charge_calc_hints(msg);
            }

            ForgeError::MissingBioMetadata => {
                self.add("Hybrid charge method requires biological metadata");
                self.add("Use 'dforge bio' for structures with residue information");
                self.add("Or use --charge=qeq for small molecules with 'dforge chem'");
            }

            ForgeError::HybridChargeAssignment {
                chain_id,
                residue_id,
                residue_name,
                detail,
            } => {
                self.add(format!(
                    "Hybrid charge assignment failed for {} {} in chain {}",
                    residue_name, residue_id, chain_id
                ));
                self.add(format!("Issue: {}", detail));
                self.add("Check that the residue has correct atom naming");
                self.add("Verify force field scheme supports this residue type");
            }

            ForgeError::MissingParameter { atom_type, detail } => {
                self.add(format!("No parameters found for atom type '{}'", atom_type));
                self.add(format!("Missing: {}", detail));
                self.add("Provide custom parameters via --params");
                self.add("Or add the missing entry to your parameter file");
            }

            ForgeError::InvalidBond { i, j, detail } => {
                self.add(format!("Invalid bond between atoms {} and {}", i, j));
                self.add(format!("Issue: {}", detail));
                self.add("Check input structure for unusual bonding");
            }

            ForgeError::EmptySystem => {
                self.add("Input system contains no atoms");
                self.add("Verify the input file is not empty or corrupt");
                self.add("Check that the correct format was detected");
            }

            ForgeError::Conversion(msg) => {
                self.add("Internal data conversion failed during parameterization");
                self.add(format!("Details: {}", msg));
            }
        }
    }

    fn collect_rule_parse_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("smarts") || msg_lower.contains("pattern") {
            self.add("Check SMARTS pattern syntax in typing rules");
        } else if msg_lower.contains("priority") {
            self.add("Verify rule priority values are valid integers");
        } else {
            self.add("Review the typing rules file structure");
        }
    }

    fn collect_atom_typing_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("element") || msg_lower.contains("unsupported") {
            self.add("Structure contains unsupported elements");
            self.add("DREIDING supports: H, C, N, O, F, Al, Si, P, S, Cl, Ga, Ge, As, Se, Br, In, Sn, Sb, Te, I, ...");
        } else if msg_lower.contains("hybridization") || msg_lower.contains("geometry") {
            self.add("Could not determine hybridization for some atoms");
            self.add("Check for unusual bonding patterns or coordination");
        } else if msg_lower.contains("no rule") || msg_lower.contains("match") {
            self.add("No typing rule matched some atoms");
            self.add("Provide custom rules via --rules for unusual cases");
        } else {
            self.add("Try providing custom typing rules with --rules");
            self.add("Check for unusual molecular structures");
        }
    }

    fn collect_charge_calc_hints(&mut self, msg: &str) {
        let msg_lower = msg.to_lowercase();

        if msg_lower.contains("converge") || msg_lower.contains("iteration") {
            self.add("QEq solver failed to converge");
            self.add("Try --charge none to skip charge calculation");
        } else if msg_lower.contains("electronegativity") || msg_lower.contains("parameter") {
            self.add("Missing QEq parameters for some atom types");
            self.add("Ensure all atom types have electronegativity values");
        } else {
            self.add("Charge calculation failed");
            self.add("Try --charge none to skip charge calculation");
        }
    }

    fn collect_fallback_hints(&mut self, err: &Error) {
        let msg = error_chain_text(err);

        if msg.contains("terminal") || msg.contains("stdin") || msg.contains("tty") {
            self.add("Input appears to be from a terminal");
            self.add("Provide input via -i/--input or pipe data to stdin");
            return;
        }

        if msg.contains("no such file") || msg.contains("not found") {
            self.add("Check that the file path is correct");
            self.add("Verify the file exists and is readable");
            return;
        }

        if msg.contains("permission denied") {
            self.add("Check file permissions with `ls -la`");
            self.add("Ensure you have the required access rights");
            return;
        }

        if msg.contains("empty") && !self.has_typed_hints {
            self.add("Input appears to be empty");
            self.add("Verify the input contains valid molecular data");
        }
    }
}

fn error_chain_text(err: &Error) -> String {
    let mut text = String::new();

    text.push_str(&err.to_string());

    let mut source = err.source();
    while let Some(cause) = source {
        text.push('\n');
        text.push_str(&cause.to_string());
        source = cause.source();
    }

    text.to_lowercase()
}
