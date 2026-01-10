use anyhow::{Context, Result, bail};

use dreid_forge::forge;
use dreid_forge::io::{
    BioReader, BioWriter, ChemWriter, Format, write_bgf, write_lammps_data, write_lammps_settings,
};

use crate::cli::BioArgs;
use crate::config::{
    build_bio_forge_config, build_clean_config, build_lammps_config, build_protonate_config,
    build_solvate_config, build_topology_config, potential_names,
};
use crate::display::{
    Context as DisplayContext, Progress, print_atom_types, print_chain_breakdown, print_parameters,
    print_structure_info,
};
use crate::io::{
    OutputSpec, create_output, infer_bio_input_format, infer_output_format, open_input,
    stdin_is_tty, stdout_is_tty,
};

const TOTAL_STEPS: u8 = 3;

pub fn run_bio(args: BioArgs, ctx: DisplayContext) -> Result<()> {
    if args.io.input.is_none() && stdin_is_tty() {
        bail!(
            "No input file specified and stdin is a terminal.\n\nUsage: dforge bio <INPUT> or pipe data via stdin."
        );
    }

    let input_format = resolve_input_format(&args)?;
    let output_specs = resolve_outputs(&args)?;

    let mut progress = Progress::new(ctx.interactive, TOTAL_STEPS);

    progress.step("Preparing structure");
    let system = read_bio_structure(&args, input_format)?;

    let prep_substeps = build_prep_substeps(&args, input_format);
    let prep_substeps_ref: Vec<&str> = prep_substeps.iter().map(|s| s.as_str()).collect();
    progress.complete_step("Preparing structure", &prep_substeps_ref);

    if ctx.interactive {
        print_structure_info(&system);
        print_chain_breakdown(&system);
    }

    progress.step("Running DREIDING parameterization");
    let forge_config =
        build_bio_forge_config(&args.charge, &args.hybrid, &args.qeq, &args.potential)?;
    let forged = forge(&system, &forge_config).context("Parameterization failed")?;

    let param_substeps = build_param_substeps(&args);
    let param_substeps_ref: Vec<&str> = param_substeps.iter().map(|s| s.as_str()).collect();
    progress.complete_step("Running DREIDING parameterization", &param_substeps_ref);

    if ctx.interactive {
        let (bond_type, angle_type, vdw_type) = potential_names(&args.potential);
        print_atom_types(&system, Some(&forged));
        print_parameters(&forged, bond_type, angle_type, vdw_type);
    }

    progress.step("Writing output");
    let lammps_config = build_lammps_config(&args.lammps);
    write_outputs(&forged, &output_specs, &lammps_config)?;

    let write_substeps = build_write_substeps(&output_specs);
    let write_substeps_ref: Vec<&str> = write_substeps.iter().map(|s| s.as_str()).collect();
    progress.complete_step("Writing output", &write_substeps_ref);

    progress.finish();

    Ok(())
}

fn build_prep_substeps(args: &BioArgs, format: Format) -> Vec<String> {
    use crate::cli::HisStrategy;

    let mut steps = Vec::new();

    let fmt_name = match format {
        Format::Pdb => "PDB",
        Format::Mmcif => "mmCIF",
        _ => "structure",
    };
    steps.push(format!("Parse {} file", fmt_name));

    let mut clean_actions = Vec::new();
    if args.clean.no_water {
        clean_actions.push("water");
    }
    if args.clean.no_ions {
        clean_actions.push("ions");
    }
    if args.clean.no_hetero {
        clean_actions.push("hetero");
    }
    if args.clean.no_hydrogens {
        clean_actions.push("hydrogens");
    }
    if !args.clean.remove.is_empty() {
        clean_actions.push("specified residues");
    }
    if !clean_actions.is_empty() {
        steps.push(format!("Remove {}", clean_actions.join(", ")));
    }

    let his_method = match args.protonation.his {
        HisStrategy::Network => "network analysis",
        HisStrategy::Hid => "HID",
        HisStrategy::Hie => "HIE",
        HisStrategy::Random => "random",
    };
    steps.push(format!(
        "Add hydrogens (pH {:.1}, HIS: {})",
        args.protonation.ph, his_method
    ));

    if args.solvation.solvate {
        steps.push(format!(
            "Solvate (margin: {:.1} Å, spacing: {:.1} Å)",
            args.solvation.box_margin, args.solvation.spacing
        ));
        if args.solvation.target_charge != 0 {
            steps.push(format!(
                "Neutralize ({:?}/{:?} ions, target: {}e)",
                args.solvation.cation, args.solvation.anion, args.solvation.target_charge
            ));
        } else {
            steps.push(format!(
                "Neutralize ({:?}/{:?} ions)",
                args.solvation.cation, args.solvation.anion
            ));
        }
    }

    steps.push(format!(
        "Detect disulfide bonds (cutoff: {:.1} Å)",
        args.topology.ss_cutoff
    ));
    if !args.topology.templates.is_empty() {
        steps.push(format!(
            "Load {} hetero template(s)",
            args.topology.templates.len()
        ));
    }
    steps.push("Build molecular topology".to_string());

    steps
}

fn build_param_substeps(args: &BioArgs) -> Vec<String> {
    use crate::cli::{AnglePotential, BondPotential, ChargeMethod, VdwPotential};
    use crate::util::convert::{
        ligand_method_display_name, nucleic_scheme_display_name, protein_scheme_display_name,
        water_scheme_display_name,
    };

    let mut steps = Vec::new();

    steps.push("Build molecular graph".to_string());

    if args.potential.rules.is_some() {
        steps.push("Assign atom types (custom rules)".to_string());
    } else {
        steps.push("Assign atom types (DREIDING rules)".to_string());
    }

    let charge_step = match args.charge.method {
        ChargeMethod::None => "Skip charge calculation".to_string(),
        ChargeMethod::Qeq => {
            format!(
                "Compute charges (QEq, total: {:.1}e)",
                args.charge.total_charge
            )
        }
        ChargeMethod::Hybrid => {
            let protein = protein_scheme_display_name(args.hybrid.protein_scheme);
            let nucleic = nucleic_scheme_display_name(args.hybrid.nucleic_scheme);
            let water = water_scheme_display_name(args.hybrid.water_scheme);
            let default_method = ligand_method_display_name(args.hybrid.default_ligand_method);
            let ligand_count = args.hybrid.ligands.len();
            if ligand_count > 0 {
                format!(
                    "Compute charges (Hybrid: {}/{}/{}, {} ligand(s), default: {})",
                    protein, nucleic, water, ligand_count, default_method
                )
            } else {
                format!(
                    "Compute charges (Hybrid: {}/{}/{}, ligand default: {})",
                    protein, nucleic, water, default_method
                )
            }
        }
    };
    steps.push(charge_step);

    let bond_type = match args.potential.bond_potential {
        BondPotential::Harmonic => "harmonic",
        BondPotential::Morse => "Morse",
    };
    let angle_type = match args.potential.angle_potential {
        AnglePotential::Cosine => "cosine",
        AnglePotential::ThetaHarmonic => "θ-harmonic",
    };
    let vdw_type = match args.potential.vdw_potential {
        VdwPotential::Lj => "LJ 12-6",
        VdwPotential::Exp6 => "exp-6",
    };

    steps.push(format!(
        "Generate potentials ({} bond, {} angle, {} vdW)",
        bond_type, angle_type, vdw_type
    ));

    if args.potential.params.is_some() {
        steps.push("Apply custom parameters".to_string());
    }

    steps.push("Compute H-bond pairs".to_string());

    steps
}

fn build_write_substeps(specs: &[OutputSpec]) -> Vec<String> {
    specs
        .iter()
        .map(|spec| {
            let path_str = spec
                .path
                .as_ref()
                .map(|p| {
                    p.file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned()
                })
                .unwrap_or_else(|| "stdout".to_string());

            match spec.format {
                Format::LammpsData => format!("Write LAMMPS data → {}", path_str),
                Format::LammpsSettings => format!("Write LAMMPS settings → {}", path_str),
                Format::Bgf => format!("Write BGF → {}", path_str),
                Format::Pdb => format!("Write PDB → {}", path_str),
                Format::Mmcif => format!("Write mmCIF → {}", path_str),
                Format::Mol2 => format!("Write MOL2 → {}", path_str),
                Format::Sdf => format!("Write SDF → {}", path_str),
            }
        })
        .collect()
}

fn resolve_input_format(args: &BioArgs) -> Result<Format> {
    if let Some(fmt) = args.input_format {
        return Ok(fmt.into());
    }

    if let Some(path) = &args.io.input {
        if let Some(fmt) = infer_bio_input_format(path) {
            return Ok(fmt);
        }
        bail!(
            "Cannot infer format from '{}'. Use --infmt to specify.",
            path.display()
        );
    }

    bail!("Reading from stdin requires --infmt");
}

fn resolve_outputs(args: &BioArgs) -> Result<Vec<OutputSpec>> {
    if args.io.output.is_empty() {
        if stdout_is_tty() {
            bail!(
                "No output file specified and stdout is a terminal.\n\nUsage: dforge bio <INPUT> -o <OUTPUT> or pipe output."
            );
        }
        let format = args
            .output_format
            .map(|f| f.into())
            .unwrap_or(Format::LammpsData);
        return Ok(vec![OutputSpec { path: None, format }]);
    }

    let mut specs = Vec::with_capacity(args.io.output.len());

    let first = &args.io.output[0];
    let first_format = if let Some(fmt) = args.output_format {
        fmt.into()
    } else if let Some(fmt) = infer_output_format(first) {
        fmt
    } else {
        bail!(
            "Cannot infer format from '{}'. Use --outfmt to specify.",
            first.display()
        );
    };
    specs.push(OutputSpec {
        path: Some(first.clone()),
        format: first_format,
    });

    for path in &args.io.output[1..] {
        let format = infer_output_format(path).ok_or_else(|| {
            anyhow::anyhow!(
                "Cannot infer format from '{}'. Use explicit extension.",
                path.display()
            )
        })?;
        specs.push(OutputSpec {
            path: Some(path.clone()),
            format,
        });
    }

    Ok(specs)
}

fn read_bio_structure(args: &BioArgs, format: Format) -> Result<dreid_forge::System> {
    let input = open_input(args.io.input.as_deref())?;

    let clean_config = build_clean_config(&args.clean);
    let protonation_config = build_protonate_config(&args.protonation);
    let solvation_config = build_solvate_config(&args.solvation);
    let topology_config = build_topology_config(&args.topology)?;

    let mut reader = BioReader::new(input, format)
        .clean(clean_config)
        .protonate(protonation_config)
        .topology(topology_config);

    if let Some(solv) = solvation_config {
        reader = reader.solvate(solv);
    }

    reader.read().context("Failed to read structure")
}

fn write_outputs(
    forged: &dreid_forge::ForgedSystem,
    specs: &[OutputSpec],
    lammps_config: &dreid_forge::io::LammpsConfig,
) -> Result<()> {
    use Format::*;

    for spec in specs {
        let mut writer = create_output(spec.path.as_deref())?;

        match spec.format {
            LammpsData => {
                write_lammps_data(&mut writer, forged, lammps_config)
                    .context("Failed to write LAMMPS data file")?;
            }
            LammpsSettings => {
                write_lammps_settings(&mut writer, forged, lammps_config)
                    .context("Failed to write LAMMPS settings file")?;
            }
            Bgf => {
                write_bgf(&mut writer, forged).context("Failed to write BGF file")?;
            }
            Pdb | Mmcif => {
                BioWriter::new(writer, spec.format)
                    .write(&forged.system)
                    .context("Failed to write structure file")?;
            }
            Mol2 | Sdf => {
                ChemWriter::new(writer, spec.format)
                    .write(&forged.system)
                    .context("Failed to write structure file")?;
            }
        }
    }

    Ok(())
}
