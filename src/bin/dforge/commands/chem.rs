use anyhow::{Context, Result, bail};

use dreid_forge::forge;
use dreid_forge::io::{ChemReader, ChemWriter, Format, write_lammps_data, write_lammps_settings};

use crate::cli::ChemArgs;
use crate::config::{build_chem_forge_config, build_lammps_config, potential_names};
use crate::display::{
    Context as DisplayContext, Progress, print_atom_types, print_parameters, print_structure_info,
};
use crate::io::{
    OutputSpec, create_output, infer_chem_input_format, infer_output_format, open_input,
    stdin_is_tty, stdout_is_tty,
};

const TOTAL_STEPS: u8 = 3;

pub fn run_chem(args: ChemArgs, ctx: DisplayContext) -> Result<()> {
    if args.io.input.is_none() && stdin_is_tty() {
        bail!(
            "No input file specified and stdin is a terminal.\n\nUsage: dforge chem <INPUT> or pipe data via stdin."
        );
    }

    let input_format = resolve_input_format(&args)?;
    let output_specs = resolve_outputs(&args)?;

    let mut progress = Progress::new(ctx.interactive, TOTAL_STEPS);

    progress.step("Reading structure");
    let system = read_chem_structure(&args, input_format)?;

    let read_substeps = build_read_substeps(input_format);
    let read_substeps_ref: Vec<&str> = read_substeps.iter().map(|s| s.as_str()).collect();
    progress.complete_step("Reading structure", &read_substeps_ref);

    if ctx.interactive {
        print_structure_info(&system);
    }

    progress.step("Running DREIDING parameterization");
    let forge_config = build_chem_forge_config(&args.charge, &args.qeq, &args.potential)?;
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

fn build_read_substeps(format: Format) -> Vec<String> {
    let fmt_name = match format {
        Format::Mol2 => "MOL2",
        Format::Sdf => "SDF",
        _ => "molecule",
    };

    vec![
        format!("Parse {} file", fmt_name),
        "Perceive bond orders".to_string(),
        "Build connectivity graph".to_string(),
    ]
}

fn build_param_substeps(args: &ChemArgs) -> Vec<String> {
    use crate::cli::{AnglePotential, BondPotential, ChargeMethod, VdwPotential};

    let mut steps = Vec::new();

    if args.potential.rules.is_some() {
        steps.push("Assign atom types (custom rules)".to_string());
    } else {
        steps.push("Assign atom types (DREIDING rules)".to_string());
    }

    let charge_step = match args.charge.method {
        ChargeMethod::None => "Skip charge calculation".to_string(),
        ChargeMethod::Qeq | ChargeMethod::Hybrid => {
            format!(
                "Compute charges (QEq, total: {:.1}e)",
                args.charge.total_charge
            )
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
                Format::Mol2 => format!("Write MOL2 → {}", path_str),
                Format::Sdf => format!("Write SDF → {}", path_str),
                _ => format!("Write {:?} → {}", spec.format, path_str),
            }
        })
        .collect()
}

fn resolve_input_format(args: &ChemArgs) -> Result<Format> {
    if let Some(fmt) = args.input_format {
        return Ok(fmt.into());
    }

    if let Some(path) = &args.io.input {
        if let Some(fmt) = infer_chem_input_format(path) {
            return Ok(fmt);
        }
        bail!(
            "Cannot infer format from '{}'. Use --infmt to specify.",
            path.display()
        );
    }

    bail!("Reading from stdin requires --infmt");
}

fn resolve_outputs(args: &ChemArgs) -> Result<Vec<OutputSpec>> {
    if args.io.output.is_empty() {
        if stdout_is_tty() {
            bail!(
                "No output file specified and stdout is a terminal.\n\nUsage: dforge chem <INPUT> -o <OUTPUT> or pipe output."
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

fn read_chem_structure(args: &ChemArgs, format: Format) -> Result<dreid_forge::System> {
    let input = open_input(args.io.input.as_deref())?;
    ChemReader::new(input, format)
        .read()
        .context("Failed to read structure")
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
            Mol2 | Sdf => {
                ChemWriter::new(writer, spec.format)
                    .write(&forged.system)
                    .context("Failed to write structure file")?;
            }
            Bgf | Pdb | Mmcif => {
                bail!(
                    "Output format {:?} requires bio metadata (use 'dforge bio' instead)",
                    spec.format
                );
            }
        }
    }

    Ok(())
}
