# DREID-Forge CLI Manual

`dforge` is the command-line interface shipped with the `dreid-forge` crate. It reads molecular structures, performs DREIDING parameterization (biomolecular preparation (only for `dforge bio`) → typing → optional charge calculation → potential generation), and writes simulation-ready outputs.

---

## Installation

### Option 1: Cargo Install (from source, needs Rust toolchain)

```bash
cargo install dreid-forge
```

This installs the `dforge` binary into your Cargo bin directory (usually `~/.cargo/bin`).

### Option 2: Precompiled Binaries

Download the latest release for your platform from the [GitHub Releases](https://github.com/caltechmsc/dreid-forge/releases).

> **Note**: Make sure to add the binary location to your `PATH` environment variable if necessary.

---

## Concepts and Data Flow

At a high level, the CLI runs:

1. **Read / Prepare Structure**
   - `dforge bio`: PDB/mmCIF with biomolecular preparation (clean → protonate → optional solvate → topology)
   - `dforge chem`: MOL2/SDF direct chemistry parsing
2. **Forge (DREIDING Parameterization)**
   - Atom typing (DREIDING rules or custom typing rules)
   - Charge assignment (`qeq` or `none`)
   - Parameter generation (bond/angle/vdW functional forms selectable)
3. **Write Outputs**
   - LAMMPS data/settings and/or structure formats (depending on command)

---

## Global CLI

Top-level usage:

```text
dforge <COMMAND>
```

Commands:

- `bio` (alias: `b`) — Parameterize biological macromolecules (PDB/mmCIF)
- `chem` (alias: `c`) — Parameterize small molecules (MOL2/SDF)

Top-level flags:

- `-h, --help` — Print help
- `-V, --version` — Print version

---

## I/O Model (shared)

Both `bio` and `chem` share these I/O flags:

| Flag                  | Type              | Meaning                                                        |
| --------------------- | ----------------- | -------------------------------------------------------------- |
| `-i, --input <FILE>`  | path              | Input file. If omitted, reads from stdin (requires `--infmt`). |
| `-o, --output <FILE>` | path (repeatable) | Output file path(s). Repeat to write multiple outputs.         |
| `-q, --quiet`         | flag              | Suppress progress output (disables interactive UI).            |

### Stdin / Stdout Safeguards

- If `--input` is omitted **and stdin is a terminal**, the command errors (prevents “waiting for input” by accident).
- If `--output` is omitted **and stdout is a terminal**, the command errors (prevents dumping structured output into your terminal).

### Input Format Inference

If `--infmt` is not provided and `--input` is a file, the format is inferred from the extension:

- Bio input inference (`dforge bio`):
  - `.pdb`, `.ent` → `pdb`
  - `.cif`, `.mmcif` → `mmcif`
- Chem input inference (`dforge chem`):
  - `.mol2` → `mol2`
  - `.sdf`, `.mol` → `sdf`

Reading from stdin requires an explicit `--infmt`.

### Output Format Inference and Multi-Output Rules

If one or more `--output` paths are provided:

- The **first** output format is:
  - `--outfmt` if provided, otherwise inferred from the first output path extension.
- Any **additional** output formats are inferred from each additional output path extension.

If no `--output` is provided (stdout output):

- Output format is `--outfmt` if provided, otherwise defaults to `lammps-data`.

Output inference mapping:

- `.data`, `.lammps` → `lammps-data`
- `.settings` → `lammps-settings`
- `.bgf` → `bgf`
- `.pdb`, `.ent` → `pdb`
- `.cif`, `.mmcif` → `mmcif`
- `.mol2` → `mol2`
- `.sdf`, `.mol` → `sdf`

---

## Command: `dforge bio`

Parameterize biological macromolecules (PDB/mmCIF).

Usage:

```text
dforge bio [OPTIONS]
```

### Format Flags

| Flag                | Values                                                                 | Meaning                                                   |
| ------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------- |
| `--infmt <FORMAT>`  | `pdb`, `mmcif`                                                         | Force input parsing format.                               |
| `--outfmt <FORMAT>` | `lammps-data`, `lammps-settings`, `bgf`, `pdb`, `mmcif`, `mol2`, `sdf` | Output format for stdout output or the first `-o` output. |

Aliases:

- `lammps-data` accepts `data`
- `lammps-settings` accepts `settings`

### Structure Cleaning

Flags:

| Flag             | Type                   | Meaning                                         |
| ---------------- | ---------------------- | ----------------------------------------------- |
| `--no-water`     | flag                   | Remove water molecules.                         |
| `--no-ions`      | flag                   | Remove ions.                                    |
| `--no-hydrogens` | flag                   | Remove existing hydrogens.                      |
| `--no-hetero`    | flag                   | Remove hetero atoms (HETATM).                   |
| `--remove <RES>` | list (comma-separated) | Remove specific residues by residue name.       |
| `--keep <RES>`   | list (comma-separated) | Keep only these residue names (whitelist mode). |

Notes:

- `--remove` and `--keep` are both accepted; if both are used, the effective behavior is determined by the underlying cleaning stage.

### Protonation

| Flag               | Type  | Default                        | Meaning                                     |
| ------------------ | ----- | ------------------------------ | ------------------------------------------- |
| `--ph <PH>`        | float | **N/A** (by each residue name) | Target pH for protonation state assignment. |
| `--his <STRATEGY>` | enum  | `network`                      | Histidine tautomer strategy.                |

`--his` values:

- `hid` — Always HID (Nδ protonated)
- `hie` — Always HIE (Nε protonated)
- `random` — Random selection
- `network` — H-bond network analysis

### Solvation

Solvation is opt-in via `--solvate`.

| Flag                  | Type  | Default | Meaning                                     |
| --------------------- | ----- | ------- | ------------------------------------------- |
| `--solvate`           | flag  | off     | Enable solvation (add water box).           |
| `--solv-margin <Å>`   | float | `10.0`  | Water box margin around solute (Å).         |
| `--solv-spacing <Å>`  | float | `3.1`   | Water molecule spacing (Å).                 |
| `--solv-cutoff <Å>`   | float | `2.4`   | Minimum solute-water distance (Å).          |
| `--solv-cation <ION>` | enum  | `na`    | Cation type for charge targeting.           |
| `--solv-anion <ION>`  | enum  | `cl`    | Anion type for charge targeting.            |
| `--solv-charge <Q>`   | int   | `0`     | Target net charge after ion addition.       |
| `--solv-seed <SEED>`  | uint  | none    | Random seed for reproducible ion placement. |

`--solv-cation` values: `na`, `k`, `mg`, `ca`, `li`, `zn`

`--solv-anion` values: `cl`, `br`, `i`, `f`

### Topology

| Flag                | Type              | Default | Meaning                               |
| ------------------- | ----------------- | ------- | ------------------------------------- |
| `--ss-cutoff <Å>`   | float             | `2.2`   | Disulfide bond detection cutoff (Å).  |
| `--template <FILE>` | path (repeatable) | none    | Hetero residue template files (MOL2). |

Template notes:

- Templates are used during topology building for non-standard residues (ligands/cofactors).
- If a template does not match the structure’s residue/atom naming, topology building may fail.

### Force Field Options

These affect the DREIDING parameterization stage.

| Flag                       | Values / Type              | Default          | Meaning                                                                                           |
| -------------------------- | -------------------------- | ---------------- | ------------------------------------------------------------------------------------------------- |
| `--charge <METHOD>`        | `qeq`, `none`              | `qeq`            | Charge calculation method.                                                                        |
| `--total-charge <Q>`       | float                      | `0.0`            | Total system charge constraint used by QEq.                                                       |
| `--bond-potential <TYPE>`  | `harmonic`, `morse`        | `harmonic`       | Bond potential functional form.                                                                   |
| `--angle-potential <TYPE>` | `cosine`, `theta-harmonic` | `theta-harmonic` | Angle potential functional form. (`theta` is accepted as an alias.)                               |
| `--vdw-potential <TYPE>`   | `lj`, `exp6`               | `lj`             | van der Waals potential functional form. (`lennard-jones` and `buckingham` are accepted aliases.) |
| `--rules <FILE>`           | path                       | none             | Custom typing rules (TOML file).                                                                  |
| `--params <FILE>`          | path                       | none             | Custom force field parameters (TOML file).                                                        |

### LAMMPS Output Options

These affect how LAMMPS outputs are generated.

| Flag                      | Type                       | Default    | Meaning                                                                         |
| ------------------------- | -------------------------- | ---------- | ------------------------------------------------------------------------------- |
| `--lmp-cutoff <Å>`        | float                      | `12.0`     | Non-bonded cutoff distance (Å).                                                 |
| `--lmp-hbond-inner <Å>`   | float                      | `10.0`     | Hydrogen bond inner cutoff (Å).                                                 |
| `--lmp-hbond-outer <Å>`   | float                      | `12.0`     | Hydrogen bond outer cutoff (Å).                                                 |
| `--lmp-hbond-angle <DEG>` | float                      | `90.0`     | Hydrogen bond angle cutoff (degrees).                                           |
| `--lmp-margin <Å>`        | float                      | `5.0`      | Box margin for non-periodic systems (Å).                                        |
| `--lmp-system <TYPE>`     | `periodic`, `non-periodic` | `periodic` | Boundary condition type. (`shrink` is accepted as an alias for `non-periodic`.) |

---

## Command: `dforge chem`

Parameterize small molecules (MOL2/SDF).

Usage:

```text
dforge chem [OPTIONS]
```

### Format Flags

| Flag                | Values                                          | Meaning                                                   |
| ------------------- | ----------------------------------------------- | --------------------------------------------------------- |
| `--infmt <FORMAT>`  | `mol2`, `sdf`                                   | Force input parsing format.                               |
| `--outfmt <FORMAT>` | `lammps-data`, `lammps-settings`, `mol2`, `sdf` | Output format for stdout output or the first `-o` output. |

Aliases:

- `lammps-data` accepts `data`
- `lammps-settings` accepts `settings`

### Force Field Options

These affect the DREIDING parameterization stage.

| Flag                       | Values / Type              | Default          | Meaning                                                                                           |
| -------------------------- | -------------------------- | ---------------- | ------------------------------------------------------------------------------------------------- |
| `--charge <METHOD>`        | `qeq`, `none`              | `qeq`            | Charge calculation method.                                                                        |
| `--total-charge <Q>`       | float                      | `0.0`            | Total system charge constraint used by QEq.                                                       |
| `--bond-potential <TYPE>`  | `harmonic`, `morse`        | `harmonic`       | Bond potential functional form.                                                                   |
| `--angle-potential <TYPE>` | `cosine`, `theta-harmonic` | `theta-harmonic` | Angle potential functional form. (`theta` is accepted as an alias.)                               |
| `--vdw-potential <TYPE>`   | `lj`, `exp6`               | `lj`             | van der Waals potential functional form. (`lennard-jones` and `buckingham` are accepted aliases.) |
| `--rules <FILE>`           | path                       | none             | Custom typing rules (TOML file).                                                                  |
| `--params <FILE>`          | path                       | none             | Custom force field parameters (TOML file).                                                        |

### LAMMPS Output Options

These affect how LAMMPS outputs are generated.

| Flag                      | Type                       | Default    | Meaning                                                                         |
| ------------------------- | -------------------------- | ---------- | ------------------------------------------------------------------------------- |
| `--lmp-cutoff <Å>`        | float                      | `12.0`     | Non-bonded cutoff distance (Å).                                                 |
| `--lmp-hbond-inner <Å>`   | float                      | `10.0`     | Hydrogen bond inner cutoff (Å).                                                 |
| `--lmp-hbond-outer <Å>`   | float                      | `12.0`     | Hydrogen bond outer cutoff (Å).                                                 |
| `--lmp-hbond-angle <DEG>` | float                      | `90.0`     | Hydrogen bond angle cutoff (degrees).                                           |
| `--lmp-margin <Å>`        | float                      | `5.0`      | Box margin for non-periodic systems (Å).                                        |
| `--lmp-system <TYPE>`     | `periodic`, `non-periodic` | `periodic` | Boundary condition type. (`shrink` is accepted as an alias for `non-periodic`.) |

### Output Format Constraints

- `dforge chem` can write: `lammps-data`, `lammps-settings`, `mol2`, `sdf`.
- Attempting to write `bgf`, `pdb`, or `mmcif` from `dforge chem` is rejected because those formats require biomolecular metadata.

---

## Exit Status

- Success: exit code `0`
- Failure: exit code `1`

---

## Workflow Examples

Please refer to the [examples directory](https://github.com/caltechmsc/dreid-forge/tree/main/examples) for end-to-end usage examples and sample input/output files.
