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
   - Charge assignment (`none`, `qeq`, or `hybrid`)
   - Parameter generation (bond/angle/vdW functional forms selectable)
3. **Write Outputs**
   - BGF and/or structure formats (depending on command)

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

- If `--input` is omitted **and stdin is a terminal**, the command errors (prevents "waiting for input" by accident).
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

- Output format is `--outfmt` if provided, otherwise defaults to `bgf`.

Output inference mapping:

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

| Flag                | Values                               | Meaning                                                   |
| ------------------- | ------------------------------------ | --------------------------------------------------------- |
| `--infmt <FORMAT>`  | `pdb`, `mmcif`                       | Force input parsing format.                               |
| `--outfmt <FORMAT>` | `bgf`, `pdb`, `mmcif`, `mol2`, `sdf` | Output format for stdout output or the first `-o` output. |

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
- If a template does not match the structure's residue/atom naming, topology building may fail.

### Charge Calculation

| Flag                 | Values / Type           | Default | Meaning                                                    |
| -------------------- | ----------------------- | ------- | ---------------------------------------------------------- |
| `--charge <METHOD>`  | `none`, `qeq`, `hybrid` | `none`  | Charge calculation method.                                 |
| `--total-charge <Q>` | float                   | `0.0`   | Total system charge constraint used by QEq/Hybrid methods. |

`--charge` values:

- `none` — All charges remain zero (default).
- `qeq` — QEq charge equilibration for all atoms.
- `hybrid` — Hybrid mode: force field charges for biomolecules, QEq for ligands (bio only).

### Hybrid Charge Options (bio only)

These options are only used when `--charge hybrid` is selected.

| Flag                               | Values / Type                               | Default      | Meaning                                         |
| ---------------------------------- | ------------------------------------------- | ------------ | ----------------------------------------------- |
| `--protein-scheme <SCHEME>`        | `amber-ffsb`, `amber-ff03`, `charmm`        | `amber-ffsb` | Protein force field charge scheme.              |
| `--nucleic-scheme <SCHEME>`        | `amber`, `charmm`                           | `amber`      | Nucleic acid force field charge scheme.         |
| `--water-scheme <SCHEME>`          | `tip3p`, `tip3p-fb`, `spc`, `spc-e`, `opc3` | `tip3p`      | Water model charge scheme.                      |
| `--ligand <CONFIG>`                | string (repeatable)                         | none         | Per-ligand configuration (see below).           |
| `--default-ligand-method <METHOD>` | `vacuum`, `embedded`                        | `embedded`   | Default QEq method for unlisted ligands.        |
| `--default-ligand-cutoff <Å>`      | float                                       | `10.0`       | Default embedded QEq environment cutoff radius. |

**Force field scheme options:**

- `--protein-scheme`:

  - `amber-ffsb` (alias: `ffsb`) — AMBER ff99SB/ff14SB/ff19SB (default)
  - `amber-ff03` (alias: `ff03`) — AMBER ff03
  - `charmm` — CHARMM22/27/36/36m

- `--nucleic-scheme`:

  - `amber` — AMBER OL15/OL21/OL24/bsc1/OL3 (default)
  - `charmm` — CHARMM C27/C36

- `--water-scheme`:
  - `tip3p` — TIP3P (default)
  - `tip3p-fb` — TIP3P-FB
  - `spc` — SPC
  - `spc-e` — SPC/E
  - `opc3` — OPC3

**Ligand configuration format:** `CHAIN:RESID[:ICODE][:METHOD[:CUTOFF]]`

- `CHAIN` — Chain identifier (e.g., `A`, `L`)
- `RESID` — Residue sequence number
- `ICODE` — Optional insertion code (single character)
- `METHOD` — Optional QEq method: `vacuum` or `embedded`
- `CUTOFF` — Optional cutoff radius for embedded method (Å)

Examples:

```bash
# Use embedded QEq for ligand at chain A, residue 500
--ligand A:500

# Use vacuum QEq for ligand at chain L, residue 1
--ligand L:1:vacuum

# Use embedded QEq with 15 Å cutoff for residue 100 with insertion code B
--ligand A:100:B:embedded:15.0
```

### QEq Solver Options

These options tune the QEq solver and apply to both `qeq` and `hybrid` charge methods.

| Flag                       | Values / Type | Default | Meaning                                          |
| -------------------------- | ------------- | ------- | ------------------------------------------------ |
| `--qeq-tolerance <TOL>`    | float         | `1e-6`  | Convergence tolerance for charge equilibration.  |
| `--qeq-max-iter <N>`       | integer       | `100`   | Maximum iterations for QEq solver.               |
| `--qeq-lambda <λ>`         | float         | `0.5`   | Orbital screening parameter λ (Rappe–Goddard).   |
| `--qeq-hydrogen-scf`       | boolean       | `true`  | Enable hydrogen SCF (nonlinear hardness update). |
| `--qeq-basis <TYPE>`       | `gto`, `sto`  | `sto`   | Basis function type for Coulomb integrals.       |
| `--qeq-damping <STRATEGY>` | string        | `auto`  | SCF damping strategy.                            |

`--qeq-damping` values:

- `none` — No damping (fastest, may not converge for difficult systems).
- `fixed:<f>` — Fixed damping factor (0 < f ≤ 1).
- `auto` — Automatic adaptive damping (default).
- `auto:<f>` — Automatic damping with initial factor f.

### Force Field Options

These affect the DREIDING parameterization stage.

| Flag                       | Values / Type              | Default          | Meaning                                                                                           |
| -------------------------- | -------------------------- | ---------------- | ------------------------------------------------------------------------------------------------- |
| `--bond-potential <TYPE>`  | `harmonic`, `morse`        | `harmonic`       | Bond potential functional form.                                                                   |
| `--angle-potential <TYPE>` | `cosine`, `theta-harmonic` | `theta-harmonic` | Angle potential functional form. (`theta` is accepted as an alias.)                               |
| `--vdw-potential <TYPE>`   | `lj`, `exp6`               | `lj`             | van der Waals potential functional form. (`lennard-jones` and `buckingham` are accepted aliases.) |
| `--rules <FILE>`           | path                       | none             | Custom typing rules (TOML file).                                                                  |
| `--params <FILE>`          | path                       | none             | Custom force field parameters (TOML file).                                                        |

---

## Command: `dforge chem`

Parameterize small molecules (MOL2/SDF).

Usage:

```text
dforge chem [OPTIONS]
```

### Format Flags

| Flag                | Values        | Meaning                                                   |
| ------------------- | ------------- | --------------------------------------------------------- |
| `--infmt <FORMAT>`  | `mol2`, `sdf` | Force input parsing format.                               |
| `--outfmt <FORMAT>` | `mol2`, `sdf` | Output format for stdout output or the first `-o` output. |

### Charge Calculation

| Flag                 | Values / Type           | Default | Meaning                                     |
| -------------------- | ----------------------- | ------- | ------------------------------------------- |
| `--charge <METHOD>`  | `none`, `qeq`, `hybrid` | `none`  | Charge calculation method.                  |
| `--total-charge <Q>` | float                   | `0.0`   | Total system charge constraint used by QEq. |

Notes:

- For `dforge chem`, `--charge hybrid` behaves the same as `--charge qeq` (no biological metadata).

### QEq Solver Options

Same as `dforge bio` — see above.

### Force Field Options

These affect the DREIDING parameterization stage.

| Flag                       | Values / Type              | Default          | Meaning                                                                                           |
| -------------------------- | -------------------------- | ---------------- | ------------------------------------------------------------------------------------------------- |
| `--bond-potential <TYPE>`  | `harmonic`, `morse`        | `harmonic`       | Bond potential functional form.                                                                   |
| `--angle-potential <TYPE>` | `cosine`, `theta-harmonic` | `theta-harmonic` | Angle potential functional form. (`theta` is accepted as an alias.)                               |
| `--vdw-potential <TYPE>`   | `lj`, `exp6`               | `lj`             | van der Waals potential functional form. (`lennard-jones` and `buckingham` are accepted aliases.) |
| `--rules <FILE>`           | path                       | none             | Custom typing rules (TOML file).                                                                  |
| `--params <FILE>`          | path                       | none             | Custom force field parameters (TOML file).                                                        |

### Output Format Constraints

- `dforge chem` can write: `mol2`, `sdf`.
- Attempting to write `bgf`, `pdb`, or `mmcif` from `dforge chem` is rejected because those formats require biomolecular metadata.

---

## Exit Status

- Success: exit code `0`
- Failure: exit code `1`

---

## Workflow Examples

Please refer to the [examples directory](https://github.com/caltechmsc/dreid-forge/tree/main/examples) for end-to-end usage examples and sample input/output files.
