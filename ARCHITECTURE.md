# DREID-Forge Architecture

This document provides a comprehensive guide to the internal design, data flow, and algorithms of the `dreid-forge` library. It is intended for contributors, integrators, and users who wish to understand the technical details behind automated DREIDING force field parameterization.

## 1. System Overview

DREID-Forge is architected as a **multi-stage transformation pipeline** that converts raw molecular geometry into simulation-ready force field parameters. The system integrates four external crates (`bio-forge`, `dreid-typer`, `cheq`, `ffcharge`) and orchestrates them through a unified API.

```mermaid
flowchart TD
    subgraph "Input Sources"
        MOL2["MOL2"]
        SDF["SDF"]
        PDB["PDB"]
        mmCIF["mmCIF"]
        TPL["MOL2 Template"]
        MEM["In-Memory"]
    end

    subgraph "I/O Layer"
        direction TB
        CR["ChemReader"]
        BR["BioReader"]
        RT["read_mol2_template"]
    end

    subgraph "bio-forge Pipeline"
        direction TB
        BF_PARSE["Parse Structure"]
        BF_CLEAN["Clean"]
        BF_REPAIR["Repair"]
        BF_HYDRO["Protonate"]
        BF_SOLV["Solvate"]
        BF_TOPO["Build Topology"]
    end

    subgraph "Model Layer"
        SYS["System"]
    end

    subgraph "Forge Pipeline"
        direction TB
        INT["IntermediateSystem"]
        TYPER["Atom Typing<br><i>(dreid-typer)</i>"]
        CHARGE["Charge Calculation<br><i>(cheq/ffcharge)</i>"]
        PGEN["Parameter Generation"]
    end

    subgraph "Output"
        FS["ForgedSystem"]
    end

    subgraph "Writers"
        direction LR
        W_BGF["BGF Writer"]
        W_LAM["LAMMPS Writer"]
        W_PDB["PDB/mmCIF Writer"]
        W_CHEM["MOL2/SDF Writer"]
    end

    MOL2 --> CR
    SDF --> CR
    PDB --> BR
    mmCIF --> BR
    TPL --> RT
    RT --> BR
    MEM --> SYS

    CR --> SYS
    BR --> BF_PARSE
    BF_PARSE --> BF_CLEAN
    BF_CLEAN --> BF_REPAIR
    BF_REPAIR --> BF_HYDRO
    BF_HYDRO --> BF_SOLV
    BF_SOLV --> BF_TOPO
    BF_TOPO --> SYS

    SYS --> INT
    INT --> TYPER
    TYPER --> CHARGE
    CHARGE --> PGEN
    PGEN --> FS

    FS --> W_BGF
    FS --> W_LAM
    SYS --> W_PDB
    SYS --> W_CHEM
```

### Core Components

| Component           | Description                                                                                                     |
| ------------------- | --------------------------------------------------------------------------------------------------------------- |
| **I/O Layer**       | Readers and writers for molecular structure file formats                                                        |
| **Model Layer**     | Neutral data structures (`System`, `Atom`, `Bond`, `BioMetadata`)                                               |
| **Forge Pipeline**  | The parameterization engine: typing → charging → parameter generation                                           |
| **External Crates** | `bio-forge` (structure prep), `dreid-typer` (atom typing), `cheq` (QEq charges), `ffcharge` (classical charges) |

---

## 2. I/O Layer Architecture

The I/O layer provides a dual-track system for reading molecular structures: **ChemReader** for small molecules and **BioReader** for biomolecules with full preparation capabilities.

```mermaid
flowchart TB
    subgraph "File Formats"
        direction LR
        F_MOL2["MOL2<br><i>.mol2</i>"]
        F_SDF["SDF/MOL<br><i>.sdf, .mol</i>"]
        F_PDB["PDB<br><i>.pdb</i>"]
        F_CIF["mmCIF<br><i>.cif</i>"]
        F_TPL["MOL2 Template<br><i>hetero residues</i>"]
    end

    subgraph "Chemistry Track"
        CR["ChemReader"]
        CR_MOL2["mol2::reader"]
        CR_SDF["sdf::reader"]
    end

    subgraph "Bio Track"
        BR["BioReader"]
        RT["read_mol2_template"]
    end

    subgraph "bio-forge Internal"
        direction TB
        BF_IO["bf::io::read_*"]
        BF_CLEAN["bf::ops::clean_structure"]
        BF_REPAIR["bf::ops::repair_structure"]
        BF_HYDRO["bf::ops::add_hydrogens"]
        BF_SOLV["bf::ops::solvate_structure"]
        BF_TOPO["bf::ops::TopologyBuilder"]
    end

    subgraph "Conversion Layer"
        CONV["util::from_bio_topology"]
    end

    subgraph "Output"
        SYS["System<br><i>atoms + bonds + metadata</i>"]
    end

    F_MOL2 --> CR
    F_SDF --> CR
    CR --> CR_MOL2
    CR --> CR_SDF
    CR_MOL2 --> SYS
    CR_SDF --> SYS

    F_PDB --> BR
    F_CIF --> BR
    F_TPL --> RT
    RT --> BR

    BR --> BF_IO
    BF_IO --> BF_CLEAN
    BF_CLEAN --> BF_REPAIR
    BF_REPAIR --> BF_HYDRO
    BF_HYDRO --> BF_SOLV
    BF_SOLV --> BF_TOPO
    BF_TOPO --> CONV
    CONV --> SYS
```

### 2.1 ChemReader: Direct Small-Molecule Parsing

`ChemReader` provides lightweight parsing for chemistry-focused formats without biological context:

- **MOL2 Reader** — Parses `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, and `@<TRIPOS>BOND` sections. Element inference uses atom type strings with fallback to atom names.
- **SDF Reader** — Parses MDL V2000 connection tables with fixed-width atom blocks and bond blocks. V3000 format is explicitly unsupported.

Both readers produce a `System` with `atoms`, `bonds`, and no biological metadata.

### 2.2 BioReader: Biomolecular Preparation Pipeline

`BioReader` wraps the `bio-forge` crate to provide a complete biomolecular preparation workflow:

```mermaid
flowchart LR
    subgraph "Configuration"
        CC["CleanConfig"]
        PC["ProtonationConfig"]
        SC["SolvateConfig"]
        TC["TopologyConfig"]
    end

    subgraph "Pipeline Stages"
        S1["1. Parse"]
        S2["2. Clean"]
        S3["3. Repair"]
        S4["4. Protonate"]
        S5["5. Solvate"]
        S6["6. Topology"]
    end

    CC --> S2
    PC --> S4
    SC --> S5
    TC --> S6

    S1 --> S2 --> S3 --> S4 --> S5 --> S6
```

| Stage     | bio-forge Function                                    | Purpose                                                |
| --------- | ----------------------------------------------------- | ------------------------------------------------------ |
| Parse     | `io::read_pdb_structure` / `io::read_mmcif_structure` | Raw structure ingestion with alias resolution          |
| Clean     | `ops::clean_structure`                                | Remove water, ions, hydrogens, or specific residues    |
| Repair    | `ops::repair_structure`                               | Reconstruct missing heavy atoms via template alignment |
| Protonate | `ops::add_hydrogens`                                  | Add hydrogens with pH-aware protonation states         |
| Solvate   | `ops::solvate_structure`                              | Add water box and counterions                          |
| Topology  | `ops::TopologyBuilder`                                | Build covalent bonds from residue templates            |

### 2.3 Writers: Output Format Support

| Writer         | Function                       | Input                       | Output                         |
| -------------- | ------------------------------ | --------------------------- | ------------------------------ |
| **ChemWriter** | `mol2::writer`, `sdf::writer`  | `System`                    | Chemical formats               |
| **BioWriter**  | `pdb::writer`, `mmcif::writer` | `System` with `BioMetadata` | Biological formats             |
| **BGF**        | `write_bgf`                    | `ForgedSystem`              | Biograf format with atom types |
| **LAMMPS**     | `write_lammps_package`         | `ForgedSystem`              | Data + settings files          |

### 2.4 Data Model Conversion

The `io::util` module provides bidirectional conversion between `dreid-forge` and `bio-forge` data models:

```mermaid
flowchart LR
    subgraph "dreid-forge Types"
        DF_SYS["System"]
        DF_ATOM["Atom"]
        DF_BOND["Bond"]
        DF_META["BioMetadata"]
    end

    subgraph "bio-forge Types"
        BF_TOPO["bf::Topology"]
        BF_STRUCT["bf::Structure"]
        BF_ATOM["bf::Atom"]
    end

    BF_TOPO -->|"from_bio_topology"| DF_SYS
    DF_SYS -->|"to_bio_topology"| BF_TOPO

    DF_ATOM <-->|"convert_element"| BF_ATOM
    DF_META <-->|"convert_residue_info"| BF_STRUCT
```

---

## 3. Forge Pipeline

The `forge()` function is the central entry point that orchestrates the complete parameterization workflow. It transforms a `System` into a `ForgedSystem` through four sequential stages.

```mermaid
flowchart TD
    subgraph "Input"
        SYS["System<br><i>atoms + bonds</i>"]
        CFG["ForgeConfig<br><i>options</i>"]
    end

    subgraph "Stage 1: Conversion"
        S1["IntermediateSystem::from_system"]
        INT["IntermediateSystem<br><i>+ neighbor lists</i>"]
    end

    subgraph "Stage 2: Atom Typing"
        S2["typer::assign_atom_types"]
        subgraph "dreid-typer"
            DT_P["Perception<br><i>rings, aromaticity, hybridization</i>"]
            DT_T["Typing Engine<br><i>rule evaluation</i>"]
            DT_B["Topology Builder<br><i>angles, dihedrals, impropers</i>"]
        end
    end

    subgraph "Stage 3: Charge Calculation"
        S3["charge::assign_charges"]
        subgraph "cheq"
            CQ_B["Build Invariant System"]
            CQ_S["SCF Iteration"]
            CQ_R["Return Charges"]
        end
    end

    subgraph "Stage 4: Parameter Generation"
        S4["paramgen::generate_parameters"]
        PG_AT["Collect Atom Types"]
        PG_BP["Generate Bond Potentials"]
        PG_AP["Generate Angle Potentials"]
        PG_DP["Generate Dihedral Potentials"]
        PG_IP["Generate Improper Potentials"]
        PG_VDW["Generate VdW Potentials"]
        PG_HB["Generate H-Bond Potentials"]
    end

    subgraph "Output"
        FS["ForgedSystem<br><i>ready for simulation</i>"]
    end

    SYS --> S1
    CFG --> S1
    S1 --> INT

    INT --> S2
    S2 --> DT_P --> DT_T --> DT_B

    DT_B --> S3
    S3 --> CQ_B --> CQ_S --> CQ_R

    CQ_R --> S4
    S4 --> PG_AT --> PG_BP --> PG_AP --> PG_DP --> PG_IP --> PG_VDW --> PG_HB

    PG_HB --> FS
```

### 3.1 Stage 1: System Conversion

`IntermediateSystem::from_system` prepares the molecular data for parameterization:

```mermaid
flowchart LR
    subgraph "Input"
        SYS["System"]
        ATOMS["atoms: Vec&lt;Atom&gt;"]
        BONDS["bonds: Vec&lt;Bond&gt;"]
    end

    subgraph "Processing"
        CONV["Convert atoms"]
        NEIGH["Build neighbor lists"]
        VALID["Validate bonds"]
    end

    subgraph "Output"
        INT["IntermediateSystem"]
        IA["IntermediateAtom[]<br><i>+ neighbors, charge=0</i>"]
        IB["IntermediateBond[]<br><i>+ physical_order=None</i>"]
    end

    ATOMS --> CONV --> IA
    BONDS --> VALID --> IB
    BONDS --> NEIGH --> IA
```

**Key operations:**

- Validates all bond indices are within bounds
- Builds bidirectional neighbor lists for each atom
- Initializes empty atom types and zero charges
- Returns `Error::EmptySystem` for empty input

### 3.2 Stage 2: Atom Typing via dreid-typer

The typing stage delegates to `dreid-typer` for DREIDING atom type assignment:

```mermaid
flowchart TD
    subgraph "dreid-forge"
        INT["IntermediateSystem"]
        BUILD["build_molecular_graph"]
        APPLY["apply_topology"]
    end

    subgraph "dreid-typer Pipeline"
        MG["MolecularGraph<br><i>connectivity only</i>"]

        subgraph "Phase 1: Perception"
            P1["Ring Detection"]
            P2["Kekulé Expansion"]
            P3["Electron Assignment"]
            P4["Aromaticity Detection"]
            P5["Resonance Analysis"]
            P6["Hybridization Inference"]
        end

        AM["AnnotatedMolecule<br><i>fully characterized</i>"]

        subgraph "Phase 2: Typing"
            TE["TyperEngine"]
            RULES["TOML Rules"]
            FP["Fixed-Point Iteration"]
        end

        TYPES["Vec&lt;String&gt; types"]

        subgraph "Phase 3: Building"
            TB["TopologyBuilder"]
            TERMS["Angles, Dihedrals, Impropers"]
        end

        MT["MolecularTopology"]
    end

    INT --> BUILD --> MG
    MG --> P1 --> P2 --> P3 --> P4 --> P5 --> P6 --> AM

    AM --> TE
    RULES --> TE
    TE --> FP --> TYPES

    AM --> TB
    TYPES --> TB
    TB --> TERMS --> MT

    MT --> APPLY --> INT
```

**Perception passes (in order):**

| Pass             | Module                      | Output                             |
| ---------------- | --------------------------- | ---------------------------------- |
| 1. Rings         | `perception::rings`         | `is_in_ring`, `smallest_ring_size` |
| 2. Kekulization  | `perception::kekulize`      | Explicit bond orders               |
| 3. Electrons     | `perception::electrons`     | `formal_charge`, `lone_pairs`      |
| 4. Aromaticity   | `perception::aromaticity`   | `is_aromatic`, `is_anti_aromatic`  |
| 5. Resonance     | `perception::resonance`     | `is_resonant`                      |
| 6. Hybridization | `perception::hybridization` | `hybridization`, `steric_number`   |

**Typing engine:**

- Evaluates TOML rules sorted by priority (descending)
- Uses fixed-point iteration for `neighbor_types` dependencies
- Converges when no atom's type changes between rounds

**Applied results:**

- `IntermediateAtom.atom_type` ← assigned type string (e.g., `"C_3"`, `"O_R"`)
- `IntermediateAtom.hybridization` ← `Hybridization` enum
- `IntermediateBond.physical_order` ← `PhysicalBondOrder` (Single/Double/Triple/Resonant)
- `IntermediateSystem.angles` ← enumerated angle terms
- `IntermediateSystem.dihedrals` ← enumerated proper dihedral terms
- `IntermediateSystem.impropers` ← enumerated improper dihedral terms

### 3.3 Stage 3: Charge Calculation

The charge module supports three methods: **None** (zero charges), **QEq** (global charge equilibration), and **Hybrid** (classical force field charges for biomolecules + QEq for ligands).

```mermaid
flowchart TD
    subgraph "Input"
        INT["IntermediateSystem"]
        CM["ChargeMethod"]
    end

    subgraph "Method Dispatch"
        NONE["None"]
        QEQ["Qeq"]
        HYB["Hybrid"]
    end

    subgraph "None Path"
        N_OUT["charges = 0.0"]
    end

    subgraph "QEq Path"
        Q_SOLVE["cheq::QEqSolver"]
        Q_OUT["Equilibrated charges"]
    end

    subgraph "Hybrid Path"
        H_CLASS["classify_atoms"]
        H_FIXED["assign_fixed_charges<br><i>(ffcharge)</i>"]
        H_LIGAND["assign_ligand_charges<br><i>(cheq)</i>"]
    end

    INT --> CM
    CM --> NONE --> N_OUT
    CM --> QEQ --> Q_SOLVE --> Q_OUT
    CM --> HYB --> H_CLASS --> H_FIXED --> H_LIGAND
```

#### 3.3.1 Global QEq Method

When `ChargeMethod::Qeq` is configured, partial charges are computed using the `cheq` crate:

```mermaid
flowchart TD
    subgraph "cheq Solver"
        PARAMS(["ElementData<br><i>χ, J, r, n</i>"])
        BUILD(["Build Invariant System"])
        MATRIX(["Coefficient Matrix<br><i>J_ii + J_ij(r) terms</i>"])
        RHS(["RHS Vector<br><i>-χ_i terms</i>"])
        SOLVE(["Solve Linear System"])
        SCF(["SCF Iteration<br><i>(if hydrogen_scf)</i>"])
    end

    subgraph "Output"
        CHARGES["charges: Vec&lt;f64&gt;"]
        MU["μ_eq: equilibrated potential"]
    end

    PARAMS --> BUILD
    BUILD --> MATRIX
    BUILD --> RHS
    MATRIX --> SOLVE
    RHS --> SOLVE
    SOLVE --> SCF --> CHARGES
    SOLVE --> MU
```

**QEq algorithm summary:**

1. **Build system:** Construct matrix $A$ and vector $b$ where:

   - Diagonal: $A_{ii} = J_i$ (atomic hardness)
   - Off-diagonal: $A_{ij} = J_{ij}(r_{ij})$ (screened Coulomb)
   - Constraint row: charge conservation

2. **Solve:** $A (q, \mu)^T = b$ yields charges and chemical potential

3. **Iterate:** If hydrogen SCF is enabled, update hydrogen hardness based on charge and re-solve

#### 3.3.2 Hybrid Charge Method

The hybrid method combines classical force field charges for biomolecules with QEq for ligands. This requires biological metadata (`BioMetadata`) to classify atoms.

```mermaid
flowchart TD
    subgraph "Atom Classification"
        META["BioMetadata"]
        CLASS["classify_atoms"]
        PROT["Protein"]
        NUC["Nucleic Acid"]
        WAT["Water"]
        ION["Ion"]
        LIG["Ligand"]
    end

    subgraph "Fixed Charges (ffcharge)"
        FF_PROT["ProteinScheme<br><i>AMBER/CHARMM</i>"]
        FF_NUC["NucleicScheme<br><i>AMBER/CHARMM</i>"]
        FF_WAT["WaterScheme<br><i>TIP3P/TIP3P-FB/SPC/OPC3</i>"]
        FF_ION["IonScheme<br><i>Formal charges</i>"]
    end

    subgraph "Ligand QEq"
        LIG_CFG["LigandChargeConfig[]"]
        VAC["Vacuum QEq<br><i>isolated</i>"]
        EMB["Embedded QEq<br><i>polarized by environment</i>"]
    end

    META --> CLASS
    CLASS --> PROT --> FF_PROT
    CLASS --> NUC --> FF_NUC
    CLASS --> WAT --> FF_WAT
    CLASS --> ION --> FF_ION
    CLASS --> LIG --> LIG_CFG
    LIG_CFG --> VAC
    LIG_CFG --> EMB
```

**Atom classification rules:**

| Category     | Source                                    | Charge Source             |
| ------------ | ----------------------------------------- | ------------------------- |
| Protein      | `StandardResidue::ALA..VAL`               | `ffcharge::ProteinScheme` |
| Nucleic Acid | `StandardResidue::A..DI`                  | `ffcharge::NucleicScheme` |
| Water        | `StandardResidue::HOH`                    | `ffcharge::WaterScheme`   |
| Ion          | `ResidueCategory::Ion`                    | `ffcharge::IonScheme`     |
| Ligand       | `ResidueCategory::Hetero` or unrecognized | QEq (vacuum or embedded)  |

**pH-aware terminal handling:**

Terminal protonation states are determined by comparing `BioMetadata.target_ph` with pKa values:

| Terminal   | pKa | Low pH            | High pH             |
| ---------- | --- | ----------------- | ------------------- |
| N-terminal | 8.0 | NH₃⁺ (protonated) | NH₂ (neutral)       |
| C-terminal | 3.1 | COOH (protonated) | COO⁻ (deprotonated) |

#### 3.3.3 Embedded QEq for Ligands

Embedded QEq polarizes the ligand's charge distribution based on the electrostatic potential from surrounding fixed-charge atoms (proteins, nucleic acids).

```mermaid
flowchart TD
    subgraph "Input"
        LIG_ATOMS["Ligand atoms"]
        FIXED["Fixed-charge atoms"]
        CUTOFF["cutoff_radius"]
    end

    subgraph "Spatial Query"
        GRID["SpatialGrid<br><i>cell_size = cutoff</i>"]
        QUERY["query_radius_multi"]
        ENV["Environment atoms"]
    end

    subgraph "cheq Solver"
        EXT["ExternalPotential<br><i>from point charges</i>"]
        SOLVE["solve_in_field"]
    end

    subgraph "Output"
        CHARGES["Polarized ligand charges"]
    end

    FIXED --> GRID
    LIG_ATOMS --> QUERY
    CUTOFF --> QUERY
    GRID --> QUERY --> ENV
    ENV --> EXT --> SOLVE --> CHARGES
```

**Spatial grid optimization:**

The `SpatialGrid` data structure provides O(1) amortized neighbor lookups:

- Space is divided into cubic cells of size equal to the cutoff radius
- Each cell stores indices of atoms within its bounds
- Range queries check the 27 neighboring cells (3×3×3 cube)
- Multi-point queries (`query_radius_multi`) efficiently find all atoms within range of any ligand atom

**Why fixed charges only as environment:**

Only atoms with pre-assigned fixed charges (proteins, nucleic acids, water, ions) are included in the external potential. Including other ligands would create mutual dependencies requiring self-consistent iteration, adding complexity without significant accuracy improvement for typical drug-protein systems

### 3.4 Stage 4: Parameter Generation

The final stage generates all force field parameters from the typed and charged system:

```mermaid
flowchart TD
    subgraph "Input"
        INT["IntermediateSystem<br><i>typed + charged</i>"]
        FFP["ForceFieldParams<br><i>from TOML</i>"]
        CFG["ForgeConfig"]
    end

    subgraph "Atom Processing"
        AT["Collect unique atom types"]
        AP["Generate AtomParam[]<br><i>charge, mass, type_idx</i>"]
    end

    subgraph "Bonded Terms"
        BP["BondPotential[]<br><i>Harmonic or Morse</i>"]
        ANG["AnglePotential[]<br><i>CosineHarmonic or ThetaHarmonic</i>"]
        DIH["DihedralPotential[]<br><i>periodic torsion</i>"]
        IMP["ImproperPotential[]<br><i>Planar or Umbrella</i>"]
    end

    subgraph "Non-Bonded Terms"
        VDW["VdwPairPotential[]<br><i>LJ or Exp-6</i>"]
        HB["HBondPotential[]<br><i>directional H-bonds</i>"]
    end

    subgraph "Output"
        FS["ForgedSystem"]
    end

    INT --> AT --> AP
    INT --> BP
    INT --> ANG
    INT --> DIH
    INT --> IMP
    FFP --> BP
    FFP --> ANG
    FFP --> DIH
    FFP --> IMP
    FFP --> VDW
    FFP --> HB
    CFG --> BP
    CFG --> ANG
    CFG --> VDW

    AP --> FS
    BP --> FS
    ANG --> FS
    DIH --> FS
    IMP --> FS
    VDW --> FS
    HB --> FS
```

---

## 4. Data Structures Reference

### 4.1 Input Types

```mermaid
classDiagram
    class System {
        +atoms: Vec~Atom~
        +bonds: Vec~Bond~
        +box_vectors: Option~[[f64; 3]; 3]~
        +bio_metadata: Option~BioMetadata~
        +atom_count() usize
        +bond_count() usize
        +is_periodic() bool
    }

    class Atom {
        +element: Element
        +position: [f64; 3]
        +new(element, position) Atom
    }

    class Bond {
        +i: usize
        +j: usize
        +order: BondOrder
        +new(i, j, order) Bond
    }

    class BioMetadata {
        +atom_info: Vec~AtomResidueInfo~
        +target_ph: Option~f64~
        +effective_ph() f64
    }

    class AtomResidueInfo {
        +atom_name: String
        +residue_name: String
        +residue_id: i32
        +chain_id: String
        +insertion_code: Option~char~
        +standard_name: Option~StandardResidue~
        +category: ResidueCategory
        +position: ResiduePosition
    }

    System "1" *-- "*" Atom
    System "1" *-- "*" Bond
    System "1" *-- "0..1" BioMetadata
    BioMetadata "1" *-- "*" AtomResidueInfo
```

### 4.2 Output Types

```mermaid
classDiagram
    class ForgedSystem {
        +system: System
        +atom_types: Vec~String~
        +atom_properties: Vec~AtomParam~
        +potentials: Potentials
    }

    class AtomParam {
        +charge: f64
        +mass: f64
        +type_idx: usize
    }

    class Potentials {
        +bonds: Vec~BondPotential~
        +angles: Vec~AnglePotential~
        +dihedrals: Vec~DihedralPotential~
        +impropers: Vec~ImproperPotential~
        +vdw_pairs: Vec~VdwPairPotential~
        +h_bonds: Vec~HBondPotential~
    }

    ForgedSystem "1" *-- "1" System
    ForgedSystem "1" *-- "*" AtomParam
    ForgedSystem "1" *-- "1" Potentials
```

### 4.3 Potential Types

```mermaid
classDiagram
    class BondPotential {
        <<enumeration>>
        Harmonic(i, j, k_force, r0)
        Morse(i, j, r0, d0, alpha)
    }

    class AnglePotential {
        <<enumeration>>
        CosineHarmonic(i, j, k, k_force, theta0)
        ThetaHarmonic(i, j, k, k_force, theta0)
    }

    class DihedralPotential {
        +i: usize
        +j: usize
        +k: usize
        +l: usize
        +v_barrier: f64
        +periodicity: i32
        +phase_offset: f64
    }

    class ImproperPotential {
        <<enumeration>>
        Planar(i, j, k, l, k_force, chi0)
        Umbrella(center, p1, p2, p3, k_force, psi0)
    }

    class VdwPairPotential {
        <<enumeration>>
        LennardJones(type1_idx, type2_idx, sigma, epsilon)
        Exponential6(type1_idx, type2_idx, a, b, c)
    }

    class HBondPotential {
        +donor_type_idx: usize
        +hydrogen_type_idx: usize
        +acceptor_type_idx: usize
        +d0: f64
        +r0: f64
    }
```

---

## 5. Algorithm Deep Dives

### 5.1 DREIDING Atom Typing Rules

The typing engine uses a priority-based rule system. Key atom types:

| Atom Type | Description        | Hybridization   | Priority |
| --------- | ------------------ | --------------- | -------- |
| `C_3`     | sp³ carbon         | SP3             | 100      |
| `C_2`     | sp² carbon         | SP2             | 200      |
| `C_R`     | Resonant carbon    | Resonant        | 400      |
| `C_1`     | sp carbon          | SP              | 300      |
| `N_3`     | sp³ nitrogen       | SP3             | 100      |
| `N_R`     | Resonant nitrogen  | Resonant        | 400      |
| `O_3`     | sp³ oxygen         | SP3             | 100      |
| `O_2`     | sp² oxygen         | SP2             | 200      |
| `O_R`     | Resonant oxygen    | Resonant        | 400      |
| `H_`      | Generic hydrogen   | —               | 1        |
| `H_HB`    | H-bonding hydrogen | neighbor O or N | 78-80    |

**Fixed-point iteration:**

```pseudo
atom_states[*] = (None, 0)

repeat:
    changed = false
    for each atom i:
        for each rule r (sorted by priority desc):
            if matches(r.conditions, atom[i], atom_states):
                if r.priority > atom_states[i].priority:
                    atom_states[i] = (r.type, r.priority)
                    changed = true
                break
until not changed or rounds > 100
```

### 5.2 QEq Charge Equilibration

The QEq method solves for partial charges that equalize the chemical potential:

**Electronegativity equalization principle:**
$$\chi_i + J_i q_i + \sum_{j \neq i} J_{ij}(r_{ij}) q_j = \mu \quad \forall i$$

**Screened Coulomb interaction:**
$$J_{ij}(r_{ij}) = \frac{14.4}{\sqrt{r_{ij}^2 + \frac{1}{\lambda^2}(r_i r_j)^{n_i + n_j}}}$$

### 5.3 Bond Potential Calculation

**Equilibrium bond length:**
$$r_0 = R_i + R_j - \delta$$

where $R_i$, $R_j$ are covalent radii and $\delta = 0.01$ Å.

**Force constant scaling by bond order:**

| Bond Order | Multiplier |
| ---------- | ---------- |
| Single     | 1.0        |
| Resonant   | 1.5        |
| Double     | 2.0        |
| Triple     | 3.0        |

**Harmonic potential:**
$$V_{bond}(r) = \frac{1}{2} k_b \cdot n \cdot (r - r_0)^2$$

**Morse potential:**
$$V_{bond}(r) = D_0 \cdot n \cdot \left[1 - e^{-\alpha(r - r_0)}\right]^2$$

where $\alpha = \sqrt{k_b / (2 D_0)}$.

### 5.4 Torsion Parameter Rules

Dihedral parameters depend on the hybridization of the central bond atoms:

```mermaid
flowchart TD
    subgraph "j-k Hybridization"
        SP3_SP3["SP3 — SP3"]
        SP2_SP2["SP2 — SP2"]
        RES_RES["Resonant — Resonant"]
        SP2_RES["SP2 — Resonant"]
        MIXED["SP2/R — SP3"]
    end

    subgraph "Parameters"
        P1["V=2.0, n=3, φ=180°"]
        P2["V=45.0, n=2, φ=180°"]
        P3["V=25.0, n=2, φ=180°"]
        P4["V=5.0, n=2, φ=180°"]
        P5["V=2.0, n=2/3/6, varies"]
    end

    SP3_SP3 --> P1
    SP2_SP2 --> P2
    RES_RES --> P3
    SP2_RES --> P4
    MIXED --> P5
```

**Special case:** When both central atoms are in the oxygen column (O, S, Se, Te), SP3-SP3 torsions use $n=2$, $\phi=90°$.

### 5.5 Van der Waals Mixing Rules

**Lennard-Jones 12-6:**
$$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}$$
$$\epsilon_{ij} = \sqrt{\epsilon_i \cdot \epsilon_j}$$

**Exponential-6:**
$$A_{ij} = \sqrt{A_i \cdot A_j}$$
$$B_{ij} = \frac{B_i + B_j}{2}$$
$$C_{ij} = \sqrt{C_i \cdot C_j}$$

### 5.6 Hydrogen Bond Detection

H-bond potentials are generated when:

1. Hydrogen has type `H_HB` (bonded to O or N)
2. An acceptor atom element is O, N, F
3. The acceptor is not the hydrogen's bonded atom

---

## 6. Error Handling Strategy

### 6.1 Error Types

```mermaid
flowchart TD
    subgraph "forge::Error"
        FE_PP["ParameterParse<br><i>TOML parsing failure</i>"]
        FE_RP["RuleParse<br><i>typing rules malformed</i>"]
        FE_AT["AtomTyping<br><i>dreid-typer failure</i>"]
        FE_CC["ChargeCalculation<br><i>cheq failure</i>"]
        FE_MB["MissingBioMetadata<br><i>hybrid requires metadata</i>"]
        FE_HC["HybridChargeAssignment<br><i>classical charge lookup failed</i>"]
        FE_MP["MissingParameter<br><i>no params for atom type</i>"]
        FE_IB["InvalidBond<br><i>bond index out of range</i>"]
        FE_ES["EmptySystem<br><i>no atoms</i>"]
        FE_CV["Conversion<br><i>internal consistency</i>"]
    end

    subgraph "io::Error"
        IE_IO["Io<br><i>file system errors</i>"]
        IE_PA["Parse<br><i>format + line + details</i>"]
        IE_UR["UnsupportedReadFormat"]
        IE_UW["UnsupportedWriteFormat"]
        IE_MM["MissingMetadata<br><i>bio_metadata required</i>"]
        IE_BF["BioForgePreparation"]
        IE_CV["Conversion"]
    end
```

### 6.2 Error Propagation

- **Eager validation:** Bond indices are checked immediately in `IntermediateSystem::from_system`
- **Typed errors:** Each error variant carries context (atom type, line number, format)
- **thiserror integration:** All errors implement `std::error::Error` with proper `Display`
- **No panics:** Library functions return `Result` types; panics only occur for internal bugs

---

## 7. Configuration Reference

### 7.1 ForgeConfig

```rust
pub struct ForgeConfig {
    pub rules: Option<String>,        // Custom typing rules (TOML)
    pub params: Option<String>,       // Custom FF params (TOML)
    pub charge_method: ChargeMethod,  // None, Qeq, or Hybrid
    pub bond_potential: BondPotentialType,   // Harmonic or Morse
    pub angle_potential: AnglePotentialType, // ThetaHarmonic or CosineHarmonic
    pub vdw_potential: VdwPotentialType,     // LennardJones or Exponential6
}
```

### 7.2 Potential Type Selection

| Config Field      | Options                           | DREIDING Default |
| ----------------- | --------------------------------- | ---------------- |
| `bond_potential`  | `Harmonic`, `Morse`               | `Harmonic`       |
| `angle_potential` | `ThetaHarmonic`, `CosineHarmonic` | `ThetaHarmonic`  |
| `vdw_potential`   | `LennardJones`, `Exponential6`    | `LennardJones`   |

### 7.3 Charge Method Configuration

```rust
pub enum ChargeMethod {
    None,              // All charges = 0.0
    Qeq(QeqConfig),    // Global QEq for all atoms
    Hybrid(HybridConfig), // Classical + QEq (requires BioMetadata)
}
```

### 7.4 QeqConfig

```rust
pub struct QeqConfig {
    pub total_charge: f64,       // Target system charge (default: 0.0)
    pub solver_options: SolverOptions,
}
```

### 7.5 HybridConfig

```rust
pub struct HybridConfig {
    pub protein_scheme: ProteinScheme,   // AMBER ff99SB/ff14SB/ff19SB, AMBER ff03, CHARMM C22/C27/C36/C36m
    pub nucleic_scheme: NucleicScheme,   // AMBER OL15/OL21/OL24/bsc1/OL3, CHARMM C27/C36
    pub water_scheme: WaterScheme,       // TIP3P, TIP3P-FB, SPC, SPC/E, OPC3
    pub ligand_configs: Vec<LigandChargeConfig>, // Per-ligand QEq configuration
    pub default_ligand_qeq: QeqConfig,   // Default QEq for unlisted ligands
}
```

### 7.6 Ligand Charge Configuration

```rust
pub struct LigandChargeConfig {
    pub selector: ResidueSelector,  // Target residue (chain_id, residue_id, insertion_code)
    pub method: LigandQeqMethod,    // Vacuum or Embedded QEq
}

pub enum LigandQeqMethod {
    Vacuum(QeqConfig),        // Isolated QEq calculation
    Embedded(EmbeddedQeqConfig), // QEq polarized by environment
}

pub struct EmbeddedQeqConfig {
    pub cutoff_radius: f64,   // Environment search radius in Å (default: 10.0)
    pub qeq: QeqConfig,       // QEq solver settings
}
```

### 7.7 Residue Selector

```rust
pub struct ResidueSelector {
    pub chain_id: String,
    pub residue_id: i32,
    pub insertion_code: Option<char>,  // None matches any insertion code
}
```

### 7.8 I/O Configurations

| Config              | Purpose                            |
| ------------------- | ---------------------------------- |
| `CleanConfig`       | Control water/ion/hydrogen removal |
| `ProtonationConfig` | pH, histidine strategy             |
| `SolvateConfig`     | Water box margin, ion types        |
| `TopologyConfig`    | Hetero templates, disulfide cutoff |
| `LammpsConfig`      | Cutoffs, boundary conditions       |

---

## 8. Dependency Architecture

```mermaid
flowchart BT
    subgraph "External Crates"
        BF["bio-forge v0.3<br><i>structure preparation</i>"]
        DT["dreid-typer v0.4<br><i>atom typing</i>"]
        CQ["cheq v0.5<br><i>QEq charges</i>"]
        FF["ffcharge v0.2<br><i>classical FF charges</i>"]
    end

    subgraph "dreid-forge Modules"
        IO["io<br><i>readers/writers</i>"]
        MODEL["model<br><i>data structures</i>"]
        FORGE["forge<br><i>parameterization</i>"]
    end

    subgraph "Shared Dependencies"
        SERDE["serde + toml"]
        THISERROR["thiserror"]
    end

    BF --> IO
    DT --> FORGE
    CQ --> FORGE
    FF --> FORGE
    SERDE --> BF
    SERDE --> DT
    SERDE --> CQ
    THISERROR --> IO
    THISERROR --> FORGE
```

### Integration Points

| Crate         | Version | Integration Module                            | Key Types Exchanged                                          |
| ------------- | ------- | --------------------------------------------- | ------------------------------------------------------------ |
| `bio-forge`   | 0.3     | `io::util`, `io::pdb`, `io::mmcif`            | `Structure`, `Topology`, `Template`                          |
| `dreid-typer` | 0.4     | `forge::typer`                                | `MolecularGraph`, `MolecularTopology`, `Hybridization`       |
| `cheq`        | 0.5     | `forge::charge::qeq`, `forge::charge::hybrid` | `QEqSolver`, `ExternalPotential`, `PointCharge`              |
| `ffcharge`    | 0.2     | `forge::charge::hybrid`                       | `ProteinScheme`, `NucleicScheme`, `WaterScheme`, `IonScheme` |

---

## 9. Performance Considerations

### 9.1 Computational Complexity

| Stage                | Complexity            | Dominant Operation         |
| -------------------- | --------------------- | -------------------------- |
| System conversion    | O(A + B)              | Neighbor list construction |
| Ring detection       | O(B × R)              | SSSR enumeration           |
| Aromaticity          | O(R)                  | π-electron counting        |
| Typing               | O(A × Rules × Rounds) | Rule matching              |
| QEq                  | O(A³)                 | Matrix solve               |
| Parameter generation | O(B + Ang + Dih)      | Term enumeration           |

Where A = atoms, B = bonds, R = rings.

### 9.2 Memory Layout

- **IntermediateSystem:** Owned data, heap-allocated vectors
- **ForgedSystem:** Clones input System, owns all generated parameters
- **Potentials:** Flat vectors with no internal pointers

### 9.3 Parallelization Opportunities

- **cheq:** Coulomb matrix construction is parallelized via `rayon`
- **dreid-typer:** Perception passes are sequential but O(A)
- **Parameter generation:** Each potential type can be generated independently
