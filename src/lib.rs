//! A high-performance, pure Rust library for automated DREIDING force field parameterization.
//! It orchestrates structure repair, topology perception, and partial atomic charge calculation
//! to produce simulation-ready inputs for both biological macromolecules and arbitrary chemical systems.
//!
//! # Features
//!
//! - **Atom typing** — Automatic assignment of DREIDING atom types based on
//!   element, hybridization, and local bonding environment
//! - **Charge calculation** — QEq charge equilibration for partial atomic charges
//! - **Parameter generation** — Bond, angle, dihedral, improper, van der Waals,
//!   and hydrogen bond potentials with multiple functional forms
//! - **Flexible I/O** — Read/write PDB, mmCIF, MOL2, SDF formats; export to
//!   BGF and LAMMPS simulation input files
//!
//! # Quick Start
//!
//! The main entry point is the [`forge`] function, which takes a [`System`] and
//! [`ForgeConfig`] and produces a fully parameterized [`ForgedSystem`]:
//!
//! ```
//! use dreid_forge::{System, Atom, Bond, Element, BondOrder};
//! use dreid_forge::{forge, ForgeConfig, ForgeError};
//!
//! // Build an ethanol molecule (C₂H₅OH)
//! let mut system = System::new();
//!
//! // Atoms: C1 (methyl), C2 (methylene), O, and hydrogens
//! system.atoms.push(Atom::new(Element::C, [-1.270,  0.248,  0.000])); // C1
//! system.atoms.push(Atom::new(Element::C, [ 0.139, -0.308,  0.000])); // C2
//! system.atoms.push(Atom::new(Element::O, [ 1.036,  0.789,  0.000])); // O
//! system.atoms.push(Atom::new(Element::H, [-1.317,  0.885,  0.883])); // H on C1
//! system.atoms.push(Atom::new(Element::H, [-1.317,  0.885, -0.883])); // H on C1
//! system.atoms.push(Atom::new(Element::H, [-2.030, -0.533,  0.000])); // H on C1
//! system.atoms.push(Atom::new(Element::H, [ 0.358, -0.920,  0.876])); // H on C2
//! system.atoms.push(Atom::new(Element::H, [ 0.358, -0.920, -0.876])); // H on C2
//! system.atoms.push(Atom::new(Element::H, [ 1.939,  0.473,  0.000])); // H on O
//!
//! // Bonds
//! system.bonds.push(Bond::new(0, 1, BondOrder::Single)); // C1-C2
//! system.bonds.push(Bond::new(1, 2, BondOrder::Single)); // C2-O
//! system.bonds.push(Bond::new(0, 3, BondOrder::Single)); // C1-H
//! system.bonds.push(Bond::new(0, 4, BondOrder::Single)); // C1-H
//! system.bonds.push(Bond::new(0, 5, BondOrder::Single)); // C1-H
//! system.bonds.push(Bond::new(1, 6, BondOrder::Single)); // C2-H
//! system.bonds.push(Bond::new(1, 7, BondOrder::Single)); // C2-H
//! system.bonds.push(Bond::new(2, 8, BondOrder::Single)); // O-H
//!
//! // Parameterize with default settings
//! let forged = forge(&system, &ForgeConfig::default())?;
//!
//! // Atom types: C_3 (sp³ carbon), O_3 (sp³ oxygen), H_, H_HB (H-bond donor)
//! assert_eq!(forged.atom_types.len(), 4);
//! assert!(forged.atom_types.contains(&"C_3".to_string()));
//! assert!(forged.atom_types.contains(&"O_3".to_string()));
//! assert!(forged.atom_types.contains(&"H_HB".to_string()));
//!
//! // Per-atom properties: charge, mass, type index
//! assert_eq!(forged.atom_properties.len(), 9);
//!
//! // Bond potentials: C-C, C-O, C-H, O-H
//! assert_eq!(forged.potentials.bonds.len(), 8);
//!
//! // Angle potentials: all angles around sp³ centers
//! assert_eq!(forged.potentials.angles.len(), 13);
//!
//! // Dihedral potentials: H-C-C-H, H-C-C-O, C-C-O-H, H-C-O-H
//! assert_eq!(forged.potentials.dihedrals.len(), 12);
//!
//! // Improper potentials: none for all-sp³ molecule
//! assert!(forged.potentials.impropers.is_empty());
//!
//! // VdW pair potentials: n(n+1)/2 = 4×5/2 = 10 pairs
//! assert_eq!(forged.potentials.vdw_pairs.len(), 10);
//!
//! // H-bond potential: O_3 as donor/acceptor with H_HB
//! assert_eq!(forged.potentials.h_bonds.len(), 1);
//! # Ok::<(), ForgeError>(())
//! ```
//!
//! # Module Organization
//!
//! - [`io`] — File I/O for molecular structures (PDB, mmCIF, MOL2, SDF, BGF, LAMMPS)
//! - [`forge`] — Main parameterization function
//! - [`ForgeConfig`] — Configuration for potential types and charge methods
//!
//! # Data Types
//!
//! ## Input Structures
//!
//! - [`System`] — Molecular system with atoms, bonds, and optional metadata
//! - [`Atom`] — Single atom with element and Cartesian coordinates
//! - [`Bond`] — Bond between two atoms with bond order
//! - [`Element`] — Chemical element (H through Og)
//! - [`BondOrder`] — Bond order (Single, Double, Triple, Aromatic)
//!
//! ## Output Structures
//!
//! - [`ForgedSystem`] — Fully parameterized system ready for simulation
//! - [`AtomParam`] — Per-atom charge, mass, and type index
//! - [`Potentials`] — Collection of all potential energy functions
//! - [`BondPotential`] — Harmonic or Morse bond stretching
//! - [`AnglePotential`] — Cosine-harmonic or theta-harmonic bending
//! - [`DihedralPotential`] — Periodic torsion potentials
//! - [`ImproperPotential`] — Planar or umbrella out-of-plane terms
//! - [`VdwPairPotential`] — Lennard-Jones or Exponential-6 dispersion
//! - [`HBondPotential`] — Directional hydrogen bond terms
//!
//! ## Configuration
//!
//! - [`ChargeMethod`] — None, QEq, or Hybrid charge equilibration
//! - [`QeqConfig`] — QEq solver settings (total charge, convergence)
//! - [`HybridConfig`] — Hybrid biological/QEq charge assignment settings
//! - [`BondPotentialType`] — Harmonic vs Morse selection
//! - [`AnglePotentialType`] — Cosine-harmonic vs theta-harmonic
//! - [`VdwPotentialType`] — Lennard-Jones vs Exponential-6
//!
//! ## Biological Metadata
//!
//! - [`BioMetadata`] — Per-atom PDB/mmCIF annotations
//! - [`AtomResidueInfo`] — Residue name, chain, sequence number
//! - [`StandardResidue`] — Standard amino acids and nucleotides
//! - [`ResidueCategory`] — Standard, hetero, or ion classification
//! - [`ResiduePosition`] — Terminal position (N/C-terminal, 5'/3'-end)

mod forge;
mod model;

pub mod io;

pub use model::atom::Atom;
pub use model::system::{Bond, System};
pub use model::types::{BondOrder, Element, ParseBondOrderError, ParseElementError};

pub use model::topology::{
    AnglePotential, AtomParam, BondPotential, DihedralPotential, ForgedSystem, HBondPotential,
    ImproperPotential, Potentials, VdwPairPotential,
};

pub use model::metadata::{
    AtomResidueBuilder, AtomResidueInfo, BioMetadata, ResidueCategory, ResiduePosition,
    StandardResidue,
};

pub use forge::{
    AnglePotentialType, BondPotentialType, ChargeMethod, EmbeddedQeqConfig, ForgeConfig,
    HybridConfig, LigandChargeConfig, LigandQeqMethod, NucleicScheme, ProteinScheme, QeqConfig,
    ResidueSelector, SolverOptions, VdwPotentialType, WaterScheme, forge,
};

pub use forge::Error as ForgeError;
