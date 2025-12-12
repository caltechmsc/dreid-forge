//! Core data structures representing molecular systems and force field topologies.
//!
//! This module provides the foundational types that flow through `dreid-forge`:
//!
//! - [`atom`] – Minimal atom representation with element and Cartesian coordinates.
//! - [`types`] – Periodic table elements and bond order classifications.
//! - [`system`] – Complete molecular systems with atoms, bonds, and optional box vectors.
//! - [`metadata`] – Biological context (residue names, chain IDs) for macromolecular structures.
//! - [`topology`] – Force field output including potentials, charges, and atom type assignments.
//!
//! The data model intentionally separates raw molecular geometry ([`System`]) from
//! parameterized topology ([`ForgedSystem`]), allowing the [`crate::forge`] pipeline
//! to transform one into the other while preserving provenance.
//!
//! [`System`]: system::System
//! [`ForgedSystem`]: topology::ForgedSystem

pub mod atom;
pub mod metadata;
pub mod system;
pub mod topology;
pub mod types;
