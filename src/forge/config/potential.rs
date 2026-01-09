//! Potential energy function type configurations.
//!
//! This module defines the available (configurable) functional forms
//! of potential energy term in the DREIDING force field.

/// Bond stretching potential function type.
///
/// Determines the functional form used for bond stretching terms
/// in the force field. Both options use the same equilibrium bond
/// length but differ in behavior far from equilibrium.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondPotentialType {
    /// Harmonic bond potential.
    #[default]
    Harmonic,

    /// Morse anharmonic potential.
    Morse,
}

/// Angle bending potential function type.
///
/// Determines the functional form used for angle bending terms.
/// The choice affects behavior especially for linear or near-linear
/// angles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AnglePotentialType {
    /// Cosine-harmonic angle potential (DREIDING original).
    CosineHarmonic,

    /// Theta-harmonic angle potential.
    #[default]
    ThetaHarmonic,
}

/// Van der Waals non-bonded potential function type.
///
/// Determines the functional form used for van der Waals (dispersion
/// and repulsion) interactions between non-bonded atom pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VdwPotentialType {
    /// Lennard-Jones 12-6 potential.
    #[default]
    LennardJones,

    /// Exponential-6 (Buckingham) potential.
    Exponential6,
}
