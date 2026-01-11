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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_potential_default() {
        assert_eq!(BondPotentialType::default(), BondPotentialType::Harmonic);
    }

    #[test]
    fn angle_potential_default() {
        assert_eq!(
            AnglePotentialType::default(),
            AnglePotentialType::ThetaHarmonic
        );
    }

    #[test]
    fn vdw_potential_default() {
        assert_eq!(VdwPotentialType::default(), VdwPotentialType::LennardJones);
    }

    #[test]
    fn potential_types_are_copy() {
        let bond = BondPotentialType::Morse;
        let bond_copy = bond;
        assert_eq!(bond, bond_copy);

        let angle = AnglePotentialType::CosineHarmonic;
        let angle_copy = angle;
        assert_eq!(angle, angle_copy);

        let vdw = VdwPotentialType::Exponential6;
        let vdw_copy = vdw;
        assert_eq!(vdw, vdw_copy);
    }
}
