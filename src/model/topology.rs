//! Force field topology and potential function definitions.
//!
//! This module defines the data structures for DREIDING force field
//! parameters and potential energy functions. After parameterization,
//! a [`ForgedSystem`] contains all information needed to run molecular
//! dynamics or energy minimization simulations.
//!
//! # Potential Functions
//!
//! The DREIDING force field supports several potential function types:
//!
//! - **Bonds**: [`BondPotential`] — Harmonic or Morse stretching
//! - **Angles**: [`AnglePotential`] — Cosine-harmonic, cosine-linear, or theta-harmonic bending
//! - **Torsions**: [`TorsionPotential`] — Periodic torsion potentials
//! - **Inversions**: [`InversionPotential`] — Out-of-plane (planar/umbrella)
//! - **Van der Waals**: [`VdwPairPotential`] — LJ 12-6 or Buckingham Exp-6
//! - **Hydrogen bonds**: [`HBondPotential`] — Directional H-bond terms
//!
//! # Output Structure
//!
//! The [`ForgedSystem`] struct combines the original molecular system
//! with computed atom types, partial charges, and all bonded/non-bonded
//! potential parameters needed for simulation.

use super::system::System;

/// Per-atom force field parameters.
///
/// Contains the computed properties for a single atom after
/// DREIDING parameterization.
///
/// # Fields
///
/// * `charge` — Partial atomic charge in elementary charge units
/// * `mass` — Atomic mass in atomic mass units (amu)
/// * `type_idx` — Index into the atom type table
#[derive(Debug, Clone, PartialEq)]
pub struct AtomParam {
    /// Partial atomic charge (e).
    pub charge: f64,
    /// Atomic mass (amu).
    pub mass: f64,
    /// Index into the atom type name table.
    pub type_idx: usize,
}

/// Bond stretching potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum BondPotential {
    /// Harmonic bond stretching potential.
    Harmonic {
        /// Atom indices (i, j).
        atoms: (usize, usize),
        /// Half force constant (kcal/mol/Å²).
        k_half: f64,
        /// Equilibrium bond length (Å).
        r0: f64,
    },
    /// Morse anharmonic bond potential.
    Morse {
        /// Atom indices (i, j).
        atoms: (usize, usize),
        /// Dissociation energy (kcal/mol).
        de: f64,
        /// Equilibrium bond length (Å).
        r0: f64,
        /// Stiffness parameter alpha (Å⁻¹).
        alpha: f64,
    },
}

/// Angle bending potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum AnglePotential {
    /// Cosine-harmonic angle potential (for θ₀ ≠ 180°).
    CosineHarmonic {
        /// Atom indices (i, j, k) where j is the central atom.
        atoms: (usize, usize, usize),
        /// Half force constant (kcal/mol/rad²).
        k_half: f64,
        /// Cosine of equilibrium angle.
        cos0: f64,
    },
    /// Cosine-linear angle potential for linear geometries (θ₀ = 180°).
    CosineLinear {
        /// Atom indices (i, j, k) where j is the central atom.
        atoms: (usize, usize, usize),
        /// Force constant (kcal/mol/rad²).
        k: f64,
    },
    /// Simple theta-harmonic angle potential.
    ThetaHarmonic {
        /// Atom indices (i, j, k) where j is the central atom.
        atoms: (usize, usize, usize),
        /// Half force constant (kcal/mol/rad²).
        k_half: f64,
        /// Equilibrium angle (radians).
        theta0: f64,
    },
}

/// Torsion (proper dihedral) potential.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionPotential {
    /// Atom indices (i, j, k, l) where j-k is the central bond.
    pub atoms: (usize, usize, usize, usize),
    /// Half barrier height (kcal/mol).
    pub v_half: f64,
    /// Periodicity (number of minima in 360°).
    pub n: u8,
    /// Precomputed cos(n·φ₀).
    pub cos_n_phi0: f64,
    /// Precomputed sin(n·φ₀).
    pub sin_n_phi0: f64,
}

/// Inversion (improper dihedral) potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum InversionPotential {
    /// Planar inversion for sp² centers.
    Planar {
        /// Atom indices (center, axis, plane1, plane2) where center is the sp² atom,
        /// axis is the out-of-plane neighbor, and plane1/plane2 are in-plane neighbors.
        atoms: (usize, usize, usize, usize),
        /// Half force constant (kcal/mol/rad²).
        c_half: f64,
    },
    /// Umbrella inversion for pyramidal centers.
    Umbrella {
        /// Atom indices (center, axis, plane1, plane2) where center is the pyramidal atom,
        /// axis is one neighbor defining the inversion axis, and plane1/plane2 are the other neighbors.
        atoms: (usize, usize, usize, usize),
        /// Half force constant (kcal/mol/rad²).
        c_half: f64,
        /// Cosine of equilibrium inversion angle.
        cos_psi0: f64,
    },
}

/// Van der Waals non-bonded pair potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum VdwPairPotential {
    /// Lennard-Jones 12-6 potential.
    LennardJones {
        /// First atom type index.
        type1_idx: usize,
        /// Second atom type index.
        type2_idx: usize,
        /// Energy well depth (kcal/mol).
        d0: f64,
        /// Squared equilibrium distance (Å²).
        r0_sq: f64,
    },
    /// Buckingham (Exponential-6) potential with energy reflection.
    Buckingham {
        /// First atom type index.
        type1_idx: usize,
        /// Second atom type index.
        type2_idx: usize,
        /// Repulsion prefactor (kcal/mol).
        a: f64,
        /// Repulsion decay (Å⁻¹).
        b: f64,
        /// Attraction coefficient (kcal·Å⁶/mol).
        c: f64,
        /// Squared distance of energy maximum (Å²).
        r_max_sq: f64,
        /// Twice the energy at maximum (kcal/mol).
        two_e_max: f64,
    },
}

/// Hydrogen bond directional potential.
#[derive(Debug, Clone, PartialEq)]
pub struct HBondPotential {
    /// Donor atom type index (D in D-H···A).
    pub donor_type_idx: usize,
    /// Hydrogen atom type index (H in D-H···A).
    pub hydrogen_type_idx: usize,
    /// Acceptor atom type index (A in D-H···A).
    pub acceptor_type_idx: usize,
    /// Energy well depth (kcal/mol).
    pub d_hb: f64,
    /// Squared equilibrium distance (Å²).
    pub r_hb_sq: f64,
}

/// Collection of all potential energy functions for a system.
///
/// Groups all bonded and non-bonded interaction parameters
/// computed during DREIDING parameterization.
#[derive(Debug, Clone, Default)]
pub struct Potentials {
    /// Bond stretching potentials.
    pub bonds: Vec<BondPotential>,
    /// Angle bending potentials.
    pub angles: Vec<AnglePotential>,
    /// Torsion (proper dihedral) potentials.
    pub torsions: Vec<TorsionPotential>,
    /// Inversion (improper dihedral) potentials.
    pub inversions: Vec<InversionPotential>,
    /// Van der Waals pair potentials between atom types.
    pub vdw_pairs: Vec<VdwPairPotential>,
    /// Hydrogen bond potentials.
    pub h_bonds: Vec<HBondPotential>,
}

/// A fully parameterized molecular system.
///
/// Contains the original [`System`] along with computed DREIDING
/// force field parameters including atom types, partial charges,
/// and all potential energy function parameters.
///
/// This is the primary output of the [`forge`](crate::forge) function
/// and contains everything needed to write simulation input files
/// for molecular dynamics packages.
///
/// # Fields
///
/// * `system` — Original molecular structure
/// * `atom_types` — DREIDING atom type names (e.g., "C_3", "O_2")
/// * `atom_properties` — Per-atom charges, masses, and type indices
/// * `potentials` — All bonded and non-bonded potential parameters
#[derive(Debug, Clone)]
pub struct ForgedSystem {
    /// The original molecular system with atoms and bonds.
    pub system: System,
    /// DREIDING atom type names, indexed by type_idx.
    pub atom_types: Vec<String>,
    /// Per-atom force field parameters.
    pub atom_properties: Vec<AtomParam>,
    /// All potential energy function parameters.
    pub potentials: Potentials,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_harmonic_construction() {
        let bond = BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 350.0,
            r0: 1.54,
        };

        match bond {
            BondPotential::Harmonic { atoms, k_half, r0 } => {
                assert_eq!(atoms.0, 0);
                assert_eq!(atoms.1, 1);
                assert_eq!(k_half, 350.0);
                assert_eq!(r0, 1.54);
            }
            _ => panic!("expected Harmonic variant"),
        }
    }

    #[test]
    fn bond_morse_construction() {
        let bond = BondPotential::Morse {
            atoms: (2, 3),
            de: 70.0,
            r0: 1.10,
            alpha: 2.0,
        };

        match bond {
            BondPotential::Morse {
                atoms,
                de,
                r0,
                alpha,
            } => {
                assert_eq!(atoms.0, 2);
                assert_eq!(atoms.1, 3);
                assert_eq!(de, 70.0);
                assert_eq!(r0, 1.10);
                assert_eq!(alpha, 2.0);
            }
            _ => panic!("expected Morse variant"),
        }
    }

    #[test]
    fn bond_partial_eq() {
        let b1 = BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 350.0,
            r0: 1.54,
        };
        let b2 = BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 350.0,
            r0: 1.54,
        };
        let b3 = BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 400.0,
            r0: 1.54,
        };

        assert_eq!(b1, b2);
        assert_ne!(b1, b3);
    }

    #[test]
    fn bond_clone() {
        let bond = BondPotential::Harmonic {
            atoms: (5, 6),
            k_half: 350.0,
            r0: 1.09,
        };
        let cloned = bond.clone();
        assert_eq!(bond, cloned);
    }

    #[test]
    fn angle_cosine_harmonic_construction() {
        let angle = AnglePotential::CosineHarmonic {
            atoms: (0, 1, 2),
            k_half: 50.0,
            cos0: 0.5,
        };

        match angle {
            AnglePotential::CosineHarmonic {
                atoms,
                k_half,
                cos0,
            } => {
                assert_eq!(atoms.0, 0);
                assert_eq!(atoms.1, 1);
                assert_eq!(atoms.2, 2);
                assert_eq!(k_half, 50.0);
                assert_eq!(cos0, 0.5);
            }
            _ => panic!("expected CosineHarmonic variant"),
        }
    }

    #[test]
    fn angle_cosine_linear_construction() {
        let angle = AnglePotential::CosineLinear {
            atoms: (6, 7, 8),
            k: 100.0,
        };

        match angle {
            AnglePotential::CosineLinear { atoms, k } => {
                assert_eq!(atoms.0, 6);
                assert_eq!(atoms.1, 7);
                assert_eq!(atoms.2, 8);
                assert_eq!(k, 100.0);
            }
            _ => panic!("expected CosineLinear variant"),
        }
    }

    #[test]
    fn angle_theta_harmonic_construction() {
        let angle = AnglePotential::ThetaHarmonic {
            atoms: (3, 4, 5),
            k_half: 50.0,
            theta0: 1.911,
        };

        match angle {
            AnglePotential::ThetaHarmonic {
                atoms,
                k_half,
                theta0,
            } => {
                assert_eq!(atoms.0, 3);
                assert_eq!(atoms.1, 4);
                assert_eq!(atoms.2, 5);
                assert_eq!(k_half, 50.0);
                assert_eq!(theta0, 1.911);
            }
            _ => panic!("expected ThetaHarmonic variant"),
        }
    }

    #[test]
    fn angle_partial_eq() {
        let a1 = AnglePotential::CosineHarmonic {
            atoms: (0, 1, 2),
            k_half: 50.0,
            cos0: 0.5,
        };
        let a2 = AnglePotential::CosineHarmonic {
            atoms: (0, 1, 2),
            k_half: 50.0,
            cos0: 0.5,
        };
        assert_eq!(a1, a2);
    }

    #[test]
    fn torsion_construction() {
        let torsion = TorsionPotential {
            atoms: (0, 1, 2, 3),
            v_half: 1.0,
            n: 3,
            cos_n_phi0: -1.0,
            sin_n_phi0: 0.0,
        };

        assert_eq!(torsion.atoms.0, 0);
        assert_eq!(torsion.atoms.1, 1);
        assert_eq!(torsion.atoms.2, 2);
        assert_eq!(torsion.atoms.3, 3);
        assert_eq!(torsion.v_half, 1.0);
        assert_eq!(torsion.n, 3);
        assert_eq!(torsion.cos_n_phi0, -1.0);
        assert_eq!(torsion.sin_n_phi0, 0.0);
    }

    #[test]
    fn torsion_partial_eq() {
        let t1 = TorsionPotential {
            atoms: (0, 1, 2, 3),
            v_half: 1.0,
            n: 3,
            cos_n_phi0: -1.0,
            sin_n_phi0: 0.0,
        };
        let t2 = TorsionPotential {
            atoms: (0, 1, 2, 3),
            v_half: 1.0,
            n: 3,
            cos_n_phi0: -1.0,
            sin_n_phi0: 0.0,
        };
        let t3 = TorsionPotential {
            atoms: (0, 1, 2, 3),
            v_half: 2.0,
            n: 3,
            cos_n_phi0: -1.0,
            sin_n_phi0: 0.0,
        };

        assert_eq!(t1, t2);
        assert_ne!(t1, t3);
    }

    #[test]
    fn torsion_clone() {
        let torsion = TorsionPotential {
            atoms: (5, 6, 7, 8),
            v_half: 22.5,
            n: 2,
            cos_n_phi0: 1.0,
            sin_n_phi0: 0.0,
        };
        let cloned = torsion.clone();
        assert_eq!(torsion, cloned);
    }

    #[test]
    fn inversion_planar_construction() {
        let inversion = InversionPotential::Planar {
            atoms: (1, 0, 2, 3),
            c_half: 20.0,
        };

        match inversion {
            InversionPotential::Planar { atoms, c_half } => {
                assert_eq!(atoms.0, 1);
                assert_eq!(atoms.1, 0);
                assert_eq!(atoms.2, 2);
                assert_eq!(atoms.3, 3);
                assert_eq!(c_half, 20.0);
            }
            _ => panic!("expected Planar variant"),
        }
    }

    #[test]
    fn inversion_umbrella_construction() {
        let inversion = InversionPotential::Umbrella {
            atoms: (5, 6, 7, 8),
            c_half: 13.33,
            cos_psi0: 0.5774,
        };

        match inversion {
            InversionPotential::Umbrella {
                atoms,
                c_half,
                cos_psi0,
            } => {
                assert_eq!(atoms.0, 5);
                assert_eq!(atoms.1, 6);
                assert_eq!(atoms.2, 7);
                assert_eq!(atoms.3, 8);
                assert_eq!(c_half, 13.33);
                assert_eq!(cos_psi0, 0.5774);
            }
            _ => panic!("expected Umbrella variant"),
        }
    }

    #[test]
    fn inversion_partial_eq() {
        let i1 = InversionPotential::Planar {
            atoms: (1, 0, 2, 3),
            c_half: 20.0,
        };
        let i2 = InversionPotential::Planar {
            atoms: (1, 0, 2, 3),
            c_half: 20.0,
        };
        assert_eq!(i1, i2);
    }

    #[test]
    fn vdw_lennard_jones_construction() {
        let vdw = VdwPairPotential::LennardJones {
            type1_idx: 0,
            type2_idx: 1,
            d0: 0.0951,
            r0_sq: 15.195,
        };

        match vdw {
            VdwPairPotential::LennardJones {
                type1_idx,
                type2_idx,
                d0,
                r0_sq,
            } => {
                assert_eq!(type1_idx, 0);
                assert_eq!(type2_idx, 1);
                assert_eq!(d0, 0.0951);
                assert_eq!(r0_sq, 15.195);
            }
            _ => panic!("expected LennardJones variant"),
        }
    }

    #[test]
    fn vdw_buckingham_construction() {
        let vdw = VdwPairPotential::Buckingham {
            type1_idx: 2,
            type2_idx: 3,
            a: 1000.0,
            b: 3.08,
            c: 500.0,
            r_max_sq: 16.0,
            two_e_max: 10.5,
        };

        match vdw {
            VdwPairPotential::Buckingham {
                type1_idx,
                type2_idx,
                a,
                b,
                c,
                r_max_sq,
                two_e_max,
            } => {
                assert_eq!(type1_idx, 2);
                assert_eq!(type2_idx, 3);
                assert_eq!(a, 1000.0);
                assert_eq!(b, 3.08);
                assert_eq!(c, 500.0);
                assert_eq!(r_max_sq, 16.0);
                assert_eq!(two_e_max, 10.5);
            }
            _ => panic!("expected Buckingham variant"),
        }
    }

    #[test]
    fn vdw_partial_eq() {
        let v1 = VdwPairPotential::LennardJones {
            type1_idx: 0,
            type2_idx: 1,
            d0: 0.0951,
            r0_sq: 15.195,
        };
        let v2 = VdwPairPotential::LennardJones {
            type1_idx: 0,
            type2_idx: 1,
            d0: 0.0951,
            r0_sq: 15.195,
        };
        let v3 = VdwPairPotential::LennardJones {
            type1_idx: 0,
            type2_idx: 1,
            d0: 0.1000,
            r0_sq: 15.195,
        };

        assert_eq!(v1, v2);
        assert_ne!(v1, v3);
    }

    #[test]
    fn hbond_construction() {
        let hbond = HBondPotential {
            donor_type_idx: 0,
            hydrogen_type_idx: 1,
            acceptor_type_idx: 2,
            d_hb: 4.0,
            r_hb_sq: 7.5625,
        };

        assert_eq!(hbond.donor_type_idx, 0);
        assert_eq!(hbond.hydrogen_type_idx, 1);
        assert_eq!(hbond.acceptor_type_idx, 2);
        assert_eq!(hbond.d_hb, 4.0);
        assert_eq!(hbond.r_hb_sq, 7.5625);
    }

    #[test]
    fn hbond_partial_eq() {
        let h1 = HBondPotential {
            donor_type_idx: 0,
            hydrogen_type_idx: 1,
            acceptor_type_idx: 2,
            d_hb: 4.0,
            r_hb_sq: 7.5625,
        };
        let h2 = HBondPotential {
            donor_type_idx: 0,
            hydrogen_type_idx: 1,
            acceptor_type_idx: 2,
            d_hb: 4.0,
            r_hb_sq: 7.5625,
        };
        let h3 = HBondPotential {
            donor_type_idx: 0,
            hydrogen_type_idx: 1,
            acceptor_type_idx: 2,
            d_hb: 5.0,
            r_hb_sq: 7.5625,
        };

        assert_eq!(h1, h2);
        assert_ne!(h1, h3);
    }

    #[test]
    fn hbond_clone() {
        let hbond = HBondPotential {
            donor_type_idx: 3,
            hydrogen_type_idx: 4,
            acceptor_type_idx: 5,
            d_hb: 7.0,
            r_hb_sq: 7.5625,
        };
        let cloned = hbond.clone();
        assert_eq!(hbond, cloned);
    }

    #[test]
    fn atom_param_construction() {
        let param = AtomParam {
            charge: -0.5,
            mass: 12.011,
            type_idx: 2,
        };

        assert_eq!(param.charge, -0.5);
        assert_eq!(param.mass, 12.011);
        assert_eq!(param.type_idx, 2);
    }

    #[test]
    fn atom_param_partial_eq() {
        let p1 = AtomParam {
            charge: 0.0,
            mass: 1.008,
            type_idx: 0,
        };
        let p2 = AtomParam {
            charge: 0.0,
            mass: 1.008,
            type_idx: 0,
        };
        let p3 = AtomParam {
            charge: 0.1,
            mass: 1.008,
            type_idx: 0,
        };

        assert_eq!(p1, p2);
        assert_ne!(p1, p3);
    }

    #[test]
    fn atom_param_clone() {
        let param = AtomParam {
            charge: -0.3,
            mass: 15.999,
            type_idx: 5,
        };
        let cloned = param.clone();
        assert_eq!(param, cloned);
    }

    #[test]
    fn potentials_default() {
        let pots = Potentials::default();

        assert_eq!(pots.bonds.len(), 0);
        assert_eq!(pots.angles.len(), 0);
        assert_eq!(pots.torsions.len(), 0);
        assert_eq!(pots.inversions.len(), 0);
        assert_eq!(pots.vdw_pairs.len(), 0);
        assert_eq!(pots.h_bonds.len(), 0);
    }

    #[test]
    fn potentials_add_elements() {
        let mut pots = Potentials::default();

        pots.bonds.push(BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 350.0,
            r0: 1.54,
        });

        pots.angles.push(AnglePotential::CosineHarmonic {
            atoms: (0, 1, 2),
            k_half: 50.0,
            cos0: -0.333,
        });

        pots.torsions.push(TorsionPotential {
            atoms: (0, 1, 2, 3),
            v_half: 1.0,
            n: 3,
            cos_n_phi0: -1.0,
            sin_n_phi0: 0.0,
        });

        assert_eq!(pots.bonds.len(), 1);
        assert_eq!(pots.angles.len(), 1);
        assert_eq!(pots.torsions.len(), 1);
    }

    #[test]
    fn potentials_clone() {
        let mut pots = Potentials::default();
        pots.bonds.push(BondPotential::Harmonic {
            atoms: (0, 1),
            k_half: 350.0,
            r0: 1.54,
        });

        let cloned = pots.clone();
        assert_eq!(pots.bonds.len(), cloned.bonds.len());
    }
}
