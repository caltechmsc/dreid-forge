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
//! - **Angles**: [`AnglePotential`] — Cosine-harmonic or theta-harmonic bending
//! - **Dihedrals**: [`DihedralPotential`] — Periodic torsion potentials
//! - **Impropers**: [`ImproperPotential`] — Out-of-plane (planar/umbrella)
//! - **Van der Waals**: [`VdwPairPotential`] — LJ 12-6 or Exp-6
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
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
        /// Force constant (kcal/mol/Å²).
        k_force: f64,
        /// Equilibrium bond length (Å).
        r0: f64,
    },
    /// Morse anharmonic bond potential.
    Morse {
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
        /// Equilibrium bond length (Å).
        r0: f64,
        /// Bond dissociation energy (kcal/mol).
        d0: f64,
        /// Morse alpha parameter (Å⁻¹).
        alpha: f64,
    },
}

/// Angle bending potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum AnglePotential {
    /// Cosine-harmonic angle potential (DREIDING default).
    CosineHarmonic {
        /// First atom index.
        i: usize,
        /// Central atom index.
        j: usize,
        /// Third atom index.
        k: usize,
        /// Force constant (kcal/mol/rad²).
        k_force: f64,
        /// Equilibrium angle (degrees).
        theta0: f64,
    },
    /// Simple theta-harmonic angle potential.
    ThetaHarmonic {
        /// First atom index.
        i: usize,
        /// Central atom index.
        j: usize,
        /// Third atom index.
        k: usize,
        /// Force constant (kcal/mol/rad²).
        k_force: f64,
        /// Equilibrium angle (degrees).
        theta0: f64,
    },
}

/// Proper dihedral (torsion) potential.
#[derive(Debug, Clone, PartialEq)]
pub struct DihedralPotential {
    /// First atom index.
    pub i: usize,
    /// Second atom index (bond axis).
    pub j: usize,
    /// Third atom index (bond axis).
    pub k: usize,
    /// Fourth atom index.
    pub l: usize,
    /// Barrier height (kcal/mol).
    pub v_barrier: f64,
    /// Periodicity (number of minima in 360°).
    pub periodicity: i32,
    /// Phase offset (degrees).
    pub phase_offset: f64,
}

/// Improper dihedral (out-of-plane) potential functions.
#[derive(Debug, Clone, PartialEq)]
pub enum ImproperPotential {
    /// Planar improper for sp² centers.
    Planar {
        /// First peripheral atom.
        i: usize,
        /// Central atom (sp² center).
        j: usize,
        /// Second peripheral atom.
        k: usize,
        /// Third peripheral atom.
        l: usize,
        /// Force constant (kcal/mol/rad²).
        k_force: f64,
        /// Equilibrium out-of-plane angle (degrees).
        chi0: f64,
    },
    /// Umbrella improper for pyramidal centers.
    Umbrella {
        /// Central atom (pyramidal center).
        center: usize,
        /// First peripheral atom.
        p1: usize,
        /// Second peripheral atom.
        p2: usize,
        /// Third peripheral atom.
        p3: usize,
        /// Force constant (kcal/mol).
        k_force: f64,
        /// Equilibrium umbrella angle (degrees).
        psi0: f64,
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
        /// LJ sigma parameter (Å).
        sigma: f64,
        /// LJ epsilon parameter (kcal/mol).
        epsilon: f64,
    },
    /// Exponential-6 potential (Buckingham-like).
    Exponential6 {
        /// First atom type index.
        type1_idx: usize,
        /// Second atom type index.
        type2_idx: usize,
        /// Repulsive prefactor A.
        a: f64,
        /// Exponential decay parameter B (Å⁻¹).
        b: f64,
        /// Attractive coefficient C (kcal·Å⁶/mol).
        c: f64,
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
    /// H-bond equilibrium energy (kcal/mol).
    pub d0: f64,
    /// Equilibrium H···A distance (Å).
    pub r0: f64,
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
    /// Proper dihedral (torsion) potentials.
    pub dihedrals: Vec<DihedralPotential>,
    /// Improper dihedral (out-of-plane) potentials.
    pub impropers: Vec<ImproperPotential>,
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
    use crate::model::system::System;

    #[test]
    fn atom_param_fields_and_clone() {
        let p = AtomParam {
            charge: -0.34,
            mass: 12.011,
            type_idx: 2,
        };
        assert!(p.charge < 0.0);
        assert_eq!(p.mass, 12.011);
        assert_eq!(p.type_idx, 2);
        let q = p.clone();
        assert_eq!(p, q);
    }

    #[test]
    fn bond_potential_variants_and_debug() {
        let h = BondPotential::Harmonic {
            i: 0,
            j: 1,
            k_force: 300.0,
            r0: 1.23,
        };
        let m = BondPotential::Morse {
            i: 1,
            j: 2,
            r0: 1.5,
            d0: 4.0,
            alpha: 2.0,
        };
        match h {
            BondPotential::Harmonic { i, j, k_force, r0 } => {
                assert_eq!(i, 0);
                assert_eq!(j, 1);
                assert!(k_force > 0.0);
                assert!(r0 > 0.0);
            }
            _ => panic!("expected Harmonic variant"),
        }
        match m {
            BondPotential::Morse {
                i,
                j,
                r0,
                d0,
                alpha,
            } => {
                assert_eq!(i, 1);
                assert_eq!(j, 2);
                assert!(d0 > 0.0);
                assert!(alpha > 0.0);
                assert!(r0 > 0.0);
            }
            _ => panic!("expected Morse variant"),
        }
        let s = format!("{:?} {:?}", h, m);
        assert!(s.contains("Harmonic"));
        assert!(s.contains("Morse"));
    }

    #[test]
    fn angle_potential_variants() {
        let a1 = AnglePotential::CosineHarmonic {
            i: 0,
            j: 1,
            k: 2,
            k_force: 50.0,
            theta0: 109.5,
        };
        let a2 = AnglePotential::ThetaHarmonic {
            i: 2,
            j: 1,
            k: 0,
            k_force: 40.0,
            theta0: 120.0,
        };
        match a1 {
            AnglePotential::CosineHarmonic {
                k_force, theta0, ..
            } => {
                assert_eq!(k_force, 50.0);
                assert_eq!(theta0, 109.5);
            }
            _ => panic!("expected CosineHarmonic"),
        }
        match a2 {
            AnglePotential::ThetaHarmonic {
                k_force, theta0, ..
            } => {
                assert_eq!(k_force, 40.0);
                assert_eq!(theta0, 120.0);
            }
            _ => panic!("expected ThetaHarmonic"),
        }
    }

    #[test]
    fn dihedral_and_improper_variants() {
        let d = DihedralPotential {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            v_barrier: 2.5,
            periodicity: 3,
            phase_offset: 180.0,
        };
        assert_eq!(d.periodicity, 3);
        assert!(d.v_barrier > 0.0);

        let imp1 = ImproperPotential::Planar {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            k_force: 10.0,
            chi0: 0.0,
        };
        let imp2 = ImproperPotential::Umbrella {
            center: 1,
            p1: 2,
            p2: 3,
            p3: 4,
            k_force: 5.0,
            psi0: 180.0,
        };
        match imp1 {
            ImproperPotential::Planar { k_force, chi0, .. } => {
                assert_eq!(k_force, 10.0);
                assert_eq!(chi0, 0.0);
            }
            _ => panic!("expected Planar"),
        }
        match imp2 {
            ImproperPotential::Umbrella {
                center, p1, p2, p3, ..
            } => {
                assert_eq!(center, 1);
                assert_eq!(p1, 2);
                assert_eq!(p2, 3);
                assert_eq!(p3, 4);
            }
            _ => panic!("expected Umbrella"),
        }
    }

    #[test]
    fn vdw_and_hbond_variants_and_potentials_container() {
        let lj = VdwPairPotential::LennardJones {
            type1_idx: 0,
            type2_idx: 1,
            sigma: 3.5,
            epsilon: 0.2,
        };
        let ex6 = VdwPairPotential::Exponential6 {
            type1_idx: 1,
            type2_idx: 2,
            a: 1000.0,
            b: 50.0,
            c: 2.0,
        };
        match lj {
            VdwPairPotential::LennardJones { sigma, epsilon, .. } => {
                assert_eq!(sigma, 3.5);
                assert_eq!(epsilon, 0.2);
            }
            _ => panic!("expected LennardJones"),
        }
        match ex6 {
            VdwPairPotential::Exponential6 { a, b, c, .. } => {
                assert!(a > 0.0);
                assert!(b > 0.0);
                assert!(c > 0.0);
            }
            _ => panic!("expected Exponential6"),
        }

        let hb = HBondPotential {
            donor_type_idx: 0,
            hydrogen_type_idx: 1,
            acceptor_type_idx: 2,
            d0: 9.5,
            r0: 2.75,
        };
        assert_eq!(hb.donor_type_idx, 0);
        assert_eq!(hb.hydrogen_type_idx, 1);
        assert_eq!(hb.acceptor_type_idx, 2);
        assert_eq!(hb.d0, 9.5);
        assert_eq!(hb.r0, 2.75);

        let mut pots = Potentials::default();
        pots.vdw_pairs.push(lj.clone());
        pots.vdw_pairs.push(ex6.clone());
        pots.h_bonds.push(hb.clone());
        assert_eq!(pots.vdw_pairs.len(), 2);
        assert_eq!(pots.h_bonds.len(), 1);
    }

    #[test]
    fn forged_system_basic_fields_and_debug() {
        let sys = System::new();
        let atom_types = vec!["C_3".to_string(), "O_2".to_string()];
        let atom_props = vec![
            AtomParam {
                charge: -0.1,
                mass: 12.0,
                type_idx: 0,
            },
            AtomParam {
                charge: -0.2,
                mass: 16.0,
                type_idx: 1,
            },
        ];
        let pots = Potentials::default();
        let fs = ForgedSystem {
            system: sys.clone(),
            atom_types: atom_types.clone(),
            atom_properties: atom_props.clone(),
            potentials: pots,
        };
        assert_eq!(fs.atom_types.len(), 2);
        assert_eq!(fs.atom_properties.len(), 2);
        let s = format!("{:?}", fs);
        assert!(s.contains("ForgedSystem"));
        assert!(s.contains("atom_types"));
        assert!(s.contains("atom_properties"));
    }
}
