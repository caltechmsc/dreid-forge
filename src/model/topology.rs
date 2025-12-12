use super::system::System;

#[derive(Debug, Clone, PartialEq)]
pub struct AtomParam {
    pub charge: f64,
    pub mass: f64,
    pub type_idx: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BondPotential {
    Harmonic {
        i: usize,
        j: usize,
        k_force: f64,
        r0: f64,
    },
    Morse {
        i: usize,
        j: usize,
        r0: f64,
        d0: f64,
        alpha: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum AnglePotential {
    CosineHarmonic {
        i: usize,
        j: usize,
        k: usize,
        k_force: f64,
        theta0: f64,
    },
    ThetaHarmonic {
        i: usize,
        j: usize,
        k: usize,
        k_force: f64,
        theta0: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct DihedralPotential {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub v_barrier: f64,
    pub periodicity: i32,
    pub phase_offset: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ImproperPotential {
    Planar {
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        k_force: f64,
        chi0: f64,
    },
    Umbrella {
        center: usize,
        p1: usize,
        p2: usize,
        p3: usize,
        k_force: f64,
        psi0: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum VdwPairPotential {
    LennardJones {
        type1_idx: usize,
        type2_idx: usize,
        sigma: f64,
        epsilon: f64,
    },
    Exponential6 {
        type1_idx: usize,
        type2_idx: usize,
        a: f64,
        b: f64,
        c: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct HBondPotential {
    pub donor_type_idx: usize,
    pub hydrogen_type_idx: usize,
    pub acceptor_type_idx: usize,
    pub d0: f64,
    pub r0: f64,
}

#[derive(Debug, Clone, Default)]
pub struct Potentials {
    pub bonds: Vec<BondPotential>,
    pub angles: Vec<AnglePotential>,
    pub dihedrals: Vec<DihedralPotential>,
    pub impropers: Vec<ImproperPotential>,
    pub vdw_pairs: Vec<VdwPairPotential>,
    pub h_bonds: Vec<HBondPotential>,
}

#[derive(Debug, Clone)]
pub struct ForgedSystem {
    pub system: System,
    pub atom_types: Vec<String>,
    pub atom_properties: Vec<AtomParam>,
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
