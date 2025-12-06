use super::system::System;

#[derive(Debug, Clone, PartialEq)]
pub struct AtomParam {
    pub charge: f64,
    pub mass: f64,
    pub type_index: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BondPotential {
    Harmonic {
        i: usize,
        j: usize,
        k: f64,
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
    pub donor_idx: usize,
    pub hydrogen_idx: usize,
    pub acceptor_idx: usize,
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
