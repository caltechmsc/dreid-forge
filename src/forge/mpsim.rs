//! MPSim compatibility adapter.
//!
//! Implements the two structural/naming conventions required by the legacy
//! MPSim EM/MM engine, driven by [`MpsimConfig`](super::config::MpsimConfig):
//!
//! - [`rename_hb_hydrogens`] rewrites emitted `H_HB` force-field type names to
//!   MPSim's `H___A`.
//! - [`normalize_terminal_hydrogens`] rebuilds protein chain termini into their
//!   neutral protonation state — stripping the extra N-terminal hydrogen
//!   (`NH₃⁺ → NH₂`) and building the C-terminal carboxyl proton
//!   (`COO⁻ → COOH`) — so that the neutral terminal charge sets applied during
//!   hybrid charge assignment stay self-consistent.
//!
//! The neutral-terminal *charge* selection itself lives in the hybrid charge
//! module; this module owns the matching topology transform. Both are gated on
//! the same [`MpsimConfig::neutral_termini`] flag.

use crate::model::atom::Atom;
use crate::model::metadata::{AtomResidueInfo, ResidueCategory, ResiduePosition, StandardResidue};
use crate::model::system::{Bond, System};
use crate::model::types::{BondOrder, Element};
use std::collections::{HashMap, HashSet};

/// DREIDING force-field type for polar (hydrogen-bond) hydrogens.
pub const HB_HYDROGEN_TYPE: &str = "H_HB";
/// MPSim force-field type replacing [`HB_HYDROGEN_TYPE`].
pub const MPSIM_HB_HYDROGEN_TYPE: &str = "H___A";

/// PDB atom name of the extra N-terminal ammonium hydrogen (`NH₃⁺`).
const N_TERM_EXTRA_HYDROGEN: &str = "H3";
/// PDB atom name of the C-terminal carboxyl hydrogen (`COOH`).
const C_TERM_HYDROGEN: &str = "HOXT";

/// O–H bond length for the carboxyl proton (Å); mirrors bio-forge.
const COOH_BOND_LENGTH: f64 = 0.97;
/// sp³ tetrahedral bond angle for hydroxyl placement (degrees).
const SP3_ANGLE_DEG: f64 = 109.5;
/// Dihedral offset used when placing the carboxyl hydrogen (degrees).
const HOXT_DIHEDRAL_DEG: f64 = 60.0;

/// Rewrites `H_HB` force-field type names to MPSim's `H___A`.
///
/// Operates only on the emitted type-name table. Because potentials reference
/// atom types by index and the rename is one-to-one, all downstream references
/// (per-atom type index, hydrogen-bond terms, BGF output) remain valid.
pub fn rename_hb_hydrogens(atom_types: &mut [String]) {
    for atom_type in atom_types.iter_mut() {
        if atom_type == HB_HYDROGEN_TYPE {
            *atom_type = MPSIM_HB_HYDROGEN_TYPE.to_string();
        }
    }
}

/// Normalizes protein chain termini to their neutral protonation state.
///
/// For every standard protein residue at a terminus:
///
/// - **N-terminal** — removes the extra `H3` hydrogen, leaving the neutral
///   amine (`–NH₂`, or `–NH` for proline).
/// - **C-terminal** — adds the `HOXT` carboxyl hydrogen (bonded to `OXT`),
///   forming the neutral acid (`–COOH`), unless it is already present.
///
/// Atom indices, bonds, and biological metadata are kept consistent. Systems
/// without biological metadata (e.g. small molecules) are left untouched, as
/// are termini that already carry the neutral hydrogen set.
pub fn normalize_terminal_hydrogens(system: &mut System) {
    if system.bio_metadata.is_none() {
        return;
    }

    let removals = plan_nterm_removals(system);
    let additions = plan_cterm_additions(system);

    if removals.is_empty() && additions.is_empty() {
        return;
    }

    let old_len = system.atoms.len();

    // Rebuild the atom / metadata arrays, dropping removed atoms and recording
    // the old-index -> new-index mapping (None for removed atoms).
    let mut remap = vec![None; old_len];
    {
        let mut new_atoms = Vec::with_capacity(old_len);
        let mut new_info = Vec::with_capacity(old_len);
        let info = &system.bio_metadata.as_ref().unwrap().atom_info;
        for i in 0..old_len {
            if removals.contains(&i) {
                continue;
            }
            remap[i] = Some(new_atoms.len());
            new_atoms.push(system.atoms[i].clone());
            new_info.push(info[i].clone());
        }
        system.atoms = new_atoms;
        system.bio_metadata.as_mut().unwrap().atom_info = new_info;
    }

    // Remap surviving bonds, dropping any that touched a removed atom.
    let mut new_bonds = Vec::with_capacity(system.bonds.len());
    for bond in &system.bonds {
        if let (Some(i), Some(j)) = (remap[bond.i], remap[bond.j]) {
            new_bonds.push(Bond::new(i, j, bond.order));
        }
    }
    system.bonds = new_bonds;

    // Append the newly built C-terminal carboxyl hydrogens.
    for add in additions {
        let Some(partner) = remap[add.partner] else {
            continue;
        };
        let new_idx = system.atoms.len();
        system.atoms.push(Atom::new(Element::H, add.position));
        system
            .bio_metadata
            .as_mut()
            .unwrap()
            .atom_info
            .push(add.info);
        system
            .bonds
            .push(Bond::new(partner, new_idx, BondOrder::Single));
    }
}

/// Collects atom indices of extra N-terminal hydrogens to remove.
fn plan_nterm_removals(system: &System) -> HashSet<usize> {
    let meta = system.bio_metadata.as_ref().unwrap();
    let mut removals = HashSet::new();

    for (idx, info) in meta.atom_info.iter().enumerate() {
        if info.position == ResiduePosition::NTerminal
            && info.category == ResidueCategory::Standard
            && is_protein_residue(info.standard_name)
            && info.atom_name == N_TERM_EXTRA_HYDROGEN
        {
            removals.insert(idx);
        }
    }

    removals
}

/// A carboxyl hydrogen to be appended, referencing its `OXT` partner by the
/// atom index in the *original* (pre-removal) system.
struct PendingHydrogen {
    position: [f64; 3],
    info: AtomResidueInfo,
    partner: usize,
}

/// Plans C-terminal `HOXT` additions for all neutral-capped protein termini.
fn plan_cterm_additions(system: &System) -> Vec<PendingHydrogen> {
    let meta = system.bio_metadata.as_ref().unwrap();

    // Group C-terminal protein residue atom indices by residue identity.
    let mut groups: HashMap<(&str, i32, Option<char>), Vec<usize>> = HashMap::new();
    for (idx, info) in meta.atom_info.iter().enumerate() {
        if info.position == ResiduePosition::CTerminal
            && info.category == ResidueCategory::Standard
            && is_protein_residue(info.standard_name)
        {
            let key = (info.chain_id.as_str(), info.residue_id, info.insertion_code);
            groups.entry(key).or_default().push(idx);
        }
    }

    let mut additions = Vec::new();
    for atoms in groups.values() {
        let find = |name: &str| {
            atoms
                .iter()
                .copied()
                .find(|&i| meta.atom_info[i].atom_name == name)
        };

        // Skip if already protonated.
        if find(C_TERM_HYDROGEN).is_some() {
            continue;
        }

        let (Some(c_idx), Some(oxt_idx)) = (find("C"), find("OXT")) else {
            continue;
        };
        let Some(ref_idx) = find("CA").or_else(|| find("O")) else {
            continue;
        };

        let position = place_hydroxyl_hydrogen(
            system.atoms[oxt_idx].position,
            system.atoms[c_idx].position,
            system.atoms[ref_idx].position,
        );

        let template = &meta.atom_info[oxt_idx];
        let info = AtomResidueInfo::builder(
            C_TERM_HYDROGEN,
            template.residue_name.clone(),
            template.residue_id,
            template.chain_id.clone(),
        )
        .insertion_code_opt(template.insertion_code)
        .standard_name(template.standard_name)
        .category(template.category)
        .position(template.position)
        .build();

        additions.push(PendingHydrogen {
            position,
            info,
            partner: oxt_idx,
        });
    }

    // Deterministic append order regardless of HashMap iteration order.
    additions.sort_by_key(|a| a.partner);
    additions
}

/// Returns `true` for the 20 standard amino-acid residues.
fn is_protein_residue(name: Option<StandardResidue>) -> bool {
    use StandardResidue::*;
    matches!(
        name,
        Some(
            ALA | ARG
                | ASN
                | ASP
                | CYS
                | GLN
                | GLU
                | GLY
                | HIS
                | ILE
                | LEU
                | LYS
                | MET
                | PHE
                | PRO
                | SER
                | THR
                | TRP
                | TYR
                | VAL
        )
    )
}

/// Places a hydroxyl hydrogen on `oxt` using sp³ tetrahedral geometry, mirroring
/// bio-forge's carboxyl-proton construction (bond length, angle, dihedral).
fn place_hydroxyl_hydrogen(oxt: [f64; 3], attached: [f64; 3], reference: [f64; 3]) -> [f64; 3] {
    let (x, y, z) = build_sp3_frame(oxt, attached, reference);

    let theta = SP3_ANGLE_DEG.to_radians();
    let phi = HOXT_DIHEDRAL_DEG.to_radians();
    let (sin_theta, cos_theta) = (theta.sin(), theta.cos());

    let local = [sin_theta * phi.cos(), sin_theta * phi.sin(), -cos_theta];
    let dir = add(
        add(scale(x, local[0]), scale(y, local[1])),
        scale(z, local[2]),
    );

    add(oxt, scale(dir, COOH_BOND_LENGTH))
}

/// Builds an orthonormal sp³ frame at `center` relative to `attached`, oriented
/// by `reference`.
fn build_sp3_frame(
    center: [f64; 3],
    attached: [f64; 3],
    reference: [f64; 3],
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let z = normalize(sub(center, attached));
    let ref_vec = normalize(sub(reference, attached));
    let x = normalize(sub(ref_vec, scale(z, dot(z, ref_vec))));
    let y = cross(z, x);
    (x, y, z)
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = dot(a, a).sqrt();
    if n < 1e-12 { a } else { scale(a, 1.0 / n) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::metadata::BioMetadata;

    #[test]
    fn renames_only_hb_hydrogens() {
        let mut types = vec![
            "C_3".to_string(),
            "H_HB".to_string(),
            "H_".to_string(),
            "O_3".to_string(),
        ];
        rename_hb_hydrogens(&mut types);
        assert_eq!(types, vec!["C_3", "H___A", "H_", "O_3"]);
    }

    fn atom(
        name: &str,
        res: StandardResidue,
        pos: ResiduePosition,
        xyz: [f64; 3],
    ) -> (Atom, AtomResidueInfo) {
        let element = if name.starts_with('H') {
            Element::H
        } else if name.starts_with('O') {
            Element::O
        } else if name.starts_with('N') {
            Element::N
        } else {
            Element::C
        };
        let info = AtomResidueInfo::builder(name, format!("{:?}", res), 1, "A")
            .standard_name(Some(res))
            .category(ResidueCategory::Standard)
            .position(pos)
            .build();
        (Atom::new(element, xyz), info)
    }

    fn build_system(entries: Vec<(Atom, AtomResidueInfo)>, bonds: Vec<Bond>) -> System {
        let mut atoms = Vec::new();
        let mut atom_info = Vec::new();
        for (a, i) in entries {
            atoms.push(a);
            atom_info.push(i);
        }
        System {
            atoms,
            bonds,
            box_vectors: None,
            bio_metadata: Some(BioMetadata {
                atom_info,
                target_ph: None,
            }),
        }
    }

    #[test]
    fn strips_nterminal_h3() {
        use ResiduePosition::NTerminal;
        use StandardResidue::ALA;
        let entries = vec![
            atom("N", ALA, NTerminal, [0.0, 0.0, 0.0]),
            atom("H1", ALA, NTerminal, [-0.5, 0.8, 0.0]),
            atom("H2", ALA, NTerminal, [-0.5, -0.8, 0.0]),
            atom("H3", ALA, NTerminal, [-0.5, 0.0, 0.8]),
            atom("CA", ALA, NTerminal, [1.4, 0.0, 0.0]),
        ];
        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(0, 2, BondOrder::Single),
            Bond::new(0, 3, BondOrder::Single),
            Bond::new(0, 4, BondOrder::Single),
        ];
        let mut system = build_system(entries, bonds);
        normalize_terminal_hydrogens(&mut system);

        let names: Vec<&str> = system
            .bio_metadata
            .as_ref()
            .unwrap()
            .atom_info
            .iter()
            .map(|i| i.atom_name.as_str())
            .collect();
        assert_eq!(names, vec!["N", "H1", "H2", "CA"]);
        assert_eq!(system.atoms.len(), 4);
        // The N–H3 bond is gone; three bonds remain, all valid indices.
        assert_eq!(system.bonds.len(), 3);
        for b in &system.bonds {
            assert!(b.i < 4 && b.j < 4);
        }
    }

    #[test]
    fn builds_cterminal_hoxt() {
        use ResiduePosition::CTerminal;
        use StandardResidue::ALA;
        let entries = vec![
            atom("CA", ALA, CTerminal, [0.0, 0.0, 0.0]),
            atom("C", ALA, CTerminal, [1.5, 0.0, 0.0]),
            atom("O", ALA, CTerminal, [2.1, 1.0, 0.0]),
            atom("OXT", ALA, CTerminal, [2.1, -1.0, 0.0]),
        ];
        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
            Bond::new(1, 3, BondOrder::Single),
        ];
        let mut system = build_system(entries, bonds);
        normalize_terminal_hydrogens(&mut system);

        let meta = system.bio_metadata.as_ref().unwrap();
        assert_eq!(system.atoms.len(), 5);
        let hoxt_idx = meta
            .atom_info
            .iter()
            .position(|i| i.atom_name == "HOXT")
            .expect("HOXT added");
        assert_eq!(system.atoms[hoxt_idx].element, Element::H);

        // Bond OXT–HOXT exists.
        let oxt_idx = meta
            .atom_info
            .iter()
            .position(|i| i.atom_name == "OXT")
            .unwrap();
        assert!(
            system
                .bonds
                .iter()
                .any(|b| (b.i == oxt_idx && b.j == hoxt_idx) || (b.j == oxt_idx && b.i == hoxt_idx))
        );

        // Placed near OXT at ~0.97 Å.
        let d = dot(
            sub(
                system.atoms[hoxt_idx].position,
                system.atoms[oxt_idx].position,
            ),
            sub(
                system.atoms[hoxt_idx].position,
                system.atoms[oxt_idx].position,
            ),
        )
        .sqrt();
        assert!((d - COOH_BOND_LENGTH).abs() < 1e-6, "HOXT bond length {d}");
    }

    #[test]
    fn idempotent_when_already_neutral() {
        use ResiduePosition::CTerminal;
        use StandardResidue::ALA;
        let entries = vec![
            atom("CA", ALA, CTerminal, [0.0, 0.0, 0.0]),
            atom("C", ALA, CTerminal, [1.5, 0.0, 0.0]),
            atom("O", ALA, CTerminal, [2.1, 1.0, 0.0]),
            atom("OXT", ALA, CTerminal, [2.1, -1.0, 0.0]),
            atom("HOXT", ALA, CTerminal, [3.0, -1.0, 0.0]),
        ];
        let bonds = vec![
            Bond::new(1, 3, BondOrder::Single),
            Bond::new(3, 4, BondOrder::Single),
        ];
        let mut system = build_system(entries, bonds);
        let before = system.atoms.len();
        normalize_terminal_hydrogens(&mut system);
        assert_eq!(system.atoms.len(), before, "no duplicate HOXT");
    }

    #[test]
    fn no_metadata_is_noop() {
        let mut system = System {
            atoms: vec![Atom::new(Element::O, [0.0, 0.0, 0.0])],
            bonds: vec![],
            box_vectors: None,
            bio_metadata: None,
        };
        normalize_terminal_hydrogens(&mut system);
        assert_eq!(system.atoms.len(), 1);
    }
}
