//! MPSim compatibility adapter.
//!
//! Implements the conventions of the legacy MPSim BGF conversion, driven by
//! [`MpsimConfig`](super::config::MpsimConfig):
//!
//! - [`rename_hb_hydrogens`] rewrites emitted `H_HB` force-field type names to
//!   MPSim's `H___A`.
//! - [`normalize_open_ends`] caps every open chain end — a real terminus and a
//!   break in the model alike — in its neutral state: `–NH₂` on the N side and
//!   an aldehyde `–CHO` on the C side.
//! - [`open_end_charges`] supplies the charges those caps carry. They are flat
//!   conventions of the converter rather than a residue library's values, and
//!   the aldehyde cap has no library entry at all because no residue ends in an
//!   aldehyde.
//! - [`relabel_ligand_residues`] forces the single residue label `RES 999` onto
//!   every ligand.
//!
//! The structural transforms run on the input system before typing, so the typer
//! reads the final connectivity. The charges are handed to charge assignment as
//! preassigned values, since the lookups they replace either do not exist or
//! would return the modern parameterization instead of the convention.

use crate::model::atom::Atom;
use crate::model::metadata::{AtomResidueInfo, ResidueCategory, StandardResidue};
use crate::model::system::{Bond, System};
use crate::model::types::{BondOrder, Element};
use std::collections::{HashMap, HashSet};

use super::config::{MPSIM_LIGAND_RESIDUE_ID, MPSIM_LIGAND_RESIDUE_NAME};
use super::intermediate::IntermediateSystem;

/// DREIDING force-field type for polar (hydrogen-bond) hydrogens.
pub const HB_HYDROGEN_TYPE: &str = "H_HB";
/// MPSim force-field type replacing [`HB_HYDROGEN_TYPE`].
pub const MPSIM_HB_HYDROGEN_TYPE: &str = "H___A";

/// PDB atom name of the carboxyl oxygen the aldehyde cap replaces.
const CARBOXYL_OXYGEN: &str = "OXT";
/// PDB atom name of the carboxyl proton removed alongside it.
const CARBOXYL_HYDROGEN: &str = "HOXT";
/// PDB atom name of the aldehyde cap hydrogen bonded to the backbone carbon.
const ALDEHYDE_HYDROGEN: &str = "HC";

/// Sigma bonds on a neutral amine nitrogen, which fixes how many hydrogens an
/// `–NH₂` cap carries once its heavy neighbours are counted — two for most
/// residues, one for proline, whose ring supplies a second heavy neighbour.
const AMINE_SIGMA_BONDS: usize = 3;

/// N–H bond length of the amine cap (Å), as written by the converter.
const AMINE_BOND_LENGTH: f64 = 1.000;
/// C–H bond length of the aldehyde cap (Å), as written by the converter.
const ALDEHYDE_BOND_LENGTH: f64 = 1.020;
/// H–C=O angle of the aldehyde cap (degrees), as written by the converter.
const ALDEHYDE_ANGLE_DEG: f64 = 120.0;
/// Tetrahedral angle separating two hydrogens built on one nitrogen (degrees).
const SP3_ANGLE_DEG: f64 = 109.5;

/// Charge the converter writes on a capped amine nitrogen.
///
/// Flat values, identical at a real N-terminus and at a break and identical
/// across every occurrence in the reference set. They replace the residue
/// library's neutral N-terminal set, which distributes the same net charge
/// differently.
const AMINE_NITROGEN_CHARGE: f64 = -0.47;
/// Charge on each hydrogen of a capped amine nitrogen.
const AMINE_HYDROGEN_CHARGE: f64 = 0.155;

/// Charge on the backbone carbon of an aldehyde cap.
///
/// The three aldehyde charges sum to zero. No residue library defines an
/// aldehyde C-terminus, so they come from the reference files.
const ALDEHYDE_CARBON_CHARGE: f64 = 0.51;
/// Charge on the carbonyl oxygen of an aldehyde cap.
const ALDEHYDE_OXYGEN_CHARGE: f64 = -0.51;
/// Charge on the cap hydrogen of an aldehyde cap.
const ALDEHYDE_HYDROGEN_CHARGE: f64 = 0.0;

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

/// Caps every open chain end in its neutral state.
///
/// An end is *open* when the backbone bond that would continue the chain is
/// absent — at a real terminus, and equally at a break where the model simply
/// stops. Both are treated identically, which is the convention this adapter
/// reproduces:
///
/// - **N side** — the nitrogen is left with exactly the hydrogens a neutral
///   amine carries. An ammonium's surplus proton is removed; a break's missing
///   proton is built.
/// - **C side** — `OXT` and its proton are removed and a cap hydrogen `HC` is
///   bonded to the backbone carbon, in the `CA–C=O` plane.
///
/// Atom indices, bonds and biological metadata are kept consistent. Systems
/// without biological metadata are left untouched, as are ends already carrying
/// the neutral hydrogen set.
pub fn normalize_open_ends(system: &mut System) {
    if system.bio_metadata.is_none() {
        return;
    }

    let (removals, additions) = {
        let structure = StructureIndex::new(system);
        let open = structure.open_ends();
        let removals = plan_removals(&structure, &open);
        let additions = plan_additions(&structure, &open, &removals);
        (removals, additions)
    };
    if removals.is_empty() && additions.is_empty() {
        return;
    }

    let remap = drop_atoms(system, &removals);
    append_hydrogens(system, additions, &remap);
}

/// Rebuilds the atom, metadata and bond arrays without the removed atoms,
/// returning the old-index to new-index mapping.
fn drop_atoms(system: &mut System, removals: &HashSet<usize>) -> Vec<Option<usize>> {
    let old_len = system.atoms.len();
    let mut remap = vec![None; old_len];

    let mut atoms = Vec::with_capacity(old_len);
    let mut info = Vec::with_capacity(old_len);
    let old_info = &system.bio_metadata.as_ref().unwrap().atom_info;
    for i in 0..old_len {
        if removals.contains(&i) {
            continue;
        }
        remap[i] = Some(atoms.len());
        atoms.push(system.atoms[i].clone());
        info.push(old_info[i].clone());
    }
    system.atoms = atoms;
    system.bio_metadata.as_mut().unwrap().atom_info = info;

    system.bonds = system
        .bonds
        .iter()
        .filter_map(|b| match (remap[b.i], remap[b.j]) {
            (Some(i), Some(j)) => Some(Bond::new(i, j, b.order)),
            _ => None,
        })
        .collect();

    remap
}

/// Appends the planned hydrogens, translating each partner index through
/// `remap`.
fn append_hydrogens(system: &mut System, additions: Vec<PendingHydrogen>, remap: &[Option<usize>]) {
    for add in additions {
        let Some(partner) = remap[add.partner] else {
            continue;
        };
        let idx = system.atoms.len();
        system.atoms.push(Atom::new(Element::H, add.position));
        system
            .bio_metadata
            .as_mut()
            .unwrap()
            .atom_info
            .push(add.info);
        system
            .bonds
            .push(Bond::new(partner, idx, BondOrder::Single));
    }
}

/// Charges carried by every capped chain end, keyed by atom index.
///
/// Covers the amine nitrogen and its hydrogens on the N side, and the carbonyl
/// carbon, its oxygen and the cap hydrogen on the C side. These are the
/// converter's flat conventions: they replace the library's neutral N-terminal
/// set, and on the C side they stand in for a lookup that cannot succeed at all.
///
/// Only an end [`normalize_open_ends`] actually reached is listed; one left
/// alone for want of a backbone atom keeps whatever the library assigns it.
pub fn open_end_charges(system: &IntermediateSystem) -> HashMap<usize, f64> {
    let Some(meta) = system.bio_metadata.as_ref() else {
        return HashMap::new();
    };

    let mut residues: HashMap<ResidueKey, Vec<usize>> = HashMap::new();
    for (idx, info) in meta.atom_info.iter().enumerate() {
        if is_protein(info) {
            residues.entry(residue_key(info)).or_default().push(idx);
        }
    }

    let name_of = |idx: usize| meta.atom_info[idx].atom_name.as_str();
    let continues = |idx: usize| {
        let own = residue_key(&meta.atom_info[idx]);
        system.atoms[idx]
            .neighbors
            .iter()
            .any(|&n| name_of(n) == "C" && residue_key(&meta.atom_info[n]) != own)
    };

    let mut out = HashMap::new();
    for atoms in residues.values() {
        let find = |wanted: &str| atoms.iter().copied().find(|&i| name_of(i) == wanted);

        // N side: capped exactly when the chain does not continue past the
        // nitrogen, whatever hydrogens it ended up carrying.
        if let Some(n) = find("N")
            && !continues(n)
        {
            out.insert(n, AMINE_NITROGEN_CHARGE);
            for &h in &system.atoms[n].neighbors {
                if name_of(h).starts_with('H') {
                    out.insert(h, AMINE_HYDROGEN_CHARGE);
                }
            }
        }

        // C side: capped exactly when the cap hydrogen is present.
        if let Some(hc) = find(ALDEHYDE_HYDROGEN) {
            out.insert(hc, ALDEHYDE_HYDROGEN_CHARGE);
            if let Some(c) = find("C") {
                out.insert(c, ALDEHYDE_CARBON_CHARGE);
            }
            if let Some(o) = find("O") {
                out.insert(o, ALDEHYDE_OXYGEN_CHARGE);
            }
        }
    }

    out
}

/// Forces the MPSim ligand label onto every ligand in the system.
///
/// Ligands are the hetero-compound residues; ions, water and the biopolymers
/// keep their own labels. Chain assignment is untouched.
pub fn relabel_ligand_residues(system: &mut System) {
    let Some(meta) = system.bio_metadata.as_mut() else {
        return;
    };

    for info in meta.atom_info.iter_mut() {
        if info.category == ResidueCategory::Hetero {
            info.residue_name = MPSIM_LIGAND_RESIDUE_NAME.to_string();
            info.residue_id = MPSIM_LIGAND_RESIDUE_ID;
            info.insertion_code = None;
        }
    }
}

/// Residue identity: chain, sequence number, insertion code.
type ResidueKey = (String, i32, Option<char>);

fn residue_key(info: &AtomResidueInfo) -> ResidueKey {
    (info.chain_id.clone(), info.residue_id, info.insertion_code)
}

fn is_protein(info: &AtomResidueInfo) -> bool {
    info.category == ResidueCategory::Standard && is_protein_residue(info.standard_name)
}

/// The residues whose chain stops, grouped by the side it stops on.
#[derive(Default)]
struct OpenEnds {
    n: HashSet<ResidueKey>,
    c: HashSet<ResidueKey>,
}

/// An index over the structure — adjacency plus each residue's atoms — built
/// once so locating the open ends, planning removals and planning additions all
/// work from the same picture.
struct StructureIndex<'a> {
    atoms: &'a [Atom],
    info: &'a [AtomResidueInfo],
    neighbors: Vec<Vec<usize>>,
    residues: HashMap<ResidueKey, Vec<usize>>,
}

impl<'a> StructureIndex<'a> {
    fn new(system: &'a System) -> Self {
        let info = &system.bio_metadata.as_ref().unwrap().atom_info;

        let mut neighbors = vec![Vec::new(); system.atoms.len()];
        for bond in &system.bonds {
            neighbors[bond.i].push(bond.j);
            neighbors[bond.j].push(bond.i);
        }

        let mut residues: HashMap<ResidueKey, Vec<usize>> = HashMap::new();
        for (idx, atom_info) in info.iter().enumerate() {
            if is_protein(atom_info) {
                residues
                    .entry(residue_key(atom_info))
                    .or_default()
                    .push(idx);
            }
        }

        Self {
            atoms: &system.atoms,
            info,
            neighbors,
            residues,
        }
    }

    fn name_of(&self, idx: usize) -> &str {
        &self.info[idx].atom_name
    }

    /// The atom called `name` in `residue`, skipping any already withdrawn.
    fn find(&self, residue: &ResidueKey, name: &str, withdrawn: &HashSet<usize>) -> Option<usize> {
        self.residues
            .get(residue)?
            .iter()
            .copied()
            .find(|&i| self.name_of(i) == name && !withdrawn.contains(&i))
    }

    /// The neighbours of `idx` split into heavy atoms and hydrogens, skipping
    /// any already withdrawn.
    fn split_neighbors(&self, idx: usize, withdrawn: &HashSet<usize>) -> (Vec<usize>, Vec<usize>) {
        self.neighbors[idx]
            .iter()
            .copied()
            .filter(|n| !withdrawn.contains(n))
            .partition(|&n| !self.name_of(n).starts_with('H'))
    }

    /// Whether the chain continues past `idx` — that is, whether `idx` is bonded
    /// to a backbone atom called `partner` in a *different* residue.
    fn chain_continues(&self, idx: usize, partner: &str) -> bool {
        let own = residue_key(&self.info[idx]);
        self.neighbors[idx]
            .iter()
            .any(|&n| self.name_of(n) == partner && residue_key(&self.info[n]) != own)
    }

    /// Locates every residue whose chain stops, on either side.
    fn open_ends(&self) -> OpenEnds {
        let mut open = OpenEnds::default();
        for (idx, info) in self.info.iter().enumerate() {
            if !is_protein(info) {
                continue;
            }
            match info.atom_name.as_str() {
                "N" if !self.chain_continues(idx, "C") => {
                    open.n.insert(residue_key(info));
                }
                "C" if !self.chain_continues(idx, "N") => {
                    open.c.insert(residue_key(info));
                }
                _ => {}
            }
        }
        open
    }

    /// An `AtomResidueInfo` for a new atom called `name`, taking its residue
    /// identity from an existing atom.
    fn derive_info(&self, from: usize, name: &str) -> AtomResidueInfo {
        let t = &self.info[from];
        AtomResidueInfo::builder(
            name,
            t.residue_name.clone(),
            t.residue_id,
            t.chain_id.clone(),
        )
        .insertion_code_opt(t.insertion_code)
        .standard_name(t.standard_name)
        .category(t.category)
        .position(t.position)
        .build()
    }
}

/// Collects the atoms the caps replace: an ammonium's surplus hydrogens on the
/// N side, and the whole carboxyl group on the C side.
fn plan_removals(structure: &StructureIndex, open: &OpenEnds) -> HashSet<usize> {
    let none = HashSet::new();
    let mut removals = HashSet::new();

    for residue in &open.c {
        for name in [CARBOXYL_OXYGEN, CARBOXYL_HYDROGEN] {
            if let Some(idx) = structure.find(residue, name, &none) {
                removals.insert(idx);
            }
        }
    }

    for residue in &open.n {
        let Some(n) = structure.find(residue, "N", &none) else {
            continue;
        };
        let (heavy, mut hydrogens) = structure.split_neighbors(n, &none);
        let wanted = AMINE_SIGMA_BONDS.saturating_sub(heavy.len());
        if hydrogens.len() <= wanted {
            continue;
        }
        // Drop the surplus highest-numbered hydrogens, so which ones go does not
        // depend on the order the input happened to list them in.
        hydrogens.sort_by(|&a, &b| structure.name_of(a).cmp(structure.name_of(b)));
        removals.extend(hydrogens.into_iter().skip(wanted));
    }

    removals
}

/// A hydrogen to be appended, referencing its partner by the atom index in the
/// *original* (pre-removal) system.
struct PendingHydrogen {
    position: [f64; 3],
    info: AtomResidueInfo,
    partner: usize,
}

/// Plans the hydrogens the caps are missing, on both sides.
fn plan_additions(
    structure: &StructureIndex,
    open: &OpenEnds,
    removals: &HashSet<usize>,
) -> Vec<PendingHydrogen> {
    let mut additions: Vec<PendingHydrogen> = open
        .c
        .iter()
        .filter_map(|r| plan_aldehyde_cap(structure, r, removals))
        .collect();
    for residue in &open.n {
        additions.extend(plan_amine_hydrogens(structure, residue, removals));
    }

    // Deterministic append order regardless of set iteration order.
    additions.sort_by(|a, b| {
        a.partner
            .cmp(&b.partner)
            .then_with(|| a.info.atom_name.cmp(&b.info.atom_name))
    });
    additions
}

/// Plans the cap hydrogen that turns an open C side into an aldehyde, or `None`
/// when it is already capped or the backbone is incomplete.
fn plan_aldehyde_cap(
    structure: &StructureIndex,
    residue: &ResidueKey,
    removals: &HashSet<usize>,
) -> Option<PendingHydrogen> {
    if structure
        .find(residue, ALDEHYDE_HYDROGEN, removals)
        .is_some()
    {
        return None;
    }
    let c = structure.find(residue, "C", removals)?;
    let o = structure.find(residue, "O", removals)?;
    let ca = structure.find(residue, "CA", removals)?;

    Some(PendingHydrogen {
        position: place_aldehyde_hydrogen(
            structure.atoms[c].position,
            structure.atoms[o].position,
            structure.atoms[ca].position,
        ),
        info: structure.derive_info(c, ALDEHYDE_HYDROGEN),
        partner: c,
    })
}

/// Plans the hydrogens an open N side is missing to reach a neutral amine.
fn plan_amine_hydrogens(
    structure: &StructureIndex,
    residue: &ResidueKey,
    removals: &HashSet<usize>,
) -> Vec<PendingHydrogen> {
    let Some(n) = structure.find(residue, "N", removals) else {
        return Vec::new();
    };
    let (heavy, hydrogens) = structure.split_neighbors(n, removals);
    let wanted = AMINE_SIGMA_BONDS.saturating_sub(heavy.len());

    // Each hydrogen is placed against the substituents already bonded plus the
    // ones planned so far, so two built on the same nitrogen cannot coincide.
    let mut planned = Vec::new();
    for slot in hydrogens.len()..wanted {
        let Some(position) = place_amine_hydrogen(structure, n, &heavy, &hydrogens, &planned)
        else {
            break;
        };
        planned.push(PendingHydrogen {
            // Name the group the way a built N-terminus is named, so the
            // hydrogen already present and the new one read as one pair.
            info: structure.derive_info(n, &format!("H{}", slot + 1)),
            position,
            partner: n,
        });
    }
    planned
}

/// Places the aldehyde cap hydrogen on `c`.
///
/// It sits in the `ca`–`c`–`o` plane, [`ALDEHYDE_ANGLE_DEG`] from the C=O bond
/// and on the far side from `ca`, at [`ALDEHYDE_BOND_LENGTH`]. A degenerate
/// triple that leaves the plane undefined falls back to the direction opposite
/// the other two substituents.
fn place_aldehyde_hydrogen(c: [f64; 3], o: [f64; 3], ca: [f64; 3]) -> [f64; 3] {
    let to_o = normalize(sub(o, c));
    let to_ca = normalize(sub(ca, c));

    let dir = match normalize_checked(cross(to_o, to_ca)) {
        Some(axis) => {
            let theta = ALDEHYDE_ANGLE_DEG.to_radians();
            let forward = rotate_about(to_o, axis, theta);
            let backward = rotate_about(to_o, axis, -theta);
            if dot(forward, to_ca) <= dot(backward, to_ca) {
                forward
            } else {
                backward
            }
        }
        None => normalize(scale(add(to_o, to_ca), -1.0)),
    };

    add(c, scale(dir, ALDEHYDE_BOND_LENGTH))
}

/// Places one more hydrogen on the amine nitrogen `n`, at
/// [`AMINE_BOND_LENGTH`] and turned away from every substituent already there,
/// including hydrogens planned but not yet appended.
fn place_amine_hydrogen(
    structure: &StructureIndex,
    n: usize,
    heavy: &[usize],
    hydrogens: &[usize],
    planned: &[PendingHydrogen],
) -> Option<[f64; 3]> {
    let origin = structure.atoms[n].position;
    let mut existing: Vec<[f64; 3]> = heavy
        .iter()
        .chain(hydrogens)
        .map(|&i| normalize(sub(structure.atoms[i].position, origin)))
        .collect();
    existing.extend(planned.iter().map(|p| normalize(sub(p.position, origin))));
    let anchor = *existing.first()?;

    // Opposite the sum of what is already bonded, then opened out to the
    // tetrahedral angle from the first substituent so two hydrogens built on one
    // nitrogen cannot coincide.
    let away = normalize_checked(scale(
        existing.iter().fold([0.0; 3], |acc, v| add(acc, *v)),
        -1.0,
    ))
    .unwrap_or_else(|| perpendicular(anchor));

    let axis = normalize_checked(cross(anchor, away)).unwrap_or_else(|| perpendicular(anchor));
    let theta = SP3_ANGLE_DEG.to_radians();
    let dir = rotate_about(anchor, axis, theta);
    let dir = if dot(dir, away) < 0.0 {
        rotate_about(anchor, axis, -theta)
    } else {
        dir
    };

    Some(add(origin, scale(dir, AMINE_BOND_LENGTH)))
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

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize_checked(v: [f64; 3]) -> Option<[f64; 3]> {
    let len = dot(v, v).sqrt();
    (len > 1e-9).then(|| scale(v, 1.0 / len))
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    normalize_checked(v).unwrap_or([1.0, 0.0, 0.0])
}

/// Any unit vector orthogonal to `v`.
fn perpendicular(v: [f64; 3]) -> [f64; 3] {
    let seed = if v[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    normalize(cross(v, seed))
}

/// Rotates `v` about the unit `axis` by `theta` radians (Rodrigues' formula).
fn rotate_about(v: [f64; 3], axis: [f64; 3], theta: f64) -> [f64; 3] {
    let (sin, cos) = theta.sin_cos();
    add(
        add(scale(v, cos), scale(cross(axis, v), sin)),
        scale(axis, dot(axis, v) * (1.0 - cos)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::metadata::{BioMetadata, ResiduePosition};

    fn angle_deg(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let u = normalize(sub(a, b));
        let v = normalize(sub(c, b));
        dot(u, v).clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// Builds a glycine dipeptide with the requested caps already present, so
    /// each transform can be exercised against a known starting point.
    ///
    /// Residue 10 opens on the N side, residue 11 on the C side, and the two are
    /// joined by a peptide bond so neither is mistaken for a break.
    fn dipeptide(n_hydrogens: usize, c_carboxyl: bool) -> System {
        let mut atoms: Vec<(String, Element, [f64; 3], i32)> = Vec::new();
        let mut push = |name: &str, e: Element, p: [f64; 3], r: i32| {
            atoms.push((name.to_string(), e, p, r));
        };

        push("N", Element::N, [0.000, 0.000, 0.000], 10);
        for i in 0..n_hydrogens {
            let a = i as f64 * 2.0;
            push(
                &format!("H{}", i + 1),
                Element::H,
                [-0.5, 0.8 * a.cos(), 0.8 * a.sin()],
                10,
            );
        }
        push("CA", Element::C, [1.450, 0.000, 0.000], 10);
        push("C", Element::C, [2.000, 1.400, 0.000], 10);
        push("O", Element::O, [1.300, 2.400, 0.000], 10);

        push("N", Element::N, [3.330, 1.500, 0.000], 11);
        push("H", Element::H, [3.800, 0.620, 0.000], 11);
        push("CA", Element::C, [4.100, 2.720, 0.000], 11);
        push("C", Element::C, [5.600, 2.500, 0.000], 11);
        push("O", Element::O, [6.200, 3.560, 0.000], 11);
        if c_carboxyl {
            push("OXT", Element::O, [6.250, 1.320, 0.000], 11);
            push("HOXT", Element::H, [7.200, 1.400, 0.000], 11);
        }

        let index = |name: &str, r: i32| {
            atoms
                .iter()
                .position(|(n, _, _, rr)| n == name && *rr == r)
                .unwrap()
        };
        let mut bonds = vec![
            Bond::new(index("N", 10), index("CA", 10), BondOrder::Single),
            Bond::new(index("CA", 10), index("C", 10), BondOrder::Single),
            Bond::new(index("C", 10), index("O", 10), BondOrder::Double),
            // The peptide bond joining the two residues.
            Bond::new(index("C", 10), index("N", 11), BondOrder::Single),
            Bond::new(index("N", 11), index("H", 11), BondOrder::Single),
            Bond::new(index("N", 11), index("CA", 11), BondOrder::Single),
            Bond::new(index("CA", 11), index("C", 11), BondOrder::Single),
            Bond::new(index("C", 11), index("O", 11), BondOrder::Double),
        ];
        for i in 0..n_hydrogens {
            bonds.push(Bond::new(
                index("N", 10),
                index(&format!("H{}", i + 1), 10),
                BondOrder::Single,
            ));
        }
        if c_carboxyl {
            bonds.push(Bond::new(
                index("C", 11),
                index("OXT", 11),
                BondOrder::Single,
            ));
            bonds.push(Bond::new(
                index("OXT", 11),
                index("HOXT", 11),
                BondOrder::Single,
            ));
        }

        let atom_info = atoms
            .iter()
            .map(|(name, _, _, r)| {
                AtomResidueInfo::builder(name, "GLY", *r, "A")
                    .standard_name(Some(StandardResidue::GLY))
                    .category(ResidueCategory::Standard)
                    .position(ResiduePosition::Internal)
                    .build()
            })
            .collect();

        let mut system = System {
            atoms: atoms.iter().map(|(_, e, p, _)| Atom::new(*e, *p)).collect(),
            bonds,
            ..Default::default()
        };
        system.bio_metadata = Some(BioMetadata {
            atom_info,
            target_ph: None,
        });
        system
    }

    fn names(system: &System, residue: i32) -> Vec<String> {
        let meta = system.bio_metadata.as_ref().unwrap();
        let mut out: Vec<String> = meta
            .atom_info
            .iter()
            .filter(|i| i.residue_id == residue)
            .map(|i| i.atom_name.clone())
            .collect();
        out.sort();
        out
    }

    fn position_of(system: &System, residue: i32, name: &str) -> [f64; 3] {
        let meta = system.bio_metadata.as_ref().unwrap();
        let idx = meta
            .atom_info
            .iter()
            .position(|i| i.residue_id == residue && i.atom_name == name)
            .unwrap();
        system.atoms[idx].position
    }

    #[test]
    fn renames_hb_hydrogen_to_mpsim_type() {
        let mut types = vec![HB_HYDROGEN_TYPE.to_string(), "H_".to_string()];
        rename_hb_hydrogens(&mut types);
        assert_eq!(types, vec![MPSIM_HB_HYDROGEN_TYPE, "H_"]);
    }

    #[test]
    fn ammonium_n_terminus_loses_its_surplus_proton() {
        let mut system = dipeptide(3, true);
        normalize_open_ends(&mut system);
        assert_eq!(names(&system, 10), vec!["C", "CA", "H1", "H2", "N", "O"]);
    }

    #[test]
    fn open_n_side_gains_the_hydrogen_it_is_missing() {
        // One hydrogen on the nitrogen is what a chain break leaves behind.
        let mut system = dipeptide(1, true);
        normalize_open_ends(&mut system);
        assert_eq!(names(&system, 10), vec!["C", "CA", "H1", "H2", "N", "O"]);
    }

    #[test]
    fn neutral_n_terminus_is_left_alone() {
        let mut system = dipeptide(2, true);
        let before = system.atoms.len();
        normalize_open_ends(&mut system);
        // OXT and its proton go, the aldehyde cap hydrogen arrives.
        assert_eq!(system.atoms.len(), before - 2 + 1);
        assert_eq!(names(&system, 10), vec!["C", "CA", "H1", "H2", "N", "O"]);
    }

    #[test]
    fn carboxyl_c_terminus_becomes_an_aldehyde() {
        let mut system = dipeptide(2, true);
        normalize_open_ends(&mut system);
        assert_eq!(names(&system, 11), vec!["C", "CA", "H", "HC", "N", "O"]);
    }

    #[test]
    fn open_c_side_without_a_carboxyl_still_gets_its_cap() {
        // A chain break leaves the backbone carbon with only CA and O.
        let mut system = dipeptide(2, false);
        normalize_open_ends(&mut system);
        assert_eq!(names(&system, 11), vec!["C", "CA", "H", "HC", "N", "O"]);
    }

    #[test]
    fn cap_hydrogen_matches_the_reference_geometry() {
        let mut system = dipeptide(2, true);
        normalize_open_ends(&mut system);

        let (c, o, ca, hc) = (
            position_of(&system, 11, "C"),
            position_of(&system, 11, "O"),
            position_of(&system, 11, "CA"),
            position_of(&system, 11, "HC"),
        );

        let bond = dot(sub(hc, c), sub(hc, c)).sqrt();
        assert!((bond - ALDEHYDE_BOND_LENGTH).abs() < 1e-6, "bond {bond}");

        let hco = angle_deg(hc, c, o);
        assert!((hco - ALDEHYDE_ANGLE_DEG).abs() < 1e-6, "H-C=O {hco}");

        // Planar, and on the far side from CA: the three angles close the circle.
        let total = hco + angle_deg(hc, c, ca) + angle_deg(o, c, ca);
        assert!((total - 360.0).abs() < 1e-6, "angle sum {total}");
    }

    #[test]
    fn built_amine_hydrogen_uses_the_reference_bond_length() {
        let mut system = dipeptide(1, true);
        normalize_open_ends(&mut system);

        let n = position_of(&system, 10, "N");
        for name in ["H1", "H2"] {
            let h = position_of(&system, 10, name);
            let bond = dot(sub(h, n), sub(h, n)).sqrt();
            // Only the built hydrogen is placed by us; the input one keeps its
            // own geometry, so check the built one is at the right distance.
            if name == "H2" {
                assert!((bond - AMINE_BOND_LENGTH).abs() < 1e-6, "N-H {bond}");
            }
        }
        // The two hydrogens must not coincide.
        let h1 = position_of(&system, 10, "H1");
        let h2 = position_of(&system, 10, "H2");
        assert!(dot(sub(h1, h2), sub(h1, h2)).sqrt() > 1.0);
    }

    #[test]
    fn interior_residues_are_never_capped() {
        let mut system = dipeptide(2, true);
        normalize_open_ends(&mut system);
        // Residue 10's C is bonded to residue 11's N, so it keeps no cap.
        assert!(!names(&system, 10).contains(&ALDEHYDE_HYDROGEN.to_string()));
        // Residue 11's N is bonded to residue 10's C, so it keeps its single H.
        assert_eq!(
            names(&system, 11)
                .iter()
                .filter(|n| n.starts_with('H'))
                .count(),
            2 // "H" and "HC"
        );
    }

    #[test]
    fn is_idempotent() {
        let mut system = dipeptide(3, true);
        normalize_open_ends(&mut system);
        let once: Vec<String> = (10..=11).flat_map(|r| names(&system, r)).collect();
        normalize_open_ends(&mut system);
        let twice: Vec<String> = (10..=11).flat_map(|r| names(&system, r)).collect();
        assert_eq!(once, twice);
    }

    #[test]
    fn relabel_renames_only_hetero_residues() {
        let mut system = dipeptide(2, true);
        {
            let info = &mut system.bio_metadata.as_mut().unwrap().atom_info;
            info[0].category = ResidueCategory::Hetero;
            info[0].residue_name = "PPI".into();
        }

        relabel_ligand_residues(&mut system);

        let info = &system.bio_metadata.as_ref().unwrap().atom_info;
        assert_eq!(info[0].residue_name, MPSIM_LIGAND_RESIDUE_NAME);
        assert_eq!(info[0].residue_id, MPSIM_LIGAND_RESIDUE_ID);
        assert_eq!(info[1].residue_name, "GLY");
    }

    #[test]
    fn systems_without_metadata_are_untouched() {
        let mut system = System {
            atoms: vec![Atom::new(Element::C, [0.0; 3])],
            ..Default::default()
        };
        normalize_open_ends(&mut system);
        relabel_ligand_residues(&mut system);
        assert_eq!(system.atoms.len(), 1);
    }
}
