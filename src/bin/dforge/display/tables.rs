use std::collections::HashMap;
use std::io::{self, Write};

use dreid_forge::{Element, ForgedSystem, System};

use crate::util::text::truncate;

const INDENT: &str = "      ";

const BOX_INNER_WIDTH: usize = 62;
const SAFE_TABLE_WIDTH: usize = BOX_INNER_WIDTH - INDENT.len();

pub fn print_structure_info(system: &System) {
    let stderr = io::stderr();
    let mut out = stderr.lock();

    let mut rows = vec![
        ("Total Atoms", format!("{}", system.atom_count())),
        ("Total Bonds", format!("{}", system.bond_count())),
    ];

    if let Some(box_vectors) = &system.box_vectors {
        let a = vec_len(&box_vectors[0]);
        let b = vec_len(&box_vectors[1]);
        let c = vec_len(&box_vectors[2]);
        rows.push(("Box (Å)", format!("{:.1} × {:.1} × {:.1}", a, b, c)));

        let (alpha, beta, gamma) = calc_angles(box_vectors);
        rows.push((
            "Angles (α β γ)",
            format!("{:.1}° {:.1}° {:.1}°", alpha, beta, gamma),
        ));
    }

    print_kv_table(&mut out, "Structure Summary", &rows);
}

pub fn print_atom_types(system: &System, forged: Option<&ForgedSystem>) {
    let stderr = io::stderr();
    let mut out = stderr.lock();

    if let Some(forged) = forged {
        print_type_distribution(&mut out, forged);
    } else {
        print_element_distribution(&mut out, system);
    }
}

fn print_type_distribution(out: &mut impl Write, forged: &ForgedSystem) {
    let mut type_counts: HashMap<usize, usize> = HashMap::new();
    for param in &forged.atom_properties {
        *type_counts.entry(param.type_idx).or_insert(0) += 1;
    }

    let total = forged.atom_properties.len();
    let mut sorted: Vec<_> = type_counts
        .into_iter()
        .map(|(idx, count)| {
            let name = forged
                .atom_types
                .get(idx)
                .cloned()
                .unwrap_or_else(|| format!("Type{}", idx));
            (name, count)
        })
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    print_distribution_table(out, "Atom Type Distribution", &sorted, total);
}

fn print_element_distribution(out: &mut impl Write, system: &System) {
    let mut element_counts: HashMap<Element, usize> = HashMap::new();
    for atom in &system.atoms {
        *element_counts.entry(atom.element).or_insert(0) += 1;
    }

    let total = system.atoms.len();
    let mut sorted: Vec<_> = element_counts
        .into_iter()
        .map(|(e, c)| (format!("{:?}", e), c))
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    print_distribution_table(out, "Element Distribution", &sorted, total);
}

fn print_distribution_table(
    out: &mut impl Write,
    title: &str,
    data: &[(String, usize)],
    total: usize,
) {
    let name_w = 10usize;
    let count_w = 8usize;
    let sep_overhead = 6;
    let dist_w = SAFE_TABLE_WIDTH.saturating_sub(name_w + count_w + sep_overhead);
    let max_bar_width = dist_w.saturating_sub(8).min(20);

    let _ = writeln!(
        out,
        "{}┌─ {} ─┐",
        INDENT,
        truncate(title, SAFE_TABLE_WIDTH - 6)
    );
    let _ = writeln!(
        out,
        "{}┌{name_line}┬{count_line}┬{dist_line}┐",
        INDENT,
        name_line = "─".repeat(name_w + 2),
        count_line = "─".repeat(count_w + 2),
        dist_line = "─".repeat(dist_w + 2)
    );
    let _ = writeln!(
        out,
        "{}│ {:<name_w$} │ {:>count_w$} │ {:<dist_w$} │",
        INDENT,
        "Type",
        "Count",
        "Distribution",
        name_w = name_w,
        count_w = count_w,
        dist_w = dist_w
    );
    let _ = writeln!(
        out,
        "{}├{name_line}┼{count_line}┼{dist_line}┤",
        INDENT,
        name_line = "─".repeat(name_w + 2),
        count_line = "─".repeat(count_w + 2),
        dist_line = "─".repeat(dist_w + 2)
    );

    for (name, count) in data.iter().take(15) {
        let pct = (*count as f64 / total as f64) * 100.0;
        let bar = make_bar(pct, max_bar_width);
        let name_s = truncate(name, name_w);
        let dist_cell = format!("{}  {:>5.1}%", bar, pct);
        let _ = writeln!(
            out,
            "{}│ {:<name_w$} │ {:>count_w$} │ {:<dist_w$} │",
            INDENT,
            name_s,
            count,
            dist_cell,
            name_w = name_w,
            count_w = count_w,
            dist_w = dist_w
        );
    }

    if data.len() > 15 {
        let _ = writeln!(
            out,
            "{}│ {:<name_w$} │ {:>count_w$} │ {:<dist_w$} │",
            INDENT,
            "...",
            "...",
            format!("({} more types)", data.len() - 15),
            name_w = name_w,
            count_w = count_w,
            dist_w = dist_w
        );
    }

    let _ = writeln!(
        out,
        "{}└{name_line}┴{count_line}┴{dist_line}┘",
        INDENT,
        name_line = "─".repeat(name_w + 2),
        count_line = "─".repeat(count_w + 2),
        dist_line = "─".repeat(dist_w + 2)
    );
}

pub fn print_parameters(forged: &ForgedSystem, bond_type: &str, angle_type: &str, vdw_type: &str) {
    let stderr = io::stderr();
    let mut out = stderr.lock();

    let unique_types: std::collections::HashSet<_> = forged.atom_types.iter().collect();

    let rows: Vec<(String, String, &str)> = vec![
        (
            "Atom Types".to_string(),
            format!("{}", unique_types.len()),
            "unique",
        ),
        (
            format!("Bonds ({})", bond_type),
            format!("{}", forged.potentials.bonds.len()),
            "terms",
        ),
        (
            format!("Angles ({})", angle_type),
            format!("{}", forged.potentials.angles.len()),
            "terms",
        ),
        (
            "Dihedrals".to_string(),
            format!("{}", forged.potentials.dihedrals.len()),
            "terms",
        ),
        (
            "Impropers".to_string(),
            format!("{}", forged.potentials.impropers.len()),
            "terms",
        ),
        (
            format!("VdW ({})", vdw_type),
            format!("{}", forged.potentials.vdw_pairs.len()),
            "pairs",
        ),
        (
            "H-Bond Pairs".to_string(),
            format!("{}", forged.potentials.h_bonds.len()),
            "pairs",
        ),
    ];

    let _ = writeln!(out, "{}┌─ Force Field Parameters ─┐", INDENT);
    let _ = writeln!(out, "{}┌─────────────────────┬────────┬────────┐", INDENT);
    let _ = writeln!(out, "{}│ Category            │  Count │ Type   │", INDENT);
    let _ = writeln!(out, "{}├─────────────────────┼────────┼────────┤", INDENT);

    for (cat, count, typ) in &rows {
        let _ = writeln!(out, "{}│ {:<19} │ {:>6} │ {:<6} │", INDENT, cat, count, typ);
    }

    let _ = writeln!(out, "{}└─────────────────────┴────────┴────────┘", INDENT);
}

pub fn print_chain_breakdown(system: &System) {
    let bio_meta = match &system.bio_metadata {
        Some(m) => m,
        None => return,
    };

    let mut chain_residues: HashMap<char, std::collections::HashSet<(i32, char)>> = HashMap::new();
    let mut chain_atoms: HashMap<char, usize> = HashMap::new();

    for info in &bio_meta.atom_info {
        chain_residues
            .entry(info.chain_id)
            .or_default()
            .insert((info.residue_id, info.insertion_code));
        *chain_atoms.entry(info.chain_id).or_insert(0) += 1;
    }

    if chain_atoms.is_empty() {
        return;
    }

    let mut sorted: Vec<_> = chain_atoms.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));

    let stderr = io::stderr();
    let mut out = stderr.lock();

    let chain_w = 7usize;
    let residues_w = 10usize;
    let sep_overhead = 6;
    let atoms_w = SAFE_TABLE_WIDTH.saturating_sub(chain_w + residues_w + sep_overhead);

    let _ = writeln!(out, "{}┌─ Chain Breakdown ─┐", INDENT);
    let _ = writeln!(
        out,
        "{}┌{c_line}┬{r_line}┬{a_line}┐",
        INDENT,
        c_line = "─".repeat(chain_w + 2),
        r_line = "─".repeat(residues_w + 2),
        a_line = "─".repeat(atoms_w + 2)
    );
    let _ = writeln!(
        out,
        "{}│ {:<chain_w$} │ {:>residues_w$} │ {:>atoms_w$} │",
        INDENT,
        "Chain",
        "Residues",
        "Atoms",
        chain_w = chain_w,
        residues_w = residues_w,
        atoms_w = atoms_w
    );
    let _ = writeln!(
        out,
        "{}├{c_line}┼{r_line}┼{a_line}┤",
        INDENT,
        c_line = "─".repeat(chain_w + 2),
        r_line = "─".repeat(residues_w + 2),
        a_line = "─".repeat(atoms_w + 2)
    );

    for (chain, atoms) in &sorted {
        let residues = chain_residues.get(chain).map(|s| s.len()).unwrap_or(0);
        let _ = writeln!(
            out,
            "{}│ {:<chain_w$} │ {:>residues_w$} │ {:>atoms_w$} │",
            INDENT,
            chain,
            residues,
            atoms,
            chain_w = chain_w,
            residues_w = residues_w,
            atoms_w = atoms_w
        );
    }

    let _ = writeln!(
        out,
        "{}└{c_line}┴{r_line}┴{a_line}┘",
        INDENT,
        c_line = "─".repeat(chain_w + 2),
        r_line = "─".repeat(residues_w + 2),
        a_line = "─".repeat(atoms_w + 2)
    );
}

fn print_kv_table(out: &mut impl Write, title: &str, rows: &[(&str, String)]) {
    let key_w = 16usize;
    let sep_overhead = 6;
    let val_w = SAFE_TABLE_WIDTH.saturating_sub(key_w + sep_overhead);

    let _ = writeln!(
        out,
        "{}┌─ {} ─┐",
        INDENT,
        truncate(title, SAFE_TABLE_WIDTH - 6)
    );
    let _ = writeln!(
        out,
        "{}┌{k_line}┬{v_line}┐",
        INDENT,
        k_line = "─".repeat(key_w + 2),
        v_line = "─".repeat(val_w + 2)
    );
    let _ = writeln!(
        out,
        "{}│ {:<key_w$} │ {:>val_w$} │",
        INDENT,
        "Metric",
        "Value",
        key_w = key_w,
        val_w = val_w
    );
    let _ = writeln!(
        out,
        "{}├{k_line}┼{v_line}┤",
        INDENT,
        k_line = "─".repeat(key_w + 2),
        v_line = "─".repeat(val_w + 2)
    );

    for (key, val) in rows {
        let _ = writeln!(
            out,
            "{}│ {:<key_w$} │ {:>val_w$} │",
            INDENT,
            truncate(key, key_w),
            truncate(val, val_w),
            key_w = key_w,
            val_w = val_w
        );
    }

    let _ = writeln!(
        out,
        "{}└{k_line}┴{v_line}┘",
        INDENT,
        k_line = "─".repeat(key_w + 2),
        v_line = "─".repeat(val_w + 2)
    );
}

fn make_bar(pct: f64, max_width: usize) -> String {
    let filled = ((pct / 100.0) * max_width as f64).round() as usize;
    let empty = max_width.saturating_sub(filled);
    format!("{}{}", "█".repeat(filled), "░".repeat(empty))
}

fn vec_len(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn calc_angles(box_vectors: &[[f64; 3]; 3]) -> (f64, f64, f64) {
    let a = &box_vectors[0];
    let b = &box_vectors[1];
    let c = &box_vectors[2];

    let len_a = vec_len(a);
    let len_b = vec_len(b);
    let len_c = vec_len(c);

    let alpha = ((b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / (len_b * len_c))
        .acos()
        .to_degrees();
    let beta = ((a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / (len_a * len_c))
        .acos()
        .to_degrees();
    let gamma = ((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (len_a * len_b))
        .acos()
        .to_degrees();

    (alpha, beta, gamma)
}
