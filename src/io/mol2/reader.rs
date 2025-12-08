use crate::io::{Format, error::Error, util};
use crate::model::{
    atom::Atom,
    system::{Bond, System},
    types::BondOrder,
};
use std::collections::HashMap;
use std::io::BufRead;

pub fn read<R: BufRead>(reader: R) -> Result<System, Error> {
    let lines = collect_lines(reader)?;

    let mol_idx = find_section(&lines, "@<TRIPOS>MOLECULE")
        .ok_or_else(|| Error::parse(Format::Mol2, 1, "missing @<TRIPOS>MOLECULE section"))?;

    let mut cursor = mol_idx + 1;
    let _name = next_data_line(&lines, &mut cursor).unwrap_or_else(|| (0, String::from("MOL2")));

    let (count_line_no, count_line) = next_data_line(&lines, &mut cursor)
        .ok_or_else(|| Error::parse(Format::Mol2, cursor + 1, "missing counts line"))?;
    let (atom_count, bond_count) = parse_counts(&count_line, count_line_no)?;

    let atom_section = find_section(&lines, "@<TRIPOS>ATOM")
        .ok_or_else(|| Error::parse(Format::Mol2, cursor + 1, "missing @<TRIPOS>ATOM section"))?;
    let bond_section = find_section(&lines, "@<TRIPOS>BOND").unwrap_or(lines.len());

    let (atoms, id_map) = parse_atoms(&lines, atom_section + 1, bond_section, atom_count)?;
    let bonds = parse_bonds(&lines, bond_section + 1, atom_count, bond_count, &id_map)?;

    Ok(System {
        atoms,
        bonds,
        box_vectors: None,
        bio_metadata: None,
    })
}

fn collect_lines<R: BufRead>(reader: R) -> Result<Vec<(usize, String)>, Error> {
    reader
        .lines()
        .enumerate()
        .map(|(i, line)| {
            line.map(|v| (i + 1, v))
                .map_err(|e| Error::Io { source: e })
        })
        .collect()
}

fn find_section(lines: &[(usize, String)], name: &str) -> Option<usize> {
    lines
        .iter()
        .position(|(_, line)| line.trim().eq_ignore_ascii_case(name))
}

fn next_data_line(lines: &[(usize, String)], cursor: &mut usize) -> Option<(usize, String)> {
    while *cursor < lines.len() {
        let (ln, content) = &lines[*cursor];
        *cursor += 1;
        let trimmed = content.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        return Some((*ln, content.clone()));
    }
    None
}

fn parse_counts(line: &str, line_no: usize) -> Result<(usize, usize), Error> {
    let parts: Vec<_> = line.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(Error::parse(
            Format::Mol2,
            line_no,
            "counts line must have at least atom and bond counts",
        ));
    }
    let atoms = parts[0]
        .parse::<usize>()
        .map_err(|_| Error::parse(Format::Mol2, line_no, "invalid atom count in counts line"))?;
    let bonds = parts[1]
        .parse::<usize>()
        .map_err(|_| Error::parse(Format::Mol2, line_no, "invalid bond count in counts line"))?;
    Ok((atoms, bonds))
}

fn parse_atoms(
    lines: &[(usize, String)],
    start: usize,
    end: usize,
    expected: usize,
) -> Result<(Vec<Atom>, HashMap<usize, usize>), Error> {
    let mut atoms = Vec::with_capacity(expected);
    let mut id_map = HashMap::new();

    for idx in 0..expected {
        let line_idx = start + idx;
        if line_idx >= end {
            return Err(Error::parse(
                Format::Mol2,
                lines.last().map(|(ln, _)| *ln).unwrap_or(0),
                "ATOM section ended before expected atom count",
            ));
        }
        let (ln, raw) = &lines[line_idx];
        let parts: Vec<_> = raw.split_whitespace().collect();
        if parts.len() < 6 {
            return Err(Error::parse(Format::Mol2, *ln, "invalid ATOM line"));
        }

        let atom_id = parts[0]
            .parse::<usize>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid atom id in ATOM line"))?;
        let x = parts[2]
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid x coordinate in ATOM line"))?;
        let y = parts[3]
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid y coordinate in ATOM line"))?;
        let z = parts[4]
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid z coordinate in ATOM line"))?;

        let element = util::guess_element_symbol(parts[5])
            .or_else(|| util::guess_element_symbol(parts[1]))
            .ok_or_else(|| Error::parse(Format::Mol2, *ln, "unable to infer element"))?;

        id_map.insert(atom_id, atoms.len());
        atoms.push(Atom::new(element, [x, y, z]));
    }

    Ok((atoms, id_map))
}

fn parse_bonds(
    lines: &[(usize, String)],
    start: usize,
    atom_count: usize,
    expected: usize,
    id_map: &HashMap<usize, usize>,
) -> Result<Vec<Bond>, Error> {
    let mut bonds = Vec::with_capacity(expected);

    for idx in 0..expected {
        let line_idx = start + idx;
        if line_idx >= lines.len() {
            return Err(Error::parse(
                Format::Mol2,
                lines.last().map(|(ln, _)| *ln).unwrap_or(0),
                "BOND section ended before expected bond count",
            ));
        }
        let (ln, raw) = &lines[line_idx];
        let parts: Vec<_> = raw.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(Error::parse(Format::Mol2, *ln, "invalid BOND line"));
        }

        let a1 = parts[1]
            .parse::<usize>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid first atom id in BOND line"))?;
        let a2 = parts[2]
            .parse::<usize>()
            .map_err(|_| Error::parse(Format::Mol2, *ln, "invalid second atom id in BOND line"))?;

        let order: BondOrder = util::bond_order_from_mol2(parts[3])
            .ok_or_else(|| Error::parse(Format::Mol2, *ln, "unsupported bond type in BOND line"))?;

        let i = *id_map
            .get(&a1)
            .ok_or_else(|| Error::parse(Format::Mol2, *ln, "bond references unknown atom id"))?;
        let j = *id_map
            .get(&a2)
            .ok_or_else(|| Error::parse(Format::Mol2, *ln, "bond references unknown atom id"))?;

        if i >= atom_count || j >= atom_count {
            return Err(Error::parse(
                Format::Mol2,
                *ln,
                "bond references atom beyond declared count",
            ));
        }

        bonds.push(Bond::new(i, j, order));
    }

    Ok(bonds)
}
