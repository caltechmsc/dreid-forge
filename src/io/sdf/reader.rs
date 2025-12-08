use crate::io::{Format, error::Error, util};
use crate::model::{
    atom::Atom,
    system::{Bond, System},
};
use std::io::BufRead;

pub fn read<R: BufRead>(reader: R) -> Result<System, Error> {
    let lines = collect_first_block(reader)?;
    if lines.len() < 4 {
        return Err(Error::parse(
            Format::Sdf,
            1,
            "SDF block must contain at least a header and counts line",
        ));
    }

    let counts_line_no = lines[3].0;
    let counts_line = &lines[3].1;
    if counts_line.contains("V3000") {
        return Err(Error::parse(
            Format::Sdf,
            counts_line_no,
            "V3000 is not supported",
        ));
    }

    let (atom_count, bond_count) = parse_counts(counts_line, counts_line_no)?;
    let atom_start = 4;
    let bond_start = atom_start + atom_count;

    if lines.len() < bond_start + bond_count {
        return Err(Error::parse(
            Format::Sdf,
            lines.last().map(|(ln, _)| *ln).unwrap_or(counts_line_no),
            "SDF block ended before atoms/bonds were fully specified",
        ));
    }

    let atoms = parse_atoms(&lines[atom_start..atom_start + atom_count])?;
    let bonds = parse_bonds(&lines[bond_start..bond_start + bond_count], atom_count)?;

    Ok(System {
        atoms,
        bonds,
        box_vectors: None,
        bio_metadata: None,
    })
}

fn collect_first_block<R: BufRead>(reader: R) -> Result<Vec<(usize, String)>, Error> {
    let mut lines = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let content = line.map_err(|e| Error::Io { source: e })?;
        let ln = i + 1;
        if content.trim() == "$$$$" && !lines.is_empty() {
            break;
        }
        lines.push((ln, content));
    }
    Ok(lines)
}

fn parse_counts(line: &str, line_no: usize) -> Result<(usize, usize), Error> {
    let tokens: Vec<_> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(Error::parse(
            Format::Sdf,
            line_no,
            "counts line must contain atom and bond counts",
        ));
    }
    let atoms = tokens[0]
        .parse::<usize>()
        .map_err(|_| Error::parse(Format::Sdf, line_no, "invalid atom count"))?;
    let bonds = tokens[1]
        .parse::<usize>()
        .map_err(|_| Error::parse(Format::Sdf, line_no, "invalid bond count"))?;
    Ok((atoms, bonds))
}

fn parse_atoms(lines: &[(usize, String)]) -> Result<Vec<Atom>, Error> {
    let mut atoms = Vec::with_capacity(lines.len());
    for (ln, raw) in lines {
        let padded = format!("{raw:<40}");
        let x = padded[0..10]
            .trim()
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid x coordinate in atom line"))?;
        let y = padded[10..20]
            .trim()
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid y coordinate in atom line"))?;
        let z = padded[20..30]
            .trim()
            .parse::<f64>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid z coordinate in atom line"))?;
        let element_token = padded[31..34].trim();
        let element = util::guess_element_symbol(element_token)
            .ok_or_else(|| Error::parse(Format::Sdf, *ln, "unable to infer element symbol"))?;
        atoms.push(Atom::new(element, [x, y, z]));
    }
    Ok(atoms)
}

fn parse_bonds(lines: &[(usize, String)], atom_count: usize) -> Result<Vec<Bond>, Error> {
    let mut bonds = Vec::with_capacity(lines.len());
    for (ln, raw) in lines {
        let tokens: Vec<_> = raw.split_whitespace().collect();
        if tokens.len() < 3 {
            return Err(Error::parse(Format::Sdf, *ln, "invalid bond line"));
        }

        let a1 = tokens[0]
            .parse::<usize>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid first atom index"))?;
        let a2 = tokens[1]
            .parse::<usize>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid second atom index"))?;
        let order_val = tokens[2]
            .parse::<i32>()
            .map_err(|_| Error::parse(Format::Sdf, *ln, "invalid bond order value"))?;

        let order = util::bond_order_from_ctfile(order_val)
            .ok_or_else(|| Error::parse(Format::Sdf, *ln, "unsupported bond order in bond line"))?;

        if a1 == 0 || a2 == 0 || a1 > atom_count || a2 > atom_count {
            return Err(Error::parse(
                Format::Sdf,
                *ln,
                "bond references atom outside declared range",
            ));
        }

        bonds.push(Bond::new(a1 - 1, a2 - 1, order));
    }
    Ok(bonds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{atom::Atom, types::Element};
    use std::io::Cursor;

    fn sample_sdf() -> String {
        r#"Methanol
        Program
        Comment
        4  3  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            1.0900    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        -0.5400    0.9350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        -0.5400   -0.9350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        1  2  1  0  0  0  0
        1  3  1  0  0  0  0
        1  4  1  0  0  0  0
        M  END
        $$$$
        "#
        .to_string()
    }

    #[test]
    fn parses_first_block_only() {
        let mut multi = sample_sdf();
        multi.push_str(&sample_sdf());
        let system = read(Cursor::new(multi)).expect("parse sdf");
        assert_eq!(system.atom_count(), 4);
        assert_eq!(system.bond_count(), 3);
        assert_eq!(system.atoms[1], Atom::new(Element::O, [1.09, 0.0, 0.0]));
    }

    #[test]
    fn errors_on_v3000() {
        let mut sdf = sample_sdf();
        sdf = sdf.replace("V2000", "V3000");
        let err = read(Cursor::new(sdf)).expect_err("V3000 should be rejected");
        match err {
            Error::Parse { format, .. } => assert_eq!(format, Format::Sdf),
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
