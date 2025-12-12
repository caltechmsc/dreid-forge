use super::error::Error;
use super::intermediate::{
    IntermediateAngle, IntermediateDihedral, IntermediateImproper, IntermediateSystem,
};
use crate::model::types::BondOrder;
use dreid_typer::{
    Element as TyperElement, GraphBondOrder, MolecularGraph, MolecularTopology, assign_topology,
    assign_topology_with_rules, rules,
};

pub fn assign_atom_types(
    system: &mut IntermediateSystem,
    rules: Option<&str>,
) -> Result<(), Error> {
    let graph = build_molecular_graph(system)?;

    let topology = if let Some(rules_toml) = rules {
        let rules = rules::parse_rules(rules_toml).map_err(|e| Error::RuleParse(e.to_string()))?;
        assign_topology_with_rules(&graph, &rules)?
    } else {
        assign_topology(&graph)?
    };

    apply_topology(system, &topology);

    Ok(())
}

fn apply_topology(system: &mut IntermediateSystem, topology: &MolecularTopology) {
    for (int_atom, topo_atom) in system.atoms.iter_mut().zip(topology.atoms.iter()) {
        int_atom.atom_type = topo_atom.atom_type.clone();
        int_atom.hybridization = topo_atom.hybridization;
    }

    for int_bond in &mut system.bonds {
        let (i, j) = (int_bond.i.min(int_bond.j), int_bond.i.max(int_bond.j));
        if let Some(topo_bond) = topology.bonds.iter().find(|b| b.atom_ids == (i, j)) {
            int_bond.physical_order = Some(topo_bond.order);
        }
    }

    system.angles = topology
        .angles
        .iter()
        .map(|a| IntermediateAngle {
            i: a.atom_ids.0,
            j: a.atom_ids.1,
            k: a.atom_ids.2,
        })
        .collect();

    system.dihedrals = topology
        .propers
        .iter()
        .map(|d| IntermediateDihedral {
            i: d.atom_ids.0,
            j: d.atom_ids.1,
            k: d.atom_ids.2,
            l: d.atom_ids.3,
        })
        .collect();

    system.impropers = topology
        .impropers
        .iter()
        .map(|imp| IntermediateImproper {
            p1: imp.atom_ids.0,
            p2: imp.atom_ids.1,
            center: imp.atom_ids.2,
            p3: imp.atom_ids.3,
        })
        .collect();
}

fn build_molecular_graph(system: &IntermediateSystem) -> Result<MolecularGraph, Error> {
    let mut graph = MolecularGraph::new();

    for atom in &system.atoms {
        let element = convert_element(atom.element)?;
        graph.add_atom(element);
    }

    for bond in &system.bonds {
        let order = bond_order_to_graph_order(bond.order);
        graph
            .add_bond(bond.i, bond.j, order)
            .map_err(|e| Error::AtomTyping(e.to_string()))?;
    }

    Ok(graph)
}

fn convert_element(elem: crate::model::types::Element) -> Result<TyperElement, Error> {
    let symbol = elem.symbol();
    symbol.parse::<TyperElement>().map_err(|_| {
        Error::Conversion(format!(
            "element symbol {} not recognized by dreid-typer",
            symbol
        ))
    })
}

fn bond_order_to_graph_order(order: BondOrder) -> GraphBondOrder {
    match order {
        BondOrder::Single => GraphBondOrder::Single,
        BondOrder::Double => GraphBondOrder::Double,
        BondOrder::Triple => GraphBondOrder::Triple,
        BondOrder::Aromatic => GraphBondOrder::Aromatic,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forge::intermediate::{IntermediateSystem, PhysicalBondOrder};
    use crate::model::atom::Atom;
    use crate::model::system::{Bond, System};
    use crate::model::types::Element;

    fn make_water() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::O, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [0.96, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.24, 0.93, 0.0]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys
    }

    fn make_benzene() -> System {
        let mut sys = System::new();
        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            let x = 1.4 * angle.cos();
            let y = 1.4 * angle.sin();
            sys.atoms.push(Atom::new(Element::C, [x, y, 0.0]));
        }
        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            let x = 2.5 * angle.cos();
            let y = 2.5 * angle.sin();
            sys.atoms.push(Atom::new(Element::H, [x, y, 0.0]));
        }
        for i in 0..6 {
            let j = (i + 1) % 6;
            sys.bonds.push(Bond::new(i, j, BondOrder::Aromatic));
        }
        for i in 0..6 {
            sys.bonds.push(Bond::new(i, i + 6, BondOrder::Single));
        }
        sys
    }

    fn make_ethane() -> System {
        let mut sys = System::new();
        sys.atoms.push(Atom::new(Element::C, [0.0, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::C, [1.54, 0.0, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [-0.36, -0.51, 0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, 1.03, 0.0]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, -0.89]));
        sys.atoms.push(Atom::new(Element::H, [1.90, -0.51, 0.89]));
        sys.bonds.push(Bond::new(0, 1, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 2, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 3, BondOrder::Single));
        sys.bonds.push(Bond::new(0, 4, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 5, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 6, BondOrder::Single));
        sys.bonds.push(Bond::new(1, 7, BondOrder::Single));
        sys
    }

    #[test]
    fn types_water_atoms_correctly() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        assert_eq!(int.atoms[0].atom_type, "O_3");
        assert_eq!(int.atoms[1].atom_type, "H_HB");
        assert_eq!(int.atoms[2].atom_type, "H_HB");
    }

    #[test]
    fn populates_water_topology() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        assert_eq!(int.angles.len(), 1);
        assert!(int.dihedrals.is_empty());
        assert!(int.impropers.is_empty());
    }

    #[test]
    fn types_ethane_atoms_correctly() {
        let ethane = make_ethane();
        let mut int = IntermediateSystem::from_system(&ethane).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        assert_eq!(int.atoms[0].atom_type, "C_3");
        assert_eq!(int.atoms[1].atom_type, "C_3");
        for i in 2..8 {
            assert_eq!(int.atoms[i].atom_type, "H_");
        }
    }

    #[test]
    fn populates_ethane_topology() {
        let ethane = make_ethane();
        let mut int = IntermediateSystem::from_system(&ethane).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        assert!(!int.angles.is_empty());
        assert_eq!(int.dihedrals.len(), 9);
        assert!(int.impropers.is_empty());
    }

    #[test]
    fn types_benzene_carbons_as_resonant() {
        let benzene = make_benzene();
        let mut int = IntermediateSystem::from_system(&benzene).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        for i in 0..6 {
            assert_eq!(int.atoms[i].atom_type, "C_R");
        }
        for i in 6..12 {
            assert_eq!(int.atoms[i].atom_type, "H_");
        }
    }

    #[test]
    fn sets_aromatic_bonds_to_resonant() {
        let benzene = make_benzene();
        let mut int = IntermediateSystem::from_system(&benzene).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        for bond in &int.bonds {
            if bond.i < 6 && bond.j < 6 {
                assert_eq!(bond.physical_order, Some(PhysicalBondOrder::Resonant));
            }
        }
    }

    #[test]
    fn generates_impropers_for_planar_centers() {
        let benzene = make_benzene();
        let mut int = IntermediateSystem::from_system(&benzene).unwrap();
        assign_atom_types(&mut int, None).unwrap();

        assert_eq!(int.impropers.len(), 6);
    }

    #[test]
    fn errors_on_invalid_custom_rules() {
        let water = make_water();
        let mut int = IntermediateSystem::from_system(&water).unwrap();

        let invalid_rules = "not valid [[[ toml rules";
        let result = assign_atom_types(&mut int, Some(invalid_rules));
        assert!(matches!(result, Err(Error::RuleParse(_))));
    }

    #[test]
    fn element_conversion_common_elements() {
        assert_eq!(convert_element(Element::C).unwrap(), TyperElement::C);
        assert_eq!(convert_element(Element::N).unwrap(), TyperElement::N);
        assert_eq!(convert_element(Element::O).unwrap(), TyperElement::O);
        assert_eq!(convert_element(Element::H).unwrap(), TyperElement::H);
        assert_eq!(convert_element(Element::S).unwrap(), TyperElement::S);
        assert_eq!(convert_element(Element::P).unwrap(), TyperElement::P);
    }

    #[test]
    fn bond_order_conversion_all_variants() {
        assert_eq!(
            bond_order_to_graph_order(BondOrder::Single),
            GraphBondOrder::Single
        );
        assert_eq!(
            bond_order_to_graph_order(BondOrder::Double),
            GraphBondOrder::Double
        );
        assert_eq!(
            bond_order_to_graph_order(BondOrder::Triple),
            GraphBondOrder::Triple
        );
        assert_eq!(
            bond_order_to_graph_order(BondOrder::Aromatic),
            GraphBondOrder::Aromatic
        );
    }
}
