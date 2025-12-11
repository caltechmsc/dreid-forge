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
