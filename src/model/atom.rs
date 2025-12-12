//! Atomic representation for molecular systems.
//!
//! Provides the [`Atom`] struct that stores the chemical identity
//! and Cartesian coordinates of an atom within a molecular structure.

use super::types::Element;

/// A single atom with element type and 3D position.
///
/// Represents the minimal information needed to define an atom
/// in a molecular system: its chemical element and Cartesian
/// coordinates in Ångströms.
///
/// # Examples
///
/// ```
/// use dreid_forge::{Atom, Element};
///
/// // Create a carbon atom at the origin
/// let carbon = Atom::new(Element::C, [0.0, 0.0, 0.0]);
///
/// assert_eq!(carbon.element, Element::C);
/// assert_eq!(carbon.position, [0.0, 0.0, 0.0]);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Chemical element of this atom.
    pub element: Element,
    /// Cartesian coordinates [x, y, z] in Ångströms.
    pub position: [f64; 3],
}

impl Atom {
    /// Creates a new atom with the specified element and position.
    ///
    /// # Arguments
    ///
    /// * `element` — Chemical element of the atom
    /// * `position` — Cartesian coordinates [x, y, z] in Ångströms
    ///
    /// # Returns
    ///
    /// A new [`Atom`] instance with the given properties.
    pub fn new(element: Element, position: [f64; 3]) -> Self {
        Self { element, position }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::Element;

    #[test]
    fn atom_new_and_fields() {
        let pos = [0.0_f64, 1.5, -2.25];
        let a = Atom::new(Element::C, pos);
        assert_eq!(a.element, Element::C);
        assert_eq!(a.position, pos);
    }

    #[test]
    fn atom_clone_and_eq() {
        let a = Atom::new(Element::O, [1.0, 2.0, 3.0]);
        let b = a.clone();
        assert_eq!(a, b);
        let mut c = b.clone();
        c.position[0] = 9.0;
        assert_ne!(a, c);
    }

    #[test]
    fn atom_debug_contains_fields() {
        let a = Atom::new(Element::Na, [0.0, 0.0, 0.0]);
        let dbg = format!("{:?}", a);
        assert!(dbg.contains("Atom"));
        assert!(dbg.contains("Na"));
        assert!(dbg.contains("position"));
    }

    #[test]
    fn atom_position_precision_preserved() {
        let pos = [1e-12_f64, -1e6, 3.141592653589793];
        let a = Atom::new(Element::H, pos);
        assert_eq!(a.position[0], pos[0]);
        assert_eq!(a.position[1], pos[1]);
        assert_eq!(a.position[2], pos[2]);
    }
}
