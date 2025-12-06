use super::types::Element;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub position: [f64; 3],
}

impl Atom {
    pub fn new(element: Element, position: [f64; 3]) -> Self {
        Self { element, position }
    }
}
