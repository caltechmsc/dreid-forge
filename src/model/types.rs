use std::fmt;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("invalid or unsupported element symbol: '{0}'")]
pub struct ParseElementError(String);

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("invalid bond order string: '{0}'")]
pub struct ParseBondOrderError(String);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Element {
    H = 1,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Ac,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Nh,
    Fl,
    Mc,
    Lv,
    Ts,
    Og = 118,
}

impl Element {
    pub fn atomic_mass(&self) -> f64 {
        match self {
            Element::H => 1.008,
            Element::He => 4.0026,
            Element::Li => 6.94,
            Element::Be => 9.0122,
            Element::B => 10.81,
            Element::C => 12.011,
            Element::N => 14.007,
            Element::O => 15.999,
            Element::F => 18.998,
            Element::Ne => 20.18,
            Element::Na => 22.99,
            Element::Mg => 24.305,
            Element::Al => 26.982,
            Element::Si => 28.085,
            Element::P => 30.974,
            Element::S => 32.06,
            Element::Cl => 35.45,
            Element::Ar => 39.948,
            Element::K => 39.098,
            Element::Ca => 40.078,
            Element::Sc => 44.956,
            Element::Ti => 47.867,
            Element::V => 50.942,
            Element::Cr => 51.996,
            Element::Mn => 54.938,
            Element::Fe => 55.845,
            Element::Co => 58.933,
            Element::Ni => 58.693,
            Element::Cu => 63.546,
            Element::Zn => 65.38,
            Element::Ga => 69.723,
            Element::Ge => 72.63,
            Element::As => 74.922,
            Element::Se => 78.971,
            Element::Br => 79.904,
            Element::Kr => 83.798,
            Element::Rb => 85.468,
            Element::Sr => 87.62,
            Element::Y => 88.906,
            Element::Zr => 91.224,
            Element::Nb => 92.906,
            Element::Mo => 95.96,
            Element::Tc => 98.0,
            Element::Ru => 101.07,
            Element::Rh => 102.91,
            Element::Pd => 106.42,
            Element::Ag => 107.87,
            Element::Cd => 112.41,
            Element::In => 114.82,
            Element::Sn => 118.71,
            Element::Sb => 121.76,
            Element::Te => 127.6,
            Element::I => 126.9,
            Element::Xe => 131.29,
            Element::Cs => 132.91,
            Element::Ba => 137.33,
            Element::La => 138.91,
            Element::Ce => 140.12,
            Element::Pr => 140.91,
            Element::Nd => 144.24,
            Element::Pm => 145.0,
            Element::Sm => 150.36,
            Element::Eu => 151.96,
            Element::Gd => 157.25,
            Element::Tb => 158.93,
            Element::Dy => 162.5,
            Element::Ho => 164.93,
            Element::Er => 167.26,
            Element::Tm => 168.93,
            Element::Yb => 173.05,
            Element::Lu => 174.97,
            Element::Hf => 178.49,
            Element::Ta => 180.95,
            Element::W => 183.84,
            Element::Re => 186.21,
            Element::Os => 190.23,
            Element::Ir => 192.22,
            Element::Pt => 195.08,
            Element::Au => 196.97,
            Element::Hg => 200.59,
            Element::Tl => 204.38,
            Element::Pb => 207.2,
            Element::Bi => 208.98,
            Element::Po => 209.0,
            Element::At => 210.0,
            Element::Rn => 222.0,
            Element::Fr => 223.0,
            Element::Ra => 226.0,
            Element::Ac => 227.0,
            Element::Th => 232.04,
            Element::Pa => 231.04,
            Element::U => 238.03,
            Element::Np => 237.0,
            Element::Pu => 244.0,
            Element::Am => 243.0,
            Element::Cm => 247.0,
            Element::Bk => 247.0,
            Element::Cf => 251.0,
            Element::Es => 252.0,
            Element::Fm => 257.0,
            Element::Md => 258.0,
            Element::No => 259.0,
            Element::Lr => 262.0,
            Element::Rf => 267.0,
            Element::Db => 270.0,
            Element::Sg => 271.0,
            Element::Bh => 270.0,
            Element::Hs => 277.0,
            Element::Mt => 276.0,
            Element::Ds => 281.0,
            Element::Rg => 280.0,
            Element::Cn => 285.0,
            Element::Nh => 284.0,
            Element::Fl => 289.0,
            Element::Mc => 288.0,
            Element::Lv => 293.0,
            Element::Ts => 294.0,
            Element::Og => 294.0,
        }
    }

    #[inline]
    pub fn atomic_number(&self) -> u8 {
        *self as u8
    }

    pub fn symbol(&self) -> &'static str {
        match self {
            Element::H => "H",
            Element::He => "He",
            Element::Li => "Li",
            Element::Be => "Be",
            Element::B => "B",
            Element::C => "C",
            Element::N => "N",
            Element::O => "O",
            Element::F => "F",
            Element::Ne => "Ne",
            Element::Na => "Na",
            Element::Mg => "Mg",
            Element::Al => "Al",
            Element::Si => "Si",
            Element::P => "P",
            Element::S => "S",
            Element::Cl => "Cl",
            Element::Ar => "Ar",
            Element::K => "K",
            Element::Ca => "Ca",
            Element::Sc => "Sc",
            Element::Ti => "Ti",
            Element::V => "V",
            Element::Cr => "Cr",
            Element::Mn => "Mn",
            Element::Fe => "Fe",
            Element::Co => "Co",
            Element::Ni => "Ni",
            Element::Cu => "Cu",
            Element::Zn => "Zn",
            Element::Ga => "Ga",
            Element::Ge => "Ge",
            Element::As => "As",
            Element::Se => "Se",
            Element::Br => "Br",
            Element::Kr => "Kr",
            Element::Rb => "Rb",
            Element::Sr => "Sr",
            Element::Y => "Y",
            Element::Zr => "Zr",
            Element::Nb => "Nb",
            Element::Mo => "Mo",
            Element::Tc => "Tc",
            Element::Ru => "Ru",
            Element::Rh => "Rh",
            Element::Pd => "Pd",
            Element::Ag => "Ag",
            Element::Cd => "Cd",
            Element::In => "In",
            Element::Sn => "Sn",
            Element::Sb => "Sb",
            Element::Te => "Te",
            Element::I => "I",
            Element::Xe => "Xe",
            Element::Cs => "Cs",
            Element::Ba => "Ba",
            Element::La => "La",
            Element::Ce => "Ce",
            Element::Pr => "Pr",
            Element::Nd => "Nd",
            Element::Pm => "Pm",
            Element::Sm => "Sm",
            Element::Eu => "Eu",
            Element::Gd => "Gd",
            Element::Tb => "Tb",
            Element::Dy => "Dy",
            Element::Ho => "Ho",
            Element::Er => "Er",
            Element::Tm => "Tm",
            Element::Yb => "Yb",
            Element::Lu => "Lu",
            Element::Hf => "Hf",
            Element::Ta => "Ta",
            Element::W => "W",
            Element::Re => "Re",
            Element::Os => "Os",
            Element::Ir => "Ir",
            Element::Pt => "Pt",
            Element::Au => "Au",
            Element::Hg => "Hg",
            Element::Tl => "Tl",
            Element::Pb => "Pb",
            Element::Bi => "Bi",
            Element::Po => "Po",
            Element::At => "At",
            Element::Rn => "Rn",
            Element::Fr => "Fr",
            Element::Ra => "Ra",
            Element::Ac => "Ac",
            Element::Th => "Th",
            Element::Pa => "Pa",
            Element::U => "U",
            Element::Np => "Np",
            Element::Pu => "Pu",
            Element::Am => "Am",
            Element::Cm => "Cm",
            Element::Bk => "Bk",
            Element::Cf => "Cf",
            Element::Es => "Es",
            Element::Fm => "Fm",
            Element::Md => "Md",
            Element::No => "No",
            Element::Lr => "Lr",
            Element::Rf => "Rf",
            Element::Db => "Db",
            Element::Sg => "Sg",
            Element::Bh => "Bh",
            Element::Hs => "Hs",
            Element::Mt => "Mt",
            Element::Ds => "Ds",
            Element::Rg => "Rg",
            Element::Cn => "Cn",
            Element::Nh => "Nh",
            Element::Fl => "Fl",
            Element::Mc => "Mc",
            Element::Lv => "Lv",
            Element::Ts => "Ts",
            Element::Og => "Og",
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

impl FromStr for Element {
    type Err = ParseElementError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Element::H),
            "He" => Ok(Element::He),
            "Li" => Ok(Element::Li),
            "Be" => Ok(Element::Be),
            "B" => Ok(Element::B),
            "C" => Ok(Element::C),
            "N" => Ok(Element::N),
            "O" => Ok(Element::O),
            "F" => Ok(Element::F),
            "Ne" => Ok(Element::Ne),
            "Na" => Ok(Element::Na),
            "Mg" => Ok(Element::Mg),
            "Al" => Ok(Element::Al),
            "Si" => Ok(Element::Si),
            "P" => Ok(Element::P),
            "S" => Ok(Element::S),
            "Cl" => Ok(Element::Cl),
            "Ar" => Ok(Element::Ar),
            "K" => Ok(Element::K),
            "Ca" => Ok(Element::Ca),
            "Sc" => Ok(Element::Sc),
            "Ti" => Ok(Element::Ti),
            "V" => Ok(Element::V),
            "Cr" => Ok(Element::Cr),
            "Mn" => Ok(Element::Mn),
            "Fe" => Ok(Element::Fe),
            "Co" => Ok(Element::Co),
            "Ni" => Ok(Element::Ni),
            "Cu" => Ok(Element::Cu),
            "Zn" => Ok(Element::Zn),
            "Ga" => Ok(Element::Ga),
            "Ge" => Ok(Element::Ge),
            "As" => Ok(Element::As),
            "Se" => Ok(Element::Se),
            "Br" => Ok(Element::Br),
            "Kr" => Ok(Element::Kr),
            "Rb" => Ok(Element::Rb),
            "Sr" => Ok(Element::Sr),
            "Y" => Ok(Element::Y),
            "Zr" => Ok(Element::Zr),
            "Nb" => Ok(Element::Nb),
            "Mo" => Ok(Element::Mo),
            "Tc" => Ok(Element::Tc),
            "Ru" => Ok(Element::Ru),
            "Rh" => Ok(Element::Rh),
            "Pd" => Ok(Element::Pd),
            "Ag" => Ok(Element::Ag),
            "Cd" => Ok(Element::Cd),
            "In" => Ok(Element::In),
            "Sn" => Ok(Element::Sn),
            "Sb" => Ok(Element::Sb),
            "Te" => Ok(Element::Te),
            "I" => Ok(Element::I),
            "Xe" => Ok(Element::Xe),
            "Cs" => Ok(Element::Cs),
            "Ba" => Ok(Element::Ba),
            "La" => Ok(Element::La),
            "Ce" => Ok(Element::Ce),
            "Pr" => Ok(Element::Pr),
            "Nd" => Ok(Element::Nd),
            "Pm" => Ok(Element::Pm),
            "Sm" => Ok(Element::Sm),
            "Eu" => Ok(Element::Eu),
            "Gd" => Ok(Element::Gd),
            "Tb" => Ok(Element::Tb),
            "Dy" => Ok(Element::Dy),
            "Ho" => Ok(Element::Ho),
            "Er" => Ok(Element::Er),
            "Tm" => Ok(Element::Tm),
            "Yb" => Ok(Element::Yb),
            "Lu" => Ok(Element::Lu),
            "Hf" => Ok(Element::Hf),
            "Ta" => Ok(Element::Ta),
            "W" => Ok(Element::W),
            "Re" => Ok(Element::Re),
            "Os" => Ok(Element::Os),
            "Ir" => Ok(Element::Ir),
            "Pt" => Ok(Element::Pt),
            "Au" => Ok(Element::Au),
            "Hg" => Ok(Element::Hg),
            "Tl" => Ok(Element::Tl),
            "Pb" => Ok(Element::Pb),
            "Bi" => Ok(Element::Bi),
            "Po" => Ok(Element::Po),
            "At" => Ok(Element::At),
            "Rn" => Ok(Element::Rn),
            "Fr" => Ok(Element::Fr),
            "Ra" => Ok(Element::Ra),
            "Ac" => Ok(Element::Ac),
            "Th" => Ok(Element::Th),
            "Pa" => Ok(Element::Pa),
            "U" => Ok(Element::U),
            "Np" => Ok(Element::Np),
            "Pu" => Ok(Element::Pu),
            "Am" => Ok(Element::Am),
            "Cm" => Ok(Element::Cm),
            "Bk" => Ok(Element::Bk),
            "Cf" => Ok(Element::Cf),
            "Es" => Ok(Element::Es),
            "Fm" => Ok(Element::Fm),
            "Md" => Ok(Element::Md),
            "No" => Ok(Element::No),
            "Lr" => Ok(Element::Lr),
            "Rf" => Ok(Element::Rf),
            "Db" => Ok(Element::Db),
            "Sg" => Ok(Element::Sg),
            "Bh" => Ok(Element::Bh),
            "Hs" => Ok(Element::Hs),
            "Mt" => Ok(Element::Mt),
            "Ds" => Ok(Element::Ds),
            "Rg" => Ok(Element::Rg),
            "Cn" => Ok(Element::Cn),
            "Nh" => Ok(Element::Nh),
            "Fl" => Ok(Element::Fl),
            "Mc" => Ok(Element::Mc),
            "Lv" => Ok(Element::Lv),
            "Ts" => Ok(Element::Ts),
            "Og" => Ok(Element::Og),
            _ => Err(ParseElementError(s.to_string())),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    pub fn value(&self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BondOrder::Single => write!(f, "Single"),
            BondOrder::Double => write!(f, "Double"),
            BondOrder::Triple => write!(f, "Triple"),
            BondOrder::Aromatic => write!(f, "Aromatic"),
        }
    }
}

impl FromStr for BondOrder {
    type Err = ParseBondOrderError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "single" | "1" => Ok(BondOrder::Single),
            "double" | "2" => Ok(BondOrder::Double),
            "triple" | "3" => Ok(BondOrder::Triple),
            "aromatic" | "ar" => Ok(BondOrder::Aromatic),
            _ => Err(ParseBondOrderError(s.to_string())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn element_from_str_valid() {
        assert_eq!(Element::from_str("H").unwrap(), Element::H);
        assert_eq!(Element::from_str("He").unwrap(), Element::He);
        assert_eq!(Element::from_str("Fe").unwrap(), Element::Fe);
        assert_eq!(Element::from_str("Og").unwrap(), Element::Og);
    }

    #[test]
    fn element_from_str_invalid_case() {
        let err = Element::from_str("h").unwrap_err();
        let s = format!("{}", err);
        assert_eq!(s, "invalid or unsupported element symbol: 'h'");
    }

    #[test]
    fn element_symbol_display_and_atomic_number() {
        let el = Element::Na;
        assert_eq!(el.symbol(), "Na");
        assert_eq!(el.to_string(), "Na");
        assert_eq!(el.atomic_number(), 11u8);
    }

    #[test]
    fn atomic_mass_values() {
        assert!(approx_eq(Element::H.atomic_mass(), 1.008, 1e-6));
        assert!(approx_eq(Element::C.atomic_mass(), 12.011, 1e-6));
        assert!(approx_eq(Element::Fe.atomic_mass(), 55.845, 1e-6));
        assert!(approx_eq(Element::Og.atomic_mass(), 294.0, 1e-6));
    }

    #[test]
    fn bondorder_from_str_variants() {
        assert_eq!(BondOrder::from_str("single").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("1").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("Double").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("2").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("triple").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("3").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("AR").unwrap(), BondOrder::Aromatic);
        assert_eq!(
            BondOrder::from_str("aromatic").unwrap(),
            BondOrder::Aromatic
        );
    }

    #[test]
    fn bondorder_from_str_invalid() {
        let err = BondOrder::from_str("quad").unwrap_err();
        let s = format!("{}", err);
        assert_eq!(s, "invalid bond order string: 'quad'");
    }

    #[test]
    fn bondorder_value_and_display() {
        assert!(approx_eq(BondOrder::Single.value(), 1.0, 1e-12));
        assert!(approx_eq(BondOrder::Double.value(), 2.0, 1e-12));
        assert!(approx_eq(BondOrder::Triple.value(), 3.0, 1e-12));
        assert!(approx_eq(BondOrder::Aromatic.value(), 1.5, 1e-12));

        assert_eq!(BondOrder::Single.to_string(), "Single");
        assert_eq!(BondOrder::Double.to_string(), "Double");
        assert_eq!(BondOrder::Triple.to_string(), "Triple");
        assert_eq!(BondOrder::Aromatic.to_string(), "Aromatic");
    }
}
