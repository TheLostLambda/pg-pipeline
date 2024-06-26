use std::fmt::{self, Display, Formatter};

use crate::FunctionalGroup;

use super::polymer_database::FunctionalGroupDescription;

impl<'p> FunctionalGroup<'p> {
    // MISSING: Only `pub(crate)` so that users of `Polychem` are forced to get functional groups from a query method
    // on `Polymer` — making it impossible for users to make up a group that doesn't exist
    #[must_use]
    pub(crate) const fn new(name: &'p str, location: &'p str) -> Self {
        Self { name, location }
    }
}

impl<'p> From<&'p FunctionalGroupDescription> for FunctionalGroup<'p> {
    fn from(value: &'p FunctionalGroupDescription) -> Self {
        FunctionalGroup {
            name: value.name.as_str(),
            location: value.location.as_str(),
        }
    }
}

impl Display for FunctionalGroup<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} at={:?}", self.name, self.location)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        assert_eq!(n_terminal.to_string(), r#""Amino" at="N-Terminal""#);
        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        assert_eq!(c_terminal.to_string(), r#""Carboxyl" at="C-Terminal""#);
    }
}
