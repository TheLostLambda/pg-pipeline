// Local Crate Imports
use crate::FoundFragment;

// Public API ==========================================================================================================

impl<'p, 'f> FoundFragment<'p, 'f, '_> {
    #[must_use]
    pub fn theoretical_precursor_name(&self) -> &'p str {
        self.theoretical_precursor.name()
    }

    #[must_use]
    pub fn theoretical_fragment_name(&self) -> &'f str {
        self.theoretical_fragment.name()
    }

    #[must_use]
    pub fn theoretical_precursor_mz(&self) -> f64 {
        self.theoretical_precursor.mz()
    }

    #[must_use]
    pub fn theoretical_fragment_mz(&self) -> f64 {
        self.theoretical_fragment.mz()
    }

    #[must_use]
    pub fn observed_precursor_mz(&self) -> f64 {
        self.observed_precursor_mz.into()
    }

    #[must_use]
    pub fn observed_fragment_mz(&self) -> f64 {
        self.observed_fragment_mz.into()
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn scan_number(&self) -> usize {
        self.scan_number
    }

    #[must_use]
    pub fn start_time(&self) -> f64 {
        self.start_time.into()
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;

    use crate::{
        NamedIon,
        ordered_floats::{Minutes, Mz},
    };

    use super::*;

    #[test]
    fn found_fragment_getters() {
        let theoretical_precursor = NamedIon::new("gm-AEJA", 942.41);
        let theoretical_fragment = NamedIon::new("J", 173.09);
        let ff = FoundFragment::new(
            &theoretical_precursor,
            &theoretical_fragment,
            Mz::from(942.4121),
            Mz::from(173.0923),
            44,
            Minutes::from(12.345),
        );

        assert_eq!(ff.theoretical_precursor_name(), "gm-AEJA");
        assert_eq!(ff.theoretical_fragment_name(), "J");
        assert_float_absolute_eq!(ff.theoretical_precursor_mz(), 942.41);
        assert_float_absolute_eq!(ff.theoretical_fragment_mz(), 173.09);
        assert_float_absolute_eq!(ff.observed_precursor_mz(), 942.4121);
        assert_float_absolute_eq!(ff.observed_fragment_mz(), 173.0923);
        assert_eq!(ff.scan_number(), 44);
        assert_float_absolute_eq!(ff.start_time(), 12.345);
    }
}
