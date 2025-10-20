// Standard Library Imports
use std::ops::RangeInclusive;

// External Crate Imports
use mzdata::mzpeaks::Tolerance;
use ordered_float::OrderedFloat;

// Local Crate Imports
use crate::ppm_window::PpmWindow;

// Public API ==========================================================================================================

// NOTE: I've tried with `f32` instead of `f64`, but performance actually got slightly worse. It's not shrinking things
// enough to meaningfully increase cache efficiency (but it did, predictably, save memory on the heap!)
pub type Mz = OrderedFloat<f64>;
pub type Minutes = OrderedFloat<f64>;

impl PpmWindow for Mz {
    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self(min_mz)..=Self(max_mz)
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::ops::{Bound, RangeBounds};

    use assert_float_eq::assert_float_absolute_eq;

    use super::*;

    #[test]
    fn mz_ppm_window() {
        let window = Mz::ppm_window(941.407_703, 7.5);
        assert!(window.contains(&Mz::from(941.407_703)));

        let Bound::Included(&OrderedFloat(start)) = window.start_bound() else {
            panic!("expected inclusive start bound");
        };
        let Bound::Included(&OrderedFloat(end)) = window.end_bound() else {
            panic!("expected inclusive end bound")
        };
        assert_float_absolute_eq!(start, 941.400_642);
        assert_float_absolute_eq!(end, 941.414_764);
    }
}
