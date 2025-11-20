// Standard Library Imports
use std::ops::RangeInclusive;

// External Crate Imports
use mzdata::mzpeaks::Tolerance;

// Local Crate Imports
use crate::{
    ordered_floats::{Minutes, Mz},
    ppm_window::PpmWindow,
};

// Public API ==========================================================================================================

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct ScanKey {
    pub mz: Mz,
    pub scan_number: usize,
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct ScanValue<V> {
    pub start_time: Minutes,
    pub value: V,
}

impl ScanKey {
    pub fn new(mz: f64, scan_number: usize) -> Self {
        Self {
            mz: Mz::from(mz),
            scan_number,
        }
    }
}

impl<V> ScanValue<V> {
    pub fn new(start_time: f64, value: V) -> Self {
        Self {
            start_time: Minutes::from(start_time),
            value,
        }
    }
}

impl PpmWindow for ScanKey {
    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self::new(min_mz, usize::MIN)..=Self::new(max_mz, usize::MAX)
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::ops::{Bound, RangeBounds};

    use assert_float_eq::assert_float_absolute_eq;

    use super::*;

    #[test]
    fn scan_key_ppm_window() {
        let window = ScanKey::ppm_window(471.711_128, 10.0);
        assert!(window.contains(&ScanKey::new(471.711_128, 42)));

        let Bound::Included(&ScanKey {
            mz: start_mz,
            scan_number: start_scan,
        }) = window.start_bound()
        else {
            panic!("expected inclusive start bound");
        };
        assert_float_absolute_eq!(start_mz.into(), 471.706_411);
        assert_eq!(start_scan, usize::MIN);

        let Bound::Included(&ScanKey {
            mz: end_mz,
            scan_number: end_scan,
        }) = window.end_bound()
        else {
            panic!("expected inclusive end bound");
        };
        assert_float_absolute_eq!(end_mz.into(), 471.715_845);
        assert_eq!(end_scan, usize::MAX);
    }
}
