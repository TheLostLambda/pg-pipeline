// Standard Library Imports
use std::collections::BTreeSet;

// External Crate Imports
use mzdata::mzpeaks::{CentroidPeak, MZPeakSetType};

// Local Crate Imports
use crate::{ordered_floats::Mz, ppm_window::PpmWindow};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Peaks(BTreeSet<Mz>);

impl From<&MZPeakSetType<CentroidPeak>> for Peaks {
    fn from(peak_set: &MZPeakSetType<CentroidPeak>) -> Self {
        Self(peak_set.iter().map(|peak| peak.mz.into()).collect())
    }
}

impl Peaks {
    pub fn find_peaks(&self, mz: f64, ppm: f64) -> impl Iterator<Item = Mz> {
        self.0.range(Mz::ppm_window(mz, ppm)).copied()
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peaks_filter_peaks() {
        let query = 942.412_793_603_516;
        let peaks = Peaks(BTreeSet::from(
            [
                942.413_696_289_063,
                942.4121,
                942.413_452_148_438,
                942.4121,
                942.4119,
                942.413_513_183_594,
            ]
            .map(Mz::from),
        ));
        let filter_peaks = |ppm| -> Vec<_> { peaks.find_peaks(query, ppm).collect() };
        assert_eq!(filter_peaks(0.0), [] as [Mz; 0]);
        assert_eq!(filter_peaks(0.75), [942.4121, 942.413_452_148_438]);
        assert_eq!(
            filter_peaks(10.0),
            [
                942.4119,
                942.4121,
                942.413_452_148_438,
                942.413_513_183_594,
                942.413_696_289_063
            ]
        );
    }
}
