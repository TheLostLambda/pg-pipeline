mod found_fragment;
mod found_precursor;
mod index;
mod peaks;

// Standard Library Imports
use std::collections::BTreeMap;

// External Crate Imports
use derive_more::Constructor;

// Local Crate Imports
use crate::{
    NamedIon,
    ms2::peaks::Peaks,
    ordered_floats::{Minutes, Mz},
    scan_kv::{ScanKey, ScanValue},
};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Index(BTreeMap<ScanKey, ScanValue<Peaks>>);

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Constructor)]
pub struct FoundPrecursor<'p, 'n> {
    theoretical: &'p NamedIon<'n>,
    observed_mz: Mz,
    scan_number: usize,
    start_time: Minutes,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Constructor)]
pub struct FoundFragment<'p, 'f, 'n> {
    theoretical_precursor: &'p NamedIon<'n>,
    theoretical_fragment: &'f NamedIon<'n>,
    observed_precursor_mz: Mz,
    observed_fragment_mz: Mz,
    scan_number: usize,
    start_time: Minutes,
}
