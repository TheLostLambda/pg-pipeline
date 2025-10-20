mod found_fragment;
mod found_precursor;
mod ms2_index;
mod named_ion;
mod ordered_floats;
mod peaks;
mod ppm_window;
mod scan_kv;

// Standard Library Imports
use std::{borrow::Cow, collections::BTreeMap};

// External Crate Imports
use derive_more::Constructor;
use thiserror::Error;

// Local Crate Imports
use crate::{
    ordered_floats::{Minutes, Mz},
    scan_kv::{ScanKey, ScanValue},
};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Ms2Index(BTreeMap<ScanKey, ScanValue>);

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

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct NamedIon<'n> {
    name: Cow<'n, str>,
    mz: Mz,
}

type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("MS2 spectrum was missing precursor ion information")]
    MissingPrecursor,
    #[error("failed to find centroided peak data")]
    UncentroidedData,
    #[error("failed to decompress gzipped bytes")]
    GzipError,
}
