mod ms2;
mod named_ion;
mod ordered_floats;
mod ppm_window;
mod scan_kv;

// Standard Library Imports
use std::borrow::Cow;

// External Crate Imports
use thiserror::Error;

// Local Crate Imports
use crate::ordered_floats::Mz;

// Public API ==========================================================================================================

pub use ms2::*;

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
