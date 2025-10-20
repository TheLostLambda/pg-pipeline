// Standard Library Imports
use std::{borrow::Cow, io::Read};

// External Crate Imports
use flate2::read::GzDecoder;
use mzdata::{
    MzMLReader,
    io::DetailLevel,
    prelude::{IonProperties, SpectrumLike},
    spectrum::MultiLayerSpectrum,
};

// Local Crate Imports
use crate::{
    FoundFragment, FoundPrecursor, Ms2Index, NamedIon, Result,
    ppm_window::PpmWindow,
    scan_kv::{ScanKey, ScanValue},
};

// Public API ==========================================================================================================

impl Ms2Index {
    pub fn from_bytes(bytes: impl AsRef<[u8]>) -> Result<Self> {
        let bytes = Self::decode_bytes_if_compressed(bytes.as_ref())?;

        // FIXME: I think `mzdata` should have a way to construct an `MzMLReader` *without* a buffer! For now, the
        // buffer capacity is `10,000` bytes, completely arbitrarily
        let spectra = MzMLReader::with_buffer_capacity_and_detail_level(
            bytes.as_ref(),
            10_000,
            DetailLevel::Lazy,
        );

        Self::from_spectra(spectra)
    }

    // PERF: This could take an `IntoParallelIterator` and use `rayon` to speed things up somewhat if you're dealing
    // with a large number of spectra. The overhead, however, isn't worth it for smaller datasets...
    pub fn from_spectra(spectra: impl IntoIterator<Item = MultiLayerSpectrum>) -> Result<Self> {
        spectra
            .into_iter()
            .filter(|spectrum| spectrum.ms_level() == 2)
            .map(|mut spectrum| {
                let precursor = spectrum
                    .precursor()
                    .ok_or("MS2 spectrum was missing precursor ion information")?
                    .mz();
                // NOTE: This is assuming that the mzML file we've been given represents a full MS run — if this file
                // is a slice of a larger run, then the scan `.id()` will differ from this!
                let scan_number = spectrum.index() + 1;
                let start_time = spectrum.start_time();

                let peaks = spectrum
                    .try_build_centroids()
                    .map_err(|_| "failed to find centroided peak data")?
                    .into();

                Ok((
                    ScanKey::new(precursor, scan_number),
                    ScanValue::new(start_time, peaks),
                ))
            })
            .collect::<Result<_>>()
            .map(Ms2Index)
    }

    // PERF: Add a `par_find_precursors()` that takes and returns a parallel iterator (if benchmarks justify it)
    pub fn find_precursors<'p, 'n: 'p>(
        &self,
        precursors: impl IntoIterator<Item = &'p NamedIon<'n>>,
        ppm_tolerance: impl Into<f64>,
    ) -> impl Iterator<Item = FoundPrecursor<'p, 'n>> {
        let ppm_tolerance = ppm_tolerance.into();
        precursors.into_iter().flat_map(move |theoretical| {
            self.filter_scans(theoretical.mz(), ppm_tolerance).map(
                |(
                    &ScanKey {
                        precursor,
                        scan_number,
                    },
                    &ScanValue { start_time, .. },
                )| {
                    FoundPrecursor::new(theoretical, precursor, scan_number, start_time)
                },
            )
        })
    }

    // PERF: Add a `par_find_fragments()` that takes and returns a parallel iterator (if benchmarks justify it)
    pub fn find_fragments<'f, 'p, 'n: 'f + 'p>(
        &self,
        fragments: impl IntoIterator<
            Item = &'p (
                NamedIon<'n>,
                impl IntoIterator<Item = &'f NamedIon<'n>> + Copy + 'p,
            ),
        >,
        ppm_tolerance: impl Into<f64>,
    ) -> impl Iterator<Item = FoundFragment<'p, 'f, 'n>> {
        let ppm_tolerance = ppm_tolerance.into();
        fragments
            .into_iter()
            .flat_map(move |(theoretical_precursor, theoretical_fragments)| {
                theoretical_fragments
                    .into_iter()
                    .flat_map(move |theoretical_fragment| {
                        self.filter_scans(theoretical_precursor.mz(), ppm_tolerance)
                            .flat_map(
                                move |(
                                    &ScanKey {
                                        precursor,
                                        scan_number,
                                    },
                                    &ScanValue {
                                        start_time,
                                        ref peaks,
                                    },
                                )| {
                                    peaks
                                        .find_peaks(theoretical_fragment.mz(), ppm_tolerance)
                                        .map(move |fragment_mz| {
                                            FoundFragment::new(
                                                theoretical_precursor,
                                                theoretical_fragment,
                                                precursor,
                                                fragment_mz,
                                                scan_number,
                                                start_time,
                                            )
                                        })
                                },
                            )
                    })
            })
    }
}

// Private Methods =====================================================================================================

impl Ms2Index {
    fn decode_bytes_if_compressed(bytes: &[u8]) -> Result<Cow<'_, [u8]>> {
        let mut gz_decoder = GzDecoder::new(bytes);

        if gz_decoder.header().is_some() {
            // NOTE: Our decompressed data should be *at least* as big as the compressed data, so we can use the size
            // of the compressed data to allocate some initial capacity
            let mut decoded_bytes = Vec::with_capacity(bytes.len());
            gz_decoder
                .read_to_end(&mut decoded_bytes)
                .map_err(|_| "failed to decompress bytes")?;
            Ok(Cow::Owned(decoded_bytes))
        } else {
            Ok(Cow::Borrowed(bytes))
        }
    }

    fn filter_scans(
        &self,
        precursor_mz: f64,
        ppm_tolerance: f64,
    ) -> impl Iterator<Item = (&ScanKey, &ScanValue)> {
        self.0
            .range(ScanKey::ppm_window(precursor_mz, ppm_tolerance))
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use mzdata::{MzMLReader, io::DetailLevel};

    use super::*;

    const MZML: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML");
    const MZML_GZ: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML.gz");

    #[test]
    fn from() {
        let mzml =
            MzMLReader::with_buffer_capacity_and_detail_level(MZML, 10_000, DetailLevel::Lazy);
        let spectra: Vec<_> = mzml.collect();
        let ms2_index = Ms2Index::from_spectra(spectra).unwrap();
        assert_debug_snapshot!(ms2_index);

        let from_bytes = Ms2Index::from_bytes(MZML).unwrap();
        assert_eq!(from_bytes, ms2_index);
        let from_gzipped_bytes = Ms2Index::from_bytes(MZML_GZ).unwrap();
        assert_eq!(from_gzipped_bytes, ms2_index);
    }

    #[test]
    fn find_precursors() {
        let ms2_index = Ms2Index::from_bytes(MZML).unwrap();
        let monomer_ions = [
            NamedIon::new("g(r)m-AEJA(+p)", 942.414_979),
            NamedIon::new("g(r)m-AEJA(+2p)", 471.711_127_5),
        ];

        let found_precursors: Vec<_> = ms2_index.find_precursors(&monomer_ions, 10).collect();
        assert_debug_snapshot!(found_precursors);
    }

    #[test]
    fn find_fragments() {
        let ms2_index = Ms2Index::from_bytes(MZML).unwrap();
        let monomer_fragments = [
            ("gm(r)-AEJA", 942.414_979),
            ("gm(r)-AEJ", 853.367_300),
            ("m(r)-AEJA", 739.335_606),
            ("gm(r)-AE", 681.282_508),
            ("m(r)-AEJ", 650.287_928),
            ("gm(r)-A", 552.239_915),
            ("gm(r)", 481.202_801),
            ("m(r)-AE", 478.203_135),
            ("AEJA", 462.219_454),
            ("EJA", 391.182_340),
            ("AEJ", 373.171_776),
            ("m(r)-A", 349.160_542),
            ("EJ", 302.134_662),
            ("m(r)", 278.123_428),
            ("JA", 262.139_747),
            ("g", 204.086_649),
            ("AE", 201.086_983),
            ("J", 173.092_069),
            ("E", 130.049_870),
            ("A", 90.054_955),
            ("A", 72.044_390),
        ]
        .map(|(name, mz)| NamedIon::new(name, mz));
        let monomer_fragments = [
            NamedIon::new("g(r)m-AEJA(+p)", 942.414_979),
            NamedIon::new("g(r)m-AEJA(+2p)", 471.711_127_5),
        ]
        .map(|named_ion| (named_ion, &monomer_fragments));

        let found_fragments: Vec<_> = ms2_index.find_fragments(&monomer_fragments, 10).collect();
        assert_debug_snapshot!(found_fragments);
    }
}
