// NOTE: Remember! I'm prototyping! Don't be afraid to write shit code, the point is to learn (and then re-write all of
// this anyways)!

use std::{borrow::Cow, collections::BTreeMap, io::Read};

use derive_more::Constructor;
use flate2::read::GzDecoder;
use mzdata::{
    io::{DetailLevel, mzml::BufferedMzMLReaderType},
    mzpeaks::CentroidPeak,
    prelude::SpectrumLike,
    spectrum::MultiLayerSpectrum,
};
use ordered_float::OrderedFloat;

use crate::{
    Error, NamedIon, Result,
    ordered_floats::{Minutes, Mz},
    ppm_window::PpmWindow,
    scan_kv::{ScanKey, ScanValue},
};

// FIXME: Move this to `src/ordered_floats.rs`
type Intensity = OrderedFloat<f32>;

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Index(BTreeMap<ScanKey, ScanValue<Intensity>>);

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Constructor)]
pub struct FoundIon<'p, 'n> {
    theoretical: &'p NamedIon<'n>,
    observed_mz: Mz,
    scan_number: usize,
    start_time: Minutes,
    intensity: Intensity,
}

// FIXME: Move this impl to a dedicated module (like the other `Found*` types)!
impl<'p> FoundIon<'p, '_> {
    #[must_use]
    pub fn theoretical_name(&self) -> &'p str {
        self.theoretical.name()
    }

    #[must_use]
    pub fn theoretical_mz(&self) -> f64 {
        self.theoretical.mz()
    }

    #[must_use]
    pub fn observed_mz(&self) -> f64 {
        self.observed_mz.into()
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

    #[must_use]
    pub fn intensity(&self) -> f32 {
        self.intensity.into()
    }
}

// FIXME: Move `from_bytes` to some module that can be shared between MS1 and MS2 indicies
impl Index {
    pub fn from_bytes(bytes: impl AsRef<[u8]>) -> Result<Self> {
        let bytes = Self::decode_bytes_if_compressed(bytes.as_ref())?;

        let spectra = BufferedMzMLReaderType::new_buffered_with_detail_level(
            bytes.as_ref(),
            DetailLevel::Lazy,
        );

        Self::from_spectra(spectra)
    }

    pub fn from_spectra(spectra: impl IntoIterator<Item = MultiLayerSpectrum>) -> Result<Self> {
        spectra
            .into_iter()
            .filter(|spectrum| spectrum.ms_level() == 1)
            .flat_map(|mut spectrum| {
                // NOTE: This is assuming that the mzML file we've been given represents a full MS run — if this file
                // is a slice of a larger run, then the scan `.id()` will differ from this!
                let scan_number = spectrum.index() + 1;
                let start_time = spectrum.start_time();

                spectrum
                    .try_build_centroids()
                    .map_err(|_| Error::UncentroidedData)
                    // FIXME: I need to find a way to properly propagate this error!
                    .unwrap()
                    .iter()
                    .map(move |&CentroidPeak { mz, intensity, .. }| {
                        Ok((
                            ScanKey::new(mz, scan_number),
                            ScanValue::new(start_time, intensity.into()),
                        ))
                    })
                    // FIXME: Get rid of this nonsense! Surely there is a way to do this without collecting?
                    .collect::<Vec<_>>()
            })
            .collect::<Result<_>>()
            .map(Index)
    }

    pub fn find_ions<'p, 'n: 'p>(
        &self,
        ions: impl IntoIterator<Item = &'p NamedIon<'n>>,
        ppm_tolerance: impl Into<f64>,
    ) -> impl Iterator<Item = FoundIon<'p, 'n>> {
        let ppm_tolerance = ppm_tolerance.into();
        ions.into_iter().flat_map(move |theoretical| {
            self.filter_scans(theoretical.mz(), ppm_tolerance).map(
                |(
                    &ScanKey {
                        mz: observed_mz,
                        scan_number,
                    },
                    &ScanValue {
                        start_time,
                        value: intensity,
                    },
                )| {
                    FoundIon::new(theoretical, observed_mz, scan_number, start_time, intensity)
                },
            )
        })
    }

    // FIXME: This can fully be pulled out into a shared module for MS1/2 indicies
    fn decode_bytes_if_compressed(bytes: &[u8]) -> Result<Cow<'_, [u8]>> {
        let mut gz_decoder = GzDecoder::new(bytes);

        if gz_decoder.header().is_some() {
            // NOTE: Our decompressed data should be *at least* as big as the compressed data, so we can use the size
            // of the compressed data to allocate some initial capacity
            let mut decoded_bytes = Vec::with_capacity(bytes.len());
            gz_decoder
                .read_to_end(&mut decoded_bytes)
                .map_err(|_| Error::GzipError)?;
            Ok(Cow::Owned(decoded_bytes))
        } else {
            Ok(Cow::Borrowed(bytes))
        }
    }

    // FIXME: This is more duplicated code from the `Ms2Index`!
    fn filter_scans(
        &self,
        precursor_mz: f64,
        ppm_tolerance: f64,
    ) -> impl Iterator<Item = (&ScanKey, &ScanValue<Intensity>)> {
        self.0
            .range(ScanKey::ppm_window(precursor_mz, ppm_tolerance))
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;

    use super::*;

    const MZML: &[u8] = include_bytes!("../../tests/data/WT (6.7–7.3 min).mzML");
    const MZML_GZ: &[u8] = include_bytes!("../../tests/data/WT (6.7–7.3 min).mzML.gz");

    #[test]
    fn from() {
        let mzml = BufferedMzMLReaderType::new_buffered_with_detail_level(MZML, DetailLevel::Lazy);
        let spectra: Vec<_> = mzml.collect();
        let ms1_index = Index::from_spectra(spectra).unwrap();
        let monomer_matches: Vec<_> = ms1_index.filter_scans(942.414_979, 20.).collect();
        assert_debug_snapshot!(monomer_matches);

        let from_bytes = Index::from_bytes(MZML).unwrap();
        assert_eq!(from_bytes, ms1_index);
        let from_gzipped_bytes = Index::from_bytes(MZML_GZ).unwrap();
        assert_eq!(from_gzipped_bytes, ms1_index);
    }

    #[test]
    fn main() {
        use std::io::Write;
        const FULL_MZML: &[u8] = include_bytes!("/home/tll/Downloads/WT.mzML");
        let monomer_isoforms = [
            ("gm-AEJA(+1p)", 942.414_979),
            ("gm-AEJA(+1p, -C1, +[13C]1)", 943.418333),
            ("gm-AEJA(+1p, -C2, +[13C]2)", 944.421688),
            ("gm-AEJA(+2p)", 471.711_128),
            ("gm-AEJA(+2p, -C1, +[13C]1)", 472.212_805),
            ("gm-AEJA(+2p, -C2, +[13C]2)", 472.714_482),
            ("gm-AEJA(+3p)", 314.809_844),
            ("gm-AEJA(+3p, -C1, +[13C]1)", 315.144_295),
            ("gm-AEJA(+3p, -C2, +[13C]2)", 315.478_747),
        ]
        .map(|(name, mz)| NamedIon::new(name, mz));
        let ms1_index = Index::from_bytes(FULL_MZML).unwrap();
        let mut csv = Vec::new();
        writeln!(
            &mut csv,
            "Ion,Theoretical M/Z (Da),Observed M/Z (Da),Scan Number,Start Time (min),Intensity"
        )
        .unwrap();
        for found_ion in ms1_index.find_ions(&monomer_isoforms, 20) {
            writeln!(
                &mut csv,
                "\"{}\",{},{},{},{},{}",
                found_ion.theoretical_name(),
                found_ion.theoretical_mz(),
                found_ion.observed_mz(),
                found_ion.scan_number(),
                found_ion.start_time(),
                found_ion.intensity(),
            )
            .unwrap();
        }
        std::fs::write(
            "/home/tll/Documents/University/PhD/Wasteland/gm-AEJA XICs.csv",
            csv,
        )
        .unwrap();

        let mean_intensity = ms1_index
            .0
            .values()
            .map(|sv| f32::from(sv.value) as f64)
            .sum::<f64>()
            / ms1_index.0.len() as f64;

        dbg!(mean_intensity);

        // FIXME: Think about moving this into the index constructor? So these stats are computed before I flatten
        // everything?
        let mut rts = ms1_index
            .0
            .iter()
            .map(|(sk, sv)| (sk.scan_number, sv.start_time))
            .collect::<Vec<_>>();

        rts.sort_unstable();
        rts.dedup_by_key(|&mut (sn, _)| sn);

        // FIXME: Will crash if there aren't at least two scans — should probably give some nice error instead...
        let mean_rt_gap = rts.windows(2).map(|w| w[1].1 - w[0].1).collect::<Vec<_>>();
        // .sum::<f64>()
        // / (rts.len() - 1) as f64;

        dbg!(mean_rt_gap.iter().max());

        let mut csv = Vec::new();
        writeln!(&mut csv, "MZ,Scan,RT,Intensity").unwrap();
        for (ScanKey { mz, scan_number }, ScanValue { start_time, value }) in ms1_index.0 {
            writeln!(&mut csv, "{mz},{scan_number},{start_time},{value}",).unwrap();
        }
        std::fs::write(
            "/home/tll/Documents/University/PhD/Wasteland/All WT MS1 Ions.csv",
            csv,
        )
        .unwrap();

        panic!()
    }
}
