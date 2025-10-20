use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use sifter::{Ms2Index, NamedIon};

const MZML: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML");
const MZML_GZ: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML.gz");

#[inline]
fn from_mzml() -> Ms2Index {
    Ms2Index::from_bytes(MZML).unwrap()
}

#[inline]
fn from_mzml_gz() -> Ms2Index {
    Ms2Index::from_bytes(MZML_GZ).unwrap()
}

#[inline]
fn find_fragments(ms2_index: &Ms2Index) {
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

    let _found_fragments: Vec<_> = ms2_index.find_fragments(&monomer_fragments, 10).collect();
}

fn build_ms2_index(c: &mut Criterion) {
    c.bench_function("from_mzml", |b| b.iter(|| black_box(from_mzml())));
    c.bench_function("from_mzml_gz", |b| b.iter(|| black_box(from_mzml_gz())));
}

fn search_ms2_index(c: &mut Criterion) {
    let ms2_index = from_mzml();
    c.bench_function("find_fragments", |b| {
        b.iter(|| black_box(find_fragments(&ms2_index)))
    });
}

criterion_group!(benches, build_ms2_index, search_ms2_index);
criterion_main!(benches);
