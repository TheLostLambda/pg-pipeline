use std::{hint::black_box, sync::LazyLock};

use gungraun::{
    Callgrind, Dhat, FlamegraphConfig, LibraryBenchmarkConfig, library_benchmark,
    library_benchmark_group, main,
};
use sifter::{Ms2Index, NamedIon};

const MZML: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML");
const MZML_GZ: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML.gz");

static MS2_INDEX: LazyLock<Ms2Index> = LazyLock::new(|| Ms2Index::from_bytes(MZML).unwrap());

static MONOMER_MS2_FRAGMENTS: LazyLock<[NamedIon<'static>; 21]> = LazyLock::new(|| {
    [
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
    .map(|(name, mz)| NamedIon::new(name, mz))
});

static MONOMER_FRAGMENTS: LazyLock<[(NamedIon<'static>, &'static [NamedIon<'static>; 21]); 2]> =
    LazyLock::new(|| {
        [
            NamedIon::new("g(r)m-AEJA(+p)", 942.414_979),
            NamedIon::new("g(r)m-AEJA(+2p)", 471.711_127_5),
        ]
        .map(|named_ion| (named_ion, &*MONOMER_MS2_FRAGMENTS))
    });

fn force_lazy_statics() -> (
    &'static Ms2Index,
    &'static [(NamedIon<'static>, &'static [NamedIon<'static>; 21]); 2],
) {
    (
        LazyLock::force(&MS2_INDEX),
        LazyLock::force(&MONOMER_FRAGMENTS),
    )
}

#[library_benchmark]
#[bench::mzml(MZML)]
#[bench::mzml_gz(MZML_GZ)]
fn from(bytes: &[u8]) -> Ms2Index {
    Ms2Index::from_bytes(bytes).unwrap()
}

#[library_benchmark]
#[bench::monomer(setup = force_lazy_statics)]
fn find_fragments((ms2_index, fragments): (&Ms2Index, &[(NamedIon, &[NamedIon; 21]); 2])) {
    black_box(ms2_index.find_fragments(fragments, 10).for_each(drop));
}

library_benchmark_group!(
    name = ms2_index;
    benchmarks = from, find_fragments
);

main!(
    config = LibraryBenchmarkConfig::default()
        .tool(Callgrind::default().flamegraph(FlamegraphConfig::default()))
        .tool(Dhat::default());
    library_benchmark_groups = ms2_index
);
