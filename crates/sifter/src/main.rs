// Standard Library Imports
use std::{
    fs,
    io::stdout,
    path::{Path, PathBuf},
};

// External Crate Imports
use anyhow::Result;
use clap::Parser;
use csv::StringRecord;
use sifter::{FoundFragment, FoundPrecursor, Ms2Index, NamedIon};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// The .mzML file containing MS/MS data to search
    mzml: PathBuf,
    /// A CSV file describing precursor ions with "Structure" and "M/Z" columns
    #[arg(short, long)]
    precursors: PathBuf,
    /// A CSV file describing fragment ions with "Structure" and "M/Z" columns
    #[arg(short, long)]
    fragments: Option<PathBuf>,
    /// The PPM tolerance to use when matching
    #[arg(short, long, default_value_t = 10.0)]
    tolerance: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let mzml_bytes = fs::read(cli.mzml)?;
    let ms2_index = Ms2Index::from_bytes(mzml_bytes)?;
    let precursors = load_mass_database(cli.precursors)?;

    let found: Found = if let Some(fragments) = cli.fragments {
        let fragments = load_mass_database(fragments)?;
        let precursors_and_fragments: Vec<_> = precursors
            .into_iter()
            .map(|precursor| (precursor, &fragments))
            .collect();

        ms2_index
            .find_fragments(&precursors_and_fragments, cli.tolerance)
            .collect()
    } else {
        ms2_index
            .find_precursors(&precursors, cli.tolerance)
            .collect()
    };

    let mut csv_writer = csv::Writer::from_writer(stdout());

    csv_writer.write_record(&found.header())?;
    for record in found.records() {
        csv_writer.write_record(record)?;
    }

    Ok(())
}

fn load_mass_database(path: impl AsRef<Path>) -> Result<Vec<NamedIon<'static>>> {
    csv::Reader::from_path(path)?
        .records()
        .map(|structure_mz_pair| {
            let structure_mz_pair = structure_mz_pair?;
            Ok(NamedIon::new(
                structure_mz_pair[0].to_owned(),
                structure_mz_pair[1].parse::<f64>()?,
            ))
        })
        .collect()
}

enum Found {
    Precursors(Vec<StringRecord>),
    Fragments(Vec<StringRecord>),
}

impl Found {
    fn header(&self) -> StringRecord {
        let header = match self {
            Self::Precursors(_) => [
                "Precursor",
                "Theoretical M/Z (Da)",
                "Observed M/Z (Da)",
                "Scan Number",
                "Start Time (min)",
            ]
            .as_slice(),
            Self::Fragments(_) => [
                "Precursor",
                "Theoretical Precursor M/Z (Da)",
                "Fragment",
                "Theoretical Fragment M/Z (Da)",
                "Observed Precursor M/Z (Da)",
                "Observed Fragment M/Z (Da)",
                "Scan Number",
                "Start Time (min)",
            ]
            .as_slice(),
        };
        StringRecord::from(header)
    }

    fn records(&self) -> impl Iterator<Item = &StringRecord> {
        let (Self::Precursors(records) | Self::Fragments(records)) = &self;
        records.iter()
    }
}

impl<'p, 'n> FromIterator<FoundPrecursor<'p, 'n>> for Found {
    fn from_iter<T: IntoIterator<Item = FoundPrecursor<'p, 'n>>>(iter: T) -> Self {
        Self::Precursors(
            iter.into_iter()
                .map(|fp| {
                    StringRecord::from(
                        [
                            fp.theoretical_name().to_owned(),
                            fp.theoretical_mz().to_string(),
                            fp.observed_mz().to_string(),
                            fp.scan_number().to_string(),
                            fp.start_time().to_string(),
                        ]
                        .as_slice(),
                    )
                })
                .collect(),
        )
    }
}

impl<'p, 'f, 'n> FromIterator<FoundFragment<'p, 'f, 'n>> for Found {
    fn from_iter<T: IntoIterator<Item = FoundFragment<'p, 'f, 'n>>>(iter: T) -> Self {
        Self::Fragments(
            iter.into_iter()
                .map(|ff| {
                    StringRecord::from(
                        [
                            ff.theoretical_precursor_name().to_owned(),
                            ff.theoretical_precursor_mz().to_string(),
                            ff.theoretical_fragment_name().to_owned(),
                            ff.theoretical_fragment_mz().to_string(),
                            ff.observed_precursor_mz().to_string(),
                            ff.observed_fragment_mz().to_string(),
                            ff.scan_number().to_string(),
                            ff.start_time().to_string(),
                        ]
                        .as_slice(),
                    )
                })
                .collect(),
        )
    }
}
