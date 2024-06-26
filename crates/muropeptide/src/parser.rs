use std::cell::RefCell;

use miette::Diagnostic;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, one_of, space0, space1},
    combinator::{cut, map, opt, recognize},
    error::ErrorKind,
    multi::{many0, many1, separated_list1},
    sequence::{delimited, pair, preceded, terminated, tuple},
    IResult,
};
use nom_miette::{map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError};
use polychem::{
    errors::PolychemError,
    parsers::{
        chemical_composition,
        errors::PolychemErrorKind,
        primitives::{count, lowercase, offset_kind, uppercase},
    },
    Count, ModificationId, Polymer, Polymerizer,
};
use thiserror::Error;

use crate::{
    AminoAcid, LateralChain, Monomer, Monosaccharide, Muropeptide, PeptideDirection,
    UnbranchedAminoAcid,
};

// FIXME: Need to think about if these should really live in another KDL config?
const PEPTIDE_BOND: &str = "Pep";
const GLYCOSIDIC_BOND: &str = "Gly";
const STEM_BOND: &str = "Stem";
const NTOC_BOND: &str = "NToC";
const CTON_BOND: &str = "CToN";
const CROSSLINK_BOND: &str = "Link";

// FIXME: A horrible hack that's needed to specify the lifetimes captured by `impl FnMut(...) -> ...` correctly. Once
// Rust 2024 is stabilized, however, this hack can be removed. Keep an eye on:
// https://github.com/rust-lang/rust/issues/117587
// FIXME: Make private again
pub trait Captures<U> {}
impl<T: ?Sized, U> Captures<U> for T {}

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!

/// Muropeptide = Monomer , { Connection , Monomer } , [ Connection ] , [ { " " }- ,
///   ( Modifications , [ { " " }- , Crosslinks ]
///   | Crosslinks , [ { " " }- , Modifications ]
///   ) ] ;
// FIXME: Very very incomplete!
pub fn muropeptide<'z, 'a, 'p, 's>(
    polymerizer: &'z Polymerizer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Muropeptide<'a, 'p>> + Captures<(&'z (), &'a (), &'p ())> {
    move |i| {
        let polymer = RefCell::new(polymerizer.new_polymer());
        // FIXME: Perhaps there is a better way to shorten that `polymer` borrow...
        let (rest, (monomer, _)) = {
            let mut parser = pair(
                monomer(&polymer),
                opt(preceded(space1, modifications(&polymer))),
            );
            parser(i)?
        };

        Ok((
            rest,
            Muropeptide {
                polymer: polymer.into_inner(),
                monomers: vec![monomer],
                connections: Vec::new(),
            },
        ))
    }
}

/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
// FIXME: Make private again
pub fn monomer<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer> + Captures<(&'c (), &'a (), &'p ())> {
    let optional_peptide = opt(preceded(char('-'), cut(peptide(polymer))));
    let glycan_and_peptide = map_res(
        pair(glycan(polymer), optional_peptide),
        |(glycan, peptide)| {
            if let Some(peptide) = peptide {
                // SAFETY: Both the `glycan` and `peptide` parsers ensure at least one residue is present, so `.last()` and
                // `.first()` will never return `None`!
                let donor = *glycan.last().unwrap();
                let acceptor = peptide.first().unwrap().residue;

                polymer
                    .borrow_mut()
                    .bond_residues(STEM_BOND, donor, acceptor)?;
                Ok(Monomer { glycan, peptide })
            } else {
                Ok(Monomer {
                    glycan,
                    peptide: Vec::new(),
                })
            }
        },
    );

    let just_peptide = map(peptide(polymer), |peptide| Monomer {
        glycan: Vec::new(),
        peptide,
    });

    // FIXME: Add a `map_res` wrapping this final parser
    alt((glycan_and_peptide, just_peptide))
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<Monosaccharide>> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = many1(monosaccharide(polymer));
    map_res(parser, |residues| {
        polymer
            .borrow_mut()
            .bond_chain(GLYCOSIDIC_BOND, &residues)?;
        Ok(residues)
    })
}

/// Peptide = { Amino Acid }- ;
fn peptide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AminoAcid>> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = many1(amino_acid(polymer));
    map_res(parser, |residues| {
        let residue_ids = residues.iter().map(|aa| aa.residue);
        polymer.borrow_mut().bond_chain(PEPTIDE_BOND, residue_ids)?;
        Ok(residues)
    })
}

// =

/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Monosaccharide> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = pair(recognize(lowercase), opt(modifications(polymer)));
    map_res(parser, |(abbr, modifications)| {
        let residue = polymer.borrow_mut().new_residue(abbr)?;
        for modification in modifications.into_iter().flatten() {
            polymer
                .borrow_mut()
                .localize_modification(modification, residue)?;
        }

        Ok(residue)
    })
}

// FIXME: Damn... This is messy... Need to sort that out!
/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = pair(unbranched_amino_acid(polymer), opt(lateral_chain(polymer)));
    map_res(parser, |(residue, lateral_chain)| {
        if let Some(LateralChain { direction, peptide }) = &lateral_chain {
            let c_to_n = || -> polychem::Result<_> {
                polymer
                    .borrow_mut()
                    .bond_residues(CTON_BOND, peptide[0], residue)?;
                // FIXME: Replace with `.clone()` then `.reverse()`?
                Ok(peptide.iter().copied().rev().collect())
            };
            let n_to_c = || -> polychem::Result<_> {
                polymer
                    .borrow_mut()
                    .bond_residues(NTOC_BOND, residue, peptide[0])?;
                // FIXME: This feels silly... Maybe better when `bond_chain` takes an `impl IntoIterator`?
                Ok(peptide.clone())
            };

            let chain: Vec<_> = match direction {
                // FIXME: This default should be moved to a configuration file and not be hard-coded!
                // FIXME: Furthermore, this should really return several possible polymers instead of picking one.
                // That might end up being difficult to do inline with this parsing stuff...
                PeptideDirection::Unspecified => c_to_n().or_else(|_| n_to_c())?,
                PeptideDirection::CToN => c_to_n()?,
                PeptideDirection::NToC => n_to_c()?,
            };
            polymer.borrow_mut().bond_chain(PEPTIDE_BOND, &chain)?;
        }
        Ok(AminoAcid {
            residue,
            lateral_chain,
        })
    })
}

// =

/// Modifications = "(" , Any Modification ,
///   { { " " } , "," , { " " } , Any Modification } , ")" ;
fn modifications<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<ModificationId>> + Captures<(&'c (), &'a (), &'p ())> {
    let separator = delimited(space0, char(','), space0);
    delimited(
        char('('),
        separated_list1(separator, any_modification(polymer)),
        char(')'),
    )
}

/// Unbranched Amino Acid = [ lowercase ] , uppercase , [ Modifications ] ;
fn unbranched_amino_acid<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<UnbranchedAminoAcid> + Captures<(&'c (), &'a (), &'p ())> {
    let abbr = recognize(preceded(opt(lowercase), uppercase));
    let parser = pair(abbr, opt(modifications(polymer)));
    map_res(parser, |(abbr, modifications)| {
        let residue = polymer.borrow_mut().new_residue(abbr)?;
        for modification in modifications.into_iter().flatten() {
            polymer
                .borrow_mut()
                .localize_modification(modification, residue)?;
        }

        Ok(residue)
    })
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Lateral Chain = "[" , Peptide Direction , { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain> + Captures<(&'c (), &'a (), &'p ())> {
    let peptide = many1(unbranched_amino_acid(polymer));
    let parser = delimited(char('['), pair(peptide_direction, peptide), char(']'));
    map(parser, |(direction, peptide)| LateralChain {
        direction,
        peptide,
    })
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Peptide Direction = [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ;
fn peptide_direction(i: &str) -> ParseResult<PeptideDirection> {
    map(opt(one_of("<>")), |c| match c {
        Some('<') => PeptideDirection::CToN,
        Some('>') => PeptideDirection::NToC,
        None => PeptideDirection::Unspecified,
        _ => unreachable!(),
    })(i)
}
// =

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

/// Any Modification = Named Modification | Offset Modification
pub fn any_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId> + Captures<(&'c (), &'a (), &'p ())> {
    alt((named_modification(polymer), offset_modification(polymer)))
}

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Named Modification = [ Multiplier ] , Identifier
pub fn named_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = pair(opt(multiplier), identifier);
    map_res(parser, |(multiplier, named_mod)| {
        polymer
            .borrow_mut()
            .new_modification(multiplier.unwrap_or_default(), named_mod)
    })
}

/// Offset Modification = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
pub fn offset_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId> + Captures<(&'c (), &'a (), &'p ())> {
    let chemical_composition = chemical_composition(polymer.borrow().atomic_db());
    let parser = tuple((offset_kind, opt(multiplier), chemical_composition));

    map_res(parser, |(kind, multiplier, composition)| {
        polymer.borrow_mut().new_offset_with_composition(
            kind,
            multiplier.unwrap_or_default(),
            composition,
        )
    })
}

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let mut parser = terminated(count, char('x'));
    // FIXME: Add error handling / reporting!
    parser(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores")]
    ExpectedIdentifier,

    // FIXME: Kill this and merge into the error below!
    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(Box<PolychemError>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] PolychemErrorKind),

    #[diagnostic(help(
        "this is an internal error that you shouldn't ever see! If you have gotten this error, \
        then please report it as a bug!"
    ))]
    #[error("internal `nom` error: {0:?}")]
    NomError(ErrorKind),

    #[diagnostic(help(
        "check the unparsed region for errors, or remove it from the rest of the muropeptide"
    ))]
    #[error("could not interpret the full input as a valid muropeptide structure")]
    Incomplete,
}

impl LabeledErrorKind for MuropeptideErrorKind {
    fn label(&self) -> Option<&'static str> {
        Some(match self {
            // FIXME: Need to add branches for passing labels through the transparent errors?
            Self::Incomplete => "input was valid up until this point",
            Self::NomError(_) => "the region that triggered this bug!",
            _ => return None,
        })
    }
}

// FIXME: Can I get rid of this?
impl<'a> FromExternalError<'a, Box<PolychemError>> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: Box<PolychemError>) -> LabeledParseError<'_, Self> {
        LabeledParseError::new(input, Self::PolychemError(e))
    }
}

impl From<ErrorKind> for MuropeptideErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::Incomplete,
            kind => Self::NomError(kind),
        }
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use polychem::{AtomicDatabase, Charged, Massive, PolymerDatabase, Polymerizer};
    use rust_decimal::Decimal;
    use rust_decimal_macros::dec;

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

    #[test]
    fn test_identifier() {
        // Valid Identifiers
        assert_eq!(identifier("Ac"), Ok(("", "Ac")));
        assert_eq!(identifier("H2O"), Ok(("", "H2O")));
        assert_eq!(identifier("Anh"), Ok(("", "Anh")));
        assert_eq!(identifier("E2E"), Ok(("", "E2E")));
        assert_eq!(identifier("no_way"), Ok(("", "no_way")));
        assert_eq!(identifier("H"), Ok(("", "H")));
        assert_eq!(identifier("p"), Ok(("", "p")));
        // Invalid Identifiers
        assert!(identifier(" H2O").is_err());
        assert!(identifier("1").is_err());
        assert!(identifier("9999").is_err());
        assert!(identifier("0").is_err());
        assert!(identifier("00145").is_err());
        assert!(identifier("+H").is_err());
        assert!(identifier("[H]").is_err());
        assert!(identifier("Øof").is_err());
        assert!(identifier("2xAc").is_err());
        assert!(identifier("-Ac").is_err());
        assert!(identifier("_Ac").is_err());
        // Multiple Identifiers
        assert_eq!(identifier("OH-p"), Ok(("-p", "OH")));
        assert_eq!(identifier("HeH 2slow"), Ok((" 2slow", "HeH")));
        assert_eq!(identifier("Gefählt"), Ok(("ählt", "Gef")));
        // This is a weird unicode 6
        assert!('𝟨'.is_numeric());
        assert!(!'𝟨'.is_ascii_digit());
        assert_eq!(identifier("C2H𝟨O"), Ok(("𝟨O", "C2H")));
    }

    #[test]
    fn test_multiplier() {
        macro_rules! assert_multiplier {
            ($input:literal, $output:literal, $count:literal) => {
                let (rest, count) = multiplier($input).unwrap();
                assert_eq!((rest, u32::from(count)), ($output, $count));
            };
        }

        // Valid Multipliers
        assert_multiplier!("1x", "", 1);
        assert_multiplier!("10x", "", 10);
        assert_multiplier!("422x", "", 422);
        assert_multiplier!("9999x", "", 9999);
        // Invalid Multipliers
        assert!(multiplier("1").is_err());
        assert!(multiplier("10").is_err());
        assert!(multiplier("422").is_err());
        assert!(multiplier("9999").is_err());
        assert!(multiplier("0").is_err());
        assert!(multiplier("01").is_err());
        assert!(multiplier("00145").is_err());
        assert!(multiplier("H").is_err());
        assert!(multiplier("p").is_err());
        assert!(multiplier("+H").is_err());
        assert!(multiplier("[H]").is_err());
        // Multiple Multipliers
        assert_multiplier!("1xOH", "OH", 1);
        assert_multiplier!("42xHeH", "HeH", 42);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_named_modification() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_named_modification = named_modification(&polymer);
        macro_rules! assert_named_modification {
            ($input:literal, $output:literal, $multiplier:literal, $name:expr, $mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_id) = named_modification(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let modification = polymer
                    .modification(parsed_id)
                    .unwrap()
                    .clone()
                    .unwrap_unlocalized();

                let multiplier = modification.multiplier();
                assert_eq!(u32::from(multiplier), $multiplier);

                let name = modification.kind().clone().unwrap_named().name();
                assert_eq!(name, $name);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
            };
        }

        // Valid Named Modifications
        assert_named_modification!("Am", "", 1, "Amidation", -0.98401558291);
        assert_named_modification!("Ac", "", 1, "O-Acetylation", 42.01056468403);
        assert_named_modification!("Poly", "", 1, "Wall Polymer Linkage", 77.95068082490);
        assert_named_modification!("DeAc", "", 1, "De-N-Acetylation", -42.01056468403);
        assert_named_modification!("Red", "", 1, "Reduced", 2.01565006446);
        assert_named_modification!("Anh", "", 1, "1,6-Anhydro", -18.01056468403);
        assert_named_modification!("1xAm", "", 1, "Amidation", -0.98401558291);
        assert_named_modification!("2xRed", "", 2, "Reduced", 4.03130012892);
        assert_named_modification!("3xAnh", "", 3, "1,6-Anhydro", -54.03169405209);
        // Invalid Named Modifications
        assert!(err_named_modification(" H2O").is_err());
        assert!(err_named_modification("1").is_err());
        assert!(err_named_modification("9999").is_err());
        assert!(err_named_modification("0").is_err());
        assert!(err_named_modification("00145").is_err());
        assert!(err_named_modification("+H").is_err());
        assert!(err_named_modification("[H]").is_err());
        assert!(err_named_modification("Øof").is_err());
        assert!(err_named_modification("-Ac").is_err());
        assert!(err_named_modification("_Ac").is_err());
        assert!(err_named_modification("+Am").is_err());
        assert!(err_named_modification("-2xAm").is_err());
        assert!(err_named_modification("(Am)").is_err());
        assert!(err_named_modification("-4xH2O").is_err());
        assert!(err_named_modification("-2p").is_err());
        assert!(err_named_modification("+C2H2O-2e").is_err());
        assert!(err_named_modification("-3xC2H2O-2e").is_err());
        assert!(err_named_modification("+NH3+p").is_err());
        assert!(err_named_modification("+2xD2O").is_err());
        assert!(err_named_modification("-2x[2H]2O").is_err());
        // Non-Existent Named Modifications
        assert!(err_named_modification("Blue").is_err());
        assert!(err_named_modification("Hydro").is_err());
        assert!(err_named_modification("1xAm2").is_err());
        assert!(err_named_modification("2xR_ed").is_err());
        // Multiple Named Modifications
        assert_named_modification!("Anh, Am", ", Am", 1, "1,6-Anhydro", -18.01056468403);
        assert_named_modification!("1xAm)JAA", ")JAA", 1, "Amidation", -0.98401558291);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_offset_modification() {
        use polychem::OffsetKind::{Add, Remove};
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_offset_modification = offset_modification(&polymer);
        macro_rules! assert_offset_modification {
            ($input:literal, $output:literal, $kind:ident, $multiplier:literal, $mass:literal, $charge:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_id) = offset_modification(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let modification = polymer
                    .modification(parsed_id)
                    .unwrap()
                    .clone()
                    .unwrap_unlocalized();

                let kind = modification.kind().clone().unwrap_offset().kind();
                assert_eq!(kind, $kind);

                let multiplier = modification.multiplier();
                assert_eq!(u32::from(multiplier), $multiplier);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
                assert_eq!(i64::from(polymer.charge()), $charge);
            };
        }

        // Valid Offset Modifications
        assert_offset_modification!("+H2O", "", Add, 1, 18.01056468403, 0);
        assert_offset_modification!("-H2O", "", Remove, 1, -18.01056468403, 0);
        assert_offset_modification!("+2xH2O", "", Add, 2, 36.02112936806, 0);
        assert_offset_modification!("-4xH2O", "", Remove, 4, -72.04225873612, 0);
        assert_offset_modification!("-2p", "", Remove, 1, -2.014552933242, -2);
        assert_offset_modification!("+H", "", Add, 1, 1.00782503223, 0);
        assert_offset_modification!("+C2H2O-2e", "", Add, 1, 42.009467524211870, 2);
        assert_offset_modification!("-3xC2H2O-2e", "", Remove, 3, -126.02840257263561, -6);
        assert_offset_modification!("+NH3+p", "", Add, 1, 18.033825567741, 1);
        assert_offset_modification!("+2xD2O", "", Add, 2, 40.04623635162, 0);
        assert_offset_modification!("-2x[2H]2O", "", Remove, 2, -40.04623635162, 0);
        assert_offset_modification!("+[37Cl]5-2p", "", Add, 1, 182.814960076758, -2);
        assert_offset_modification!("-NH2[99Tc]", "", Remove, 1, -114.92497486889, 0);
        // Invalid Offset Modifications
        assert!(err_offset_modification(" ").is_err());
        assert!(err_offset_modification("H2O").is_err());
        assert!(err_offset_modification("(-H2O)").is_err());
        assert!(err_offset_modification("+0xH2O").is_err());
        assert!(err_offset_modification("2xH2O").is_err());
        assert!(err_offset_modification("-2x3xH2O").is_err());
        assert!(err_offset_modification("-2x+H2O").is_err());
        assert!(err_offset_modification("+2[2H]").is_err());
        assert!(err_offset_modification("-[H+p]O").is_err());
        assert!(err_offset_modification("+NH2[100Tc]").is_err());
        // Multiple Offset Modifications
        assert_offset_modification!("+[37Cl]5-2p10", "10", Add, 1, 182.814960076758, -2);
        assert_offset_modification!("+[2H]2O*H2O", "*H2O", Add, 1, 20.02311817581, 0);
        assert_offset_modification!("+NH2{100Tc", "{100Tc", Add, 1, 16.01872406889, 0);
        assert_offset_modification!("+C11H12N2O2 H2O", " H2O", Add, 1, 204.08987763476, 0);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_modifications() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_modifications = modifications(&polymer);
        macro_rules! assert_modifications {
            ($input:literal, $output:literal, $count:literal, $mass:literal, $charge:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = modifications(&polymer)($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(parsed_ids.len(), $count);

                let polymer = polymer.borrow();
                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
                assert_eq!(i64::from(polymer.charge()), $charge);
            };
        }
        // Valid Modifications
        assert_modifications!("(-H2O)", "", 1, -18.01056468403, 0);
        assert_modifications!("(+2xH2O)", "", 1, 36.02112936806, 0);
        assert_modifications!("(-2p)", "", 1, -2.014552933242, -2);
        assert_modifications!("(-3xC2H2O-2e)", "", 1, -126.02840257263561, -6);
        assert_modifications!("(+[37Cl]5-2p)", "", 1, 182.814960076758, -2);
        assert_modifications!("(Red)", "", 1, 2.01565006446, 0);
        assert_modifications!("(Anh)", "", 1, -18.01056468403, 0);
        assert_modifications!("(1xAm)", "", 1, -0.98401558291, 0);
        assert_modifications!("(2xRed)", "", 1, 4.03130012892, 0);
        assert_modifications!("(-OH, +NH2)", "", 2, -0.98401558291, 0);
        assert_modifications!("(Anh, +H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh,+H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh   ,+H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh  ,  +H2O)", "", 2, 0, 0);
        assert_modifications!("(2xAnh, +3xH2O)", "", 2, 18.01056468403, 0);
        assert_modifications!("(Anh, Anh, +3xH2O)", "", 3, 18.01056468403, 0);
        assert_modifications!("(-H2, +Ca)", "", 2, 37.94694079854, 0);
        // NOTE: There is a super small mass defect (13.6 eV, or ~1e-8 u) stored in the binding energy between a proton
        // and electon — that's why this result is slightly different from the one above!
        assert_modifications!("(-2p, +Ca-2e)", "", 2, 37.946940769939870, 0);
        assert_modifications!("(+2p, -2p, +Ca-2e)", "", 3, 39.961493703181870, 2);
        // Invalid Modifications
        assert!(err_modifications(" ").is_err());
        assert!(err_modifications("H2O").is_err());
        assert!(err_modifications("(-H2O").is_err());
        assert!(err_modifications("(+0xH2O)").is_err());
        assert!(err_modifications("(2xH2O)").is_err());
        assert!(err_modifications("(-2x3xH2O)").is_err());
        assert!(err_modifications("(-2x+H2O)").is_err());
        assert!(err_modifications("(+2[2H])").is_err());
        assert!(err_modifications("(-[H+p]O)").is_err());
        assert!(err_modifications("(+NH2[100Tc])").is_err());
        assert!(err_modifications("( H2O)").is_err());
        assert!(err_modifications("(1)").is_err());
        assert!(err_modifications("(9999)").is_err());
        assert!(err_modifications("(0)").is_err());
        assert!(err_modifications("(00145)").is_err());
        assert!(err_modifications("([H])").is_err());
        assert!(err_modifications("(Øof)").is_err());
        assert!(err_modifications("(-Ac)").is_err());
        assert!(err_modifications("(_Ac)").is_err());
        assert!(err_modifications("(+Am)").is_err());
        assert!(err_modifications("(-2xAm)").is_err());
        assert!(err_modifications("((Am))").is_err());
        assert!(err_modifications("(Anh +H2O)").is_err());
        assert!(err_modifications("(Anh; +H2O)").is_err());
        // Non-Existent Modifications
        assert!(err_modifications("(Blue)").is_err());
        assert!(err_modifications("(Hydro)").is_err());
        assert!(err_modifications("(1xAm2)").is_err());
        assert!(err_modifications("(2xR_ed)").is_err());
        // Multiple Modifications
        assert_modifications!("(+[37Cl]5-2p)10", "10", 1, 182.814960076758, -2);
        assert_modifications!("(+[2H]2O)*H2O", "*H2O", 1, 20.02311817581, 0);
        assert_modifications!("(+NH2){100Tc", "{100Tc", 1, 16.01872406889, 0);
        assert_modifications!("(+C11H12N2O2) H2O", " H2O", 1, 204.08987763476, 0);
        assert_modifications!(
            "(2xAnh, +3xH2O)AA=gm-AEJA",
            "AA=gm-AEJA",
            2,
            18.01056468403,
            0
        );
    }

    // FIXME: Unfininshed! Needs modification support — same with unbranched_amino_acid!
    #[test]
    fn test_monosaccharide() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut monosaccharide = monosaccharide(&polymer);
        macro_rules! assert_monosaccharide_name {
            ($input:literal, $output:literal, $name:literal) => {
                let (rest, id) = monosaccharide($input).unwrap();
                assert_eq!(
                    (rest, polymer.borrow().residue(id).unwrap().name()),
                    ($output, $name)
                );
            };
        }

        // Valid Monosaccharides
        assert_monosaccharide_name!("g", "", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m", "", "N-Acetylmuramic Acid");
        // Invalid Monosaccharides
        assert!(monosaccharide("P").is_err());
        assert!(monosaccharide("EP").is_err());
        assert!(monosaccharide("1h").is_err());
        assert!(monosaccharide("+m").is_err());
        assert!(monosaccharide("-g").is_err());
        assert!(monosaccharide("[h]").is_err());
        // Non-Existent Monosaccharides
        assert!(monosaccharide("s").is_err());
        assert!(monosaccharide("f").is_err());
        // Multiple Monosaccharides
        assert_monosaccharide_name!("gm", "m", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m-A", "-A", "N-Acetylmuramic Acid");
    }

    // FIXME: Unfininshed! Needs modification support — same with monosaccharide!
    #[test]
    fn test_unbranched_amino_acid() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut unbranched_amino_acid = unbranched_amino_acid(&polymer);
        macro_rules! assert_unbranched_aa_name {
            ($input:literal, $output:literal, $name:literal) => {
                let (rest, id) = unbranched_amino_acid($input).unwrap();
                assert_eq!(
                    (rest, polymer.borrow().residue(id).unwrap().name()),
                    ($output, $name)
                );
            };
        }

        // Valid Unbranched Amino Acids
        assert_unbranched_aa_name!("A", "", "Alanine");
        assert_unbranched_aa_name!("E", "", "Glutamic Acid");
        assert_unbranched_aa_name!("J", "", "Diaminopimelic Acid");
        assert_unbranched_aa_name!("yE", "", "γ-Glutamate");
        assert_unbranched_aa_name!("eK", "", "ε-Lysine");
        // Invalid Unbranched Amino Acids
        assert!(unbranched_amino_acid("p").is_err());
        assert!(unbranched_amino_acid("eP").is_err());
        assert!(unbranched_amino_acid("1H").is_err());
        assert!(unbranched_amino_acid("+M").is_err());
        assert!(unbranched_amino_acid("-G").is_err());
        assert!(unbranched_amino_acid("[H]").is_err());
        // Non-Existent Unbranched Amino Acids
        assert!(unbranched_amino_acid("iA").is_err());
        assert!(unbranched_amino_acid("yK").is_err());
        // Multiple Unbranched Amino Acids
        assert_unbranched_aa_name!("AEJA", "EJA", "Alanine");
        assert_unbranched_aa_name!("EJA", "JA", "Glutamic Acid");
        assert_unbranched_aa_name!("JA", "A", "Diaminopimelic Acid");
        assert_unbranched_aa_name!("yEJA", "JA", "γ-Glutamate");
        assert_unbranched_aa_name!("eK[GGGGG]", "[GGGGG]", "ε-Lysine");
    }

    // FIXME: Add modification testing!
    // FIXME: Add lateral chain testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_peptide() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_peptide = peptide(&polymer);
        macro_rules! assert_chain_residues_and_masses {
            ($input:literal, $output:literal, $residues:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = peptide(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let parsed_ids: Vec<_> = parsed_ids
                    .into_iter()
                    .map(|id| polymer.residue(id.residue).unwrap().name())
                    .collect();
                let residues = Vec::from($residues);
                assert_eq!(parsed_ids, residues);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Peptides
        assert_chain_residues_and_masses!(
            "AEJA",
            "",
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_chain_residues_and_masses!(
            "AyEJA",
            "",
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_chain_residues_and_masses!(
            "AE",
            "",
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_chain_residues_and_masses!(
            "A",
            "",
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        // Invalid Peptides
        assert!(err_peptide("y").is_err());
        assert!(err_peptide("yrE").is_err());
        assert!(err_peptide("-AEJA").is_err());
        assert!(err_peptide("[GGGGG]").is_err());
        assert!(err_peptide("gm-AEJA").is_err());
        assert!(err_peptide("(Am)").is_err());
        // Non-Existent Peptide Residues
        assert!(err_peptide("AEJiA").is_err());
        assert!(err_peptide("AQyK").is_err());
        // Multiple Peptides
        assert_chain_residues_and_masses!(
            "AE=gm-AEJ",
            "=gm-AEJ",
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_chain_residues_and_masses!(
            "AeeK",
            "eeK",
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
    }

    // FIXME: Add modification testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_glycan() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_glycan = glycan(&polymer);
        macro_rules! assert_chain_residues_and_masses {
            ($input:literal, $output:literal, $residues:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = glycan(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let parsed_ids: Vec<_> = parsed_ids
                    .into_iter()
                    .map(|id| polymer.residue(id).unwrap().name())
                    .collect();
                let residues = Vec::from($residues);
                assert_eq!(parsed_ids, residues);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Glycans
        assert_chain_residues_and_masses!(
            "gmgm",
            "",
            [
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid",
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid"
            ],
            974.37031350523,
            974.91222678113779720
        );
        assert_chain_residues_and_masses!(
            "gm",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            496.19043909463,
            496.46375660678381490
        );
        assert_chain_residues_and_masses!(
            "g",
            "",
            ["N-Acetylglucosamine"],
            221.08993720530,
            221.20813124207411765
        );
        assert_chain_residues_and_masses!(
            "m",
            "",
            ["N-Acetylmuramic Acid"],
            293.11106657336,
            293.27091179713952985
        );
        // Invalid Glycans
        assert!(err_glycan("Y").is_err());
        assert!(err_glycan("Ygm").is_err());
        assert!(err_glycan("-AEJA").is_err());
        assert!(err_glycan("[GGGGG]").is_err());
        assert!(err_glycan("EA=gm-AEJA").is_err());
        assert!(err_glycan("(Am)").is_err());
        // Non-Existent Glycan Residues
        assert!(err_glycan("y").is_err());
        assert!(err_glycan("fp").is_err());
        // Multiple Glycans
        assert_chain_residues_and_masses!(
            "gm-AEJ",
            "-AEJ",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            496.19043909463,
            496.46375660678381490
        );
        assert_chain_residues_and_masses!("xAJgmK", "AJgmK", ["Unknown Monosaccharide"], 0.0, 0.0);
    }

    // FIXME: Add modification testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    #[allow(clippy::too_many_lines)]
    fn test_monomer() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_monomer = monomer(&polymer);
        macro_rules! assert_monomer_residues_and_masses {
            ($input:literal, $output:literal, $glycan:expr, $peptide:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, Monomer { glycan, peptide }) = monomer(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let glycan: Vec<_> = glycan
                    .into_iter()
                    .map(|id| polymer.residue(id).unwrap().name())
                    .collect();
                let peptide: Vec<_> = peptide
                    .into_iter()
                    .map(|id| polymer.residue(id.residue).unwrap().name())
                    .collect();
                assert_eq!(glycan, Vec::<&str>::from($glycan));
                assert_eq!(peptide, Vec::<&str>::from($peptide));

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Monomers
        assert_monomer_residues_and_masses!(
            "gmgm",
            "",
            [
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid",
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid"
            ],
            [],
            974.37031350523,
            974.91222678113779720
        );
        assert_monomer_residues_and_masses!(
            "gm",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            [],
            496.19043909463,
            496.46375660678381490
        );
        assert_monomer_residues_and_masses!(
            "g",
            "",
            ["N-Acetylglucosamine"],
            [],
            221.08993720530,
            221.20813124207411765
        );
        assert_monomer_residues_and_masses!(
            "m",
            "",
            ["N-Acetylmuramic Acid"],
            [],
            293.11106657336,
            293.27091179713952985
        );
        assert_monomer_residues_and_masses!(
            "AEJA",
            "",
            [],
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_monomer_residues_and_masses!(
            "AyEJA",
            "",
            [],
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_monomer_residues_and_masses!(
            "AE",
            "",
            [],
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_monomer_residues_and_masses!(
            "A",
            "",
            [],
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        assert_monomer_residues_and_masses!(
            "gm-AEJA",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            939.39205200801,
            939.91604006741105325
        );
        assert_monomer_residues_and_masses!(
            "gm-AyEJA",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            939.39205200801,
            939.91604006741105325
        );
        // Invalid Monomers
        assert!(err_monomer("-AEJA").is_err());
        assert!(err_monomer("[GGGGG]").is_err());
        assert!(err_monomer("(Am)").is_err());
        // Non-Existent Monomer Residues & Bonds
        assert!(err_monomer("y").is_err());
        assert!(err_monomer("fp").is_err());
        assert!(err_monomer("AEJiA").is_err());
        assert!(err_monomer("AQyK").is_err());
        assert!(err_monomer("g-A").is_err());
        // Multiple Monomers
        assert_monomer_residues_and_masses!(
            "gm,AEJ",
            ",AEJ",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            [],
            496.19043909463,
            496.46375660678381490
        );
        assert_monomer_residues_and_masses!(
            "xAJgmK",
            "AJgmK",
            ["Unknown Monosaccharide"],
            [],
            0.0,
            0.0
        );
        assert_monomer_residues_and_masses!(
            "AE=gm-AEJ",
            "=gm-AEJ",
            [],
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_monomer_residues_and_masses!(
            "AeeK",
            "eeK",
            [],
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        assert_monomer_residues_and_masses!(
            "gm-AE=gm-AEJA",
            "=gm-AEJA",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "Glutamic Acid"],
            696.27014596853,
            696.65595894949984270
        );
    }

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
