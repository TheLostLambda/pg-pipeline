use nom::{branch::alt, character::complete::char, sequence::terminated, Parser};
use nom_miette::{wrap_err, LabeledParseError};

use crate::{Count, ModificationId, Polymer};

use super::{
    count,
    errors::{ParseResult, PolychemErrorKind},
};

// FIXME: The errors for these parsers need to be tested and improved!
// FIXME: Ensure that users of these parsers *don't* need to use nom-miette!

/// Any Modification = Named Modification | Offset Modification
pub fn any<'a, 'p, 's, K>(
    polymer: &mut Polymer<'a, 'p>,
    identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind> + From<nom::error::ErrorKind>,
{
    alt((named(polymer, identifier), offset::<K>(polymer)))
}

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Named Modification = [ Multiplier ] , Identifier
pub fn named<'a, 'p, 's, K>(
    _polymer: &mut Polymer<'a, 'p>,
    _identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind> + From<nom::error::ErrorKind>,
{
    |_| todo!()
}

/// Offset Modification = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
pub fn offset<'a, 's, K>(
    _polymer: &mut Polymer<'a, '_>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind>,
{
    |_| todo!()
}

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let parser = terminated(count, char('x'));
    wrap_err(parser, PolychemErrorKind::ExpectedMultiplier)(i)
}

#[cfg(test)]
mod tests {
    use nom::{
        bytes::complete::tag,
        character::complete::{alpha1, alphanumeric1},
        combinator::recognize,
        multi::many0,
        sequence::pair,
    };
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use super::*;
    use crate::{AtomicDatabase, PolymerDatabase};

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    /// Identifier = letter , { letter | digit | "_" } ;
    fn identifier(i: &str) -> ParseResult<&str> {
        recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))))(i)
    }

    #[test]
    fn test_multiplier() {
        macro_rules! assert_multiplier {
            ($input:literal, $output:literal, $count:expr) => {
                assert_eq!(
                    multiplier($input),
                    Ok(($output, Count::new($count).unwrap()))
                );
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
        let mut named_modification = named(todo!(), identifier);
        macro_rules! assert_offset_mass {
            ($input:literal, $output:literal, $mass:expr) => {
                let (rest, modification) = named_modification($input).unwrap();
                assert_eq!(rest, $output);
                // TODO: Need to get the monoisotopic mass of the modification, extracted from the `Polymer`...
                assert_eq!(dec!(0), $mass);
                // FIXME: This also needs to test charge!
            };
        }
        // Valid Named Modifications
        assert_offset_mass!("Am", "", dec!(-0.98401558291));
        assert_offset_mass!("Ac", "", dec!(42.01056468403));
        assert_offset_mass!("Poly", "", dec!(77.95068082490));
        assert_offset_mass!("DeAc", "", dec!(-42.01056468403));
        assert_offset_mass!("Red", "", dec!(2.01565006446));
        assert_offset_mass!("Anh", "", dec!(-18.01056468403));
        assert_offset_mass!("1xAm", "", dec!(-0.98401558291));
        assert_offset_mass!("2xRed", "", dec!(4.03130012892));
        assert_offset_mass!("3xAnh", "", dec!(-54.03169405209));
        // Invalid Named Modifications
        assert!(named_modification(" H2O").is_err());
        assert!(named_modification("1").is_err());
        assert!(named_modification("9999").is_err());
        assert!(named_modification("0").is_err());
        assert!(named_modification("00145").is_err());
        assert!(named_modification("+H").is_err());
        assert!(named_modification("[H]").is_err());
        assert!(named_modification("Øof").is_err());
        assert!(named_modification("-Ac").is_err());
        assert!(named_modification("_Ac").is_err());
        assert!(named_modification("+Am").is_err());
        assert!(named_modification("-2xAm").is_err());
        assert!(named_modification("(Am)").is_err());
        assert!(named_modification("-4xH2O").is_err());
        assert!(named_modification("-2p").is_err());
        assert!(named_modification("+C2H2O-2e").is_err());
        assert!(named_modification("-3xC2H2O-2e").is_err());
        assert!(named_modification("+NH3+p").is_err());
        assert!(named_modification("+2xD2O").is_err());
        assert!(named_modification("-2x[2H]2O").is_err());
        // Non-Existent Named Modifications
        assert!(named_modification("Blue").is_err());
        assert!(named_modification("Hydro").is_err());
        assert!(named_modification("1xAm2").is_err());
        assert!(named_modification("2xR_ed").is_err());
        // Multiple Named Modifications
        assert_offset_mass!("Anh, Am", ", Am", dec!(-18.01056468403));
        assert_offset_mass!("1xAm)JAA", ")JAA", dec!(-0.98401558291));
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_offset_modification() {
        let mut offset_modification = offset::<PolychemErrorKind>(todo!());
        macro_rules! assert_offset_mz {
            ($input:literal, $output:literal, $mass:expr, $charge:literal) => {
                let (rest, modification) = offset_modification($input).unwrap();
                assert_eq!(rest, $output);
                // TODO: Need to get modification from `Polymer`, then check mass and charge!
                assert_eq!(dec!(0), $mass);
                assert_eq!(0, $charge);
            };
        }
        // Valid Offset Modifications
        assert_offset_mz!("+H2O", "", dec!(18.01056468403), 0);
        assert_offset_mz!("-H2O", "", dec!(-18.01056468403), 0);
        assert_offset_mz!("+2xH2O", "", dec!(36.02112936806), 0);
        assert_offset_mz!("-4xH2O", "", dec!(-72.04225873612), 0);
        assert_offset_mz!("-2p", "", dec!(-2.014552933242), -2);
        assert_offset_mz!("+H", "", dec!(1.00782503223), 0);
        assert_offset_mz!("+C2H2O-2e", "", dec!(42.009467524211870), 2);
        assert_offset_mz!("-3xC2H2O-2e", "", dec!(-126.02840257263561), -6);
        assert_offset_mz!("+NH3+p", "", dec!(18.033825567741), 1);
        assert_offset_mz!("+2xD2O", "", dec!(40.04623635162), 0);
        assert_offset_mz!("-2x[2H]2O", "", dec!(-40.04623635162), 0);
        assert_offset_mz!("+[37Cl]5-2p", "", dec!(182.814960076758), -2);
        assert_offset_mz!("-NH2[99Tc]", "", dec!(-114.92497486889), 0);
        // Invalid Offset Modifications
        assert!(offset_modification(" ").is_err());
        assert!(offset_modification("H2O").is_err());
        assert!(offset_modification("(-H2O)").is_err());
        assert!(offset_modification("+0xH2O").is_err());
        assert!(offset_modification("2xH2O").is_err());
        assert!(offset_modification("-2x3xH2O").is_err());
        assert!(offset_modification("-2x+H2O").is_err());
        assert!(offset_modification("+2[2H]").is_err());
        assert!(offset_modification("-[H+p]O").is_err());
        assert!(offset_modification("+NH2[100Tc]").is_err());
        // Multiple Offset Modifications
        assert_offset_mz!("+[37Cl]5-2p10", "10", dec!(182.814960076758), -2);
        assert_offset_mz!("+[2H]2O*H2O", "*H2O", dec!(20.02311817581), 0);
        assert_offset_mz!("+NH2{100Tc", "{100Tc", dec!(16.01872406889), 0);
        assert_offset_mz!("+C11H12N2O2 H2O", " H2O", dec!(204.08987763476), 0);
    }
}
