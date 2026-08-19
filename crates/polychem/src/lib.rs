//! An abstraction for building chemically validated polymers
pub mod atoms;
pub mod errors;
pub mod moieties;
pub mod parsers;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use std::num::NonZero;

use derive_more::{
    Add, AddAssign, Display, From, Into, IsVariant, Neg, Sub, SubAssign, Sum, Unwrap,
};
use polymers::polymerizer_state::PolymerizerState;
use serde::Serialize;

// External Crate Imports
use ahash::{HashMap, HashSet};
use rust_decimal::Decimal;

// FIXME: Work on what's publicly exported / part of the API! — maybe create a prelude?
pub use atoms::atomic_database::AtomicDatabase;
pub use errors::Result;
pub use moieties::polymer_database::PolymerDatabase;
pub use polymers::polymerizer::Polymerizer;

// FIXME: I've exported a lot of things that previously weren't exported! Make sure that all of that new public API has
// as many traits automatically derived as possible!

// Core Data Types =====================================================================================================

// NOTE: For the types in this module, 'a lifetimes indicate references to the AtomicDatabase, whilst 'p lifetimes
// indicate references to the PolymerDatabase
#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct Polymer<'a, 'p> {
    // DESIGN: The `PolymerizerState` struct keeps track of fields that are only useful when adding new components to
    // the current polymer. Abstracting those fields out into a `PolymerizerState` helps to separate concerns, and makes
    // it clearer that the `residues`, `modifications`, and `bonds` fields are, by themselves, a complete description of
    // the `Polymer`s structure. A discarded alternative was two versions of the `Polymer` struct ("complete" and
    // "non-complete" versions), which would have introduced additional complexity to the user API, and the free group
    // index lost by discarding the `PolymerizerState` would have needed to be reconstructed during fragmentation
    // anyways for quick group lookup. Having two structs just for the ID counter and database references didn't seem
    // worth it. This decision accepts the downside of increasing the size of the `Polymer` struct.
    #[serde(skip)]
    polymerizer_state: PolymerizerState<'a, 'p>,
    // NOTE: Whilst the `Polymerizer` struct is defined elsewhere, the following fields are part of the core polymer
    // representation and are therefore defined in this file, with `impl`s added in separate modules as needed.
    residues: HashMap<ResidueId, Residue<'a, 'p>>,
    modifications: HashMap<ModificationId, ModificationInfo<'a, 'p>>,
    bonds: HashMap<BondId, BondInfo<'a, 'p>>,
}

// ---------------------------------------------------------------------------------------------------------------------

// FIXME: Pass through Display implementations for all Id newtypes!
// FIXME: Add tests that fail if `Default` is implemented for `Id`s!
// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Serialize)]
pub struct ResidueId(Id);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct Residue<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    composition: &'p ChemicalComposition<'a>,
    functional_groups: HashMap<FunctionalGroup<'p>, GroupState>,
    offset_modifications: HashSet<ModificationId>,
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Serialize)]
pub struct ModificationId(Id);

#[derive(Clone, Eq, PartialEq, Debug, IsVariant, Unwrap, Serialize)]
pub enum ModificationInfo<'a, 'p> {
    Named(NamedMod<'a, 'p>, ResidueGroup<'p>),
    Offset(Modification<OffsetMod<'a>>, ResidueId),
    Unlocalized(AnyModification<'a, 'p>),
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Serialize)]
pub struct BondId(Id);

// FIXME: Perhaps I should consider changing these `*Info` structs to have named fields? Is the donor -> acceptor order
// obvious enough for users?
#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct BondInfo<'a, 'p>(pub ResidueGroup<'p>, pub Bond<'a, 'p>, pub ResidueGroup<'p>);

// ---------------------------------------------------------------------------------------------------------------------

// NOTE: This underlying `Id` type is just a synonym since it's private to this crate — there is no way for users to
// directly provide or modify `Id`s, so I don't need to worry about a newtype wrapper here
type Id = u64;

#[derive(Clone, Eq, PartialEq, Debug, Default, Serialize)]
pub struct ChemicalComposition<'a> {
    chemical_formula: Vec<(Element<'a>, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle<'a>)>,
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct FunctionalGroup<'p> {
    name: &'p str,
    location: &'p str,
}

#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, IsVariant, Serialize,
)]
pub enum GroupState {
    #[default]
    Free,
    Modified(ModificationId),
    Donor(BondId),
    Acceptor(BondId),
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct NamedMod<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

// FIXME: Should I be using named fields here?
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug, Serialize)]
pub struct ResidueGroup<'p>(pub ResidueId, pub FunctionalGroup<'p>);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct OffsetMod<'a> {
    kind: OffsetKind,
    composition: ChemicalComposition<'a>,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct Modification<K> {
    multiplier: Count,
    kind: K,
}

// FIXME: Remove this? Or add a comment saying why it should exist?
pub type AnyModification<'a, 'p> = Modification<AnyMod<'a, 'p>>;

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct Bond<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Element<'a> {
    symbol: &'a str,
    name: &'a str,
    mass_number: Option<MassNumber>,
    isotopes: &'a HashMap<MassNumber, Isotope>,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Into, Serialize)]
pub struct Count(NonZero<u32>);

// MISSING: No `Default` — this *should* be user constructable, but there is no sensible default here
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub enum OffsetKind {
    Add,
    Remove,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Particle<'a> {
    symbol: &'a str,
    name: &'a str,
    mass: &'a Mass,
    charge: &'a Charge,
}

#[derive(Clone, Eq, PartialEq, Debug, IsVariant, Unwrap, Serialize)]
pub enum AnyMod<'a, 'p> {
    Named(NamedMod<'a, 'p>),
    Offset(OffsetMod<'a>),
}

// ---------------------------------------------------------------------------------------------------------------------

// MISSING: No `Default` — this *should* be user constructable, but there is no sensible default here
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, From, Into, Serialize,
)]
pub struct MassNumber(NonZero<u32>);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Isotope {
    relative_mass: Mass,
    abundance: Option<Abundance>,
}

// NOTE: `Mass` is *private* to this crate, and must be converted to either a `MonoisotopicMass` or `AverageMass`
// before being passed to the user. `Mass` should *not* show up in public API!
#[derive(Copy, Clone, Eq, PartialEq, Debug, Neg, Add, Sum, Serialize)]
struct Mass(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Display,
    Into,
    Neg,
    Add,
    Sub,
    Sum,
    AddAssign,
    SubAssign,
    Serialize,
)]
pub struct Charge(i64);

// ---------------------------------------------------------------------------------------------------------------------

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Into, Serialize)]
struct Abundance(Decimal);

// =====================================================================================================================

pub trait Massive {
    fn monoisotopic_mass(&self) -> MonoisotopicMass;
    fn average_mass(&self) -> AverageMass;
}

pub trait Charged {
    fn charge(&self) -> Charge;
}

// FIXME: Not super sold on that trait name...
pub trait ChargedParticle: Massive + Charged {
    fn monoisotopic_mz(&self) -> Option<MonoisotopicMz> {
        let mass = self.monoisotopic_mass();
        let charge = self.charge().abs();
        (charge.0 != 0).then(|| mass / charge)
    }

    fn average_mz(&self) -> Option<AverageMz> {
        let mass = self.average_mass();
        let charge = self.charge().abs();
        (charge.0 != 0).then(|| mass / charge)
    }
}

// ---------------------------------------------------------------------------------------------------------------------

// DESIGN: These types are all "duplicated" instead of providing something like `Monoisotopic<T>` and `Average<T>`,
// since writing something like `Monoisotopic<Mass>` would leak the private `Mass` struct into the API, and would
// prevent me from hiding that implementation detail from the user. We're explicitly accepting the repetition here...
// MISSING: No `Default` — should not be constructable by the user
#[derive(
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Display,
    Into,
    Neg,
    Add,
    Sub,
    Sum,
    AddAssign,
    SubAssign,
    Serialize,
)]
pub struct MonoisotopicMass(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Display,
    Into,
    Neg,
    Add,
    Sub,
    Sum,
    AddAssign,
    SubAssign,
    Serialize,
)]
pub struct AverageMass(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Display,
    Into,
    Add,
    Sub,
    Sum,
    AddAssign,
    SubAssign,
    Serialize,
)]
pub struct MonoisotopicMz(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Display,
    Into,
    Add,
    Sub,
    Sum,
    AddAssign,
    SubAssign,
    Serialize,
)]
pub struct AverageMz(Decimal);

// =====================================================================================================================

// Blanket impls

macro_rules! massive_ref_impls {
    ($($ref_type:ty),+ $(,)?) => {
        $(
            impl<T: Massive> Massive for $ref_type {
                fn monoisotopic_mass(&self) -> MonoisotopicMass {
                    (**self).monoisotopic_mass()
                }

                fn average_mass(&self) -> AverageMass {
                    (**self).average_mass()
                }
            }
        )+
    };
}

massive_ref_impls!(&T, &mut T, Box<T>);

macro_rules! charged_ref_impls {
    ($($ref_type:ty),+ $(,)?) => {
        $(
            impl<T: Charged> Charged for $ref_type {
                fn charge(&self) -> Charge {
                    (**self).charge()
                }
            }
        )+
    };
}

charged_ref_impls!(&T, &mut T, Box<T>);

impl<T: Massive + Charged> ChargedParticle for T {}

#[cfg(test)]
mod tests {
    use rust_decimal_macros::dec;

    use super::*;

    #[derive(Clone)]
    struct Chonky;

    impl Massive for Chonky {
        fn monoisotopic_mass(&self) -> MonoisotopicMass {
            MonoisotopicMass(dec!(42))
        }

        fn average_mass(&self) -> AverageMass {
            AverageMass(dec!(42.42))
        }
    }

    #[test]
    fn massive_blanket_impls() {
        let big_boi = Chonky;
        let chonky_ref = &big_boi;
        let chonky_mut = &mut big_boi.clone();
        let chonky_box = Box::new(big_boi.clone());
        assert_eq!(big_boi.monoisotopic_mass(), chonky_ref.monoisotopic_mass());
        assert_eq!(
            chonky_ref.monoisotopic_mass(),
            chonky_mut.monoisotopic_mass(),
        );
        assert_eq!(
            chonky_mut.monoisotopic_mass(),
            chonky_box.monoisotopic_mass(),
        );
        assert_eq!(big_boi.average_mass(), chonky_ref.average_mass());
        assert_eq!(chonky_ref.average_mass(), chonky_mut.average_mass());
        assert_eq!(chonky_mut.average_mass(), chonky_box.average_mass());
    }

    #[derive(Clone)]
    struct Pikachu;

    impl Charged for Pikachu {
        fn charge(&self) -> Charge {
            Charge(9001)
        }
    }

    #[test]
    fn charged_blanket_impls() {
        let yellow_rat = Pikachu;
        let pikachu_ref = &yellow_rat;
        let pikachu_mut = &mut yellow_rat.clone();
        let pikachu_box = Box::new(yellow_rat.clone());
        assert_eq!(yellow_rat.charge(), pikachu_ref.charge());
        assert_eq!(pikachu_ref.charge(), pikachu_mut.charge());
        assert_eq!(pikachu_mut.charge(), pikachu_box.charge());
    }
}
