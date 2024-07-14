use muropeptide::{Muropeptide, POLYMERIZER};
use polychem::{ChargedParticle, Massive};
use smithereens::Dissociable;
use std::fmt::Write;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
struct Peptidoglycan(Muropeptide<'static, 'static>);

#[wasm_bindgen]
impl Peptidoglycan {
    #[wasm_bindgen(constructor)]
    pub fn new(structure: &str) -> Result<Peptidoglycan, String> {
        // NOTE: This ensures the panic hook is set before any other shim code can be run!
        console_error_panic_hook::set_once();
        Muropeptide::new(&POLYMERIZER, structure)
            .map(Self)
            .map_err(|e| e.to_string())
    }

    pub fn monoisotopic_mass(&self) -> String {
        let mass = self.0.monoisotopic_mass();
        format!("{mass:.6}")
    }

    pub fn fragment(&self) -> String {
        let mut fragments: Vec<_> = self
            .0
            .fragment(None)
            .map(|fragment| (fragment.to_string(), fragment.monoisotopic_mz().unwrap()))
            .collect();
        fragments.sort_unstable_by(|(s1, mz1), (s2, mz2)| {
            mz1.cmp(mz2).then_with(|| s1.cmp(s2)).reverse()
        });
        fragments.dedup();

        let mut csv = String::new();
        for (structure, mz) in fragments {
            writeln!(&mut csv, r#""{structure}",{mz:.6}"#).unwrap();
        }
        csv
    }
}