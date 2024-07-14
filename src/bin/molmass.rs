use miette::{Diagnostic, GraphicalReportHandler, GraphicalTheme};
use once_cell::sync::Lazy;
use polychem::{AtomicDatabase, Charged, ChargedParticle, ChemicalComposition, Massive, Result};
use rustyline::DefaultEditor;
use std::fmt::Write;

static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

fn main() {
    let mut rl = DefaultEditor::new().unwrap();
    while let Ok(formula) = rl.readline("Molecule: ") {
        rl.add_history_entry(&formula).unwrap();
        match molecule_info(&formula) {
            Ok(info) => print!("{info}"),
            Err(diagnostic) => render_error(*diagnostic),
        }
    }
}

fn molecule_info(formula: &str) -> Result<String> {
    let mut buf = String::new();
    let molecule = ChemicalComposition::new(&DB, formula)?;

    let mono_mass = molecule.monoisotopic_mass();
    let avg_mass = molecule.average_mass();
    let charge = molecule.charge();

    writeln!(buf, "Monoisotopic Mass: {mono_mass:.6}").unwrap();
    writeln!(buf, "Average Mass: {avg_mass:.4}").unwrap();
    writeln!(buf, "Charge: {charge}").unwrap();

    if i64::from(charge) != 0 {
        let mono_mz = molecule.monoisotopic_mz().unwrap();
        let avg_mz = molecule.average_mz().unwrap();
        writeln!(buf, "Monoisotopic m/z: {mono_mz:.6}").unwrap();
        writeln!(buf, "Average m/z: {avg_mz:.4}").unwrap();
    }

    writeln!(buf).unwrap();

    Ok(buf)
}

fn render_error(diagnostic: impl Into<Box<dyn Diagnostic + 'static>>) {
    let mut buf = String::new();
    GraphicalReportHandler::new_themed(GraphicalTheme::unicode())
        .render_report(&mut buf, diagnostic.into().as_ref())
        .unwrap();
    println!("{buf}");
}
