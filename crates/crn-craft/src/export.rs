//! Explicit elementary stochastic-count exporter. No target runtime is linked.
//!
//! [`export_gillespy2`] translates a validated model into Python source for
//! GillesPy2 [`GILLESPY2_VERSION`] as a `crn-craft.export/0.1` report. Only
//! mass action with total reactant order 0, 1 or 2 and no population limits
//! is supported; anything else fails with `CRN-E050` rather than being
//! dropped, and a homodimer coefficient that cannot be halved exactly fails
//! with `CRN-E051`. The source is for reference cross-checks; it carries no
//! budgets, seeds or time grid. See
//! [docs/reference-export.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference-export.md).
use crate::{
    contracts::{Rate, Scalar},
    diagnostic::{Diagnostic, Severity},
    parse_model, sha256,
    validate::validate_model,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The GillesPy2 version the generated source requires; it raises `RuntimeError` on any other.
pub const GILLESPY2_VERSION: &str = "1.8.3";

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of an export report.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ExportVersion {
    #[doc = concat!("`crn-craft.export/0.1`.")]
    #[serde(rename = "crn-craft.export/0.1")]
    V01,
}
/// Generated reference source for one target and how its names map back to the model.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PythonExport {
    /// Target package, `gillespy2`.
    pub target: String,
    /// Exact target version the source checks for at run time.
    pub target_version: String,
    /// Python module defining `build_model()`, with one frozen rate constant per channel and
    /// the target volume fixed at 1.
    pub source: String,
    /// Model species name to target name (`S0`, `S1`, … in declaration order).
    pub species_aliases: BTreeMap<String, String>,
    /// Model reaction name to target name (`R0`, `R1`, … in declaration order).
    pub reaction_aliases: BTreeMap<String, String>,
    /// Semantics the translation does not preserve, for readers.
    pub notes: Vec<String>,
}
/// Result of `export`: success flag, model hash, the generated source and diagnostics.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportReport {
    #[doc = concat!("Always `crn-craft.export/0.1`.")]
    pub schema: ExportVersion,
    /// True when source was generated; then no error diagnostic is present.
    pub success: bool,
    /// Always false: export runs no simulation.
    pub simulation_available: bool,
    /// SHA-256 of the exact model bytes; null only when the model file could not be read.
    pub model_sha256: Option<String>,
    /// The generated source on success, otherwise null; no partial source is emitted.
    pub artifact: Option<PythonExport>,
    /// Validation findings, plus `CRN-E050` and `CRN-E051` for unsupported features.
    pub diagnostics: Vec<Diagnostic>,
}
impl ExportReport {
    /// A failed report carrying one diagnostic, such as a read error (`CRN-E030`).
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: ExportVersion::V01,
            success: false,
            simulation_available: false,
            model_sha256: None,
            artifact: None,
            diagnostics: vec![error],
        }
    }
}
fn unsupported(path: &str, message: &str) -> Diagnostic {
    Diagnostic::error(
        "CRN-E050",
        path,
        message,
        "This exporter accepts elementary mass action without population limits; keep unsupported semantics in the source model.",
    )
}
/// Exports a model from its exact bytes; see the [module docs](self) for
/// the supported subset. Failures give `success: false` and no artifact.
///
/// ```
/// let model = br#"{"schema": "crn-craft.model/0.1", "name": "decay",
///   "convention": "stochastic_counts", "omega": 1.0, "parameters": {"k": 0.1},
///   "species": [{"name": "A", "initial_count": 10}],
///   "reactions": [{"name": "death", "reactants": {"A": 1}, "products": {},
///                  "rate": {"kind": "mass_action", "coefficient": "k"}}]}"#;
/// let report = crn_craft::export::export_gillespy2(model);
/// assert!(report.success);
/// let artifact = report.artifact.unwrap();
/// assert!(artifact.source.contains("def build_model():"));
/// assert_eq!(artifact.species_aliases["A"], "S0");
/// ```
pub fn export_gillespy2(model_bytes: &[u8]) -> ExportReport {
    let mut report = ExportReport {
        schema: ExportVersion::V01,
        success: false,
        simulation_available: false,
        model_sha256: Some(sha256(model_bytes)),
        artifact: None,
        diagnostics: Vec::new(),
    };
    let model = match parse_model(model_bytes) {
        Ok(model) => model,
        Err(error) => {
            report.diagnostics.push(error);
            return report;
        }
    };
    report.diagnostics = validate_model(&model);
    if report
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error)
    {
        return report;
    }
    for (i, species) in model.species.iter().enumerate() {
        if species.population_limit.is_some() {
            report.diagnostics.push(unsupported(
                &format!("/species/{i}/population_limit"),
                "GillesPy2 cannot preserve the CRN-craft population safety stop",
            ));
        }
        if species.initial_count > 9_007_199_254_740_992 {
            report.diagnostics.push(unsupported(
                &format!("/species/{i}/initial_count"),
                "initial counts above 2^53 cannot be represented exactly by the target solver",
            ));
        }
    }
    let mut coefficients = Vec::new();
    for (i, reaction) in model.reactions.iter().enumerate() {
        let Rate::MassAction { coefficient } = &reaction.rate else {
            report.diagnostics.push(unsupported(
                &format!("/reactions/{i}/rate"),
                "custom propensity export is not implemented",
            ));
            continue;
        };
        if reaction.reactants.values().sum::<i64>() > 2 {
            report.diagnostics.push(unsupported(
                &format!("/reactions/{i}/reactants"),
                "GillesPy2 elementary mass action supports total reactant order at most two",
            ));
            continue;
        }
        let coefficient = match coefficient {
            Scalar::Literal(x) => *x,
            Scalar::Parameter(name) => model.parameters[name],
        };
        // GillesPy2's homodimer uses k*x*(x-1)/volume, without the /2 factorial.
        let translated = if reaction.reactants.values().any(|&v| v == 2) {
            coefficient / 2.0
        } else {
            coefficient
        };
        if reaction.reactants.values().any(|&v| v == 2) && translated * 2.0 != coefficient {
            report.diagnostics.push(Diagnostic::error("CRN-E051", &format!("/reactions/{i}/rate"),
                "homodimer coefficient cannot be halved exactly in f64",
                "Rescale the reference fixture; do not silently round or disable the exported channel."));
        }
        coefficients.push(translated);
    }
    if report
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error)
    {
        return report;
    }
    let aliases: BTreeMap<_, _> = model
        .species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.clone(), format!("S{i}")))
        .collect();
    let indices: BTreeMap<_, _> = model
        .species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();
    let mut source = format!(
        "# Generated by CRN-craft. Stochastic-count semantics; reference SSA only.\n# Model SHA-256: {}\nimport gillespy2 as gp\n\ndef build_model():\n    if gp.__version__ != {GILLESPY2_VERSION:?}:\n        raise RuntimeError('Requires GillesPy2 {GILLESPY2_VERSION}')\n    model = gp.Model(name={:?}, population=True, volume=1.0)\n    species = [\n",
        sha256(model_bytes),
        model.name
    );
    for (i, s) in model.species.iter().enumerate() {
        source.push_str(&format!(
            "        gp.Species(name='S{i}', initial_value={}, mode='discrete'),\n",
            s.initial_count
        ));
    }
    source.push_str("    ]\n    model.add_species(species)\n");
    let side = |terms: &BTreeMap<String, i64>| {
        let text: Vec<_> = terms
            .iter()
            .map(|(name, c)| format!("species[{}]: {c}", indices[name.as_str()]))
            .collect();
        format!("{{{}}}", text.join(", "))
    };
    for (i, reaction) in model.reactions.iter().enumerate() {
        source.push_str(&format!("    rate = gp.Parameter(name='K{i}', expression={:?})\n    model.add_parameter(rate)\n    model.add_reaction(gp.Reaction(name='R{i}', reactants={}, products={}, rate=rate))\n",
            coefficients[i], side(&reaction.reactants), side(&reaction.products)));
    }
    source.push_str("    return model\n");
    report.artifact = Some(PythonExport { target: "gillespy2".into(), target_version: GILLESPY2_VERSION.into(), source,
        species_aliases: aliases, reaction_aliases: model.reactions.iter().enumerate().map(|(i, r)| (r.name.clone(), format!("R{i}"))).collect(),
        notes: vec![
            "Target volume is fixed at 1: source coefficients already use stochastic counts; Omega does not rescale them.".into(),
            "Homodimer rate constants are divided by 2 to preserve the binomial factorial. Parameters are resolved into one constant per channel.".into(),
            "Use a stochastic solver. Target ODE semantics, live budgets, i64 range/underflow diagnostics and cross-runtime bit identity are not preserved.".into(),
            "The target uses floating-point populations; even supported initial counts may exceed its exact range during a reference run.".into(),
        ] });
    report.success = true;
    report
}
