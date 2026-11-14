//! Stochastic-count exporter to GillesPy2 source. No target runtime is linked.
//!
//! [`export_gillespy2`] translates a model of any version the engine reads
//! (a network `crn-craft.model/0.1`, converted first, a rules model or a
//! `crn-craft.model/1.0`) into Python source for GillesPy2
//! [`GILLESPY2_VERSION`] as a `crn-craft.export/0.1` report. It works on the
//! model's flat network ([`FlatView::with_every_type`]): one GillesPy2
//! species per species label and one reaction per channel.
//!
//! - A channel whose rule has a constant `rate` (one that reads no
//!   observable) and that takes at most two molecules exports as GillesPy2
//!   mass action, with the rate times the patterns' `factor` as its
//!   constant; a homodimer's constant is halved, and a constant that cannot
//!   be halved exactly in f64 fails with `CRN-E051`.
//! - Any other channel exports as a custom propensity string: the rule's
//!   expression with parameter values inlined, each observable as the
//!   weighted sum of its species, `^` as `**`, times the falling factorial
//!   of the reactant counts (`rate` rules) or times an availability factor
//!   `x / (x + 1e-300)` per molecule taken (whole-rule `propensity`).
//!
//! What GillesPy2 cannot evaluate in both its Python and C++ solvers fails
//! with `CRN-E050` rather than being dropped: comparisons, logic and
//! functions in an expression, observables reduced other than by `sum`,
//! more than 64 molecules of one species in a custom propensity, a
//! negative constant rate, and a `stop` condition. The source is for
//! reference cross-checks; it carries no budgets, seeds or time grid. See
//! [docs/reference/gillespy2-export.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/gillespy2-export.md).
use crate::{
    diagnostic::{Diagnostic, Severity, pointer},
    rules::FlatView,
    sha256,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The GillesPy2 version the generated source requires; it raises `RuntimeError` on any other.
pub const GILLESPY2_VERSION: &str = "1.8.3";

/// Most molecules of one species a custom propensity writes out as factors.
pub const MAX_CUSTOM_ORDER: u64 = 64;

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
    /// Python module defining `build_model()`, with one frozen rate constant per mass-action
    /// reaction, a custom propensity string for every other reaction, and the target volume
    /// fixed at 1.
    pub source: String,
    /// Species label to target name (`S0`, `S1`, … in the order of the flat view; for a
    /// network model its species in declaration order).
    pub species_aliases: BTreeMap<String, String>,
    /// Rule name to target name (`R0`, `R1`, … one per channel of the flat view; for a
    /// network model one per reaction). The channels of a rule with several reactant tuples
    /// are named `rule[0]`, `rule[1]`, ….
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
    /// The errors of reading and lowering the model, plus `CRN-E050` and `CRN-E051` for
    /// unsupported features.
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
fn unsupported(path: &str, message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E050",
        path,
        message,
        "This exporter accepts what GillesPy2 can evaluate: mass action, and custom propensities of + - * / and ^ over parameters and summed observables; keep unsupported semantics in the source model.",
    )
}

/// How one channel is written.
enum Law {
    /// GillesPy2 mass action with this rate constant.
    MassAction(f64),
    /// A custom propensity string.
    Custom(String),
}

/// A number as a Python literal that reads back to the same `f64`.
fn literal(x: f64) -> String {
    if x < 0.0 {
        format!("(-{:?})", -x)
    } else {
        format!("{x:?}")
    }
}

/// A rule's rate or propensity as a Python expression over the target
/// species, or the diagnostic of its first unsupported construct.
fn expression(view: &FlatView, rule: usize, path: &str) -> Result<String, Diagnostic> {
    view.python_rate(rule).map_err(|construct| {
        unsupported(
            path,
            format!("a custom propensity cannot contain {construct}"),
        )
    })
}

/// The custom propensity of a channel.
fn custom(view: &FlatView, channel: usize) -> Result<String, Diagnostic> {
    let rule = view.rule(channel);
    let whole = view.whole(channel);
    let base = format!("/rules/{rule}");
    let path = pointer(&base, if whole { "propensity" } else { "rate" });
    // Operations come parenthesized, so the expression is one factor.
    let mut factors = vec![expression(view, rule, &path)?];
    if !whole && view.factor(channel) != 1.0 {
        factors.push(literal(view.factor(channel)));
    }
    for &(species, count) in view.reactants(channel) {
        if count > MAX_CUSTOM_ORDER {
            return Err(unsupported(
                &pointer(&base, "reactants"),
                format!(
                    "a custom propensity writes at most {MAX_CUSTOM_ORDER} molecules of one species"
                ),
            ));
        }
        let name = format!("S{species}");
        let shifted = |j: u64| {
            if j == 0 {
                name.clone()
            } else {
                format!("({name} - {j}.0)")
            }
        };
        if whole {
            // 1 while at least `count` molecules are present, else 0.
            for j in 0..count {
                let x = shifted(j);
                factors.push(format!("({x}/({x} + 1e-300))"));
            }
        } else if count == 1 {
            factors.push(name.clone());
        } else {
            // The number of ways to choose `count` molecules: n (n - 1) … / count!.
            let falling: Vec<String> = (0..count).map(shifted).collect();
            let factorial: f64 = (1..=count).map(|k| k as f64).product();
            factors.push(format!("(({}) / {factorial:?})", falling.join(" * ")));
        }
    }
    Ok(factors.join(" * "))
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
    let view = match FlatView::of_any_model(model_bytes) {
        Ok(view) => view,
        Err(diagnostics) => {
            report.diagnostics = diagnostics;
            return report;
        }
    };
    if view.has_stop() {
        report.diagnostics.push(unsupported(
            "/stop",
            "GillesPy2 cannot preserve the model's stop condition",
        ));
    }
    let mut laws = Vec::with_capacity(view.channels());
    for c in 0..view.channels() {
        let rule = view.rule(c);
        let order: u64 = view.reactants(c).iter().map(|&(_, k)| k).sum();
        let constant = view.rate_constant(rule).filter(|_| !view.whole(c));
        match constant {
            Some(rate) if rate < 0.0 || !rate.is_finite() => report.diagnostics.push(unsupported(
                &format!("/rules/{rule}/rate"),
                "a negative constant rate cannot be exported",
            )),
            Some(rate) if order <= 2 => {
                let coefficient = rate * view.factor(c);
                // GillesPy2's homodimer uses k*x*(x-1)/volume, without the /2 factorial.
                let homodimer = view.reactants(c).iter().any(|&(_, k)| k == 2);
                let translated = if homodimer {
                    coefficient / 2.0
                } else {
                    coefficient
                };
                if homodimer && translated * 2.0 != coefficient {
                    report.diagnostics.push(Diagnostic::error("CRN-E051", &format!("/rules/{rule}/rate"),
                        "homodimer coefficient cannot be halved exactly in f64",
                        "Rescale the reference fixture; do not silently round or disable the exported channel."));
                }
                laws.push(Law::MassAction(translated));
            }
            _ => match custom(&view, c) {
                Ok(text) => laws.push(Law::Custom(text)),
                Err(diagnostic) => report.diagnostics.push(diagnostic),
            },
        }
    }
    if report
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error)
    {
        return report;
    }
    let mut source = format!(
        "# Generated by CRN-craft. Stochastic-count semantics; reference SSA only.\n# Model SHA-256: {}\nimport gillespy2 as gp\n\ndef build_model():\n    if gp.__version__ != {GILLESPY2_VERSION:?}:\n        raise RuntimeError('Requires GillesPy2 {GILLESPY2_VERSION}')\n    model = gp.Model(name={:?}, population=True, volume=1.0)\n    species = [\n",
        sha256(model_bytes),
        view.model().name
    );
    for (i, count) in view.initial().iter().enumerate() {
        source.push_str(&format!(
            "        gp.Species(name='S{i}', initial_value={count}, mode='discrete'),\n"
        ));
    }
    source.push_str("    ]\n    model.add_species(species)\n");
    // Terms in label order, as 0.6.0 listed them by species name.
    let side = |terms: &[(usize, u64)]| {
        let mut terms = terms.to_vec();
        terms.sort_by(|a, b| view.species()[a.0].cmp(&view.species()[b.0]));
        let text: Vec<_> = terms
            .iter()
            .map(|(s, c)| format!("species[{s}]: {c}"))
            .collect();
        format!("{{{}}}", text.join(", "))
    };
    for (i, law) in laws.iter().enumerate() {
        let (reactants, products) = (side(view.reactants(i)), side(view.products(i)));
        match law {
            Law::MassAction(coefficient) => source.push_str(&format!("    rate = gp.Parameter(name='K{i}', expression={coefficient:?})\n    model.add_parameter(rate)\n    model.add_reaction(gp.Reaction(name='R{i}', reactants={reactants}, products={products}, rate=rate))\n")),
            Law::Custom(text) => source.push_str(&format!("    model.add_reaction(gp.Reaction(name='R{i}', reactants={reactants}, products={products}, propensity_function={text:?}))\n")),
        }
    }
    source.push_str("    return model\n");
    let species_aliases = view
        .species()
        .iter()
        .enumerate()
        .map(|(i, label)| (label.clone(), format!("S{i}")))
        .collect();
    let mut per_rule = vec![0_usize; view.rules()];
    for c in 0..view.channels() {
        per_rule[view.rule(c)] += 1;
    }
    let mut seen = vec![0_usize; view.rules()];
    let reaction_aliases = (0..view.channels())
        .map(|c| {
            let rule = view.rule(c);
            let name = view.rule_name(rule);
            let key = if per_rule[rule] == 1 {
                name.to_string()
            } else {
                format!("{name}[{}]", seen[rule])
            };
            seen[rule] += 1;
            (key, format!("R{c}"))
        })
        .collect();
    let mut notes = vec![
        "Target volume is fixed at 1: source coefficients already use stochastic counts; Omega does not rescale them.".into(),
        "Homodimer rate constants are divided by 2 to preserve the binomial factorial. Parameters are resolved into one constant per channel.".into(),
        "Use a stochastic solver. Target ODE semantics, live budgets, i64 range/underflow diagnostics and cross-runtime bit identity are not preserved.".into(),
        "The target uses floating-point populations; even supported initial counts may exceed its exact range during a reference run.".into(),
    ];
    if laws.iter().any(|law| matches!(law, Law::Custom(_))) {
        notes.push("Custom propensities inline parameter values, write observables as sums of species, reactant combinations as falling factorials, and the availability of a whole-rule propensity's reactants as x/(x + 1e-300) factors; GillesPy2 evaluates them with Python or C++ floating point (a division by a zero count raises ZeroDivisionError in Python).".into());
    }
    report.artifact = Some(PythonExport {
        target: "gillespy2".into(),
        target_version: GILLESPY2_VERSION.into(),
        source,
        species_aliases,
        reaction_aliases,
        notes,
    });
    report.success = true;
    report
}
