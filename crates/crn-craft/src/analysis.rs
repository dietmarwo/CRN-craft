//! Versioned, JSON-friendly network analysis of a closed model.
//!
//! [`analyze_bytes`] (and [`crate::analyze_files`]) reads a model of any
//! version the engine reads (a network `crn-craft.model/0.1`, converted
//! first by [`crate::rules::convert::from_network`], a rules model or a
//! `crn-craft.model/1.0`), lowers it to its flat network
//! ([`FlatView::with_every_type`]), evaluates every propensity at the
//! initial counts or at a supplied `crn-craft.state/0.1` document, and
//! reports the structure, dependencies, linear conservation laws, the
//! non-negative conservation cone and the complex graph as
//! `crn-craft.analysis/1.0`. It runs no simulation.
//!
//! The species are the view's species labels (for a converted network model
//! the network's species, in declaration order); the reactions are the
//! view's channels (for a converted network model, one per reaction, in
//! declaration order). A model that has no flat view (ATP, continuous
//! variables, products that create new attribute tuples) is refused with
//! `CRN-R017`.
//!
//! ```
//! use crn_craft::analyze_bytes;
//!
//! let model = br#"{"schema": "crn-craft.model/0.1", "name": "isomer",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
//!   "species": [{"name": "A", "initial_count": 7}, {"name": "B", "initial_count": 3}],
//!   "reactions": [
//!     {"name": "forward", "reactants": {"A": 1}, "products": {"B": 1},
//!      "rate": {"kind": "mass_action", "coefficient": 2.0}},
//!     {"name": "back", "reactants": {"B": 1}, "products": {"A": 1},
//!      "rate": {"kind": "mass_action", "coefficient": 1.0}}]}"#;
//! let state = br#"{"schema": "crn-craft.state/0.1", "counts": {"A": 1, "B": 9}}"#;
//! let report = analyze_bytes(model, Some(state));
//! assert!(report.success, "{:?}", report.diagnostics);
//! let analysis = report.analysis.unwrap();
//! assert_eq!(analysis.propensities, vec![2.0, 9.0]);
//! // A + B is conserved, so each count is at most the initial total 10.
//! let law = &analysis.conservation.basis[0];
//! assert_eq!(law.coefficients, ["1", "1"]);
//! assert_eq!(law.initial_total, "10");
//! assert_eq!(analysis.initial_class_upper_bounds, vec![Some("10".into()), Some("10".into())]);
//! assert_eq!(analysis.linkage.deficiency, Some(0));
//! ```
use crate::{
    cone::{ConeReport, ConeStatus, analyze_cone},
    conservation::{ConservationReport, analyze_conservation},
    diagnostic::{Diagnostic, pointer},
    linkage::{LinkageReport, analyze_linkage},
    rules::{FlatView, model::MAX_COUNT},
    sha256, strict_json,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of a named-state document.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum StateVersion {
    #[doc = concat!("`crn-craft.state/0.1`.")]
    #[serde(rename = "crn-craft.state/0.1")]
    V01,
}
/// Schema tag of the analysis report.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum AnalysisVersion {
    #[doc = concat!("`crn-craft.analysis/1.0`.")]
    #[serde(rename = "crn-craft.analysis/1.0")]
    V10,
}

/// A complete state given by species label, the optional input of `analyze`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NamedState {
    #[doc = concat!("Always `crn-craft.state/0.1`.")]
    pub schema: StateVersion,
    /// Count of every species, by label, with no missing or extra labels. Counts are
    /// non-negative and at most 2^53 (`CRN-E041`).
    pub counts: BTreeMap<String, i64>,
}

/// Sparse entry. Indices refer to the report's `species` list.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Term {
    /// Zero-based species index.
    pub species: usize,
    /// Stoichiometric coefficient: positive on reactant and product sides, signed in net changes.
    pub coefficient: i64,
}

impl Term {
    /// Terms of a (species, count) list, in the list's order.
    pub(crate) fn list<T: Copy + TryInto<i64>>(entries: &[(usize, T)]) -> Vec<Self> {
        entries
            .iter()
            .map(|&(species, count)| Self {
                species,
                coefficient: count.try_into().unwrap_or(i64::MAX),
            })
            .collect()
    }
}

/// Indexed structure of one reaction (a channel of the flat view), with terms sorted by
/// species index.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReactionStructure {
    /// The name of the reaction's rule; the channels of a rule with several reactant
    /// tuples share it.
    pub name: String,
    /// Consumed species and coefficients; they also gate the channel.
    pub reactants: Vec<Term>,
    /// Produced species and coefficients.
    pub products: Vec<Term>,
    /// One stoichiometric matrix column: products minus reactants, zeros omitted.
    pub net_changes: Vec<Term>,
    /// Conservative reads, including reactant gates and the species of every observable
    /// the rate reads.
    pub reads: Vec<usize>,
}

/// Structure, rates and conservation analysis of a successfully lowered model.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NetworkAnalysis {
    /// Species labels of the flat view (for a network model its species names in
    /// declaration order); this order indexes every other per-species list.
    pub species: Vec<String>,
    /// Initial counts from the model.
    pub initial_counts: Vec<i64>,
    /// The state at which propensities were evaluated: the supplied state, or the initial
    /// counts. It need not lie in the initial conservation class.
    pub evaluation_counts: Vec<i64>,
    /// Indexed reactions (the channels of the flat view) in rule declaration order.
    pub reactions: Vec<ReactionStructure>,
    /// Sorted propensity invalidations for each fired reaction index.
    pub dependencies: Vec<Vec<usize>>,
    /// Propensity of each reaction at `evaluation_counts`, in events per unit time.
    pub propensities: Vec<f64>,
    #[doc = concat!("Linear conservation basis and its bounds.")]
    pub conservation: ConservationReport,
    #[doc = concat!("Extreme rays of the non-negative conservation cone and their bounds.")]
    pub nonnegative_conservation: ConeReport,
    #[doc = concat!("Complex graph, linkage classes and deficiency.")]
    pub linkage: LinkageReport,
    /// Cone bounds when exact, otherwise sound linear-basis bounds.
    /// Null means not proved, not unbounded.
    pub initial_class_upper_bounds: Vec<Option<String>>,
}

/// Result of `analyze`: success flag, input hashes, the analysis and every diagnostic.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnalysisReport {
    #[doc = concat!("Always `crn-craft.analysis/1.0`.")]
    pub schema: AnalysisVersion,
    /// True when the model was lowered and every propensity evaluated. Conservation or
    /// cone limits leave it true and show up in their status fields.
    pub success: bool,
    /// Always false: analysis runs no simulation.
    pub simulation_available: bool,
    /// SHA-256 of the exact model bytes; null only when the model file could not be read.
    pub model_sha256: Option<String>,
    /// SHA-256 of the exact state bytes when a state was supplied.
    pub state_sha256: Option<String>,
    /// The analysis when `success` is true, otherwise null.
    pub analysis: Option<NetworkAnalysis>,
    /// On success the structural warnings (`CRN-W001` dropped for species with a proven
    /// bound, `CRN-W002`, `CRN-W003`) and the conservation and cone warnings; on failure
    /// the diagnostics of the failed step.
    pub diagnostics: Vec<Diagnostic>,
}

impl AnalysisReport {
    /// A failed report carrying one diagnostic, such as a read error (`CRN-E030`).
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: AnalysisVersion::V10,
            success: false,
            simulation_available: false,
            model_sha256: None,
            state_sha256: None,
            analysis: None,
            diagnostics: vec![error],
        }
    }
}

/// Analyzes a model from its exact bytes, optionally at the state in
/// `state_bytes` (a `crn-craft.state/0.1` document); see the
/// [module docs](self). Never panics on bad input: failures are reported
/// with `success: false` and their diagnostics.
pub fn analyze_bytes(model_bytes: &[u8], state_bytes: Option<&[u8]>) -> AnalysisReport {
    let mut report = AnalysisReport {
        schema: AnalysisVersion::V10,
        success: false,
        simulation_available: false,
        model_sha256: Some(sha256(model_bytes)),
        state_sha256: state_bytes.map(sha256),
        analysis: None,
        diagnostics: Vec::new(),
    };
    let result = FlatView::of_any_model(model_bytes).and_then(|view| {
        let state = match state_bytes {
            Some(bytes) => {
                let state: NamedState = strict_json::parse(bytes).map_err(|d| vec![d])?;
                named_state(&view, &state.counts).map_err(|d| vec![d])?
            }
            None => view.initial().to_vec(),
        };
        let propensities = view.propensities(&state).map_err(|d| vec![d])?;
        Ok(analyze_view(&view, &state, propensities))
    });
    match result {
        Ok((analysis, diagnostics)) => {
            report.analysis = Some(analysis);
            report.diagnostics = diagnostics;
            report.success = true;
        }
        Err(diagnostics) => report.diagnostics = diagnostics,
    }
    report
}

/// The analysis of a view at `state`, with its warnings.
fn analyze_view(
    view: &FlatView,
    state: &[u64],
    propensities: Vec<f64>,
) -> (NetworkAnalysis, Vec<Diagnostic>) {
    let conservation = analyze_conservation(view);
    let nonnegative_conservation = analyze_cone(view);
    let linkage = analyze_linkage(view, conservation.rank);
    let bounds = if nonnegative_conservation.status == ConeStatus::Exact {
        nonnegative_conservation.initial_class_upper_bounds.clone()
    } else {
        conservation.initial_class_upper_bounds.clone()
    };
    let bounded_paths: Vec<String> = bounds
        .iter()
        .enumerate()
        .filter(|(_, bound)| bound.is_some())
        .map(|(i, _)| format!("/species/{i}"))
        .collect();
    let mut diagnostics: Vec<Diagnostic> = warnings(view)
        .into_iter()
        .filter(|d| d.code != "CRN-W001" || !bounded_paths.contains(&d.path))
        .collect();
    diagnostics.extend(conservation.diagnostics.iter().cloned());
    diagnostics.extend(nonnegative_conservation.diagnostics.iter().cloned());
    // Counts are at most 2^53.
    let counts = |list: &[u64]| -> Vec<i64> { list.iter().map(|&n| n as i64).collect() };
    let reactions = (0..view.channels())
        .map(|c| ReactionStructure {
            name: view.rule_name(view.rule(c)).to_string(),
            reactants: Term::list(view.reactants(c)),
            products: Term::list(view.products(c)),
            net_changes: Term::list(view.changes(c)),
            reads: view.reads(c).to_vec(),
        })
        .collect();
    let analysis = NetworkAnalysis {
        species: view.species().to_vec(),
        initial_counts: counts(view.initial()),
        evaluation_counts: counts(state),
        reactions,
        dependencies: (0..view.channels())
            .map(|c| view.dependencies(c).to_vec())
            .collect(),
        propensities,
        conservation,
        nonnegative_conservation,
        linkage,
        initial_class_upper_bounds: bounds,
    };
    (analysis, diagnostics)
}

/// The structural warnings of the network engine, on the view's lists:
/// per reaction `CRN-W003` (no net change, at `/reactions/i`) and
/// `CRN-W002` (a rate that is the constant 0, at `/reactions/i`), then per
/// species `CRN-W001` (net produced by some reaction and net consumed by
/// none, at `/species/i`).
fn warnings(view: &FlatView) -> Vec<Diagnostic> {
    let n = view.species().len();
    let (mut produced, mut consumed) = (vec![false; n], vec![false; n]);
    let mut out = Vec::new();
    for c in 0..view.channels() {
        let path = format!("/reactions/{c}");
        if view.changes(c).is_empty() {
            out.push(Diagnostic::warning(
                "CRN-W003",
                &path,
                "reaction has no population effect",
                "Remove it unless event counting has a stated purpose.",
            ));
        }
        if view.rate_constant(view.rule(c)) == Some(0.0) {
            out.push(Diagnostic::warning(
                "CRN-W002",
                &path,
                "reaction rate is zero",
                "Confirm that this disabled reaction is intentional.",
            ));
        }
        for &(species, change) in view.changes(c) {
            if change > 0 {
                produced[species] = true;
            } else {
                consumed[species] = true;
            }
        }
    }
    for species in 0..n {
        if produced[species] && !consumed[species] {
            out.push(Diagnostic::warning(
                "CRN-W001",
                &format!("/species/{species}"),
                "species has a net-producing reaction and no net-consuming reaction",
                "Review growth and budgets; this structural warning does not prove unboundedness.",
            ));
        }
    }
    out
}

/// A state by label in view order; `CRN-E041` for a missing or unknown
/// label or a count outside 0..=2^53.
fn named_state(view: &FlatView, counts: &BTreeMap<String, i64>) -> Result<Vec<u64>, Diagnostic> {
    let invalid = |label: &str, message: &str, hint: &str| {
        Diagnostic::error("CRN-E041", &pointer("/counts", label), message, hint)
    };
    let mut state = Vec::with_capacity(view.species().len());
    for label in view.species() {
        let count = counts.get(label).ok_or_else(|| {
            invalid(
                label,
                "missing species count",
                "Provide every species label exactly once.",
            )
        })?;
        state.push(*count);
    }
    if let Some(label) = counts.keys().find(|label| view.index(label).is_none()) {
        return Err(invalid(
            label,
            "unknown species in state",
            "Use the species labels of the analysis report.",
        ));
    }
    state
        .iter()
        .zip(view.species())
        .map(|(&count, label)| {
            u64::try_from(count)
                .ok()
                .filter(|&n| n <= MAX_COUNT)
                .ok_or_else(|| {
                    invalid(
                        label,
                        "invalid count",
                        "Use non-negative counts of at most 2^53.",
                    )
                })
        })
        .collect()
}
