//! Versioned, JSON-friendly analysis and named-state input.
//!
//! [`analyze_bytes`] (and [`crate::analyze_files`]) parse and compile a
//! model, evaluate every propensity at the initial counts or at a supplied
//! `crn-craft.state/0.1` document, and report the structure, dependencies,
//! linear conservation laws, the non-negative conservation cone and the
//! complex graph as `crn-craft.analysis/0.2`. It runs no simulation. The
//! legacy 0.1 report is in [`crate::analysis_v1`].
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
    diagnostic::Diagnostic,
    linkage::{LinkageReport, analyze_linkage},
    network::{CompiledNetwork, ReactionStructure},
    parse_model, sha256, strict_json,
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
/// Schema tag of the current analysis report.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum AnalysisVersion {
    #[doc = concat!("`crn-craft.analysis/0.2`.")]
    #[serde(rename = "crn-craft.analysis/0.2")]
    V02,
}

/// A complete state given by species name, the optional input of `analyze`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NamedState {
    #[doc = concat!("Always `crn-craft.state/0.1`.")]
    pub schema: StateVersion,
    /// Count of every declared species, by name, with no missing or extra names. Counts are
    /// non-negative and within each population limit (`CRN-E041`).
    pub counts: BTreeMap<String, i64>,
}

/// Structure, rates and conservation analysis of a successfully compiled model.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NetworkAnalysis {
    /// Species names in declaration order; this order indexes every other per-species list.
    pub species: Vec<String>,
    /// Initial counts from the model.
    pub initial_counts: Vec<i64>,
    /// The state at which propensities were evaluated: the supplied state, or the initial
    /// counts. It need not lie in the initial conservation class.
    pub evaluation_counts: Vec<i64>,
    /// Indexed reactions in declaration order.
    pub reactions: Vec<ReactionStructure>,
    /// Sorted propensity invalidations for each fired reaction index.
    pub dependencies: Vec<Vec<usize>>,
    /// Propensity of each reaction at `evaluation_counts`, in events per second.
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
    #[doc = concat!("Always `crn-craft.analysis/0.2`.")]
    pub schema: AnalysisVersion,
    /// True when the model compiled and every propensity evaluated. Conservation or cone
    /// limits leave it true and show up in their status fields.
    pub success: bool,
    /// Always false: analysis runs no simulation.
    pub simulation_available: bool,
    /// SHA-256 of the exact model bytes; null only when the model file could not be read.
    pub model_sha256: Option<String>,
    /// SHA-256 of the exact state bytes when a state was supplied.
    pub state_sha256: Option<String>,
    /// The analysis when `success` is true, otherwise null.
    pub analysis: Option<NetworkAnalysis>,
    /// On success the model warnings (`CRN-W001` dropped for species with a proven bound)
    /// and the conservation and cone warnings; on failure the diagnostics of the failed step.
    pub diagnostics: Vec<Diagnostic>,
}

impl AnalysisReport {
    /// A failed report carrying one diagnostic, such as a read error (`CRN-E030`).
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: AnalysisVersion::V02,
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
        schema: AnalysisVersion::V02,
        success: false,
        simulation_available: false,
        model_sha256: Some(sha256(model_bytes)),
        state_sha256: state_bytes.map(sha256),
        analysis: None,
        diagnostics: Vec::new(),
    };
    let result = (|| {
        let model = parse_model(model_bytes).map_err(|d| vec![d])?;
        let network = CompiledNetwork::compile(&model)?;
        let state = if let Some(bytes) = state_bytes {
            let state: NamedState = strict_json::parse(bytes).map_err(|d| vec![d])?;
            network.named_state(&state.counts).map_err(|d| vec![d])?
        } else {
            network.initial_counts().to_vec()
        };
        let propensities = network.propensities(&state).map_err(|d| vec![d])?;
        let conservation = analyze_conservation(&network);
        let nonnegative_conservation = analyze_cone(&network);
        let linkage = analyze_linkage(&network, conservation.rank);
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
        report.diagnostics = network
            .warnings()
            .iter()
            .filter(|d| d.code != "CRN-W001" || !bounded_paths.contains(&d.path))
            .cloned()
            .collect();
        report
            .diagnostics
            .extend(conservation.diagnostics.iter().cloned());
        report
            .diagnostics
            .extend(nonnegative_conservation.diagnostics.iter().cloned());
        Ok(NetworkAnalysis {
            species: network.species().to_vec(),
            initial_counts: network.initial_counts().to_vec(),
            evaluation_counts: state,
            reactions: network.reactions().to_vec(),
            dependencies: network.dependencies().to_vec(),
            propensities,
            conservation,
            nonnegative_conservation,
            linkage,
            initial_class_upper_bounds: bounds,
        })
    })();
    match result {
        Ok(analysis) => {
            report.analysis = Some(analysis);
            report.success = true;
        }
        Err(diagnostics) => {
            report.diagnostics = diagnostics;
        }
    }
    report
}
