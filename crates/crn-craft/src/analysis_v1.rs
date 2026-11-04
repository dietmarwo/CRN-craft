//! Versioned, JSON-friendly analysis and named-state input.
//!
//! Legacy `crn-craft.analysis/0.1` (`analyze-v1`), kept unchanged for existing
//! consumers. It reports structure, propensities and the linear conservation
//! basis, without the cone, linkage and effective bounds of the current
//! [`crate::analysis`] report, and drops `CRN-W001` only for species bounded
//! by the basis. New code should use [`crate::analyze_bytes`].
use crate::{
    conservation::{ConservationReport, analyze_conservation},
    diagnostic::Diagnostic,
    network::{CompiledNetwork, ReactionStructure},
    parse_model, sha256, strict_json,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of a named-state document, as the legacy analysis reads it.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum StateVersion {
    #[doc = concat!("`crn-craft.state/0.1`.")]
    #[serde(rename = "crn-craft.state/0.1")]
    V01,
}
/// Schema tag of the legacy analysis report.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum AnalysisVersion {
    #[doc = concat!("`crn-craft.analysis/0.1`.")]
    #[serde(rename = "crn-craft.analysis/0.1")]
    V01,
}

/// A complete state given by species name; the same document as [`crate::analysis::NamedState`].
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NamedState {
    #[doc = concat!("Always `crn-craft.state/0.1`.")]
    pub schema: StateVersion,
    /// Count of every declared species, by name, with no missing or extra names.
    pub counts: BTreeMap<String, i64>,
}

/// Structure, rates and linear conservation of a successfully compiled model (legacy 0.1).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NetworkAnalysis {
    /// Species names in declaration order; this order indexes every other per-species list.
    pub species: Vec<String>,
    /// Initial counts from the model.
    pub initial_counts: Vec<i64>,
    /// The state at which propensities were evaluated: the supplied state, or the initial
    /// counts.
    pub evaluation_counts: Vec<i64>,
    /// Indexed reactions in declaration order.
    pub reactions: Vec<ReactionStructure>,
    /// Sorted propensity invalidations for each fired reaction index.
    pub dependencies: Vec<Vec<usize>>,
    /// Propensity of each reaction at `evaluation_counts`, in events per second.
    pub propensities: Vec<f64>,
    #[doc = concat!("Linear conservation basis and its bounds.")]
    pub conservation: ConservationReport,
}

/// Result of `analyze-v1`: success flag, input hashes, the legacy analysis and every diagnostic.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnalysisReport {
    #[doc = concat!("Always `crn-craft.analysis/0.1`.")]
    pub schema: AnalysisVersion,
    /// True when the model compiled and every propensity evaluated.
    pub success: bool,
    /// Always false: analysis runs no simulation.
    pub simulation_available: bool,
    /// SHA-256 of the exact model bytes; null only when the model file could not be read.
    pub model_sha256: Option<String>,
    /// SHA-256 of the exact state bytes when a state was supplied.
    pub state_sha256: Option<String>,
    /// The analysis when `success` is true, otherwise null.
    pub analysis: Option<NetworkAnalysis>,
    /// Model and conservation warnings on success, the diagnostics of the failed step otherwise.
    pub diagnostics: Vec<Diagnostic>,
}

impl AnalysisReport {
    /// A failed report carrying one diagnostic, such as a read error (`CRN-E030`).
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: AnalysisVersion::V01,
            success: false,
            simulation_available: false,
            model_sha256: None,
            state_sha256: None,
            analysis: None,
            diagnostics: vec![error],
        }
    }
}

/// Legacy analysis of a model from its exact bytes, optionally at a
/// `crn-craft.state/0.1` state. Failures give `success: false`.
pub fn analyze_bytes(model_bytes: &[u8], state_bytes: Option<&[u8]>) -> AnalysisReport {
    let mut report = AnalysisReport {
        schema: AnalysisVersion::V01,
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
        let bounded_paths: Vec<String> = conservation
            .initial_class_upper_bounds
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
        Ok(NetworkAnalysis {
            species: network.species().to_vec(),
            initial_counts: network.initial_counts().to_vec(),
            evaluation_counts: state,
            reactions: network.reactions().to_vec(),
            dependencies: network.dependencies().to_vec(),
            propensities,
            conservation,
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
