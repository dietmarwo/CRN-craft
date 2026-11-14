#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![recursion_limit = "256"]
//! Exact and approximate stochastic simulation, finite state projection and
//! analysis of chemical reaction networks, from fixed species lists to
//! dynamic, attributed species.
//!
//! `crn-craft` is the engine of [CRN-craft](https://github.com/dietmarwo/CRN-craft).
//! Models, requests and results are versioned JSON documents whose schemas
//! are generated from the Rust types in this crate. Inputs are read strictly
//! (duplicate keys and unknown fields are errors), bound to their exact bytes
//! by SHA-256 and validated with stable diagnostic codes. Every run states its
//! budgets and seed and accounts for every replicate. Completed results
//! replay exactly from the same input bytes, method and seed, whatever the
//! thread count.
//!
//! # Status
//!
//! The engine is [`rules`]: models are `crn-craft.model/1.0` (and the older
//! `crn-craft.rules-model` versions, which keep replaying byte for byte),
//! requests `crn-craft.run/1.0` (and `crn-craft.rules-run`). Exact methods
//! are Static Direct with four event selectors (and `auto`), the tree method,
//! next reaction, the mtDNA path and thinning for the ATP block; tau leaping,
//! binomial steps and the hybrid method are explicitly approximate. Beyond
//! that:
//!
//! - [`fsp`] reports its truncation and leakage terms separately; they are
//!   estimates, not a certified floating-point error bound.
//! - Generated event loops and kernels ([`rules::generate_rules_loop`]) are
//!   opt-in; benchmark them before adopting them.
//! - [`kinetics`] serves the optional deterministic (ODE) integration, whose
//!   solutions are not exact means of nonlinear stochastic models.
//! - A network `crn-craft.model/0.1` is read through the converter
//!   ([`rules::convert::from_network`]) by check, analysis, export, the
//!   kinetics and `crn-craft convert`.
//!
//! A successful numerical run is not biological validation.
//!
//! # Quickstart
//!
//! Check a model, then simulate an ensemble with an exact method:
//!
//! ```
//! use crn_craft::rules::{Control, Termination, check_model, solve};
//! use crn_craft::sha256;
//!
//! // Immigration–death: X is made at rate 4 and each molecule decays at rate 0.2.
//! let model = br#"{
//!   "schema": "crn-craft.model/1.0", "name": "immigration_death",
//!   "parameters": {"birth": 4.0, "death": 0.2},
//!   "types": [{"name": "X"}], "initial": [],
//!   "observables": {"X": {"type": "X"}},
//!   "rules": [
//!     {"name": "birth", "reactants": [], "rate": "birth", "products": [{"type": "X"}]},
//!     {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "death",
//!      "products": []}]}"#;
//! let report = check_model(model);
//! assert!(report.valid, "{:?}", report.diagnostics);
//!
//! // The request binds the exact model bytes by SHA-256.
//! let request = format!(r#"{{
//!   "schema": "crn-craft.run/1.0",
//!   "model": {{"path": "model.json", "sha256": "{}"}},
//!   "method": {{"kind": "static_direct"}},
//!   "time": {{"end": 10.0, "report_times": [0.0, 5.0, 10.0]}},
//!   "ensemble": {{"replicates": 32, "root_seed": 42}},
//!   "outputs": {{"statistics": ["mean", "variance"], "keep_paths": 0, "keep_instances": 0}},
//!   "budget": {{"max_events_per_replicate": 100000,
//!              "max_sampling_steps_per_replicate": 1000000, "max_instances": 1000,
//!              "threads": 1, "max_output_bytes": 1000000}}}}"#, sha256(model));
//!
//! let result = solve(model, request.as_bytes(), &Control::default()).unwrap();
//! assert_eq!(result.termination, Termination::Completed);
//! let statistics = result.statistics.expect("every replicate completed");
//! // mean[report time][observable]; the exact mean at t = 10 is 20 (1 − e⁻²) ≈ 17.3.
//! let mean = statistics.mean.expect("the request asks for the mean");
//! assert!((mean[2][0] - 17.3).abs() < 3.0);
//! ```
//!
//! # Entry points
//!
//! | Task | API | Contracts (`crn-craft.` prefix omitted) |
//! |---|---|---|
//! | Check models and requests | [`check_any`], [`check_any_files`], [`rules::check_model`] | any model, `run/1.0`, rules requests → `check/1.0`, `rules-check/0.1` |
//! | Simulate | [`rules::prepare`], [`rules::execute`], [`rules::solve`], [`rules::solve_files`] | `model/1.0`, rules models; `run/1.0`, `rules-run/0.1`–`0.6` → `result/1.0`, `rules-result/0.1`–`0.3` |
//! | Convert network models | [`rules::convert::from_network`] | `model/0.1` → `model/1.0` |
//! | Analyze structure and conservation | [`analyze_bytes`], [`analyze_files`], [`rules::FlatView`] | `model/0.1`, rules models, `model/1.0`, `state/0.1` → `analysis/1.0` |
//! | Sealed archives | [`archive::v1::run_files`], [`archive::v1::write`], [`archive::v1::verify`] | rules requests, `run/1.0` → `archive/1.0`, `archive-check/1.0` |
//! | Finite state projection | [`fsp::v1::solve`], [`fsp::v1::solve_files`], [`rules::FlatView`] | `run/1.0` with `fsp` → `fsp-result/1.0` |
//! | Reference export | [`export::export_gillespy2`], [`export_gillespy2_file`] | `model/0.1`, rules models, `model/1.0` → `export/0.1` |
//! | Deterministic kinetics and Jacobians | [`kinetics::Kinetics`], [`rules::FlatView`] | `model/0.1`, rules models, `model/1.0` |
//! | Discovery | [`describe`], [`schema`], [`rules::schema`] | `describe/1.0` |
//!
//! Lower-level building blocks are [`strict_json`], [`diagnostic`],
//! [`sampling`], [`conservation`], [`cone`] and [`linkage`].
//!
//! # Errors, stops and cancellation
//!
//! Errors are [`diagnostic::Diagnostic`] values with a stable code (`CRN-E…`,
//! `CRN-R…`, warnings `CRN-W…`) and a JSON Pointer to the offending input,
//! listed in
//! [docs/reference/diagnostics.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/diagnostics.md).
//! Report-style functions ([`check_any`], [`analyze_bytes`],
//! [`export::export_gillespy2`], [`archive::v1::run_files`],
//! [`archive::v1::verify`]) never fail; they return a report with a `valid`
//! or `success` flag. In stochastic simulations, event, sampling, instance
//! and wall-time limits and cancellation are typed outcomes in each
//! replicate's record, not errors; an exceeded output budget is the run-level
//! outcome `output_limit`; and ensemble statistics are withheld unless every
//! replicate completed. FSP has no replicates and fails with `CRN-E070`
//! instead. Cancel a running call through [`rules::Control`].
//!
//! # Related crates and documentation
//!
//! - [`crn-craft-cli`](https://crates.io/crates/crn-craft-cli): the
//!   `crn-craft` command-line tool.
//! - [`crn-craft-rules`](https://crates.io/crates/crn-craft-rules): the
//!   former rules-engine crate, now a deprecated re-export of [`rules`].
//! - [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/): Python
//!   bindings with prepared runs, importers and trajectory analysis.
//! - Repository documentation: the
//!   [user guide](https://github.com/dietmarwo/CRN-craft/blob/main/docs/guide/01-install.md),
//!   [architecture](https://github.com/dietmarwo/CRN-craft/blob/main/docs/developer/architecture.md),
//!   [model language](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/model-language.md),
//!   [input contracts](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/contracts.md),
//!   [result and archive contract](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/archives.md)
//!   and the [technical report](https://github.com/dietmarwo/CRN-craft/blob/main/docs/evidence/technical-report.md)
//!   on numerical methods and validation.

// Contributor note: doc comments on types that derive `JsonSchema` become
// "description" entries in the generated schemas (`schemas/*.schema.json`,
// `describe`). schemars moves a documented unit enum variant into its own
// `oneOf` branch and wraps a documented field whose type is a schema
// definition in `allOf`, which changes a delivered schema's shape. Such
// items are documented with `#[doc = concat!("…")]`, which rustdoc reads and
// schemars ignores. After editing docs on schema types, regenerate the
// snapshots with `crn-craft describe KIND > schemas/KIND.schema.json` and
// check that only descriptions changed.
pub mod analysis;
pub mod archive;
pub mod cone;
pub mod conservation;
pub mod diagnostic;
pub mod export;
pub mod fsp;
pub mod kinetics;
pub mod linkage;
pub mod rules;
pub mod sampling;
pub mod ssa;
pub mod strict_json;
pub mod unified;

pub use analysis::{AnalysisReport, analyze_bytes};
pub use unified::{check_any, check_any_files, describe};

use diagnostic::Diagnostic;
use schemars::schema_for;
use serde::Deserialize;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::Path;

/// Lowercase hexadecimal SHA-256 of `bytes`: the digest every `sha256` field
/// and model reference uses. Hash the exact file bytes, never re-serialized
/// JSON.
pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// The top-level `schema` string of a JSON document, if it has one; read
/// leniently, for choosing a diagnostic or a report version.
pub(crate) fn schema_field(bytes: &[u8]) -> Option<String> {
    #[derive(Deserialize)]
    struct Peek {
        schema: Option<Value>,
    }
    serde_json::from_slice::<Peek>(bytes)
        .ok()?
        .schema?
        .as_str()
        .map(str::to_string)
}

/// Reads a model file (a network `crn-craft.model/0.1`, a rules model or a
/// `crn-craft.model/1.0`) and an optional `crn-craft.state/0.1` file, then
/// runs [`analyze_bytes`]. Each must be a regular file of at most 16 MiB;
/// read failures give a failed report with `CRN-E030` or `CRN-E004`.
pub fn analyze_files(model_path: &Path, state_path: Option<&Path>) -> AnalysisReport {
    let model = match read_document(model_path) {
        Ok(bytes) => bytes,
        Err(error) => return AnalysisReport::failure(error),
    };
    let state = match state_path.map(read_document).transpose() {
        Ok(bytes) => bytes,
        Err(error) => return AnalysisReport::failure(error),
    };
    analyze_bytes(&model, state.as_deref())
}

/// Reads a model file (a network `crn-craft.model/0.1`, a rules model or a
/// `crn-craft.model/1.0`) and runs [`export::export_gillespy2`]; read
/// failures give a failed report with `CRN-E030` or `CRN-E004`.
pub fn export_gillespy2_file(model_path: &Path) -> export::ExportReport {
    match read_document(model_path) {
        Ok(bytes) => export::export_gillespy2(&bytes),
        Err(error) => export::ExportReport::failure(error),
    }
}

fn read_document(path: &Path) -> Result<Vec<u8>, Diagnostic> {
    let io_error = |error: std::io::Error| {
        Diagnostic::error(
            "CRN-E030",
            "",
            format!("{}: {error}", path.display()),
            "Provide a readable regular JSON file.",
        )
    };
    if !std::fs::metadata(path).map_err(io_error)?.is_file() {
        return Err(Diagnostic::error(
            "CRN-E030",
            "",
            format!("{} is not a regular file", path.display()),
            "Provide a regular JSON file.",
        ));
    }
    let file = std::fs::File::open(path).map_err(io_error)?;
    if !file.metadata().map_err(io_error)?.is_file() {
        return Err(Diagnostic::error(
            "CRN-E030",
            "",
            format!("{} is not a regular file", path.display()),
            "Provide a regular JSON file.",
        ));
    }
    let mut bytes = Vec::new();
    file.take(strict_json::MAX_DOCUMENT_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(io_error)?;
    if bytes.len() > strict_json::MAX_DOCUMENT_BYTES {
        return Err(Diagnostic::error(
            "CRN-E004",
            "",
            "document exceeds 16 MiB",
            "Reduce the input file.",
        ));
    }
    Ok(bytes)
}

/// The generated JSON Schema of one contract, as `crn-craft describe KIND`
/// prints it, or `None` for an unknown kind.
///
/// Kinds: `model` (the network `crn-craft.model/0.1`, read through the
/// converter), `archive-1.0`, `archive-check-1.0`, `state`, `analysis-1.0`,
/// `export` and `fsp-result-1.0`. The `rules-*` kinds and the other 1.0
/// contracts (`model-1.0`, `run-1.0`, ...) come from [`rules::schema`].
///
/// ```
/// let schema = crn_craft::schema("model").unwrap();
/// assert_eq!(schema["title"], "Model");
/// assert_eq!(schema["properties"]["schema"]["$ref"], "#/definitions/ModelVersion");
/// assert!(crn_craft::schema("rules-model").is_none());
/// assert!(crn_craft::schema("run").is_none()); // the network run/0.1 is gone
/// ```
pub fn schema(kind: &str) -> Option<Value> {
    match kind {
        "model" => Some(rules::convert::network_schema()),
        "archive-1.0" => serde_json::to_value(schema_for!(archive::v1::Manifest)).ok(),
        "archive-check-1.0" => serde_json::to_value(schema_for!(archive::v1::ArchiveReport)).ok(),
        "state" => serde_json::to_value(schema_for!(analysis::NamedState)).ok(),
        "analysis-1.0" => serde_json::to_value(schema_for!(AnalysisReport)).ok(),
        "export" => serde_json::to_value(schema_for!(export::ExportReport)).ok(),
        "fsp-result-1.0" => serde_json::to_value(schema_for!(fsp::v1::FspResult)).ok(),
        _ => None,
    }
}
