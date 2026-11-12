#![forbid(unsafe_code)]
#![warn(missing_docs)]
//! Exact stochastic simulation, finite state projection and analysis of
//! well-mixed chemical reaction networks with fixed species and reactions.
//!
//! `crn-craft` is the fixed-network engine of
//! [CRN-craft](https://github.com/dietmarwo/CRN-craft). Models, requests and
//! results are versioned JSON documents whose schemas are generated from the
//! Rust types in this crate. Inputs are read strictly (duplicate keys and
//! unknown fields are errors), bound to their exact bytes by SHA-256 and
//! validated with stable diagnostic codes. Every run states its budgets and
//! seed and accounts for every replicate. Completed results replay exactly
//! from the same input bytes, sampler and seed, whatever the thread count.
//!
//! # Status
//!
//! Implemented: strict contracts and validation, network analysis, Direct
//! SSA with sealed archives, the alternative exact methods of [`ssa`]
//! (Direct, cached Direct, DirectCR, SortingDirect, RSSACR and next
//! reaction), finite state projection and reusable parameter templates.
//! Beyond that:
//!
//! - [`tau`] leaping is experimental and explicitly approximate; its best
//!   assessment is `approximate_ensemble`.
//! - [`fsp`] reports its truncation and leakage terms separately; they are
//!   estimates, not a certified floating-point error bound.
//! - Generated propensity kernels ([`codegen`]) are opt-in; benchmark them
//!   before adopting them.
//! - [`kinetics`] serves the optional deterministic (ODE) integration, whose
//!   solutions are not exact means of nonlinear stochastic models.
//! - [`analysis_v1`] and the reserved [`contracts::RunResult`] envelope are
//!   legacy contracts kept unchanged for existing consumers.
//!
//! A successful numerical run is not biological validation. Dynamic,
//! attributed species are the separate, experimental
//! [`crn-craft-rules`](https://crates.io/crates/crn-craft-rules) crate.
//!
//! # Quickstart
//!
//! Check a model, then simulate an ensemble with an exact method:
//!
//! ```
//! use crn_craft::{check_bytes, sha256, simulation::Control, ssa};
//!
//! // Immigration–death: X is made at rate 4 and each molecule decays at rate 0.2.
//! let model = br#"{
//!   "schema": "crn-craft.model/0.1", "name": "immigration_death",
//!   "convention": "stochastic_counts", "omega": 1.0,
//!   "parameters": {"birth": 4.0, "death": 0.2},
//!   "species": [{"name": "X", "initial_count": 0}],
//!   "reactions": [
//!     {"name": "birth", "reactants": {}, "products": {"X": 1},
//!      "rate": {"kind": "mass_action", "coefficient": "birth"}},
//!     {"name": "death", "reactants": {"X": 1}, "products": {},
//!      "rate": {"kind": "mass_action", "coefficient": "death"}}]}"#;
//! let report = check_bytes(model, None);
//! assert!(report.valid, "{:?}", report.diagnostics);
//!
//! // The request binds the exact model bytes by SHA-256.
//! let request = format!(r#"{{
//!   "schema": "crn-craft.ssa-run/0.4",
//!   "model": {{"path": "model.json", "sha256": "{}"}},
//!   "method": {{"kind": "ssa", "aggregator": "direct", "waiting_time": "inverse"}},
//!   "time": {{"end": 10.0, "save_times": [0.0, 5.0, 10.0]}},
//!   "ensemble": {{"replicates": 32, "root_seed": 42}},
//!   "outputs": {{"statistics": ["mean", "variance"], "keep_trajectories": 0}},
//!   "budget": {{"max_events_per_replicate": 100000,
//!              "max_sampling_steps_per_replicate": 1000000,
//!              "threads": 1, "max_output_bytes": 1000000}}}}"#, sha256(model));
//!
//! let result = ssa::solve(model, request.as_bytes(), &Control::default()).unwrap();
//! assert_eq!(result.termination, ssa::Termination::Completed);
//! let mean = result.statistics.expect("every replicate completed").mean.unwrap();
//! // mean[save time][species]; the exact mean at t = 10 is 20 (1 − e⁻²) ≈ 17.3.
//! assert!((mean[2][0] - 17.3).abs() < 3.0);
//! ```
//!
//! # Entry points
//!
//! | Task | API | Contracts (`crn-craft.` prefix omitted) |
//! |---|---|---|
//! | Read and validate inputs | [`parse_model`], [`parse_run`], [`check_bytes`], [`check_files`] | `model/0.1`, `run/0.1` → `check/0.1` |
//! | Analyze structure and conservation | [`analyze_bytes`], [`analyze_files`], [`CompiledNetwork`] | `state/0.1` → `analysis/0.2` |
//! | Direct SSA with sealed archives | [`archive::simulate_files`], [`archive::verify_archive`], [`simulation::simulate_bytes`] | `run/0.1` → `result/0.2`, `archive/0.1`, `archive-check/0.1` |
//! | Alternative exact SSA | [`ssa::prepare`], [`ssa::execute`], [`ssa::solve`], [`ssa::solve_files`] | `ssa-run/0.1`–`0.4` → `ssa-result/0.1`–`0.4` |
//! | Finite state projection | [`fsp::solve`] | `run/0.1` with `fsp` → `fsp-result/0.1` |
//! | Finite state projection of a closed rules model | [`fsp::v1::solve`], [`fsp::v1::solve_files`], [`rules::FlatView`] | `run/1.0` with `fsp` → `fsp-result/1.0` |
//! | Approximate tau leaping (experimental) | [`tau::solve`], [`tau::solve_files`] | `tau-request/0.1` → `tau-result/0.1` |
//! | Parameter sweeps and optimization | [`parameters::ParameterTemplate`], [`codegen::mass_action_source`], [`ssa::execute_with_kernel`] | `ssa-run` |
//! | Reference export | [`export::export_gillespy2`], [`export_gillespy2_file`] | `export/0.1` |
//! | Deterministic kinetics and Jacobians | [`kinetics::Kinetics`] | — |
//! | Discovery | [`describe`], [`schema`] | `describe/0.10` |
//!
//! Lower-level building blocks are [`strict_json`], [`validate`],
//! [`diagnostic`], [`sampling`], [`conservation`], [`cone`] and [`linkage`].
//!
//! # Errors, stops and cancellation
//!
//! Errors are [`diagnostic::Diagnostic`] values with a stable code (`CRN-E…`,
//! warnings `CRN-W…`) and a JSON Pointer to the offending input, listed in
//! [docs/diagnostics.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/diagnostics.md).
//! Report-style functions ([`check_bytes`], [`analyze_bytes`],
//! [`export::export_gillespy2`], [`archive::simulate_files`],
//! [`archive::verify_archive`]) never fail; they return a report with a
//! `valid` or `success` flag. In stochastic simulations, event, sampling,
//! population and wall-time limits and cancellation are typed outcomes in
//! each replicate's record, not errors; an exceeded output budget is the
//! run-level outcome `output_limit`; and ensemble statistics are withheld
//! unless every replicate completed. FSP has no replicates and fails with
//! `CRN-E070` instead. Cancel a running call through [`simulation::Control`].
//!
//! # Related crates and documentation
//!
//! - [`crn-craft-rules`](https://crates.io/crates/crn-craft-rules): rules
//!   engine for attributed species whose identities and numbers change.
//! - [`crn-craft-cli`](https://crates.io/crates/crn-craft-cli): the
//!   `crn-craft` command-line tool for both engines.
//! - [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/): Python
//!   bindings with prepared runs, importers and trajectory analysis.
//! - Repository documentation: the
//!   [tutorial](https://github.com/dietmarwo/CRN-craft/blob/main/docs/tutorial.md),
//!   [architecture](https://github.com/dietmarwo/CRN-craft/blob/main/docs/architecture.md),
//!   [input contracts](https://github.com/dietmarwo/CRN-craft/blob/main/docs/contracts.md),
//!   [Direct SSA](https://github.com/dietmarwo/CRN-craft/blob/main/docs/direct-ssa.md),
//!   [alternative SSA](https://github.com/dietmarwo/CRN-craft/blob/main/docs/alternative-ssa.md),
//!   [result and archive contract](https://github.com/dietmarwo/CRN-craft/blob/main/docs/archive-contract.md),
//!   [network API](https://github.com/dietmarwo/CRN-craft/blob/main/docs/network-api.md),
//!   [reusable parameters](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reusable-parameters.md)
//!   and the [technical report](https://github.com/dietmarwo/CRN-craft/blob/main/docs/technical-report.md)
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
pub mod analysis_v1;
pub mod archive;
pub mod codegen;
pub mod cone;
pub mod conservation;
pub mod contracts;
pub mod diagnostic;
pub mod export;
pub mod fsp;
pub mod kinetics;
pub mod linkage;
pub mod network;
#[doc(hidden)]
pub mod next_reaction;
pub mod parameters;
mod rates;
pub mod rules;
pub mod sampling;
pub mod simulation;
pub mod ssa;
pub mod strict_json;
pub mod tau;
pub mod validate;

pub use analysis::{AnalysisReport, analyze_bytes};
pub use network::CompiledNetwork;

use contracts::{Model, Run, RunResult};
use diagnostic::{Diagnostic, Severity};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::Path;

/// Result of `check`: whether a model, and optionally a run, passed parsing and validation.
#[derive(Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckReport {
    /// Always `crn-craft.check/0.1`.
    pub schema: String,
    /// True when no diagnostic has severity `error`.
    pub valid: bool,
    /// Always `structural_and_domain_validation`: no rate is evaluated and nothing runs.
    pub scope: String,
    /// Always false: a check report is not a simulation result.
    pub execution_available: bool,
    /// SHA-256 of the exact model bytes; null only when an input file could not be read.
    pub model_sha256: Option<String>,
    /// SHA-256 of the exact run bytes when a run was checked.
    pub run_sha256: Option<String>,
    /// Parse errors (the first per document) and every validation error and warning.
    pub diagnostics: Vec<Diagnostic>,
}

impl CheckReport {
    fn new() -> Self {
        Self {
            schema: "crn-craft.check/0.1".into(),
            valid: false,
            scope: "structural_and_domain_validation".into(),
            execution_available: false,
            model_sha256: None,
            run_sha256: None,
            diagnostics: Vec::new(),
        }
    }
    /// An invalid report carrying one diagnostic, such as a read error (`CRN-E030`).
    pub fn from_error(diagnostic: Diagnostic) -> Self {
        let mut report = Self::new();
        report.diagnostics.push(diagnostic);
        report
    }
    fn finish(&mut self) {
        self.valid = !self
            .diagnostics
            .iter()
            .any(|item| item.severity == Severity::Error);
    }
}

/// Lowercase hexadecimal SHA-256 of `bytes`: the digest every `sha256` field
/// and model reference uses. Hash the exact file bytes, never re-serialized
/// JSON.
pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// Strictly parses a `crn-craft.model/0.1` document.
///
/// This checks syntax and shape only, not semantics: use [`check_bytes`], or
/// [`CompiledNetwork::compile`], which always validates.
///
/// # Errors
///
/// The [`strict_json`] codes `CRN-E001` to `CRN-E004` and `CRN-E015`, and
/// `CRN-E032` at `/schema` for a document of the rules engine
/// (`crn-craft.model/1.0` or a `crn-craft.rules-*` contract), which the
/// network commands do not read.
///
/// ```
/// let duplicate = br#"{"schema": "crn-craft.model/0.1", "name": "m",
///   "convention": "stochastic_counts", "omega": 1.0, "parameters": {"k": 1.0, "k": 2.0},
///   "species": [{"name": "A", "initial_count": 1}], "reactions": []}"#;
/// let error = crn_craft::parse_model(duplicate).unwrap_err();
/// assert_eq!((error.code.as_str(), error.path.as_str()), ("CRN-E002", "/parameters/k"));
/// ```
pub fn parse_model(bytes: &[u8]) -> Result<Model, Diagnostic> {
    parse_network(bytes, "`crn-craft.model/0.1` network models")
}
/// Strictly parses a `crn-craft.run/0.1` document; like [`parse_model`], not
/// semantic validation.
///
/// # Errors
///
/// The [`strict_json`] codes `CRN-E001` to `CRN-E004` and `CRN-E015`, and
/// `CRN-E032` for a request of the rules engine (`crn-craft.run/1.0` or
/// `crn-craft.rules-run`).
pub fn parse_run(bytes: &[u8]) -> Result<Run, Diagnostic> {
    parse_network(bytes, "`crn-craft.run/0.1` requests")
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

/// Strictly parses a document of a network contract. A shape error
/// (`CRN-E003`) in a document whose `schema` names a rules-engine contract,
/// a `crn-craft.rules-*` or a 1.0 contract, becomes `CRN-E032` at `/schema`:
/// such a document never parses as a network contract, and the strict
/// parser would report the first field it does not know.
pub(crate) fn parse_network<T: serde::de::DeserializeOwned>(
    bytes: &[u8],
    reads: &str,
) -> Result<T, Diagnostic> {
    strict_json::parse(bytes).map_err(|error| {
        if error.code != "CRN-E003" {
            return error;
        }
        match schema_field(bytes) {
            Some(schema) if schema.starts_with("crn-craft.rules-") || schema.ends_with("/1.0") => {
                Diagnostic::error(
                    "CRN-E032",
                    "/schema",
                    format!("`{schema}` is a rules-engine contract; this command reads {reads}"),
                    "Check rules-engine models with `crn-craft rules-check` and run them with `crn-craft rules`; `crn-craft convert` writes a network model as a crn-craft.model/1.0 model.",
                )
            }
            _ => error,
        }
    })
}

/// Reads a model file and an optional `crn-craft.state/0.1` file, then runs
/// [`analyze_bytes`]. Each must be a regular file of at most 16 MiB; read
/// failures give a failed report with `CRN-E030` or `CRN-E004`.
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

/// Legacy M1a output and warning semantics, retained for existing consumers.
///
/// Like [`analyze_files`], but returns the `crn-craft.analysis/0.1` report
/// of [`analysis_v1::analyze_bytes`].
pub fn analyze_legacy_files(
    model_path: &Path,
    state_path: Option<&Path>,
) -> analysis_v1::AnalysisReport {
    let model = match read_document(model_path) {
        Ok(bytes) => bytes,
        Err(error) => return analysis_v1::AnalysisReport::failure(error),
    };
    let state = match state_path.map(read_document).transpose() {
        Ok(bytes) => bytes,
        Err(error) => return analysis_v1::AnalysisReport::failure(error),
    };
    analysis_v1::analyze_bytes(&model, state.as_deref())
}

/// Reads a model file and runs [`export::export_gillespy2`]; read failures
/// give a failed report with `CRN-E030` or `CRN-E004`.
pub fn export_gillespy2_file(model_path: &Path) -> export::ExportReport {
    match read_document(model_path) {
        Ok(bytes) => export::export_gillespy2(&bytes),
        Err(error) => export::ExportReport::failure(error),
    }
}

/// Validate in-memory input bytes. Only check_files verifies model.path identity.
///
/// Parses the model and, if given, the run, then applies
/// [`validate::validate_model`] and [`validate::validate_run`], which also
/// requires the run's `model.sha256` to match `model_bytes`. Never fails:
/// read the report's `valid` flag and diagnostics.
///
/// ```
/// let model = br#"{"schema": "crn-craft.model/0.1", "name": "typo",
///   "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
///   "species": [{"name": "A", "initial_count": 1}],
///   "reactions": [{"name": "convert", "reactants": {"A": 1}, "products": {"B": 1},
///                  "rate": {"kind": "mass_action", "coefficient": 1.0}}]}"#;
/// let report = crn_craft::check_bytes(model, None);
/// assert!(!report.valid);
/// let error = &report.diagnostics[0];
/// assert_eq!((error.code.as_str(), error.path.as_str()), ("CRN-E012", "/reactions/0/products/B"));
/// ```
pub fn check_bytes(model_bytes: &[u8], run_bytes: Option<&[u8]>) -> CheckReport {
    let mut report = CheckReport::new();
    report.model_sha256 = Some(sha256(model_bytes));
    report.run_sha256 = run_bytes.map(sha256);
    let model = parse_model(model_bytes);
    if let Ok(model) = &model {
        report.diagnostics.extend(validate::validate_model(model));
    }
    if let Some(bytes) = run_bytes {
        match parse_run(bytes) {
            Ok(run) => {
                if let Ok(model) = &model {
                    report.diagnostics.extend(validate::validate_run(
                        model,
                        &run,
                        report.model_sha256.as_deref().unwrap_or_default(),
                    ));
                }
            }
            Err(error) => report.diagnostics.push(error),
        }
    }
    if let Err(error) = model {
        report.diagnostics.insert(0, error);
    }
    report.finish();
    report
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

/// Resolve model.path relative to RUN and require it to denote the supplied MODEL.
/// A symlink is allowed only if its canonical target stays inside RUN's directory.
///
/// Reads both files (regular files of at most 16 MiB; otherwise `CRN-E030`
/// or `CRN-E004`), runs [`check_bytes`] and adds `CRN-E022` at `/model/path`
/// when the reference does not resolve to `model_path`. The check does not
/// guard against a filesystem that changes concurrently.
pub fn check_files(model_path: &Path, run_path: Option<&Path>) -> CheckReport {
    let model_bytes = match read_document(model_path) {
        Ok(x) => x,
        Err(e) => return CheckReport::from_error(e),
    };
    let run_bytes = match run_path.map(read_document).transpose() {
        Ok(x) => x,
        Err(e) => return CheckReport::from_error(e),
    };
    let mut report = check_bytes(&model_bytes, run_bytes.as_deref());
    if let (Some(path), Some(bytes)) = (run_path, run_bytes.as_deref()) {
        if let Ok(run) = parse_run(bytes) {
            if validate::valid_relative_path(&run.model.path) {
                let parent = path
                    .parent()
                    .filter(|x| !x.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                let resolved = (
                    parent.canonicalize(),
                    parent.join(&run.model.path).canonicalize(),
                    model_path.canonicalize(),
                );
                let matches = match resolved {
                    (Ok(base), Ok(referenced), Ok(supplied)) => {
                        referenced.starts_with(&base) && referenced == supplied
                    }
                    _ => false,
                };
                if !matches {
                    report.diagnostics.push(Diagnostic::error("CRN-E022", "/model/path", "model reference does not resolve to the supplied model inside the run directory", "Place the model below the run directory and pass that same file to check."));
                }
            }
        }
    }
    report.finish();
    report
}

/// The generated JSON Schema of one contract, as `crn-craft describe KIND`
/// prints it, or `None` for an unknown kind.
///
/// Kinds: `model`, `run`, `result`, `result-v1`, `archive`, `archive-check`,
/// `check`, `state`, `analysis`, `analysis-v1`, `export`, `tau-request`,
/// `tau-result`, `fsp-result`, `fsp-result-1.0`, `ssa-run` and
/// `ssa-result`. The `rules-*` kinds and the other 1.0 contracts
/// (`model-1.0`, `run-1.0`, ...) come from [`rules::schema`].
///
/// ```
/// let schema = crn_craft::schema("model").unwrap();
/// assert_eq!(schema["title"], "Model");
/// assert_eq!(schema["properties"]["schema"]["$ref"], "#/definitions/ModelVersion");
/// assert!(crn_craft::schema("rules-model").is_none());
/// ```
pub fn schema(kind: &str) -> Option<Value> {
    match kind {
        "model" => serde_json::to_value(schema_for!(Model)).ok(),
        "run" => serde_json::to_value(schema_for!(Run)).ok(),
        "result" => serde_json::to_value(schema_for!(simulation::SimulationResult)).ok(),
        "result-v1" => serde_json::to_value(schema_for!(RunResult)).ok(),
        "archive" => serde_json::to_value(schema_for!(archive::Manifest)).ok(),
        "archive-check" => serde_json::to_value(schema_for!(archive::ArchiveReport)).ok(),
        "check" => serde_json::to_value(schema_for!(CheckReport)).ok(),
        "state" => serde_json::to_value(schema_for!(analysis::NamedState)).ok(),
        "analysis" => serde_json::to_value(schema_for!(AnalysisReport)).ok(),
        "analysis-v1" => serde_json::to_value(schema_for!(analysis_v1::AnalysisReport)).ok(),
        "export" => serde_json::to_value(schema_for!(export::ExportReport)).ok(),
        "tau-request" => serde_json::to_value(schema_for!(tau::Request)).ok(),
        "tau-result" => serde_json::to_value(schema_for!(tau::TauResult)).ok(),
        "fsp-result" => serde_json::to_value(schema_for!(fsp::FspResult)).ok(),
        "fsp-result-1.0" => serde_json::to_value(schema_for!(fsp::v1::FspResult)).ok(),
        "ssa-run" => serde_json::to_value(schema_for!(ssa::Request)).ok(),
        "ssa-result" => serde_json::to_value(schema_for!(ssa::SsaResult)).ok(),
        _ => None,
    }
}

/// The `crn-craft.describe/0.10` discovery document.
///
/// It lists the version, commands, methods with their status, sampler
/// identities, the accepted `ssa-run` versions and their `fast` defaults,
/// rate laws, exporters, optional integrations, input and execution limits,
/// and the schemas of the eleven original contracts; [`schema`] returns the
/// others. The native CLI reports `crn-craft.describe/0.13`, which adds the
/// rules engine.
///
/// ```
/// let description = crn_craft::describe();
/// assert_eq!(description["schema"], "crn-craft.describe/0.10");
/// assert_eq!(description["version"], env!("CARGO_PKG_VERSION"));
/// assert!(description["commands"].as_array().unwrap().iter().any(|c| c == "ssa"));
/// assert_eq!(description["limits"]["species"], 4096);
/// ```
pub fn describe() -> Value {
    json!({
        "schema": "crn-craft.describe/0.10", "version": env!("CARGO_PKG_VERSION"), "milestone": "SSA-acceleration",
        "commands": ["describe", "check", "analyze", "analyze-v1", "export", "simulate", "verify", "simulate-cached", "fsp", "tau", "ssa"], "simulation_available": true,
        "aggregators": [{"name": "direct", "status": "implemented", "implementation_milestone": "M2"}],
        "methods": [{"name": "ssa", "status": "implemented"}, {"name": "fsp", "status": "implemented"}, {"name": "tau", "status": "experimental_approximation"}],
        "sampler": sampling::SAMPLER,
        "implementation_variants": ["direct", "dependency_cached_direct"],
        "alternative_exact_ssa": {"request":"crn-craft.ssa-run/0.3", "result":"crn-craft.ssa-result/0.3", "accepted_requests":["crn-craft.ssa-run/0.1","crn-craft.ssa-run/0.2","crn-craft.ssa-run/0.3","crn-craft.ssa-run/0.4"], "next_reaction_request":"crn-craft.ssa-run/0.4", "next_reaction_result":"crn-craft.ssa-result/0.4", "aggregators":["direct","cached_direct","direct_cr","sorting_direct","rssa_cr","next_reaction"], "waiting_times":["inverse","ziggurat"], "fast":{"default":true, "default_by_request":{"crn-craft.ssa-run/0.1":"unavailable","crn-craft.ssa-run/0.2":false,"crn-craft.ssa-run/0.3":true,"crn-craft.ssa-run/0.4":true}, "next_reaction_policy":"checked only; method.fast defaults to false and must not be true", "checked_policy":"method.fast=false: ChaCha8 words, lost-increment and 53-bit resolution checks, stop polls before every random word", "emulates":"JumpProcesses.jl 10.0.1: Xoshiro256++, running totals without loss or resolution checks, batched stop polls"}, "archive_format":null},
        "optional_integrations": {"python":"bindings/python", "jobs":"integrations/jobs", "ode_bdf_sensitivities":"integrations/diffsol"},
        "rate_laws": ["mass_action", "hill_activation", "hill_repression", "michaelis_menten"],
        "reference_exporters": [{"target":"gillespy2", "version": export::GILLESPY2_VERSION, "subset":"elementary_mass_action", "maximum_reactant_order":2}],
        "rate_evaluation_available": true, "observables": ["python_trajectory_autocorrelation_spectrum", "python_interval_censored_first_passage"],
        "model_cards": ["docs/model-cards/rate-laws.md", "docs/model-cards/direct.md"],
        "limits": {"document_bytes": strict_json::MAX_DOCUMENT_BYTES, "json_depth": strict_json::MAX_JSON_DEPTH,
            "species": validate::MAX_SPECIES, "parameters": validate::MAX_PARAMETERS,
            "reactions": validate::MAX_REACTIONS, "save_times": validate::MAX_SAVE_TIMES,
            "dependency_edges": network::MAX_DEPENDENCY_EDGES,
            "conservation_species": conservation::MAX_CONSERVATION_SPECIES,
            "conservation_reactions": conservation::MAX_CONSERVATION_REACTIONS,
            "cone_species": cone::MAX_CONE_SPECIES, "cone_reactions": cone::MAX_CONE_REACTIONS,
            "cone_rays": cone::MAX_CONE_RAYS, "cone_work": cone::MAX_CONE_WORK,
            "simulation_replicates": simulation::MAX_REPLICATES, "simulation_workers": simulation::MAX_WORKERS,
            "simulation_workspace_bytes": simulation::MAX_WORKSPACE_BYTES, "archive_bytes": simulation::MAX_RESULT_BYTES},
        "schemas": {"model": schema("model"), "run": schema("run"), "result": schema("result"), "result-v1":schema("result-v1"), "archive":schema("archive"), "archive-check":schema("archive-check"), "check": schema("check"), "state": schema("state"), "analysis": schema("analysis"), "analysis-v1": schema("analysis-v1"), "export": schema("export")}
    })
}
