#![forbid(unsafe_code)]
//! CRN-craft M1a: strict inputs, compiled networks, propensities and analysis.
//! There is deliberately no time integration or simulation API yet.
pub mod analysis;
pub mod conservation;
pub mod contracts;
pub mod diagnostic;
pub mod network;
mod rates;
pub mod strict_json;
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

#[derive(Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckReport {
    pub schema: String,
    pub valid: bool,
    pub scope: String,
    pub execution_available: bool,
    pub model_sha256: Option<String>,
    pub run_sha256: Option<String>,
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

pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

pub fn parse_model(bytes: &[u8]) -> Result<Model, Diagnostic> {
    strict_json::parse(bytes)
}
pub fn parse_run(bytes: &[u8]) -> Result<Run, Diagnostic> {
    strict_json::parse(bytes)
}

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

/// Validate in-memory input bytes. Only check_files verifies model.path identity.
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

pub fn schema(kind: &str) -> Option<Value> {
    match kind {
        "model" => serde_json::to_value(schema_for!(Model)).ok(),
        "run" => serde_json::to_value(schema_for!(Run)).ok(),
        "result" => serde_json::to_value(schema_for!(RunResult)).ok(),
        "check" => serde_json::to_value(schema_for!(CheckReport)).ok(),
        "state" => serde_json::to_value(schema_for!(analysis::NamedState)).ok(),
        "analysis" => serde_json::to_value(schema_for!(AnalysisReport)).ok(),
        _ => None,
    }
}

pub fn describe() -> Value {
    json!({
        "schema": "crn-craft.describe/0.2", "version": env!("CARGO_PKG_VERSION"), "milestone": "M1a",
        "commands": ["describe", "check", "analyze"], "simulation_available": false,
        "aggregators": [{"name": "direct", "status": "contract_only", "implementation_milestone": "M2"}],
        "methods": [{"name": "ssa", "status": "contract_only"}, {"name": "fsp", "status": "contract_only"}],
        "rate_laws": ["mass_action", "hill_activation", "hill_repression", "michaelis_menten"],
        "rate_evaluation_available": true, "observables": [],
        "model_cards": ["docs/model-cards/rate-laws.md", "docs/model-cards/direct.md"],
        "limits": {"document_bytes": strict_json::MAX_DOCUMENT_BYTES, "json_depth": strict_json::MAX_JSON_DEPTH,
            "species": validate::MAX_SPECIES, "parameters": validate::MAX_PARAMETERS,
            "reactions": validate::MAX_REACTIONS, "save_times": validate::MAX_SAVE_TIMES,
            "dependency_edges": network::MAX_DEPENDENCY_EDGES,
            "conservation_species": conservation::MAX_CONSERVATION_SPECIES,
            "conservation_reactions": conservation::MAX_CONSERVATION_REACTIONS},
        "schemas": {"model": schema("model"), "run": schema("run"), "result": schema("result"), "check": schema("check"), "state": schema("state"), "analysis": schema("analysis")}
    })
}
