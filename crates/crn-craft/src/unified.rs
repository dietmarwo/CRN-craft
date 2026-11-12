//! The unified command set of release 0.7.0: discovery and checks.
//!
//! `crn-craft describe` prints [`describe`] (`crn-craft.describe/1.0`) and
//! `crn-craft check MODEL [REQUEST]` (alias `rules-check`) prints
//! [`check_any_files`]. Python's `crn_craft.describe()` and
//! `crn_craft.check()` return the same documents.

use crate::diagnostic::Diagnostic;
use crate::rules::{self, CheckVersion, RulesCheck};
use crate::sha256;
use serde_json::{Value, json};
use std::path::Path;

/// The milestone `crn-craft --version` and [`describe`] report.
pub const MILESTONE: &str = "unified-engine";

/// The commands of the unified command line, in `--help` order.
pub const COMMANDS: [&str; 7] = [
    "describe", "check", "run", "verify", "analyze", "convert", "export",
];

/// The `describe KIND` kinds of the 1.0 contracts and the current reports.
pub const SCHEMAS: [&str; 10] = [
    "model-1.0",
    "run-1.0",
    "result-1.0",
    "check-1.0",
    "fsp-result-1.0",
    "analysis-1.0",
    "archive-1.0",
    "archive-check-1.0",
    "export",
    "state",
];

/// The `describe KIND` kinds of older contracts that are still read.
pub const LEGACY_SCHEMAS: [&str; 5] = [
    "rules-model",
    "rules-run",
    "rules-result",
    "rules-check",
    "model",
];

/// The discovery document `crn-craft.describe/1.0`: what `crn-craft
/// describe` prints and Python's `crn_craft.describe()` returns.
///
/// It lists the commands, their aliases, the 1.0 contracts and the schema
/// kinds of `describe KIND` first: methods, Static Direct selectors,
/// quantiles, finite state projection, exporters and limits. The older
/// contracts that are still read (rules models and requests, the network
/// `crn-craft.model/0.1` through the converter) are under `legacy`, with the
/// rules engine's 0.6.0 discovery entry as `legacy.rules`.
///
/// ```
/// let description = crn_craft::describe();
/// assert_eq!(description["schema"], "crn-craft.describe/1.0");
/// assert_eq!(description["crate"], "crn-craft");
/// assert_eq!(description["contracts"]["model"], "crn-craft.model/1.0");
/// assert_eq!(description["aliases"]["rules"], "run");
/// assert!(!description["commands"].as_array().unwrap().iter().any(|c| c == "simulate"));
/// ```
#[must_use]
pub fn describe() -> Value {
    let entry = rules::describe();
    let unified = &entry["contracts_1_0"];
    let mut methods = entry["methods"].clone();
    if let Some(list) = methods.as_array_mut() {
        list.push("fsp".into());
    }
    json!({
        "schema": "crn-craft.describe/1.0",
        "crate": "crn-craft",
        "version": env!("CARGO_PKG_VERSION"),
        "milestone": MILESTONE,
        "commands": COMMANDS,
        "aliases": {"rules": "run", "rules-check": "check"},
        "python": "crn_craft",
        "contracts": {
            "model": "crn-craft.model/1.0",
            "request": "crn-craft.run/1.0",
            "result": "crn-craft.result/1.0",
            "check": "crn-craft.check/1.0",
            "fsp_result": "crn-craft.fsp-result/1.0",
            "analysis": "crn-craft.analysis/1.0",
            "archive": "crn-craft.archive/1.0",
            "archive_check": "crn-craft.archive-check/1.0",
            "export": "crn-craft.export/0.1",
            "state": "crn-craft.state/0.1",
            "import": "crn-craft.import/1.0",
            "report": "crn-craft.report/1.0"
        },
        "schemas": SCHEMAS,
        "methods": methods,
        "approximate_methods": entry["approximate_methods"],
        "default_method": "tree",
        "selectors": unified["selectors"],
        "quantiles": {
            "request": "outputs.quantiles",
            "result": "statistics.quantiles",
            "max": rules::MAX_QUANTILES,
            "definition": "numpy-linear-type-7",
            "python_default": [0.025, 0.5, 0.975]
        },
        "fsp": {
            "request": "crn-craft.run/1.0 with method fsp",
            "result": "crn-craft.fsp-result/1.0",
            "max_states": crate::fsp::MAX_STATES,
            "max_edges": crate::fsp::MAX_EDGES,
            "max_work": crate::fsp::MAX_WORK
        },
        "archives": {
            "command": "run MODEL REQUEST --archive NEW_DIRECTORY",
            "verify": "verify DIRECTORY",
            "manifest": "crn-craft.archive/1.0",
            "report": "crn-craft.archive-check/1.0"
        },
        "reference_exporters": [{"target": "gillespy2", "version": crate::export::GILLESPY2_VERSION, "subset": "mass_action_and_custom_propensities", "maximum_mass_action_order": 2, "maximum_custom_order": crate::export::MAX_CUSTOM_ORDER}],
        "limits": {
            "document_bytes": crate::strict_json::MAX_DOCUMENT_BYTES,
            "json_depth": crate::strict_json::MAX_JSON_DEPTH,
            "types": unified["limits"]["types"],
            "parameters": unified["limits"]["parameters"],
            "observables": unified["limits"]["observables"],
            "rules": unified["limits"]["rules"],
            "attributes": entry["limits"]["attributes"],
            "products": entry["limits"]["products"],
            "pattern_count": rules::model::MAX_PATTERN_COUNT,
            "initial_entries": entry["limits"]["initial_entries"],
            "expression_bytes": entry["limits"]["expression_bytes"],
            "replicates": rules::MAX_REPLICATES,
            "report_times": rules::MAX_REPORT_TIMES,
            "threads": rules::MAX_WORKERS,
            "result_bytes": rules::MAX_RESULT_BYTES,
            "workspace_bytes": rules::MAX_WORKSPACE_BYTES
        },
        "sampler": entry["sampler"],
        "reference": "docs/rules-language.md",
        "legacy": {
            "models": [
                "crn-craft.rules-model/0.1", "crn-craft.rules-model/0.2",
                "crn-craft.rules-model/0.3", "crn-craft.rules-model/0.4",
                "crn-craft.model/0.1"
            ],
            "network_model": "crn-craft.model/0.1 is read through the converter by check, analyze, convert and export; run reads rules and 1.0 models only",
            "requests": entry["accepted_requests"],
            "results": ["crn-craft.rules-result/0.1", "crn-craft.rules-result/0.2", "crn-craft.rules-result/0.3"],
            "checks": ["crn-craft.rules-check/0.1"],
            "schemas": LEGACY_SCHEMAS,
            "rules": entry
        }
    })
}

fn report(model: &[u8], request: Option<&[u8]>, diagnostics: Vec<Diagnostic>) -> RulesCheck {
    RulesCheck {
        schema: CheckVersion::V10,
        valid: diagnostics.is_empty(),
        model_sha256: Some(sha256(model)),
        run_sha256: request.map(sha256),
        diagnostics,
    }
}

/// Checks a model, and optionally a request for it, without running
/// anything: `crn-craft check MODEL [REQUEST]` on bytes.
///
/// - A rules-model alone: [`rules::check_model`] (`rules-check/0.1`, as in
///   release 0.6.0; `check/1.0` for a `crn-craft.model/1.0` document).
/// - A network `crn-craft.model/0.1` alone: `check/1.0` of the conversion,
///   valid when [`rules::convert::from_network`] converts it (its
///   validation errors and `CRN-R016` refusals otherwise).
/// - Any model with a request: `check/1.0` with `run_sha256`, the model's
///   check and then the errors `crn-craft run` reports before running
///   ([`rules::prepare`], or [`crate::fsp::v1::check`] for an `fsp`
///   request). A network model with a request gets `CRN-E033` (convert it
///   first): requests never name a network model.
///
/// ```
/// let model = br#"{"schema": "crn-craft.rules-model/0.1", "name": "decay",
///   "parameters": {"k": 0.5}, "types": [{"name": "A"}],
///   "initial": [{"type": "A", "count": 10}], "observables": {"A": {"type": "A"}},
///   "rules": [{"name": "decay", "reactants": [{"label": "a", "type": "A"}],
///              "rate": "k", "products": []}]}"#;
/// let report = serde_json::to_value(crn_craft::check_any(model, None)).unwrap();
/// assert_eq!(report["schema"], "crn-craft.rules-check/0.1");
/// assert_eq!(report["valid"], true);
/// ```
#[must_use]
pub fn check_any(model: &[u8], request: Option<&[u8]>) -> RulesCheck {
    let network = crate::schema_field(model).as_deref() == Some(rules::convert::NETWORK_SCHEMA);
    match request {
        None if network => {
            let diagnostics = match rules::convert::from_network(model) {
                Ok(converted) => rules::check_model(&converted).diagnostics,
                Err(diagnostics) => diagnostics,
            };
            report(model, None, diagnostics)
        }
        None => rules::check_model(model),
        Some(run) => {
            let checked = rules::check_model(model);
            let diagnostics = if !checked.valid {
                checked.diagnostics
            } else if crate::fsp::v1::is_fsp_request(run) {
                crate::fsp::v1::check(model, run)
            } else {
                rules::prepare(model, run).err().unwrap_or_default()
            };
            report(model, Some(run), diagnostics)
        }
    }
}

/// [`check_any`] for files: `crn-craft check MODEL [REQUEST]`.
///
/// Each file must be a regular file of at most 16 MiB (`CRN-E030` or
/// `CRN-E004` otherwise, in a `check/1.0` report). With a request, its
/// `model.path` must name `model_path` below the request's directory, as
/// `crn-craft run` requires (`CRN-R020` at `/model/path`).
#[must_use]
pub fn check_any_files(model_path: &Path, request_path: Option<&Path>) -> RulesCheck {
    let failure = |model: Option<&[u8]>, error| RulesCheck {
        schema: CheckVersion::V10,
        valid: false,
        model_sha256: model.map(sha256),
        run_sha256: None,
        diagnostics: vec![error],
    };
    let model = match rules::read(model_path) {
        Ok(bytes) => bytes,
        Err(error) => return failure(None, error),
    };
    let Some(request_path) = request_path else {
        return check_any(&model, None);
    };
    let request = match rules::read(request_path) {
        Ok(bytes) => bytes,
        Err(error) => return failure(Some(&model), error),
    };
    let mut report = check_any(&model, Some(&request));
    let reference = serde_json::from_slice::<Value>(&request)
        .ok()
        .and_then(|v| {
            v.pointer("/model/path")
                .and_then(Value::as_str)
                .map(String::from)
        });
    if let Some(relative) = reference {
        if !rules::references_model(model_path, request_path, &relative) {
            report.diagnostics.push(Diagnostic::error(
                "CRN-R020",
                "/model/path",
                "request does not identify the supplied model below its directory",
                "Use the referenced model within the request directory.",
            ));
            report.valid = false;
        }
    }
    report
}
