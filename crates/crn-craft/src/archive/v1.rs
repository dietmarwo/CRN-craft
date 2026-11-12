//! Sealed archives of rules runs: `crn-craft.archive/1.0`.
//!
//! [`run_files`] (the `crn-craft rules MODEL REQUEST --archive NEW_DIRECTORY`
//! command) and [`write()`] (the byte API behind Python's
//! `crn_craft.rules.archive`) run a rules request, a simulation request of
//! any rules-run version or a `crn-craft.run/1.0` request whose method is
//! `fsp`, and write a new directory with five files:
//!
//! | File | Content |
//! | --- | --- |
//! | `model.json` | the exact model bytes |
//! | `request.json` | the exact request bytes, with its original `model.path` |
//! | `result.json` | the result as `crn-craft rules MODEL REQUEST` prints it: compact JSON with sorted keys and a trailing newline |
//! | `manifest.json` | `crn-craft.archive/1.0`: length and SHA-256 of the three payloads |
//! | `COMPLETE` | the SHA-256 of the manifest bytes and a newline, written last |
//!
//! [`verify`] checks a finished directory without simulating: exactly these
//! regular files, the marker against the manifest, lengths and digests, that
//! the inputs still prepare, that the result parses strictly and is the
//! engine's own serialization, the hashes inside the result and the
//! relations between the request and the result (identity, replicate
//! accounting, terminations, kept paths, statistics). Both return a
//! `crn-craft.archive-check/1.0` [`ArchiveReport`]; archive failures are
//! `CRN-E065`. The 0.1 archives of the network engine stay readable through
//! [`super::verify_any`], which `crn-craft verify DIRECTORY` calls.
//!
//! An archive is the sealed, verifiable record of one run for replay. A
//! complete archive can hold an incomplete simulation: read `archive_complete`
//! and `simulation_complete` together. Integrity checks are not signatures.
//!
//! ```
//! use crn_craft::{archive::v1, rules::Control, sha256};
//!
//! let dir = std::env::temp_dir().join(format!("crn-craft-doc-archive-v1-{}", std::process::id()));
//! let _ = std::fs::remove_dir_all(&dir);
//! std::fs::create_dir_all(&dir).unwrap();
//! let model = br#"{"schema": "crn-craft.model/1.0", "name": "decay",
//!   "parameters": {"k": 1.0}, "types": [{"name": "A"}],
//!   "initial": [{"type": "A", "count": 5}], "observables": {"A": {"type": "A"}},
//!   "rules": [{"name": "death", "reactants": [{"label": "a", "type": "A"}],
//!              "rate": "k", "products": []}]}"#;
//! let request = format!(r#"{{"schema": "crn-craft.run/1.0",
//!   "model": {{"path": "model.json", "sha256": "{}"}},
//!   "method": {{"kind": "static_direct"}},
//!   "time": {{"end": 1.0, "report_times": [0.0, 1.0]}},
//!   "ensemble": {{"replicates": 4, "root_seed": 3}},
//!   "outputs": {{"statistics": ["mean"], "keep_paths": 4, "keep_instances": 0}},
//!   "budget": {{"max_events_per_replicate": 100, "max_sampling_steps_per_replicate": 10000,
//!              "max_instances": 10, "threads": 1, "max_output_bytes": 100000}}}}"#,
//!   sha256(model));
//! std::fs::write(dir.join("model.json"), model).unwrap();
//! std::fs::write(dir.join("request.json"), &request).unwrap();
//!
//! let out = dir.join("archive");
//! let written = v1::run_files(&dir.join("model.json"), &dir.join("request.json"), &out,
//!     &Control::default());
//! assert!(written.valid && written.archive_complete && written.simulation_complete,
//!     "{:?}", written.diagnostics);
//! let checked = v1::verify(&out);
//! assert!(checked.valid && checked.simulation_complete, "{:?}", checked.diagnostics);
//! assert_eq!(checked.model_sha256.as_deref(), Some(sha256(model).as_str()));
//! std::fs::remove_dir_all(&dir).unwrap();
//! ```
use super::{encoded, error, io_error, read_regular, write_new};
use crate::{
    contracts::Statistic,
    diagnostic::Diagnostic,
    fsp::v1 as fsp,
    rules::{self, Control, FlatView, RulesResult, RunRequest, Termination},
    sha256, strict_json,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    path::Path,
};

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// The three payload files, in manifest order.
const PAYLOAD_NAMES: [&str; 3] = ["model.json", "request.json", "result.json"];
/// Every file of a finished archive.
const ALL_NAMES: [&str; 5] = [
    "model.json",
    "request.json",
    "result.json",
    "manifest.json",
    "COMPLETE",
];
/// Largest manifest file.
const MANIFEST_BYTES: u64 = 8192;
/// Largest completion marker: 64 hexadecimal digits and a newline.
const MARKER_BYTES: u64 = 65;
/// Largest input document, as for every contract.
const INPUT_BYTES: u64 = strict_json::MAX_DOCUMENT_BYTES as u64;
/// The manifest version string, also read leniently by [`is_archive`].
const ARCHIVE_1_0: &str = "crn-craft.archive/1.0";
/// Largest deterministic work of a finite state projection result.
const MAX_FSP_WORK: u64 = 200_000_000;

/// Schema tag of an archive manifest, `crn-craft.archive/1.0`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ArchiveVersion {
    #[doc = concat!("`crn-craft.archive/1.0`.")]
    #[serde(rename = "crn-craft.archive/1.0")]
    V10,
}

/// Length and digest of one payload file.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ArchiveFile {
    /// File name inside the archive directory: `model.json`, `request.json` or
    /// `result.json`.
    pub path: String,
    /// Exact file length in bytes.
    pub bytes: u64,
    /// SHA-256 of the file bytes, lowercase hexadecimal.
    pub sha256: String,
}

/// `manifest.json` of a `crn-craft.archive/1.0` archive: the three payload
/// files, each listed exactly once.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    #[doc = concat!("Always `crn-craft.archive/1.0`.")]
    pub schema: ArchiveVersion,
    /// The three payload files, written in the order model, request, result.
    pub files: Vec<ArchiveFile>,
}

/// Schema tag of an archive report, `crn-craft.archive-check/1.0`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum CheckVersion {
    #[doc = concat!("`crn-craft.archive-check/1.0`.")]
    #[serde(rename = "crn-craft.archive-check/1.0")]
    V10,
}

/// Report of writing or verifying a `crn-craft.archive/1.0` archive: archive
/// integrity and the outcome of the archived run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ArchiveReport {
    #[doc = concat!("Always `crn-craft.archive-check/1.0`.")]
    pub schema: CheckVersion,
    /// True when the archive was written, or verified, without error.
    pub valid: bool,
    /// True when every file, the manifest and the completion marker are
    /// present and consistent. It says nothing about whether the run
    /// finished.
    pub archive_complete: bool,
    /// True when the archived result is complete: a simulation result whose
    /// termination is `completed`, or a finite state projection.
    pub simulation_complete: bool,
    /// The archived result's `schema`, for example `crn-craft.result/1.0` or
    /// `crn-craft.fsp-result/1.0`; null when the archive is invalid.
    pub result: Option<String>,
    /// The archived simulation result's ensemble termination; null for a
    /// finite state projection and when the archive is invalid.
    pub termination: Option<Termination>,
    /// SHA-256 of the archived model bytes; null when the archive is invalid.
    pub model_sha256: Option<String>,
    /// SHA-256 of the archived request bytes; null when the archive is invalid.
    pub run_sha256: Option<String>,
    /// Why the archive is invalid; empty for a valid archive.
    pub diagnostics: Vec<Diagnostic>,
}

impl ArchiveReport {
    /// An invalid report carrying `diagnostics`.
    #[must_use]
    pub fn failure(diagnostics: Vec<Diagnostic>) -> Self {
        Self {
            schema: CheckVersion::V10,
            valid: false,
            archive_complete: false,
            simulation_complete: false,
            result: None,
            termination: None,
            model_sha256: None,
            run_sha256: None,
            diagnostics,
        }
    }
}

/// Runs a request file against the model it references and archives the
/// run in `output`, a new directory whose parent exists: the
/// `crn-craft rules MODEL REQUEST --archive NEW_DIRECTORY` command.
///
/// The request's `model.path`, resolved from the request's directory, must
/// name `model_path` inside that directory, as for `crn-craft rules`
/// (`CRN-R020` at `/model/path` otherwise). Then as [`write()`].
///
/// Never panics on bad input; failures are an invalid report.
pub fn run_files(
    model_path: &Path,
    request_path: &Path,
    output: &Path,
    control: &Control,
) -> ArchiveReport {
    let action = || -> Result<ArchiveReport, Vec<Diagnostic>> {
        let model = rules::read(model_path).map_err(|e| vec![e])?;
        let request = rules::read(request_path).map_err(|e| vec![e])?;
        let relative = if fsp::is_fsp_request(&request) {
            strict_json::parse::<fsp::FspRequest>(&request)
                .map_err(|e| vec![e])?
                .model
                .path
        } else {
            strict_json::parse::<RunRequest>(&request)
                .map_err(|e| vec![e])?
                .model
                .path
        };
        if !rules::references_model(model_path, request_path, &relative) {
            return Err(vec![Diagnostic::error(
                "CRN-R020",
                "/model/path",
                "rules request does not identify the supplied model below its directory",
                "Use the referenced model within the request directory.",
            )]);
        }
        archive(&model, &request, output, control)
    };
    action().unwrap_or_else(ArchiveReport::failure)
}

/// Runs `request` against `model` (exact bytes; the request's `model.path`
/// is not resolved) and archives the run in `output`, a new directory whose
/// parent exists.
///
/// A simulation request is prepared first: an input error (strict parsing
/// `CRN-E001`–`CRN-E004`, `CRN-R…` request and model errors) creates no
/// directory. Then the directory is created (`CRN-E065` if it exists), the
/// model and request bytes are written, the run executes and the result,
/// the manifest and last the completion marker are written. A finite state
/// projection is solved before the directory is created, so a failed solve
/// (`CRN-E070`) leaves nothing. A run that fails even after trimming
/// (`CRN-R045`) or a write failure leaves a directory without `COMPLETE`,
/// which is not an archive. Each file is created exclusively and
/// synchronized; on Unix the directory is synchronized after the marker.
///
/// Never panics on bad input; failures are an invalid report.
pub fn write(model: &[u8], request: &[u8], output: &Path, control: &Control) -> ArchiveReport {
    archive(model, request, output, control).unwrap_or_else(ArchiveReport::failure)
}

fn archive(
    model: &[u8],
    request: &[u8],
    output: &Path,
    control: &Control,
) -> Result<ArchiveReport, Vec<Diagnostic>> {
    let one = |e: Diagnostic| vec![e];
    if fsp::is_fsp_request(request) {
        let result = fsp::solve(model, request, control)?;
        let run: fsp::FspRequest = strict_json::parse(request).map_err(one)?;
        let bytes = printed(&result, output_limit(run.budget.max_output_bytes)).map_err(one)?;
        fs::create_dir(output).map_err(|e| one(io_error(e)))?;
        write_new(&output.join("model.json"), model).map_err(one)?;
        write_new(&output.join("request.json"), request).map_err(one)?;
        write_new(&output.join("result.json"), &bytes).map_err(one)?;
        seal(output, [model, request, &bytes]).map_err(one)?;
        return Ok(success(&result.schema, None, model, request));
    }
    let prepared = rules::prepare(model, request)?;
    let run: RunRequest = strict_json::parse(request).map_err(one)?;
    fs::create_dir(output).map_err(|e| one(io_error(e)))?;
    write_new(&output.join("model.json"), model).map_err(one)?;
    write_new(&output.join("request.json"), request).map_err(one)?;
    let result = rules::execute(&prepared, control).map_err(one)?;
    let bytes = printed(&result, output_limit(run.budget.max_output_bytes)).map_err(one)?;
    write_new(&output.join("result.json"), &bytes).map_err(one)?;
    seal(output, [model, request, &bytes]).map_err(one)?;
    Ok(success(
        &result.schema,
        Some(result.termination),
        model,
        request,
    ))
}

/// The size limit of a result: the request's budget, capped at 1 GiB.
fn output_limit(max_output_bytes: u64) -> u64 {
    max_output_bytes.min(rules::MAX_RESULT_BYTES)
}

/// A result as `crn-craft rules` prints it: compact JSON of its value (keys
/// sorted) and a newline, within `limit` bytes.
fn printed<T: Serialize>(result: &T, limit: u64) -> Result<Vec<u8>, Diagnostic> {
    let value = serde_json::to_value(result).map_err(|e| error(e.to_string()))?;
    encoded(&value, limit)
}

/// Writes the manifest of the three payloads, then the completion marker.
fn seal(output: &Path, payloads: [&[u8]; 3]) -> Result<(), Diagnostic> {
    let files = PAYLOAD_NAMES
        .iter()
        .zip(payloads)
        .map(|(name, bytes)| ArchiveFile {
            path: (*name).into(),
            bytes: bytes.len() as u64,
            sha256: sha256(bytes),
        })
        .collect();
    let manifest = encoded(
        &Manifest {
            schema: ArchiveVersion::V10,
            files,
        },
        MANIFEST_BYTES,
    )?;
    let marker = format!("{}\n", sha256(&manifest));
    write_new(&output.join("manifest.json"), &manifest)?;
    write_new(&output.join("COMPLETE"), marker.as_bytes())?;
    #[cfg(unix)]
    File::open(output)
        .and_then(|f| f.sync_all())
        .map_err(io_error)?;
    Ok(())
}

/// The report of a valid archive of `result_schema`.
fn success(
    result_schema: &impl Serialize,
    termination: Option<Termination>,
    model: &[u8],
    request: &[u8],
) -> ArchiveReport {
    ArchiveReport {
        schema: CheckVersion::V10,
        valid: true,
        archive_complete: true,
        simulation_complete: termination.is_none_or(|t| t == Termination::Completed),
        result: serde_json::to_value(result_schema)
            .ok()
            .and_then(|v| v.as_str().map(String::from)),
        termination,
        model_sha256: Some(sha256(model)),
        run_sha256: Some(sha256(request)),
        diagnostics: Vec::new(),
    }
}

/// Whether `path` holds a `crn-craft.archive/1.0` archive rather than a 0.1
/// archive of the network engine: its `manifest.json`, read leniently, names
/// `crn-craft.archive/1.0`, or, without a readable manifest version, it has
/// a `request.json` and no `run.json`.
#[must_use]
pub fn is_archive(path: &Path) -> bool {
    let version = read_regular(&path.join("manifest.json"), MANIFEST_BYTES)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<Value>(&bytes).ok())
        .and_then(|value| {
            value
                .get("schema")
                .and_then(Value::as_str)
                .map(String::from)
        });
    match version {
        Some(version) => version == ARCHIVE_1_0,
        None => {
            fs::symlink_metadata(path.join("request.json")).is_ok()
                && fs::symlink_metadata(path.join("run.json")).is_err()
        }
    }
}

/// Integrity and relational checks of a `crn-craft.archive/1.0` directory,
/// not cryptographic authenticity or proof of stochastic correctness. It
/// never runs a simulation, never executes supplied code and never follows
/// manifest paths.
///
/// It requires exactly the five regular files (no symlinks, directories or
/// extra files) within their size limits, a marker matching the manifest,
/// a strict manifest naming the three payloads once each, matching lengths
/// and digests, inputs that still prepare, a result within the request's
/// output budget that parses strictly and is byte for byte the engine's
/// serialization of itself, and a result consistent with its inputs:
///
/// - identity: result version, both hashes, sampler, seed, types, rules,
///   observables, report times, replicate count, method and continuous
///   columns as the engine reports them for these inputs (taken from a
///   cancelled execution, which simulates nothing);
/// - every replicate record: id, times, budgets, rule events summing to the
///   events, unstarted and complete records, diagnostics of halted ones;
/// - the aggregate counts, the ensemble termination and the ensemble
///   diagnostics;
/// - kept paths (which replicates, their lengths and, when the first report
///   time is 0, their initial row), final instances and observations;
/// - statistics: present exactly for a complete ensemble that asked for
///   them, the requested estimators and quantile probabilities, dimensions,
///   finite values, standard errors equal to the square root of variance
///   over samples, quantiles nondecreasing in the probability, and, when
///   every path is kept, the mean, variance, standard error and quantiles
///   recomputed from the paths, bit for bit.
///
/// A finite state projection result is checked for its hashes, the flat
/// view's species labels and initial state, the bounds, the report times
/// and the dimensions of every snapshot.
///
/// A valid archive of an incomplete simulation reports `valid: true` and
/// `simulation_complete: false`. Failures are `CRN-E065` with the first
/// inconsistency, or the strict parser's codes.
#[must_use]
pub fn verify(path: &Path) -> ArchiveReport {
    let action = || -> Result<ArchiveReport, Diagnostic> {
        if !fs::symlink_metadata(path).map_err(io_error)?.is_dir() {
            return Err(error("archive root is not a real directory"));
        }
        let mut names = BTreeSet::new();
        let mut oversized = false;
        for entry in fs::read_dir(path).map_err(io_error)? {
            let entry = entry.map_err(io_error)?;
            let name = entry
                .file_name()
                .into_string()
                .map_err(|_| error("non-UTF8 archive name"))?;
            if !ALL_NAMES.contains(&name.as_str())
                || !entry.file_type().map_err(io_error)?.is_file()
            {
                return Err(error("unexpected file, directory or symlink in archive"));
            }
            oversized |= entry.metadata().map_err(io_error)?.len() > limit_of(&name);
            names.insert(name);
        }
        if names.len() != ALL_NAMES.len() || oversized {
            return Err(error("incomplete or oversized archive"));
        }
        let manifest_bytes = read_regular(&path.join("manifest.json"), MANIFEST_BYTES)?;
        let marker = read_regular(&path.join("COMPLETE"), MARKER_BYTES)?;
        if marker != format!("{}\n", sha256(&manifest_bytes)).as_bytes() {
            return Err(error("completion marker does not match manifest"));
        }
        let manifest: Manifest = strict_json::parse(&manifest_bytes)?;
        let paths: BTreeSet<_> = manifest.files.iter().map(|f| f.path.as_str()).collect();
        if manifest.files.len() != PAYLOAD_NAMES.len()
            || paths != PAYLOAD_NAMES.into_iter().collect()
        {
            return Err(error(
                "manifest paths are not exactly the three unique payload names",
            ));
        }
        let mut payloads = BTreeMap::new();
        for entry in manifest.files {
            let bytes = read_regular(&path.join(&entry.path), limit_of(&entry.path))?;
            if entry.bytes != bytes.len() as u64 || entry.sha256 != sha256(&bytes) {
                return Err(error("payload checksum or length mismatch"));
            }
            payloads.insert(entry.path, bytes);
        }
        let (model, request, result) = (
            &payloads["model.json"],
            &payloads["request.json"],
            &payloads["result.json"],
        );
        if fsp::is_fsp_request(request) {
            verify_fsp(model, request, result)
        } else {
            verify_simulation(model, request, result)
        }
    };
    action().unwrap_or_else(|e| ArchiveReport::failure(vec![e]))
}

/// The size limit of one archive file.
fn limit_of(name: &str) -> u64 {
    match name {
        "result.json" => rules::MAX_RESULT_BYTES,
        "manifest.json" => MANIFEST_BYTES,
        "COMPLETE" => MARKER_BYTES,
        _ => INPUT_BYTES,
    }
}

fn invalid_inputs() -> Diagnostic {
    error("invalid archived inputs")
}

/// Parses an archived result and requires its exact bytes to be the
/// engine's serialization of it, within the request's output budget.
fn parse_result<T: Serialize + serde::de::DeserializeOwned>(
    bytes: &[u8],
    max_output_bytes: u64,
) -> Result<T, Diagnostic> {
    if bytes.len() as u64 > output_limit(max_output_bytes) {
        return Err(error("result exceeds the request's output budget"));
    }
    let result: T = strict_json::parse_bounded(bytes, rules::MAX_RESULT_BYTES as usize)?;
    if printed(&result, u64::MAX)? != bytes {
        return Err(error("result bytes are not the engine's serialization"));
    }
    Ok(result)
}

fn verify_simulation(
    model: &[u8],
    request: &[u8],
    bytes: &[u8],
) -> Result<ArchiveReport, Diagnostic> {
    let prepared = rules::prepare(model, request).map_err(|_| invalid_inputs())?;
    let run: RunRequest = strict_json::parse(request)?;
    let result: RulesResult = parse_result(bytes, run.budget.max_output_bytes)?;
    // What the engine reports for these inputs, from an execution that is
    // cancelled before it simulates anything.
    let cancelled = Control::default();
    cancelled.cancel();
    let header = rules::execute(&prepared, &cancelled).map_err(|_| invalid_inputs())?;
    check_identity(&result, &header, model, request)?;
    let records = check_records(&result, &run)?;
    check_ensemble(&result, &header, &records)?;
    check_kept(&result, &run, &records, model)?;
    check_statistics(&result, &run, &records)?;
    Ok(success(
        &result.schema,
        Some(result.termination),
        model,
        request,
    ))
}

fn same(a: &impl Serialize, b: &impl Serialize) -> bool {
    matches!((serde_json::to_value(a), serde_json::to_value(b)), (Ok(a), Ok(b)) if a == b)
}

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

fn check_identity(
    result: &RulesResult,
    header: &RulesResult,
    model: &[u8],
    request: &[u8],
) -> Result<(), Diagnostic> {
    if result.schema != header.schema
        || result.model_sha256 != sha256(model)
        || result.run_sha256 != sha256(request)
        || result.model_sha256 != header.model_sha256
        || result.run_sha256 != header.run_sha256
        || result.sampler != header.sampler
        || result.root_seed != header.root_seed
        || !same(&result.types, &header.types)
        || result.rules != header.rules
        || result.observables != header.observables
        || !same_bits(&result.report_times, &header.report_times)
        || result.requested_replicates != header.requested_replicates
        || result.replicates.len() != header.replicates.len()
        || result.method != header.method
        || result.continuous != header.continuous
    {
        return Err(error("result identity or dimensions differ from inputs"));
    }
    Ok(())
}

/// What the records establish: which replicates started and how many
/// report times each filled, and whether the ensemble is complete.
struct Records {
    started: Vec<bool>,
    saved: Vec<usize>,
    complete: bool,
}

fn check_records(result: &RulesResult, run: &RunRequest) -> Result<Records, Diagnostic> {
    let rows = run.time.report_times.len();
    let end = run.time.end;
    let mut completed = 0;
    let mut stopped = 0;
    let mut events = 0_u64;
    let mut aggregate = Termination::Completed;
    for (id, r) in result.replicates.iter().enumerate() {
        let fired = r
            .rule_events
            .iter()
            .try_fold(0_u64, |total, &n| total.checked_add(n));
        if r.id != id as u64
            || !r.stop_time.is_finite()
            || r.stop_time < 0.0
            || r.stop_time > end
            || r.events > run.budget.max_events_per_replicate
            || r.saved_points > rows
            || r.rule_events.len() != result.rules.len()
            || fired != Some(r.events)
            || r.termination == Termination::OutputLimit
            || (!r.started
                && (r.events != 0
                    || r.saved_points != 0
                    || r.rng_words != 0
                    || r.peak_instances != 0
                    || r.created_instances != 0
                    || r.stop_time != 0.0
                    || r.diagnostic.is_some()
                    || !matches!(
                        r.termination,
                        Termination::Cancelled | Termination::WallTimeLimit
                    )))
        {
            return Err(error("invalid replicate accounting"));
        }
        if r.termination.complete() {
            if !r.started
                || r.saved_points != rows
                || r.diagnostic.is_some()
                || (r.termination == Termination::Completed && r.stop_time != end)
            {
                return Err(error("incomplete replicate marked complete"));
            }
            completed += 1;
            if r.termination == Termination::Stopped {
                stopped += 1;
            }
        } else if r.started && r.diagnostic.is_none() {
            return Err(error("failed or halted replicate has no diagnostic"));
        }
        events = events
            .checked_add(r.events)
            .ok_or_else(|| error("event accounting overflow"))?;
        aggregate = aggregate.max(if r.termination == Termination::Stopped {
            Termination::Completed
        } else {
            r.termination
        });
    }
    if result.completed_replicates != completed
        || result.stopped_replicates != stopped
        || result.events != events
    {
        return Err(error("aggregate accounting mismatch"));
    }
    let output_limited = result.termination == Termination::OutputLimit;
    if !output_limited && result.termination != aggregate {
        return Err(error("termination evidence mismatch"));
    }
    Ok(Records {
        started: result.replicates.iter().map(|r| r.started).collect(),
        saved: result.replicates.iter().map(|r| r.saved_points).collect(),
        complete: completed == result.requested_replicates,
    })
}

/// The ensemble diagnostics: `CRN-R050` for an incomplete ensemble, the
/// `CRN-R051` warnings of these inputs and `CRN-R045` for a trimmed result,
/// in the engine's order; a trimmed result keeps no numerical payload.
fn check_ensemble(
    result: &RulesResult,
    header: &RulesResult,
    records: &Records,
) -> Result<(), Diagnostic> {
    let output_limited = result.termination == Termination::OutputLimit;
    let incomplete = header.diagnostics.iter().find(|d| d.code == "CRN-R050");
    let mut expected: Vec<&Diagnostic> = Vec::new();
    if !records.complete {
        expected.extend(incomplete);
    }
    expected.extend(header.diagnostics.iter().filter(|d| d.code != "CRN-R050"));
    let trimmed = result.diagnostics.last().filter(|d| d.code == "CRN-R045");
    if output_limited {
        expected.extend(trimmed);
    }
    if (!records.complete && incomplete.is_none())
        || (output_limited
            && (trimmed.is_none()
                || result.statistics.is_some()
                || !result.paths.is_empty()
                || !result.observations.is_empty()
                || !result.final_instances.is_empty()))
        || !same(&result.diagnostics, &expected)
    {
        return Err(error("termination evidence mismatch"));
    }
    Ok(())
}

fn check_kept(
    result: &RulesResult,
    run: &RunRequest,
    records: &Records,
    model: &[u8],
) -> Result<(), Diagnostic> {
    let output_limited = result.termination == Termination::OutputLimit;
    let kept = |limit: u64| -> Vec<u64> {
        if output_limited {
            return Vec::new();
        }
        (0..result.requested_replicates.min(limit))
            .filter(|&id| records.started[id as usize])
            .collect()
    };
    let columns = result.observables.len();
    if result
        .paths
        .iter()
        .map(|p| p.id)
        .ne(kept(run.outputs.keep_paths))
        || result.paths.iter().any(|p| {
            p.values.len() != records.saved[p.id as usize]
                || p.values.iter().any(|row| row.len() != columns)
        })
    {
        return Err(error("kept path identity or length mismatch"));
    }
    if run.time.report_times[0] == 0.0 && result.paths.iter().any(|p| !p.values.is_empty()) {
        if let Some(initial) = initial_row(model, run)? {
            if result
                .paths
                .iter()
                .filter_map(|p| p.values.first())
                .any(|row| !same_bits(row, &initial))
            {
                return Err(error("kept path starts at the wrong initial state"));
            }
        }
    }
    if result
        .final_instances
        .iter()
        .map(|f| f.id)
        .ne(kept(run.outputs.keep_instances))
        || result.final_instances.iter().any(|f| {
            f.instances.iter().any(|i| {
                i.count == 0
                    || result.types.get(i.type_index).is_none_or(|t| {
                        t.attributes.len() != i.attributes.len()
                            || i.attributes.iter().any(|a| !a.is_finite())
                    })
            })
        })
    {
        return Err(error("final instance identity or content mismatch"));
    }
    let mut previous = None;
    for observation in &result.observations {
        let id = observation.id;
        if previous.is_some_and(|p| p >= id)
            || id >= result.requested_replicates
            || !records.started[id as usize]
            || observation.trackers.len() != run.trackers().len()
        {
            return Err(error("observations do not follow the started replicates"));
        }
        previous = Some(id);
    }
    Ok(())
}

/// The first reported row of a replicate of these inputs: the initial
/// observable values, from one replicate that may fire at most one event,
/// reported at time 0 only. `None` when that replicate reports no row.
fn initial_row(model: &[u8], run: &RunRequest) -> Result<Option<Vec<f64>>, Diagnostic> {
    let mut probe = run.clone();
    probe.ensemble.replicates = 1;
    probe.time.report_times.truncate(1);
    probe.outputs.keep_paths = 1;
    probe.outputs.keep_instances = 0;
    probe.budget.max_events_per_replicate = 1;
    probe.budget.threads = 1;
    probe.budget.wall_time_s = None;
    let bytes = serde_json::to_vec(&probe).map_err(|_| invalid_inputs())?;
    let prepared = rules::prepare(model, &bytes).map_err(|_| invalid_inputs())?;
    let probed = rules::execute(&prepared, &Control::default()).map_err(|_| invalid_inputs())?;
    Ok(probed
        .paths
        .into_iter()
        .next()
        .and_then(|p| p.values.into_iter().next()))
}

#[allow(clippy::cast_precision_loss)]
fn check_statistics(
    result: &RulesResult,
    run: &RunRequest,
    records: &Records,
) -> Result<(), Diagnostic> {
    let output_limited = result.termination == Termination::OutputLimit;
    let probabilities = run.outputs.quantiles();
    let wanted = !run.outputs.statistics.is_empty() || !probabilities.is_empty();
    if result.statistics.is_some() != (records.complete && wanted && !output_limited) {
        return Err(error(
            "statistics present for incomplete work or absent for completed work",
        ));
    }
    let Some(stats) = &result.statistics else {
        return Ok(());
    };
    let n = result.requested_replicates;
    let mean = run.outputs.statistics.contains(&Statistic::Mean);
    let variance = run.outputs.statistics.contains(&Statistic::Variance) && n > 1;
    let quantiles = stats.quantiles.as_deref().unwrap_or(&[]);
    if stats.samples != n
        || stats.mean.is_some() != mean
        || stats.sample_variance.is_some() != variance
        || stats.mean_standard_error.is_some() != variance
        || stats.quantiles.is_some() != !probabilities.is_empty()
        || !same_bits(
            &quantiles.iter().map(|q| q.p).collect::<Vec<_>>(),
            probabilities,
        )
    {
        return Err(error("statistic estimator/sample count mismatch"));
    }
    let columns = result.observables.len();
    // Without observables every matrix is empty (see `statistics`).
    let rows = if columns == 0 {
        0
    } else {
        run.time.report_times.len()
    };
    let flat = |m: &Vec<Vec<f64>>| -> Option<Vec<f64>> {
        (m.len() == rows && m.iter().all(|r| r.len() == columns))
            .then(|| m.iter().flatten().copied().collect())
    };
    let matrices: Vec<Option<Vec<f64>>> = [
        &stats.mean,
        &stats.sample_variance,
        &stats.mean_standard_error,
    ]
    .into_iter()
    .map(|m| m.as_ref().map(flat).unwrap_or(Some(Vec::new())))
    .chain(quantiles.iter().map(|q| flat(&q.values)))
    .collect();
    if matrices
        .iter()
        .any(|m| m.as_ref().is_none_or(|v| v.iter().any(|x| !x.is_finite())))
    {
        return Err(error("invalid statistic dimensions or values"));
    }
    let matrices: Vec<Vec<f64>> = matrices.into_iter().flatten().collect();
    let (mean_values, variances, errors) = (&matrices[0], &matrices[1], &matrices[2]);
    let samples = n as f64;
    if variances.iter().chain(errors).any(|x| *x < 0.0)
        || (variance
            && !same_bits(
                errors,
                &variances
                    .iter()
                    .map(|v| libm::sqrt(v / samples))
                    .collect::<Vec<_>>(),
            ))
        || matrices[3..]
            .windows(2)
            .any(|w| w[0].iter().zip(&w[1]).any(|(a, b)| a > b))
    {
        return Err(error("inconsistent statistics"));
    }
    if result.paths.len() as u64 != n {
        return Ok(());
    }
    // Every path is kept: recompute in replicate order, as the engine folds.
    let cells = rows * columns;
    let mut means = vec![0.0; cells];
    let mut m2 = vec![0.0; cells];
    let mut buffer = Vec::new();
    for (k, path) in result.paths.iter().enumerate() {
        let count = (k + 1) as f64;
        for ((m, s), &x) in means
            .iter_mut()
            .zip(&mut m2)
            .zip(path.values.iter().flatten())
        {
            let delta = x - *m;
            *m += delta / count;
            *s += delta * (x - *m);
        }
        if !probabilities.is_empty() {
            buffer.extend(path.values.iter().flatten());
        }
    }
    let sample_variances: Vec<f64> = m2.iter().map(|s| s / n.saturating_sub(1) as f64).collect();
    let recomputed = if probabilities.is_empty() {
        Vec::new()
    } else {
        rules::quantile::statistics(
            &buffer,
            usize::try_from(n).map_err(|_| invalid_inputs())?,
            rows,
            columns,
            probabilities,
        )
    };
    if (mean && !same_bits(mean_values, &means))
        || (variance && !same_bits(variances, &sample_variances))
        || quantiles.iter().zip(&recomputed).any(|(a, b)| {
            !same_bits(
                &a.values.iter().flatten().copied().collect::<Vec<_>>(),
                &b.values.iter().flatten().copied().collect::<Vec<_>>(),
            )
        })
    {
        return Err(error("statistics disagree with the complete kept ensemble"));
    }
    Ok(())
}

fn verify_fsp(model: &[u8], request: &[u8], bytes: &[u8]) -> Result<ArchiveReport, Diagnostic> {
    let run: fsp::FspRequest = strict_json::parse(request).map_err(|_| invalid_inputs())?;
    let view = FlatView::from_bytes(model).map_err(|_| invalid_inputs())?;
    if run.model.sha256 != sha256(model) {
        return Err(invalid_inputs());
    }
    let result: fsp::FspResult = parse_result(bytes, run.budget.max_output_bytes)?;
    let fsp::FspMethod::Fsp { bounds, .. } = &run.method;
    let species = view.species();
    let limits: Option<Vec<u64>> = species.iter().map(|s| bounds.get(s).copied()).collect();
    let n = species.len();
    let states = result.states.len();
    let consistent = limits.is_some_and(|limits| {
        result
            .states
            .iter()
            .all(|s| s.len() == n && s.iter().zip(&limits).all(|(x, max)| x <= max))
    });
    if result.model_sha256 != sha256(model)
        || result.run_sha256 != sha256(request)
        || result.species != species
        || result.states.first().map(Vec::as_slice) != Some(view.initial())
        || !consistent
        || !same_bits(
            &result.snapshots.iter().map(|s| s.time).collect::<Vec<_>>(),
            &run.time.report_times,
        )
        || result.snapshots.iter().any(|s| {
            s.probabilities.len() != states
                || s.probabilities
                    .iter()
                    .any(|p| !p.is_finite() || *p < 0.0 || *p > 1.0)
                || s.projected_first_moment.len() != n
                || s.projected_second_raw_moment.len() != n
        })
        || result.sparse_edges > states.saturating_mul(states.saturating_add(1))
        || result.work > MAX_FSP_WORK
        || !result.uniformization_rate.is_finite()
        || result.uniformization_rate < 0.0
        || result.numerical_error_certified
    {
        return Err(error("result identity or dimensions differ from inputs"));
    }
    Ok(success(&result.schema, None, model, request))
}
