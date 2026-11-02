//! Exact stochastic simulation of dynamic, attributed species.
//!
//! A model declares species *types* with numeric attributes. An instance is a
//! type plus one value per attribute and a copy number; instances appear when
//! a rule produces new attribute values and disappear when their count reaches
//! zero. Rules are schemata over zero, one or two reactant patterns with
//! product-form propensities (a global rate times local factors), simulated
//! exactly with Gillespie's direct method over sum trees.
//!
//! Contracts: `crn-craft.rules-model/0.1`, `crn-craft.rules-run/0.1`,
//! `crn-craft.rules-result/0.1`. See `docs/rules-language.md`.

mod engine;
mod expr;
pub mod model;
mod tree;

use crn_craft::contracts::{Ensemble, ModelRef, Statistic};
use crn_craft::diagnostic::Diagnostic;
use rand_chacha::ChaCha8Rng;
use rand_core::SeedableRng;
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

pub use model::Model;

pub const SAMPLER: &str = "chacha8-u53/rules-direct/1";
pub const MAX_REPLICATES: u64 = 1_000_000;
pub const MAX_REPORT_TIMES: usize = 100_000;
pub const MAX_WORKERS: u32 = 64;
pub const MAX_RESULT_BYTES: u64 = 1 << 30;

/// Per-replicate ChaCha8 stream in the rules domain: a SHA-256 key of the
/// domain and root seed, with the replicate index as stream id.
pub fn stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    let mut digest = Sha256::new();
    digest.update(b"crn-craft.rules/0.1");
    digest.update(root_seed.to_le_bytes());
    let mut rng = ChaCha8Rng::from_seed(digest.finalize().into());
    rng.set_stream(replicate);
    rng
}

/// Cooperative cancellation, by flag or by the existence of a file.
#[derive(Debug, Default)]
pub struct Control {
    cancelled: AtomicBool,
    pub cancel_file: Option<PathBuf>,
}
impl Control {
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    pub fn is_cancelled(&self) -> bool {
        if self.cancelled.load(Ordering::Relaxed) {
            return true;
        }
        if self.cancel_file.as_ref().is_some_and(|p| p.exists()) {
            self.cancel();
            return true;
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Run request

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum RunVersion {
    #[serde(rename = "crn-craft.rules-run/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Time {
    /// Simulated time horizon, positive and finite.
    pub end: f64,
    /// Strictly increasing times in [0, end] at which observables are
    /// reported (right-continuous: events at a report time are included).
    pub report_times: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Outputs {
    pub statistics: Vec<Statistic>,
    /// Observable series of replicates 0..keep_paths.
    pub keep_paths: u64,
    /// Final instances of replicates 0..keep_instances.
    pub keep_instances: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    pub max_events_per_replicate: u64,
    /// Every random word, including rejected pair draws and product draws.
    pub max_sampling_steps_per_replicate: u64,
    /// Distinct instances alive at once, over all types.
    pub max_instances: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    #[schemars(range(min = 1, max = 64))]
    pub threads: u32,
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunRequest {
    pub schema: RunVersion,
    pub model: ModelRef,
    pub time: Time,
    pub ensemble: Ensemble,
    pub outputs: Outputs,
    pub budget: Budget,
}

// ---------------------------------------------------------------------------
// Result

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.rules-result/0.1")]
    V01,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    Completed,
    /// The stop condition held; the state is held to the end of the grid.
    Stopped,
    EventLimit,
    SamplingLimit,
    InstanceLimit,
    WallTimeLimit,
    Cancelled,
    Failed,
    OutputLimit,
}
impl Termination {
    pub fn complete(self) -> bool {
        matches!(self, Self::Completed | Self::Stopped)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TypeHeader {
    pub name: String,
    pub attributes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct InstanceRecord {
    /// Index into `types`.
    pub type_index: usize,
    /// One value per attribute of the type, in declared order.
    pub attributes: Vec<f64>,
    pub count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub id: u64,
    pub started: bool,
    pub termination: Termination,
    pub events: u64,
    /// End of the horizon, the time the stop condition fired, or the time a
    /// limit or failure ended the path.
    pub stop_time: f64,
    pub saved_points: usize,
    /// Events per rule, in rule order.
    pub rule_events: Vec<u64>,
    pub rng_words: u64,
    /// Rejected draws of the same molecule twice for identical patterns.
    pub rejected_pairs: u64,
    pub peak_instances: u64,
    /// Instances created by events (excluding the initial state).
    pub created_instances: u64,
    pub diagnostic: Option<Diagnostic>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PathSeries {
    pub id: u64,
    /// report time × observable.
    pub values: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FinalInstances {
    pub id: u64,
    pub instances: Vec<InstanceRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Statistics {
    pub samples: u64,
    /// report time × observable.
    pub mean: Option<Vec<Vec<f64>>>,
    pub sample_variance: Option<Vec<Vec<f64>>>,
    pub mean_standard_error: Option<Vec<Vec<f64>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RulesResult {
    pub schema: ResultVersion,
    pub model_sha256: String,
    pub run_sha256: String,
    pub sampler: String,
    pub root_seed: u64,
    pub types: Vec<TypeHeader>,
    pub rules: Vec<String>,
    pub observables: Vec<String>,
    pub report_times: Vec<f64>,
    /// `completed` when every replicate completed or stopped; otherwise the
    /// most severe replicate termination.
    pub termination: Termination,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub stopped_replicates: u64,
    pub events: u64,
    pub replicates: Vec<Record>,
    pub paths: Vec<PathSeries>,
    pub final_instances: Vec<FinalInstances>,
    /// Withheld unless every replicate completed or stopped.
    pub statistics: Option<Statistics>,
    pub diagnostics: Vec<Diagnostic>,
}

// ---------------------------------------------------------------------------

pub struct Prepared {
    compiled: model::Compiled,
    request: RunRequest,
    model_sha256: String,
    run_sha256: String,
}

fn request_error(path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
    Diagnostic::error("CRN-R020", path, message, hint)
}

fn validate_request(request: &RunRequest, model_sha256: &str) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if !crn_craft::validate::valid_relative_path(&request.model.path) {
        out.push(request_error("/model/path", "non-portable or non-relative model path", "Use a slash-separated path below the request directory without '.', '..', empty components, drive prefixes or backslashes."));
    }
    if request.model.sha256 != model_sha256 {
        out.push(request_error(
            "/model/sha256",
            "model SHA-256 does not match the supplied model",
            "Recompute the hash of the exact model bytes.",
        ));
    }
    let time = &request.time;
    if !(time.end.is_finite() && time.end > 0.0) {
        out.push(request_error(
            "/time/end",
            "end must be positive and finite",
            "Use a positive finite horizon.",
        ));
    }
    if time.report_times.is_empty() || time.report_times.len() > MAX_REPORT_TIMES {
        out.push(request_error(
            "/time/report_times",
            "between 1 and 100000 report times are required",
            "Give a sorted report grid.",
        ));
    }
    for (i, &t) in time.report_times.iter().enumerate() {
        let ordered = i == 0 || t > time.report_times[i - 1];
        if !(t.is_finite() && t >= 0.0 && t <= time.end && ordered) {
            out.push(request_error(
                &format!("/time/report_times/{i}"),
                "report times must be strictly increasing within [0, end]",
                "Sort the grid and keep it inside the horizon.",
            ));
            break;
        }
    }
    let replicates = request.ensemble.replicates;
    if replicates == 0 || replicates > MAX_REPLICATES {
        out.push(request_error(
            "/ensemble/replicates",
            "replicates must be between 1 and 1000000",
            "Use a supported ensemble size.",
        ));
    }
    if request.outputs.keep_paths > replicates || request.outputs.keep_instances > replicates {
        out.push(request_error(
            "/outputs",
            "cannot keep more replicates than requested",
            "Keep at most `replicates` paths and instance sets.",
        ));
    }
    let mut seen = std::collections::BTreeSet::new();
    if !request
        .outputs
        .statistics
        .iter()
        .all(|s| seen.insert(s.clone()))
    {
        out.push(request_error(
            "/outputs/statistics",
            "duplicate statistic",
            "List each statistic once.",
        ));
    }
    let budget = &request.budget;
    if budget.threads == 0 || budget.threads > MAX_WORKERS {
        out.push(request_error(
            "/budget/threads",
            "threads must be between 1 and 64",
            "Use a supported worker count.",
        ));
    }
    if budget.max_output_bytes == 0 {
        out.push(request_error(
            "/budget/max_output_bytes",
            "output budget must be positive",
            "Give an explicit output budget.",
        ));
    }
    if budget.max_instances == 0 {
        out.push(request_error(
            "/budget/max_instances",
            "instance budget must be positive",
            "Give an explicit instance budget.",
        ));
    }
    if budget
        .wall_time_s
        .is_some_and(|w| !(w.is_finite() && w > 0.0))
    {
        out.push(request_error(
            "/budget/wall_time_s",
            "wall time must be positive and finite",
            "Omit it or give a positive number of seconds.",
        ));
    }
    out
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum CheckVersion {
    #[serde(rename = "crn-craft.rules-check/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RulesCheck {
    pub schema: CheckVersion,
    pub valid: bool,
    /// SHA-256 of the model bytes, when they could be read.
    pub model_sha256: Option<String>,
    pub diagnostics: Vec<Diagnostic>,
}

/// Parses and checks a model; every expression is compiled and name-checked.
pub fn check_model(model_bytes: &[u8]) -> RulesCheck {
    let diagnostics = crn_craft::strict_json::parse::<Model>(model_bytes)
        .map_err(|e| vec![e])
        .and_then(|model| model::compile(&model).map(|_| ()))
        .err()
        .unwrap_or_default();
    RulesCheck {
        schema: CheckVersion::V01,
        valid: diagnostics.is_empty(),
        model_sha256: Some(crn_craft::sha256(model_bytes)),
        diagnostics,
    }
}

pub fn prepare(model_bytes: &[u8], request_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let model: Model = crn_craft::strict_json::parse(model_bytes).map_err(|e| vec![e])?;
    let request: RunRequest = crn_craft::strict_json::parse(request_bytes).map_err(|e| vec![e])?;
    let model_sha256 = crn_craft::sha256(model_bytes);
    let mut diagnostics = validate_request(&request, &model_sha256);
    let compiled = match model::compile(&model) {
        Ok(c) => Some(c),
        Err(mut d) => {
            diagnostics.append(&mut d);
            None
        }
    };
    match compiled {
        Some(compiled) if diagnostics.is_empty() => Ok(Prepared {
            compiled,
            request,
            model_sha256,
            run_sha256: crn_craft::sha256(request_bytes),
        }),
        _ => Err(diagnostics),
    }
}

fn record(id: u64, path: &engine::PathOutput) -> Record {
    Record {
        id,
        started: true,
        termination: path.termination,
        events: path.events,
        stop_time: path.stop_time,
        saved_points: path.saved_points,
        rule_events: path.rule_events.clone(),
        rng_words: path.rng_words,
        rejected_pairs: path.rejected_pairs,
        peak_instances: path.peak_instances,
        created_instances: path.created_instances,
        diagnostic: path.diagnostic.clone(),
    }
}

fn fits_json(value: &impl Serialize, limit: u64) -> bool {
    struct Counter {
        used: u64,
        limit: u64,
    }
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            let n = bytes.len() as u64;
            if n > self.limit.saturating_sub(self.used) {
                return Err(std::io::Error::other("output limit"));
            }
            self.used += n;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    // The budget covers the compact JSON and a trailing newline.
    let mut counter = Counter {
        used: 0,
        limit: limit.saturating_sub(1),
    };
    serde_json::to_writer(&mut counter, value).is_ok()
}

#[allow(clippy::cast_precision_loss)]
fn statistics(
    wanted: &[Statistic],
    paths: &[engine::PathOutput],
    rows: usize,
    columns: usize,
) -> Statistics {
    let n = paths.len();
    let chunk = |flat: Vec<f64>| -> Vec<Vec<f64>> {
        flat.chunks(columns.max(1))
            .take(rows)
            .map(<[f64]>::to_vec)
            .collect()
    };
    let cells = rows * columns;
    let mut mean = vec![0.0; cells];
    for p in paths {
        for (m, v) in mean.iter_mut().zip(&p.values) {
            *m += v;
        }
    }
    for m in &mut mean {
        *m /= n as f64;
    }
    let mut m2 = vec![0.0; cells];
    for p in paths {
        for ((s, v), m) in m2.iter_mut().zip(&p.values).zip(&mean) {
            *s += (v - m) * (v - m);
        }
    }
    let variance: Vec<f64> = m2
        .iter()
        .map(|s| if n > 1 { s / (n - 1) as f64 } else { f64::NAN })
        .collect();
    let has_variance = wanted.contains(&Statistic::Variance) && n > 1;
    Statistics {
        samples: n as u64,
        mean: wanted.contains(&Statistic::Mean).then(|| chunk(mean)),
        mean_standard_error: has_variance
            .then(|| chunk(variance.iter().map(|v| libm::sqrt(v / n as f64)).collect())),
        sample_variance: has_variance.then(|| chunk(variance)),
    }
}

/// Runs the ensemble. Replicates are independent keyed streams, so results do
/// not depend on the thread count.
#[allow(clippy::missing_panics_doc)]
pub fn execute(prepared: &Prepared, control: &Control) -> Result<RulesResult, Diagnostic> {
    let model = &prepared.compiled;
    let request = &prepared.request;
    let start = Instant::now();
    let deadline = request
        .budget
        .wall_time_s
        .map(|s| start + Duration::from_secs_f64(s.min(1e9)));
    let limits = engine::Limits {
        max_events: request.budget.max_events_per_replicate,
        max_words: request.budget.max_sampling_steps_per_replicate,
        max_instances: request.budget.max_instances,
        deadline,
    };
    let n = request.ensemble.replicates;
    let keep_instances = request.outputs.keep_instances;
    let next = AtomicU64::new(0);
    let halt = AtomicBool::new(false);
    let finished: Mutex<Vec<(u64, engine::PathOutput)>> = Mutex::new(Vec::new());
    let workers = u64::from(request.budget.threads).min(n).max(1);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| {
                let mut local = Vec::new();
                loop {
                    if halt.load(Ordering::Relaxed) {
                        break;
                    }
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    if id >= n {
                        break;
                    }
                    let path = engine::Path {
                        model,
                        report_times: &request.time.report_times,
                        end: request.time.end,
                        limits: &limits,
                        control,
                        keep_instances: id < keep_instances,
                    };
                    let output = path.run(&mut stream(request.ensemble.root_seed, id));
                    if matches!(
                        output.termination,
                        Termination::Cancelled | Termination::WallTimeLimit
                    ) {
                        halt.store(true, Ordering::Relaxed);
                    }
                    local.push((id, output));
                }
                finished
                    .lock()
                    .expect("no poisoned workers")
                    .append(&mut local);
            });
        }
    });
    let mut finished = finished.into_inner().expect("no poisoned workers");
    finished.sort_by_key(|(id, _)| *id);
    let control_stop = if control.is_cancelled() {
        Termination::Cancelled
    } else {
        Termination::WallTimeLimit
    };
    let mut replicates = Vec::with_capacity(usize::try_from(n).unwrap_or(0));
    let mut done = finished.iter().peekable();
    for id in 0..n {
        if let Some((_, path)) = done.next_if(|(i, _)| *i == id) {
            replicates.push(record(id, path));
        } else {
            replicates.push(Record {
                id,
                started: false,
                termination: control_stop,
                events: 0,
                stop_time: 0.0,
                saved_points: 0,
                rule_events: vec![0; model.rules.len()],
                rng_words: 0,
                rejected_pairs: 0,
                peak_instances: 0,
                created_instances: 0,
                diagnostic: None,
            });
        }
    }
    let complete = replicates.iter().all(|r| r.termination.complete());
    let termination = replicates
        .iter()
        .map(|r| {
            if r.termination == Termination::Stopped {
                Termination::Completed
            } else {
                r.termination
            }
        })
        .max()
        .unwrap_or(Termination::Failed);
    let columns = model.observables.len();
    let rows = request.time.report_times.len();
    let to_rows = |values: &[f64]| -> Vec<Vec<f64>> {
        if columns == 0 {
            vec![Vec::new(); rows]
        } else {
            values.chunks(columns).map(<[f64]>::to_vec).collect()
        }
    };
    let paths = finished
        .iter()
        .filter(|(id, _)| *id < request.outputs.keep_paths)
        .map(|(id, p)| PathSeries {
            id: *id,
            values: to_rows(&p.values[..p.saved_points * columns]),
        })
        .collect();
    let final_instances = finished
        .iter()
        .filter_map(|(id, p)| {
            p.instances.as_ref().map(|instances| FinalInstances {
                id: *id,
                instances: instances.clone(),
            })
        })
        .collect();
    let all: Vec<engine::PathOutput> = finished.into_iter().map(|(_, p)| p).collect();
    let statistics = (complete && !request.outputs.statistics.is_empty())
        .then(|| statistics(&request.outputs.statistics, &all, rows, columns));
    let mut diagnostics = Vec::new();
    if !complete {
        diagnostics.push(Diagnostic::warning(
            "CRN-R050",
            "/replicates",
            "ensemble incomplete; statistics withheld",
            "See the per-replicate terminations and diagnostics.",
        ));
    }
    let mut result = RulesResult {
        schema: ResultVersion::V01,
        model_sha256: prepared.model_sha256.clone(),
        run_sha256: prepared.run_sha256.clone(),
        sampler: SAMPLER.into(),
        root_seed: request.ensemble.root_seed,
        types: model
            .types
            .iter()
            .map(|t| TypeHeader {
                name: t.name.clone(),
                attributes: t.attributes.clone(),
            })
            .collect(),
        rules: model.rules.iter().map(|r| r.name.clone()).collect(),
        observables: model.observables.iter().map(|o| o.name.clone()).collect(),
        report_times: request.time.report_times.clone(),
        termination,
        requested_replicates: n,
        completed_replicates: replicates
            .iter()
            .filter(|r| r.termination.complete())
            .count() as u64,
        stopped_replicates: replicates
            .iter()
            .filter(|r| r.termination == Termination::Stopped)
            .count() as u64,
        events: replicates.iter().map(|r| r.events).sum(),
        replicates,
        paths,
        final_instances,
        statistics,
        diagnostics,
    };
    let limit = request.budget.max_output_bytes.min(MAX_RESULT_BYTES);
    if !fits_json(&result, limit) {
        result.termination = Termination::OutputLimit;
        result.statistics = None;
        result.paths.clear();
        result.final_instances.clear();
        result.diagnostics.push(Diagnostic::error(
            "CRN-R045",
            "/outputs",
            "rules output exceeds its byte budget",
            "Increase max_output_bytes or keep fewer paths; statistics, paths and instances have been removed.",
        ));
        if !fits_json(&result, limit) {
            return Err(Diagnostic::error(
                "CRN-R045",
                "/budget/max_output_bytes",
                "replicate records exceed the output budget",
                "Increase max_output_bytes.",
            ));
        }
    }
    Ok(result)
}

pub fn solve(
    model_bytes: &[u8],
    request_bytes: &[u8],
    control: &Control,
) -> Result<RulesResult, Vec<Diagnostic>> {
    execute(&prepare(model_bytes, request_bytes)?, control).map_err(|e| vec![e])
}

fn read(path: &Path) -> Result<Vec<u8>, Diagnostic> {
    let error = |message: String| {
        Diagnostic::error(
            "CRN-E030",
            "",
            message,
            "Provide a readable regular JSON file.",
        )
    };
    let meta = std::fs::metadata(path).map_err(|e| error(format!("{}: {e}", path.display())))?;
    if !meta.is_file() {
        return Err(error(format!("{} is not a regular file", path.display())));
    }
    if meta.len() > crn_craft::strict_json::MAX_DOCUMENT_BYTES as u64 {
        return Err(Diagnostic::error(
            "CRN-E004",
            "",
            format!("{} exceeds the document size limit", path.display()),
            "Keep input documents below 16 MiB.",
        ));
    }
    std::fs::read(path).map_err(|e| error(format!("{}: {e}", path.display())))
}

/// Runs a request file against the model it references, which must be the
/// supplied model file below the request's directory.
pub fn solve_files(
    model_path: &Path,
    request_path: &Path,
    control: &Control,
) -> Result<RulesResult, Vec<Diagnostic>> {
    let model = read(model_path).map_err(|e| vec![e])?;
    let bytes = read(request_path).map_err(|e| vec![e])?;
    let request: RunRequest = crn_craft::strict_json::parse(&bytes).map_err(|e| vec![e])?;
    let parent = request_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let resolved = (
        parent.canonicalize(),
        parent.join(&request.model.path).canonicalize(),
        model_path.canonicalize(),
    );
    if !crn_craft::validate::valid_relative_path(&request.model.path)
        || !matches!(resolved, (Ok(ref base), Ok(ref reference), Ok(ref supplied)) if reference.starts_with(base) && reference == supplied)
    {
        return Err(vec![request_error(
            "/model/path",
            "rules request does not identify the supplied model below its directory",
            "Use the referenced model within the request directory.",
        )]);
    }
    solve(&model, &bytes, control)
}

pub fn check_model_file(model_path: &Path) -> RulesCheck {
    match read(model_path) {
        Ok(bytes) => check_model(&bytes),
        Err(e) => RulesCheck {
            schema: CheckVersion::V01,
            valid: false,
            model_sha256: None,
            diagnostics: vec![e],
        },
    }
}

/// JSON Schemas generated from the Rust contract types.
pub fn schema(kind: &str) -> Option<Value> {
    match kind {
        "rules-model" => serde_json::to_value(schema_for!(Model)).ok(),
        "rules-check" => serde_json::to_value(schema_for!(RulesCheck)).ok(),
        "rules-run" => serde_json::to_value(schema_for!(RunRequest)).ok(),
        "rules-result" => serde_json::to_value(schema_for!(RulesResult)).ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
