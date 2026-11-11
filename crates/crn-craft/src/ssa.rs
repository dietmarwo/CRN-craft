//! Versioned alternative exact SSA methods. Legacy Direct contracts are unchanged.
//!
//! The `ssa` entry point: a `crn-craft.ssa-run/0.1`–`0.4` [`Request`] selects
//! an [`Algorithm`] (`direct`, `cached_direct`, `direct_cr`,
//! `sorting_direct`, `rssa_cr`, and `next_reaction` from 0.4), a
//! [`WaitingTime`] sampler and, from 0.2, the `fast` numerical policy.
//! [`prepare`] validates the request against the model's exact bytes,
//! [`execute`] runs it and returns the matching `ssa-result` version as an
//! [`SsaResult`], a standalone JSON result rather than an archive.
//! [`solve`] combines both, and [`solve_files`] also checks that the
//! request's `model.path` names the model file.
//!
//! `fast` is not allowed in 0.1, defaults to false in 0.2 and to true in 0.3
//! and 0.4, except for `next_reaction`, which runs the checked policy only.
//! The checked policy uses ChaCha8 words and rejects lost increments and
//! channels below 53-bit resolution; `fast` emulates JumpProcesses.jl with
//! Xoshiro256++ words and running totals without those checks. Each
//! combination has its own sampler identity, recorded in
//! [`SsaResult::sampler`]. Every replicate also has a random-word budget
//! (`max_sampling_steps_per_replicate`); exhausting it stops the replicate
//! with `sampling_limit` and `CRN-E067`. See
//! [docs/alternative-ssa.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/alternative-ssa.md).
//!
//! ```
//! use crn_craft::{sha256, simulation::Control, ssa};
//!
//! let model = br#"{"schema": "crn-craft.model/0.1", "name": "birth_death",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
//!   "species": [{"name": "X", "initial_count": 0}],
//!   "reactions": [
//!     {"name": "make", "reactants": {}, "products": {"X": 1},
//!      "rate": {"kind": "mass_action", "coefficient": 5.0}},
//!     {"name": "decay", "reactants": {"X": 1}, "products": {},
//!      "rate": {"kind": "mass_action", "coefficient": 1.0}}]}"#;
//! let request = format!(r#"{{"schema": "crn-craft.ssa-run/0.3",
//!   "model": {{"path": "model.json", "sha256": "{}"}},
//!   "method": {{"kind": "ssa", "aggregator": "rssa_cr", "waiting_time": "inverse", "fast": false}},
//!   "time": {{"end": 2.0, "save_times": [0.0, 1.0, 2.0]}},
//!   "ensemble": {{"replicates": 8, "root_seed": 11}},
//!   "outputs": {{"statistics": ["mean"], "keep_trajectories": 0}},
//!   "budget": {{"max_events_per_replicate": 10000, "max_sampling_steps_per_replicate": 1000000,
//!              "threads": 1, "max_output_bytes": 100000}}}}"#, sha256(model));
//!
//! let prepared = ssa::prepare(model, request.as_bytes()).unwrap();
//! let result = ssa::execute(&prepared, &Control::default()).unwrap();
//! assert_eq!(result.termination, ssa::Termination::Completed);
//! assert_eq!(result.sampler, "chacha8-u53/rssa_cr/inverse/2");
//! // A prepared request is immutable: running it again replays the same streams.
//! let again = ssa::execute(&prepared, &Control::default()).unwrap();
//! assert_eq!(result.replicates[3].events, again.replicates[3].events);
//! ```
mod bounds;
mod selectors;
mod ziggurat;
pub use ziggurat::exponential_from_words;
mod ziggurat_tables;

use crate::{
    contracts,
    diagnostic::Diagnostic,
    sampling,
    simulation::{
        self, Assessment, Control, PathResult, Poller, Replicate, Statistics, Stop, Trajectory,
    },
};
use rand_core::RngCore;
use schemars::JsonSchema;
use selectors::{Composition, RateTree};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, path::Path, sync::Arc};

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of an alternative-SSA request; every version is still accepted.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum RequestVersion {
    /// The first version; `method.fast` is not allowed and the checked policy runs.
    #[serde(rename = "crn-craft.ssa-run/0.1")]
    V01,
    /// Adds the optional `method.fast` switch, default false.
    #[serde(rename = "crn-craft.ssa-run/0.2")]
    V02,
    /// Same fields as 0.2; `method.fast` defaults to true. Set it to false for
    /// the checked exact policy.
    #[serde(rename = "crn-craft.ssa-run/0.3")]
    V03,
    /// Same fields and defaults as 0.3, plus the `next_reaction` aggregator.
    /// Next reaction supports only the checked policy with inverse waiting
    /// times, so its `method.fast` defaults to false and must not be true.
    #[serde(rename = "crn-craft.ssa-run/0.4")]
    V04,
}
/// Schema tag of an alternative-SSA result; it matches the request version.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    /// Answers an ssa-run/0.1 request; `method` has no `fast` field.
    #[serde(rename = "crn-craft.ssa-result/0.1")]
    V01,
    /// Answers an ssa-run/0.2 request; `method` may carry `fast: true`.
    #[serde(rename = "crn-craft.ssa-result/0.2")]
    V02,
    /// Answers an ssa-run/0.3 request; `method.fast` is always present.
    #[serde(rename = "crn-craft.ssa-result/0.3")]
    V03,
    /// Answers an ssa-run/0.4 request; `method.fast` is always present.
    #[serde(rename = "crn-craft.ssa-result/0.4")]
    V04,
}
/// Exact event-selection algorithm. All sample the same law; they differ in cost and random-word use.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    #[doc = concat!("`direct`: full ordered propensity sum and binary search over prefix sums.")]
    Direct,
    #[doc = concat!("`cached_direct`: recomputes only dependent propensities and resumes the prefix")]
    #[doc = concat!("sums at the first changed channel.")]
    CachedDirect,
    #[doc = concat!("`direct_cr`: dependency updates, balanced totals and power-of-two")]
    #[doc = concat!("composition–rejection.")]
    DirectCr,
    #[doc = concat!("`sorting_direct`: dependency updates and a linear search whose order moves a")]
    #[doc = concat!("fired channel one place forward.")]
    SortingDirect,
    #[doc = concat!("`rssa_cr`: propensity bounds over a population bracket, composition–rejection")]
    #[doc = concat!("and thinning; bounds refresh only when the state leaves its bracket.")]
    RssaCr,
    #[doc = concat!("`next_reaction`: an indexed heap of absolute per-channel clocks; ssa-run/0.4")]
    #[doc = concat!("only, inverse waits and checked policy only.")]
    NextReaction,
}
impl Algorithm {
    fn name(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::CachedDirect => "cached_direct",
            Self::DirectCr => "direct_cr",
            Self::SortingDirect => "sorting_direct",
            Self::RssaCr => "rssa_cr",
            Self::NextReaction => "next_reaction",
        }
    }
}
/// Exponential waiting-time sampler.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WaitingTime {
    #[doc = concat!("`inverse`: −ln(U)/total with an open uniform U, as in the legacy Direct")]
    #[doc = concat!("sampler.")]
    Inverse,
    #[doc = concat!("`ziggurat`: the 256-layer exponential Ziggurat adapted from Julia's Random;")]
    #[doc = concat!("it consumes the stream differently.")]
    Ziggurat,
}
/// Method of an alternative-SSA request, selected by `kind`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    /// `ssa`: an exact stochastic simulation algorithm with a waiting-time sampler and policy.
    Ssa {
        #[doc = concat!("Event-selection algorithm.")]
        aggregator: Algorithm,
        #[doc = concat!("Exponential sampler; `next_reaction` requires `inverse`.")]
        waiting_time: WaitingTime,
        /// Not allowed in ssa-run/0.1; default false in 0.2 and true in 0.3
        /// and 0.4, except for next reaction, which is checked only.
        /// Emulates JumpProcesses.jl: Xoshiro256++ words, running totals
        /// updated by `+= new - old` without loss or 53-bit resolution
        /// checks, and stop polls per batch of events instead of per random
        /// word. Selection stays exact for the running totals; paths differ
        /// from the checked policy and replay only with the same setting.
        #[serde(
            default,
            deserialize_with = "present_bool",
            skip_serializing_if = "Option::is_none"
        )]
        #[schemars(with = "bool")]
        fast: Option<bool>,
    },
}
/// A present `fast` must be a boolean; `null` is rejected, not read as absent.
fn present_bool<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<bool>, D::Error> {
    bool::deserialize(d).map(Some)
}
impl Method {
    fn algorithm(self) -> Algorithm {
        let Self::Ssa { aggregator, .. } = self;
        aggregator
    }
    /// The policy in force: an explicit value, else the version's default.
    fn fast_for(self, version: RequestVersion) -> bool {
        let Self::Ssa {
            fast, aggregator, ..
        } = self;
        fast.unwrap_or(match version {
            RequestVersion::V01 | RequestVersion::V02 => false,
            RequestVersion::V03 => true,
            RequestVersion::V04 => aggregator != Algorithm::NextReaction,
        })
    }
    /// The resolved method as results echo it: 0.1/0.2 results carry `fast`
    /// only when true (unchanged bytes); 0.3 and 0.4 results always carry it.
    fn resolved(self, version: RequestVersion) -> Self {
        let fast = self.fast_for(version);
        let Self::Ssa {
            aggregator,
            waiting_time,
            ..
        } = self;
        Self::Ssa {
            aggregator,
            waiting_time,
            fast: if matches!(version, RequestVersion::V03 | RequestVersion::V04) {
                Some(fast)
            } else {
                fast.then_some(true)
            },
        }
    }
    fn fast(self) -> bool {
        let Self::Ssa { fast, .. } = self;
        fast.unwrap_or(false)
    }
}
/// Resource limits of an alternative-SSA request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    /// Most events per replicate, ≥ 1; a replicate that reaches it stops with `event_limit`.
    pub max_events_per_replicate: u64,
    /// Includes every random word drawn by every rejection/tail loop.
    /// A replicate that exhausts it stops with `sampling_limit` and `CRN-E067`.
    pub max_sampling_steps_per_replicate: u64,
    /// Optional wall-clock limit in seconds, finite and positive, polled cooperatively.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    /// Worker allowance, 1–1024; at most 64 workers and no more than the replicates are used.
    #[schemars(range(min = 1, max = 1024))]
    pub threads: u32,
    /// Byte budget of the compact result JSON with its trailing newline, also capped at 256 MiB.
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}
/// Alternative exact SSA request, bound to one model file by relative path and SHA-256.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    #[doc = concat!("Request version, `crn-craft.ssa-run/0.1` to `0.4`.")]
    pub schema: RequestVersion,
    #[doc = concat!("The model file this request belongs to.")]
    pub model: contracts::ModelRef,
    #[doc = concat!("Algorithm, waiting-time sampler and policy.")]
    pub method: Method,
    #[doc = concat!("Horizon and save grid, in seconds.")]
    pub time: contracts::TimeGrid,
    #[doc = concat!("Replicate count and root seed.")]
    pub ensemble: contracts::Ensemble,
    #[doc = concat!("Statistics and retained trajectories to report.")]
    pub outputs: contracts::Outputs,
    #[doc = concat!("Event, sampling, time, thread and output limits.")]
    pub budget: Budget,
}
impl Request {
    pub(crate) fn legacy(&self) -> contracts::Run {
        contracts::Run {
            schema: contracts::RunVersion::V01,
            model: self.model.clone(),
            method: contracts::Method::Ssa {
                aggregator: contracts::Aggregator::Direct,
            },
            time: self.time.clone(),
            ensemble: Some(self.ensemble.clone()),
            outputs: self.outputs.clone(),
            budget: contracts::Budget {
                max_events_per_replicate: Some(self.budget.max_events_per_replicate),
                wall_time_s: self.budget.wall_time_s,
                threads: self.budget.threads,
                max_output_bytes: self.budget.max_output_bytes,
            },
        }
    }
}
/// Work counters of a replicate or a whole run; for completed work they depend only on the
/// inputs and seed.
#[derive(Debug, Default, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Work {
    /// 64-bit random words drawn, the quantity `max_sampling_steps_per_replicate` limits.
    pub rng_words: u64,
    /// Channel selection attempts, including those composition–rejection rejected.
    pub proposals: u64,
    /// Propensity evaluations by the scheduler.
    pub propensity_evaluations: u64,
    /// `rssa_cr` propensity-bound refreshes.
    pub bound_updates: u64,
    /// Times an `rssa_cr` bound fell back to the propensity at the current state.
    pub point_fallbacks: u64,
    /// `rssa_cr` thinning rejections: the clock advanced but no reaction fired.
    pub rejected_candidates: u64,
}
impl Work {
    fn add(&mut self, rhs: &Self) {
        self.rng_words += rhs.rng_words;
        self.proposals += rhs.proposals;
        self.propensity_evaluations += rhs.propensity_evaluations;
        self.bound_updates += rhs.bound_updates;
        self.point_fallbacks += rhs.point_fallbacks;
        self.rejected_candidates += rhs.rejected_candidates;
    }
}
/// How a replicate, or the whole run, ended.
///
/// The declaration order is the aggregation priority: a run reports the
/// largest replicate value, or `output_limit` when retained output was removed.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    #[doc = concat!("`completed`: reached the horizon, or became absorbing and held its state.")]
    Completed,
    #[doc = concat!("`event_limit`: executed `max_events_per_replicate` events.")]
    EventLimit,
    #[doc = concat!("`sampling_limit`: drew `max_sampling_steps_per_replicate` random words")]
    #[doc = concat!("(`CRN-E067`).")]
    SamplingLimit,
    #[doc = concat!("`population_limit`: the next event would have exceeded a population limit and")]
    #[doc = concat!("was not applied.")]
    PopulationLimit,
    #[doc = concat!("`wall_time_limit`: `wall_time_s` elapsed.")]
    WallTimeLimit,
    #[doc = concat!("`cancelled`: cancellation was requested.")]
    Cancelled,
    #[doc = concat!("`failed`: a numerical or state-update check failed.")]
    Failed,
    #[doc = concat!("`output_limit`: the result did not fit `max_output_bytes`, so trajectories and")]
    #[doc = concat!("statistics were removed (`CRN-E064`). Run level only.")]
    OutputLimit,
}
fn termination(stop: Stop, diagnostic: Option<&Diagnostic>) -> Termination {
    match stop {
        Stop::Completed => Termination::Completed,
        Stop::EventLimit if diagnostic.is_some_and(|d| d.code == "CRN-E067") => {
            Termination::SamplingLimit
        }
        Stop::EventLimit => Termination::EventLimit,
        Stop::PopulationLimit => Termination::PopulationLimit,
        Stop::WallTimeLimit => Termination::WallTimeLimit,
        Stop::Cancelled => Termination::Cancelled,
        Stop::Failed => Termination::Failed,
        Stop::OutputLimit => Termination::OutputLimit,
    }
}
/// Accounting of one replicate; every requested id has exactly one record.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Record {
    /// Replicate id, which also selects its random stream.
    pub id: u64,
    /// False for a replicate that a stop prevented from starting.
    pub started: bool,
    #[doc = concat!("How the replicate ended.")]
    pub termination: Termination,
    /// Events executed; rejected candidates and refused updates are not counted.
    pub events: u64,
    /// The latest certified time of the path in seconds; the horizon when completed.
    pub stop_time: f64,
    /// Save times for which a state was recorded, a prefix of the save grid.
    pub saved_points: usize,
    /// The diagnostic behind a `failed`, `population_limit` or `sampling_limit` stop.
    pub diagnostic: Option<Diagnostic>,
    #[doc = concat!("Work counters of this replicate.")]
    pub work: Work,
}
/// Alternative exact SSA result: method, identities, replicate accounting, work, trajectories
/// and statistics.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SsaResult {
    #[doc = concat!("Matches the request version.")]
    pub schema: ResultVersion,
    #[doc = concat!("The method that ran. From 0.3 `fast` is always present; in 0.1 and 0.2 it")]
    #[doc = concat!("appears only when true.")]
    pub method: Method,
    /// SHA-256 of the exact model bytes.
    pub model_sha256: String,
    /// SHA-256 of the exact request bytes.
    pub run_sha256: String,
    /// Sampler identity: `chacha8-u53/ALGORITHM/WAIT/N` for the checked policy, with N = 2 for
    /// `direct_cr` and `rssa_cr` and 1 otherwise, and `xoshiro256pp-u53/ALGORITHM/WAIT/fast-1`
    /// for `fast`; WAIT is `inverse` or `ziggurat256-julia`. Replay needs the same identity,
    /// seed and input bytes.
    pub sampler: String,
    /// Root seed from the request.
    pub root_seed: u64,
    /// Species names in the order used by every count row.
    pub species: Vec<String>,
    /// The request's save grid in seconds.
    pub save_times: Vec<f64>,
    #[doc = concat!("Aggregate stop: the highest-priority replicate termination, or `output_limit`.")]
    pub termination: Termination,
    #[doc = concat!("What the result may be used for; `exact_ensemble` only when every replicate")]
    #[doc = concat!("completed.")]
    pub assessment: Assessment,
    /// Replicates the request asked for.
    pub requested_replicates: u64,
    /// Replicates that completed.
    pub completed_replicates: u64,
    /// Events executed over all replicates.
    pub events: u64,
    #[doc = concat!("Work counters summed over all replicates.")]
    pub work: Work,
    /// One record per requested replicate, in id order.
    pub replicates: Vec<Record>,
    /// Saved states of the first `keep_trajectories` replicates; empty when output was removed.
    pub trajectories: Vec<Trajectory>,
    /// Requested statistics, present only when every replicate completed and output fit.
    pub statistics: Option<Statistics>,
    /// Model warnings, plus `CRN-E064` when output was removed.
    pub diagnostics: Vec<Diagnostic>,
}

/// A validated, admitted alternative-SSA execution for [`execute`] or
/// [`execute_with_kernel`]. Immutable: each execution restarts the same
/// streams, and concurrent executions own separate state.
pub struct Prepared {
    pub(crate) core: simulation::Prepared,
    pub(crate) request: Request,
    pub(crate) readers: Arc<Vec<Vec<usize>>>,
    bound_groups: Vec<usize>,
}

/// Validates a request against the exact model bytes and admits it.
///
/// `model.path` is only syntax-checked; use [`solve_files`] to check that it
/// names the model file.
///
/// # Errors
///
/// Strict-parse errors of the request; `CRN-E020` for `fast` in 0.1 or an
/// invalid `next_reaction` combination (at the offending field); every
/// [`crate::check_bytes`] error of the model and the run fields; and
/// `CRN-E060` when the workspace or work accounting exceeds its limits.
pub fn prepare(model: &[u8], request_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let mut request: Request =
        crate::parse_network(request_bytes, "`crn-craft.ssa-run` requests").map_err(|e| vec![e])?;
    if request.method.fast_for(request.schema) && request.schema == RequestVersion::V01 {
        return Err(vec![Diagnostic::error(
            "CRN-E020",
            "/method/fast",
            "fast requires crn-craft.ssa-run/0.2 or later",
            "Use crn-craft.ssa-run/0.3 (fast by default), or remove fast for the checked exact policy.",
        )]);
    }
    if let Some(error) = next_reaction_error(&request) {
        return Err(vec![error]);
    }
    // From here on `fast` is explicit wherever it applies.
    request.method = request.method.resolved(request.schema);
    let legacy = serde_json::to_vec(&request.legacy())
        .map_err(|e| vec![simulation::execution_error(e.to_string())])?;
    let core = simulation::prepare(model, &legacy)?;
    finish_prepare(core, request, request_bytes, model.len(), None)
}

/// Next reaction is new in ssa-run/0.4 and runs the checked policy with
/// inverse waiting times only.
fn next_reaction_error(request: &Request) -> Option<Diagnostic> {
    let Method::Ssa { waiting_time, .. } = request.method;
    if request.method.algorithm() != Algorithm::NextReaction {
        return None;
    }
    let (path, message, hint) = if request.schema != RequestVersion::V04 {
        (
            "/method/aggregator",
            "next_reaction requires crn-craft.ssa-run/0.4",
            "Select crn-craft.ssa-run/0.4 or another aggregator.",
        )
    } else if request.method.fast_for(request.schema) {
        (
            "/method/fast",
            "next_reaction supports only the checked policy",
            "Remove method.fast or set it to false.",
        )
    } else if waiting_time != WaitingTime::Inverse {
        (
            "/method/waiting_time",
            "next_reaction requires inverse waiting times",
            "Set method.waiting_time to inverse.",
        )
    } else {
        return None;
    };
    Some(Diagnostic::error("CRN-E020", path, message, hint))
}

pub(crate) fn finish_prepare(
    mut core: simulation::Prepared,
    request: Request,
    request_bytes: &[u8],
    model_len: usize,
    shared_readers: Option<Arc<Vec<Vec<usize>>>>,
) -> Result<Prepared, Vec<Diagnostic>> {
    // Parameter binding reaches this point without `prepare`.
    if let Some(error) = next_reaction_error(&request) {
        return Err(vec![error]);
    }
    // Additional selector and bounded-range buffers, conservatively reserved.
    let reactions = core.network.reactions().len() as u64;
    let species = core.network.species().len() as u64;
    let extra =
        (core.workers as u128) * (reactions as u128 * 384 + species as u128 * 128 + 200_000);
    let reads: usize = core.network.reactions().iter().map(|r| r.reads.len()).sum();
    let workspace = core.workspace as u128 + extra + reads as u128 * 8 + species as u128 * 24;
    let max_work = (u128::from(request.budget.max_events_per_replicate) + 1)
        * (u128::from(reactions) * 260 + u128::from(species) + 1)
        * u128::from(request.ensemble.replicates);
    let sampling_work = (u128::from(request.budget.max_sampling_steps_per_replicate) + 1)
        * u128::from(request.ensemble.replicates);
    if workspace > u128::from(simulation::MAX_WORKSPACE_BYTES)
        || max_work + sampling_work > u128::from(u64::MAX)
    {
        return Err(vec![simulation::execution_error(
            "alternative SSA workspace or work accounting exceeds limits",
        )]);
    }
    core.workspace = workspace as u64;
    core.run_hash = crate::sha256(request_bytes);
    core.input_bytes = (model_len + request_bytes.len()) as u64;
    let readers = shared_readers.unwrap_or_else(|| {
        let mut readers = vec![Vec::new(); species as usize];
        for (rx, structure) in core.network.reactions().iter().enumerate() {
            for &s in &structure.reads {
                readers[s].push(rx);
            }
        }
        Arc::new(readers)
    });
    let mut group_keys = BTreeMap::new();
    let mut bound_groups = vec![usize::MAX; reactions as usize];
    for (rx, slot) in bound_groups.iter_mut().enumerate() {
        if let Some(key) = core.network.compiled_rate(rx).bounds_key() {
            let next = group_keys.len();
            *slot = *group_keys.entry(key).or_insert(next);
        }
    }
    Ok(Prepared {
        core,
        request,
        readers,
        bound_groups,
    })
}

fn numeric(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E061",
        "/simulation",
        message,
        "Rescale the rates; the exact sampler does not silently discard probability.",
    )
}
enum Bail {
    Stopped(Stop),
    SamplingLimit,
    Numerical(Diagnostic),
}
impl From<Diagnostic> for Bail {
    fn from(value: Diagnostic) -> Self {
        Self::Numerical(value)
    }
}
struct Draws<'a, 'b> {
    limit: u64,
    work: &'a mut Work,
    poller: &'a mut Poller<'b>,
    /// `fast` polls for stops per batch of event-loop iterations and on
    /// rejection retries, not before every random word.
    fast: bool,
}
impl Draws<'_, '_> {
    fn poll(&mut self) -> Result<(), Bail> {
        if let Some(s) = self.poller.poll() {
            Err(Bail::Stopped(s))
        } else {
            Ok(())
        }
    }
    /// Before redrawing in a rejection loop (exact runs already polled per word).
    fn retry(&mut self) -> Result<(), Bail> {
        if self.fast { self.poll() } else { Ok(()) }
    }
    fn word(&mut self, rng: &mut impl RngCore) -> Result<u64, Bail> {
        if !self.fast {
            self.poll()?;
        }
        if self.work.rng_words >= self.limit {
            return Err(Bail::SamplingLimit);
        }
        self.work.rng_words += 1;
        Ok(rng.next_u64())
    }
    fn closed(&mut self, rng: &mut impl RngCore) -> Result<f64, Bail> {
        Ok(sampling::uniform_closed_open(self.word(rng)?))
    }
    fn open(&mut self, rng: &mut impl RngCore) -> Result<f64, Bail> {
        loop {
            let u = self.closed(rng)?;
            if u > 0.0 {
                return Ok(u);
            }
            self.retry()?;
        }
    }
    fn index(&mut self, rng: &mut impl RngCore, n: usize) -> Result<usize, Bail> {
        let n = n as u64;
        let cutoff = n.wrapping_neg() % n;
        loop {
            let product = u128::from(self.word(rng)?) * u128::from(n);
            if product as u64 >= cutoff {
                return Ok((product >> 64) as usize);
            }
            self.retry()?;
        }
    }
}

/// Model-specific propensity evaluators must match the checked reference in
/// every valid state, including gates, rounded arithmetic and diagnostics.
/// Generated kernels are guarded by the fixed-topology fingerprint on execute.
pub trait PropensityKernel: Sync {
    /// The [`crate::CompiledNetwork::topology_sha256`] of the model the kernel was written for.
    fn topology_sha256(&self) -> &str;
    /// The propensity of `reaction` at `state`, identical to
    /// [`crate::CompiledNetwork::propensity_in_scheduler_state`], including its
    /// errors.
    fn evaluate(
        &self,
        network: &crate::CompiledNetwork,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic>;
}
struct DynamicKernel;
impl PropensityKernel for DynamicKernel {
    fn topology_sha256(&self) -> &str {
        ""
    }
    #[inline]
    fn evaluate(
        &self,
        network: &crate::CompiledNetwork,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic> {
        network.propensity_in_valid_state(reaction, state)
    }
}

#[derive(Default)]
struct Workspace {
    state: Vec<i64>,
    rates: Vec<f64>,
    prefix: Vec<f64>,
    minima: Vec<f64>,
    tree: RateTree,
    composition: Composition,
    order: Vec<usize>,
    bounds: bounds::Bounds,
}

fn path<const MODE: u8, const ZIG: bool, const FAST: bool, K: PropensityKernel>(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    space: &mut Workspace,
    kernel: &K,
) -> (PathResult, Work) {
    let seed = prepared.core.run.ensemble.as_ref().unwrap().root_seed;
    if FAST {
        run_path::<MODE, ZIG, FAST, _, _>(
            prepared,
            id,
            poller,
            space,
            sampling::fast_stream(seed, id),
            kernel,
        )
    } else {
        run_path::<MODE, ZIG, FAST, _, _>(
            prepared,
            id,
            poller,
            space,
            sampling::stream(seed, id),
            kernel,
        )
    }
}

/// Stop polls in `fast` runs: one per this many event-loop iterations.
const FAST_POLL_STRIDE: u64 = 256;

fn run_path<const MODE: u8, const ZIG: bool, const FAST: bool, R: RngCore, K: PropensityKernel>(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    space: &mut Workspace,
    mut rng: R,
    kernel: &K,
) -> (PathResult, Work) {
    let network = &prepared.core.network;
    let run = &prepared.core.run;
    let mut record = Replicate {
        id,
        started: false,
        termination: Stop::Completed,
        events: 0,
        stop_time: 0.0,
        saved_points: 0,
        diagnostic: None,
    };
    let mut work = Work::default();
    let mut draws = Draws {
        limit: prepared.request.budget.max_sampling_steps_per_replicate,
        work: &mut work,
        poller,
        fast: FAST,
    };
    let state = &mut space.state;
    state.clear();
    state.extend_from_slice(network.initial_counts());
    let n = network.reactions().len();
    space.rates.resize(n, 0.0);
    space.prefix.resize(n, 0.0);
    space.minima.resize(n, 0.0);
    if MODE >= 2 && MODE != 5 {
        space.tree.reset_fast(n, FAST);
        if MODE == 2 || MODE == 4 {
            space.composition.reset(n);
        }
        if MODE == 3 {
            space.order.clear();
            space.order.extend(0..n);
        }
    }
    if MODE == 4 {
        space.bounds.reset(network, state);
    }
    let mut counts = Vec::with_capacity(run.time.save_times.len() * state.len());
    let mut time = 0.0;
    let mut save = 0;
    let mut previous: Option<usize> = None;
    let mut total = 0.0;
    let mut bounds_initialized = false;
    let mut absorption_checked = false;
    let mut absorbed = false;
    let mut iterations: u64 = 0;
    let mut clocks = crate::next_reaction::Clocks::default();
    if MODE == 5 {
        clocks.reset(n);
    }
    let outcome = (|| -> Result<Stop, Bail> {
        draws.poll()?;
        record.started = true;
        loop {
            while save < run.time.save_times.len() && run.time.save_times[save] <= time {
                draws.poll()?;
                counts.extend_from_slice(state);
                save += 1;
            }
            if time == run.time.end {
                return Ok(Stop::Completed);
            }
            if !FAST || iterations % FAST_POLL_STRIDE == 0 {
                draws.poll()?;
            }
            iterations += 1;
            if record.events >= run.budget.max_events_per_replicate.unwrap() {
                return Ok(Stop::EventLimit);
            }
            let (next, channel) = if MODE == 5 {
                let mut update = |rx, fired| -> Result<(), Bail> {
                    draws.poll()?;
                    let a = kernel.evaluate(network, rx, state)?;
                    draws.work.propensity_evaluations += 1;
                    let fresh = if clocks.needs_draw(rx, a, fired) {
                        Some(-libm::log(draws.open(&mut rng)?))
                    } else {
                        None
                    };
                    if !clocks.update(rx, a, time, fresh) {
                        return Err(numeric("next reaction clock failed to advance").into());
                    }
                    Ok(())
                };
                if let Some(rx) = previous {
                    update(rx, true)?;
                    for &other in &network.dependencies()[rx] {
                        if other != rx {
                            update(other, false)?;
                        }
                    }
                } else {
                    for rx in 0..n {
                        update(rx, false)?;
                    }
                }
                let candidate = clocks.next();
                let until = candidate.map_or(f64::INFINITY, |(t, _)| t);
                while save < run.time.save_times.len() && run.time.save_times[save] < until {
                    draws.poll()?;
                    counts.extend_from_slice(state);
                    save += 1;
                }
                if until > run.time.end {
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                candidate.expect("finite clock")
            } else {
                if MODE == 4 {
                    if !bounds_initialized || previous.is_some() {
                        space.bounds.refresh(
                            network,
                            &prepared.readers,
                            &prepared.bound_groups,
                            state,
                            previous.take(),
                            &mut space.tree,
                            &mut space.composition,
                            &mut draws,
                        )?;
                        bounds_initialized = true;
                    }
                    if !absorption_checked && space.bounds.positive_lower == 0 {
                        absorbed = true;
                        for rx in 0..n {
                            draws.poll()?;
                            draws.work.propensity_evaluations += 1;
                            if kernel.evaluate(network, rx, state)? > 0.0 {
                                absorbed = false;
                                break;
                            }
                        }
                    }
                    absorption_checked = true;
                } else if MODE > 0 && previous.is_some() {
                    let affected = &network.dependencies()[previous.unwrap()];
                    for (k, &rx) in affected.iter().enumerate() {
                        if k % 64 == 0 && (!FAST || k > 0) {
                            draws.poll()?;
                        }
                        let a = kernel.evaluate(network, rx, state)?;
                        draws.work.propensity_evaluations += 1;
                        space.rates[rx] = a;
                        if MODE >= 2 {
                            space
                                .tree
                                .set(rx, a, if a > 0.0 { a } else { f64::INFINITY });
                            if MODE == 2 {
                                space.composition.update(rx, a)?;
                            }
                        }
                    }
                    if MODE == 1 {
                        if let Some(&first) = affected.first() {
                            total = if FAST {
                                sampling::fast_prefix(&space.rates, &mut space.prefix, first)?
                            } else {
                                sampling::resume_total(
                                    &space.rates,
                                    &mut space.prefix,
                                    &mut space.minima,
                                    first,
                                )?
                            };
                        }
                    }
                } else {
                    for rx in 0..n {
                        if rx % 64 == 0 {
                            draws.poll()?;
                        }
                        let a = kernel.evaluate(network, rx, state)?;
                        draws.work.propensity_evaluations += 1;
                        space.rates[rx] = a;
                        if MODE >= 2 {
                            space
                                .tree
                                .set(rx, a, if a > 0.0 { a } else { f64::INFINITY });
                            if MODE == 2 {
                                space.composition.update(rx, a)?;
                            }
                        }
                    }
                    if FAST && MODE <= 1 {
                        total = sampling::fast_prefix(&space.rates, &mut space.prefix, 0)?;
                    } else if MODE == 0 {
                        total = sampling::total_with_prefix(&space.rates, &mut space.prefix)?;
                    } else if MODE == 1 {
                        total = sampling::resume_total(
                            &space.rates,
                            &mut space.prefix,
                            &mut space.minima,
                            0,
                        )?;
                    }
                }
                if MODE >= 2 {
                    total = if MODE == 4 && absorbed {
                        0.0
                    } else {
                        space.tree.total()?
                    };
                }
                if total == 0.0 {
                    while save < run.time.save_times.len() {
                        draws.poll()?;
                        counts.extend_from_slice(state);
                        save += 1;
                    }
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                let wait = if ZIG {
                    ziggurat::wait(total, &mut rng, &mut draws)?
                } else {
                    sampling::waiting_time(total, draws.open(&mut rng)?)?
                };
                let next = sampling::event_time(time, wait)?;
                while save < run.time.save_times.len() && run.time.save_times[save] < next {
                    draws.poll()?;
                    counts.extend_from_slice(state);
                    save += 1;
                }
                if next > run.time.end {
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                let channel = if MODE == 2 || MODE == 4 {
                    space
                        .composition
                        .sample(&space.tree, &mut rng, &mut draws)?
                } else if MODE == 3 {
                    draws.work.proposals += 1;
                    selectors::sorting_sample(
                        &mut space.order,
                        &space.rates,
                        total,
                        &mut rng,
                        &mut draws,
                    )?
                } else {
                    draws.work.proposals += 1;
                    sampling::select_from_prefix(&space.prefix, total, draws.closed(&mut rng)?)?
                };
                (next, channel)
            };
            if MODE == 4 {
                let upper = space.tree.get(channel);
                let threshold = draws.closed(&mut rng)? * upper;
                if threshold >= upper {
                    return Err(numeric("RSSA acceptance rounded to its bound").into());
                }
                if threshold >= space.bounds.lower[channel] {
                    let a = kernel.evaluate(network, channel, state)?;
                    draws.work.propensity_evaluations += 1;
                    if a < space.bounds.lower[channel] || a > upper {
                        return Err(numeric(
                            "RSSA propensity escaped its certified integer-box bounds",
                        )
                        .into());
                    }
                    if threshold >= a {
                        time = next;
                        draws.work.rejected_candidates += 1;
                        continue;
                    }
                }
            }
            if let Err(error) = network.fire_enabled_reaction(channel, state) {
                if error.code == "CRN-E044" && error.message.contains("safety limit") {
                    record.diagnostic = Some(error);
                    return Ok(Stop::PopulationLimit);
                }
                return Err(error.into());
            }
            previous = Some(channel);
            absorption_checked = false;
            time = next;
            record.events += 1;
        }
    })();
    match outcome {
        Ok(stop) | Err(Bail::Stopped(stop)) => record.termination = stop,
        Err(Bail::SamplingLimit) => {
            record.termination = Stop::EventLimit;
            record.diagnostic = Some(Diagnostic::error(
                "CRN-E067",
                "/budget/max_sampling_steps_per_replicate",
                "sampling work limit reached",
                "Increase the explicit sampling-step budget; incomplete ensembles withhold statistics.",
            ));
        }
        Err(Bail::Numerical(error)) => {
            record.termination = Stop::Failed;
            record.diagnostic = Some(error);
        }
    }
    record.stop_time = time.max(save.checked_sub(1).map_or(0.0, |i| run.time.save_times[i]));
    record.saved_points = save;
    (PathResult { record, counts }, work)
}

/// Runs a prepared request with the built-in propensity evaluator.
///
/// # Errors
///
/// `CRN-E060` when a worker cannot start or even the replicate accounting
/// exceeds `max_output_bytes`. Stops and numerical failures are reported
/// per replicate.
pub fn execute(prepared: &Prepared, control: &Control) -> Result<SsaResult, Diagnostic> {
    execute_inner(prepared, control, &DynamicKernel)
}

/// Execute an optional trusted Rust kernel against its exact topology. As with
/// custom Rust simulation code, callers must establish evaluator correctness;
/// this fingerprint prevents accidental use with a different network.
///
/// # Errors
///
/// `CRN-E020` at `/kernel` when the kernel's fingerprint differs from the
/// prepared network's, otherwise as [`execute`].
pub fn execute_with_kernel<K: PropensityKernel>(
    prepared: &Prepared,
    control: &Control,
    kernel: &K,
) -> Result<SsaResult, Diagnostic> {
    if kernel.topology_sha256() != prepared.core.network.topology_sha256() {
        return Err(Diagnostic::error(
            "CRN-E020",
            "/kernel",
            "generated kernel topology mismatch",
            "Regenerate the kernel from this parameter template's model.",
        ));
    }
    execute_inner(prepared, control, kernel)
}

fn execute_inner<K: PropensityKernel>(
    prepared: &Prepared,
    control: &Control,
    kernel: &K,
) -> Result<SsaResult, Diagnostic> {
    let algorithm = prepared.request.method.algorithm();
    let Method::Ssa { waiting_time, .. } = prepared.request.method;
    // Resolved at preparation: explicit for 0.3, present only when true before.
    let fast = prepared.request.method.fast();
    let wait = match waiting_time {
        WaitingTime::Inverse => "inverse",
        WaitingTime::Ziggurat => "ziggurat256-julia",
    };
    let sampler = if fast {
        format!("xoshiro256pp-u53/{}/{wait}/fast-1", algorithm.name())
    } else {
        format!(
            "chacha8-u53/{}/{wait}/{}",
            algorithm.name(),
            // Version 2: composition draws member and acceptance from one word.
            if matches!(algorithm, Algorithm::DirectCr | Algorithm::RssaCr) {
                2
            } else {
                1
            }
        )
    };
    macro_rules! run {
        ($mode:literal, $zig:literal, $id:expr, $poller:expr, $work:expr) => {
            if fast {
                path::<$mode, $zig, true, _>(prepared, $id, $poller, $work, kernel)
            } else {
                path::<$mode, $zig, false, _>(prepared, $id, $poller, $work, kernel)
            }
        };
    }
    let (simulation, metrics) = simulation::execute_paths::<Workspace, Work, _>(
        &prepared.core,
        control,
        &sampler,
        None,
        |id, poller, work| match (algorithm, waiting_time) {
            (Algorithm::Direct, WaitingTime::Inverse) => run!(0, false, id, poller, work),
            (Algorithm::CachedDirect, WaitingTime::Inverse) => run!(1, false, id, poller, work),
            (Algorithm::DirectCr, WaitingTime::Inverse) => run!(2, false, id, poller, work),
            (Algorithm::SortingDirect, WaitingTime::Inverse) => run!(3, false, id, poller, work),
            (Algorithm::RssaCr, WaitingTime::Inverse) => run!(4, false, id, poller, work),
            (Algorithm::Direct, WaitingTime::Ziggurat) => run!(0, true, id, poller, work),
            (Algorithm::CachedDirect, WaitingTime::Ziggurat) => run!(1, true, id, poller, work),
            (Algorithm::DirectCr, WaitingTime::Ziggurat) => run!(2, true, id, poller, work),
            (Algorithm::SortingDirect, WaitingTime::Ziggurat) => run!(3, true, id, poller, work),
            (Algorithm::RssaCr, WaitingTime::Ziggurat) => run!(4, true, id, poller, work),
            (Algorithm::NextReaction, _) => {
                path::<5, false, false, _>(prepared, id, poller, work, kernel)
            }
        },
    )?;
    let r = simulation.result;
    let mut work = Work::default();
    let replicates: Vec<_> = r
        .replicates
        .into_iter()
        .zip(metrics)
        .map(|(r, w)| {
            work.add(&w);
            Record {
                id: r.id,
                started: r.started,
                termination: termination(r.termination, r.diagnostic.as_ref()),
                events: r.events,
                stop_time: r.stop_time,
                saved_points: r.saved_points,
                diagnostic: r.diagnostic,
                work: w,
            }
        })
        .collect();
    let mut stop = replicates
        .iter()
        .map(|r| r.termination)
        .max()
        .unwrap_or(Termination::Failed);
    if r.termination == Stop::OutputLimit {
        stop = Termination::OutputLimit;
    }
    let mut result = SsaResult {
        schema: match prepared.request.schema {
            RequestVersion::V01 => ResultVersion::V01,
            RequestVersion::V02 => ResultVersion::V02,
            RequestVersion::V03 => ResultVersion::V03,
            RequestVersion::V04 => ResultVersion::V04,
        },
        method: prepared.request.method,
        model_sha256: r.model_sha256,
        run_sha256: r.run_sha256,
        sampler: r.sampler,
        root_seed: r.root_seed,
        species: r.species,
        save_times: r.save_times,
        termination: stop,
        assessment: r.assessment,
        requested_replicates: r.requested_replicates,
        completed_replicates: r.completed_replicates,
        events: r.events,
        work,
        replicates,
        trajectories: r.trajectories,
        statistics: r.statistics,
        diagnostics: r.diagnostics,
    };
    let limit = prepared
        .request
        .budget
        .max_output_bytes
        .min(simulation::MAX_RESULT_BYTES);
    // The budget covers this compact result JSON and its trailing newline.
    if !simulation::fits_json(&result, limit) {
        result.termination = Termination::OutputLimit;
        result.assessment = Assessment::Incomplete;
        result.statistics = None;
        result.trajectories.clear();
        result.diagnostics.push(Diagnostic::error(
            "CRN-E064",
            "/outputs",
            "alternative SSA output exceeds its byte budget",
            "Increase max_output_bytes; statistics and trajectories have been removed.",
        ));
        if !simulation::fits_json(&result, limit) {
            return Err(simulation::execution_error(
                "SSA evidence exceeds output budget",
            ));
        }
    }
    Ok(result)
}
/// [`prepare`] then [`execute`] on in-memory bytes.
///
/// # Errors
///
/// The errors of both steps.
pub fn solve(
    model: &[u8],
    request: &[u8],
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    execute(&prepare(model, request)?, control).map_err(|e| vec![e])
}
/// Reads both files, requires the request's `model.path`, resolved from the
/// request's directory, to name `model_path` inside that directory, then
/// calls [`solve`].
///
/// # Errors
///
/// `CRN-E030` or `CRN-E004` for unreadable or oversized files, strict-parse
/// errors, `CRN-E022` when the reference does not name the model file, and
/// the errors of [`solve`].
pub fn solve_files(
    model_path: &Path,
    request_path: &Path,
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    let model = crate::read_document(model_path).map_err(|e| vec![e])?;
    let bytes = crate::read_document(request_path).map_err(|e| vec![e])?;
    let request: Request =
        crate::parse_network(&bytes, "`crn-craft.ssa-run` requests").map_err(|e| vec![e])?;
    let parent = request_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let resolved = (
        parent.canonicalize(),
        parent.join(&request.model.path).canonicalize(),
        model_path.canonicalize(),
    );
    if !crate::validate::valid_relative_path(&request.model.path)
        || !matches!(resolved,(Ok(ref base),Ok(ref reference),Ok(ref supplied)) if reference.starts_with(base) && reference==supplied)
    {
        return Err(vec![Diagnostic::error(
            "CRN-E022",
            "/model/path",
            "SSA request does not identify the supplied model below its directory",
            "Use the referenced model within the request directory.",
        )]);
    }
    solve(&model, &bytes, control)
}
