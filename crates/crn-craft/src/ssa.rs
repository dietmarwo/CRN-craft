//! Versioned alternative exact SSA methods. Legacy Direct contracts are unchanged.
mod bounds;
mod selectors;
mod ziggurat;
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
use std::path::Path;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum RequestVersion {
    #[serde(rename = "crn-craft.ssa-run/0.1")]
    V01,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.ssa-result/0.1")]
    V01,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    Direct,
    CachedDirect,
    DirectCr,
    SortingDirect,
    RssaCr,
}
impl Algorithm {
    fn name(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::CachedDirect => "cached_direct",
            Self::DirectCr => "direct_cr",
            Self::SortingDirect => "sorting_direct",
            Self::RssaCr => "rssa_cr",
        }
    }
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WaitingTime {
    Inverse,
    Ziggurat,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    Ssa {
        aggregator: Algorithm,
        waiting_time: WaitingTime,
    },
}
impl Method {
    fn algorithm(self) -> Algorithm {
        let Self::Ssa { aggregator, .. } = self;
        aggregator
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    pub max_events_per_replicate: u64,
    /// Includes every random word drawn by every rejection/tail loop.
    pub max_sampling_steps_per_replicate: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    #[schemars(range(min = 1, max = 1024))]
    pub threads: u32,
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub schema: RequestVersion,
    pub model: contracts::ModelRef,
    pub method: Method,
    pub time: contracts::TimeGrid,
    pub ensemble: contracts::Ensemble,
    pub outputs: contracts::Outputs,
    pub budget: Budget,
}
impl Request {
    fn legacy(&self) -> contracts::Run {
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
#[derive(Debug, Default, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Work {
    pub rng_words: u64,
    pub proposals: u64,
    pub propensity_evaluations: u64,
    pub bound_updates: u64,
    pub point_fallbacks: u64,
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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    Completed,
    EventLimit,
    SamplingLimit,
    PopulationLimit,
    WallTimeLimit,
    Cancelled,
    Failed,
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
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub id: u64,
    pub started: bool,
    pub termination: Termination,
    pub events: u64,
    pub stop_time: f64,
    pub saved_points: usize,
    pub diagnostic: Option<Diagnostic>,
    pub work: Work,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SsaResult {
    pub schema: ResultVersion,
    pub method: Method,
    pub model_sha256: String,
    pub run_sha256: String,
    pub sampler: String,
    pub root_seed: u64,
    pub species: Vec<String>,
    pub save_times: Vec<f64>,
    pub termination: Termination,
    pub assessment: Assessment,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub events: u64,
    pub work: Work,
    pub replicates: Vec<Record>,
    pub trajectories: Vec<Trajectory>,
    pub statistics: Option<Statistics>,
    pub diagnostics: Vec<Diagnostic>,
}

pub struct Prepared {
    core: simulation::Prepared,
    request: Request,
    readers: Vec<Vec<usize>>,
}

pub fn prepare(model: &[u8], request_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let request: Request = crate::strict_json::parse(request_bytes).map_err(|e| vec![e])?;
    let legacy = serde_json::to_vec(&request.legacy())
        .map_err(|e| vec![simulation::execution_error(e.to_string())])?;
    let mut core = simulation::prepare(model, &legacy)?;
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
    core.input_bytes = (model.len() + request_bytes.len()) as u64;
    let mut readers = vec![Vec::new(); species as usize];
    for (rx, structure) in core.network.reactions().iter().enumerate() {
        for &s in &structure.reads {
            readers[s].push(rx);
        }
    }
    Ok(Prepared {
        core,
        request,
        readers,
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
}
impl Draws<'_, '_> {
    #[inline]
    fn poll(&mut self) -> Result<(), Bail> {
        if let Some(s) = self.poller.poll() {
            Err(Bail::Stopped(s))
        } else {
            Ok(())
        }
    }
    /// Callers poll once per event-loop iteration, per rejection retry and
    /// every 64 inner steps; a single word is bounded by the work limit only.
    #[inline]
    fn word(&mut self, rng: &mut impl RngCore) -> Result<u64, Bail> {
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
            self.poll()?;
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
            self.poll()?;
        }
    }
}

#[derive(Default)]
struct Workspace {
    state: Vec<i64>,
    rates: Vec<f64>,
    prefix: Vec<f64>,
    floor: f64,
    tree: RateTree,
    composition: Composition,
    order: Vec<usize>,
    bounds: bounds::Bounds,
}

fn path<const MODE: u8, const ZIG: bool>(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    space: &mut Workspace,
) -> (PathResult, Work) {
    let network = &prepared.core.network;
    let run = &prepared.core.run;
    let mut rng = sampling::stream(run.ensemble.as_ref().unwrap().root_seed, id);
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
    };
    let state = &mut space.state;
    state.clear();
    state.extend_from_slice(network.initial_counts());
    let n = network.reactions().len();
    space.rates.resize(n, 0.0);
    space.prefix.resize(n, 0.0);
    space.floor = f64::INFINITY;
    if MODE >= 2 {
        space.tree.reset(n);
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
            draws.poll()?;
            if record.events >= run.budget.max_events_per_replicate.unwrap() {
                return Ok(Stop::EventLimit);
            }
            if MODE == 4 {
                if !bounds_initialized || previous.is_some() {
                    space.bounds.refresh(
                        network,
                        &prepared.readers,
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
                        if network.propensity_in_valid_state(rx, state)? > 0.0 {
                            absorbed = false;
                            break;
                        }
                    }
                }
                absorption_checked = true;
            } else if MODE > 0 && previous.is_some() {
                let affected = &network.dependencies()[previous.unwrap()];
                for (k, &rx) in affected.iter().enumerate() {
                    if k % 64 == 0 {
                        draws.poll()?;
                    }
                    let a = network.propensity_in_valid_state(rx, state)?;
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
                        total = sampling::resume_total(
                            &space.rates,
                            &mut space.prefix,
                            &mut space.floor,
                            first,
                        )?;
                    }
                }
            } else {
                for rx in 0..n {
                    if rx % 64 == 0 {
                        draws.poll()?;
                    }
                    let a = network.propensity_in_valid_state(rx, state)?;
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
                if MODE == 0 {
                    total = sampling::total_with_prefix(&space.rates, &mut space.prefix)?;
                }
                if MODE == 1 {
                    total = sampling::resume_total(
                        &space.rates,
                        &mut space.prefix,
                        &mut space.floor,
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
            if MODE == 4 {
                let upper = space.tree.get(channel);
                let threshold = draws.closed(&mut rng)? * upper;
                if threshold >= upper {
                    return Err(numeric("RSSA acceptance rounded to its bound").into());
                }
                if threshold >= space.bounds.lower[channel] {
                    let a = network.propensity_in_valid_state(channel, state)?;
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

pub fn execute(prepared: &Prepared, control: &Control) -> Result<SsaResult, Diagnostic> {
    let algorithm = prepared.request.method.algorithm();
    let Method::Ssa { waiting_time, .. } = prepared.request.method;
    let sampler = format!(
        "chacha8-u53/{}/{}/{}",
        algorithm.name(),
        match waiting_time {
            WaitingTime::Inverse => "inverse",
            WaitingTime::Ziggurat => "ziggurat256-julia",
        },
        // Version 2: composition draws member and acceptance from one word.
        if matches!(algorithm, Algorithm::DirectCr | Algorithm::RssaCr) {
            2
        } else {
            1
        }
    );
    let (simulation, metrics) = simulation::execute_paths::<Workspace, Work, _>(
        &prepared.core,
        control,
        &sampler,
        None,
        |id, poller, work| match (algorithm, waiting_time) {
            (Algorithm::Direct, WaitingTime::Inverse) => {
                path::<0, false>(prepared, id, poller, work)
            }
            (Algorithm::CachedDirect, WaitingTime::Inverse) => {
                path::<1, false>(prepared, id, poller, work)
            }
            (Algorithm::DirectCr, WaitingTime::Inverse) => {
                path::<2, false>(prepared, id, poller, work)
            }
            (Algorithm::SortingDirect, WaitingTime::Inverse) => {
                path::<3, false>(prepared, id, poller, work)
            }
            (Algorithm::Direct, WaitingTime::Ziggurat) => {
                path::<0, true>(prepared, id, poller, work)
            }
            (Algorithm::CachedDirect, WaitingTime::Ziggurat) => {
                path::<1, true>(prepared, id, poller, work)
            }
            (Algorithm::DirectCr, WaitingTime::Ziggurat) => {
                path::<2, true>(prepared, id, poller, work)
            }
            (Algorithm::SortingDirect, WaitingTime::Ziggurat) => {
                path::<3, true>(prepared, id, poller, work)
            }
            (Algorithm::RssaCr, WaitingTime::Inverse) => {
                path::<4, false>(prepared, id, poller, work)
            }
            (Algorithm::RssaCr, WaitingTime::Ziggurat) => {
                path::<4, true>(prepared, id, poller, work)
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
        schema: ResultVersion::V01,
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
pub fn solve(
    model: &[u8],
    request: &[u8],
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    execute(&prepare(model, request)?, control).map_err(|e| vec![e])
}
pub fn solve_files(
    model_path: &Path,
    request_path: &Path,
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    let model = crate::read_document(model_path).map_err(|e| vec![e])?;
    let bytes = crate::read_document(request_path).map_err(|e| vec![e])?;
    let request: Request = crate::strict_json::parse(&bytes).map_err(|e| vec![e])?;
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
