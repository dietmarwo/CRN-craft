//! Direct SSA. Fixed replicate streams and ordered reduction, bounded queues.
use crate::{
    CompiledNetwork, check_bytes,
    contracts::{Method, Run, Statistic},
    diagnostic::Diagnostic,
    parse_model, parse_run, sampling, sha256,
};
use rand_core::RngCore;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
    time::Instant,
};

pub const MAX_REPLICATES: u64 = 1_000_000;
pub const MAX_WORKSPACE_BYTES: u64 = 256 * 1024 * 1024;
pub const MAX_RESULT_BYTES: u64 = 256 * 1024 * 1024;
pub const MAX_WORKERS: usize = 64;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.result/0.2")]
    V02,
}
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Stop {
    Completed,
    EventLimit,
    PopulationLimit,
    WallTimeLimit,
    Cancelled,
    Failed,
    OutputLimit,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Assessment {
    ExactEnsemble,
    Incomplete,
    FailedChecks,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Replicate {
    pub id: u64,
    pub started: bool,
    pub termination: Stop,
    pub events: u64,
    pub stop_time: f64,
    pub saved_points: usize,
    pub diagnostic: Option<Diagnostic>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Trajectory {
    pub id: u64,
    pub counts: Vec<Vec<i64>>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Statistics {
    pub samples: u64,
    pub mean: Option<Vec<Vec<f64>>>,
    pub sample_variance: Option<Vec<Vec<f64>>>,
    pub mean_standard_error: Option<Vec<Vec<f64>>>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SimulationResult {
    pub schema: ResultVersion,
    pub model_sha256: String,
    pub run_sha256: String,
    pub sampler: String,
    pub root_seed: u64,
    pub species: Vec<String>,
    pub save_times: Vec<f64>,
    pub termination: Stop,
    pub assessment: Assessment,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub events: u64,
    pub replicates: Vec<Replicate>,
    pub trajectories: Vec<Trajectory>,
    pub statistics: Option<Statistics>,
    pub diagnostics: Vec<Diagnostic>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    pub crate_version: String,
    pub rustc: String,
    pub target_os: String,
    pub target_arch: String,
    pub cargo_lock_sha256: String,
    pub executable_sha256: Option<String>,
    pub workers: usize,
    pub elapsed_seconds: f64,
    pub reserved_workspace_bytes: u64,
    pub workspace_cap_bytes: u64,
}
pub struct Simulation {
    pub result: SimulationResult,
    pub provenance: Provenance,
}

/// Cooperative API cancellation; CLI additionally accepts a caller-owned file.
#[derive(Default)]
pub struct Control {
    cancelled: AtomicBool,
    pub cancel_file: Option<PathBuf>,
}
impl Control {
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    fn poll(&self, start: Instant, wall: Option<f64>) -> Option<Stop> {
        if self.cancelled.load(Ordering::Relaxed)
            || self.cancel_file.as_ref().is_some_and(|p| p.exists())
        {
            return Some(Stop::Cancelled);
        }
        if wall.is_some_and(|seconds| start.elapsed().as_secs_f64() >= seconds) {
            return Some(Stop::WallTimeLimit);
        }
        None
    }
}

pub(crate) fn execution_error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E060",
        "/simulation",
        message,
        "Use Direct SSA and reduce replicas, saved output or worker count to fit the published execution limits.",
    )
}

#[derive(Debug)]
pub struct Prepared {
    pub(crate) network: CompiledNetwork,
    pub(crate) run: Run,
    pub(crate) model_hash: String,
    pub(crate) run_hash: String,
    pub(crate) workers: usize,
    pub(crate) workspace: u64,
    input_bytes: u64,
}

pub fn prepare(model_bytes: &[u8], run_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let report = check_bytes(model_bytes, Some(run_bytes));
    if !report.valid {
        return Err(report.diagnostics);
    }
    let model = parse_model(model_bytes).map_err(|e| vec![e])?;
    let run = parse_run(run_bytes).map_err(|e| vec![e])?;
    if !matches!(run.method, Method::Ssa { .. }) {
        return Err(vec![execution_error("FSP execution is not implemented")]);
    }
    let count = run
        .ensemble
        .as_ref()
        .expect("validated ensemble")
        .replicates;
    let workers = (run.budget.threads as usize)
        .min(MAX_WORKERS)
        .min(count as usize);
    let cells = (model.species.len() as u128) * run.time.save_times.len() as u128;
    // Conservative allocation reservation, not a claim to cap allocator/OS RSS.
    let workspace = (cells * 8 + run.time.save_times.len() as u128 * 24)
        * (3 * workers as u128 + u128::from(run.outputs.keep_trajectories))
        + cells * 24
        + u128::from(count) * 1024
        + workers as u128 * (2 * 1024 * 1024 + model.reactions.len() as u128 * 8)
        + 32 * 1024 * 1024;
    let evidence_reserve = model_bytes.len() as u128
        + run_bytes.len() as u128
        + 16_384
        + u128::from(count) * 1024
        + model.species.len() as u128 * 80
        + run.time.save_times.len() as u128 * 32;
    if count > MAX_REPLICATES
        || workspace > u128::from(MAX_WORKSPACE_BYTES)
        || evidence_reserve > u128::from(run.budget.max_output_bytes.min(MAX_RESULT_BYTES))
    {
        return Err(vec![execution_error(format!(
            "execution reservation exceeds limits: {workspace} workspace bytes, {evidence_reserve} evidence bytes"
        ))]);
    }
    // Ensure event accounting cannot wrap, even if every replicate reaches its cap.
    if count
        .checked_mul(run.budget.max_events_per_replicate.unwrap())
        .is_none()
    {
        return Err(vec![execution_error("aggregate event cap exceeds u64")]);
    }
    let network = CompiledNetwork::compile(&model)?;
    Ok(Prepared {
        network,
        run,
        model_hash: sha256(model_bytes),
        run_hash: sha256(run_bytes),
        workers,
        workspace: workspace as u64,
        input_bytes: (model_bytes.len() + run_bytes.len()) as u64,
    })
}

struct PathResult {
    record: Replicate,
    counts: Vec<Vec<i64>>,
}

fn trajectory(prepared: &Prepared, id: u64, control: &Control, start: Instant) -> PathResult {
    let run = &prepared.run;
    let mut rng = sampling::stream(run.ensemble.as_ref().unwrap().root_seed, id);
    path_with_rng(&prepared.network, run, id, control, start, &mut rng)
}

fn path_with_rng(
    network: &CompiledNetwork,
    run: &Run,
    id: u64,
    control: &Control,
    start: Instant,
    rng: &mut impl RngCore,
) -> PathResult {
    let mut record = Replicate {
        id,
        started: false,
        termination: Stop::Completed,
        events: 0,
        stop_time: 0.0,
        saved_points: 0,
        diagnostic: None,
    };
    let mut counts = Vec::new();
    let mut state = network.initial_counts().to_vec();
    let mut time = 0.0;
    let mut save = 0;
    let mut rates = Vec::with_capacity(network.reactions().len());
    let outcome = (|| -> Result<Stop, Diagnostic> {
        if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
            return Ok(stop);
        }
        record.started = true;
        loop {
            while save < run.time.save_times.len() && run.time.save_times[save] <= time {
                if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                    return Ok(stop);
                }
                counts.push(state.clone());
                save += 1;
            }
            if time == run.time.end {
                return Ok(Stop::Completed);
            }
            if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                return Ok(stop);
            }
            if record.events >= run.budget.max_events_per_replicate.unwrap() {
                return Ok(Stop::EventLimit);
            }
            rates.clear();
            for reaction in 0..network.reactions().len() {
                if reaction % 64 == 0 {
                    if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                        return Ok(stop);
                    }
                }
                rates.push(network.propensity_in_valid_state(reaction, &state)?);
            }
            let total = sampling::total_propensity(&rates)?;
            if total == 0.0 {
                // Absorption consumes no RNG and retains the state to the horizon.
                while save < run.time.save_times.len() {
                    if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                        return Ok(stop);
                    }
                    counts.push(state.clone());
                    save += 1;
                }
                time = run.time.end;
                return Ok(Stop::Completed);
            }
            let wait = sampling::waiting_time(total, sampling::uniform_open(rng))?;
            let next = sampling::event_time(time, wait)?;
            // Strict inequality: an event exactly at a save time comes first.
            while save < run.time.save_times.len() && run.time.save_times[save] < next {
                if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                    return Ok(stop);
                }
                counts.push(state.clone());
                save += 1;
            }
            if next > run.time.end {
                time = run.time.end;
                return Ok(Stop::Completed);
            }
            if let Some(stop) = control.poll(start, run.budget.wall_time_s) {
                return Ok(stop);
            }
            let channel = sampling::select_channel(
                &rates,
                total,
                sampling::uniform_closed_open(rng.next_u64()),
            )?;
            if let Err(error) = network.fire_reaction(channel, &mut state) {
                if error.code == "CRN-E044" && error.message.contains("safety limit") {
                    record.diagnostic = Some(error);
                    return Ok(Stop::PopulationLimit);
                }
                return Err(error);
            }
            time = next;
            record.events += 1;
        }
    })();
    match outcome {
        Ok(stop) => record.termination = stop,
        Err(error) => {
            record.termination = Stop::Failed;
            record.diagnostic = Some(error);
        }
    }
    record.stop_time = time.max(
        counts
            .len()
            .checked_sub(1)
            .map_or(0.0, |i| run.time.save_times[i]),
    );
    record.saved_points = counts.len();
    PathResult { record, counts }
}

pub fn execute(prepared: &Prepared, control: &Control) -> Result<Simulation, Diagnostic> {
    let start = Instant::now();
    let run = &prepared.run;
    let ensemble = run.ensemble.as_ref().unwrap();
    let n_species = prepared.network.species().len();
    let cells = n_species * run.time.save_times.len();
    let mut mean = vec![0.0; cells];
    let mut m2 = vec![0.0; cells];
    let mut completed = 0;
    let mut records = Vec::with_capacity(ensemble.replicates as usize);
    let mut kept = Vec::new();
    std::thread::scope(|scope| -> Result<(), Diagnostic> {
        let mut receivers = Vec::new();
        for worker in 0..prepared.workers {
            let (tx, rx) = mpsc::sync_channel(1);
            receivers.push(rx);
            std::thread::Builder::new()
                .name(format!("crn-{worker}"))
                .spawn_scoped(scope, move || {
                    for id in (worker as u64..ensemble.replicates).step_by(prepared.workers) {
                        if tx.send(trajectory(prepared, id, control, start)).is_err() {
                            break;
                        }
                    }
                })
                .map_err(|e| execution_error(format!("cannot start worker: {e}")))?;
        }
        // Per-worker bounded queues prevent an arbitrarily large reorder buffer.
        for id in 0..ensemble.replicates {
            let path = receivers[id as usize % prepared.workers]
                .recv()
                .map_err(|e| execution_error(format!("worker did not return a trajectory: {e}")))?;
            if path.record.termination == Stop::Completed {
                completed += 1;
                for (i, count) in path.counts.iter().flatten().enumerate() {
                    let delta = *count as f64 - mean[i];
                    mean[i] += delta / completed as f64;
                    m2[i] += delta * (*count as f64 - mean[i]);
                }
            }
            if id < run.outputs.keep_trajectories {
                kept.push(Trajectory {
                    id,
                    counts: path.counts,
                });
            }
            records.push(path.record);
        }
        Ok(())
    })?;
    let termination = records
        .iter()
        .map(|r| r.termination)
        .max()
        .unwrap_or(Stop::Failed);
    let assessment = if termination == Stop::Completed {
        Assessment::ExactEnsemble
    } else if termination == Stop::Failed {
        Assessment::FailedChecks
    } else {
        Assessment::Incomplete
    };
    let statistics = if completed == ensemble.replicates && !run.outputs.statistics.is_empty() {
        let chunks = |values: Vec<f64>| values.chunks(n_species).map(<[f64]>::to_vec).collect();
        let want_mean = run.outputs.statistics.contains(&Statistic::Mean);
        let variance: Vec<_> = m2
            .iter()
            .map(|x| x / (completed - 1).max(1) as f64)
            .collect();
        Some(Statistics {
            samples: completed,
            mean_standard_error: (want_mean && completed > 1).then(|| {
                chunks(
                    variance
                        .iter()
                        .map(|v| libm::sqrt(v / completed as f64))
                        .collect(),
                )
            }),
            mean: want_mean.then(|| chunks(mean)),
            sample_variance: run
                .outputs
                .statistics
                .contains(&Statistic::Variance)
                .then(|| chunks(variance)),
        })
    } else {
        None
    };
    let mut simulation = Simulation {
        result: SimulationResult {
            schema: ResultVersion::V02,
            model_sha256: prepared.model_hash.clone(),
            run_sha256: prepared.run_hash.clone(),
            sampler: sampling::SAMPLER.into(),
            root_seed: ensemble.root_seed,
            species: prepared.network.species().to_vec(),
            save_times: run.time.save_times.clone(),
            termination,
            assessment,
            requested_replicates: ensemble.replicates,
            completed_replicates: completed,
            events: records.iter().map(|r| r.events).sum(),
            replicates: records,
            trajectories: kept,
            statistics,
            diagnostics: prepared.network.warnings().to_vec(),
        },
        provenance: Provenance {
            crate_version: env!("CARGO_PKG_VERSION").into(),
            rustc: env!("CRN_RUSTC_VERSION").into(),
            target_os: std::env::consts::OS.into(),
            target_arch: std::env::consts::ARCH.into(),
            cargo_lock_sha256: sha256(include_bytes!("../../../Cargo.lock")),
            executable_sha256: None,
            workers: prepared.workers,
            elapsed_seconds: start.elapsed().as_secs_f64(),
            reserved_workspace_bytes: prepared.workspace,
            workspace_cap_bytes: MAX_WORKSPACE_BYTES,
        },
    };
    constrain_output(
        &mut simulation.result,
        run.budget
            .max_output_bytes
            .min(MAX_RESULT_BYTES)
            .saturating_sub(prepared.input_bytes + 16_384),
    )?;
    simulation.provenance.elapsed_seconds = start.elapsed().as_secs_f64();
    Ok(simulation)
}

pub(crate) fn constrain_output(
    result: &mut SimulationResult,
    limit: u64,
) -> Result<(), Diagnostic> {
    struct Counter {
        used: u64,
        limit: u64,
    }
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() as u64 > self.limit.saturating_sub(self.used) {
                return Err(std::io::Error::other("output limit"));
            }
            self.used += bytes.len() as u64;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let fits =
        |value: &SimulationResult| serde_json::to_writer(Counter { used: 1, limit }, value).is_ok();
    if fits(result) {
        return Ok(());
    }
    result.termination = Stop::OutputLimit;
    result.assessment = Assessment::Incomplete;
    result.statistics = None;
    result.trajectories.clear();
    let diagnostic = Diagnostic::error(
        "CRN-E064",
        "/outputs",
        "serialized numerical payload exceeded the output byte budget",
        "Increase max_output_bytes; trajectories/statistics were removed, replicate accounting was retained.",
    );
    result.diagnostics.push(diagnostic.clone());
    if fits(result) {
        Ok(())
    } else {
        Err(diagnostic)
    }
}

pub fn simulate_bytes(
    model: &[u8],
    run: &[u8],
    control: &Control,
) -> Result<Simulation, Vec<Diagnostic>> {
    execute(&prepare(model, run)?, control).map_err(|e| vec![e])
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::Error;
    use serde_json::json;

    struct Words {
        values: Vec<u64>,
        index: usize,
    }
    impl RngCore for Words {
        fn next_u32(&mut self) -> u32 {
            panic!("unexpected 32-bit draw")
        }
        fn next_u64(&mut self) -> u64 {
            let value = self.values[self.index];
            self.index += 1;
            value
        }
        fn fill_bytes(&mut self, _: &mut [u8]) {
            panic!("unexpected byte draw")
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), Error> {
            panic!("unexpected byte draw")
        }
    }
    #[test]
    fn open_uniform_rejects_zero_and_never_returns_one() {
        let mut rng = Words {
            values: vec![0, u64::MAX],
            index: 0,
        };
        assert_eq!(
            sampling::uniform_open(&mut rng),
            1.0 - sampling::UNIFORM_STEP
        );
        assert_eq!(rng.index, 2);
    }
    #[test]
    fn event_exactly_at_horizon_is_included_and_consumes_two_words() {
        let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"boundary","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":0}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1.0}}]})).unwrap();
        let end = libm::log(2.0);
        let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
            "ensemble":{"replicates":1,"root_seed":0},"time":{"end":end,"save_times":[0.0,end]},"outputs":{"statistics":["mean"],"keep_trajectories":1},
            "budget":{"max_events_per_replicate":1,"threads":1,"max_output_bytes":100000}})).unwrap();
        let prepared = prepare(&model, &run).unwrap();
        let mut rng = Words {
            values: vec![1_u64 << 63, 0],
            index: 0,
        };
        let path = path_with_rng(
            &prepared.network,
            &prepared.run,
            0,
            &Control::default(),
            Instant::now(),
            &mut rng,
        );
        assert_eq!(path.record.termination, Stop::Completed);
        assert_eq!(path.record.events, 1);
        assert_eq!(path.counts, vec![vec![0], vec![1]]);
        assert_eq!(rng.index, 2);
    }
    #[test]
    fn encoded_budget_checks_exact_boundary_without_allocating_payload() {
        let mut r = SimulationResult {
            schema: ResultVersion::V02,
            model_sha256: "a".repeat(64),
            run_sha256: "b".repeat(64),
            sampler: sampling::SAMPLER.into(),
            root_seed: 1,
            species: vec!["X".into()],
            save_times: vec![0.0, 1.0],
            termination: Stop::Completed,
            assessment: Assessment::ExactEnsemble,
            requested_replicates: 0,
            completed_replicates: 0,
            events: 0,
            replicates: vec![],
            trajectories: vec![],
            statistics: None,
            diagnostics: vec![],
        };
        let exact = serde_json::to_vec(&r).unwrap().len() as u64 + 1;
        constrain_output(&mut r, exact).unwrap();
        assert_eq!(r.termination, Stop::Completed);
        assert!(constrain_output(&mut r, exact - 1).is_err());
        assert_eq!(r.termination, Stop::OutputLimit);
    }
}
