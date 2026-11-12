//! Finite state projection (Munsky & Khammash 2006), row-oriented generator.
//! Positive sparse uniformization, split so lambda*dt <= 16. No normalization.
//!
//! The numerical core enumerates the states reachable from the initial
//! counts inside an inclusive `bounds` box (at most `max_states`, and never
//! more than [`MAX_STATES`]), routes every transition that leaves the box to
//! one absorbing sink, and propagates the probability distribution to each
//! report time by uniformization. Each snapshot reports the retained
//! probabilities, the probability in the sink and the error terms
//! separately; nothing is renormalized, and `numerical_error_certified` is
//! always false. Capacity, work, output and interruption failures are
//! `CRN-E070` and return no partial result.
//!
//! Two entry points share that core:
//!
//! - [`v1::solve`] reads a `crn-craft.run/1.0` request whose method is
//!   `fsp` ([`v1::FspRequest`]) for any closed rules model
//!   (`crn-craft.model/1.0` or a `crn-craft.rules-model`), through its
//!   [`FlatView`](crate::rules::FlatView), and returns
//!   `crn-craft.fsp-result/1.0` ([`v1::FspResult`]) with the flat species
//!   labels. A state where the model's `stop` condition holds is absorbing,
//!   as a stopped replicate holds its state. This is what `crn-craft rules
//!   MODEL REQUEST` runs for such a request.
//! - [`solve`] reads a network `run/0.1` request with `method.kind = "fsp"`
//!   for a `crn-craft.model/0.1` network model and returns
//!   `crn-craft.fsp-result/0.1` ([`FspResult`]); this is the `crn-craft fsp`
//!   command, kept until the network engine is removed.
//!
//! ```
//! use crn_craft::{fsp, sha256, simulation::Control};
//!
//! let model = br#"{"schema": "crn-craft.model/0.1", "name": "decay",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
//!   "species": [{"name": "A", "initial_count": 2}],
//!   "reactions": [{"name": "death", "reactants": {"A": 1}, "products": {},
//!                  "rate": {"kind": "mass_action", "coefficient": 1.0}}]}"#;
//! let run = format!(r#"{{"schema": "crn-craft.run/0.1",
//!   "model": {{"path": "model.json", "sha256": "{}"}},
//!   "method": {{"kind": "fsp", "bounds": {{"A": 2}}, "max_states": 10, "tolerance": 1e-12}},
//!   "time": {{"end": 1.0, "save_times": [0.0, 1.0]}},
//!   "outputs": {{"statistics": ["mean"], "keep_trajectories": 0}},
//!   "budget": {{"threads": 1, "max_output_bytes": 100000}}}}"#, sha256(model));
//! let result = fsp::solve(model, run.as_bytes(), &Control::default()).unwrap();
//! assert_eq!(result.states, vec![vec![2], vec![1], vec![0]]);
//! let p_two = result.snapshots[1].probabilities[0];
//! assert!((p_two - (-2.0f64).exp()).abs() < 1e-12); // both molecules survive: e^{-2}
//! assert_eq!(result.snapshots[1].physical_leakage_estimate, 0.0);
//! ```
use crate::{
    CompiledNetwork, check_bytes,
    contracts::Method,
    diagnostic::Diagnostic,
    parse_model, parse_run, sha256,
    simulation::{Control, Poller},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, time::Instant};

/// Most projected states, whatever `max_states` says (`CRN-E070` beyond).
pub const MAX_STATES: usize = 100_000;
/// Most sparse generator edges (`CRN-E070` beyond).
pub const MAX_EDGES: usize = 2_000_000;
/// Most deterministic work units for enumeration and propagation (`CRN-E070` beyond).
pub const MAX_WORK: u64 = 200_000_000;
const MAX_BYTES: u128 = 256 * 1024 * 1024;

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of an FSP result.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum FspVersion {
    #[doc = concat!("`crn-craft.fsp-result/0.1`.")]
    #[serde(rename = "crn-craft.fsp-result/0.1")]
    V01,
}
/// The projected distribution at one save time and its separately reported error terms.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FspSnapshot {
    /// Save time in seconds.
    pub time: f64,
    /// Index order matches states, excludes the absorbing escaping-mass sink.
    pub probabilities: Vec<f64>,
    /// Probability in the sink: mass that left the bounds box by this time.
    pub physical_leakage_estimate: f64,
    /// Unnormalized moments of the retained law, never of the full unbounded process.
    pub projected_first_moment: Vec<f64>,
    /// Σ p(x)·x² per species over the retained states, unnormalized like the first moment.
    pub projected_second_raw_moment: Vec<f64>,
    /// Cumulative discarded positive Poisson weight, excluding roundoff.
    pub series_tail_upper_bound: f64,
    /// Probability mass dropped by truncating the uniformization series so far.
    pub discarded_weight_estimate: f64,
    /// 1 − (retained probability + sink probability): the mass missing from the vector.
    pub mass_defect: f64,
    /// Mass-balance residual only; not a rigorous floating-point error bound.
    /// It is |mass_defect − discarded_weight_estimate|.
    pub roundoff_mass_residual: f64,
}
/// Finite state projection result: the enumerated states and one snapshot per save time.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FspResult {
    #[doc = concat!("Always `crn-craft.fsp-result/0.1`.")]
    pub schema: FspVersion,
    /// SHA-256 of the exact model bytes.
    pub model_sha256: String,
    /// SHA-256 of the exact run bytes.
    pub run_sha256: String,
    /// Species names in the order of each state's counts.
    pub species: Vec<String>,
    /// Projected states in discovery order; state 0 is the initial state.
    pub states: Vec<Vec<i64>>,
    /// One snapshot per save time, in order, starting at time 0.
    pub snapshots: Vec<FspSnapshot>,
    /// Uniformization rate λ: the largest total exit rate of any projected state, per second.
    pub uniformization_rate: f64,
    /// Nonzero off-diagonal generator entries, transitions to the sink included.
    pub sparse_edges: usize,
    /// Deterministic work units spent, at most 200000000.
    pub work: u64,
    /// Always false: the reported error terms are estimates, not a certified bound.
    pub numerical_error_certified: bool,
}
fn error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E070",
        "/fsp",
        message,
        "Reduce the projection/work or revise bounds/tolerance. No boundary is reflected or normalized.",
    )
}
/// Work accounting and the interruption poll of one solve. `poll` returns
/// the reason of a stop.
struct Work<'a> {
    poll: &'a mut dyn FnMut() -> Option<String>,
    used: u64,
}
impl Work<'_> {
    fn add(&mut self, amount: usize) -> Result<(), Diagnostic> {
        self.used = self.used.saturating_add(amount as u64);
        if self.used > MAX_WORK {
            return Err(error("FSP work limit exceeded"));
        }
        if let Some(stop) = (self.poll)() {
            return Err(error(format!("FSP interrupted: {stop}")));
        }
        Ok(())
    }
}

/// The Poisson weight omitted after term `first - 1`: the sum of the
/// weights from index `first` on, starting with `first_weight`, summed
/// explicitly rather than as `1 - kept weight`, which cancels to roundoff.
/// The series stops only once `mu / (first + 1) < 1`, so the terms decrease.
fn omitted_weight(first_weight: f64, mu: f64, first: u64) -> f64 {
    let mut term = first_weight;
    let mut sum = 0.0;
    let mut index = first;
    while term > 0.0 && term > sum * f64::EPSILON {
        sum += term;
        index += 1;
        term *= mu / index as f64;
    }
    sum
}

/// Compact JSON length plus the trailing newline the CLI writes.
fn serialized_len(result: &impl Serialize) -> u64 {
    struct Counter(u64);
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 += bytes.len() as u64;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut counter = Counter(1);
    serde_json::to_writer(&mut counter, result).expect("finite FSP result");
    counter.0
}

/// Normalized rows: the diagonal weight and the off-diagonal (target, weight)
/// edges of each projected state; `usize::MAX` targets the escaping sink.
type Rows = [(f64, Vec<(usize, f64)>)];

/// The same generator by column: the incoming (source, weight) entries of
/// each target, the sink last, each in ascending source order with the
/// diagonal at its own position. Gathering a column adds its contributions
/// in the order the row-oriented product scattered them.
struct Columns {
    start: Vec<usize>,
    entries: Vec<(usize, f64)>,
}

fn columns(rows: &Rows) -> Columns {
    let n = rows.len();
    let target = |j: usize| if j == usize::MAX { n } else { j };
    let mut start = vec![0; n + 2];
    for (i, (_, edges)) in rows.iter().enumerate() {
        start[i + 1] += 1;
        for &(j, _) in edges {
            start[target(j) + 1] += 1;
        }
    }
    for c in 1..start.len() {
        start[c] += start[c - 1];
    }
    let mut fill = start.clone();
    let mut entries = vec![(0, 0.0); start[n + 1]];
    for (i, (diagonal, edges)) in rows.iter().enumerate() {
        entries[fill[i]] = (i, *diagonal);
        fill[i] += 1;
        for &(j, weight) in edges {
            let c = target(j);
            entries[fill[c]] = (i, weight);
            fill[c] += 1;
        }
    }
    Columns { start, entries }
}

/// One series term, `next = term · P`, with the sink (index n) keeping its
/// mass, fused with the series accumulation `result += weight * next`. Each
/// entry of `next` is the same ascending-source sum as the row-oriented
/// product computed, so results are bitwise unchanged.
fn multiply_accumulate(
    columns: &Columns,
    term: &[f64],
    next: &mut [f64],
    result: &mut [f64],
    weight: f64,
) {
    let n = columns.start.len() - 2;
    for (c, ((out, total), bounds)) in next
        .iter_mut()
        .zip(result.iter_mut())
        .zip(columns.start.windows(2))
        .enumerate()
    {
        let mut sum = if c == n { term[n] } else { 0.0 };
        for &(i, rate) in &columns.entries[bounds[0]..bounds[1]] {
            sum += term[i] * rate;
        }
        *out = sum;
        *total += weight * sum;
    }
}

/// The network the numerical core projects: a fixed species list, channels
/// with their net changes, and the propensities of a state.
trait Generator {
    /// The initial state.
    fn initial(&self) -> Vec<i64>;
    /// Net changes of one firing of each channel, as (species, change).
    fn changes(&self) -> Vec<Vec<(usize, i64)>>;
    /// Writes every channel's propensity at `state` into `rates` and
    /// returns `false`; or, for an absorbing state (a `stop` condition
    /// holds), writes zeros without evaluating any channel and returns
    /// `true`.
    fn rates(&mut self, state: &[i64], rates: &mut Vec<f64>) -> Result<bool, Diagnostic>;
    /// Checks a newly reached state before it is projected.
    fn admit(&self, state: &[i64]) -> Result<(), Diagnostic>;
}

/// The settings of one projection.
struct Settings<'a> {
    /// Inclusive upper bound per species, in state order.
    bounds: &'a [i64],
    /// Most states to enumerate.
    capacity: usize,
    /// Report times, strictly increasing and nonnegative.
    times: &'a [f64],
    /// Total Poisson-series truncation budget.
    tolerance: f64,
    /// Output budget of the early reservation check.
    max_output_bytes: u64,
}

/// What the core computes: states, snapshots and generator facts.
struct Projection {
    states: Vec<Vec<i64>>,
    snapshots: Vec<FspSnapshot>,
    uniformization_rate: f64,
    sparse_edges: usize,
}

/// The numerical core: enumeration, uniformization and error accounting.
#[allow(clippy::too_many_lines)]
fn project<G: Generator>(
    network: &mut G,
    settings: &Settings<'_>,
    work: &mut Work<'_>,
) -> Result<Projection, Diagnostic> {
    let bounds = settings.bounds;
    let changes = network.changes();
    let mut states = vec![network.initial()];
    let mut index = BTreeMap::from([(states[0].clone(), 0)]);
    let mut rows = Vec::new();
    let mut edge_count = 0;
    let mut lambda: f64 = 0.0;
    let mut rates = Vec::with_capacity(changes.len());
    let mut i = 0;
    while i < states.len() {
        work.add(1)?;
        network.rates(&states[i], &mut rates)?;
        let mut row = BTreeMap::new();
        for (r, &rate) in rates.iter().enumerate() {
            work.add(1)?;
            if rate == 0.0 {
                continue;
            }
            let mut next = states[i].clone();
            let mut escaped = false;
            for &(species, coefficient) in &changes[r] {
                let value = i128::from(next[species]) + i128::from(coefficient);
                if value < 0 {
                    return Err(error("negative reachable count"));
                }
                if value > i128::from(bounds[species]) {
                    escaped = true;
                    break;
                }
                next[species] = value as i64;
            }
            let target = if escaped {
                usize::MAX
            } else if let Some(&j) = index.get(&next) {
                j
            } else {
                if states.len() >= settings.capacity {
                    return Err(error(
                        "reachable-state cap exceeded; no incomplete generator is solved",
                    ));
                }
                network.admit(&next)?;
                let j = states.len();
                index.insert(next.clone(), j);
                states.push(next);
                j
            };
            // A physical null reaction contributes no generator exit or sink.
            if target != i {
                *row.entry(target).or_insert(0.0) += rate;
            }
        }
        let exit: f64 = row.values().sum();
        if !exit.is_finite() {
            return Err(error("generator exit rate overflows"));
        }
        lambda = lambda.max(exit);
        edge_count += row.len();
        if edge_count > MAX_EDGES {
            return Err(error("FSP sparse-edge cap exceeded"));
        }
        // Retained snapshots: probabilities, two moment vectors and the record.
        let memory = states.len() as u128 * (bounds.len() as u128 * 32 + 128)
            + edge_count as u128 * 64
            + states.len() as u128 * (settings.times.len() as u128 * 8 + 64)
            + settings.times.len() as u128 * (bounds.len() as u128 * 16 + 256);
        if memory > MAX_BYTES {
            return Err(error("FSP workspace reservation exceeds 256 MiB"));
        }
        rows.push((exit, row.into_iter().collect::<Vec<_>>()));
        i += 1;
    }
    let n = states.len();
    // An early estimate; the serialized result is checked exactly at the end.
    let output_reserve = n as u128 * (bounds.len() as u128 * 24 + 32)
        + n as u128 * settings.times.len() as u128 * 32
        + settings.times.len() as u128 * (bounds.len() as u128 * 48 + 400)
        + 16384;
    if output_reserve > u128::from(settings.max_output_bytes) {
        return Err(error("FSP output reservation exceeds budget"));
    }
    // Uniformization steps of each interval from 0 through the report times.
    let total_steps_f = settings
        .times
        .iter()
        .scan(0.0, |previous, &time| {
            let dt = time - *previous;
            *previous = time;
            Some(dt)
        })
        .filter(|&dt| dt > 0.0)
        .map(|dt| (dt * lambda / 16.0).ceil().max(1.0))
        .sum::<f64>();
    if !total_steps_f.is_finite() || total_steps_f > MAX_WORK as f64 {
        return Err(error("FSP uniformization step limit exceeded"));
    }
    let tail_per_step = settings.tolerance / total_steps_f;
    if tail_per_step < 1e-18 {
        return Err(error(
            "requested series tolerance is below the supported per-step range",
        ));
    }
    // Normalize each immutable row once; retain the exact per-entry division.
    if lambda > 0.0 {
        for (diagonal, edges) in &mut rows {
            *diagonal = 1.0 - *diagonal / lambda;
            for (_, weight) in edges {
                *weight /= lambda;
                if *weight == 0.0 {
                    return Err(error(
                        "positive FSP transition underflows during normalization",
                    ));
                }
            }
        }
    }
    let columns = columns(&rows);
    drop(rows);
    let mut term = vec![0.0; n + 1];
    let mut next = vec![0.0; n + 1];
    let mut result = vec![0.0; n + 1];
    let mut p = vec![0.0; n + 1];
    p[0] = 1.0;
    let mut snapshots = Vec::new();
    let mut tail_bound = 0.0;
    let mut discarded = 0.0;
    let mut previous = 0.0;
    for &time in settings.times {
        let dt = time - previous;
        if lambda > 0.0 && dt > 0.0 {
            let steps = (dt * lambda / 16.0).ceil().max(1.0) as u64;
            let mu = lambda * (dt / steps as f64);
            for _ in 0..steps {
                term.clone_from(&p);
                let mut weight = libm::exp(-mu);
                for (out, value) in result.iter_mut().zip(&term) {
                    *out = value * weight;
                }
                let mut k = 0;
                loop {
                    work.add(n + edge_count)?;
                    let next_weight = weight * mu / (k + 1) as f64;
                    let ratio = mu / (k + 2) as f64;
                    if ratio < 1.0 && next_weight / (1.0 - ratio) <= tail_per_step {
                        tail_bound += next_weight / (1.0 - ratio);
                        discarded += omitted_weight(next_weight, mu, k + 1) * p.iter().sum::<f64>();
                        break;
                    }
                    multiply_accumulate(&columns, &term, &mut next, &mut result, next_weight);
                    std::mem::swap(&mut term, &mut next);
                    weight = next_weight;
                    k += 1;
                }
                std::mem::swap(&mut p, &mut result);
            }
        }
        // Moments and checks visit every state once per species.
        work.add(n * (2 * bounds.len() + 1))?;
        let mass: f64 = p.iter().sum();
        if p.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(error("invalid probability after uniformization"));
        }
        snapshots.push(FspSnapshot {
            time,
            probabilities: p[..n].to_vec(),
            physical_leakage_estimate: p[n],
            projected_first_moment: (0..bounds.len())
                .map(|j| states.iter().zip(&p).map(|(x, q)| x[j] as f64 * q).sum())
                .collect(),
            projected_second_raw_moment: (0..bounds.len())
                .map(|j| {
                    states
                        .iter()
                        .zip(&p)
                        .map(|(x, q)| (x[j] as f64).powi(2) * q)
                        .sum()
                })
                .collect(),
            series_tail_upper_bound: tail_bound,
            discarded_weight_estimate: discarded,
            mass_defect: 1.0 - mass,
            roundoff_mass_residual: (1.0 - mass - discarded).abs(),
        });
        previous = time;
    }
    Ok(Projection {
        states,
        snapshots,
        uniformization_rate: lambda,
        sparse_edges: edge_count,
    })
}

/// A compiled network model as the core's generator.
struct Network<'a>(&'a CompiledNetwork);

impl Generator for Network<'_> {
    fn initial(&self) -> Vec<i64> {
        self.0.initial_counts().to_vec()
    }
    fn changes(&self) -> Vec<Vec<(usize, i64)>> {
        self.0
            .reactions()
            .iter()
            .map(|r| {
                r.net_changes
                    .iter()
                    .map(|t| (t.species, t.coefficient))
                    .collect()
            })
            .collect()
    }
    fn rates(&mut self, state: &[i64], rates: &mut Vec<f64>) -> Result<bool, Diagnostic> {
        *rates = self.0.propensities(state)?;
        Ok(false)
    }
    fn admit(&self, state: &[i64]) -> Result<(), Diagnostic> {
        self.0.validate_state(state)
    }
}

/// Solves an FSP run against the exact model bytes; see the [module docs](self).
///
/// `model.path` is only syntax-checked. The run's wall-time limit and the
/// [`Control`] are polled during the computation; a stop fails the solve.
///
/// # Errors
///
/// The [`crate::check_bytes`] errors of model and run (including
/// `CRN-E023` for invalid bounds, capacity or tolerance), and `CRN-E070` for
/// an SSA run, exceeded state, edge, work, memory or output limits, a
/// tolerance below 10⁻¹⁸ per step, a transition that underflows, invalid
/// probabilities or an interruption. A reachable state above a population
/// limit (`CRN-E041`) or a rate that cannot be evaluated (`CRN-E042`) also
/// fails the solve.
pub fn solve(
    model_bytes: &[u8],
    run_bytes: &[u8],
    control: &Control,
) -> Result<FspResult, Vec<Diagnostic>> {
    let checked = check_bytes(model_bytes, Some(run_bytes));
    if !checked.valid {
        return Err(checked.diagnostics);
    }
    let model = parse_model(model_bytes).map_err(|e| vec![e])?;
    let run = parse_run(run_bytes).map_err(|e| vec![e])?;
    let Method::Fsp {
        bounds,
        max_states,
        tolerance,
    } = &run.method
    else {
        return Err(vec![error("FSP requires an FSP run")]);
    };
    let network = CompiledNetwork::compile(&model)?;
    let mut poller = Poller::new(control, Instant::now(), run.budget.wall_time_s);
    let mut poll = || poller.poll().map(|stop| format!("{stop:?}"));
    let mut work = Work {
        poll: &mut poll,
        used: 0,
    };
    let action = (|| -> Result<FspResult, Diagnostic> {
        let bounds: Vec<_> = network.species().iter().map(|s| bounds[s]).collect();
        let settings = Settings {
            bounds: &bounds,
            capacity: (*max_states).min(MAX_STATES as u64) as usize,
            times: &run.time.save_times,
            tolerance: *tolerance,
            max_output_bytes: run.budget.max_output_bytes,
        };
        let projection = project(&mut Network(&network), &settings, &mut work)?;
        let result = FspResult {
            schema: FspVersion::V01,
            model_sha256: sha256(model_bytes),
            run_sha256: sha256(run_bytes),
            species: network.species().to_vec(),
            states: projection.states,
            snapshots: projection.snapshots,
            uniformization_rate: projection.uniformization_rate,
            sparse_edges: projection.sparse_edges,
            work: work.used,
            numerical_error_certified: false,
        };
        if serialized_len(&result) > run.budget.max_output_bytes {
            return Err(error("FSP result exceeds max_output_bytes"));
        }
        Ok(result)
    })();
    action.map_err(|e| vec![e])
}

pub mod v1 {
    //! Finite state projection of a closed rules model: the `fsp` method of
    //! `crn-craft.run/1.0` and its result `crn-craft.fsp-result/1.0`.
    //!
    //! The request ([`FspRequest`]) binds a model by hash like a simulation
    //! request and has `time`, a `method` of kind `fsp` and a `budget`; it
    //! has no ensemble, outputs, trackers or event log. `bounds` gives an
    //! inclusive upper count for every species label of the model's
    //! [`FlatView`] (`A`, or `mtDNA{wt=1}` for an attributed type). The
    //! result ([`FspResult`]) has the fields of `crn-craft.fsp-result/0.1`,
    //! with the labels as `species` and one snapshot per report time.
    //!
    //! ```
    //! use crn_craft::fsp::v1;
    //! use crn_craft::rules::Control;
    //! use serde_json::json;
    //!
    //! // Immigration-death from 0 copies: Poisson with mean 10 (1 - e^-t).
    //! let model = serde_json::to_vec(&json!({
    //!     "schema": "crn-craft.model/1.0", "name": "immigration_death",
    //!     "parameters": {"k_in": 10.0, "k_out": 1.0},
    //!     "types": [{"name": "X"}], "initial": [],
    //!     "observables": {"X": {"type": "X"}},
    //!     "rules": [
    //!         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
    //!         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
    //!     ]
    //! }))
    //! .unwrap();
    //! let request = serde_json::to_vec(&json!({
    //!     "schema": "crn-craft.run/1.0",
    //!     "model": {"path": "model.json", "sha256": crn_craft::sha256(&model)},
    //!     "time": {"end": 2.0, "report_times": [2.0]},
    //!     "method": {"kind": "fsp", "bounds": {"X": 80}, "max_states": 100, "tolerance": 1e-12},
    //!     "budget": {"max_output_bytes": 1000000}
    //! }))
    //! .unwrap();
    //! let result = v1::solve(&model, &request, &Control::default()).unwrap();
    //! assert_eq!(result.species, ["X"]);
    //! let mean = 10.0 * (1.0 - (-2.0f64).exp());
    //! let p0 = result.snapshots[0].probabilities[0]; // state 0 is the initial state, X = 0
    //! assert!((p0 - (-mean).exp()).abs() < 1e-12);
    //! assert!((result.snapshots[0].projected_first_moment[0] - mean).abs() < 1e-9);
    //! ```
    use super::{Generator, Settings, Work, error, project, serialized_len};
    use crate::contracts::ModelRef;
    use crate::diagnostic::{Diagnostic, pointer};
    use crate::rules::flat_view::Evaluator;
    use crate::rules::{Control, FlatView, MAX_REPORT_TIMES, MAX_RESULT_BYTES, RunVersion, Time};
    use crate::sha256;
    use schemars::{JsonSchema, schema_for};
    use serde::{Deserialize, Serialize};
    use serde_json::Value;
    use std::collections::BTreeMap;
    use std::path::Path;
    use std::time::Instant;

    /// A `crn-craft.run/1.0` request whose method is `fsp`: compute the
    /// probability distribution of a closed model's counts at the report
    /// times by finite state projection.
    ///
    /// Unknown fields and duplicate keys are rejected; a request without
    /// `ensemble`, `outputs`, `trackers` and `event_log`.
    #[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
    #[serde(deny_unknown_fields)]
    pub struct FspRequest {
        // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
        #[doc = concat!("Always `crn-craft.run/1.0` (`CRN-R020` at `/schema` otherwise).")]
        pub schema: RunVersion,
        #[doc = concat!(
            "The model by portable relative path below the request's directory ",
            "and by the SHA-256 of its exact bytes."
        )]
        pub model: ModelRef,
        #[doc = concat!(
            "Horizon and report grid: one snapshot per report time. The ",
            "distribution starts at time 0 from the initial counts."
        )]
        pub time: Time,
        #[doc = concat!("The `fsp` method and its truncation.")]
        pub method: FspMethod,
        #[doc = concat!("Time and output limits.")]
        pub budget: FspBudget,
    }

    /// The method of an FSP request, tagged by `kind`.
    #[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
    #[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
    pub enum FspMethod {
        /// Finite state projection on the box of the bounds: the states
        /// reachable from the initial counts inside it, with every transition
        /// that leaves it moved to an absorbing sink. Result
        /// `crn-craft.fsp-result/1.0`.
        Fsp {
            /// Inclusive upper count bound for every species label of the
            /// model's flat view (`A`, or `T{a=1,b=0.5}` for an attributed
            /// type), each at least the initial count. Probability that leaves
            /// the box moves to an absorbing sink.
            bounds: BTreeMap<String, u64>,
            /// Most reachable states to enumerate, ≥ 1; execution also caps it
            /// at 100000.
            #[schemars(range(min = 1))]
            max_states: u64,
            /// Total Poisson-series truncation budget over the whole run,
            /// 0 < tolerance < 1. It does not bound floating-point error.
            tolerance: f64,
        },
    }

    /// Resource limits of an FSP request.
    #[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
    #[serde(deny_unknown_fields)]
    pub struct FspBudget {
        /// Size limit of the result as compact JSON plus a trailing newline,
        /// positive; execution never uses more than 1 GiB.
        #[schemars(range(min = 1))]
        pub max_output_bytes: u64,
        /// Optional wall-clock limit in seconds, positive and finite. It is
        /// polled cooperatively; reaching it fails the solve (`CRN-E070`).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub wall_time_s: Option<f64>,
    }

    /// Schema tag of a `crn-craft.fsp-result/1.0` result.
    #[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
    pub enum FspVersion {
        #[doc = concat!("`crn-craft.fsp-result/1.0`.")]
        #[serde(rename = "crn-craft.fsp-result/1.0")]
        V10,
    }

    /// The projected distribution at one report time and its separately
    /// reported error terms.
    #[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
    #[serde(deny_unknown_fields)]
    pub struct FspSnapshot {
        /// Report time.
        pub time: f64,
        /// Index order matches states, excludes the absorbing escaping-mass sink.
        pub probabilities: Vec<f64>,
        /// Probability in the sink: mass that left the bounds box by this time.
        pub physical_leakage_estimate: f64,
        /// Σ p(x)·x per species over the retained states: unnormalized moments
        /// of the retained law, never of the full unbounded process.
        pub projected_first_moment: Vec<f64>,
        /// Σ p(x)·x² per species over the retained states, unnormalized like the first moment.
        pub projected_second_raw_moment: Vec<f64>,
        /// Cumulative discarded positive Poisson weight, excluding roundoff.
        pub series_tail_upper_bound: f64,
        /// Probability mass dropped by truncating the uniformization series so far.
        pub discarded_weight_estimate: f64,
        /// 1 − (retained probability + sink probability): the mass missing from the vector.
        pub mass_defect: f64,
        /// Mass-balance residual only; not a rigorous floating-point error bound.
        /// It is |mass_defect − discarded_weight_estimate|.
        pub roundoff_mass_residual: f64,
    }

    impl From<super::FspSnapshot> for FspSnapshot {
        fn from(s: super::FspSnapshot) -> Self {
            Self {
                time: s.time,
                probabilities: s.probabilities,
                physical_leakage_estimate: s.physical_leakage_estimate,
                projected_first_moment: s.projected_first_moment,
                projected_second_raw_moment: s.projected_second_raw_moment,
                series_tail_upper_bound: s.series_tail_upper_bound,
                discarded_weight_estimate: s.discarded_weight_estimate,
                mass_defect: s.mass_defect,
                roundoff_mass_residual: s.roundoff_mass_residual,
            }
        }
    }

    /// The result of an FSP request, `crn-craft.fsp-result/1.0`: the
    /// enumerated states and one snapshot per report time.
    #[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
    #[serde(deny_unknown_fields)]
    pub struct FspResult {
        #[doc = concat!("Always `crn-craft.fsp-result/1.0`.")]
        pub schema: FspVersion,
        /// SHA-256 of the exact model bytes.
        pub model_sha256: String,
        /// SHA-256 of the exact request bytes.
        pub run_sha256: String,
        /// Species labels of the model's flat view, in the order of each
        /// state's counts: by type in declaration order, `A` for a type
        /// without attributes and `T{a=1,b=0.5}` otherwise.
        pub species: Vec<String>,
        /// Projected states in discovery order; state 0 is the initial state.
        pub states: Vec<Vec<u64>>,
        /// One snapshot per report time, in order.
        pub snapshots: Vec<FspSnapshot>,
        /// Uniformization rate λ: the largest total exit rate of any projected
        /// state, per unit of time.
        pub uniformization_rate: f64,
        /// Nonzero off-diagonal generator entries, transitions to the sink included.
        pub sparse_edges: usize,
        /// Deterministic work units spent, at most 200000000.
        pub work: u64,
        /// Always false: the reported error terms are estimates, not a certified bound.
        pub numerical_error_certified: bool,
    }

    /// Whether request bytes ask for the `fsp` method: a JSON object whose
    /// `method.kind` is `"fsp"`. It reads the document leniently, so
    /// [`solve`] still reports any strict-parsing error.
    #[must_use]
    pub fn is_fsp_request(request_bytes: &[u8]) -> bool {
        serde_json::from_slice::<Value>(request_bytes)
            .ok()
            .and_then(|v| v.pointer("/method/kind").map(|k| k == "fsp"))
            .unwrap_or(false)
    }

    /// [`is_fsp_request`] for a request file; `false` when it cannot be read.
    #[must_use]
    pub fn is_fsp_request_file(request_path: &Path) -> bool {
        crate::rules::read(request_path).is_ok_and(|bytes| is_fsp_request(&bytes))
    }

    /// The refusal of an FSP request where a simulation request is expected
    /// (`prepare`): `CRN-R020` at `/method`.
    pub(crate) fn not_a_simulation() -> Diagnostic {
        request_error(
            "/method",
            "the fsp method computes a distribution, not an ensemble",
            "Solve FSP requests with crn_craft::fsp::v1::solve, `crn-craft rules MODEL REQUEST` or crn_craft.rules.fsp.",
        )
    }

    fn request_error(path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
        Diagnostic::error("CRN-R020", path, message, hint)
    }

    /// The checks that need no model: version, binding, time, method and
    /// budget. The messages of shared fields are those of simulation
    /// requests.
    fn validate(request: &FspRequest, model_sha256: &str) -> Vec<Diagnostic> {
        let mut out = Vec::new();
        if request.schema != RunVersion::V10 {
            out.push(request_error(
                "/schema",
                "the fsp method requires crn-craft.run/1.0",
                "Select the crn-craft.run/1.0 request contract.",
            ));
        }
        if !crate::validate::valid_relative_path(&request.model.path) {
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
        let FspMethod::Fsp {
            max_states,
            tolerance,
            ..
        } = &request.method;
        if *max_states == 0 {
            out.push(request_error(
                "/method/max_states",
                "max_states must be at least 1",
                "Give the most states to enumerate, at most 100000 are used.",
            ));
        }
        if !(tolerance.is_finite() && *tolerance > 0.0 && *tolerance < 1.0) {
            out.push(request_error(
                "/method/tolerance",
                "tolerance must be in (0, 1)",
                "Give the total series truncation budget, for example 1e-12.",
            ));
        }
        if request.budget.max_output_bytes == 0 {
            out.push(request_error(
                "/budget/max_output_bytes",
                "output budget must be positive",
                "Give an explicit output budget.",
            ));
        }
        if request
            .budget
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

    /// The bound of every species, in view order, or the bound errors.
    fn bounds(
        view: &FlatView,
        bounds: &BTreeMap<String, u64>,
    ) -> Result<Vec<i64>, Vec<Diagnostic>> {
        let mut out = Vec::new();
        for label in bounds.keys() {
            if view.index(label).is_none() {
                let path = pointer("/method/bounds", label);
                let never = view.model().types.iter().any(|t| &t.name == label);
                out.push(if never {
                    request_error(
                        &path,
                        format!("species type `{label}` never holds molecules in this model, so it has no species label"),
                        "Leave out types that have no initial instance and that no rule produces.",
                    )
                } else {
                    request_error(
                        &path,
                        format!("unknown species label `{label}`"),
                        "Bounds name the species labels of the model's flat view: the type name, or Type{attribute=value,...} for an attributed type.",
                    )
                });
            }
        }
        let mut result = Vec::with_capacity(view.species().len());
        for (label, &initial) in view.species().iter().zip(view.initial()) {
            match bounds.get(label) {
                Some(&bound) if bound >= initial => {
                    result.push(i64::try_from(bound).unwrap_or(i64::MAX));
                }
                found => out.push(request_error(
                    &pointer("/method/bounds", label),
                    if found.is_some() {
                        "FSP bound below the initial count"
                    } else {
                        "missing FSP bound"
                    },
                    "Give every species a truncation bound at least its initial count; escaping probability goes to a sink.",
                )),
            }
        }
        if out.is_empty() { Ok(result) } else { Err(out) }
    }

    /// A flat view as the core's generator.
    struct View<'a> {
        view: &'a FlatView,
        evaluator: Evaluator<'a>,
        counts: Vec<u64>,
    }

    impl Generator for View<'_> {
        fn initial(&self) -> Vec<i64> {
            // Counts are at most 2^53.
            self.view.initial().iter().map(|&n| n as i64).collect()
        }
        fn changes(&self) -> Vec<Vec<(usize, i64)>> {
            (0..self.view.channels())
                .map(|c| self.view.changes(c).to_vec())
                .collect()
        }
        fn rates(&mut self, state: &[i64], rates: &mut Vec<f64>) -> Result<bool, Diagnostic> {
            for (out, &n) in self.counts.iter_mut().zip(state) {
                *out = n as u64;
            }
            self.evaluator.load(&self.counts)?;
            rates.clear();
            if self.evaluator.stops()? {
                rates.resize(self.view.channels(), 0.0);
                return Ok(true);
            }
            for c in 0..self.view.channels() {
                rates.push(self.evaluator.channel(c)?);
            }
            Ok(false)
        }
        fn admit(&self, _state: &[i64]) -> Result<(), Diagnostic> {
            Ok(())
        }
    }

    /// Polls a rules [`Control`] and the wall clock: the flag at every poll,
    /// the cancel file and the clock at the first poll and every 1024th.
    struct Poll<'a> {
        control: &'a Control,
        start: Instant,
        wall: Option<f64>,
        countdown: u32,
    }

    impl Poll<'_> {
        fn poll(&mut self) -> Option<String> {
            if self.control.cancel_flagged() {
                return Some("cancelled".into());
            }
            if self.countdown > 0 {
                self.countdown -= 1;
                return None;
            }
            self.countdown = 1023;
            if self.control.is_cancelled() {
                return Some("cancelled".into());
            }
            if self
                .wall
                .is_some_and(|seconds| self.start.elapsed().as_secs_f64() >= seconds)
            {
                return Some("wall_time_limit".into());
            }
            None
        }
    }

    /// Solves an FSP request against the exact model bytes; see the
    /// [module docs](self).
    ///
    /// The model is any closed rules model (`crn-craft.model/1.0` or a
    /// `crn-craft.rules-model` version). States are the model's
    /// [`FlatView`] counts, and every propensity is the one Static Direct
    /// uses in that state. A state where the model's `stop` condition holds
    /// is absorbing. The request's wall-time limit and the [`Control`] are
    /// polled during the computation; a stop fails the solve.
    ///
    /// # Errors
    ///
    /// `CRN-E001` to `CRN-E004` for documents that do not parse; the model's
    /// validation diagnostics; `CRN-R017` for a model without a flat view;
    /// `CRN-R020` for the request (a version other than run/1.0 at
    /// `/schema`, the model binding, the time grid, `max_states`,
    /// `tolerance`, the budget, and at `/method/bounds/LABEL` an unknown
    /// label, a missing bound or a bound below the initial count); `CRN-E070`
    /// at `/fsp` for exceeded state, edge, work, memory or output limits, a
    /// tolerance below 10⁻¹⁸ per step, a transition that underflows, invalid
    /// probabilities or an interruption; and `CRN-R030` for a propensity,
    /// observable or stop condition that cannot be evaluated in a reachable
    /// state.
    pub fn solve(
        model_bytes: &[u8],
        request_bytes: &[u8],
        control: &Control,
    ) -> Result<FspResult, Vec<Diagnostic>> {
        let model = crate::rules::flat_view::parse(model_bytes).map_err(|e| vec![e])?;
        let request: FspRequest = crate::strict_json::parse(request_bytes).map_err(|e| vec![e])?;
        let model_sha256 = sha256(model_bytes);
        let mut diagnostics = validate(&request, &model_sha256);
        let view = match FlatView::new(&model) {
            Ok(view) => Some(view),
            Err(mut d) => {
                diagnostics.append(&mut d);
                None
            }
        };
        if !diagnostics.is_empty() {
            return Err(diagnostics);
        }
        let view = view.expect("a view when valid");
        let FspMethod::Fsp {
            bounds: requested,
            max_states,
            tolerance,
        } = &request.method;
        let bounds = bounds(&view, requested)?;
        let max_output_bytes = request.budget.max_output_bytes.min(MAX_RESULT_BYTES);
        let mut poller = Poll {
            control,
            start: Instant::now(),
            wall: request.budget.wall_time_s,
            countdown: 0,
        };
        let mut poll = || poller.poll();
        let mut work = Work {
            poll: &mut poll,
            used: 0,
        };
        let settings = Settings {
            bounds: &bounds,
            capacity: (*max_states).min(super::MAX_STATES as u64) as usize,
            times: &request.time.report_times,
            tolerance: *tolerance,
            max_output_bytes,
        };
        let mut generator = View {
            view: &view,
            evaluator: view.evaluator(),
            counts: vec![0; view.species().len()],
        };
        let projection = project(&mut generator, &settings, &mut work).map_err(|e| vec![e])?;
        let result = FspResult {
            schema: FspVersion::V10,
            model_sha256,
            run_sha256: sha256(request_bytes),
            species: view.species().to_vec(),
            states: projection
                .states
                .into_iter()
                .map(|s| s.into_iter().map(|n| n as u64).collect())
                .collect(),
            snapshots: projection.snapshots.into_iter().map(Into::into).collect(),
            uniformization_rate: projection.uniformization_rate,
            sparse_edges: projection.sparse_edges,
            work: work.used,
            numerical_error_certified: false,
        };
        if serialized_len(&result) > max_output_bytes {
            return Err(vec![error("FSP result exceeds max_output_bytes")]);
        }
        Ok(result)
    }

    /// Solves a request file against the model it references, which must be
    /// the supplied model file below the request's directory, as
    /// [`crate::rules::solve_files`] requires. `crn-craft rules MODEL
    /// REQUEST` calls it for an FSP request ([`is_fsp_request_file`]).
    ///
    /// # Errors
    ///
    /// `CRN-E030` for a file that cannot be read, `CRN-E004` above 16 MiB,
    /// `CRN-E001` to `CRN-E004` for a request that does not parse, `CRN-R020`
    /// at `/model/path` when the request does not reference `model_path`,
    /// and otherwise the errors of [`solve`].
    pub fn solve_files(
        model_path: &Path,
        request_path: &Path,
        control: &Control,
    ) -> Result<FspResult, Vec<Diagnostic>> {
        let model = crate::rules::read(model_path).map_err(|e| vec![e])?;
        let bytes = crate::rules::read(request_path).map_err(|e| vec![e])?;
        let request: FspRequest = crate::strict_json::parse(&bytes).map_err(|e| vec![e])?;
        if !crate::rules::references_model(model_path, request_path, &request.model.path) {
            return Err(vec![request_error(
                "/model/path",
                "rules request does not identify the supplied model below its directory",
                "Use the referenced model within the request directory.",
            )]);
        }
        solve(&model, &bytes, control)
    }

    /// The `run-1.0` schema: the simulation request ([`crate::rules::RunRequest`],
    /// as `rules::schema` generated it) and [`FspRequest`] as alternatives.
    pub(crate) fn unified_run_schema(mut simulation: Value) -> Option<Value> {
        let fsp = serde_json::to_value(schema_for!(FspRequest)).ok()?;
        let object = simulation.as_object_mut()?;
        let uri = object.remove("$schema")?;
        let Value::Object(mut definitions) = object.remove("definitions")? else {
            return None;
        };
        object.remove("title");
        let Value::Object(mut fsp) = fsp else {
            return None;
        };
        let Some(Value::Object(fsp_definitions)) = fsp.remove("definitions") else {
            return None;
        };
        // Shared definitions (the model reference, the time grid and the
        // filtered request version) keep the simulation request's form.
        for (name, definition) in fsp_definitions {
            definitions.entry(name).or_insert(definition);
        }
        fsp.remove("$schema");
        fsp.remove("title");
        definitions.insert("RunRequest".into(), simulation);
        definitions.insert("FspRequest".into(), Value::Object(fsp));
        Some(serde_json::json!({
            "$schema": uri,
            "title": "Run",
            "description": RUN_1_0,
            "oneOf": [
                {"$ref": "#/definitions/RunRequest"},
                {"$ref": "#/definitions/FspRequest"}
            ],
            "definitions": definitions
        }))
    }

    const RUN_1_0: &str = "A run request: `crn-craft.run/1.0`, the unified request contract. A simulation request (`RunRequest`) has the fields and defaults of `crn-craft.rules-run/0.6` plus `outputs.quantiles` and returns `crn-craft.result/1.0`. A request whose method is `fsp` (`FspRequest`) computes the distribution of a closed model's counts by finite state projection and returns `crn-craft.fsp-result/1.0`.\n\nBoth bind one model, `crn-craft.model/1.0` or any `crn-craft.rules-model` version, by hash. Unknown fields and duplicate keys are rejected.";
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The row-oriented scatter product used before P17 and P18.
    fn reference(rows: &Rows, term: &[f64], next: &mut [f64]) {
        let n = rows.len();
        next.fill(0.0);
        next[n] = term[n];
        for (i, (diagonal, edges)) in rows.iter().enumerate() {
            next[i] += term[i] * diagonal;
            for &(j, rate) in edges {
                next[if j == usize::MAX { n } else { j }] += term[i] * rate;
            }
        }
    }

    fn random_rows(n: usize, seed: &mut u64) -> Vec<(f64, Vec<(usize, f64)>)> {
        let mut uniform = || {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*seed >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n)
            .map(|_| {
                let mut edges: Vec<(usize, f64)> = (0..4)
                    .map(|_| {
                        let j = (uniform() * (n + 1) as f64) as usize;
                        (if j >= n { usize::MAX } else { j }, uniform() * 0.2)
                    })
                    .collect();
                edges.sort_by_key(|&(j, _)| j);
                edges.dedup_by_key(|&mut (j, _)| j);
                (uniform() * 0.2, edges)
            })
            .collect()
    }

    #[test]
    fn the_gathered_product_and_fused_sum_match_the_row_product_bit_for_bit() {
        let mut seed = 20261101;
        for n in [1, 2, 17, 300] {
            let rows = random_rows(n, &mut seed);
            let columns = columns(&rows);
            for density in [0.0, 0.1, 0.5, 1.0] {
                let mut term = vec![0.0; n + 1];
                for (i, value) in term.iter_mut().enumerate() {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let u = (seed >> 11) as f64 / (1u64 << 53) as f64;
                    if u < density {
                        // Include subnormal and tiny values.
                        *value = if i % 7 == 0 {
                            5e-324 * (i + 1) as f64
                        } else {
                            u * 1e-3
                        };
                    }
                }
                let weight = 0.37;
                let (mut fast, mut slow) = (vec![1.0; n + 1], vec![2.0; n + 1]);
                let mut fused: Vec<f64> = (0..=n).map(|i| i as f64 * 1e-4).collect();
                let mut separate = fused.clone();
                multiply_accumulate(&columns, &term, &mut fast, &mut fused, weight);
                reference(&rows, &term, &mut slow);
                for (out, value) in separate.iter_mut().zip(&slow) {
                    *out += weight * value;
                }
                let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(&fast), bits(&slow), "n {n}, density {density}");
                assert_eq!(bits(&fused), bits(&separate), "n {n}, density {density}");
            }
        }
    }
}
