//! Experimental adaptive explicit Poisson leaps with rejection and Direct fallback.
//! Approximation is explicit; no clipping, reflecting boundaries or exact assessment.
//!
//! [`solve`] reads a `crn-craft.tau-request/0.1` [`Request`]: an accuracy
//! `epsilon`, a leap budget and a complete `run/0.1` SSA run as the
//! reference for model binding, time grid, ensemble, outputs and budgets. Each
//! step chooses a leap that bounds the expected relative change of every
//! species, draws Poisson firing counts, halves the step after a proposal
//! that would make a count negative or consume more molecules than present,
//! and takes exact Direct steps when a reactant is nearly exhausted or a
//! step has been halved 20 times. Replicates run in
//! order on one thread and use their own streams ([`crate::sampling::tau_stream`]),
//! so a tau ensemble is not coupled to an exact one with the same seed. The
//! result is a `crn-craft.tau-result/0.1` [`TauResult`] with assessment
//! `approximate_ensemble` at best. Request and numerical failures are
//! `CRN-E085`.
//!
//! ```
//! use crn_craft::{sha256, simulation::{Control, Stop}, tau};
//!
//! let model = br#"{"schema": "crn-craft.model/0.1", "name": "birth_death",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
//!   "species": [{"name": "X", "initial_count": 30}],
//!   "reactions": [
//!     {"name": "make", "reactants": {}, "products": {"X": 1},
//!      "rate": {"kind": "mass_action", "coefficient": 20.0}},
//!     {"name": "decay", "reactants": {"X": 1}, "products": {},
//!      "rate": {"kind": "mass_action", "coefficient": 0.5}}]}"#;
//! let request = format!(r#"{{"schema": "crn-craft.tau-request/0.1", "epsilon": 0.05,
//!   "max_leaps_per_replicate": 10000,
//!   "reference_run": {{"schema": "crn-craft.run/0.1",
//!     "model": {{"path": "model.json", "sha256": "{}"}},
//!     "method": {{"kind": "ssa", "aggregator": "direct"}},
//!     "time": {{"end": 2.0, "save_times": [0.0, 1.0, 2.0]}},
//!     "ensemble": {{"replicates": 8, "root_seed": 91}},
//!     "outputs": {{"statistics": ["mean"], "keep_trajectories": 0}},
//!     "budget": {{"max_events_per_replicate": 100000, "threads": 1,
//!                "max_output_bytes": 100000}}}}}}"#, sha256(model));
//! let result = tau::solve(model, request.as_bytes(), &Control::default()).unwrap();
//! assert_eq!(result.termination, Stop::Completed);
//! // Leaping is approximate: a completed ensemble is `approximate_ensemble`, never exact.
//! assert!(result.assessment.starts_with("approximate"));
//! assert!(result.leaps > 0);
//! ```
use crate::{
    CompiledNetwork,
    contracts::{Run, Statistic},
    diagnostic::Diagnostic,
    sampling,
    simulation::{self, Control, Poller, Replicate, Statistics, Stop, Trajectory},
};
use rand_core::RngCore;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::time::Instant;
// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of a tau-leaping request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum RequestVersion {
    #[doc = concat!("`crn-craft.tau-request/0.1`.")]
    #[serde(rename = "crn-craft.tau-request/0.1")]
    V01,
}
/// Experimental approximate Poisson-leaping request around an exact SSA reference run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    #[doc = concat!("Always `crn-craft.tau-request/0.1`.")]
    pub schema: RequestVersion,
    /// Accuracy parameter in [0.001, 0.2] that bounds the expected relative change of species
    /// counts, and so of rates, within one leap. Smaller is more accurate and slower.
    pub epsilon: f64,
    /// Most leaps and fallback steps per replicate, 1–10000000; exceeding it fails the
    /// replicate.
    pub max_leaps_per_replicate: u64,
    #[doc = concat!("A complete `run/0.1` SSA run supplying the model reference, time grid,")]
    #[doc = concat!("ensemble, outputs and budgets. Its event cap also bounds the events per")]
    #[doc = concat!("replicate that leaps may fire.")]
    pub reference_run: Run,
}
/// Schema tag of a tau-leaping result.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ResultVersion {
    #[doc = concat!("`crn-craft.tau-result/0.1`.")]
    #[serde(rename = "crn-craft.tau-result/0.1")]
    V01,
}
/// Approximate Poisson-leaping result. It is never an exact ensemble.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TauResult {
    #[doc = concat!("Always `crn-craft.tau-result/0.1`.")]
    pub schema: ResultVersion,
    /// Algorithm identity, including its version.
    pub algorithm: String,
    /// `approximate_ensemble` when every replicate completed, otherwise `incomplete`.
    pub assessment: String,
    /// SHA-256 of the exact model bytes.
    pub model_sha256: String,
    /// SHA-256 of the exact tau request bytes.
    pub request_sha256: String,
    /// The request's accuracy parameter.
    pub epsilon: f64,
    #[doc = concat!("Aggregate stop: the highest-priority replicate stop.")]
    pub termination: Stop,
    /// Species names in the order used by every count row.
    pub species: Vec<String>,
    /// The reference run's save grid in seconds.
    pub save_times: Vec<f64>,
    /// Replicates the reference run requested.
    pub requested_replicates: u64,
    /// Reaction firings over all replicates, leaps and exact steps together.
    pub events: u64,
    /// Accepted leaps and exact fallback steps over all replicates.
    pub leaps: u64,
    /// Leap proposals rejected because they would make a count negative or consume more
    /// molecules than present.
    pub rejected_leaps: u64,
    /// One record per requested replicate, in id order.
    pub replicates: Vec<Replicate>,
    /// Saved states of the first `keep_trajectories` replicates.
    pub trajectories: Vec<Trajectory>,
    /// Requested statistics, present only when every replicate completed.
    pub statistics: Option<Statistics>,
}
fn error(s: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E085",
        "/tau",
        s,
        "Use epsilon in [0.001,0.2], finite budgets and non-negative checked populations; this method is approximate.",
    )
}
/// Largest Poisson mean a single leap may draw for one channel.
const POISSON_MEAN_LIMIT: f64 = 1e7;
/// Algorithm identity recorded in [`TauResult::algorithm`].
pub const ALGORITHM: &str = "adaptive_explicit_poisson_cgp_reject_halve_direct_fallback/0.3";

/// Cao, Gillespie and Petzold (2006), eq. 27: the factor g_i for a species
/// whose highest-order reaction as a reactant has total `order` and needs
/// `multiplicity` molecules of it. Orders above three use the order itself.
pub(crate) fn cgp_g(order: i64, multiplicity: i64, count: i64) -> f64 {
    let x = count as f64;
    match (order, multiplicity) {
        (..=1, _) => 1.0,
        (2, 2) => 2.0 + 1.0 / (x - 1.0),
        (2, _) => 2.0,
        (3, 3) => 3.0 + 1.0 / (x - 1.0) + 2.0 / (x - 2.0),
        (3, 2) => 1.5 * (2.0 + 1.0 / (x - 1.0)),
        (3, _) => 3.0,
        _ => order as f64,
    }
}

/// Per species: the highest total reactant order of a reaction consuming it,
/// and the largest multiplicity of the species among reactions of that order.
fn highest_orders(network: &CompiledNetwork) -> Vec<(i64, i64)> {
    let mut orders = vec![(0, 0); network.species().len()];
    for reaction in network.reactions() {
        let order: i64 = reaction.reactants.iter().map(|t| t.coefficient).sum();
        for term in &reaction.reactants {
            let entry = &mut orders[term.species];
            match order.cmp(&entry.0) {
                std::cmp::Ordering::Greater => *entry = (order, term.coefficient),
                std::cmp::Ordering::Equal => entry.1 = entry.1.max(term.coefficient),
                std::cmp::Ordering::Less => {}
            }
        }
    }
    orders
}

/// Leap means above this use transformed rejection; at or below it the
/// multiplication method draws one Poisson variate directly.
const PTRS_MEAN: f64 = 16.0;

/// A Poisson variate. Means up to 16 use the multiplication method
/// (product of uniforms against exp(-mean)). Larger means use Hörmann's
/// transformed rejection with squeeze, PTRS (Hörmann 1993, "The transformed
/// rejection method for generating Poisson random variables", Insurance:
/// Mathematics and Economics 12, 39–45), whose expected number of uniform
/// pairs is about 1.1 for every mean, instead of about `mean` uniforms.
fn poisson(mean: f64, rng: &mut impl RngCore, poller: &mut Poller) -> Result<u64, Diagnostic> {
    if !mean.is_finite() || !(0.0..=POISSON_MEAN_LIMIT).contains(&mean) {
        return Err(error("Poisson mean exceeds supported work range"));
    }
    if mean == 0. {
        return Ok(0);
    }
    if let Some(stop) = poller.poll() {
        return Err(error(format!("interrupted Poisson sampling: {stop:?}")));
    }
    poisson_with(mean, || Ok(sampling::uniform_open(rng)))
}

/// Shared Poisson kernel. Callers validate the mean and supply open uniforms
/// from their own budgeted sampler.
#[doc(hidden)]
pub fn poisson_with<E>(mean: f64, mut open: impl FnMut() -> Result<f64, E>) -> Result<u64, E> {
    if mean == 0.0 {
        return Ok(0);
    }
    if mean > PTRS_MEAN {
        return ptrs(mean, open);
    }
    let threshold = libm::exp(-mean);
    let mut count = 0;
    let mut product = open()?;
    while product > threshold {
        count += 1;
        product *= open()?;
    }
    Ok(count)
}

/// Hörmann's PTRS for means of at least 10, with his published constants.
/// Uniforms are drawn from (0, 1), so `us` is positive and every logarithm is
/// finite. Each pair of uniforms is accepted by the squeeze, rejected, or
/// accepted by the exact log-density comparison.
fn ptrs<E>(mean: f64, mut open: impl FnMut() -> Result<f64, E>) -> Result<u64, E> {
    let root = mean.sqrt();
    let b = 0.931 + 2.53 * root;
    let a = -0.059 + 0.02483 * b;
    let inverse_alpha = 1.1239 + 1.1328 / (b - 3.4);
    let v_r = 0.9277 - 3.6224 / (b - 2.0);
    let log_mean = mean.ln();
    loop {
        let u = open()? - 0.5;
        let v = open()?;
        let us = 0.5 - u.abs();
        let k = ((2.0 * a / us + b) * u + mean + 0.43).floor();
        if us >= 0.07 && v <= v_r {
            return Ok(k as u64);
        }
        if k < 0.0 || (us < 0.013 && v > us) {
            continue;
        }
        if v.ln() + inverse_alpha.ln() - (a / (us * us) + b).ln()
            <= -mean + k * log_mean - libm::lgamma(k + 1.0)
        {
            return Ok(k as u64);
        }
    }
}
fn path(
    network: &CompiledNetwork,
    run: &Run,
    epsilon: f64,
    max_leaps: u64,
    id: u64,
    poller: &mut Poller,
) -> (Replicate, Vec<Vec<i64>>, u64, u64) {
    let mut record = Replicate {
        id,
        started: false,
        termination: Stop::Completed,
        events: 0,
        stop_time: 0.,
        saved_points: 0,
        diagnostic: None,
    };
    let mut rng = sampling::tau_stream(run.ensemble.as_ref().unwrap().root_seed, id);
    let orders = highest_orders(network);
    let mut state = network.initial_counts().to_vec();
    let mut counts = Vec::new();
    let mut time = 0.;
    let mut leaps = 0;
    let mut rejected = 0;
    let result = (|| -> Result<Stop, Diagnostic> {
        if let Some(stop) = poller.poll() {
            return Ok(stop);
        }
        record.started = true;
        for &target in &run.time.save_times {
            while time < target {
                if let Some(stop) = poller.poll() {
                    return Ok(stop);
                }
                if record.events >= run.budget.max_events_per_replicate.unwrap() {
                    return Ok(Stop::EventLimit);
                }
                if leaps >= max_leaps {
                    return Err(error("leap budget exceeded"));
                }
                let rates = network.propensities(&state)?;
                let total = sampling::total_propensity(&rates)?;
                if total == 0. {
                    time = target;
                    break;
                }
                let critical = network.reactions().iter().zip(&rates).any(|(rx, &a)| {
                    a > 0.
                        && rx
                            .reactants
                            .iter()
                            .any(|term| state[term.species] / term.coefficient <= 10)
                });
                if critical {
                    let next = sampling::event_time(
                        time,
                        sampling::waiting_time(total, sampling::uniform_open(&mut rng))?,
                    )?;
                    if next > target {
                        time = target;
                        break;
                    }
                    let channel = sampling::select_channel(
                        &rates,
                        total,
                        sampling::uniform_closed_open(rng.next_u64()),
                    )?;
                    network.fire_reaction(channel, &mut state)?;
                    record.events += 1;
                    time = next;
                    leaps += 1;
                    continue;
                }
                let mut drift = vec![0.; state.len()];
                let mut variance = vec![0.; state.len()];
                for (rx, &a) in network.reactions().iter().zip(&rates) {
                    for term in &rx.net_changes {
                        let c = term.coefficient as f64;
                        drift[term.species] += c * a;
                        variance[term.species] += c * c * a;
                    }
                }
                let mut dt = target - time;
                for i in 0..state.len() {
                    // No critical reactant here, so every reactant count exceeds ten.
                    let (order, multiplicity) = orders[i];
                    let g = cgp_g(order, multiplicity, state[i]);
                    let scale = (epsilon * state[i] as f64 / g).max(1.);
                    if drift[i] != 0. {
                        dt = dt.min(scale / drift[i].abs());
                    }
                    if variance[i] > 0. {
                        dt = dt.min(scale * scale / variance[i]);
                    }
                }
                // Keep every channel's leap mean within the sampler's range.
                let largest = rates.iter().copied().fold(0.0, f64::max);
                dt = dt.min(POISSON_MEAN_LIMIT / largest);
                for attempt in 0..=20 {
                    if let Some(stop) = poller.poll() {
                        return Ok(stop);
                    }
                    if !dt.is_finite() || dt <= 0. || time + dt <= time {
                        return Err(error("tau clock cannot advance"));
                    }
                    if attempt == 20 {
                        let next = sampling::event_time(
                            time,
                            sampling::waiting_time(total, sampling::uniform_open(&mut rng))?,
                        )?;
                        if next > target {
                            time = target;
                            break;
                        }
                        let channel = sampling::select_channel(
                            &rates,
                            total,
                            sampling::uniform_closed_open(rng.next_u64()),
                        )?;
                        network.fire_reaction(channel, &mut state)?;
                        record.events += 1;
                        time = next;
                        leaps += 1;
                        break;
                    }
                    let mut proposed: Vec<i128> = state.iter().map(|x| i128::from(*x)).collect();
                    let mut events = 0_u64;
                    let mut consumption = vec![0_i128; state.len()];
                    for (r, &a) in rates.iter().enumerate() {
                        let k = poisson(a * dt, &mut rng, poller)?;
                        events = events
                            .checked_add(k)
                            .ok_or_else(|| error("tau event count overflow"))?;
                        for term in &network.reactions()[r].reactants {
                            consumption[term.species] +=
                                i128::from(term.coefficient) * i128::from(k);
                        }
                        for term in &network.reactions()[r].net_changes {
                            proposed[term.species] += i128::from(term.coefficient) * i128::from(k);
                        }
                    }
                    if proposed.iter().any(|x| *x < 0)
                        || consumption
                            .iter()
                            .zip(&state)
                            .any(|(used, available)| *used > i128::from(*available))
                    {
                        rejected += 1;
                        dt *= 0.5;
                        continue;
                    }
                    if proposed.iter().any(|x| *x > i128::from(i64::MAX)) {
                        return Err(error("tau population overflow"));
                    }
                    let proposed: Vec<_> = proposed.into_iter().map(|x| x as i64).collect();
                    if let Err(e) = network.validate_state(&proposed) {
                        record.diagnostic = Some(e);
                        return Ok(Stop::PopulationLimit);
                    }
                    if events > run.budget.max_events_per_replicate.unwrap() - record.events {
                        return Ok(Stop::EventLimit);
                    }
                    state = proposed;
                    record.events += events;
                    leaps += 1;
                    time = (time + dt).min(target);
                    break;
                }
            }
            counts.push(state.clone());
            record.stop_time = target;
        }
        Ok(Stop::Completed)
    })();
    match result {
        Ok(stop) => record.termination = stop,
        Err(e) => {
            record.termination = poller.poll().unwrap_or_else(|| {
                if e.code == "CRN-E044" && e.message.contains("safety limit") {
                    Stop::PopulationLimit
                } else {
                    Stop::Failed
                }
            });
            record.diagnostic = Some(e);
        }
    }
    record.stop_time = record.stop_time.max(time);
    record.saved_points = counts.len();
    (record, counts, leaps, rejected)
}
/// Runs a tau-leaping request against the exact model bytes; see the
/// [module docs](self). `reference_run.model.path` is only syntax-checked;
/// use [`solve_files`] to check that it names the model file.
///
/// # Errors
///
/// Strict-parse errors; `CRN-E085` for `epsilon` outside [0.001, 0.2], a
/// leap budget outside 1–10000000 or a result above `max_output_bytes`; the
/// preparation errors of [`crate::simulation::prepare`] for the reference
/// run. Stops and numerical failures during a replicate are reported in its
/// record.
pub fn solve(
    model: &[u8],
    request_bytes: &[u8],
    control: &Control,
) -> Result<TauResult, Vec<Diagnostic>> {
    let request: Request =
        crate::parse_network(request_bytes, "`crn-craft.tau-request/0.1` requests")
            .map_err(|e| vec![e])?;
    if !request.epsilon.is_finite()
        || !(0.001..=0.2).contains(&request.epsilon)
        || request.max_leaps_per_replicate == 0
        || request.max_leaps_per_replicate > 10_000_000
    {
        return Err(vec![error("invalid epsilon or leap budget")]);
    }
    let bytes = serde_json::to_vec(&request.reference_run)
        .map_err(|_| vec![error("invalid reference run")])?;
    let prepared = simulation::prepare(model, &bytes)?;
    let run = &prepared.run;
    let ensemble = run.ensemble.as_ref().unwrap();
    let cells = prepared.network.species().len() * run.time.save_times.len();
    let n = prepared.network.species().len();
    let mut mean = vec![0.; cells];
    let mut m2 = vec![0.; cells];
    let mut records = Vec::new();
    let mut kept = Vec::new();
    let mut total_leaps = 0;
    let mut rejected = 0;
    let mut poller = Poller::new(control, Instant::now(), run.budget.wall_time_s);
    for id in 0..ensemble.replicates {
        let (record, counts, leaps, rejections) = path(
            &prepared.network,
            run,
            request.epsilon,
            request.max_leaps_per_replicate,
            id,
            &mut poller,
        );
        total_leaps += leaps;
        rejected += rejections;
        if record.termination == Stop::Completed {
            for (i, &x) in counts.iter().flatten().enumerate() {
                let delta = x as f64 - mean[i];
                mean[i] += delta / (id + 1) as f64;
                m2[i] += delta * (x as f64 - mean[i]);
            }
        }
        if id < run.outputs.keep_trajectories {
            kept.push(Trajectory { id, counts });
        }
        records.push(record);
    }
    let termination = records
        .iter()
        .map(|r| r.termination)
        .max()
        .unwrap_or(Stop::Failed);
    let chunks = |x: Vec<f64>| x.chunks(n).map(<[f64]>::to_vec).collect();
    let variance: Vec<_> = m2
        .iter()
        .map(|x| x / (ensemble.replicates - 1).max(1) as f64)
        .collect();
    let statistics =
        (termination == Stop::Completed && !run.outputs.statistics.is_empty()).then(|| {
            Statistics {
                samples: ensemble.replicates,
                mean: run
                    .outputs
                    .statistics
                    .contains(&Statistic::Mean)
                    .then(|| chunks(mean)),
                sample_variance: run
                    .outputs
                    .statistics
                    .contains(&Statistic::Variance)
                    .then(|| chunks(variance.clone())),
                mean_standard_error: (run.outputs.statistics.contains(&Statistic::Mean)
                    && ensemble.replicates > 1)
                    .then(|| {
                        chunks(
                            variance
                                .iter()
                                .map(|v| libm::sqrt(v / ensemble.replicates as f64))
                                .collect(),
                        )
                    }),
            }
        });
    let result = TauResult {
        schema: ResultVersion::V01,
        algorithm: ALGORITHM.into(),
        assessment: if termination == Stop::Completed {
            "approximate_ensemble"
        } else {
            "incomplete"
        }
        .into(),
        model_sha256: crate::sha256(model),
        request_sha256: crate::sha256(request_bytes),
        epsilon: request.epsilon,
        termination,
        species: prepared.network.species().to_vec(),
        save_times: run.time.save_times.clone(),
        requested_replicates: ensemble.replicates,
        events: records.iter().map(|r| r.events).sum(),
        leaps: total_leaps,
        rejected_leaps: rejected,
        replicates: records,
        trajectories: kept,
        statistics,
    };
    let size = serde_json::to_vec(&result)
        .map_err(|_| vec![error("non-finite tau output")])?
        .len() as u64;
    if size + 1
        > run
            .budget
            .max_output_bytes
            .min(simulation::MAX_RESULT_BYTES)
    {
        return Err(vec![error("tau output byte budget exceeded")]);
    }
    Ok(result)
}

/// File API checks that the embedded reference identifies the supplied model.
///
/// `reference_run.model.path` is resolved from the request's directory and
/// must name `model_path` inside it (`CRN-E022` at
/// `/reference_run/model/path`); then [`solve`] runs.
///
/// # Errors
///
/// `CRN-E030` or `CRN-E004` for unreadable or oversized files, `CRN-E022`,
/// and the errors of [`solve`].
pub fn solve_files(
    model_path: &std::path::Path,
    request_path: &std::path::Path,
    control: &Control,
) -> Result<TauResult, Vec<Diagnostic>> {
    let model = crate::read_document(model_path).map_err(|e| vec![e])?;
    let bytes = crate::read_document(request_path).map_err(|e| vec![e])?;
    let request: Request = crate::parse_network(&bytes, "`crn-craft.tau-request/0.1` requests")
        .map_err(|e| vec![e])?;
    let parent = request_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(std::path::Path::new("."));
    let resolved = (
        parent.canonicalize(),
        parent
            .join(&request.reference_run.model.path)
            .canonicalize(),
        model_path.canonicalize(),
    );
    if !crate::validate::valid_relative_path(&request.reference_run.model.path)
        || !matches!(resolved,(Ok(ref base),Ok(ref reference),Ok(ref supplied)) if reference.starts_with(base) && reference==supplied)
    {
        return Err(vec![Diagnostic::error(
            "CRN-E022",
            "/reference_run/model/path",
            "tau reference does not identify the supplied model below the request directory",
            "Use the referenced model within the request directory.",
        )]);
    }
    solve(&model, &bytes, control)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cgp_factors_follow_the_published_cases() {
        assert_eq!(cgp_g(1, 1, 50), 1.0);
        assert_eq!(cgp_g(2, 1, 50), 2.0);
        assert_eq!(cgp_g(2, 2, 11), 2.0 + 1.0 / 10.0);
        assert_eq!(cgp_g(3, 1, 50), 3.0);
        assert_eq!(cgp_g(3, 2, 11), 1.5 * (2.0 + 1.0 / 10.0));
        assert_eq!(cgp_g(3, 3, 12), 3.0 + 1.0 / 11.0 + 2.0 / 10.0);
        assert_eq!(cgp_g(4, 1, 50), 4.0);
        assert_eq!(cgp_g(0, 0, 50), 1.0);
    }

    #[test]
    fn tau_keys_differ_from_direct_keys() {
        for seed in [0, 1, 4242] {
            for replicate in [0, 7] {
                assert_ne!(
                    sampling::tau_stream(seed, replicate).next_u64(),
                    sampling::stream(seed, replicate).next_u64()
                );
            }
        }
    }

    /// Declared in docs/campaign/performance-p11-p20.md (P19, part a): 200,000
    /// draws at each of five means from fixed streams, a two-sided DKW band at
    /// family alpha 0.001 (Bonferroni over the means) against the exact
    /// Poisson CDF, and every sample mean within five standard errors.
    #[test]
    fn poisson_draws_follow_the_poisson_law_on_both_sides_of_the_switch() {
        let control = Control::default();
        let mut poller = Poller::new(&control, Instant::now(), None);
        let draws = 200_000u32;
        let alpha = 0.001 / 5.0;
        let band = ((2.0f64 / alpha).ln() / (2.0 * f64::from(draws))).sqrt();
        for (stream, mean) in [16.5f64, 30.0, 100.0, 1000.0, 1e5].into_iter().enumerate() {
            let mut rng = sampling::tau_stream(20261101, stream as u64);
            let mut counts = std::collections::BTreeMap::<u64, u32>::new();
            let mut sum = 0.0;
            for _ in 0..draws {
                let k = poisson(mean, &mut rng, &mut poller).unwrap();
                *counts.entry(k).or_default() += 1;
                sum += k as f64;
            }
            let sample_mean = sum / f64::from(draws);
            let standard_error = (mean / f64::from(draws)).sqrt();
            assert!(
                (sample_mean - mean).abs() <= 5.0 * standard_error,
                "mean {mean}: sample mean {sample_mean}"
            );
            // Exact CDF from the log pmf; mass beyond 40 standard deviations is
            // below 1e-300 and irrelevant at this band.
            let spread = 40.0 * mean.sqrt();
            let low = (mean - spread).max(0.0).floor() as u64;
            let high = (mean + spread).ceil() as u64;
            let (first, last) = (
                *counts.keys().next().unwrap(),
                *counts.keys().next_back().unwrap(),
            );
            assert!(
                first >= low && last <= high,
                "draws outside 40 sd at {mean}"
            );
            let (mut cdf, mut empirical, mut worst) = (0.0f64, 0u64, 0.0f64);
            for k in low..=high {
                let x = k as f64;
                cdf += (-mean + x * mean.ln() - libm::lgamma(x + 1.0)).exp();
                empirical += u64::from(counts.get(&k).copied().unwrap_or(0));
                let gap = (empirical as f64 / f64::from(draws) - cdf).abs();
                worst = worst.max(gap);
            }
            assert!(worst <= band, "mean {mean}: DKW distance {worst} > {band}");
        }
    }

    #[test]
    fn small_means_keep_the_multiplication_method_and_its_words() {
        let control = Control::default();
        let mut poller = Poller::new(&control, Instant::now(), None);
        for mean in [0.5, 3.0, 16.0] {
            let (mut a, mut b) = (sampling::tau_stream(7, 1), sampling::tau_stream(7, 1));
            for _ in 0..1000 {
                let k = poisson(mean, &mut a, &mut poller).unwrap();
                // The pre-P19 single-chunk multiplication method.
                let threshold = libm::exp(-mean);
                let mut expected = 0;
                let mut product = sampling::uniform_open(&mut b);
                while product > threshold {
                    expected += 1;
                    product *= sampling::uniform_open(&mut b);
                }
                assert_eq!(k, expected);
            }
            assert_eq!(a.next_u64(), b.next_u64(), "same words consumed at {mean}");
        }
    }

    /// P19's point: a large leap mean costs a few uniforms, not about `mean`.
    #[test]
    fn large_means_draw_a_few_uniforms_per_variate() {
        struct Counting<R>(R, u64);
        impl<R: RngCore> RngCore for Counting<R> {
            fn next_u32(&mut self) -> u32 {
                self.1 += 1;
                self.0.next_u32()
            }
            fn next_u64(&mut self) -> u64 {
                self.1 += 1;
                self.0.next_u64()
            }
            fn fill_bytes(&mut self, bytes: &mut [u8]) {
                self.1 += 1;
                self.0.fill_bytes(bytes)
            }
            fn try_fill_bytes(&mut self, bytes: &mut [u8]) -> Result<(), rand_core::Error> {
                self.1 += 1;
                self.0.try_fill_bytes(bytes)
            }
        }
        let control = Control::default();
        let mut poller = Poller::new(&control, Instant::now(), None);
        for mean in [17.0, 1e3, 1e5, 1e7] {
            let mut rng = Counting(sampling::tau_stream(3, 0), 0);
            for _ in 0..10_000 {
                poisson(mean, &mut rng, &mut poller).unwrap();
            }
            assert!(
                rng.1 <= 30_000,
                "{} words for 10,000 draws at {mean}",
                rng.1
            );
        }
    }
}
