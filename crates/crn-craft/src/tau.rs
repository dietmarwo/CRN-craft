//! Experimental adaptive explicit Poisson leaps with rejection and Direct fallback.
//! Approximation is explicit; no clipping, reflecting boundaries or exact assessment.
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
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum RequestVersion {
    #[serde(rename = "crn-craft.tau-request/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub schema: RequestVersion,
    pub epsilon: f64,
    pub max_leaps_per_replicate: u64,
    pub reference_run: Run,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.tau-result/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TauResult {
    pub schema: ResultVersion,
    pub algorithm: String,
    pub assessment: String,
    pub model_sha256: String,
    pub request_sha256: String,
    pub epsilon: f64,
    pub termination: Stop,
    pub species: Vec<String>,
    pub save_times: Vec<f64>,
    pub requested_replicates: u64,
    pub events: u64,
    pub leaps: u64,
    pub rejected_leaps: u64,
    pub replicates: Vec<Replicate>,
    pub trajectories: Vec<Trajectory>,
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
pub const ALGORITHM: &str = "adaptive_explicit_poisson_cgp_reject_halve_direct_fallback/0.2";

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
            if order > entry.0 {
                *entry = (order, term.coefficient);
            } else if order == entry.0 {
                entry.1 = entry.1.max(term.coefficient);
            }
        }
    }
    orders
}

fn poisson(mean: f64, rng: &mut impl RngCore, poller: &mut Poller) -> Result<u64, Diagnostic> {
    if !mean.is_finite() || !(0.0..=POISSON_MEAN_LIMIT).contains(&mean) {
        return Err(error("Poisson mean exceeds supported work range"));
    }
    if mean == 0. {
        return Ok(0);
    }
    let chunks = (mean / 16.).ceil() as u64;
    let lambda = mean / chunks as f64;
    let threshold = libm::exp(-lambda);
    let mut count = 0;
    for _ in 0..chunks {
        if let Some(stop) = poller.poll() {
            return Err(error(format!("interrupted Poisson sampling: {stop:?}")));
        }
        let mut product = sampling::uniform_open(rng);
        while product > threshold {
            count += 1;
            product *= sampling::uniform_open(rng);
        }
    }
    Ok(count)
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
pub fn solve(
    model: &[u8],
    request_bytes: &[u8],
    control: &Control,
) -> Result<TauResult, Vec<Diagnostic>> {
    let request: Request = crate::strict_json::parse(request_bytes).map_err(|e| vec![e])?;
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
pub fn solve_files(
    model_path: &std::path::Path,
    request_path: &std::path::Path,
    control: &Control,
) -> Result<TauResult, Vec<Diagnostic>> {
    let model = crate::read_document(model_path).map_err(|e| vec![e])?;
    let bytes = crate::read_document(request_path).map_err(|e| vec![e])?;
    let request: Request = crate::strict_json::parse(&bytes).map_err(|e| vec![e])?;
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
}
