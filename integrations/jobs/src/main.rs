#![forbid(unsafe_code)]
//! Register fallible CRN applications in the existing fcmaes worker, no new protocol.
//!
//! Each candidate binds its parameter values to a prepared rules run
//! ([`rules::Prepared::bind_parameters`]) and runs it with Static Direct. The
//! model input is a `crn-craft.model/1.0` document (convert a network model
//! with `crn-craft convert`) or any rules model; the fit-data columns follow
//! the result's observables, sorted by name.
use crn_craft::rules::{self, Control, PathSeries, Prepared, RulesResult, Termination};
use fcmaes_job_worker::{
    Adapter,
    registry::{Example, Instance, Registry, Scalar},
    simulation::SimContext,
};
use optimization_job_api::Result;
use optimization_job_api::worker::Worker;
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    collections::BTreeMap,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
    },
    time::Duration,
};

#[derive(Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
enum DataVersion {
    #[serde(rename = "crn-craft.fit-data/0.1")]
    V01,
}
#[derive(Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Data {
    schema: DataVersion,
    times: Vec<f64>,
    mean: Vec<Vec<f64>>,
    variance: Vec<Vec<f64>>,
}

/// The run of one seed: Static Direct over the data's report times, every
/// path kept for the circuit score. `threads` is the evaluation's thread
/// allowance; budgets are enforced by the job context, not by the request.
fn request(
    model: &[u8],
    times: &[f64],
    replicates: u64,
    seed: u64,
    keep: bool,
    max_events: u64,
    threads: usize,
) -> Vec<u8> {
    serde_json::to_vec(&json!({"schema":"crn-craft.run/1.0",
        "model":{"path":"model.json","sha256":crn_craft::sha256(model)},
        "method":{"kind":"static_direct"},
        "ensemble":{"replicates":replicates,"root_seed":seed},
        "time":{"end":times.last(),"report_times":times},
        "outputs":{"statistics":["mean","variance"],"keep_paths":if keep {replicates} else {0},"keep_instances":0},
        "budget":{"threads":threads,"max_events_per_replicate":max_events,
                  "max_sampling_steps_per_replicate":u64::MAX,"max_instances":u64::MAX,
                  "max_output_bytes":268435456}}))
    .expect("serializable request")
}

/// The prepared runs of every seed for each thread allowance, built once.
struct Templates {
    model: Vec<u8>,
    times: Vec<f64>,
    replicates: u64,
    seeds: Vec<u64>,
    keep: bool,
    max_events: u64,
    by_threads: Mutex<BTreeMap<usize, Arc<Vec<Prepared>>>>,
}
impl Templates {
    fn get(&self, threads: usize) -> std::result::Result<Arc<Vec<Prepared>>, String> {
        let mut cache = self.by_threads.lock().map_err(|e| e.to_string())?;
        if let Some(found) = cache.get(&threads) {
            return Ok(Arc::clone(found));
        }
        let prepared = self
            .seeds
            .iter()
            .map(|&seed| {
                let run = request(
                    &self.model,
                    &self.times,
                    self.replicates,
                    seed,
                    self.keep,
                    self.max_events,
                    threads,
                );
                rules::prepare(&self.model, &run).map_err(|e| serde_json::to_string(&e).unwrap())
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let prepared = Arc::new(prepared);
        cache.insert(threads, Arc::clone(&prepared));
        Ok(prepared)
    }
}

/// Executes a bound candidate, cancelling it when the job is cancelled or
/// reaches its deadline.
fn execute(prepared: &Prepared, ctx: &SimContext) -> std::result::Result<RulesResult, String> {
    let control = Control::default();
    let done = AtomicBool::new(false);
    std::thread::scope(|scope| {
        scope.spawn(|| {
            while !done.load(Ordering::Acquire) {
                if ctx.is_cancelled() {
                    control.cancel();
                    break;
                }
                std::thread::sleep(Duration::from_millis(2));
            }
        });
        let result = rules::execute(prepared, &control);
        done.store(true, Ordering::Release);
        result.map_err(|e| e.message)
    })
}

fn simulation(
    templates: &Templates,
    x: &[f64],
    names: &[String],
    ctx: &SimContext,
) -> std::result::Result<Vec<RulesResult>, String> {
    let mut values = BTreeMap::new();
    for (name, log_rate) in names.iter().zip(x) {
        let rate = log_rate.exp();
        if !rate.is_finite() || rate <= 0.0 {
            return Err("invalid log-rate".into());
        }
        values.insert(name.clone(), rate);
    }
    let prepared = templates.get(ctx.threads().saturating_sub(1).clamp(1, 64))?;
    let mut results = Vec::new();
    for template in prepared.iter() {
        ctx.checkpoint()?;
        let candidate = template
            .bind_parameters(&values)
            .map_err(|e| serde_json::to_string(&e).unwrap())?;
        let result = execute(&candidate, ctx)?;
        if result.termination != Termination::Completed {
            return Err(format!(
                "incomplete CRN objective: {:?}",
                result.termination
            ));
        }
        results.push(result);
    }
    Ok(results)
}
/// Mean over kept paths of 1 / (1 + variance of the amplitude observable over report times
/// after `burn_in`, which start at index `first`). A mean, so the penalty does not grow with
/// the replicate count.
fn amplitude_penalty(
    paths: &[PathSeries],
    observable: usize,
    first: usize,
    replicates: u64,
) -> std::result::Result<f64, String> {
    if paths.len() as u64 != replicates || paths.is_empty() {
        return Err("circuit score requires every replicate path".into());
    }
    let mut total = 0.0;
    for path in paths {
        let x: Vec<_> = path.values[first..].iter().map(|s| s[observable]).collect();
        let mean = x.iter().sum::<f64>() / x.len() as f64;
        let variance = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64;
        total += 1.0 / (1.0 + variance);
    }
    Ok(total / paths.len() as f64)
}
fn application(instance: &Instance, circuit: bool) -> Result<Scalar> {
    let bytes = instance.input("model")?;
    let network = serde_json::from_slice::<serde_json::Value>(&bytes)
        .is_ok_and(|v| v["schema"] == "crn-craft.model/0.1");
    if network {
        return Err("crn-craft.model/0.1 network models are not accepted: convert the model with `crn-craft convert` (crn-craft.model/1.0) and order the fit-data columns like the result's observables, sorted by name".into());
    }
    let model: rules::Model = crn_craft::strict_json::parse(&bytes).map_err(|e| e.message)?;
    let names: Vec<String> = instance.parameter(
        "parameter_names",
        model.parameters.keys().cloned().collect::<Vec<_>>(),
    )?;
    instance.require_dimension(names.len()..=names.len())?;
    if names
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
        != names.len()
        || names
            .iter()
            .any(|n| model.parameters.get(n).is_some_and(|v| *v <= 0.0))
        || names.is_empty()
        || names.iter().any(|n| !model.parameters.contains_key(n))
    {
        return Err("parameter_names must name positive model parameters".into());
    }
    // Fixed seeds across candidates, always explicit: search and held-out jobs must name their
    // own seed sets. SimContext.seed is point-derived and deliberately unused.
    let seeds: Vec<u64> = instance.required("seed_set")?;
    let replicates: u64 = instance.parameter("replicates", 64_u64)?;
    let max_events: u64 = instance.parameter("max_events", 1_000_000_u64)?;
    if seeds
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
        != seeds.len()
        || seeds.is_empty()
        || seeds.len() > 64
        || !(2..=10000).contains(&replicates)
        || max_events == 0
    {
        return Err("invalid seeds/replicates/event budget".into());
    }
    let data: Data =
        crn_craft::strict_json::parse(&instance.input("data")?).map_err(|e| e.message)?;
    let _schema = &data.schema;
    if data.times.len() < 2
        || data.times[0] != 0.0
        || data
            .times
            .windows(2)
            .any(|w| !w[1].is_finite() || w[1] <= w[0])
    {
        return Err("invalid fit-data dimensions/domain/grid".into());
    }
    let templates = Templates {
        model: bytes,
        times: data.times.clone(),
        replicates,
        seeds,
        keep: circuit,
        max_events,
        by_threads: Mutex::new(BTreeMap::new()),
    };
    // The model and the runs are checked here, once; the observables (result
    // columns) come from a cancelled run, which simulates nothing.
    let prepared = templates.get(1)?;
    let cancelled = Control::default();
    cancelled.cancel();
    let observables = rules::execute(&prepared[0], &cancelled)
        .map_err(|e| e.message)?
        .observables;
    if data.mean.len() != data.times.len()
        || data.variance.len() != data.times.len()
        || data
            .mean
            .iter()
            .chain(&data.variance)
            .any(|r| r.len() != observables.len() || r.iter().any(|v| !v.is_finite() || *v < 0.0))
    {
        return Err("invalid fit-data dimensions/domain/grid".into());
    }
    // Circuit amplitude: a named observable, over report times strictly after a declared
    // burn-in, so the deterministic initial state never contributes.
    let amplitude = if circuit {
        let name: String = instance.required("amplitude_species")?;
        let burn_in: f64 = instance.parameter("burn_in", 0.0)?;
        let observable = observables
            .iter()
            .position(|o| *o == name)
            .ok_or("amplitude_species must name a model observable")?;
        let first = data.times.partition_point(|t| *t <= burn_in);
        if !burn_in.is_finite() || burn_in < 0.0 || data.times.len() - first < 2 {
            return Err("burn_in must leave at least two save times after it".into());
        }
        Some((observable, first))
    } else {
        None
    };
    let lower = vec![-12.0; names.len()];
    let upper = vec![12.0; names.len()];
    let objective = move |x: &[f64], ctx: &SimContext| -> std::result::Result<f64, String> {
        let results = simulation(&templates, x, &names, ctx)?;
        let mut loss = 0.0;
        for result in &results {
            let stats = result.statistics.as_ref().ok_or("statistics missing")?;
            let (mean, variance) = stats
                .mean
                .as_ref()
                .zip(stats.sample_variance.as_ref())
                .ok_or("statistics missing")?;
            for ((actual, target), (variance, target_v)) in
                mean.iter().flatten().zip(data.mean.iter().flatten()).zip(
                    variance
                        .iter()
                        .flatten()
                        .zip(data.variance.iter().flatten()),
                )
            {
                loss += ((actual - target) / target.abs().max(1.0)).powi(2)
                    + ((variance - target_v) / target_v.abs().max(1.0)).powi(2);
            }
            if let Some((observable, first)) = amplitude {
                // Fixed-topology temporal shape loss plus a per-path amplitude floor.
                // An ODE screen is never allowed to remove the stochastic evaluation.
                loss += amplitude_penalty(&result.paths, observable, first, replicates)?;
            }
        }
        let value = loss / results.len() as f64;
        if !value.is_finite() {
            return Err("non-finite CRN loss".into());
        }
        Ok(value)
    };
    Ok(Scalar::simulation(objective).domain(lower, upper))
}
fn main() {
    let mut registry = Registry::new();
    for (id, circuit) in [
        ("app.crn.fit-moments", false),
        ("app.crn.circuit-score", true),
    ] {
        let mut about = json!({"simulation":true,"variables":"natural logarithms of parameter_names","inputs":{"model":"crn-craft.model/1.0 (or rules model) exact bytes; convert network models with crn-craft convert","data":"crn-craft.fit-data/0.1; columns in the order of the result's observables, sorted by name"},"engine":"Static Direct (crn-craft.run/1.0); each candidate binds its values to a prepared run","seed_policy":"required explicit seed_set, fixed across every candidate; separate held-out job with its own seed_set required","failure_policy":"incomplete simulation is a failed evaluation","objective":"mean over seed_set of the normalized mean/variance residual sum over all report times and observables"});
        if circuit {
            about["objective"] = json!(
                "mean over seed_set of the normalized mean/variance residual sum plus the mean over all replicate paths of 1/(1+variance) of amplitude_species over report times after burn_in"
            );
            about["parameter_notes"] = json!({"amplitude_species":"required model observable name (a species name for a converted network model)","burn_in":"report times t <= burn_in are excluded from the amplitude penalty; default 0, which excludes the initial state"});
        }
        // Version 2: seed_set is required; circuit amplitude is a per-path mean of a named species
        // after burn_in. Version 3: crn-craft.model/1.0 inputs run on the rules engine; fit-data
        // columns and amplitude_species follow the result's observables.
        registry.scalar(id, "3", about, move |instance| {
            application(instance, circuit)
        });
        let data = json!({"schema":"crn-craft.fit-data/0.1","times":[0.,1.,2.],"mean":[[0.],[2.065016384560419],[2.992888064407119]],"variance":[[0.],[2.065016384560419],[2.992888064407119]]});
        if circuit {
            registry.example(
                id,
                Example::new(2)
                    .bounds(vec![4.0, -1.0], vec![6.0, 1.0])
                    .parameter("parameter_names", vec!["alpha", "gamma_m"])
                    .parameter("replicates", 8)
                    .parameter("seed_set", vec![20261101_u64])
                    .parameter("amplitude_species", "m1")
                    .parameter("burn_in", 0.0)
                    .input_text(
                        "model",
                        "model.json",
                        include_str!("../examples/circuit.model.json"),
                    )
                    .input_text(
                        "data",
                        "data.json",
                        include_str!("../examples/circuit.data.json"),
                    ),
            );
        } else {
            registry.example(
                id,
                Example::new(2)
                    .bounds(vec![-1., -2.], vec![2., 1.])
                    .parameter("replicates", 64)
                    .parameter("seed_set", vec![20261101_u64])
                    .input_text(
                        "model",
                        "model.json",
                        include_str!("../examples/fit.model.json"),
                    )
                    .input_text("data", "data.json", data.to_string()),
            );
        }
    }
    optimization_job_api::worker::run(DiscoveredWorker(Adapter::new(registry)));
}

struct DiscoveredWorker(Adapter);
impl Worker for DiscoveredWorker {
    fn describe(&self) -> serde_json::Value {
        let mut report = self.0.describe();
        report["input_schemas"] = json!({"crn-craft.model/1.0":rules::schema("model-1.0").unwrap(),"crn-craft.fit-data/0.1":schema_for!(Data)});
        report
    }
    fn prepare(
        &self,
        spec: &mut optimization_job_api::JobSpec,
        base: &std::path::Path,
    ) -> Result<Box<dyn optimization_job_api::worker::Prepared>> {
        self.0.prepare(spec, base)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn path(id: u64, values: &[[f64; 2]]) -> PathSeries {
        PathSeries {
            id,
            values: values.iter().map(|c| c.to_vec()).collect(),
        }
    }
    #[test]
    fn amplitude_penalty_is_a_per_path_mean_of_the_named_observable_after_burn_in() {
        // Observable 1 after index 1: [2, 4] has variance 1 (penalty 1/2); [5, 5] has 0 (penalty 1).
        let paths = [
            path(0, &[[9., 0.], [9., 2.], [9., 4.]]),
            path(1, &[[9., 100.], [9., 5.], [9., 5.]]),
        ];
        assert_eq!(amplitude_penalty(&paths, 1, 1, 2), Ok(0.75));
        // Doubling the replicate count with the same paths leaves the mean unchanged.
        let doubled: Vec<_> = paths.iter().chain(&paths).cloned().collect();
        assert_eq!(amplitude_penalty(&doubled, 1, 1, 4), Ok(0.75));
        // Observable 0 is constant: no amplitude on any path.
        assert_eq!(amplitude_penalty(&paths, 0, 1, 2), Ok(1.0));
        // Including the initial state would credit the 100 -> 5 transient as amplitude.
        assert!(amplitude_penalty(&paths, 1, 0, 2).unwrap() < 0.75);
        assert!(amplitude_penalty(&paths, 1, 1, 3).is_err());
    }

    #[test]
    fn the_examples_are_converted_models_with_observable_columns() {
        let fit = include_bytes!("../examples/fit.model.json");
        let circuit = include_bytes!("../examples/circuit.model.json");
        for model in [&fit[..], &circuit[..]] {
            let value: serde_json::Value = serde_json::from_slice(model).unwrap();
            assert_eq!(value["schema"], "crn-craft.model/1.0");
            assert!(rules::check_model(model).valid);
        }
        let data: Data =
            serde_json::from_str(include_str!("../examples/circuit.data.json")).unwrap();
        let run = request(circuit, &data.times, 4, 1, true, 1000, 1);
        let prepared = rules::prepare(circuit, &run).unwrap();
        let cancelled = Control::default();
        cancelled.cancel();
        let observables = rules::execute(&prepared, &cancelled).unwrap().observables;
        assert_eq!(observables, ["m1", "m2", "m3", "p1", "p2", "p3"]);
        assert!(data.mean.iter().all(|r| r.len() == observables.len()));
    }

    #[test]
    fn bound_candidates_run_the_values_and_match_a_fresh_preparation() {
        let model = include_bytes!("../examples/fit.model.json");
        let run = request(model, &[0.0, 1.0, 2.0], 64, 20261101, false, 1_000_000, 2);
        let template = rules::prepare(model, &run).unwrap();
        let values = BTreeMap::from([("birth".to_string(), 5.0), ("death".to_string(), 0.5)]);
        let bound = rules::execute(
            &template.bind_parameters(&values).unwrap(),
            &Control::default(),
        )
        .unwrap();
        let mut edited: serde_json::Value = serde_json::from_slice(model).unwrap();
        edited["parameters"] = json!({"birth": 5.0, "death": 0.5});
        let edited = serde_json::to_vec(&edited).unwrap();
        let fresh = rules::execute(
            &rules::prepare(
                &edited,
                &request(&edited, &[0.0, 1.0, 2.0], 64, 20261101, false, 1_000_000, 1),
            )
            .unwrap(),
            &Control::default(),
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&bound.statistics).unwrap(),
            serde_json::to_value(&fresh.statistics).unwrap()
        );
        // The stationary mean is birth / death = 10; at t = 2 it is 10 (1 - e^-1).
        let mean = bound.statistics.unwrap().mean.unwrap()[2][0];
        assert!(
            (mean - 10.0 * (1.0 - (-1.0f64).exp())).abs() < 1.0,
            "{mean}"
        );
    }
}
