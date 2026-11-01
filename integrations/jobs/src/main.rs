#![forbid(unsafe_code)]
//! Register fallible CRN applications in the existing fcmaes worker, no new protocol.
use crn_craft::{
    contracts::Model,
    simulation::{self, Control, Stop},
};
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
#[allow(clippy::too_many_arguments)]
fn simulation(
    model: &Model,
    times: &[f64],
    replicates: u64,
    seeds: &[u64],
    x: &[f64],
    names: &[String],
    ctx: &SimContext,
    keep: bool,
    max_events: u64,
) -> std::result::Result<Vec<simulation::SimulationResult>, String> {
    let mut model = model.clone();
    for (name, log_rate) in names.iter().zip(x) {
        let rate = log_rate.exp();
        if !rate.is_finite() || rate <= 0.0 {
            return Err("invalid log-rate".into());
        }
        model.parameters.insert(name.clone(), rate);
    }
    let bytes = serde_json::to_vec(&model).map_err(|e| e.to_string())?;
    let mut results = Vec::new();
    for &seed in seeds {
        ctx.checkpoint()?;
        let request = json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(&bytes)},"method":{"kind":"ssa","aggregator":"direct"},"ensemble":{"replicates":replicates,"root_seed":seed},"time":{"end":times.last(),"save_times":times},"outputs":{"statistics":["mean","variance"],"keep_trajectories":if keep {replicates} else {0}},"budget":{"threads":ctx.threads().saturating_sub(1).clamp(1,64),"max_events_per_replicate":max_events,"max_output_bytes":268435456,"wall_time_s":ctx.remaining_seconds().max(1e-12)}});
        let run = serde_json::to_vec(&request).unwrap();
        let prepared =
            simulation::prepare(&bytes, &run).map_err(|e| serde_json::to_string(&e).unwrap())?;
        let context = ctx.clone();
        let mut control = Control::default();
        control.external_stop = Some(Box::new(move || context.is_cancelled()));
        let result = simulation::execute_cached(&prepared, &control)
            .map_err(|e| e.message)?
            .result;
        if result.termination != Stop::Completed {
            return Err(format!(
                "incomplete CRN objective: {:?}",
                result.termination
            ));
        }
        results.push(result);
    }
    Ok(results)
}
fn application(instance: &Instance, circuit: bool) -> Result<Scalar> {
    let model: Model = crn_craft::parse_model(&instance.input("model")?).map_err(|e| e.message)?;
    crn_craft::CompiledNetwork::compile(&model).map_err(|e| serde_json::to_string(&e).unwrap())?;
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
    // Fixed seeds across candidates. SimContext.seed is point-derived and deliberately unused.
    let seeds: Vec<u64> = instance.parameter("seed_set", vec![20261101_u64])?;
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
        || data.mean.len() != data.times.len()
        || data.variance.len() != data.times.len()
        || data
            .mean
            .iter()
            .chain(&data.variance)
            .any(|r| r.len() != model.species.len() || r.iter().any(|v| !v.is_finite() || *v < 0.0))
    {
        return Err("invalid fit-data dimensions/domain/grid".into());
    }
    let lower = vec![-12.0; names.len()];
    let upper = vec![12.0; names.len()];
    let objective = move |x: &[f64], ctx: &SimContext| -> std::result::Result<f64, String> {
        let results = simulation(
            &model,
            &data.times,
            replicates,
            &seeds,
            x,
            &names,
            ctx,
            circuit,
            max_events,
        )?;
        let mut loss = 0.0;
        for result in &results {
            let stats = result.statistics.as_ref().ok_or("statistics missing")?;
            for ((actual, target), (variance, target_v)) in stats
                .mean
                .as_ref()
                .unwrap()
                .iter()
                .flatten()
                .zip(data.mean.iter().flatten())
                .zip(
                    stats
                        .sample_variance
                        .as_ref()
                        .unwrap()
                        .iter()
                        .flatten()
                        .zip(data.variance.iter().flatten()),
                )
            {
                loss += ((actual - target) / target.abs().max(1.0)).powi(2)
                    + ((variance - target_v) / target_v.abs().max(1.0)).powi(2);
            }
            if circuit {
                // Fixed-topology temporal shape loss, per-trajectory amplitude floor.
                // An ODE screen is never allowed to remove the stochastic evaluation.
                for path in &result.trajectories {
                    let x: Vec<_> = path.counts.iter().map(|s| s[0] as f64).collect();
                    let mean = x.iter().sum::<f64>() / x.len() as f64;
                    let variance =
                        x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64;
                    loss += 1.0 / (1.0 + variance);
                }
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
        registry.scalar(id,"1",json!({"simulation":true,"variables":"natural logarithms of parameter_names","inputs":{"model":"crn-craft.model/0.1 exact bytes","data":"crn-craft.fit-data/0.1"},"seed_policy":"explicit fixed seed_set across every candidate; separate held-out job required","failure_policy":"incomplete simulation is a failed evaluation"}),move |instance|application(instance,circuit));
        let model = json!({"schema":"crn-craft.model/0.1","name":"ImmigrationDeathFit","convention":"stochastic_counts","omega":1.0,"parameters":{"birth":3.0,"death":0.8},"species":[{"name":"X","initial_count":0}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":"birth"}},{"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":"death"}}]});
        let data = json!({"schema":"crn-craft.fit-data/0.1","times":[0.,1.,2.],"mean":[[0.],[2.065016384560419],[2.992888064407119]],"variance":[[0.],[2.065016384560419],[2.992888064407119]]});
        if circuit {
            registry.example(
                id,
                Example::new(2)
                    .bounds(vec![4.0, -1.0], vec![6.0, 1.0])
                    .parameter("parameter_names", vec!["alpha", "gamma_m"])
                    .parameter("replicates", 8)
                    .parameter("seed_set", vec![20261101_u64])
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
                    .input_text("model", "model.json", model.to_string())
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
        report["input_schemas"] = json!({"crn-craft.model/0.1":crn_craft::schema("model").unwrap(),"crn-craft.fit-data/0.1":schema_for!(Data)});
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
