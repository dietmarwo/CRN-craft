#![forbid(unsafe_code)]
#![recursion_limit = "512"]
//! Isolated benchmark worker. The numerical library has no optimizer/rebop dependency.
use crn_craft::{
    contracts::Model,
    parameters::{BoundParameters, ParameterTemplate},
    simulation::Control,
    ssa,
};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    io::{self, BufRead, Write},
    time::Instant,
};

#[cfg(feature = "generated")]
mod crn_binding12;
#[cfg(feature = "generated")]
mod crn_birth_death;
#[cfg(feature = "generated")]
mod crn_guard;
#[cfg(feature = "generated")]
mod crn_sparse33;
#[cfg(feature = "rebop-macro")]
mod rebop_binding12;
#[cfg(feature = "rebop-macro")]
mod rebop_birth_death;
#[cfg(feature = "rebop-macro")]
mod rebop_sparse33;

#[cfg(feature = "rebop-macro")]
struct Moments {
    n: u64,
    mean: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
}
#[cfg(feature = "rebop-macro")]
impl Moments {
    fn new(times: usize, species: usize) -> Self {
        Self {
            n: 0,
            mean: vec![vec![0.; species]; times],
            m2: vec![vec![0.; species]; times],
        }
    }
    fn add(&mut self, states: &[Vec<f64>]) {
        self.n += 1;
        for ((mean, m2), state) in self.mean.iter_mut().zip(&mut self.m2).zip(states) {
            for ((mean, m2), x) in mean.iter_mut().zip(m2).zip(state) {
                let delta = x - *mean;
                *mean += delta / self.n as f64;
                *m2 += delta * (x - *mean);
            }
        }
    }
    fn finish(self) -> Value {
        let variance: Vec<Vec<_>> = self
            .m2
            .into_iter()
            .map(|r| r.into_iter().map(|v| v / (self.n - 1) as f64).collect())
            .collect();
        json!({"mean":self.mean,"variance":variance})
    }
}
struct State {
    engine: String,
    case: String,
    model: Model,
    names: Vec<String>,
    times: Vec<f64>,
    requests: BTreeMap<String, ssa::Request>,
    templates: BTreeMap<String, ParameterTemplate>,
}
fn error(value: impl std::fmt::Debug) -> String {
    format!("{value:?}")
}
fn diagnostics(value: Vec<crn_craft::diagnostic::Diagnostic>) -> String {
    serde_json::to_string(&value).unwrap()
}

impl State {
    fn new(v: &Value) -> Result<Self, String> {
        let model_bytes = serde_json::to_vec(&v["model"]).map_err(error)?;
        let model = crn_craft::parse_model(&model_bytes).map_err(error)?;
        let times: Vec<f64> = serde_json::from_value(v["times"].clone()).map_err(error)?;
        let engine = v["engine"].as_str().ok_or("engine")?.to_owned();
        let case = v["case"].as_str().ok_or("case")?.to_owned();
        let mut requests = BTreeMap::new();
        let mut templates = BTreeMap::new();
        for (phase, paths, seed) in [
            (
                "train",
                v["paths"].as_u64().unwrap_or(32),
                v["seed"].as_u64().unwrap_or(77201),
            ),
            (
                "heldout",
                v["heldout_paths"].as_u64().unwrap_or(4096),
                880431,
            ),
        ] {
            let request = serde_json::to_vec(&json!({"schema":if v["fast"].as_bool().unwrap_or(true) {"crn-craft.ssa-run/0.2"} else {"crn-craft.ssa-run/0.1"},
                "model":{"path":"model.json","sha256":crn_craft::sha256(&model_bytes)},
                "method":{"kind":"ssa","aggregator":v["algorithm"].as_str().unwrap_or("cached_direct"),"waiting_time":v["wait"].as_str().unwrap_or("ziggurat"),"fast":v["fast"].as_bool().unwrap_or(true)},
                "time":{"end":times.last(),"save_times":times},"ensemble":{"replicates":paths,"root_seed":seed},
                "outputs":{"statistics":["mean","variance"],"keep_trajectories":v["keep"].as_u64().unwrap_or(0)},
                "budget":{"threads":1,"max_events_per_replicate":1000000,"max_sampling_steps_per_replicate":1000000000,"max_output_bytes":64000000}})).map_err(error)?;
            if engine != "rebop_macro" {
                templates.insert(
                    phase.into(),
                    ParameterTemplate::prepare(&model_bytes, &request).map_err(diagnostics)?,
                );
            }
            requests.insert(
                phase.into(),
                crn_craft::strict_json::parse(&request).map_err(error)?,
            );
        }
        let names = model.parameters.keys().cloned().collect();
        Ok(Self {
            engine,
            case,
            model,
            names,
            times,
            requests,
            templates,
        })
    }
    fn evaluate(&self, v: &Value) -> Result<Value, String> {
        let values: Vec<f64> = serde_json::from_value(v["values"].clone()).map_err(error)?;
        if values.len() != self.names.len() || values.iter().any(|x| !x.is_finite() || *x < 0.) {
            return Err("invalid benchmark values".into());
        }
        let phase = v["phase"].as_str().unwrap_or("train");
        let request = self.requests.get(phase).ok_or("phase")?;
        if self.engine == "rebop_macro" {
            #[cfg(feature = "rebop-macro")]
            {
                let seeds: Vec<u64> =
                    serde_json::from_value(v["path_seeds"].clone()).map_err(error)?;
                let paths = request.ensemble.replicates;
                if seeds.len() != paths as usize {
                    return Err("seed count".into());
                }
                let start = Instant::now();
                let moments = match self.case.as_str() {
                    "birth_death" => rebop_birth_death::run(&values, &self.times, paths, &seeds)?,
                    "sparse33" => rebop_sparse33::run(&values, &self.times, paths, &seeds)?,
                    "binding12" => rebop_binding12::run(&values, &self.times, paths, &seeds)?,
                    _ => return Err("case".into()),
                };
                return Ok(
                    json!({"moments":moments,"simulation_seconds":start.elapsed().as_secs_f64(),"bind_seconds":0.}),
                );
            }
            #[cfg(not(feature = "rebop-macro"))]
            return Err("rebop-macro feature not built".into());
        }
        let start = Instant::now();
        let bound = if self.engine == "fresh" {
            let mut model = self.model.clone();
            for (name, value) in self.names.iter().zip(&values) {
                model.parameters.insert(name.clone(), *value);
            }
            let model_bytes = serde_json::to_vec(&model).map_err(error)?;
            let mut request = request.clone();
            request.model.sha256 = crn_craft::sha256(&model_bytes);
            let request_bytes = serde_json::to_vec(&request).map_err(error)?;
            let prepared = ssa::prepare(&model_bytes, &request_bytes).map_err(diagnostics)?;
            BoundParameters {
                prepared,
                model_bytes,
                request_bytes,
            }
        } else {
            self.templates
                .get(phase)
                .ok_or("template")?
                .bind(&values, None)
                .map_err(diagnostics)?
        };
        let binding_seconds = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let control = Control::default();
        let result = if self.engine == "generated" {
            #[cfg(feature = "generated")]
            {
                match self.case.as_str() {
                    "birth_death" => ssa::execute_with_kernel(
                        &bound.prepared,
                        &control,
                        &crn_birth_death::GeneratedKernel,
                    ),
                    "sparse33" => ssa::execute_with_kernel(
                        &bound.prepared,
                        &control,
                        &crn_sparse33::GeneratedKernel,
                    ),
                    "binding12" => ssa::execute_with_kernel(
                        &bound.prepared,
                        &control,
                        &crn_binding12::GeneratedKernel,
                    ),
                    _ => return Err("case".into()),
                }
            }
            #[cfg(not(feature = "generated"))]
            {
                return Err("generated feature not built".into());
            }
        } else if self.engine == "fresh" || self.engine == "reuse" {
            ssa::execute(&bound.prepared, &control)
        } else {
            return Err("engine".into());
        }
        .map_err(error)?;
        let simulation_seconds = start.elapsed().as_secs_f64();
        if result.termination != ssa::Termination::Completed {
            return Err(format!("incomplete: {:?}", result.termination));
        }
        let stats = result.statistics.as_ref().ok_or("missing statistics")?;
        let mut reply = json!({"moments":{"mean":stats.mean,"variance":stats.sample_variance},"events":result.events,
            "bind_seconds":binding_seconds,"simulation_seconds":simulation_seconds,
            "result_sha256":crn_craft::sha256(&serde_json::to_vec(&result).map_err(error)?)});
        if v["full"].as_bool().unwrap_or(false) {
            reply["result"] = serde_json::to_value(&result).map_err(error)?;
        }
        Ok(reply)
    }
}

#[cfg(feature = "generated")]
fn check_kernel<K: ssa::PropensityKernel>(model: &Model, kernel: &K) -> Result<usize, String> {
    let mut checks = 0;
    for coefficient in [0., -0., f64::from_bits(1), f64::MIN_POSITIVE, 1., f64::MAX] {
        let mut model = model.clone();
        for value in model.parameters.values_mut() {
            *value = coefficient;
        }
        let net = crn_craft::CompiledNetwork::compile(&model).map_err(diagnostics)?;
        for offset in [0, 1, 2, 3, 5, 31, i64::MAX - 64] {
            let state: Vec<_> = (0..net.species().len())
                .map(|i| offset + i as i64)
                .collect();
            for rx in 0..net.reactions().len() {
                let a = net.propensity(rx, &state);
                let b = kernel.evaluate(&net, rx, &state);
                match (a, b) {
                    (Ok(a), Ok(b)) if a.to_bits() == b.to_bits() => (),
                    (Err(a), Err(b))
                        if serde_json::to_value(&a).unwrap()
                            == serde_json::to_value(&b).unwrap() => {}
                    mismatch => return Err(format!("kernel evaluation mismatch: {mismatch:?}")),
                }
                checks += 1;
            }
        }
    }
    Ok(checks)
}
// Internal benchmark requests are strict too: malformed optional fields never
// silently select a default or hide a misspelled setting.
fn validate_fields(value: &Value) -> Result<(), String> {
    let object = value.as_object().ok_or("request must be an object")?;
    let command = value["command"].as_str().ok_or("command")?;
    let allowed: &[&str] = match command {
        "init" => &[
            "command",
            "engine",
            "case",
            "model",
            "times",
            "paths",
            "heldout_paths",
            "seed",
            "keep",
            "algorithm",
            "wait",
            "fast",
        ],
        "eval" => &["command", "values", "phase", "path_seeds", "full"],
        "validate_kernel" => &["command"],
        _ => return Err("unknown command".into()),
    };
    for key in object.keys() {
        if !allowed.contains(&key.as_str()) {
            return Err(format!("unknown field {key:?}"));
        }
    }
    for key in ["paths", "heldout_paths", "seed", "keep"] {
        if object.contains_key(key) && value[key].as_u64().is_none() {
            return Err(format!("{key} must be an unsigned integer"));
        }
    }
    for key in ["fast", "full"] {
        if object.contains_key(key) && value[key].as_bool().is_none() {
            return Err(format!("{key} must be boolean"));
        }
    }
    for key in ["algorithm", "wait", "phase"] {
        if object.contains_key(key) && value[key].as_str().is_none() {
            return Err(format!("{key} must be a string"));
        }
    }
    Ok(())
}

fn main() {
    let mut state: Option<State> = None;
    let mut stdout = io::stdout().lock();
    for line in io::stdin().lock().lines() {
        let answer = (|| -> Result<Value, String> {
            let line = line.map_err(error)?;
            let value: Value = crn_craft::strict_json::parse(line.as_bytes()).map_err(error)?;
            validate_fields(&value)?;
            match value["command"].as_str().ok_or("command")? {
                "init" => {
                    let start = Instant::now();
                    state = Some(State::new(&value)?);
                    Ok(
                        json!({"ready":true,"case":state.as_ref().unwrap().case,"save_points":state.as_ref().unwrap().times.len(),"initialization_seconds":start.elapsed().as_secs_f64()}),
                    )
                }
                "eval" => state.as_ref().ok_or("not initialized")?.evaluate(&value),
                "validate_kernel" => {
                    #[cfg(feature = "generated")]
                    {
                        let state = state.as_ref().ok_or("not initialized")?;
                        let checks = match state.case.as_str() {
                            "birth_death" => {
                                check_kernel(&state.model, &crn_birth_death::GeneratedKernel)?
                            }
                            "sparse33" => {
                                check_kernel(&state.model, &crn_sparse33::GeneratedKernel)?
                            }
                            "binding12" => {
                                check_kernel(&state.model, &crn_binding12::GeneratedKernel)?
                            }
                            "guard" => check_kernel(&state.model, &crn_guard::GeneratedKernel)?,
                            _ => return Err("case".into()),
                        };
                        Ok(json!({"checks":checks,"passed":true}))
                    }
                    #[cfg(not(feature = "generated"))]
                    {
                        Err("generated feature not built".into())
                    }
                }
                _ => Err("unknown command".into()),
            }
        })();
        let reply = match answer {
            Ok(value) => json!({"ok":true,"result":value}),
            Err(message) => json!({"ok":false,"error":message}),
        };
        if writeln!(stdout, "{reply}").is_err() || stdout.flush().is_err() {
            break;
        }
    }
}
