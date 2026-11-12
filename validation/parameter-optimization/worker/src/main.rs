#![forbid(unsafe_code)]
#![recursion_limit = "512"]
//! Isolated benchmark worker. The numerical library has no optimizer/rebop dependency.
//!
//! The CRN routes run the rules engine's Static Direct on the model converted
//! to `crn-craft.model/1.0`: `fresh` prepares every candidate from its model
//! bytes, `reuse` binds the candidate's values to a prepared run
//! (`Prepared::bind_parameters`). Generated rules kernels and loops compile
//! the parameter values in and refuse rebound candidates, so no generated
//! route serves a parameter search; both routes are interpreted.
use crn_craft::rules::{self, Control, Prepared, RulesResult, RunRequest, Termination};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    io::{self, BufRead, Write},
    time::Instant,
};

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
    /// The rules model (a network model converted to `crn-craft.model/1.0`).
    model: rules::Model,
    names: Vec<String>,
    times: Vec<f64>,
    requests: BTreeMap<String, RunRequest>,
    templates: BTreeMap<String, Prepared>,
    /// Result column of each type, in type declaration order: the network
    /// model's species order for a converted model.
    columns: Option<Vec<usize>>,
}
fn error(value: impl std::fmt::Debug) -> String {
    format!("{value:?}")
}
fn diagnostics(value: Vec<crn_craft::diagnostic::Diagnostic>) -> String {
    serde_json::to_string(&value).unwrap()
}

impl State {
    fn new(v: &Value) -> Result<Self, String> {
        let supplied = serde_json::to_vec(&v["model"]).map_err(error)?;
        let model_bytes = if v["model"]["schema"] == "crn-craft.model/0.1" {
            rules::convert::from_network(&supplied).map_err(diagnostics)?
        } else {
            supplied
        };
        let model: rules::Model = crn_craft::strict_json::parse(&model_bytes).map_err(error)?;
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
            let request = serde_json::to_vec(&json!({"schema":"crn-craft.run/1.0",
                "model":{"path":"model.json","sha256":crn_craft::sha256(&model_bytes)},
                "method":{"kind":"static_direct","waiting_time":v["wait"].as_str().unwrap_or("ziggurat")},
                "time":{"end":times.last(),"report_times":times},"ensemble":{"replicates":paths,"root_seed":seed},
                "outputs":{"statistics":["mean","variance"],"keep_paths":v["keep"].as_u64().unwrap_or(0),"keep_instances":0},
                "budget":{"threads":1,"max_events_per_replicate":1000000,"max_sampling_steps_per_replicate":1000000000,
                          "max_instances":1000000,"max_output_bytes":64000000}})).map_err(error)?;
            if engine != "rebop_macro" {
                templates.insert(
                    phase.into(),
                    rules::prepare(&model_bytes, &request).map_err(diagnostics)?,
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
            columns: None,
        })
    }
    /// The result's columns in type declaration order when every type has an
    /// observable of its name (a converted network model), else result order.
    fn columns(&self, result: &RulesResult) -> Vec<usize> {
        let by_type: Option<Vec<usize>> = self
            .model
            .types
            .iter()
            .map(|t| result.observables.iter().position(|o| *o == t.name))
            .collect();
        match by_type {
            Some(columns) if columns.len() == result.observables.len() => columns,
            _ => (0..result.observables.len()).collect(),
        }
    }
    fn evaluate(&mut self, v: &Value) -> Result<Value, String> {
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
            rules::prepare(&model_bytes, &request_bytes).map_err(diagnostics)?
        } else if self.engine == "reuse" {
            let overrides: BTreeMap<String, f64> = self
                .names
                .iter()
                .cloned()
                .zip(values.iter().copied())
                .collect();
            self.templates
                .get(phase)
                .ok_or("template")?
                .bind_parameters(&overrides)
                .map_err(diagnostics)?
        } else {
            return Err("engine".into());
        };
        let binding_seconds = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let result = rules::execute(&bound, &Control::default()).map_err(error)?;
        let simulation_seconds = start.elapsed().as_secs_f64();
        if result.termination != Termination::Completed {
            return Err(format!("incomplete: {:?}", result.termination));
        }
        let columns = match &self.columns {
            Some(columns) => columns.clone(),
            None => {
                let columns = self.columns(&result);
                self.columns = Some(columns.clone());
                columns
            }
        };
        let stats = result.statistics.as_ref().ok_or("missing statistics")?;
        let order = |rows: &Option<Vec<Vec<f64>>>| -> Result<Vec<Vec<f64>>, String> {
            Ok(rows
                .as_ref()
                .ok_or("missing statistics")?
                .iter()
                .map(|row| columns.iter().map(|&c| row[c]).collect())
                .collect())
        };
        let mut reply = json!({"moments":{"mean":order(&stats.mean)?,"variance":order(&stats.sample_variance)?},
            "events":result.events,"bind_seconds":binding_seconds,"simulation_seconds":simulation_seconds,
            "result_sha256":crn_craft::sha256(&serde_json::to_vec(&result).map_err(error)?)});
        if v["full"].as_bool().unwrap_or(false) {
            reply["result"] = serde_json::to_value(&result).map_err(error)?;
        }
        Ok(reply)
    }
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
            "wait",
        ],
        "eval" => &["command", "values", "phase", "path_seeds", "full"],
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
    if object.contains_key("full") && value["full"].as_bool().is_none() {
        return Err("full must be boolean".into());
    }
    for key in ["wait", "phase"] {
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
                "eval" => state.as_mut().ok_or("not initialized")?.evaluate(&value),
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
