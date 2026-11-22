//! Concentration-only BDF, with a persistent batch evaluator for fitting.
use crn_craft::{concentration::Concentration, strict_json};
use diffsol::{NalgebraLU, NalgebraMat, OdeBuilder, OdeSolverMethod, OdeSolverStopReason};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, Write};

#[derive(Deserialize, JsonSchema)]
enum ModelVersion {
    #[serde(rename = "crn-craft.concentration-model/0.1")]
    V01,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Model {
    schema: ModelVersion,
    states: Vec<String>,
    parameters: Vec<String>,
    initial: Vec<String>,
    #[serde(default)]
    assignments: Vec<(String, String)>,
    derivatives: Vec<String>,
    #[serde(default)]
    observations: Vec<String>,
}

#[derive(Deserialize, JsonSchema)]
enum RequestVersion {
    #[serde(rename = "crn-craft.concentration-request/0.1")]
    V01,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Batch {
    schema: RequestVersion,
    cases: Vec<Case>,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Case {
    parameters: Vec<f64>,
    times: Vec<f64>,
    rtol: f64,
    atol: f64,
}

#[derive(Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum Evaluation {
    Values(Trajectory),
    Error(String),
}

#[derive(Serialize, JsonSchema)]
struct Trajectory {
    states: Vec<Vec<f64>>,
    observations: Vec<Vec<f64>>,
}

#[derive(Serialize, JsonSchema)]
enum ResultVersion {
    #[serde(rename = "crn-craft.concentration-result/0.1")]
    V01,
}

#[derive(Serialize, JsonSchema)]
struct BatchResult {
    schema: ResultVersion,
    cases: Vec<Evaluation>,
}

fn solve(model: &Concentration, case: &Case) -> Result<Trajectory, String> {
    if !(2..=201).contains(&case.times.len())
        || case.times[0] != 0.0
        || case
            .times
            .windows(2)
            .any(|w| !w[1].is_finite() || w[1] <= w[0])
        || !(1e-13..=1e-2).contains(&case.rtol)
        || !(1e-15..=1e-2).contains(&case.atol)
    {
        return Err("invalid concentration times or tolerances".into());
    }
    let initial = model.initial(&case.parameters)?;
    let n = initial.len();
    let problem = OdeBuilder::<NalgebraMat<f64>>::new()
        .p(case.parameters.clone())
        .rtol(case.rtol)
        .atol([case.atol])
        .rhs_implicit(
            |y, p, t, out| model.rhs(t, y, p, out),
            |y, p, t, v, out| model.jacobian_product(t, y, p, v, out),
        )
        .init(move |_, _, y| y.copy_from_slice(&initial), n)
        .build()
        .map_err(|e| e.to_string())?;
    let mut solver = problem
        .bdf::<NalgebraLU<f64>>()
        .map_err(|e| e.to_string())?;
    solver
        .set_stop_time(*case.times.last().expect("validated times"))
        .map_err(|e| e.to_string())?;
    let mut states = vec![(0..n).map(|i| solver.state().y[i]).collect::<Vec<_>>()];
    let mut reached = false;
    for _ in 0..20_000 {
        let status = solver.step().map_err(|e| e.to_string())?;
        while states.len() < case.times.len() && case.times[states.len()] <= solver.state().t {
            let y = solver
                .interpolate(case.times[states.len()])
                .map_err(|e| e.to_string())?;
            states.push((0..n).map(|i| y[i]).collect());
        }
        if matches!(status, OdeSolverStopReason::TstopReached) {
            reached = true;
            break;
        }
    }
    if !reached || states.len() != case.times.len() {
        return Err("incomplete concentration solve within the step limit".into());
    }
    let scales: Vec<_> = (0..n)
        .map(|i| states.iter().map(|row| row[i]).fold(0.0, f64::max))
        .collect();
    let mut observations = Vec::with_capacity(case.times.len());
    for (values, &time) in states.iter().zip(&case.times) {
        if values
            .iter()
            .zip(&scales)
            .any(|(v, scale)| !v.is_finite() || *v < -(case.atol + case.rtol * scale))
        {
            return Err(format!("invalid reported concentration at {time}"));
        }
        observations.push(model.observe(time, values, &case.parameters)?);
    }
    Ok(Trajectory {
        states,
        observations,
    })
}

pub(crate) fn command(args: &[String]) -> Result<(), String> {
    if args.len() == 1 && args[0] == "describe" {
        println!(
            "{}",
            serde_json::json!({
                "model_schema": schema_for!(Model),
                "request_schema": schema_for!(Batch),
                "result_schema": schema_for!(BatchResult)
            })
        );
        return Ok(());
    }
    if args.len() != 1 {
        return Err("usage: crn-diffsol concentrations MODEL (read JSON lines on stdin)".into());
    }
    let bytes = std::fs::read(&args[0]).map_err(|e| e.to_string())?;
    if bytes.len() > 1 << 20 {
        return Err("concentration model is too large".into());
    }
    let model: Model = strict_json::parse(&bytes).map_err(|e| e.message)?;
    let _ = model.schema;
    let prepared = Concentration::new(
        &model.states,
        &model.parameters,
        &model.initial,
        &model.assignments,
        &model.derivatives,
        &model.observations,
    )?;
    let mut out = std::io::stdout().lock();
    for line in std::io::stdin().lock().lines() {
        let line = line.map_err(|e| e.to_string())?;
        let reply = match strict_json::parse::<Batch>(line.as_bytes()) {
            Ok(batch) if line.len() <= 1 << 20 && (1..=32).contains(&batch.cases.len()) => {
                let _ = batch.schema;
                let cases = batch
                    .cases
                    .iter()
                    .map(|case| match solve(&prepared, case) {
                        Ok(value) => Evaluation::Values(value),
                        Err(error) => Evaluation::Error(error),
                    })
                    .collect();
                serde_json::to_value(BatchResult {
                    schema: ResultVersion::V01,
                    cases,
                })
                .map_err(|e| e.to_string())?
            }
            Ok(_) => serde_json::json!({"error": "invalid concentration batch size"}),
            Err(error) => serde_json::json!({"error": error.message}),
        };
        writeln!(out, "{reply}").map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}
