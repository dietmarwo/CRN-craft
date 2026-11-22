//! Concentration-only BDF, with a persistent batch evaluator for fitting.
use crn_craft::{concentration::Concentration, strict_json};
use diffsol::{NalgebraLU, NalgebraMat, OdeBuilder, OdeSolverMethod, OdeSolverStopReason};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, Read, Write};
use std::panic::{AssertUnwindSafe, catch_unwind};

/// Largest request line in bytes, newline excluded.
const LINE_LIMIT: usize = 1 << 20;

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

/// A refused case or line: a stable diagnostic code and a free message.
/// `CRN-E001`-`CRN-E003`: the line is not a valid request document;
/// `CRN-E004`: line or batch size; `CRN-E015`: times, tolerances or
/// parameters out of their domain; `CRN-E080`: the integration failed, a
/// reported state went negative beyond the tolerance, or a readout is
/// not finite.
#[derive(Serialize, JsonSchema)]
struct Failure {
    code: String,
    message: String,
}

impl Failure {
    fn new(code: &str, message: impl Into<String>) -> Self {
        Self {
            code: code.into(),
            message: message.into(),
        }
    }
}

#[derive(Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum Evaluation {
    Values(Trajectory),
    Error(Failure),
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

/// One result line: `cases` in request order, or `error` for a line that
/// was refused as a whole.
#[derive(Serialize, JsonSchema)]
struct BatchResult {
    schema: ResultVersion,
    #[serde(skip_serializing_if = "Option::is_none")]
    cases: Option<Vec<Evaluation>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<Failure>,
}

fn solve(model: &Concentration, case: &Case) -> Result<Trajectory, Failure> {
    if !(2..=201).contains(&case.times.len())
        || case.times[0] != 0.0
        || case
            .times
            .windows(2)
            .any(|w| !w[1].is_finite() || w[1] <= w[0])
        || !(1e-13..=1e-2).contains(&case.rtol)
        || !(1e-15..=1e-2).contains(&case.atol)
    {
        return Err(Failure::new(
            "CRN-E015",
            "invalid concentration times or tolerances",
        ));
    }
    let initial = model
        .initial(&case.parameters)
        .map_err(|e| Failure::new("CRN-E015", e))?;
    integrate(model, case, initial).map_err(|e| Failure::new("CRN-E080", e))
}

fn integrate(model: &Concentration, case: &Case, initial: Vec<f64>) -> Result<Trajectory, String> {
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

fn reply(model: &Concentration, line: &[u8]) -> BatchResult {
    let refuse = |error| BatchResult {
        schema: ResultVersion::V01,
        cases: None,
        error: Some(error),
    };
    if line.len() > LINE_LIMIT {
        return refuse(Failure::new(
            "CRN-E004",
            format!("request line exceeds {LINE_LIMIT} bytes"),
        ));
    }
    let batch = match strict_json::parse::<Batch>(line) {
        Ok(batch) => batch,
        Err(error) => return refuse(Failure::new(&error.code, error.message)),
    };
    let _ = batch.schema;
    if !(1..=32).contains(&batch.cases.len()) {
        return refuse(Failure::new("CRN-E004", "a request holds 1 to 32 cases"));
    }
    let cases = batch
        .cases
        .iter()
        .map(|case| {
            // A panic inside the integrator fails this case only.
            match catch_unwind(AssertUnwindSafe(|| solve(model, case))) {
                Ok(Ok(value)) => Evaluation::Values(value),
                Ok(Err(error)) => Evaluation::Error(error),
                Err(_) => Evaluation::Error(Failure::new("CRN-E080", "internal ODE solver panic")),
            }
        })
        .collect();
    BatchResult {
        schema: ResultVersion::V01,
        cases: Some(cases),
        error: None,
    }
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
    let mut input = std::io::stdin().lock();
    let mut line = Vec::new();
    loop {
        line.clear();
        let read = input
            .by_ref()
            .take(LINE_LIMIT as u64 + 1)
            .read_until(b'\n', &mut line)
            .map_err(|e| e.to_string())?;
        if read == 0 {
            break;
        }
        let complete = line.last() == Some(&b'\n');
        if complete {
            line.pop();
        } else if line.len() > LINE_LIMIT {
            // Discard the rest of an oversized line.
            let mut rest = Vec::new();
            input
                .read_until(b'\n', &mut rest)
                .map_err(|e| e.to_string())?;
        }
        let reply = reply(&prepared, &line);
        let text = serde_json::to_string(&reply).map_err(|e| e.to_string())?;
        writeln!(out, "{text}").map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}
