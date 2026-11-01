#![forbid(unsafe_code)]
use crn_craft::kinetics::Kinetics;
use diffsol::{OdeBuilder, OdeSolverStopReason, SensitivitiesOdeSolverMethod};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
enum RequestVersion {
    #[serde(rename = "crn-craft.ode-request/0.1")]
    V01,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Request {
    schema: RequestVersion,
    times: Vec<f64>,
    rtol: f64,
    atol: f64,
}
#[derive(Serialize, JsonSchema)]
enum OutputVersion {
    #[serde(rename = "crn-craft.ode-result/0.1")]
    V01,
}
#[derive(Serialize, JsonSchema)]
struct Output {
    schema: OutputVersion,
    solver: &'static str,
    convention: &'static str,
    model_sha256: String,
    times: Vec<f64>,
    parameters: Vec<String>,
    values: Vec<Vec<f64>>,
    sensitivities: Vec<Vec<Vec<f64>>>,
}
fn main() {
    if let Err(e) = run() {
        eprintln!(
            "{}",
            serde_json::json!({"success":false,"code":"CRN-E080","message":e.to_string()})
        );
        std::process::exit(1);
    }
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() == 2 && args[1] == "describe" {
        println!(
            "{}",
            serde_json::json!({"schema":"crn-craft.ode-discovery/0.1","request_schema":schema_for!(Request),"result_schema":schema_for!(Output)})
        );
        return Ok(());
    }
    if args.len() != 3 {
        return Err("usage: crn-diffsol MODEL REQUEST".into());
    }
    let bytes = std::fs::read(&args[1])?;
    let model = crn_craft::parse_model(&bytes).map_err(|e| e.message)?;
    let request: Request =
        crn_craft::strict_json::parse(&std::fs::read(&args[2])?).map_err(|e| e.message)?;
    let _schema = &request.schema;
    if request.times.len() < 2
        || request.times.len() > 10000
        || request.times[0] != 0.
        || request
            .times
            .windows(2)
            .any(|w| !w[1].is_finite() || w[1] <= w[0])
        || !request.rtol.is_finite()
        || request.rtol <= 0.
        || request.rtol > 0.1
        || !request.atol.is_finite()
        || request.atol <= 0.
        || request.atol > 0.1
    {
        return Err("invalid ODE request".into());
    }
    let kinetics =
        Arc::new(Kinetics::compile(&model).map_err(|e| serde_json::to_string(&e).unwrap())?);
    let parameters = kinetics.parameters.clone();
    let p: Vec<_> = parameters.iter().map(|s| model.parameters[s]).collect();
    if p.is_empty() {
        return Err("sensitivity path requires at least one parameter".into());
    }
    let n = model.species.len();
    let bad = Arc::new(AtomicBool::new(false));
    let rhs = Arc::clone(&kinetics);
    let jac = Arc::clone(&kinetics);
    let sens = Arc::clone(&kinetics);
    let br = Arc::clone(&bad);
    let bj = Arc::clone(&bad);
    let bs = Arc::clone(&bad);
    let initial: Vec<_> = model
        .species
        .iter()
        .map(|s| s.initial_count as f64)
        .collect();
    let problem = OdeBuilder::<diffsol::NalgebraMat<f64>>::new()
        .p(p)
        .rtol(request.rtol)
        .atol([request.atol])
        .sens_rtol(request.rtol)
        .sens_atol([request.atol])
        .use_coloring(false)
        .rhs_sens_implicit(
            move |y, p, _t, out| match rhs.evaluate(y, p) {
                Ok(v) => out.copy_from_slice(&v.rhs),
                Err(_) => {
                    br.store(true, Ordering::Relaxed);
                    out.fill(f64::NAN);
                }
            },
            move |y, p, _t, v, out| match jac.evaluate(y, p) {
                Ok(d) => {
                    for (i, row) in d.state_jacobian.iter().enumerate() {
                        out[i] = row.iter().zip(v).map(|(a, b)| a * b).sum();
                    }
                }
                Err(_) => {
                    bj.store(true, Ordering::Relaxed);
                    out.fill(f64::NAN);
                }
            },
            move |y, p, _t, v, out| match sens.evaluate(y, p) {
                Ok(d) => {
                    for (i, row) in d.parameter_jacobian.iter().enumerate() {
                        out[i] = row.iter().zip(v).map(|(a, b)| a * b).sum();
                    }
                }
                Err(_) => {
                    bs.store(true, Ordering::Relaxed);
                    out.fill(f64::NAN);
                }
            },
        )
        .init_sens(
            move |_, _, y| y.copy_from_slice(&initial),
            |_, _, _, y| y.fill(0.),
            n,
        )
        .build()?;
    let mut solver = problem.bdf_sens::<diffsol::NalgebraLU<f64>>()?;
    let (values, sens, stop) = solver.solve_dense_sensitivities(&request.times)?;
    if bad.load(Ordering::Relaxed) || stop != OdeSolverStopReason::TstopReached {
        return Err("ODE domain error or incomplete solve".into());
    }
    let output = Output {
        schema: OutputVersion::V01,
        solver: "diffsol/0.18.0 BDF",
        convention: "deterministic_counts_macroscopic_not_exact_stochastic_means",
        model_sha256: crn_craft::sha256(&bytes),
        times: request.times.clone(),
        parameters,
        values: (0..request.times.len())
            .map(|j| (0..n).map(|i| values[(i, j)]).collect())
            .collect(),
        sensitivities: sens
            .iter()
            .map(|matrix| {
                (0..request.times.len())
                    .map(|j| (0..n).map(|i| matrix[(i, j)]).collect())
                    .collect()
            })
            .collect(),
    };
    println!("{}", serde_json::to_string(&output)?);
    Ok(())
}
