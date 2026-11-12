#![forbid(unsafe_code)]
//! Optional BDF solve of the deterministic count kinetics with forward sensitivities.
//!
//! The model is a network `crn-craft.model/0.1`, a rules model or a `crn-craft.model/1.0`; the
//! kinetics read it through its flat view (`Kinetics::from_bytes`). Population limits exist in
//! network models only.
//!
//! State-domain rule (docs/campaign/M7.md):
//! - Rates and parameter derivatives are evaluated at `max(y, 0)`, and the state Jacobian is the
//!   Jacobian of that extension (a negative component has no effect). Every channel that consumes
//!   a species has a zero rate when that species is zero (a rate rule's propensity contains
//!   y^k / k!, and the kinetics accept a whole-rule propensity that takes molecules only when it
//!   is 0 with each taken species, as Michaelis-Menten is; Hill laws are production-only), so the
//!   exact trajectory from non-negative counts stays non-negative and is unchanged by this
//!   extension. It only gives BDF predictor/Newton trial states a defined value.
//! - The tolerance of species i is one solver error weight at its scale, `atol + rtol * m_i`,
//!   where `m_i` is the largest value of species i seen so far (accepted states, root-search
//!   points and reported values). An accepted state or reported value below `-(atol + rtol * m_i)`
//!   is CRN-E080. A value above `population_limit + atol + rtol * population_limit` is CRN-E044,
//!   the SSA safety-stop code. Crossings between accepted steps are located with diffsol root
//!   functions; reported values are checked again. A failed solve returns no values, and
//!   reported values are never clipped.
use crn_craft::{diagnostic::Diagnostic, kinetics::Kinetics};
use diffsol::{OdeBuilder, OdeSolverMethod, OdeSolverStopReason, SensitivitiesOdeSolverMethod};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    sync::{Arc, Mutex},
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
const HINT: &str = "Use supported deterministic kinetics and finite non-negative states/parameters; ODE trajectories are not exact nonlinear SSA means.";
fn failure(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error("CRN-E080", "/ode", message, HINT)
}
/// Message of the last panic, recorded by the hook so that `main` can report it as a diagnostic.
static PANIC: Mutex<Option<String>> = Mutex::new(None);
fn main() {
    std::panic::set_hook(Box::new(|info| {
        let payload = info
            .payload()
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| info.payload().downcast_ref::<String>().cloned())
            .unwrap_or_default();
        let place = info
            .location()
            .map(|l| {
                let file = std::path::Path::new(l.file()).file_name();
                let file = file.map_or(l.file().into(), |f| f.to_string_lossy());
                format!(" at {file}:{}", l.line())
            })
            .unwrap_or_default();
        if let Ok(mut slot) = PANIC.lock() {
            *slot = Some(format!("{payload}{place}"));
        }
    }));
    let error = match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => return,
        Ok(Err(error)) => error,
        Err(_) => failure(format!(
            "internal ODE solver panic: {}",
            PANIC
                .lock()
                .ok()
                .and_then(|mut slot| slot.take())
                .unwrap_or_default()
        )),
    };
    eprintln!(
        "{}",
        serde_json::json!({"success":false,"code":error.code,"message":error.message})
    );
    std::process::exit(1);
}
/// Diagnostic for a state that leaves the accepted band. `limits` holds (species index, limit)
/// of the species with a population limit, in species order.
fn domain_error(
    names: &[String],
    limits: &[(usize, i64)],
    root: usize,
    value: f64,
    bound: f64,
    t: f64,
) -> Diagnostic {
    let n = names.len();
    if root < n {
        let name = &names[root];
        return failure(format!(
            "species {name} reaches {value:e} below -{bound:e} (atol + rtol * its largest value) at t={t:e}: outside the non-negative domain beyond the solver tolerance"
        ));
    }
    let (index, limit) = limits[root - n];
    Diagnostic::error(
        "CRN-E044",
        &format!("/species/{index}/population_limit"),
        format!(
            "deterministic trajectory of species {} reaches {value:e}, above population_limit {limit} by more than {bound:e} (atol + rtol * limit), at t={t:e}; ODE output withheld",
            names[index]
        ),
        "Stop or explicitly revise the safety limit; the ODE never clips a state at a limit.",
    )
}

/// Whether the document names the network model contract (population limits live there only).
fn is_network_model(bytes: &[u8]) -> bool {
    serde_json::from_slice::<serde_json::Value>(bytes)
        .ok()
        .and_then(|v| {
            v.get("schema")?
                .as_str()
                .map(|s| s == "crn-craft.model/0.1")
        })
        .unwrap_or(false)
}
fn run() -> Result<(), Diagnostic> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() == 2 && args[1] == "describe" {
        println!(
            "{}",
            serde_json::json!({"schema":"crn-craft.ode-discovery/0.1","request_schema":schema_for!(Request),"result_schema":schema_for!(Output)})
        );
        return Ok(());
    }
    if args.len() != 3 {
        return Err(failure("usage: crn-diffsol MODEL REQUEST"));
    }
    let read = |path: &str| std::fs::read(path).map_err(|e| failure(format!("{path}: {e}")));
    let bytes = read(&args[1])?;
    let network = if is_network_model(&bytes) {
        Some(crn_craft::parse_model(&bytes).map_err(|e| failure(e.message))?)
    } else {
        None
    };
    let request: Request =
        crn_craft::strict_json::parse(&read(&args[2])?).map_err(|e| failure(e.message))?;
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
        return Err(failure("invalid ODE request"));
    }
    let kinetics = Arc::new(
        Kinetics::from_bytes(&bytes).map_err(|e| failure(serde_json::to_string(&e).unwrap()))?,
    );
    let parameters = kinetics.parameters.clone();
    let p = kinetics.parameter_values().to_vec();
    if p.is_empty() {
        return Err(failure("sensitivity path requires at least one parameter"));
    }
    let n = kinetics.species().len();
    // A network model keeps its species order and names in the kinetics (names the rules
    // language reserves are renamed there); diagnostics use the model's own names.
    let names: Vec<String> = match &network {
        Some(model) => model.species.iter().map(|s| s.name.clone()).collect(),
        None => kinetics.species().to_vec(),
    };
    if names.len() != n {
        return Err(failure("species of the model and of its kinetics differ"));
    }
    let named_limits: Vec<(usize, i64)> = network
        .iter()
        .flat_map(|model| model.species.iter().enumerate())
        .filter_map(|(i, s)| s.population_limit.map(|l| (i, l)))
        .collect();
    let (atol, rtol) = (request.atol, request.rtol);
    // (species index, limit, accepted excess above the limit)
    let limits: Vec<(usize, f64, f64)> = named_limits
        .iter()
        .map(|&(i, l)| (i, l as f64, atol + rtol * l as f64))
        .collect();
    let roots = n + limits.len();
    let root_limits = limits.clone();
    let jac_limits = limits.clone();
    // Largest value of each species seen by the root function (initial, accepted and search states).
    let scale = Arc::new(Mutex::new(vec![0.0_f64; n]));
    let root_scale = Arc::clone(&scale);
    // The first evaluation failure is kept; any failure fails the solve even if diffsol recovers.
    let first: Arc<Mutex<Option<Diagnostic>>> = Arc::new(Mutex::new(None));
    let evaluator = |first: &Arc<Mutex<Option<Diagnostic>>>| {
        let kinetics = Arc::clone(&kinetics);
        let first = Arc::clone(first);
        move |y: &[f64], p: &[f64]| {
            let clamped: Vec<f64> = y.iter().map(|&v| if v < 0. { 0. } else { v }).collect();
            kinetics.evaluate(&clamped, p).map_err(|d| {
                let mut slot = first.lock().unwrap();
                if slot.is_none() {
                    *slot = Some(d);
                }
            })
        }
    };
    let rhs = evaluator(&first);
    let jac = evaluator(&first);
    let sens = evaluator(&first);
    let initial: Vec<_> = kinetics.initial().iter().map(|&c| c as f64).collect();
    let problem = OdeBuilder::<diffsol::NalgebraMat<f64>>::new()
        .p(p)
        .rtol(request.rtol)
        .atol([request.atol])
        .sens_rtol(request.rtol)
        .sens_atol([request.atol])
        .use_coloring(false)
        .rhs_sens_implicit(
            move |y, p, _t, out| match rhs(y, p) {
                Ok(v) => out.copy_from_slice(&v.rhs),
                Err(()) => out.fill(f64::NAN),
            },
            move |y, p, _t, v, out| match jac(y, p) {
                // Jacobian of the solved extension f(max(y, 0)): a negative component has no effect.
                Ok(d) => {
                    for (i, row) in d.state_jacobian.iter().enumerate() {
                        out[i] = row
                            .iter()
                            .zip(v.iter().zip(y))
                            .map(|(a, (b, y))| if *y < 0. { 0. } else { a * b })
                            .sum();
                    }
                }
                Err(()) => out.fill(f64::NAN),
            },
            move |y, p, _t, v, out| match sens(y, p) {
                Ok(d) => {
                    for (i, row) in d.parameter_jacobian.iter().enumerate() {
                        out[i] = row.iter().zip(v).map(|(a, b)| a * b).sum();
                    }
                }
                Err(()) => out.fill(f64::NAN),
            },
        )
        .root_sens_implicit(
            move |y, _p, _t, g| {
                let mut scale = root_scale.lock().unwrap();
                for i in 0..n {
                    scale[i] = scale[i].max(y[i]);
                    g[i] = y[i] + atol + rtol * scale[i];
                }
                for (k, &(i, limit, excess)) in root_limits.iter().enumerate() {
                    g[n + k] = limit + excess - y[i];
                }
            },
            move |_y, _p, _t, v, g| {
                g[..n].copy_from_slice(v);
                for (k, &(i, _, _)) in jac_limits.iter().enumerate() {
                    g[n + k] = -v[i];
                }
            },
            |_y, _p, _t, _v, g| g.fill(0.),
            roots,
        )
        .init_sens(
            move |_, _, y| y.copy_from_slice(&initial),
            |_, _, _, y| y.fill(0.),
            n,
        )
        .build()
        .map_err(|e| failure(format!("ODE problem: {e}")))?;
    let mut solver = problem
        .bdf_sens::<diffsol::NalgebraLU<f64>>()
        .map_err(|e| failure(format!("ODE solver: {e}")))?;
    let solved = solver.solve_dense_sensitivities(&request.times);
    // An evaluation diagnostic explains a later solver error, so it is reported first.
    if let Some(diagnostic) = first.lock().unwrap().take() {
        return Err(diagnostic);
    }
    let (values, sens, stop) = solved.map_err(|e| failure(e.to_string()))?;
    match stop {
        OdeSolverStopReason::TstopReached => {}
        OdeSolverStopReason::RootFound(t, root) => {
            let state = solver.interpolate(t).map_err(|e| failure(e.to_string()))?;
            let (i, bound) = if root < n {
                (root, atol + rtol * scale.lock().unwrap()[root])
            } else {
                (limits[root - n].0, limits[root - n].2)
            };
            return Err(domain_error(
                &names,
                &named_limits,
                root,
                state[i],
                bound,
                t,
            ));
        }
        OdeSolverStopReason::InternalTimestep => return Err(failure("incomplete ODE solve")),
    }
    let times = request.times.len();
    let mut scale = scale.lock().unwrap().clone();
    for (i, m) in scale.iter_mut().enumerate() {
        *m = (0..times).map(|j| values[(i, j)]).fold(*m, f64::max);
    }
    for j in 0..times {
        for i in 0..n {
            let value = values[(i, j)];
            let bound = atol + rtol * scale[i];
            if !value.is_finite() || value < -bound {
                return Err(domain_error(
                    &names,
                    &named_limits,
                    i,
                    value,
                    bound,
                    request.times[j],
                ));
            }
        }
        for (k, &(i, limit, excess)) in limits.iter().enumerate() {
            if values[(i, j)] > limit + excess {
                let t = request.times[j];
                return Err(domain_error(
                    &names,
                    &named_limits,
                    n + k,
                    values[(i, j)],
                    excess,
                    t,
                ));
            }
        }
    }
    if sens
        .iter()
        .any(|m| (0..times).any(|j| (0..n).any(|i| !m[(i, j)].is_finite())))
    {
        return Err(failure("non-finite forward sensitivity"));
    }
    let output = Output {
        schema: OutputVersion::V01,
        solver: "diffsol/0.18.0 BDF",
        convention: "deterministic_counts_macroscopic_not_exact_stochastic_means",
        model_sha256: crn_craft::sha256(&bytes),
        times: request.times.clone(),
        parameters,
        values: (0..times)
            .map(|j| (0..n).map(|i| values[(i, j)]).collect())
            .collect(),
        sensitivities: sens
            .iter()
            .map(|matrix| {
                (0..times)
                    .map(|j| (0..n).map(|i| matrix[(i, j)]).collect())
                    .collect()
            })
            .collect(),
    };
    let text = serde_json::to_string(&output).map_err(|e| failure(e.to_string()))?;
    println!("{text}");
    Ok(())
}
