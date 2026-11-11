//! Prints the 0.6.0 `kinetics` right-hand side and Jacobians of each model
//! file named on the command line as one JSON document on stdout.
//!
//! Used by `validation/merge/references.py` (task M0c of the 0.7.0 engine
//! merge) to record the ODE references before the network engine is ported.
//! Each accepted model is evaluated with its own parameter values at two
//! points: the initial state and a probe state whose species `i` holds its
//! initial count plus `i + 1`, so every rate and derivative term is active.
//! A model that `Kinetics::compile` refuses is reported with its diagnostics.
//! Numbers are printed in shortest round-trip form, so they parse back to
//! the same `f64` values.
//!
//! This program builds only against the 0.6.0 crate API (`parse_model`,
//! `kinetics::Kinetics`); it is a recording tool, not part of the workspace.
use crn_craft::{kinetics::Kinetics, parse_model, sha256};
use serde_json::{Value, json};

fn point(kinetics: &Kinetics, label: &str, y: &[f64], p: &[f64]) -> Value {
    match kinetics.evaluate(y, p) {
        Ok(d) => json!({
            "label": label,
            "y": y,
            "rhs": d.rhs,
            "state_jacobian": d.state_jacobian,
            "parameter_jacobian": d.parameter_jacobian,
        }),
        Err(diagnostic) => json!({"label": label, "y": y, "diagnostic": diagnostic}),
    }
}

fn entry(path: &str) -> Value {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let head = json!({"model": path, "model_sha256": sha256(&bytes)});
    let body = match parse_model(&bytes) {
        Err(diagnostic) => json!({"accepted": false, "diagnostics": [diagnostic]}),
        Ok(model) => match Kinetics::compile(&model) {
            Err(diagnostics) => json!({"accepted": false, "diagnostics": diagnostics}),
            Ok(kinetics) => {
                let p: Vec<f64> = kinetics
                    .parameters
                    .iter()
                    .map(|name| model.parameters[name])
                    .collect();
                let initial: Vec<f64> = model
                    .species
                    .iter()
                    .map(|s| s.initial_count as f64)
                    .collect();
                let probe: Vec<f64> = model
                    .species
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.initial_count + i as i64 + 1) as f64)
                    .collect();
                let species: Vec<&str> = model.species.iter().map(|s| s.name.as_str()).collect();
                json!({
                    "accepted": true,
                    "species": species,
                    "parameters": kinetics.parameters,
                    "parameter_values": p,
                    "points": [
                        point(&kinetics, "initial", &initial, &p),
                        point(&kinetics, "probe", &probe, &p),
                    ],
                })
            }
        },
    };
    let mut merged = head;
    for (key, value) in body.as_object().expect("object") {
        merged[key] = value.clone();
    }
    merged
}

fn main() {
    let models: Vec<Value> = std::env::args().skip(1).map(|path| entry(&path)).collect();
    let document = json!({
        "schema": "crn-craft.merge-kinetics-reference/0.1",
        "crn_craft": crn_craft::describe()["version"],
        "models": models,
    });
    println!(
        "{}",
        serde_json::to_string(&document).expect("finite values")
    );
}
