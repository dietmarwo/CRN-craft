//! Independent continuous variables (rules-model/0.4) and the shared
//! `hybrid` solution: the contract, determinism across threads, bounds and
//! failures. The statistical gates are in `validation/check_hybrid.py`.
use crn_craft_rules::{Control, RulesResult, Termination, check_model, execute, prepare};
use serde_json::{Value, json};

const V04: &str = "crn-craft.rules-model/0.4";

/// Gompertz mortality: dm/dt = g m drives the per-molecule death rate m of
/// `n` molecules; `x` is a dependent variable that reads `m`.
fn gompertz(independent: bool, n: u64) -> Value {
    json!({
        "schema": V04,
        "name": "gompertz",
        "parameters": {"g": 0.5, "m0": 0.1},
        "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": n}],
        "observables": {"A": {"type": "A"}},
        "rules": [{"name": "death", "reactants": [{"label": "a", "type": "A"}],
                   "rate": "m", "products": []}],
        "continuous": [
            {"name": "m", "initial": "m0", "derivative": "g * m", "independent": independent},
            {"name": "x", "initial": 0, "derivative": "m * A - x"}
        ]
    })
}

fn request(model: &[u8], threads: u32, end: f64, replicates: u64) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": "crn-craft.rules-run/0.5",
        "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": [0.0, end / 4.0, end / 2.0, end]},
        "ensemble": {"replicates": replicates, "root_seed": 31},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates,
                    "keep_instances": 0},
        "budget": {"threads": threads, "max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 10_000_000u64,
                   "max_instances": 10, "max_output_bytes": 1u64 << 26},
        "method": {"kind": "hybrid", "rtol": 1e-8, "atol": 1e-10}
    }))
    .unwrap()
}

fn run(model: &Value, threads: u32, end: f64, replicates: u64) -> RulesResult {
    let bytes = serde_json::to_vec(model).unwrap();
    let prepared = prepare(&bytes, &request(&bytes, threads, end, replicates)).unwrap();
    execute(&prepared, &Control::default()).unwrap()
}

fn first(model: &Value) -> Option<(String, String)> {
    check_model(&serde_json::to_vec(model).unwrap())
        .diagnostics
        .first()
        .map(|d| (d.code.clone(), d.path.clone()))
}

fn pair(code: &str, path: &str) -> Option<(String, String)> {
    Some((code.into(), path.into()))
}

#[test]
fn the_contract_accepts_independent_variables_in_0_4_only() {
    let model = gompertz(true, 10);
    assert_eq!(first(&model), None);
    // `false` is the default and also valid.
    assert_eq!(first(&gompertz(false, 10)), None);
    let mut old = model.clone();
    old["schema"] = json!("crn-craft.rules-model/0.3");
    assert_eq!(first(&old), pair("CRN-R015", "/continuous/0/independent"));
    old["continuous"][0]["independent"] = json!(false);
    assert_eq!(first(&old), pair("CRN-R015", "/continuous/0/independent"));
    let mut null = model.clone();
    null["continuous"][0]["independent"] = Value::Null;
    assert_eq!(first(&null).map(|p| p.0), Some("CRN-E003".into()));
    let mut text = model.clone();
    text["continuous"][0]["independent"] = json!("yes");
    assert_eq!(first(&text).map(|p| p.0), Some("CRN-E003".into()));
    // A 0.4 model without the field is a 0.3 model under a new version.
    let mut plain = model.clone();
    plain["continuous"][0]
        .as_object_mut()
        .unwrap()
        .remove("independent");
    assert_eq!(first(&plain), None);
}

#[test]
fn an_independent_derivative_reads_only_parameters_and_independent_variables() {
    let read = |derivative: &str| {
        let mut model = gompertz(true, 10);
        model["continuous"][0]["derivative"] = json!(derivative);
        model["continuous"]
            .as_array_mut()
            .unwrap()
            .push(json!({"name": "z", "initial": 1, "derivative": "-z", "independent": true}));
        first(&model)
    };
    assert_eq!(read("g * m * z"), None);
    assert_eq!(read("g * m"), None);
    // An observable, a dependent variable, and both in a guarded branch.
    for derivative in ["g * m * A", "g * m + x", "if(z > 2, A, 0)"] {
        assert_eq!(
            read(derivative),
            pair("CRN-R012", "/continuous/0/derivative"),
            "{derivative}"
        );
    }
    let report = {
        let mut model = gompertz(true, 10);
        model["continuous"][0]["derivative"] = json!("g * x");
        check_model(&serde_json::to_vec(&model).unwrap())
    };
    let message = &report.diagnostics[0].message;
    assert!(
        message.contains("dependent continuous variable `x`"),
        "{message}"
    );
    // A dependent variable may read an independent one.
    let mut model = gompertz(true, 10);
    model["continuous"][1]["derivative"] = json!("m - x");
    assert_eq!(first(&model), None);
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let model = gompertz(true, 40);
    let mut one = serde_json::to_value(run(&model, 1, 6.0, 24)).unwrap();
    let mut three = serde_json::to_value(run(&model, 3, 6.0, 24)).unwrap();
    // run_sha256 hashes the request, which names the thread count.
    one["run_sha256"] = Value::Null;
    three["run_sha256"] = Value::Null;
    assert_eq!(one, three);
    // And a second execution of the same prepared run is identical.
    assert_eq!(
        serde_json::to_value(run(&model, 3, 6.0, 24)).unwrap()["paths"],
        three["paths"]
    );
}

#[test]
fn independent_values_are_shared_and_follow_the_ode() {
    let result = run(&gompertz(true, 40), 2, 4.0, 12);
    assert_eq!(result.continuous, ["m", "x"]);
    let column = result.observables.iter().position(|o| o == "m").unwrap();
    for path in &result.paths {
        for (row, &t) in result.report_times.iter().enumerate() {
            let m = path.values[row][column];
            let exact = 0.1 * (0.5 * t).exp();
            assert!((m - exact).abs() <= 1e-6 * exact, "{t}: {m} vs {exact}");
            assert_eq!(m.to_bits(), result.paths[0].values[row][column].to_bits());
        }
    }
    // The dependent variable and the counts still vary by replicate.
    let x = result.observables.iter().position(|o| o == "x").unwrap();
    let last = result.report_times.len() - 1;
    assert!(
        result
            .paths
            .iter()
            .any(|p| p.values[last][x] != result.paths[0].values[last][x])
    );
}

#[test]
fn an_independent_stop_bound_ends_every_replicate_at_the_same_time() {
    let mut model = gompertz(true, 1000);
    model["parameters"]["m0"] = json!(0.001);
    model["continuous"][0] = json!({"name": "m", "initial": 3, "derivative": "-1",
                                    "lower_bound": 0, "on_bound": "stop",
                                    "independent": true});
    model["rules"][0]["rate"] = json!("m0 * m");
    let result = run(&model, 2, 5.0, 16);
    let m = result.observables.iter().position(|o| o == "m").unwrap();
    for (record, path) in result.replicates.iter().zip(&result.paths) {
        assert_eq!(record.termination, Termination::Stopped);
        assert!(
            (record.stop_time - 3.0).abs() <= 1e-9,
            "{}",
            record.stop_time
        );
        // Reports: 0, 1.25, 2.5 before the stop, then the bound.
        assert!((path.values[2][m] - 0.5).abs() <= 1e-9);
        assert_eq!(path.values[3][m], 0.0);
    }
    let first = result.replicates[0].stop_time.to_bits();
    assert!(
        result
            .replicates
            .iter()
            .all(|r| r.stop_time.to_bits() == first)
    );
}

#[test]
fn an_independent_clamp_bound_holds_the_variable() {
    let mut model = gompertz(true, 50);
    model["continuous"][0] = json!({"name": "m", "initial": 2, "derivative": "-1",
                                    "lower_bound": 0, "on_bound": "clamp",
                                    "independent": true});
    let result = run(&model, 2, 8.0, 8);
    let m = result.observables.iter().position(|o| o == "m").unwrap();
    for (record, path) in result.replicates.iter().zip(&result.paths) {
        assert_eq!(record.termination, Termination::Completed);
        assert!((path.values[1][m] - 0.0).abs() <= 1e-9);
        assert_eq!(path.values[2][m], 0.0);
        assert_eq!(path.values[3][m], 0.0);
    }
}

#[test]
fn a_failed_independent_solve_fails_every_replicate_that_reaches_it() {
    // m' = m^2 from 1 blows up at t = 1.
    let mut model = gompertz(true, 3);
    model["continuous"][0]["derivative"] = json!("m * m");
    model["parameters"]["m0"] = json!(1.0);
    model["rules"][0]["rate"] = json!("0 * m + 0.01");
    let result = run(&model, 2, 4.0, 6);
    let first = &result.replicates[0];
    for record in &result.replicates {
        assert_eq!(record.termination, Termination::Failed);
        let diagnostic = serde_json::to_value(record.diagnostic.as_ref().unwrap()).unwrap();
        assert_eq!(
            diagnostic,
            serde_json::to_value(first.diagnostic.as_ref().unwrap()).unwrap()
        );
        assert_eq!(record.stop_time.to_bits(), first.stop_time.to_bits());
        // The step underflows within the requested accuracy of t = 1.
        assert!(
            (record.stop_time - 1.0).abs() < 1e-6,
            "{}",
            record.stop_time
        );
        assert_eq!(record.events, 0);
    }
}
