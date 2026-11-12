//! End-to-end checks of the crn-diffsol binary against closed forms and its failure contract.
use serde_json::{Value, json};
use std::process::Command;

const ENZYME: &str = include_str!("../../../models/legacy/enzyme-mechanism.model.json");
const MM: &str = include_str!("../../../models/legacy/michaelis-menten.model.json");
const IMMIGRATION: &str = include_str!("../../../models/legacy/immigration-death.model.json");
const DIMER: &str = include_str!("../fixtures/dimerization.model.json");
const HILL: &str = include_str!("../fixtures/hill-constant.model.json");

/// Runs the binary; returns the exit code and the stdout result or the one-line stderr diagnostic.
fn solve(key: &str, model: &Value, times: &[f64], rtol: f64, atol: f64) -> (i32, Value) {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("crn-diffsol-{key}"));
    std::fs::create_dir_all(&dir).unwrap();
    let (m, r) = (dir.join("model.json"), dir.join("request.json"));
    std::fs::write(&m, serde_json::to_vec(model).unwrap()).unwrap();
    let request =
        json!({"schema":"crn-craft.ode-request/0.1","times":times,"rtol":rtol,"atol":atol});
    std::fs::write(&r, serde_json::to_vec(&request).unwrap()).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_crn-diffsol"))
        .arg(&m)
        .arg(&r)
        .output()
        .unwrap();
    let code = out.status.code().expect("exit code");
    if code == 0 {
        assert!(out.stderr.is_empty());
        (0, serde_json::from_slice(&out.stdout).unwrap())
    } else {
        assert!(out.stdout.is_empty());
        let text = String::from_utf8(out.stderr).unwrap();
        assert_eq!(text.lines().count(), 1, "{text}");
        let diagnostic: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(diagnostic["success"], false);
        (code, diagnostic)
    }
}
fn parse(text: &str) -> Value {
    serde_json::from_str(text).unwrap()
}
fn column(result: &Value, i: usize) -> Vec<f64> {
    let values = result["values"].as_array().unwrap();
    values.iter().map(|row| row[i].as_f64().unwrap()).collect()
}
fn sensitivity(result: &Value, parameter: &str, i: usize) -> Vec<f64> {
    let names = result["parameters"].as_array().unwrap();
    let j = names.iter().position(|p| p == parameter).unwrap();
    let rows = result["sensitivities"][j].as_array().unwrap();
    rows.iter().map(|row| row[i].as_f64().unwrap()).collect()
}
fn decay(a0: i64) -> Value {
    json!({"schema":"crn-craft.model/0.1","name":"decay","convention":"stochastic_counts","omega":1.0,
        "parameters":{"k":1.0},"species":[{"name":"A","initial_count":a0}],
        "reactions":[{"name":"d","reactants":{"A":1},"products":{},"rate":{"kind":"mass_action","coefficient":"k"}}]})
}

#[test]
fn decay_through_zero_matches_the_exact_exponential() {
    // Declared tolerance: one solver error weight at the species scale, atol + rtol * A0.
    let a0 = 10.0;
    for (h, times) in [
        vec![0., 20., 40.],
        vec![0., 100., 1000.],
        vec![0., 1000., 10000.],
    ]
    .iter()
    .enumerate()
    {
        for (rtol, atol) in [(1e-10, 1e-12), (1e-8, 1e-10), (1e-6, 1e-8), (1e-4, 1e-6)] {
            let key = format!("decay-{h}-{rtol:e}");
            let (code, result) = solve(&key, &decay(10), times, rtol, atol);
            assert_eq!(code, 0, "{times:?} {rtol} {atol}: {result}");
            for (t, a) in times.iter().zip(column(&result, 0)) {
                let exact = a0 * (-t).exp();
                assert!(
                    (a - exact).abs() <= atol + rtol * a0,
                    "t={t} rtol={rtol}: {a} vs {exact}"
                );
            }
        }
    }
    // The forward sensitivity dA/dk = -t A0 exp(-k t) of the same solve.
    let times = [0., 1., 20., 40.];
    let (_, result) = solve("decay-sens", &decay(10), &times, 1e-10, 1e-12);
    for (t, s) in times.iter().zip(sensitivity(&result, "k", 0)) {
        let exact = -t * a0 * (-t).exp();
        assert!((s - exact).abs() <= 1e-7, "t={t}: {s} vs {exact}");
    }
}

#[test]
fn long_horizon_enzyme_and_michaelis_menten_complete_and_conserve() {
    for (model, horizon) in [(ENZYME, 1000.), (ENZYME, 5000.), (MM, 200.), (MM, 1000.)] {
        let model = parse(model);
        let times = [0., horizon / 2., horizon];
        let (rtol, atol) = (1e-10, 1e-12);
        let (code, result) = solve(&format!("long-{horizon}"), &model, &times, rtol, atol);
        assert_eq!(code, 0, "{horizon}: {result}");
        let names: Vec<_> = model["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["name"].as_str().unwrap())
            .collect();
        let get = |row: &Value, name: &str| {
            row[names.iter().position(|n| *n == name).unwrap()]
                .as_f64()
                .unwrap()
        };
        for row in result["values"].as_array().unwrap() {
            let substrate = get(row, "S") + get(row, "P");
            let total = if names.contains(&"ES") {
                assert!((get(row, "E") + get(row, "ES") - 10.).abs() <= 1e-8);
                substrate + get(row, "ES")
            } else {
                substrate
            };
            assert!((total - 50.).abs() <= 1e-8, "{row}");
            // Every species is non-negative within one error weight of its scale (50).
            for v in row.as_array().unwrap() {
                assert!(v.as_f64().unwrap() >= -(atol + rtol * 50.), "{row}");
            }
        }
        let last = &result["values"][2];
        assert!((get(last, "P") - 50.).abs() <= 1e-6, "{last}");
    }
}

#[test]
fn dimerization_uses_the_macroscopic_factorial_convention() {
    // c * A^2 / 2! per dimerization event and two A per event: dA/dt = -c A^2.
    let (c, a0) = (0.01, 100.0);
    let times = [0., 0.5, 1., 2., 10.];
    let (code, result) = solve("dimer", &parse(DIMER), &times, 1e-10, 1e-12);
    assert_eq!(code, 0, "{result}");
    let a = column(&result, 0);
    let b = column(&result, 1);
    assert!((a[2] - 50.).abs() <= 1e-7, "A(1) = {}", a[2]);
    for (k, t) in times.iter().enumerate() {
        let exact = a0 / (1. + c * a0 * t);
        assert!((a[k] - exact).abs() <= 1e-7 * exact.max(1.), "t={t}");
        assert!((b[k] - (a0 - exact) / 2.).abs() <= 1e-7 * a0, "t={t}");
        let ds = -a0 * a0 * t / (1. + c * a0 * t).powi(2);
        let s = sensitivity(&result, "c", 0)[k];
        assert!(
            (s - ds).abs() <= 1e-5 * ds.abs().max(1.),
            "t={t}: {s} vs {ds}"
        );
        let sb = sensitivity(&result, "c", 1)[k];
        assert!((sb + ds / 2.).abs() <= 1e-5 * ds.abs().max(1.), "t={t}");
    }
}

#[test]
fn hill_activation_and_repression_match_closed_forms() {
    let model = parse(HILL);
    let p = &model["parameters"];
    let get = |k: &str| p[k].as_f64().unwrap();
    let (v, b, k, n, g, r) = (get("v"), get("b"), get("K"), get("n"), get("g"), 20.0);
    let f = 1. / (1. + (k / r).powf(n));
    let df_dk = -n / k * f * (1. - f);
    let df_dn = f * (1. - f) * (r / k).ln();
    let times = [0., 0.5, 2., 10.];
    let (code, result) = solve("hill", &model, &times, 1e-10, 1e-12);
    assert_eq!(code, 0, "{result}");
    // Species order R, X (activated), Y (repressed); sign is +1 for activation, -1 for repression.
    for (i, sign) in [(1usize, 1.0), (2, -1.0)] {
        let fraction = if sign > 0. { f } else { 1. - f };
        let a = b + v * fraction;
        for (j, &t) in times.iter().enumerate() {
            let shape = -(-g * t).exp_m1() / g;
            let expected = [
                ("value", a * shape, column(&result, i)[j]),
                ("v", fraction * shape, sensitivity(&result, "v", i)[j]),
                ("b", shape, sensitivity(&result, "b", i)[j]),
                (
                    "K",
                    sign * v * df_dk * shape,
                    sensitivity(&result, "K", i)[j],
                ),
                (
                    "n",
                    sign * v * df_dn * shape,
                    sensitivity(&result, "n", i)[j],
                ),
                (
                    "g",
                    -a / g * shape + a * t * (-g * t).exp() / g,
                    sensitivity(&result, "g", i)[j],
                ),
            ];
            for (what, exact, actual) in expected {
                assert!(
                    (actual - exact).abs() <= 1e-6 * exact.abs().max(1.),
                    "species {i} {what} t={t}: {actual} vs {exact}"
                );
            }
        }
    }
}

#[test]
fn population_limit_is_a_safety_stop_with_a_diagnostic() {
    let mut model = parse(IMMIGRATION);
    model["parameters"]["death"] = json!(0.01);
    model["species"][0]["population_limit"] = json!(5);
    let (code, d) = solve("limit", &model, &[0., 1., 2.], 1e-8, 1e-10);
    assert_eq!(code, 1);
    assert_eq!(d["code"], "CRN-E044", "{d}");
    assert!(
        d["message"]
            .as_str()
            .unwrap()
            .contains("population_limit 5")
    );
    // A limit that is never reached leaves the solve untouched.
    model["species"][0]["population_limit"] = json!(1000);
    let (code, _) = solve("limit-ok", &model, &[0., 1., 2.], 1e-8, 1e-10);
    assert_eq!(code, 0);
    // A transient between save times is found by the root function: B(t) = 100 t exp(-t) peaks at
    // 36.8 near t=1, while B(0)=0 and B(10)=0.45 are both below the limit 30.
    let chain = json!({"schema":"crn-craft.model/0.1","name":"chain","convention":"stochastic_counts","omega":1.0,
        "parameters":{"k":1.0},"species":[{"name":"A","initial_count":100},{"name":"B","initial_count":0,"population_limit":30}],
        "reactions":[{"name":"ab","reactants":{"A":1},"products":{"B":1},"rate":{"kind":"mass_action","coefficient":"k"}},
                     {"name":"b","reactants":{"B":1},"products":{},"rate":{"kind":"mass_action","coefficient":"k"}}]});
    let (code, d) = solve("limit-transient", &chain, &[0., 10.], 1e-8, 1e-10);
    assert_eq!((code, d["code"].as_str()), (1, Some("CRN-E044")), "{d}");
    assert!(d["message"].as_str().unwrap().contains("species B"));
}

#[test]
fn huge_times_and_evaluation_failures_return_diagnostics() {
    let (code, d) = solve("huge", &parse(IMMIGRATION), &[0., f64::MAX], 1e-8, 1e-10);
    assert_eq!(code, 1, "{d}");
    assert_eq!(d["code"], "CRN-E080");
    // The kinetics diagnostic is reported instead of the solver failure it causes.
    let singular = json!({"schema":"crn-craft.model/0.1","name":"singular","convention":"stochastic_counts","omega":1.0,
        "parameters":{"v":2.0},"species":[{"name":"X","initial_count":0},{"name":"Y","initial_count":1}],
        "reactions":[{"name":"h","reactants":{},"products":{"Y":1},
            "rate":{"kind":"hill_activation","regulator":"X","vmax":"v","basal":0.0,"k":5.0,"n":0.5}}]});
    let (code, d) = solve("singular", &singular, &[0., 1.], 1e-8, 1e-10);
    assert_eq!(code, 1);
    assert!(
        d["message"].as_str().unwrap().contains("singular at zero"),
        "{d}"
    );
}

#[test]
fn rules_models_solve_like_their_network_models() {
    // The converted model/1.0 has the network's species, parameters and kinetics, so the BDF
    // solve gives the same values and sensitivities bit for bit.
    for (key, text) in [("dimer", DIMER), ("hill", HILL), ("enzyme", ENZYME)] {
        let converted = crn_craft::rules::convert::from_network(text.as_bytes()).unwrap();
        let converted: Value = serde_json::from_slice(&converted).unwrap();
        assert_eq!(converted["schema"], "crn-craft.model/1.0");
        let times = [0., 0.5, 2., 10.];
        let (code, network) = solve(
            &format!("{key}-network"),
            &parse(text),
            &times,
            1e-10,
            1e-12,
        );
        assert_eq!(code, 0, "{network}");
        let (code, rules) = solve(&format!("{key}-rules"), &converted, &times, 1e-10, 1e-12);
        assert_eq!(code, 0, "{rules}");
        for field in ["parameters", "values", "sensitivities", "times", "solver"] {
            assert_eq!(network[field], rules[field], "{key} {field}");
        }
        assert_ne!(network["model_sha256"], rules["model_sha256"]);
    }
}

#[test]
fn expression_rates_of_a_model_1_0_match_the_closed_form() {
    // dA/dt = -k exp(-B / 10) A with B held at 5: A(t) = 10 exp(-k e^-0.5 t), and
    // dA/dk = -t e^-0.5 A(t).
    let model = json!({"schema":"crn-craft.model/1.0","name":"damped","parameters":{"k":0.8},
        "types":[{"name":"A"},{"name":"B"}],
        "initial":[{"type":"A","count":10},{"type":"B","count":5}],
        "observables":{"A":{"type":"A"},"B":{"type":"B"}},
        "rules":[{"name":"decay","reactants":[{"label":"x","type":"A"}],
                  "rate":"k * exp(-B / 10)","products":[]}]});
    let times = [0., 0.5, 2., 10.];
    let (code, result) = solve("damped", &model, &times, 1e-10, 1e-12);
    assert_eq!(code, 0, "{result}");
    let rate = 0.8 * (-0.5_f64).exp();
    for (j, &t) in times.iter().enumerate() {
        let a = 10. * (-rate * t).exp();
        assert!((column(&result, 0)[j] - a).abs() <= 1e-8, "t={t}");
        assert_eq!(column(&result, 1)[j], 5.0);
        let s = -t * (-0.5_f64).exp() * a;
        assert!((sensitivity(&result, "k", 0)[j] - s).abs() <= 1e-7, "t={t}");
    }
    // A whole-rule propensity that does not vanish with its reactant is refused.
    let mut gated = model.clone();
    gated["rules"][0] = json!({"name":"leak","reactants":[{"label":"x","type":"A"}],
        "propensity":"k + B","products":[]});
    let (code, d) = solve("gated", &gated, &[0., 1.], 1e-8, 1e-10);
    assert_eq!((code, d["code"].as_str()), (1, Some("CRN-E080")), "{d}");
    assert!(
        d["message"].as_str().unwrap().contains("availability gate"),
        "{d}"
    );
}
