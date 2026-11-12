//! Finite state projection through `crn-craft.run/1.0` with method `fsp`
//! (M5a). The cases of the network engine's `fsp` tests of 0.6.0 run on
//! converted models; the network `fsp` command was removed in 0.7.0, and
//! `validation/merge/check_fsp_references.py` compares the 1.0 path with
//! its stored 0.6.0 results.
use crn_craft::fsp::v1;
use crn_craft::rules::{Control as RulesControl, convert};
use serde_json::{Value, json};

/// The decay network model of the 0.6.0 tests: four molecules of `A` that
/// die at rate `k`.
fn fixture(k: f64) -> Vec<u8> {
    serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"Decay","convention":"stochastic_counts","omega":1.0,"parameters":{},"species":[{"name":"A","initial_count":4}],"reactions":[{"name":"death","reactants":{"A":1},"products":{},"rate":{"kind":"mass_action","coefficient":k}}]})).unwrap()
}

/// A `run/1.0` FSP request for model bytes.
fn request(model: &[u8], bounds: Value, end: f64, report_times: &[f64]) -> Value {
    json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": report_times},
        "method": {"kind": "fsp", "bounds": bounds, "max_states": 100, "tolerance": 1e-14},
        "budget": {"max_output_bytes": 1000000}
    })
}

/// The converted fixture and its request: the 0.6.0 network run's settings
/// with report times equal to its save times.
fn converted(bound: i64, k: f64, end: f64) -> (Vec<u8>, Value) {
    let network = fixture(k);
    let model = convert::from_network(&network).unwrap();
    let request = request(&model, json!({"A": bound}), end, &[0.0, end]);
    (model, request)
}

fn solve_v1(
    model: &[u8],
    request: &Value,
) -> Result<v1::FspResult, Vec<crn_craft::diagnostic::Diagnostic>> {
    v1::solve(
        model,
        &serde_json::to_vec(request).unwrap(),
        &RulesControl::default(),
    )
}

/// A converted network model with one reaction replaced.
fn with_reaction(reactants: Value, products: Value, initial: i64) -> Vec<u8> {
    let m = fixture(2.0);
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["reactions"][0]["reactants"] = reactants;
    model["reactions"][0]["products"] = products;
    model["species"][0]["initial_count"] = json!(initial);
    convert::from_network(&serde_json::to_vec(&model).unwrap()).unwrap()
}

#[test]
fn v1_binomial_decay_and_large_lambda_t() {
    for end in [0.5, 1000.0] {
        let (m, r) = converted(4, 2.0, end);
        let f = solve_v1(&m, &r).unwrap();
        assert_eq!(f.schema, v1::FspVersion::V10);
        assert_eq!(f.species, ["A"]);
        let q = (-2.0 * end).exp();
        for (state, &p) in f.states.iter().zip(&f.snapshots[1].probabilities) {
            let k = state[0] as i32;
            let binomial = [1.0, 4.0, 6.0, 4.0, 1.0][k as usize];
            let expected = binomial * q.powi(k) * (1.0 - q).powi(4 - k);
            assert!((p - expected).abs() < 1e-12, "{state:?} {p} {expected}");
        }
        assert_eq!(f.snapshots[1].physical_leakage_estimate, 0.0);
        assert!(f.snapshots[1].series_tail_upper_bound <= 1.01e-14);
        assert!(f.snapshots[1].series_tail_upper_bound > 0.0);
        assert!(!f.numerical_error_certified);
    }
}

#[test]
fn v1_leakage_is_not_reflected_or_double_counted() {
    let m = with_reaction(json!({}), json!({"A": 1}), 4);
    let r = request(&m, json!({"A": 4}), 1.0, &[0.0, 1.0]);
    let f = solve_v1(&m, &r).unwrap();
    assert_eq!(f.states.len(), 1);
    let s = &f.snapshots[1];
    assert!((s.probabilities[0] - (-2.0_f64).exp()).abs() < 1e-12);
    assert!((s.physical_leakage_estimate - (1.0 - (-2.0_f64).exp())).abs() < 1e-12);
    assert!(s.mass_defect.abs() < 1e-12);
}

#[test]
fn v1_capacity_cancellation_bounds_and_null_reaction() {
    let (m, r) = converted(4, 1.0, 1.0);
    let mut capped = r.clone();
    capped["method"]["max_states"] = json!(2);
    let errors = solve_v1(&m, &capped).unwrap_err();
    assert_eq!(errors[0].code, "CRN-E070");
    assert!(errors[0].message.contains("cap"), "{}", errors[0].message);
    let control = RulesControl::default();
    control.cancel();
    let errors = v1::solve(&m, &serde_json::to_vec(&r).unwrap(), &control).unwrap_err();
    assert_eq!(errors[0].code, "CRN-E070");
    assert!(errors[0].message.contains("interrupted"));
    let mut below = r.clone();
    below["method"]["bounds"]["A"] = json!(3);
    let errors = solve_v1(&m, &below).unwrap_err();
    assert_eq!(
        (errors[0].code.as_str(), errors[0].path.as_str()),
        ("CRN-R020", "/method/bounds/A")
    );
    let m = with_reaction(json!({"A": 1}), json!({"A": 1}), 4);
    let r = request(&m, json!({"A": 4}), 1.0, &[0.0, 1.0]);
    let f = solve_v1(&m, &r).unwrap();
    assert_eq!(f.uniformization_rate, 0.0);
    assert_eq!(f.states.len(), 1);
    assert_eq!(f.snapshots[1].probabilities, vec![1.0]);
}

#[test]
fn v1_omitted_weight_is_summed_not_cancelled_and_stays_below_the_tail_bound() {
    let (m, r) = converted(4, 2.0, 1000.0);
    let f = solve_v1(&m, &r).unwrap();
    for s in &f.snapshots {
        assert!(s.discarded_weight_estimate >= 0.0);
        assert!(s.discarded_weight_estimate <= s.series_tail_upper_bound);
        assert!(
            s.roundoff_mass_residual < 1e-13,
            "{}",
            s.roundoff_mass_residual
        );
    }
}

#[test]
fn v1_the_tail_bound_bounds_the_missing_mass_at_a_loose_tolerance() {
    let (m, mut r) = converted(4, 2.0, 3.0);
    r["method"]["tolerance"] = json!(0.5);
    let f = solve_v1(&m, &r).unwrap();
    let last = f.snapshots.last().unwrap();
    assert!(
        last.mass_defect > 1e-6,
        "tolerance 0.5 should lose visible mass"
    );
    assert!(last.mass_defect <= last.series_tail_upper_bound + 1e-12);
    assert!(last.discarded_weight_estimate <= last.series_tail_upper_bound);
    assert!((last.discarded_weight_estimate - last.mass_defect).abs() < 1e-12);
}

#[test]
fn v1_many_snapshots_respect_the_output_budget() {
    // One reachable state (a null reaction); each snapshot record
    // serializes to about 240 bytes.
    let m = with_reaction(json!({"A": 1}), json!({"A": 1}), 1);
    let saves: Vec<f64> = (0..20000).map(|i| f64::from(i) * 1e-4).collect();
    let mut r = request(&m, json!({"A": 4}), saves[saves.len() - 1], &saves);
    r["method"]["tolerance"] = json!(1e-10);
    r["budget"]["max_output_bytes"] = json!(2_000_000);
    let errors = solve_v1(&m, &r).unwrap_err();
    assert_eq!(errors[0].code, "CRN-E070");
    assert!(
        errors[0].message.contains("output"),
        "{}",
        errors[0].message
    );
    r["budget"]["max_output_bytes"] = json!(100_000_000);
    let f = solve_v1(&m, &r).unwrap();
    let bytes = serde_json::to_vec(&f).unwrap().len();
    assert!(bytes > 2_000_000, "fixture output is only {bytes} bytes");
}

#[test]
fn v1_stop_states_are_absorbing() {
    // Pure birth at rate 2 that stops at A >= 3: below 3 the Poisson law,
    // at 3 the rest of the mass.
    let model = json!({
        "schema": "crn-craft.model/1.0", "name": "birth",
        "parameters": {"k": 2.0},
        "types": [{"name": "A"}], "initial": [],
        "observables": {"A": {"type": "A"}},
        "rules": [{"name": "birth", "reactants": [], "rate": "k", "products": [{"type": "A"}]}],
        "stop": "A >= 3"
    });
    let m = serde_json::to_vec(&model).unwrap();
    let r = request(&m, json!({"A": 10}), 1.5, &[0.5, 1.5]);
    let f = solve_v1(&m, &r).unwrap();
    assert_eq!(f.states, [[0], [1], [2], [3]]);
    for s in &f.snapshots {
        let mu = 2.0 * s.time;
        let poisson: Vec<f64> = (0..3)
            .map(|j| (-mu).exp() * mu.powi(j) / [1.0, 1.0, 2.0][j as usize])
            .collect();
        for (p, q) in s.probabilities.iter().zip(&poisson) {
            assert!((p - q).abs() < 1e-12);
        }
        let rest = 1.0 - poisson.iter().sum::<f64>();
        assert!((s.probabilities[3] - rest).abs() < 1e-12);
        assert_eq!(s.physical_leakage_estimate, 0.0);
    }
    // A stop that holds initially leaves the state where it is.
    let mut stopped = model.clone();
    stopped["stop"] = json!("A >= 0");
    let m = serde_json::to_vec(&stopped).unwrap();
    let f = solve_v1(&m, &request(&m, json!({"A": 10}), 1.0, &[1.0])).unwrap();
    assert_eq!(f.states, [[0]]);
    assert_eq!(f.snapshots[0].probabilities, [1.0]);
}

#[test]
fn v1_attributed_models_use_flat_labels() {
    // Each wt = 0 copy converts at rate 0.7 on its own: the number converted
    // by t is Binomial(2, 1 - e^(-0.7 t)).
    let model = json!({
        "schema": "crn-craft.rules-model/0.4", "name": "convert",
        "parameters": {"k": 0.7},
        "types": [{"name": "M", "attributes": ["wt"]}],
        "initial": [{"type": "M", "attributes": {"wt": 0}, "count": 2},
                    {"type": "M", "attributes": {"wt": 1}, "count": 1}],
        "observables": {"M": {"type": "M"}},
        "rules": [{"name": "convert", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 0"}],
                   "rate": "k", "products": [{"type": "M", "attributes": {"wt": 1}}]}]
    });
    let m = serde_json::to_vec(&model).unwrap();
    let r = request(&m, json!({"M{wt=0}": 2, "M{wt=1}": 3}), 2.0, &[0.0, 2.0]);
    let f = solve_v1(&m, &r).unwrap();
    assert_eq!(f.species, ["M{wt=0}", "M{wt=1}"]);
    assert_eq!(f.states, [[2, 1], [1, 2], [0, 3]]);
    let q = 1.0 - (-1.4f64).exp();
    let expected = [(1.0 - q).powi(2), 2.0 * q * (1.0 - q), q * q];
    for (p, e) in f.snapshots[1].probabilities.iter().zip(expected) {
        assert!((p - e).abs() < 1e-12, "{p} {e}");
    }
    let mean = f.snapshots[1].projected_first_moment[1];
    assert!((mean - (1.0 + 2.0 * q)).abs() < 1e-12);
}

#[test]
fn v1_requests_are_checked_with_coded_diagnostics() {
    let (m, r) = converted(4, 1.0, 1.0);
    let refused = |request: &Value| {
        let errors = solve_v1(&m, request).unwrap_err();
        (errors[0].code.clone(), errors[0].path.clone())
    };
    let mut old = r.clone();
    old["schema"] = json!("crn-craft.rules-run/0.6");
    assert_eq!(refused(&old), ("CRN-R020".into(), "/schema".into()));
    let mut missing = r.clone();
    missing["method"]["bounds"] = json!({});
    assert_eq!(
        refused(&missing),
        ("CRN-R020".into(), "/method/bounds/A".into())
    );
    let mut unknown = r.clone();
    unknown["method"]["bounds"]["B"] = json!(3);
    assert_eq!(
        refused(&unknown),
        ("CRN-R020".into(), "/method/bounds/B".into())
    );
    let mut negative = r.clone();
    negative["method"]["bounds"]["A"] = json!(-1);
    assert_eq!(refused(&negative).0, "CRN-E003");
    let mut states = r.clone();
    states["method"]["max_states"] = json!(0);
    assert_eq!(
        refused(&states),
        ("CRN-R020".into(), "/method/max_states".into())
    );
    for tolerance in [0.0, 1.0] {
        let mut t = r.clone();
        t["method"]["tolerance"] = json!(tolerance);
        assert_eq!(refused(&t), ("CRN-R020".into(), "/method/tolerance".into()));
    }
    let mut times = r.clone();
    times["time"]["report_times"] = json!([0.5, 0.25]);
    assert_eq!(
        refused(&times),
        ("CRN-R020".into(), "/time/report_times/1".into())
    );
    let mut budget = r.clone();
    budget["budget"]["max_output_bytes"] = json!(0);
    assert_eq!(
        refused(&budget),
        ("CRN-R020".into(), "/budget/max_output_bytes".into())
    );
    budget["budget"] = json!({"max_output_bytes": 1000, "wall_time_s": -1.0});
    assert_eq!(
        refused(&budget),
        ("CRN-R020".into(), "/budget/wall_time_s".into())
    );
    let mut hash = r.clone();
    hash["model"]["sha256"] = json!(crn_craft::sha256(b"other"));
    assert_eq!(refused(&hash), ("CRN-R020".into(), "/model/sha256".into()));
    // Simulation fields are unknown in an FSP request.
    let mut ensemble = r.clone();
    ensemble["ensemble"] = json!({"replicates": 1, "root_seed": 0});
    assert_eq!(refused(&ensemble).0, "CRN-E003");
    // A model without a flat view.
    let growing = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "growing",
        "types": [{"name": "P", "attributes": ["n"]}],
        "initial": [{"type": "P", "attributes": {"n": 1}, "count": 1}],
        "observables": {"P": {"type": "P"}},
        "rules": [{"name": "grow", "reactants": [{"label": "x", "type": "P"}], "rate": "1",
                   "products": [{"type": "P", "copy": "x", "attributes": {"n": "x.n + 1"}}]}]
    }))
    .unwrap();
    let errors = solve_v1(&growing, &request(&growing, json!({}), 1.0, &[1.0])).unwrap_err();
    assert_eq!(
        (errors[0].code.as_str(), errors[0].path.as_str()),
        ("CRN-R017", "/rules")
    );
    // A type that never holds molecules has no label.
    let never = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "never",
        "types": [{"name": "A"}, {"name": "Z"}],
        "initial": [{"type": "A", "count": 2}],
        "observables": {"A": {"type": "A"}},
        "rules": [{"name": "death", "reactants": [{"label": "x", "type": "A"}], "rate": "1", "products": []},
                  {"name": "unused", "reactants": [{"label": "z", "type": "Z"}], "rate": "1", "products": []}]
    }))
    .unwrap();
    let errors = solve_v1(
        &never,
        &request(&never, json!({"A": 2, "Z": 1}), 1.0, &[1.0]),
    )
    .unwrap_err();
    assert_eq!(errors[0].path, "/method/bounds/Z");
    assert!(
        errors[0].message.contains("never holds molecules"),
        "{}",
        errors[0].message
    );
    assert!(solve_v1(&never, &request(&never, json!({"A": 2}), 1.0, &[1.0])).is_ok());
}

#[test]
fn v1_fsp_requests_are_recognized_and_refused_by_prepare() {
    let (m, r) = converted(4, 1.0, 1.0);
    let bytes = serde_json::to_vec(&r).unwrap();
    assert!(v1::is_fsp_request(&bytes));
    assert!(!v1::is_fsp_request(
        b"{\"method\": {\"kind\": \"static_direct\"}}"
    ));
    assert!(!v1::is_fsp_request(b"not json"));
    let errors = match crn_craft::rules::prepare(&m, &bytes) {
        Err(errors) => errors,
        Ok(_) => panic!("an FSP request prepared as a simulation"),
    };
    assert_eq!(
        (errors[0].code.as_str(), errors[0].path.as_str()),
        ("CRN-R020", "/method")
    );
}

#[test]
fn v1_solve_files_resolves_the_model_and_honors_a_cancel_file() {
    let folder = std::env::temp_dir().join(format!("crn-fsp-v1-{}", std::process::id()));
    std::fs::create_dir_all(&folder).unwrap();
    let (m, r) = converted(4, 1.0, 1.0);
    std::fs::write(folder.join("model.json"), &m).unwrap();
    std::fs::write(folder.join("request.json"), serde_json::to_vec(&r).unwrap()).unwrap();
    std::fs::write(folder.join("other.json"), &m).unwrap();
    let control = RulesControl::default();
    let f = v1::solve_files(
        &folder.join("model.json"),
        &folder.join("request.json"),
        &control,
    )
    .unwrap();
    assert_eq!(f.states.len(), 5);
    assert!(v1::is_fsp_request_file(&folder.join("request.json")));
    assert!(!v1::is_fsp_request_file(&folder.join("missing.json")));
    let errors = v1::solve_files(
        &folder.join("other.json"),
        &folder.join("request.json"),
        &control,
    )
    .unwrap_err();
    assert_eq!(errors[0].path, "/model/path");
    std::fs::write(folder.join("cancel"), b"").unwrap();
    let mut control = RulesControl::default();
    control.cancel_file = Some(folder.join("cancel"));
    let errors = v1::solve_files(
        &folder.join("model.json"),
        &folder.join("request.json"),
        &control,
    )
    .unwrap_err();
    assert!(
        errors[0].message.contains("interrupted: cancelled"),
        "{}",
        errors[0].message
    );
    std::fs::remove_dir_all(&folder).unwrap();
}
