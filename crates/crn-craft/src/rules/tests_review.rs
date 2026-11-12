//! Regression tests from the review of the mtDNA A0–A6 patch.
use super::*;
use serde_json::{Value, json};

const THREE_CLASS: &str = include_str!("../../../../models/rules/mtdna-three-class.rules.json");

fn three_class() -> Value {
    serde_json::from_str(THREE_CLASS).unwrap()
}

fn request_v(schema: &str, model: &[u8], replicates: u64, end: f64, times: &[f64]) -> Value {
    json!({
        "schema": schema,
        "model": {"path": "m.rules.json", "sha256": crate::sha256(model)},
        "time": {"end": end, "report_times": times},
        "ensemble": {"replicates": replicates, "root_seed": 23},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 10_000_000u64,
                   "max_sampling_steps_per_replicate": 100_000_000u64,
                   "max_instances": 10_000u64, "threads": 2, "max_output_bytes": 1u64 << 30}
    })
}

fn prepare_values(model: &Value, request: &Value) -> Result<Prepared, Vec<Diagnostic>> {
    prepare(
        &serde_json::to_vec(model).unwrap(),
        &serde_json::to_vec(request).unwrap(),
    )
}

/// Prepare `model` with a request edited by `edit`; the request hash is kept.
fn prepare_with(
    model: &Value,
    schema: &str,
    edit: impl FnOnce(&mut Value),
) -> Result<Prepared, Vec<Diagnostic>> {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut req = request_v(schema, &bytes, 8, 50.0, &[0.0, 25.0, 50.0]);
    edit(&mut req);
    prepare_values(model, &req)
}

fn first_code(result: Result<Prepared, Vec<Diagnostic>>) -> Option<(String, String)> {
    result.err().map(|d| (d[0].code.clone(), d[0].path.clone()))
}

const V01: &str = "crn-craft.rules-run/0.1";
const V02: &str = "crn-craft.rules-run/0.2";

#[test]
fn method_objects_reject_unknown_fields() {
    let model = three_class();
    for method in [
        json!({"kind": "static_direct", "epsilon": 0.5}),
        json!({"kind": "tree", "bogus": 1}),
        json!({"kind": "next_reaction", "step": 1.0}),
        json!({"kind": "thinning", "x": 0}),
    ] {
        let result = prepare_with(&model, V02, |r| r["method"] = method.clone());
        assert!(result.is_err(), "accepted {method}");
    }
    assert!(
        prepare_with(&model, V02, |r| r["method"] =
            json!({"kind": "static_direct"}))
        .is_ok()
    );
}

#[test]
fn version_01_contracts_reject_new_keys_even_when_empty() {
    let model = three_class();
    for (key, value) in [
        ("trackers", json!([])),
        ("method", Value::Null),
        ("event_log", Value::Null),
    ] {
        let result = prepare_with(&model, V01, |r| r[key] = value.clone());
        assert!(result.is_err(), "rules-run/0.1 accepted {key}");
    }
    // A present value must not be null in later versions either.
    assert!(prepare_with(&model, V02, |r| r["method"] = Value::Null).is_err());
    let mut with_atp = three_class();
    with_atp["atp"] = Value::Null;
    assert!(!check_model(&serde_json::to_vec(&with_atp).unwrap()).valid);
}

#[test]
fn method_and_log_diagnostics_name_the_right_field() {
    let model = three_class();
    let thinning = prepare_with(&model, V02, |r| r["method"] = json!({"kind": "thinning"}));
    let diagnostics = thinning.err().unwrap();
    assert_eq!(diagnostics[0].path, "/method");
    assert!(
        diagnostics[0].message.contains("thinning requires an ATP"),
        "{}",
        diagnostics[0].message
    );
    let log = prepare_with(&model, V02, |r| {
        r["event_log"] = json!({"max_entries": 4, "rules": ["no_such_rule"]});
    });
    assert_eq!(
        first_code(log),
        Some(("CRN-R020".into(), "/event_log".into()))
    );
}

const V03: &str = "crn-craft.rules-run/0.3";
const SPECIALIZED: &str = "chacha8-u53-ziggurat/rules-mtdna-direct/1";

fn with_parameters(edit: &[(&str, f64)]) -> Value {
    let mut model = three_class();
    for (name, value) in edit {
        model["parameters"][*name] = json!(*value);
    }
    model
}

fn solve_values(model: &Value, request: &Value) -> RulesResult {
    solve(
        &serde_json::to_vec(model).unwrap(),
        &serde_json::to_vec(request).unwrap(),
        &Control::default(),
    )
    .unwrap()
}

#[test]
fn mtdna_infinite_wait_runs_to_the_horizon_like_static_direct() {
    // A hazard so small that the unit exponential divided by it overflows.
    let text = THREE_CLASS.replace("max(0, K - N) / D * (sm)", "max(0, K - N) / D * (sw)");
    let mut model: Value = serde_json::from_str(&text).unwrap();
    for (name, value) in [("K", 1.0), ("w0", 1.0), ("sw", 0.0), ("delta", 2.3e-308)] {
        model["parameters"][name] = json!(value);
    }
    let bytes = serde_json::to_vec(&model).unwrap();
    for method in ["mtdna_direct", "static_direct"] {
        let mut req = request_v(V03, &bytes, 2000, 1.0, &[0.0, 1.0]);
        req["method"] = json!({"kind": method});
        let result = solve_values(&model, &req);
        if method == "mtdna_direct" {
            assert_eq!(result.sampler, SPECIALIZED);
        }
        assert!(
            result
                .replicates
                .iter()
                .all(|r| r.termination == Termination::Completed),
            "{method} failed a path"
        );
    }
}

#[test]
fn mtdna_certification_covers_realistic_parameters_and_explains_fallback() {
    for edit in [
        vec![],
        vec![("K", 1e5), ("w0", 1e5)],
        vec![("fracAdv", 1e-6)],
        vec![("fracAdv", 0.0)],
        vec![("mu", 1e-7)],
    ] {
        let model = with_parameters(&edit);
        let prepared = prepare_with(&model, V03, |r| {
            r["method"] = json!({"kind": "mtdna_direct"})
        })
        .unwrap();
        assert!(prepared.mtdna.is_some(), "not certified: {edit:?}");
    }
    // A stop written in another form is not recognized; Static Direct runs and
    // the result says why.
    let mut model = three_class();
    model["stop"] = json!("mutants >= N * 0.9 && N > 0");
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request_v(V03, &bytes, 4, 10.0, &[0.0, 10.0]);
    req["method"] = json!({"kind": "mtdna_direct"});
    let result = solve_values(&model, &req);
    assert_eq!(result.sampler, "chacha8-u53/rules-static-direct/1");
    let warning = result
        .diagnostics
        .iter()
        .find(|d| d.code == "CRN-R051")
        .expect("fallback warning");
    assert!(warning.message.contains("stop"), "{}", warning.message);
}

struct Continue;
impl observe::EventObserver for Continue {
    fn on_event(&self, _: &observe::Event<'_>) -> observe::ObserverAction {
        observe::ObserverAction::Continue
    }
}

fn tracked_request(model: &Value, replicates: u64, threads: u32) -> Value {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut req = request_v(
        V03,
        &bytes,
        replicates,
        2000.0,
        &[0.0, 500.0, 1000.0, 2000.0],
    );
    req["method"] = json!({"kind": "mtdna_direct"});
    req["budget"]["threads"] = json!(threads);
    req["trackers"] = json!([
        {"observable": "wild_type"},
        {"observable": "neutral_mutants", "keep_episodes": 2},
        {"observable": "advantaged_mutants"}
    ]);
    req
}

#[test]
fn mtdna_continue_observer_and_threads_preserve_results() {
    let mut model = with_parameters(&[("mu", 3e-3)]);
    model["stop"] = json!("N > 0 && mutants >= 0.9 * N");
    let bytes = serde_json::to_vec(&model).unwrap();
    let req = tracked_request(&model, 24, 1);
    let prepared = prepare(&bytes, &serde_json::to_vec(&req).unwrap()).unwrap();
    assert!(prepared.mtdna.is_some());
    let plain = execute(&prepared, &Control::default()).unwrap();
    let observed = execute_with_observer(&prepared, &Control::default(), &Continue).unwrap();
    assert!(
        serde_json::to_vec(&plain).unwrap() == serde_json::to_vec(&observed).unwrap(),
        "a Continue observer changed the result"
    );
    let mut many = solve_values(&model, &tracked_request(&model, 24, 3));
    let mut one = plain;
    one.run_sha256.clear();
    many.run_sha256.clear();
    assert!(
        serde_json::to_vec(&one).unwrap() == serde_json::to_vec(&many).unwrap(),
        "worker count changed the result"
    );
}

#[test]
fn mtdna_count_trackers_match_the_general_observation_path() {
    // An event log forces the general per-event observer; counts, crossings
    // and extrema must agree exactly, integrals to rounding.
    let model = with_parameters(&[("mu", 3e-3)]);
    let fast = solve_values(&model, &tracked_request(&model, 16, 2));
    let mut req = tracked_request(&model, 16, 2);
    req["event_log"] = json!({"max_entries": 0});
    let general = solve_values(&model, &req);
    assert_eq!(fast.sampler, SPECIALIZED);
    assert_eq!(general.sampler, SPECIALIZED);
    for (a, b) in fast.observations.iter().zip(&general.observations) {
        for (x, y) in a.trackers.iter().zip(&b.trackers) {
            assert_eq!(x.appearances, y.appearances);
            assert_eq!(x.extinctions, y.extinctions);
            assert_eq!(x.minimum.to_bits(), y.minimum.to_bits());
            assert_eq!(x.maximum.to_bits(), y.maximum.to_bits());
            assert!((x.integral - y.integral).abs() <= 1e-12 * x.integral.abs().max(1.0));
        }
    }
}

#[test]
fn mtdna_stop_at_a_report_time_holds_the_state() {
    let mut model = with_parameters(&[("mu", 3e-2)]);
    model["stop"] = json!("N > 0 && mutants >= 0.9 * N");
    let bytes = serde_json::to_vec(&model).unwrap();
    for method in ["mtdna_direct", "static_direct"] {
        let mut req = request_v(V03, &bytes, 32, 20000.0, &[0.0, 5000.0, 10000.0, 20000.0]);
        req["method"] = json!({"kind": method});
        let result = solve_values(&model, &req);
        for r in &result.replicates {
            if r.termination == Termination::Stopped {
                let path = result.paths.iter().find(|p| p.id == r.id).unwrap();
                let last = path.values.last().unwrap();
                // Every report after the stop repeats the stopped state.
                for (row, &t) in path.values.iter().zip(&[0.0, 5000.0, 10000.0, 20000.0]) {
                    if t >= r.stop_time {
                        assert_eq!(row, last, "{method}");
                    }
                }
            }
        }
        assert!(
            result
                .replicates
                .iter()
                .any(|r| r.termination == Termination::Stopped),
            "{method}: no takeover in the fixture"
        );
    }
}

/// Replays fixed words, then repeats the last one.
struct Scripted(Vec<u64>, usize);
impl rand_core::RngCore for Scripted {
    fn next_u32(&mut self) -> u32 {
        self.next_u64() as u32
    }
    fn next_u64(&mut self) -> u64 {
        let w = self.0[self.1.min(self.0.len() - 1)];
        self.1 += 1;
        w
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        rand_core::impls::fill_bytes_via_next(self, dest);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

#[test]
fn a_zero_uniform_does_not_fail_the_static_wait() {
    // Word 0 maps to u = 0, which gave -log(1 - u) = 0 and a failed path.
    let mut rng = Scripted(vec![0, 1 << 63, 1 << 62], 0);
    let mut words = engine::Words::for_test(&mut rng, 16);
    let drawn = flat::direct(&[1.0, 2.0], 0.0, &mut words, &mut engine::Clock::strict());
    let Ok(Some((next, _))) = drawn else {
        panic!("zero uniform failed the draw");
    };
    assert!(next > 0.0 && next.is_finite());
}

#[test]
fn static_resolution_failures_use_the_resolution_code() {
    let mut rng = Scripted(vec![1 << 63], 0);
    let mut words = engine::Words::for_test(&mut rng, 16);
    let Err(engine::Halt::Numerical(d)) =
        flat::direct(&[1e9, 1e-9], 0.0, &mut words, &mut engine::Clock::strict())
    else {
        panic!("tiny channel accepted");
    };
    assert_eq!(d.code, "CRN-R031");
}

#[test]
fn replicates_that_fail_before_tracking_carry_no_partial_observation() {
    let mut model = three_class();
    model["initial"]
        .as_array_mut()
        .unwrap()
        .push(json!({"type": "mtDNA", "attributes": {"clone": 1, "s": "sw", "wt": 0}, "count": 1}));
    let bytes = serde_json::to_vec(&model).unwrap();
    for method in ["tree", "static_direct"] {
        let mut req = request_v(V02, &bytes, 3, 10.0, &[0.0, 10.0]);
        req["method"] = json!({"kind": method});
        req["budget"]["max_instances"] = json!(1);
        req["trackers"] = json!([{"observable": "N"}, {"observable": "mutants"}]);
        let result = solve_values(&model, &req);
        assert!(
            result
                .replicates
                .iter()
                .all(|r| r.termination == Termination::InstanceLimit),
            "{method}"
        );
        assert!(
            result.observations.iter().all(|o| o.trackers.len() == 2),
            "{method}: partial tracker list"
        );
    }
}

const ATP: &str = include_str!("../../../../models/rules/mtdna-three-class-atp.rules.json");

fn atp_run(edit: &[(&str, f64)], replicates: u64, end: f64) -> RulesResult {
    let mut model: Value = serde_json::from_str(ATP).unwrap();
    for (key, value) in edit {
        model["atp"][*key] = json!(*value);
    }
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request_v(V02, &bytes, replicates, end, &[0.0, end / 2.0, end]);
    req["method"] = json!({"kind": "thinning"});
    solve_values(&model, &req)
}

#[test]
fn atp_starting_at_zero_with_positive_forcing_rises_instead_of_stopping() {
    // 1000 wild types: f * wt - v2 * N = 0.1 * 1000 - 0.01 * 1000 = 90 > 0.
    let result = atp_run(&[("initial", 0.0)], 4, 10.0);
    for r in &result.replicates {
        assert_eq!(
            r.termination,
            Termination::Completed,
            "stopped at {}",
            r.stop_time
        );
        assert!(r.events > 0);
    }
    let atp = result.observables.len() - 1;
    assert!(result.paths[0].values[1][atp] > 0.0);
}

#[test]
fn atp_zero_time_survives_an_overflowing_ratio() {
    let flow = atp::Flow {
        atp: 1e300,
        source: -1e-300,
        decay: 1.0,
    };
    // ln(1 + v1 * A / -s) / v1 with the ratio far beyond f64::MAX.
    let expected = 1e300f64.ln() + 1.0f64.ln() - 1e-300f64.ln();
    assert!((flow.zero_time() - expected).abs() < 1e-9 * expected);
    let rising = atp::Flow {
        atp: 0.0,
        source: 1.0,
        decay: 0.5,
    };
    assert_eq!(rising.zero_time(), f64::INFINITY);
}

#[test]
fn atp_depletion_far_from_time_zero_reaches_the_zero_crossing() {
    // Falling ATP with a tiny scale: replication is off, N decays, and the
    // maintenance use 1e-5 * N empties A0 = 3 near t = 700. The
    // half-denominator cap then shrinks below the clock resolution.
    let edit = [
        ("initial", 3.0),
        ("production", 0.0),
        ("decay", 0.0),
        ("maintenance", 1e-5),
        ("scale", 1e-20),
    ];
    let result = atp_run(&edit, 4, 5000.0);
    for r in &result.replicates {
        assert_eq!(
            r.termination,
            Termination::Stopped,
            "{:?}",
            r.diagnostic.as_ref().map(|d| &d.message)
        );
        assert!(r.stop_time > 300.0);
    }
}

fn leap_run(model: &Value, method: Value, replicates: u64, end: f64, times: &[f64]) -> RulesResult {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut req = request_v(V02, &bytes, replicates, end, times);
    req["method"] = method;
    solve_values(model, &req)
}

#[test]
fn tau_leaps_respect_the_vacancy_capacity() {
    // Exact SSA never exceeds K; species-only step control let leaps overfill
    // the vacancy because the rate depends on K - N, not on each count.
    let model = three_class();
    let times: Vec<f64> = (0..=100).map(|i| f64::from(i) * 10.0).collect();
    let result = leap_run(
        &model,
        json!({"kind": "tau", "epsilon": 0.01}),
        8,
        1000.0,
        &times,
    );
    let n = result.observables.iter().position(|o| o == "N").unwrap();
    let worst = result
        .paths
        .iter()
        .flat_map(|p| p.values.iter().map(|row| row[n]))
        .fold(0.0, f64::max);
    assert!(worst <= 1000.0, "tau reached N = {worst} above K = 1000");
}

#[test]
fn binomial_steps_use_competing_risk_probabilities() {
    // Pure degradation at rate 1 over one unit step: E[n] = 1000 / e.
    let model = with_parameters(&[("delta", 1.0), ("sw", 0.0), ("sm", 0.0)]);
    let result = leap_run(
        &model,
        json!({"kind": "binomial", "step": 1.0}),
        64,
        1.0,
        &[0.0, 1.0],
    );
    let n = result.observables.iter().position(|o| o == "N").unwrap();
    let mean = result.paths.iter().map(|p| p.values[1][n]).sum::<f64>() / 64.0;
    let exact = 1000.0 * (-1.0f64).exp();
    assert!((mean - exact).abs() < 10.0, "mean {mean}, expected {exact}");
}

#[test]
fn binomial_steps_do_not_depend_on_the_report_grid() {
    let model = with_parameters(&[("mu", 3e-3)]);
    let method = json!({"kind": "binomial", "step": 1.0});
    let aligned = leap_run(&model, method.clone(), 6, 4.0, &[0.0, 1.0, 2.0, 3.0, 4.0]);
    let offset = leap_run(&model, method, 6, 4.0, &[0.0, 0.5, 1.5, 2.5, 4.0]);
    for (a, b) in aligned.paths.iter().zip(&offset.paths) {
        assert_eq!(a.values[4], b.values[4], "final state depends on the grid");
        // Reports between lattice points hold the last step's state.
        assert_eq!(b.values[1], a.values[0]);
        assert_eq!(b.values[2], a.values[1]);
        assert_eq!(b.values[3], a.values[2]);
    }
}

#[test]
fn the_cancel_flag_check_does_not_read_the_cancel_file() {
    let dir = std::env::temp_dir().join(format!("crn-flag-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("cancel");
    std::fs::write(&file, b"").unwrap();
    let control = Control {
        cancel_file: Some(file.clone()),
        ..Control::default()
    };
    assert!(!control.cancel_flagged());
    assert!(control.is_cancelled());
    assert!(control.cancel_flagged());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_disabled_observer_is_not_called_and_report_only_values_stay_fresh() {
    // An observer that declares ENABLED = false was still called whenever
    // trackers existed, and could stop a run whose report-only observables
    // had not been refreshed.
    struct QuietStop;
    impl observe::EventObserver for QuietStop {
        const ENABLED: bool = false;
        fn on_event(&self, _: &observe::Event<'_>) -> observe::ObserverAction {
            observe::ObserverAction::Stop
        }
    }
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "q",
        "types": [{"name": "X"}], "initial": [],
        "observables": {"X": {"type": "X"}, "W": {"type": "X", "weight": "2.5"}},
        "rules": [{"name": "birth", "reactants": [], "rate": "5", "products": [{"type": "X"}]}]
    });
    let bytes = serde_json::to_vec(&model).unwrap();
    for method in ["tree", "static_direct", "next_reaction"] {
        let mut req = request_v(V02, &bytes, 2, 2.0, &[0.0, 1.0, 2.0]);
        req["method"] = json!({"kind": method});
        req["trackers"] = json!([{"observable": "X"}]);
        let prepared = prepare(&bytes, &serde_json::to_vec(&req).unwrap()).unwrap();
        let plain = execute(&prepared, &Control::default()).unwrap();
        let quiet = execute_with_observer(&prepared, &Control::default(), &QuietStop).unwrap();
        assert!(
            serde_json::to_vec(&plain).unwrap() == serde_json::to_vec(&quiet).unwrap(),
            "{method}: a disabled observer changed the run"
        );
        let last = quiet.paths[0].values.last().unwrap();
        assert_eq!(last[1] * 2.5, last[0], "{method}: held state {last:?}");
    }
}

#[test]
fn sparse_trackers_without_an_eligible_tracker_keep_the_default_identity() {
    let model = three_class();
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request_v("crn-craft.rules-run/0.4", &bytes, 4, 50.0, &[0.0, 50.0]);
    req["method"] = json!({"kind": "static_direct", "sparse_trackers": true});
    // A first-passage predicate makes the only tracker ineligible.
    req["trackers"] = json!([{"observable": "mutants", "first_passage": "mutants >= 2"}]);
    let result = solve_values(&model, &req);
    assert_eq!(result.sampler, "chacha8-u53/rules-static-direct/1");
    let warning = result
        .diagnostics
        .iter()
        .find(|d| d.code == "CRN-R051")
        .expect("ineligible request not reported");
    assert!(
        warning.message.contains("sparse_trackers"),
        "{}",
        warning.message
    );
}
