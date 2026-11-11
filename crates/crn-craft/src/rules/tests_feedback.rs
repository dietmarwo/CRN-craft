//! Regression tests from the first user feedback on the 0.2.0 interface.
use super::*;
use serde_json::{Value, json};

const THREE_CLASS: &str = include_str!("../../../../models/rules/mtdna-three-class.rules.json");

/// The default `max_output_bytes` of the Python helper `rules.ensemble`.
const PYTHON_DEFAULT_OUTPUT: u64 = 64 << 20;

/// The request `crn_craft.rules.ensemble` builds for `model` with an event
/// log, as in the user manual's chapter on logs: 25,920 hours, two report
/// times, two threads and the Python default budgets.
fn python_request(model: &[u8], method: &str, replicates: u64, log: &Value) -> Value {
    json!({
        "schema": if method == "mtdna_direct" {
            "crn-craft.rules-run/0.3"
        } else {
            "crn-craft.rules-run/0.2"
        },
        "model": {"path": "model.rules.json", "sha256": crate::sha256(model)},
        "time": {"end": 25_920.0, "report_times": [0.0, 25_920.0]},
        "ensemble": {"replicates": replicates, "root_seed": 4},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"threads": 2, "max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 8_000_000u64,
                   "max_instances": 100_000u64, "max_output_bytes": PYTHON_DEFAULT_OUTPUT},
        "method": {"kind": method},
        "event_log": log
    })
}

fn prepare_json(model: &[u8], request: &Value) -> Result<Prepared, Vec<Diagnostic>> {
    prepare(model, &serde_json::to_vec(request).unwrap())
}

/// The reported request: the three-class model with a bounded log of the two
/// mutation rules, at the replicate counts the user manual documents (the
/// log chapter's 100, the first ensemble's 200, the sweep's 400 and the
/// tracker chapter's 500), with the Python default output budget.
#[test]
fn reported_event_log_fits_the_python_default_budget() {
    let model = THREE_CLASS.as_bytes();
    let log = json!({"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]});
    for method in ["static_direct", "mtdna_direct"] {
        for replicates in [100, 200, 400, 500] {
            let request = python_request(model, method, replicates, &log);
            if let Err(diagnostics) = prepare_json(model, &request) {
                panic!("{method} with {replicates} replicates: {diagnostics:?}");
            }
        }
    }
}

/// The reservation's bound and the measured compact JSON of the retained log
/// of every replicate, for logs filled to `max_entries`.
fn log_bound_and_measured(model: &[u8], method: Value, entries: usize) -> (u128, Vec<usize>) {
    let mut request = python_request(model, "tree", 3, &json!({"max_entries": entries}));
    request["method"] = method;
    request["time"] = json!({"end": 2_000.0, "report_times": [0.0, 2_000.0]});
    let prepared = prepare_json(model, &request).unwrap();
    let sizes = observation_sizes(
        &prepared.compiled,
        prepared.request.method,
        prepared.flat.as_ref(),
    );
    let reserved = prepared.observations.reservation(&sizes);
    let result = execute(&prepared, &Control::default()).unwrap();
    assert_eq!(result.termination, Termination::Completed);
    let measured = result
        .observations
        .iter()
        .map(|o| {
            assert_eq!(o.log.len(), entries, "log not filled");
            serde_json::to_vec(&o.log).unwrap().len()
        })
        .collect();
    (reserved.entries * reserved.entry, measured)
}

/// The tighter bound still covers what a full log holds: exact steps of the
/// three-class model, leaps, and the clone model on `tree`, whose new
/// identities and drawn weights are attributes, plus a type whose drawn
/// attributes print at full precision.
#[test]
fn reservation_covers_a_full_log() {
    let three = THREE_CLASS.as_bytes();
    let clones = include_bytes!("../../../../models/rules/mtdna-fixed.rules.json");
    let noisy = serde_json::to_vec(&json!({
        "schema": "crn-craft.rules-model/0.1", "name": "noisy",
        "parameters": {"k": 1.0},
        "types": [{"name": "A", "attributes": ["u", "v", "w"]}],
        "initial": [{"type": "A", "attributes": {"u": 0.1, "v": -0.3, "w": 1e-300}, "count": 1}],
        "observables": {"A": {"type": "A"}, "spread": {"type": "A", "weight": "x.u * 3.3"}},
        "rules": [{"name": "redraw", "reactants": [{"label": "a", "type": "A"}], "rate": "k",
                   "products": [{"type": "A", "attributes": {
                       "u": "normal(0, 1) * 1e-7", "v": "-uniform(0, 1)", "w": "normal(0, 1) * 1e300"}}]}]
    }))
    .unwrap();
    let cases: [(&[u8], Value); 6] = [
        (three, json!({"kind": "static_direct"})),
        (three, json!({"kind": "next_reaction"})),
        (three, json!({"kind": "tree"})),
        (three, json!({"kind": "tau", "epsilon": 0.03})),
        (clones, json!({"kind": "tree"})),
        (&noisy, json!({"kind": "tree"})),
    ];
    for (model, method) in cases {
        let (bound, measured) = log_bound_and_measured(model, method.clone(), 50);
        for bytes in measured {
            assert!(
                bytes as u128 <= bound,
                "{method}: {bytes} log bytes above the bound {bound}"
            );
        }
    }
}

/// The diagnostic itemizes the reservation and names the budget that fits.
#[test]
fn reservation_error_names_bytes_and_fix() {
    let model = THREE_CLASS.as_bytes();
    let log = json!({"max_entries": 50, "rules": ["mutation_mt"]});
    let mut request = python_request(model, "static_direct", 20_000, &log);
    let diagnostics = prepare_json(model, &request).err().unwrap();
    let error = &diagnostics[0];
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-R020", "/budget")
    );
    // Three participants of three attributes and six observables:
    // 160 + 6 * 25 + 3 * (62 + 3 * 25) = 721 bytes per entry.
    let total = 20_000u128 * 50 * 721;
    assert!(
        error
            .message
            .contains("20000 replicates x 36050 bytes (50 log entries x 721 bytes)")
            && error.message.contains(&format!("= {total} bytes"))
            && error
                .message
                .contains(&format!("max_output_bytes = {}", 64u64 << 20)),
        "{}",
        error.message
    );
    assert!(
        error.hint.contains(&format!("at least {total}")),
        "{}",
        error.hint
    );
    // That budget is enough.
    request["budget"]["max_output_bytes"] = json!(u64::try_from(total).unwrap());
    assert!(prepare_json(model, &request).is_ok());
    request["budget"]["max_output_bytes"] = json!(u64::try_from(total).unwrap() - 1);
    assert!(prepare_json(model, &request).is_err());
}

// ---------------------------------------------------------------------------
// I4: whole-rule propensities (rules-model/0.3)

const V03: &str = "crn-craft.rules-model/0.3";

/// Attribute-free species `A` with birth and death, written per reactant
/// tuple or as whole-rule propensities.
fn birth_death(whole: bool, initial: u64) -> Value {
    let rate = |per_tuple: &str, total: &str| {
        if whole {
            json!({"propensity": total})
        } else {
            json!({"rate": per_tuple})
        }
    };
    let rule = |name: &str, reactants: Value, rate: Value, products: Value| {
        let mut rule = json!({"name": name, "reactants": reactants, "products": products});
        rule.as_object_mut()
            .unwrap()
            .extend(rate.as_object().unwrap().clone());
        rule
    };
    let a = json!([{"label": "a", "type": "A"}]);
    json!({
        "schema": V03, "name": "birth_death",
        "parameters": {"b": 0.9, "d": 0.3},
        "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": initial}],
        "observables": {"A": {"type": "A"}},
        "rules": [
            rule("birth", a.clone(), rate("b", "b * A"), json!([{"type": "A"}, {"type": "A"}])),
            rule("death", a, rate("d", "d * A"), json!([])),
        ]
    })
}

fn run_request(model: &[u8], method: Value, replicates: u64) -> Value {
    // Track the first observable of the model.
    let parsed: Value = serde_json::from_slice(model).unwrap();
    let tracked = parsed["observables"]
        .as_object()
        .unwrap()
        .keys()
        .next()
        .unwrap()
        .clone();
    json!({
        "schema": "crn-craft.rules-run/0.4",
        "model": {"path": "m.rules.json", "sha256": crate::sha256(model)},
        "time": {"end": 4.0, "report_times": [0.0, 1.0, 2.0, 4.0]},
        "ensemble": {"replicates": replicates, "root_seed": 11},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates.min(3),
                    "keep_instances": replicates.min(2)},
        "budget": {"max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 10_000_000u64,
                   "max_instances": 1000u64, "threads": 2, "max_output_bytes": 1u64 << 28},
        "method": method,
        "trackers": [{"observable": tracked, "keep_episodes": 2}],
        "event_log": {"max_entries": 5}
    })
}

fn run_model(
    model: &Value,
    method: Value,
    replicates: u64,
) -> Result<RulesResult, Vec<Diagnostic>> {
    let bytes = serde_json::to_vec(model).unwrap();
    let request = run_request(&bytes, method, replicates);
    let prepared = prepare_json(&bytes, &request)?;
    Ok(execute(&prepared, &Control::default()).unwrap())
}

/// The numbers of a result, without the hashes of its documents.
fn numbers(result: &RulesResult) -> Value {
    let mut value = serde_json::to_value(result).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("model_sha256");
    object.remove("run_sha256");
    value
}

const METHODS: [&str; 6] = [
    r#"{"kind": "tree"}"#,
    r#"{"kind": "static_direct"}"#,
    r#"{"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true}"#,
    r#"{"kind": "next_reaction"}"#,
    r#"{"kind": "tau", "epsilon": 0.03}"#,
    r#"{"kind": "binomial", "step": 0.05}"#,
];

/// For one reactant, `propensity: "d * A"` and `rate: "d"` give the same
/// floating-point propensity, so every method reproduces the mass-action
/// run exactly, path by path.
#[test]
fn every_method_runs_whole_rule_propensities() {
    for method in METHODS {
        let method: Value = serde_json::from_str(method).unwrap();
        let tuple = run_model(&birth_death(false, 20), method.clone(), 64).unwrap();
        let whole = run_model(&birth_death(true, 20), method.clone(), 64).unwrap();
        assert!(numbers(&tuple) == numbers(&whole), "{method}");
        assert!(whole.events > 0 && whole.termination == Termination::Completed);
    }
}

/// The certified three-class kernel does not apply to propensity rules; the
/// run falls back to Static Direct and says why.
#[test]
fn mtdna_direct_falls_back_for_propensities() {
    let model = json!({
        "schema": V03, "name": "three_species",
        "parameters": {"K": 50, "delta": 0.05, "mu": 0.01},
        "types": [{"name": "W"}, {"name": "M"}],
        "initial": [{"type": "W", "count": 50}],
        "observables": {"W": {"type": "W"}, "M": {"type": "M"}},
        "rules": [
            {"name": "death_w", "reactants": [{"label": "a", "type": "W"}],
             "propensity": "delta * W", "products": []},
            {"name": "birth_w", "reactants": [{"label": "a", "type": "W"}],
             "propensity": "max(0, K - W - M) / (W + M) * (1 - mu) * W",
             "products": [{"type": "W"}, {"type": "W"}]},
            {"name": "mutation", "reactants": [{"label": "a", "type": "W"}],
             "propensity": "max(0, K - W - M) / (W + M) * mu * W",
             "products": [{"type": "W"}, {"type": "M"}]},
            {"name": "death_m", "reactants": [{"label": "a", "type": "M"}],
             "rate": "delta", "products": []},
            {"name": "birth_m", "reactants": [{"label": "a", "type": "M"}],
             "rate": "max(0, K - W - M) / (W + M)", "products": [{"type": "M"}, {"type": "M"}]}
        ]
    });
    let result = run_model(&model, json!({"kind": "mtdna_direct"}), 8).unwrap();
    assert_eq!(result.sampler, "chacha8-u53/rules-static-direct/1");
    let warning = &result.diagnostics[0];
    assert_eq!(warning.code, "CRN-R051");
    assert!(
        warning.message.contains("whole-rule propensities"),
        "{}",
        warning.message
    );
}

/// A whole-rule propensity applies only while each reactant pattern can take
/// distinct molecules: two of `A` for `A + A`, one of each for `A + B`.
#[test]
fn whole_rule_propensities_need_enough_molecules() {
    let model = |a: u64, b: u64| {
        json!({
            "schema": V03, "name": "pairs",
            "parameters": {"c": 5.0},
            "types": [{"name": "A"}, {"name": "B"}],
            "initial": [{"type": "A", "count": a}, {"type": "B", "count": b}],
            "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
            "rules": [
                {"name": "pair", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "A"}],
                 "propensity": "c", "products": []},
                {"name": "cross", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "B"}],
                 "propensity": "c", "products": [{"type": "A"}]}
            ]
        })
    };
    for method in [
        json!({"kind": "tree"}),
        json!({"kind": "static_direct"}),
        json!({"kind": "next_reaction"}),
    ] {
        let events = |a, b| {
            let result = run_model(&model(a, b), method.clone(), 16).unwrap();
            result
                .replicates
                .iter()
                .map(|r| r.rule_events.clone())
                .collect::<Vec<_>>()
        };
        // One A, no B: nothing can fire.
        assert!(events(1, 0).iter().all(|e| e == &[0, 0]), "{method}");
        // Three A: one pair fires and leaves a single A.
        assert!(events(3, 0).iter().all(|e| e == &[1, 0]), "{method}");
        // One A and one B: the cross rule fires (B is consumed), then
        // neither can fire again.
        assert!(events(1, 1).iter().all(|e| e == &[0, 1]), "{method}");
    }
}

/// On `tree`, a firing picks distinct molecules uniformly, whatever their
/// attributes; the static methods refuse such types.
#[test]
fn whole_rule_firings_pick_molecules_uniformly() {
    let model = |pair: bool| {
        let reactants = if pair {
            json!([{"label": "x", "type": "X"}, {"label": "y", "type": "X"}])
        } else {
            json!([{"label": "x", "type": "X"}])
        };
        json!({
            "schema": V03, "name": "pick",
            "parameters": {},
            "types": [{"name": "X", "attributes": ["v"]}],
            "initial": [{"type": "X", "attributes": {"v": 1}, "count": 1},
                        {"type": "X", "attributes": {"v": 2}, "count": 3}],
            "observables": {"X": {"type": "X"}},
            "rules": [{"name": "take", "reactants": reactants, "propensity": "1", "products": []}]
        })
    };
    let replicates = 4000u64;
    for (pair, expected) in [(false, 0.25), (true, 0.5)] {
        let bytes = serde_json::to_vec(&model(pair)).unwrap();
        let mut request = run_request(&bytes, json!({"kind": "tree"}), replicates);
        request["time"] = json!({"end": 100.0, "report_times": [100.0]});
        request["event_log"] = json!({"max_entries": 1});
        request["trackers"] = json!([]);
        let result = execute(
            &prepare_json(&bytes, &request).unwrap(),
            &Control::default(),
        )
        .unwrap();
        // The share of first firings that took the single v = 1 molecule:
        // 1 of 4 molecules, or 3 of the 6 pairs.
        let hits = result
            .observations
            .iter()
            .filter(|o| o.log[0].reactants.iter().any(|r| r.attributes == [1.0]))
            .count();
        #[allow(clippy::cast_precision_loss)]
        let share = hits as f64 / replicates as f64;
        let sd = (expected * (1.0 - expected) / 4000.0f64).sqrt();
        assert!((share - expected).abs() < 5.0 * sd, "pair={pair}: {share}");
        let refused = prepare_json(
            &bytes,
            &run_request(&bytes, json!({"kind": "static_direct"}), 2),
        )
        .err()
        .unwrap();
        assert_eq!(
            (refused[0].code.as_str(), refused[0].path.as_str()),
            ("CRN-R020", "/method")
        );
        assert!(refused[0].message.contains("several instance tuples"));
    }
}

/// Contract errors of the new field, with stable codes and paths.
#[test]
fn propensity_contract_errors() {
    let check = |model: &Value| {
        let report = check_model(&serde_json::to_vec(model).unwrap());
        report
            .diagnostics
            .first()
            .map(|d| (d.code.clone(), d.path.clone(), d.message.clone()))
    };
    let base = birth_death(true, 5);
    let edit = |f: &dyn Fn(&mut Value)| {
        let mut m = base.clone();
        f(&mut m);
        check(&m)
    };
    assert_eq!(check(&base), None);
    let older = edit(&|m| m["schema"] = json!("crn-craft.rules-model/0.2")).unwrap();
    assert_eq!(
        (older.0.as_str(), older.1.as_str()),
        ("CRN-R015", "/rules/0/propensity")
    );
    let missing = edit(&|m| {
        m["schema"] = json!("crn-craft.rules-model/0.1");
        m["rules"][0].as_object_mut().unwrap().remove("propensity");
    })
    .unwrap();
    assert_eq!(
        missing,
        (
            "CRN-E003".into(),
            "/rules/0".into(),
            "missing field `rate`".into()
        )
    );
    let both = edit(&|m| m["rules"][0]["rate"] = json!("b")).unwrap();
    assert_eq!((both.0.as_str(), both.1.as_str()), ("CRN-R013", "/rules/0"));
    let neither = edit(&|m| {
        m["rules"][0].as_object_mut().unwrap().remove("propensity");
    })
    .unwrap();
    assert_eq!(
        (neither.0.as_str(), neither.1.as_str()),
        ("CRN-R013", "/rules/0")
    );
    for key in ["where", "factor"] {
        let filtered = edit(&|m| m["rules"][1]["reactants"][0][key] = json!("1")).unwrap();
        assert_eq!(
            (filtered.0.as_str(), filtered.1.clone()),
            ("CRN-R013", format!("/rules/1/reactants/0/{key}"))
        );
    }
    let label = edit(&|m| m["rules"][1]["propensity"] = json!("d * a.size")).unwrap();
    assert_eq!(
        (label.0.as_str(), label.1.as_str()),
        ("CRN-R012", "/rules/1/propensity")
    );
    let null = edit(&|m| m["rules"][1]["propensity"] = Value::Null).unwrap();
    assert_eq!(null.0, "CRN-E003");
}

/// A negative or undefined propensity fails the replicate with its pointer.
#[test]
fn invalid_propensities_fail_the_replicate() {
    let mut model = birth_death(true, 5);
    model["rules"][1]["propensity"] = json!("d * A - 2");
    for method in [json!({"kind": "tree"}), json!({"kind": "static_direct"})] {
        let result = run_model(&model, method.clone(), 2).unwrap();
        let record = &result.replicates[0];
        assert_eq!(record.termination, Termination::Failed, "{method}");
        let diagnostic = record.diagnostic.as_ref().unwrap();
        assert_eq!(
            (diagnostic.code.as_str(), diagnostic.path.as_str()),
            ("CRN-R030", "/rules/1/propensity"),
            "{method}"
        );
    }
}

/// Generated rates kernels refuse models with whole-rule propensities;
/// generated event loops emit them (`tests/generated_loop.rs` replays them).
#[test]
fn generated_kernels_refuse_propensities_and_loops_emit_them() {
    let bytes = serde_json::to_vec(&birth_death(true, 5)).unwrap();
    let mut request = run_request(&bytes, json!({"kind": "static_direct"}), 2);
    request["schema"] = json!("crn-craft.rules-run/0.2");
    let prepared = prepare_json(&bytes, &request).unwrap();
    let error = generate_rules_kernel(&prepared).err().unwrap();
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-R020", "/method")
    );
    assert!(
        error.message.contains("whole-rule propensities"),
        "{}",
        error.message
    );
    let source = generate_rules_loop(&prepared).unwrap();
    assert!(source.contains("f64::from(u8::from(counts["), "{source}");
}

/// rules-model/0.3 keeps the ATP block, and thinning samples propensity
/// rules with it.
#[test]
fn thinning_runs_propensities_in_atp_models() {
    let mut model: Value = serde_json::from_str(include_str!(
        "../../../../models/rules/mtdna-three-class-atp.rules.json"
    ))
    .unwrap();
    model["schema"] = json!(V03);
    model["types"]
        .as_array_mut()
        .unwrap()
        .push(json!({"name": "Marker"}));
    model["observables"]["markers"] = json!({"type": "Marker"});
    model["rules"].as_array_mut().unwrap().push(json!({
        "name": "mark", "reactants": [], "propensity": "0.5", "products": [{"type": "Marker"}]
    }));
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut request = run_request(&bytes, json!({"kind": "thinning"}), 8);
    request["time"] = json!({"end": 200.0, "report_times": [0.0, 200.0]});
    request["trackers"] = json!([]);
    let result = execute(
        &prepare_json(&bytes, &request).unwrap(),
        &Control::default(),
    )
    .unwrap();
    assert_eq!(result.sampler, "chacha8-u53/rules-atp-thinning/1");
    let mean = result.statistics.unwrap().mean.unwrap();
    let markers = result
        .observables
        .iter()
        .position(|o| o == "markers")
        .unwrap();
    // Poisson with mean 0.5 * 200 = 100 per cell.
    assert!(
        (mean[1][markers] - 100.0).abs() < 15.0,
        "{}",
        mean[1][markers]
    );
}

// ---------------------------------------------------------------------------
// I5: continuous variables and the hybrid method

const ATP_LEGACY: &str = include_str!("../../../../models/rules/mtdna-three-class-atp.rules.json");
const ATP_CONTINUOUS: &str =
    include_str!("../../../../models/rules/mtdna-three-class-atp-continuous.rules.json");

/// One variable `x` with `derivative`, and molecules `A` that die at rate
/// `death` per molecule (an expression that may read `x`).
fn continuous_model(derivative: &str, x0: f64, bound: Option<&str>, death: &str, a0: u64) -> Value {
    let mut variable = json!({"name": "x", "initial": x0, "derivative": derivative});
    if let Some(action) = bound {
        variable["lower_bound"] = json!(0);
        variable["on_bound"] = json!(action);
    }
    json!({
        "schema": V03, "name": "continuous",
        "parameters": {"k": 1.0},
        "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": a0}],
        "observables": {"A": {"type": "A"}},
        "rules": [{"name": "death", "reactants": [{"label": "a", "type": "A"}],
                   "rate": death, "products": []}],
        "continuous": [variable]
    })
}

fn hybrid_request(model: &[u8], replicates: u64, end: f64, times: &[f64]) -> Value {
    json!({
        "schema": "crn-craft.rules-run/0.5",
        "model": {"path": "m.rules.json", "sha256": crate::sha256(model)},
        "time": {"end": end, "report_times": times},
        "ensemble": {"replicates": replicates, "root_seed": 5},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 10_000_000u64,
                   "max_instances": 100u64, "threads": 2, "max_output_bytes": 1u64 << 28},
        "method": {"kind": "hybrid", "rtol": 1e-8, "atol": 1e-10}
    })
}

fn hybrid_run(model: &Value, replicates: u64, end: f64, times: &[f64]) -> RulesResult {
    let bytes = serde_json::to_vec(model).unwrap();
    let request = hybrid_request(&bytes, replicates, end, times);
    execute(
        &prepare_json(&bytes, &request).unwrap(),
        &Control::default(),
    )
    .unwrap()
}

fn first_error(model: &Value, request: &Value) -> (String, String, String) {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut request = request.clone();
    request["model"]["sha256"] = json!(crate::sha256(&bytes));
    let d = prepare_json(&bytes, &request).err().unwrap();
    (d[0].code.clone(), d[0].path.clone(), d[0].message.clone())
}

/// Between events the variables follow the ODE: reports read the dense
/// output, and the result records the method and the continuous columns.
#[test]
fn hybrid_integrates_between_events() {
    // No molecule, so no event: x' = -k x and a logistic second variable.
    let mut model = continuous_model("-k * x", 2.0, None, "k", 0);
    model["continuous"].as_array_mut().unwrap().push(json!({
        "name": "z", "initial": 0.1, "derivative": "z * (1 - z)"
    }));
    let times: Vec<f64> = (0..=10).map(f64::from).collect();
    let result = hybrid_run(&model, 2, 10.0, &times);
    assert_eq!(result.schema, ResultVersion::V03);
    assert_eq!(
        result.sampler,
        "chacha8-u53/rules-approximate-hybrid-dopri54/1"
    );
    assert_eq!(result.observables, ["A", "x", "z"]);
    assert_eq!(result.continuous, ["x", "z"]);
    assert_eq!(
        result.method,
        Some(Method::Hybrid {
            rtol: 1e-8,
            atol: 1e-10
        })
    );
    for (row, &t) in result.paths[0].values.iter().zip(&times) {
        let x = 2.0 * libm::exp(-t);
        let z = 1.0 / (1.0 + 9.0 * libm::exp(-t));
        assert!((row[1] - x).abs() <= 1e-6 * x + 1e-9, "x({t}) = {}", row[1]);
        assert!((row[2] - z).abs() <= 1e-6 * z + 1e-9, "z({t}) = {}", row[2]);
    }
    let json = serde_json::to_value(&result).unwrap();
    assert_eq!(
        json["method"],
        json!({"kind": "hybrid", "rtol": 1e-8, "atol": 1e-10})
    );
}

/// With `on_bound: "stop"` the replicate stops where the variable reaches
/// the bound and reports the bound afterwards; with `"clamp"` it holds the
/// variable there and completes.
#[test]
fn hybrid_lower_bounds_stop_or_clamp() {
    let times = [0.0, 1.0, 2.5, 4.0];
    let stop = hybrid_run(
        &continuous_model("-k", 2.0, Some("stop"), "0.1", 3),
        4,
        4.0,
        &times,
    );
    for record in &stop.replicates {
        assert_eq!(record.termination, Termination::Stopped);
        assert!(
            (record.stop_time - 2.0).abs() < 1e-12,
            "{}",
            record.stop_time
        );
    }
    let x = |result: &RulesResult, row: usize| result.paths[0].values[row][1];
    assert!((x(&stop, 1) - 1.0).abs() < 1e-12);
    assert_eq!((x(&stop, 2), x(&stop, 3)), (0.0, 0.0));
    // A derivative that turns positive at the bound continues the path.
    let rising = hybrid_run(
        &continuous_model("k - A", 0.0, Some("stop"), "1", 0),
        2,
        4.0,
        &times,
    );
    assert!(
        rising
            .replicates
            .iter()
            .all(|r| r.termination == Termination::Completed)
    );
    assert!((x(&rising, 3) - 4.0).abs() < 1e-9);
    // At the bound with a negative derivative at time 0: stopped at once.
    let flat = hybrid_run(
        &continuous_model("-k", 0.0, Some("stop"), "1", 1),
        2,
        4.0,
        &times,
    );
    assert!(
        flat.replicates
            .iter()
            .all(|r| r.termination == Termination::Stopped && r.stop_time == 0.0)
    );
    let clamp = hybrid_run(
        &continuous_model("-k", 2.0, Some("clamp"), "0.1", 3),
        4,
        4.0,
        &times,
    );
    assert!(
        clamp
            .replicates
            .iter()
            .all(|r| r.termination == Termination::Completed)
    );
    assert!((x(&clamp, 1) - 1.0).abs() < 1e-12);
    assert_eq!((x(&clamp, 2), x(&clamp, 3)), (0.0, 0.0));
}

/// The ensemble does not depend on the thread count, and the legacy ATP
/// results keep their version and shape.
#[test]
fn hybrid_replays_and_legacy_results_keep_their_shape() {
    let model: Value = serde_json::from_str(ATP_CONTINUOUS).unwrap();
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut request = hybrid_request(&bytes, 6, 500.0, &[0.0, 250.0, 500.0]);
    request["method"] = json!({"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9});
    let run = |threads: u32| {
        let mut request = request.clone();
        request["budget"]["threads"] = json!(threads);
        let result = execute(
            &prepare_json(&bytes, &request).unwrap(),
            &Control::default(),
        )
        .unwrap();
        numbers(&result)
    };
    assert!(run(1) == run(2), "thread-dependent hybrid result");
    // A rules-run/0.4 thinning result carries neither `method` nor `continuous`.
    let legacy = ATP_LEGACY.as_bytes();
    let mut old = run_request(legacy, json!({"kind": "thinning"}), 2);
    old["time"] = json!({"end": 10.0, "report_times": [0.0, 10.0]});
    let result = execute(&prepare_json(legacy, &old).unwrap(), &Control::default()).unwrap();
    let json = serde_json::to_value(&result).unwrap();
    assert_eq!(json["schema"], "crn-craft.rules-result/0.2");
    assert!(json.get("method").is_none() && json.get("continuous").is_none());
    // The same thinning request as rules-run/0.5 records both.
    old["schema"] = json!("crn-craft.rules-run/0.5");
    let result = execute(&prepare_json(legacy, &old).unwrap(), &Control::default()).unwrap();
    assert_eq!(result.schema, ResultVersion::V03);
    assert_eq!(result.continuous, ["ATP"]);
    assert_eq!(result.method, Some(Method::Thinning {}));
}

/// Model contract errors for continuous variables.
#[test]
fn continuous_contract_errors() {
    let base = continuous_model("-k * x", 1.0, Some("stop"), "k * x", 1);
    let check = |edit: &dyn Fn(&mut Value)| {
        let mut m = base.clone();
        edit(&mut m);
        let report = check_model(&serde_json::to_vec(&m).unwrap());
        report
            .diagnostics
            .first()
            .map(|d| (d.code.clone(), d.path.clone()))
    };
    let pair = |code: &str, path: &str| Some((code.to_string(), path.to_string()));
    assert_eq!(check(&|_| {}), None);
    assert_eq!(
        check(&|m| m["schema"] = json!("crn-craft.rules-model/0.2")),
        pair("CRN-R015", "/continuous")
    );
    assert_eq!(
        check(&|m| m["continuous"][0]["name"] = json!("k")),
        pair("CRN-R010", "/continuous/0/name")
    );
    assert_eq!(
        check(&|m| m["continuous"][0]["name"] = json!("A")),
        pair("CRN-R010", "/continuous/0/name")
    );
    assert_eq!(
        check(&|m| m["continuous"][0]["initial"] = json!(-1)),
        pair("CRN-R014", "/continuous/0/initial")
    );
    assert_eq!(
        check(&|m| {
            m["continuous"][0]
                .as_object_mut()
                .unwrap()
                .remove("on_bound");
        }),
        pair("CRN-R013", "/continuous/0")
    );
    assert_eq!(
        check(&|m| m["continuous"][0]["derivative"] = json!("-y")),
        pair("CRN-R012", "/continuous/0/derivative")
    );
    assert_eq!(
        check(&|m| m["stop"] = json!("x < 0.5")),
        pair("CRN-R012", "/stop")
    );
    assert_eq!(
        check(&|m| m["observables"]["B"] = json!({"type": "A", "weight": "x"})),
        pair("CRN-R012", "/observables/B/weight")
    );
    assert_eq!(
        check(&|m| m["continuous"][0]["on_bound"] = Value::Null).map(|p| p.0),
        Some("CRN-E003".into())
    );
    let mut atp: Value = serde_json::from_str(ATP_LEGACY).unwrap();
    atp["schema"] = json!(V03);
    atp["continuous"] = json!([{"name": "x", "initial": 1, "derivative": "0"}]);
    let report = check_model(&serde_json::to_vec(&atp).unwrap());
    assert_eq!(report.diagnostics[0].code, "CRN-R015");
    assert_eq!(report.diagnostics[0].path, "/continuous");
}

/// Every other method refuses continuous variables, and hybrid refuses what
/// it cannot run, with coded diagnostics.
#[test]
fn methods_refuse_what_they_cannot_run() {
    let model = continuous_model("-k * x", 1.0, None, "k * x", 3);
    let bytes = serde_json::to_vec(&model).unwrap();
    let base = hybrid_request(&bytes, 2, 1.0, &[0.0, 1.0]);
    for method in METHODS
        .iter()
        .map(|m| serde_json::from_str::<Value>(m).unwrap())
        .chain([json!({"kind": "mtdna_direct"}), json!({"kind": "thinning"})])
    {
        let mut request = base.clone();
        request["method"] = method.clone();
        let (code, path, message) = first_error(&model, &request);
        assert_eq!(
            (code.as_str(), path.as_str()),
            ("CRN-R020", "/method"),
            "{method}"
        );
        assert!(
            message.contains("require the hybrid method"),
            "{method}: {message}"
        );
    }
    let mut request = base.clone();
    request["schema"] = json!("crn-craft.rules-run/0.4");
    assert_eq!(first_error(&model, &request).1, "/schema");
    for (key, value, path) in [
        ("rtol", 0.5, "/method/rtol"),
        ("rtol", 1e-13, "/method/rtol"),
        ("atol", 0.0, "/method/atol"),
        ("atol", 0.1, "/method/atol"),
    ] {
        let mut request = base.clone();
        request["method"][key] = json!(value);
        assert_eq!(first_error(&model, &request).1, path, "{key}={value}");
    }
    let mut request = base.clone();
    request["trackers"] = json!([{"observable": "x"}]);
    let (code, path, message) = first_error(&model, &request);
    assert_eq!((code.as_str(), path.as_str()), ("CRN-R020", "/trackers"));
    assert!(message.contains("x is continuous"), "{message}");
    // The legacy ATP block needs thinning.
    let atp: Value = serde_json::from_str(ATP_LEGACY).unwrap();
    let (code, path, message) = first_error(&atp, &base);
    assert_eq!((code.as_str(), path.as_str()), ("CRN-R020", "/method"));
    assert!(message.contains("thinning"), "{message}");
    // Dynamic instances do not lower to fixed tuples.
    let clones: Value = serde_json::from_str(include_str!(
        "../../../../models/rules/mtdna-fixed.rules.json"
    ))
    .unwrap();
    assert_eq!(first_error(&clones, &base).1, "/method");
    // Generated kernels and loops need Static Direct.
    let prepared = prepare_json(&bytes, &base).unwrap();
    assert_eq!(
        generate_rules_kernel(&prepared).err().unwrap().code,
        "CRN-R020"
    );
    assert_eq!(
        generate_rules_loop(&prepared).err().unwrap().code,
        "CRN-R020"
    );
}

/// Without continuous variables `hybrid` runs a jump process; a derivative
/// that is not finite fails the replicate with its pointer.
#[test]
fn hybrid_runs_jump_processes_and_reports_bad_derivatives() {
    let model = birth_death(true, 20);
    let bytes = serde_json::to_vec(&model).unwrap();
    let request = hybrid_request(&bytes, 50, 2.0, &[0.0, 2.0]);
    let result = execute(
        &prepare_json(&bytes, &request).unwrap(),
        &Control::default(),
    )
    .unwrap();
    assert_eq!(result.termination, Termination::Completed);
    assert!(result.continuous.is_empty() && result.events > 0);
    let bad = continuous_model("log(x - 2)", 1.0, None, "k", 1);
    let result = hybrid_run(&bad, 1, 1.0, &[0.0, 1.0]);
    let diagnostic = result.replicates[0].diagnostic.as_ref().unwrap();
    assert_eq!(
        (diagnostic.code.as_str(), diagnostic.path.as_str()),
        ("CRN-R030", "/continuous/0/derivative")
    );
}
