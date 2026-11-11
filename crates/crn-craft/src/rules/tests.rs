use super::*;
use serde_json::json;

fn request(model: &[u8], replicates: u64, end: f64, times: &[f64]) -> Value {
    json!({
        "schema": "crn-craft.rules-run/0.1",
        "model": {"path": "m.rules.json", "sha256": crate::sha256(model)},
        "time": {"end": end, "report_times": times},
        "ensemble": {"replicates": replicates, "root_seed": 11},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates.min(3), "keep_instances": replicates.min(3)},
        "budget": {"max_events_per_replicate": 1_000_000u64, "max_sampling_steps_per_replicate": 10_000_000u64,
                   "max_instances": 10_000u64, "threads": 2, "max_output_bytes": 100_000_000u64}
    })
}

fn run_with(
    model: &Value,
    edit: impl FnOnce(&mut Value),
    control: &Control,
) -> Result<RulesResult, Vec<Diagnostic>> {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut req = request(&bytes, 50, 2.0, &[0.0, 1.0, 2.0]);
    edit(&mut req);
    solve(&bytes, &serde_json::to_vec(&req).unwrap(), control)
}

fn run(model: &Value) -> RulesResult {
    run_with(model, |_| {}, &Control::default()).unwrap()
}

fn codes(model: &Value) -> Vec<String> {
    match check_model(&serde_json::to_vec(model).unwrap()) {
        RulesCheck { valid: true, .. } => Vec::new(),
        RulesCheck { diagnostics, .. } => diagnostics.into_iter().map(|d| d.code).collect(),
    }
}

fn aggregation(kernel: bool) -> Value {
    let factor = |l: &str| {
        if kernel {
            json!(format!("{l}.size"))
        } else {
            Value::Null
        }
    };
    let mut a = json!({"label": "a", "type": "X"});
    let mut b = json!({"label": "b", "type": "X"});
    if kernel {
        a["factor"] = factor("a");
        b["factor"] = factor("b");
    }
    json!({
        "schema": "crn-craft.rules-model/0.1", "name": "agg",
        "parameters": {"k": 1.0},
        "types": [{"name": "X", "attributes": ["size"]}],
        "initial": [{"type": "X", "attributes": {"size": 1}, "count": 12}],
        "observables": {"clusters": {"type": "X"}, "mass": {"type": "X", "weight": "x.size"},
                        "largest": {"type": "X", "weight": "x.size", "reduce": "max"},
                        "smallest": {"type": "X", "weight": "x.size", "reduce": "min"},
                        "sizes": {"type": "X", "reduce": "distinct"}},
        "rules": [{"name": "merge", "reactants": [a, b], "rate": "k",
                   "products": [{"type": "X", "attributes": {"size": "a.size + b.size"}}]}]
    })
}

fn clones() -> Value {
    json!({
        "schema": "crn-craft.rules-model/0.1", "name": "clones",
        "parameters": {"r": 1.0},
        "types": [{"name": "C", "attributes": ["id", "parent"]}],
        "initial": [{"type": "C", "attributes": {"id": 7, "parent": 0}, "count": 1}],
        "observables": {"cells": {"type": "C"}, "lineages": {"type": "C", "reduce": "distinct"}},
        "rules": [{"name": "divide", "reactants": [{"label": "x", "type": "C"}], "rate": "r",
                   "products": [{"type": "C", "copy": "x"},
                                {"type": "C", "attributes": {"id": "new_id()", "parent": "x.id"}}]}]
    })
}

#[test]
fn aggregation_conserves_mass_and_reports_reductions() {
    for kernel in [false, true] {
        let result = run(&aggregation(kernel));
        assert_eq!(result.termination, Termination::Completed);
        assert_eq!(
            result.observables,
            ["clusters", "largest", "mass", "sizes", "smallest"]
        );
        for path in &result.paths {
            for row in &path.values {
                assert_eq!(row[2], 12.0);
                assert!(row[1] >= row[4] && row[3] <= row[0]);
            }
            assert_eq!(path.values[0], vec![12.0, 1.0, 12.0, 1.0, 1.0]);
        }
        for f in &result.final_instances {
            let mass: f64 = f
                .instances
                .iter()
                .map(|i| i.attributes[0] * i.count as f64)
                .sum();
            assert_eq!(mass, 12.0);
            // Equal attributes are merged into one instance.
            let sizes: Vec<_> = f.instances.iter().map(|i| i.attributes[0]).collect();
            assert!(sizes.windows(2).all(|w| w[0] < w[1]));
        }
        assert!(result.replicates.iter().all(|r| r.events <= 11));
    }
}

#[test]
fn new_id_creates_distinct_instances() {
    let result = run(&clones());
    // id 7 is the largest initial attribute, so new identities start at 8.
    for f in &result.final_instances {
        let ids: Vec<_> = f.instances.iter().map(|i| i.attributes[0]).collect();
        assert!(
            f.instances
                .iter()
                .all(|i| i.count == 1 || i.attributes[0] == 7.0)
        );
        assert!(ids.iter().skip(1).all(|&id| id >= 8.0));
        let mut sorted = ids.clone();
        sorted.dedup();
        assert_eq!(sorted.len(), ids.len());
    }
    for (path, record) in result.paths.iter().zip(&result.replicates) {
        let last = path.values.last().unwrap();
        // Every division creates a lineage; the founder keeps one cell.
        assert_eq!(last[1], last[0]);
        assert_eq!(last[0] as u64, 1 + record.events);
        assert_eq!(record.created_instances, record.events);
    }
}

#[test]
fn replay_is_independent_of_threads() {
    let model = serde_json::to_vec(&aggregation(true)).unwrap();
    let mut outputs = Vec::new();
    for threads in [1, 3, 8] {
        let mut req = request(&model, 200, 0.5, &[0.1, 0.5]);
        req["budget"]["threads"] = json!(threads);
        req["outputs"]["keep_paths"] = json!(200);
        let mut result = solve(
            &model,
            &serde_json::to_vec(&req).unwrap(),
            &Control::default(),
        )
        .unwrap();
        result.run_sha256.clear();
        outputs.push(serde_json::to_string(&result).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0], outputs[2]);
}

#[test]
fn stop_condition_holds_the_state() {
    let mut model = clones();
    model["stop"] = json!("cells >= 5");
    let result = run_with(
        &model,
        |r| {
            r["time"]["end"] = json!(100.0);
            r["time"]["report_times"] = json!([0.0, 100.0]);
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(result.termination, Termination::Completed);
    assert_eq!(result.stopped_replicates, 50);
    for (path, record) in result.paths.iter().zip(&result.replicates) {
        assert_eq!(record.termination, Termination::Stopped);
        assert!(record.stop_time < 100.0);
        assert_eq!(record.events, 4);
        assert_eq!(path.values.last().unwrap()[0], 5.0);
    }
    assert!(result.statistics.is_some());
    // A stop condition true at t = 0 stops before any event.
    model["stop"] = json!("cells >= 1");
    let result = run(&model);
    assert!(
        result
            .replicates
            .iter()
            .all(|r| r.events == 0 && r.stop_time == 0.0)
    );
}

#[test]
fn budgets_and_cancellation_withhold_statistics() {
    let model = clones();
    let events = run_with(
        &model,
        |r| {
            r["budget"]["max_events_per_replicate"] = json!(2);
            r["time"]["end"] = json!(50.0);
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(events.termination, Termination::EventLimit);
    assert!(events.statistics.is_none());
    assert!(
        events
            .replicates
            .iter()
            .all(|r| r.events == 2 && r.diagnostic.as_ref().unwrap().code == "CRN-R040")
    );
    let words = run_with(
        &model,
        |r| {
            r["budget"]["max_sampling_steps_per_replicate"] = json!(5);
            r["time"]["end"] = json!(50.0);
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(words.termination, Termination::SamplingLimit);
    assert!(words.replicates.iter().all(|r| r.rng_words <= 5));
    let instances = run_with(
        &model,
        |r| {
            r["budget"]["max_instances"] = json!(3);
            r["time"]["end"] = json!(50.0);
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(instances.termination, Termination::InstanceLimit);
    assert!(instances.replicates.iter().all(|r| r.peak_instances == 4));
    let full = run_with(
        &model,
        |r| r["outputs"]["keep_paths"] = json!(50),
        &Control::default(),
    )
    .unwrap();
    let mut stripped = full.clone();
    stripped.paths.clear();
    stripped.final_instances.clear();
    stripped.statistics = None;
    let middle = (serde_json::to_vec(&full).unwrap().len()
        + serde_json::to_vec(&stripped).unwrap().len())
        / 2;
    let output = run_with(
        &model,
        |r| {
            r["outputs"]["keep_paths"] = json!(50);
            r["budget"]["max_output_bytes"] = json!(middle);
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(output.termination, Termination::OutputLimit);
    assert!(output.paths.is_empty() && output.statistics.is_none());
    assert!(
        run_with(
            &model,
            |r| r["budget"]["max_output_bytes"] = json!(100),
            &Control::default()
        )
        .is_err()
    );
    let control = Control::default();
    control.cancel();
    let cancelled = run_with(&model, |_| {}, &control).unwrap();
    assert_eq!(cancelled.termination, Termination::Cancelled);
    assert!(cancelled.statistics.is_none());
    assert_eq!(cancelled.completed_replicates, 0);
}

#[test]
fn extinction_does_not_evaluate_undefined_rates() {
    // D = 0 after extinction; the replication rate K/D is never evaluated
    // because no reactant remains.
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "extinct",
        "parameters": {"K": 2, "delta": 50.0},
        "types": [{"name": "M", "attributes": []}],
        "initial": [{"type": "M", "count": 2}],
        "observables": {"N": {"type": "M"}},
        "rules": [
            {"name": "die", "reactants": [{"label": "x", "type": "M"}], "rate": "delta", "products": []},
            {"name": "copy", "reactants": [{"label": "x", "type": "M"}], "rate": "max(0, K - N) / N",
             "products": [{"type": "M", "copy": "x"}, {"type": "M", "copy": "x"}]}]
    });
    let result = run(&model);
    assert_eq!(result.termination, Termination::Completed);
    assert!(result.paths.iter().all(|p| p.values[2][0] == 0.0));
}

#[test]
fn numerical_failures_are_reported() {
    let mut model = clones();
    model["rules"][0]["rate"] = json!("r - 2");
    let result = run(&model);
    assert_eq!(result.termination, Termination::Failed);
    assert_eq!(
        result.replicates[0].diagnostic.as_ref().unwrap().code,
        "CRN-R030"
    );
    // A weight that 53-bit uniforms cannot resolve.
    let mut model = aggregation(true);
    model["initial"] = json!([{"type": "X", "attributes": {"size": 1e-20}, "count": 1},
                             {"type": "X", "attributes": {"size": 1}, "count": 3}]);
    let result = run(&model);
    assert_eq!(
        result.replicates[0].diagnostic.as_ref().unwrap().code,
        "CRN-R031"
    );
}

#[test]
fn strict_model_input() {
    let base = clones();
    assert!(codes(&base).is_empty());
    let with = |path: &str, value: Value| {
        let mut m = base.clone();
        let (parent, key) = path.rsplit_once('/').unwrap();
        match m.pointer_mut(parent).unwrap() {
            Value::Object(object) => {
                object.insert(key.to_string(), value);
            }
            other => *other.pointer_mut(&format!("/{key}")).unwrap() = value,
        }
        codes(&m)
    };
    let mut unknown = base.clone();
    unknown["rules"][0]["speed"] = json!(1);
    assert_eq!(codes(&unknown), ["CRN-E003"]);
    let duplicate = br#"{"schema":"crn-craft.rules-model/0.1","name":"a","name":"b","types":[],"initial":[],"rules":[]}"#;
    assert!(!check_model(duplicate).valid);
    assert_eq!(with("/types/0/name", json!("2C"))[0], "CRN-R010");
    assert_eq!(
        with("/rules/0/rate", json!("cells / max(1)"))[0],
        "CRN-R012"
    );
    assert_eq!(with("/rules/0/rate", json!("x.id"))[0], "CRN-R012");
    assert_eq!(with("/rules/0/rate", json!("uniform(0, 1)"))[0], "CRN-R012");
    assert_eq!(
        with("/rules/0/reactants/0/factor", json!("cells"))[0],
        "CRN-R012"
    );
    assert_eq!(with("/rules/0/products/0/copy", json!("y"))[0], "CRN-R013");
    assert_eq!(
        with("/rules/0/products/1/attributes", json!({"id": 1}))[0],
        "CRN-R013"
    );
    let extra = json!({"id": "new_id()", "parent": "x.id", "size": 1});
    assert_eq!(with("/rules/0/products/1/attributes", extra)[0], "CRN-R013");
    assert_eq!(with("/initial/0/count", json!(1.5))[0], "CRN-R014");
    assert_eq!(with("/initial/0/count", json!(-1))[0], "CRN-R014");
    assert_eq!(with("/initial/0/type", json!("D"))[0], "CRN-R011");
    assert_eq!(with("/observables/cells/type", json!("D"))[0], "CRN-R011");
    assert_eq!(
        with("/parameters", json!({"cells": 1, "r": 1}))[0],
        "CRN-R010"
    );
    assert_eq!(
        with("/parameters", json!({"min": 1, "r": 1}))[0],
        "CRN-R010"
    );
    let three = json!([{"label": "a", "type": "C"}, {"label": "b", "type": "C"}, {"label": "c", "type": "C"}]);
    assert_eq!(with("/rules/0/reactants", three)[0], "CRN-R013");
    let differ =
        json!([{"label": "a", "type": "C"}, {"label": "b", "type": "C", "where": "b.id > 3"}]);
    assert!(with("/rules/0/reactants", differ).contains(&"CRN-R013".to_string()));
    let same_label = json!([{"label": "a", "type": "C"}, {"label": "a", "type": "C"}]);
    assert_eq!(with("/rules/0/reactants", same_label)[0], "CRN-R010");
}

#[test]
fn strict_request_input() {
    let model = serde_json::to_vec(&clones()).unwrap();
    let check = |edit: &dyn Fn(&mut Value)| {
        let mut req = request(&model, 5, 1.0, &[0.0, 1.0]);
        edit(&mut req);
        prepare(&model, &serde_json::to_vec(&req).unwrap())
            .err()
            .map(|d| d[0].code.clone())
    };
    assert_eq!(check(&|_| {}), None);
    assert_eq!(
        check(&|r| r["model"]["sha256"] = json!("0".repeat(64))).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["time"]["report_times"] = json!([1.0, 0.5])).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["time"]["report_times"] = json!([2.0])).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["time"]["end"] = json!(0.0)).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["outputs"]["keep_paths"] = json!(6)).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["budget"]["threads"] = json!(0)).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["model"]["path"] = json!("../m.json")).unwrap(),
        "CRN-R020"
    );
    assert_eq!(
        check(&|r| r["budget"]["extra"] = json!(1)).unwrap(),
        "CRN-E003"
    );
}

#[test]
fn pair_needs_two_molecules() {
    // A single molecule has no pair: the rule never fires.
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "pair",
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 1}],
        "observables": {"X": {"type": "X"}},
        "rules": [{"name": "annihilate", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "X"}],
                   "rate": "1", "products": []}]
    });
    let result = run(&model);
    assert!(result.replicates.iter().all(|r| r.events == 0));
    let mut two = model;
    two["initial"][0]["count"] = json!(2);
    let result = run_with(
        &two,
        |r| r["time"]["end"] = json!(100.0),
        &Control::default(),
    )
    .unwrap();
    assert!(result.replicates.iter().all(|r| r.events == 1));
}

#[test]
fn cross_type_rules() {
    // A + B -> C at rate 1 per pair; stoichiometry and attribute passing.
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "cross",
        "types": [{"name": "A", "attributes": ["v"]}, {"name": "B", "attributes": ["w"]}, {"name": "C", "attributes": ["sum"]}],
        "initial": [{"type": "A", "attributes": {"v": 1}, "count": 3}, {"type": "B", "attributes": {"w": 10}, "count": 2}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}, "C": {"type": "C", "weight": "x.sum"}},
        "rules": [{"name": "bind", "reactants": [{"label": "a", "type": "A"}, {"label": "b", "type": "B"}],
                   "rate": "1", "products": [{"type": "C", "attributes": {"sum": "a.v + b.w"}}]}]
    });
    let result = run_with(
        &model,
        |r| r["time"]["end"] = json!(100.0),
        &Control::default(),
    )
    .unwrap();
    for p in &result.paths {
        assert_eq!(p.values[2], vec![1.0, 0.0, 22.0]);
    }
}

#[test]
fn schemas_are_generated() {
    for kind in ["rules-model", "rules-check", "rules-run", "rules-result"] {
        assert!(schema(kind).is_some(), "{kind}");
    }
    assert!(schema("model").is_none());
}

#[test]
fn review_regressions() {
    // Large attributes without new_id() run; new_id() near 2^53 is checked.
    let big = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "big",
        "types": [{"name": "P", "attributes": ["mass"]}],
        "initial": [{"type": "P", "attributes": {"mass": 1e16}, "count": 3}],
        "observables": {"P": {"type": "P"}},
        "rules": [{"name": "decay", "reactants": [{"label": "p", "type": "P"}], "rate": "1", "products": []}]
    });
    assert_eq!(run(&big).termination, Termination::Completed);
    let mut ids = clones();
    ids["initial"][0]["attributes"]["id"] = json!(9_007_199_254_740_991.0_f64);
    assert_eq!(codes(&ids), ["CRN-R014"]);
    ids["initial"][0]["attributes"]["id"] = json!(9_007_199_254_740_990.0_f64);
    let result = run_with(
        &ids,
        |r| r["time"]["end"] = json!(20.0),
        &Control::default(),
    )
    .unwrap();
    let failed = result.replicates.iter().find(|r| r.events >= 1).unwrap();
    assert_eq!(failed.events, 1);
    assert_eq!(failed.diagnostic.as_ref().unwrap().code, "CRN-R030");
    // Overflowing pair weights and observable totals fail instead of
    // becoming 0 or null.
    let mut pairs = aggregation(true);
    pairs["initial"] = json!([{"type": "X", "attributes": {"size": 1e200}, "count": 2}]);
    pairs["parameters"]["k"] = json!(1e-300);
    let result = run(&pairs);
    assert_eq!(result.termination, Termination::Failed);
    assert_eq!(
        result.replicates[0].diagnostic.as_ref().unwrap().code,
        "CRN-R030"
    );
    let mut totals = aggregation(false);
    totals["initial"] = json!([{"type": "X", "attributes": {"size": 1e308}, "count": 1},
                              {"type": "X", "attributes": {"size": 1.5e308}, "count": 1}]);
    let result = run(&totals);
    assert_eq!(
        result.replicates[0].diagnostic.as_ref().unwrap().code,
        "CRN-R030"
    );
    // An undefined condition is an error, not true.
    let mut nan = clones();
    nan["rules"][0]["reactants"][0]["where"] = json!("x.parent / x.parent");
    let result = run(&nan);
    assert_eq!(
        result.replicates[0].diagnostic.as_ref().unwrap().code,
        "CRN-R030"
    );
    // choice() evaluates only the chosen branch: one identity per division.
    let mut choice = clones();
    choice["rules"][0]["products"][1]["attributes"]["id"] = json!("choice(1, new_id(), new_id())");
    let result = run_with(
        &choice,
        |r| r["time"]["end"] = json!(3.0),
        &Control::default(),
    )
    .unwrap();
    for f in &result.final_instances {
        let mut ids: Vec<_> = f.instances.iter().map(|i| i.attributes[0]).collect();
        ids.sort_by(f64::total_cmp);
        let expected: Vec<f64> = (0..ids.len())
            .map(|k| if k == 0 { 7.0 } else { 7.0 + k as f64 })
            .collect();
        assert_eq!(ids, expected);
    }
}

#[test]
fn diagnostics_paths_and_admission() {
    let mut model = clones();
    model["rules"][0]["reactants"] =
        json!([{"label": "a", "type": "C"}, {"label": "a", "type": "C"}]);
    let report = check_model(&serde_json::to_vec(&model).unwrap());
    assert_eq!(report.diagnostics[0].path, "/rules/0/reactants/1/label");
    let mut unnamed = clones();
    unnamed["name"] = json!("");
    assert_eq!(codes(&unnamed), ["CRN-R010"]);
    let mut spaced = clones();
    spaced["rules"][0]["rate"] = json!("r\u{00a0}* 1");
    assert_eq!(codes(&spaced), ["CRN-R012"]);
    // Oversized requests are rejected before running.
    let bytes = serde_json::to_vec(&clones()).unwrap();
    let mut req = request(&bytes, 1_000_000, 1.0, &[1.0]);
    req["outputs"]["keep_paths"] = json!(0);
    req["outputs"]["keep_instances"] = json!(0);
    req["budget"]["max_output_bytes"] = json!(1000);
    let errors = prepare(&bytes, &serde_json::to_vec(&req).unwrap())
        .err()
        .unwrap();
    assert!(errors.iter().any(|d| d.path == "/budget/max_output_bytes"));
    let times: Vec<f64> = (0..100_000).map(|k| k as f64 * 1e-5).collect();
    let mut req = request(&bytes, 100_000, 1.0, &times);
    req["outputs"]["keep_paths"] = json!(100_000);
    req["outputs"]["keep_instances"] = json!(0);
    let errors = prepare(&bytes, &serde_json::to_vec(&req).unwrap())
        .err()
        .unwrap();
    assert!(
        errors
            .iter()
            .any(|d| d.code == "CRN-R020" && d.path == "/budget")
    );
}

#[test]
fn pair_selection_is_exact_for_heterogeneous_factors() {
    // Instances A (factor 1, two copies) and B (factor 10, one copy): pair
    // weights AA 1, AB 2 × 10 = 20; first merge is AA with probability 1/21.
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "pairs",
        "types": [{"name": "X", "attributes": ["h"]}],
        "initial": [{"type": "X", "attributes": {"h": 1}, "count": 2}, {"type": "X", "attributes": {"h": 10}, "count": 1}],
        "observables": {"n": {"type": "X"}},
        "stop": "n <= 2",
        "rules": [{"name": "merge", "reactants": [{"label": "a", "type": "X", "factor": "a.h"}, {"label": "b", "type": "X", "factor": "b.h"}],
                   "rate": "1", "products": [{"type": "X", "attributes": {"h": "a.h + b.h + 100"}}]}]
    });
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request(&bytes, 42_000, 0.001, &[0.001]);
    req["time"]["end"] = json!(1.0);
    req["time"]["report_times"] = json!([1e-9, 1.0]);
    req["outputs"]["keep_paths"] = json!(0);
    req["outputs"]["keep_instances"] = json!(42_000);
    req["budget"]["max_output_bytes"] = json!(200_000_000u64);
    req["budget"]["max_instances"] = json!(3);
    let result = solve(
        &bytes,
        &serde_json::to_vec(&req).unwrap(),
        &Control::default(),
    )
    .unwrap();
    // The run stops after the first merge; AA leaves h = 102 behind.
    let aa = result
        .final_instances
        .iter()
        .filter(|f| f.instances.iter().any(|i| i.attributes[0] == 102.0))
        .count();
    let p = aa as f64 / 42_000.0;
    // 1/21 = 0.0476; 5 standard errors is 0.0052.
    assert!((p - 1.0 / 21.0).abs() < 0.0052, "{p}");
}

#[test]
fn generic_sparse_tracking_is_independent_of_continue_observer() {
    struct Continue;
    impl EventObserver for Continue {
        fn on_event(&self, _: &observe::Event<'_>) -> observe::ObserverAction {
            observe::ObserverAction::Continue
        }
    }
    let model = include_bytes!("../../../../models/rules/mtdna-three-class-cme.rules.json");
    let mut run = request(model, 16, 10.0, &[0.0, 1.0, 10.0]);
    run["schema"] = json!("crn-craft.rules-run/0.4");
    run["method"] =
        json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true});
    run["trackers"] =
        json!([{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]);
    let p = prepare(model, &serde_json::to_vec(&run).unwrap()).unwrap();
    let control = Control::default();
    let normal = execute(&p, &control).unwrap();
    let observed = execute_with_observer(&p, &control, &Continue).unwrap();
    assert!(normal.events > 0);
    assert_eq!(
        serde_json::to_vec(&normal).unwrap(),
        serde_json::to_vec(&observed).unwrap()
    );
    run["schema"] = json!("crn-craft.rules-run/0.3");
    assert!(prepare(model, &serde_json::to_vec(&run).unwrap()).is_err());
}

#[test]
fn generic_ziggurat_infinite_wait_completes_horizon() {
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "tiny",
        "types": [{"name": "X"}], "initial": [],
        "observables": {"X": {"type": "X"}},
        "rules": [{"name": "birth", "reactants": [], "rate": "5e-324", "products": [{"type": "X"}]}]
    });
    let result = run_with(
        &model,
        |run| {
            run["schema"] = json!("crn-craft.rules-run/0.4");
            run["method"] = json!({"kind": "static_direct", "waiting_time": "ziggurat"});
        },
        &Control::default(),
    )
    .unwrap();
    assert_eq!(result.termination, Termination::Completed);
    assert_eq!(result.events, 0);
    assert_eq!(result.paths[0].values, vec![vec![0.0]; 3]);
}

#[test]
fn unchanged_rate_is_rechecked_when_total_increases() {
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "stiff_after_birth",
        "types": [{"name": "X"}], "initial": [],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {"name": "birth", "reactants": [], "rate": "1", "products": [{"type": "X"}]},
            {"name": "rare", "reactants": [], "rate": "2e-16", "products": []},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "2", "products": []}
        ]
    });
    let result = run_with(
        &model,
        |run| {
            run["schema"] = json!("crn-craft.rules-run/0.2");
            run["method"] = json!({"kind": "static_direct"});
        },
        &Control::default(),
    )
    .unwrap();
    assert!(result.replicates.iter().any(|rep| {
        rep.events == 1
            && rep
                .diagnostic
                .as_ref()
                .is_some_and(|d| d.code == "CRN-R031")
    }));
}
