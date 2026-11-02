use crn_craft::{
    simulation::{self, Control},
    ssa,
};
use serde_json::{Value, json};

fn inputs(algorithm: &str, workers: u32) -> (Vec<u8>, Value) {
    let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"Alternatives",
        "convention":"stochastic_counts","omega":1.0,"parameters":{},
        "species":[{"name":"X","initial_count":3}],"reactions":[
            {"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":0.7}},
            {"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":1.0}}]})).unwrap();
    let run = json!({"schema":"crn-craft.ssa-run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(&model)},
        "method":{"kind":"ssa","aggregator":algorithm,"waiting_time":"inverse"},
        "time":{"end":10.0,"save_times":[0.0,1.0,3.0,10.0]},"ensemble":{"replicates":64,"root_seed":12345},
        "outputs":{"statistics":["mean","variance"],"keep_trajectories":64},
        "budget":{"threads":workers,"max_events_per_replicate":10000,"max_sampling_steps_per_replicate":100000,"max_output_bytes":10000000}});
    (model, run)
}
fn solve(m: &[u8], r: &Value) -> Value {
    serde_json::to_value(
        ssa::solve(m, &serde_json::to_vec(r).unwrap(), &Control::default()).unwrap(),
    )
    .unwrap()
}

#[test]
fn each_alternative_replays_across_threads_and_repeated_prepared_solves() {
    for algorithm in [
        "direct",
        "cached_direct",
        "direct_cr",
        "sorting_direct",
        "rssa_cr",
    ] {
        let (m, mut r) = inputs(algorithm, 1);
        for wait in ["inverse", "ziggurat"] {
            r["method"]["waiting_time"] = json!(wait);
            r["budget"]["threads"] = json!(1);
            let mut one = solve(&m, &r);
            assert_eq!(one["termination"], "completed");
            let prepared = ssa::prepare(&m, &serde_json::to_vec(&r).unwrap()).unwrap();
            assert_eq!(
                one,
                serde_json::to_value(ssa::execute(&prepared, &Control::default()).unwrap())
                    .unwrap()
            );
            assert_eq!(
                one,
                serde_json::to_value(ssa::execute(&prepared, &Control::default()).unwrap())
                    .unwrap()
            );
            r["budget"]["threads"] = json!(3);
            let mut many = solve(&m, &r);
            one.as_object_mut().unwrap().remove("run_sha256");
            many.as_object_mut().unwrap().remove("run_sha256");
            assert_eq!(one, many, "{algorithm}/{wait}");
        }
    }
}
#[test]
fn inverse_direct_variants_retain_legacy_paths() {
    for (algorithm, cached) in [("direct", false), ("cached_direct", true)] {
        let (m, r) = inputs(algorithm, 1);
        let new = solve(&m, &r);
        let mut legacy = r.clone();
        legacy["schema"] = json!("crn-craft.run/0.1");
        legacy["method"] = json!({"kind":"ssa","aggregator":"direct"});
        legacy["budget"]
            .as_object_mut()
            .unwrap()
            .remove("max_sampling_steps_per_replicate");
        let p = simulation::prepare(&m, &serde_json::to_vec(&legacy).unwrap()).unwrap();
        let old = if cached {
            simulation::execute_cached(&p, &Control::default())
        } else {
            simulation::execute(&p, &Control::default())
        }
        .unwrap();
        let old = serde_json::to_value(old.result).unwrap();
        for key in ["events", "trajectories", "statistics", "save_times"] {
            assert_eq!(new[key], old[key], "{algorithm}/{key}");
        }
    }
}
#[test]
fn sampling_limit_is_explicit_and_withholds_statistics() {
    for algorithm in [
        "direct",
        "cached_direct",
        "direct_cr",
        "sorting_direct",
        "rssa_cr",
    ] {
        let (m, mut r) = inputs(algorithm, 1);
        r["budget"]["max_sampling_steps_per_replicate"] = json!(1);
        let result = solve(&m, &r);
        assert_eq!(result["termination"], "sampling_limit");
        assert!(result["statistics"].is_null());
        assert_eq!(result["assessment"], "incomplete");
        for p in result["replicates"].as_array().unwrap() {
            assert!(p["work"]["rng_words"].as_u64().unwrap() <= 1);
        }
    }
}
#[test]
fn cancellation_and_strict_request_parsing_are_preserved() {
    let (m, mut r) = inputs("direct_cr", 2);
    let control = Control::default();
    control.cancel();
    let result = ssa::solve(&m, &serde_json::to_vec(&r).unwrap(), &control).unwrap();
    assert_eq!(result.termination, ssa::Termination::Cancelled);
    assert_eq!(result.work.rng_words, 0);
    r["unknown"] = json!(true);
    assert!(ssa::prepare(&m, &serde_json::to_vec(&r).unwrap()).is_err());
    let text = serde_json::to_string(&r).unwrap();
    let duplicated = format!("{{\"schema\":\"crn-craft.ssa-run/0.1\",{}", &text[1..]);
    assert!(ssa::prepare(&m, duplicated.as_bytes()).is_err());
}
#[test]
fn absorbing_network_consumes_no_sampling_budget() {
    let (m, mut r) = inputs("direct_cr", 1);
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["reactions"] = json!([]);
    let m = serde_json::to_vec(&model).unwrap();
    r["model"]["sha256"] = json!(crn_craft::sha256(&m));
    r["budget"]["max_sampling_steps_per_replicate"] = json!(0);
    let result = solve(&m, &r);
    assert_eq!(result["termination"], "completed");
    assert_eq!(result["work"]["rng_words"], 0);
}

#[test]
fn rssacr_detects_absorption_even_when_its_box_has_enabled_neighbors() {
    for wait in ["inverse", "ziggurat"] {
        let (m, mut r) = inputs("rssa_cr", 1);
        let mut model: Value = serde_json::from_slice(&m).unwrap();
        model["species"][0]["initial_count"] = json!(1);
        model["reactions"] = json!([{"name":"pair","reactants":{"X":2},"products":{},"rate":{"kind":"mass_action","coefficient":1.0}}]);
        let m = serde_json::to_vec(&model).unwrap();
        r["model"]["sha256"] = json!(crn_craft::sha256(&m));
        r["method"]["waiting_time"] = json!(wait);
        r["budget"]["max_sampling_steps_per_replicate"] = json!(0);
        let result = solve(&m, &r);
        assert_eq!(result["termination"], "completed");
        assert_eq!(result["work"]["rng_words"], 0);
    }
}

#[test]
fn rssacr_cooperatively_cancels_during_rejected_work() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    for wait in ["inverse", "ziggurat"] {
        let (m, mut r) = inputs("rssa_cr", 1);
        let mut model: Value = serde_json::from_slice(&m).unwrap();
        model["species"][0]["initial_count"] = json!(1);
        model["reactions"] = json!([{"name":"null","reactants":{"X":1},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1.0}}]);
        let m = serde_json::to_vec(&model).unwrap();
        r["model"]["sha256"] = json!(crn_craft::sha256(&m));
        r["method"]["waiting_time"] = json!(wait);
        r["ensemble"]["replicates"] = json!(1);
        r["outputs"]["keep_trajectories"] = json!(1);
        r["outputs"]["statistics"] = json!(["mean"]);
        r["time"] = json!({"end":10000.0,"save_times":[0.0,10000.0]});
        let calls = AtomicUsize::new(0);
        let mut control = Control::default();
        control.external_stop = Some(Box::new(move || calls.fetch_add(1, Ordering::SeqCst) >= 2));
        let result = ssa::solve(&m, &serde_json::to_vec(&r).unwrap(), &control).unwrap();
        assert_eq!(result.termination, ssa::Termination::Cancelled);
        assert!(result.work.rejected_candidates > 0);
        assert!(result.work.rng_words > 0 && result.work.rng_words < 100000);
        assert!(result.statistics.is_none());
    }
}

#[test]
fn the_output_budget_covers_exactly_the_compact_result_json() {
    // Eight-digit counts make the result larger than the request's minimum
    // numeric payload estimate (8 bytes per retained value).
    let (_, mut r) = inputs("direct_cr", 1);
    let m=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"Budget",
        "convention":"stochastic_counts","omega":1.0,"parameters":{},
        "species":[{"name":"X","initial_count":10_000_000}],"reactions":[
            {"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1.0}},
            {"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":1e-7}}]})).unwrap();
    r["model"]["sha256"] = json!(crn_craft::sha256(&m));
    let saves: Vec<f64> = (0..=1000).map(|i| f64::from(i) * 0.01).collect();
    r["time"]["save_times"] = json!(saves);
    let full = solve(&m, &r);
    assert_eq!(full["termination"], "completed");
    // The request digest has a fixed length, so the size does not depend on the budget.
    let size = serde_json::to_vec(&full).unwrap().len() as u64 + 1;
    r["budget"]["max_output_bytes"] = json!(size);
    let exact = solve(&m, &r);
    assert_eq!(exact["termination"], "completed", "budget {size}");
    assert_eq!(exact["trajectories"], full["trajectories"]);
    r["budget"]["max_output_bytes"] = json!(size - 1);
    let short = solve(&m, &r);
    assert_eq!(short["termination"], "output_limit");
    assert!(short["statistics"].is_null());
}

#[test]
fn fast_requires_ssa_run_0_2_and_replays_across_threads() {
    for algorithm in [
        "direct",
        "cached_direct",
        "direct_cr",
        "sorting_direct",
        "rssa_cr",
    ] {
        let (m, mut r) = inputs(algorithm, 1);
        r["method"]["fast"] = json!(true);
        let rejected = ssa::solve(&m, &serde_json::to_vec(&r).unwrap(), &Control::default());
        assert!(rejected.is_err(), "{algorithm}: fast under ssa-run/0.1");
        r["schema"] = json!("crn-craft.ssa-run/0.2");
        for wait in ["inverse", "ziggurat"] {
            r["method"]["waiting_time"] = json!(wait);
            r["budget"]["threads"] = json!(1);
            let mut one = solve(&m, &r);
            assert_eq!(one["termination"], "completed");
            assert_eq!(one["schema"], "crn-craft.ssa-result/0.2");
            assert_eq!(one["method"]["fast"], true);
            assert!(
                one["sampler"]
                    .as_str()
                    .unwrap()
                    .starts_with("xoshiro256pp-u53/")
            );
            r["budget"]["threads"] = json!(3);
            let mut many = solve(&m, &r);
            one.as_object_mut().unwrap().remove("run_sha256");
            many.as_object_mut().unwrap().remove("run_sha256");
            assert_eq!(one, many, "{algorithm}/{wait}");
            // fast=false under 0.2 is the default exact policy with a 0.2 envelope.
            r["method"]["fast"] = json!(false);
            r["budget"]["threads"] = json!(1);
            let exact = solve(&m, &r);
            assert_eq!(exact["schema"], "crn-craft.ssa-result/0.2");
            assert!(exact["method"].get("fast").is_none());
            assert!(
                exact["sampler"]
                    .as_str()
                    .unwrap()
                    .starts_with("chacha8-u53/")
            );
            r["method"]["fast"] = json!(true);
        }
    }
}

#[test]
fn fast_runs_still_honor_cancellation_and_sampling_budgets() {
    let (m, mut r) = inputs("direct_cr", 1);
    r["schema"] = json!("crn-craft.ssa-run/0.2");
    r["method"]["fast"] = json!(true);
    r["budget"]["max_sampling_steps_per_replicate"] = json!(5);
    let limited = solve(&m, &r);
    assert_eq!(limited["termination"], "sampling_limit");
    assert!(limited["statistics"].is_null());
    r["budget"]["max_sampling_steps_per_replicate"] = json!(100000);
    let prepared = ssa::prepare(&m, &serde_json::to_vec(&r).unwrap()).unwrap();
    let control = Control::default();
    control.cancel();
    let cancelled = ssa::execute(&prepared, &control).unwrap();
    assert_eq!(
        serde_json::to_value(&cancelled).unwrap()["termination"],
        "cancelled"
    );
}
