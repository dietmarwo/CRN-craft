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
