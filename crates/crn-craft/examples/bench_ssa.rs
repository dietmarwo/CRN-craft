//! Timing harness for exact SSA engines on one prepared request.
//!
//! cargo run --release --example bench_ssa -- MODEL.json HORIZON PATHS METHOD WAIT REPEATS [fast]
//!
//! A trailing `fast` sends an `ssa-run/0.2` request with `method.fast: true`.
//! METHOD is `legacy_direct`, `legacy_cached` or an `ssa-run/0.1` aggregator
//! (`direct`, `cached_direct`, `direct_cr`, `sorting_direct`, `rssa_cr`); WAIT is
//! `inverse` or `ziggurat` (ignored by the legacy methods). Prints one JSON line
//! with every timed run (seconds, excluding preparation), the event count, a
//! SHA-256 of the whole result and a numerical SHA-256 that ignores work
//! counters (random words, evaluations, bound updates). Neither may change
//! between repeats.
use crn_craft::{sha256, simulation, ssa};
use serde_json::json;
use std::{hint::black_box, time::Instant};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let model = std::fs::read(&args[1]).expect("model");
    let horizon: f64 = args[2].parse().expect("horizon");
    let paths: u64 = args[3].parse().expect("paths");
    let method = args[4].as_str();
    let wait = args[5].as_str();
    let repeats: usize = args[6].parse().expect("repeats");
    let fast = args.get(7).is_some_and(|a| a == "fast");
    let reference = json!({"path": "model.json", "sha256": sha256(&model)});
    let time = json!({"end": horizon, "save_times": [0.0, horizon]});
    let ensemble = json!({"replicates": paths, "root_seed": 4242});
    let outputs = json!({"statistics": ["mean", "variance"], "keep_trajectories": 0});
    let control = simulation::Control::default();
    let mut times = Vec::new();
    let mut digests = Vec::new();
    let mut numerical = Vec::new();
    let mut events = 0;
    if let Some(variant) = method.strip_prefix("legacy_") {
        let run = serde_json::to_vec(&json!({
            "schema": "crn-craft.run/0.1", "model": reference,
            "method": {"kind": "ssa", "aggregator": "direct"},
            "time": time, "ensemble": ensemble, "outputs": outputs,
            "budget": {"max_events_per_replicate": 1_000_000_000u64, "threads": 1,
                       "max_output_bytes": 100_000_000u64}}))
        .unwrap();
        let prepared = simulation::prepare(&model, &run).expect("prepare");
        for repeat in 0..=repeats {
            let start = Instant::now();
            let result = if variant == "cached" {
                simulation::execute_cached(&prepared, &control)
            } else {
                simulation::execute(&prepared, &control)
            }
            .expect("execute");
            let elapsed = start.elapsed().as_secs_f64();
            events = result.result.events;
            let digest = sha256(&serde_json::to_vec(black_box(&result.result)).unwrap());
            if repeat > 0 {
                times.push(elapsed);
                numerical.push(digest.clone());
                digests.push(digest);
            }
        }
    } else {
        let (schema, method_spec) = if fast {
            (
                "crn-craft.ssa-run/0.2",
                json!({"kind": "ssa", "aggregator": method, "waiting_time": wait, "fast": true}),
            )
        } else {
            (
                "crn-craft.ssa-run/0.1",
                json!({"kind": "ssa", "aggregator": method, "waiting_time": wait}),
            )
        };
        let request = serde_json::to_vec(&json!({
            "schema": schema, "model": reference,
            "method": method_spec,
            "time": time, "ensemble": ensemble, "outputs": outputs,
            "budget": {"max_events_per_replicate": 1_000_000_000u64,
                       "max_sampling_steps_per_replicate": 10_000_000_000u64,
                       "threads": 1, "max_output_bytes": 100_000_000u64}}))
        .unwrap();
        let prepared = ssa::prepare(&model, &request).expect("prepare");
        for repeat in 0..=repeats {
            let start = Instant::now();
            let result = ssa::execute(&prepared, &control).expect("execute");
            let elapsed = start.elapsed().as_secs_f64();
            events = result.events;
            let digest = sha256(&serde_json::to_vec(black_box(&result)).unwrap());
            let mut value = serde_json::to_value(&result).unwrap();
            value["work"] = json!(null);
            for record in value["replicates"].as_array_mut().unwrap() {
                record["work"] = json!(null);
            }
            if repeat > 0 {
                times.push(elapsed);
                numerical.push(sha256(&serde_json::to_vec(&value).unwrap()));
                digests.push(digest);
            }
        }
    }
    digests.dedup();
    numerical.dedup();
    assert_eq!(digests.len(), 1, "results changed between repeats");
    assert_eq!(numerical.len(), 1, "results changed between repeats");
    println!(
        "{}",
        json!({"method": method, "wait": wait, "events": events, "seconds": times,
               "digest": digests[0], "numerical_digest": numerical[0]})
    );
}
