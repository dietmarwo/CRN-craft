use crn_craft::{
    sampling, sha256,
    simulation::{self, Assessment, Control, Stop},
};
use rand_core::RngCore;
use serde_json::{Value, json};

fn model() -> Value {
    json!({"schema":"crn-craft.model/0.1","name":"ssa","convention":"stochastic_counts","omega":7.0,"parameters":{},
        "species":[{"name":"X","initial_count":0}],
        "reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":3.0}},
            {"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":0.8}}]})
}
fn run(model: &[u8]) -> Value {
    json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(model)},
        "method":{"kind":"ssa","aggregator":"direct"},"time":{"end":2.0,"save_times":[0.0,0.25,1.0,2.0]},
        "ensemble":{"replicates":32,"root_seed":42},"outputs":{"statistics":["mean","variance"],"keep_trajectories":32},
        "budget":{"max_events_per_replicate":10000,"threads":1,"max_output_bytes":4000000}})
}
fn simulate(model: &Value, edit: impl FnOnce(&mut Value)) -> simulation::SimulationResult {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut run = run(&bytes);
    edit(&mut run);
    simulation::simulate_bytes(
        &bytes,
        &serde_json::to_vec(&run).unwrap(),
        &Control::default(),
    )
    .unwrap()
    .result
}

#[test]
fn chacha8_matches_independent_scalar_known_answers_across_blocks_and_streams() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/rng-vectors.json")).unwrap();
    for row in fixture["streams"].as_array().unwrap() {
        let mut rng = sampling::stream(
            row["root"].as_u64().unwrap(),
            row["replicate"].as_u64().unwrap(),
        );
        for word in row["words"].as_array().unwrap() {
            assert_eq!(rng.next_u64(), word.as_u64().unwrap());
        }
    }
}

#[test]
fn sampler_endpoints_waits_selection_and_precision_failures() {
    assert_eq!(sampling::uniform_closed_open(0), 0.0);
    assert_eq!(
        sampling::uniform_closed_open(u64::MAX),
        1.0 - sampling::UNIFORM_STEP
    );
    assert_eq!(
        sampling::waiting_time(2.0, 0.5).unwrap(),
        libm::log(2.0) / 2.0
    );
    for u in [0.0, 1.0, -1.0, f64::NAN] {
        assert!(sampling::waiting_time(1.0, u).is_err());
    }
    assert_eq!(
        sampling::select_channel(&[0.0, 1.0, 0.0, 3.0], 4.0, 0.0).unwrap(),
        1
    );
    assert_eq!(
        sampling::select_channel(&[0.0, 1.0, 0.0, 3.0], 4.0, 0.25).unwrap(),
        3
    );
    assert_eq!(
        sampling::select_channel(&[1.0, 3.0], 4.0, 1.0 - sampling::UNIFORM_STEP).unwrap(),
        1
    );
    assert!(sampling::select_channel(&[0.0], 1.0, 0.5).is_err());
    assert!(sampling::total_propensity(&[1e100, 1.0]).is_err());
    assert!(sampling::total_propensity(&[1.0, 1e100]).is_err());
    assert!(sampling::total_propensity(&[f64::MAX, f64::MAX]).is_err());
    assert!(sampling::total_propensity(&[f64::NAN]).is_err());
    assert!(sampling::event_time(9_007_199_254_740_992.0, 0.5).is_err());
    assert!(sampling::event_time(f64::MAX, f64::MAX).is_err());
    assert!(sampling::waiting_time(f64::from_bits(1), 0.5).is_err());
}

#[test]
fn fixed_work_payload_is_identical_at_one_two_and_eight_threads() {
    let first = simulate(&model(), |_| {});
    assert_eq!(first.termination, Stop::Completed);
    assert_eq!(first.assessment, Assessment::ExactEnsemble);
    for threads in [2, 8] {
        let other = simulate(&model(), |r| r["budget"]["threads"] = json!(threads));
        let mut lhs = serde_json::to_value(&first).unwrap();
        let mut rhs = serde_json::to_value(&other).unwrap();
        // Raw run bytes bind the thread allowance; numerical payload is separate.
        lhs.as_object_mut().unwrap().remove("run_sha256");
        rhs.as_object_mut().unwrap().remove("run_sha256");
        assert_eq!(lhs, rhs);
    }
    assert_ne!(first.trajectories[0].counts, first.trajectories[1].counts);
}

#[test]
fn statistics_match_retained_integer_paths_and_unbiased_variance() {
    let result = simulate(&model(), |_| {});
    let stats = result.statistics.unwrap();
    for time in 0..4 {
        let values: Vec<f64> = result
            .trajectories
            .iter()
            .map(|p| p.counts[time][0] as f64)
            .collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance =
            values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
        assert!((stats.mean.as_ref().unwrap()[time][0] - mean).abs() < 1e-13);
        assert!((stats.sample_variance.as_ref().unwrap()[time][0] - variance).abs() < 1e-12);
        assert!(
            (stats.mean_standard_error.as_ref().unwrap()[time][0] - libm::sqrt(variance / 32.0))
                .abs()
                < 1e-13
        );
    }
}

#[test]
fn save_density_does_not_change_events_or_common_states() {
    let sparse = simulate(&model(), |_| {});
    let dense = simulate(&model(), |r| {
        r["time"]["save_times"] = json!([0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    });
    assert_eq!(sparse.events, dense.events);
    for (a, b) in sparse.trajectories.iter().zip(dense.trajectories) {
        for (i, j) in [0, 2, 5, 7].into_iter().enumerate() {
            assert_eq!(a.counts[i], b.counts[j]);
        }
    }
}

#[test]
fn absorption_zero_rates_and_single_replicate_se() {
    let mut m = model();
    m["reactions"] = json!([]);
    let r = simulate(&m, |r| {
        r["ensemble"]["replicates"] = json!(1);
        r["outputs"]["keep_trajectories"] = json!(1);
        r["outputs"]["statistics"] = json!(["mean"]);
    });
    assert_eq!(r.events, 0);
    assert_eq!(r.termination, Stop::Completed);
    assert_eq!(r.trajectories[0].counts, vec![vec![0]; 4]);
    assert!(r.statistics.unwrap().mean_standard_error.is_none());
}

#[test]
fn event_cap_partial_paths_and_no_conditional_statistics() {
    let r = simulate(&model(), |r| {
        r["budget"]["max_events_per_replicate"] = json!(1)
    });
    assert_eq!(r.termination, Stop::EventLimit);
    assert_eq!(r.assessment, Assessment::Incomplete);
    assert!(r.statistics.is_none());
    assert!(r.replicates.iter().all(|r| r.events <= 1));
    assert!(
        r.replicates
            .iter()
            .filter(|r| r.termination == Stop::EventLimit)
            .all(|r| r.events == 1)
    );
}

#[test]
fn mixed_completion_never_reports_conditional_ensemble_statistics() {
    let mut m = model();
    m["reactions"].as_array_mut().unwrap().truncate(1);
    m["reactions"][0]["rate"]["coefficient"] = json!(0.5);
    let r = simulate(&m, |r| r["budget"]["max_events_per_replicate"] = json!(1));
    assert!(r.completed_replicates > 0 && r.completed_replicates < r.requested_replicates);
    assert!(r.statistics.is_none());
    assert_eq!(r.assessment, Assessment::Incomplete);
}

#[test]
fn cancellation_and_wall_time_account_for_every_unstarted_replicate() {
    let bytes = serde_json::to_vec(&model()).unwrap();
    let r = run(&bytes);
    let control = Control::default();
    control.cancel();
    let cancelled = simulation::simulate_bytes(&bytes, &serde_json::to_vec(&r).unwrap(), &control)
        .unwrap()
        .result;
    assert_eq!(cancelled.termination, Stop::Cancelled);
    assert_eq!(cancelled.replicates.len(), 32);
    assert!(
        cancelled
            .replicates
            .iter()
            .all(|r| !r.started && r.events == 0)
    );
    assert!(cancelled.statistics.is_none());
    let expired = simulate(&model(), |r| r["budget"]["wall_time_s"] = json!(1e-300));
    assert_eq!(expired.termination, Stop::WallTimeLimit);
    assert!(expired.replicates.iter().all(|r| !r.started));
}

#[test]
fn population_stop_and_integer_overflow_preserve_old_counts() {
    let mut m = model();
    m["species"][0]["population_limit"] = json!(0);
    let limited = simulate(&m, |_| {});
    assert_eq!(limited.termination, Stop::PopulationLimit);
    assert!(limited.statistics.is_none());
    assert_eq!(limited.events, 0);
    assert!(
        limited
            .trajectories
            .iter()
            .flat_map(|p| &p.counts)
            .all(|c| c[0] == 0)
    );
    m["species"][0] = json!({"name":"X","initial_count":i64::MAX});
    m["reactions"].as_array_mut().unwrap().truncate(1);
    let overflow = simulate(&m, |_| {});
    assert_eq!(overflow.termination, Stop::Failed);
    assert!(
        overflow
            .replicates
            .iter()
            .filter_map(|r| r.diagnostic.as_ref())
            .all(|d| d.code == "CRN-E044")
    );
}

#[test]
fn underflow_and_dynamic_range_errors_are_not_disabled_channels() {
    let mut m = model();
    m["species"][0]["initial_count"] = json!(1);
    m["reactions"].as_array_mut().unwrap().truncate(1);
    m["reactions"][0]["rate"] = json!({"kind":"hill_activation","regulator":"X","vmax":1.0,"basal":0.0,"k":1e300,"n":100.0});
    let r = simulate(&m, |_| {});
    assert_eq!(r.termination, Stop::Failed);
    assert!(
        r.replicates
            .iter()
            .all(|r| r.diagnostic.as_ref().unwrap().code == "CRN-E042")
    );
    let mut m = model();
    m["species"][0]["initial_count"] = json!(1);
    m["reactions"][0]["rate"]["coefficient"] = json!(1e100);
    let r = simulate(&m, |_| {});
    assert!(
        r.replicates
            .iter()
            .all(|r| r.diagnostic.as_ref().unwrap().code == "CRN-E061")
    );
}

#[test]
fn null_reactions_count_as_events_and_respect_event_cap() {
    let mut m = model();
    m["reactions"] = json!([{"name":"null","reactants":{},"products":{},"rate":{"kind":"mass_action","coefficient":1e6}}]);
    let r = simulate(&m, |r| r["budget"]["max_events_per_replicate"] = json!(2));
    assert_eq!(r.events, 64);
    assert_eq!(r.termination, Stop::EventLimit);
}

#[test]
fn unreasonable_memory_work_and_fsp_requests_fail_before_execution() {
    let bytes = serde_json::to_vec(&model()).unwrap();
    for edit in [0, 1, 2] {
        let mut r = run(&bytes);
        match edit {
            0 => {
                r["ensemble"]["replicates"] = json!(u64::MAX);
                r["outputs"]["keep_trajectories"] = json!(0);
            }
            1 => r["budget"]["max_events_per_replicate"] = json!(u64::MAX),
            _ => {
                r["method"] =
                    json!({"kind":"fsp","bounds":{"X":20},"max_states":21,"tolerance":1e-8});
                r.as_object_mut().unwrap().remove("ensemble");
                r["budget"]
                    .as_object_mut()
                    .unwrap()
                    .remove("max_events_per_replicate");
                r["outputs"]["keep_trajectories"] = json!(0);
            }
        }
        let errors = simulation::prepare(&bytes, &serde_json::to_vec(&r).unwrap()).unwrap_err();
        assert!(errors.iter().any(|e| e.code == "CRN-E060"), "{errors:?}");
    }
}

#[test]
fn live_serialized_output_cap_keeps_ledger_and_discards_statistics() {
    let mut m = model();
    m["species"][0]["initial_count"] = json!(i64::MAX);
    m["reactions"] = json!([]);
    let r = simulate(&m, |r| {
        r["time"]["save_times"] =
            json!((0..200).map(|i| i as f64 * 2.0 / 199.0).collect::<Vec<_>>());
        r["budget"]["max_output_bytes"] = json!(100000);
    });
    assert_eq!(r.termination, Stop::OutputLimit);
    assert_eq!(r.completed_replicates, 32);
    assert_eq!(r.replicates.len(), 32);
    assert!(r.statistics.is_none() && r.trajectories.is_empty());
}

fn cancel_dir(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("crn-ssa-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir(&dir).unwrap();
    dir
}

#[test]
fn existing_cancel_file_stops_every_worker_before_any_replicate_starts() {
    let dir = cancel_dir("existing");
    std::fs::write(dir.join("STOP"), b"").unwrap();
    let bytes = serde_json::to_vec(&model()).unwrap();
    let mut r = run(&bytes);
    r["budget"]["threads"] = json!(4);
    let mut control = Control::default();
    control.cancel_file = Some(dir.join("STOP"));
    let result = simulation::simulate_bytes(&bytes, &serde_json::to_vec(&r).unwrap(), &control)
        .unwrap()
        .result;
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(result.termination, Stop::Cancelled);
    assert_eq!(result.replicates.len(), 32);
    assert!(
        result
            .replicates
            .iter()
            .all(|r| !r.started && r.events == 0)
    );
}

#[test]
fn cancel_file_created_during_a_run_stops_all_workers() {
    let dir = cancel_dir("live");
    let mut m = model();
    m["reactions"][0]["rate"]["coefficient"] = json!(1000.0);
    let bytes = serde_json::to_vec(&m).unwrap();
    let mut r = run(&bytes);
    // About 10^9 events per replicate: only cancellation can end this run
    // soon. The wall-time limit only bounds the test if cancellation breaks.
    r["time"] = json!({"end": 1e6, "save_times": [0.0, 1e6]});
    r["budget"]["wall_time_s"] = json!(20.0);
    r["ensemble"]["replicates"] = json!(8);
    r["outputs"] = json!({"statistics": ["mean"], "keep_trajectories": 0});
    r["budget"]["threads"] = json!(4);
    r["budget"]["max_events_per_replicate"] = json!(10_000_000_000_u64);
    let mut control = Control::default();
    control.cancel_file = Some(dir.join("STOP"));
    let start = std::time::Instant::now();
    let result = std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(50));
            std::fs::write(dir.join("STOP"), b"").unwrap();
        });
        simulation::simulate_bytes(&bytes, &serde_json::to_vec(&r).unwrap(), &control)
            .unwrap()
            .result
    });
    let elapsed = start.elapsed();
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(result.termination, Stop::Cancelled);
    assert!(result.events > 0);
    assert!(
        result
            .replicates
            .iter()
            .all(|r| r.termination == Stop::Cancelled),
        "{:?}",
        result.replicates
    );
    assert!(elapsed.as_secs() < 15, "cancellation took {elapsed:?}");
}

#[test]
fn lost_channel_and_underflowed_wait_fail_at_the_sampler_boundary() {
    // 1 + 2^-53 rounds to 1, yet the tiny channel's share is exactly the
    // 53-bit step, so only the lost-increment check can reject it.
    let lost = sampling::total_propensity(&[1.0, sampling::UNIFORM_STEP]).unwrap_err();
    assert_eq!(lost.code, "CRN-E061");
    assert!(
        lost.message.contains("loses a positive channel"),
        "{lost:?}"
    );
    // -ln(1 - 2^-53) / f64::MAX underflows to zero.
    let wait = sampling::waiting_time(f64::MAX, 1.0 - sampling::UNIFORM_STEP).unwrap_err();
    assert_eq!(wait.code, "CRN-E062");
}

#[test]
fn retained_trajectories_are_exactly_the_first_requested_ids() {
    let r = simulate(&model(), |r| r["outputs"]["keep_trajectories"] = json!(2));
    assert_eq!(r.replicates.len(), 32);
    assert_eq!(
        r.trajectories.iter().map(|t| t.id).collect::<Vec<_>>(),
        [0, 1]
    );
}
