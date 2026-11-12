//! M8a: boundaries that only the unified engine's own tests can pin.
//!
//! The mutation gate (`validation/mutation_check.py`) re-pointed from the old
//! network engine to the unified engine's code. These tests cover the
//! behavior its mutants change, which the old engine's tests covered and the
//! unified tests did not: event-cap boundaries, the numerical-resolution
//! refusals, polling during a long run, and the sampler identity of the
//! approximate methods.
use crn_craft::rules::{Control, RulesResult, Termination, execute, prepare};
use serde_json::{Value, json};
use std::sync::Arc;
use std::time::Duration;

/// `A -> B` at rate `k`, `n` molecules of A at the start. Exactly `n`
/// events fire, whatever the stream, once the horizon is far enough.
fn decay(n: u64, k: f64) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "decay",
        "parameters": {"k": k},
        "types": [{"name": "A"}, {"name": "B"}],
        "initial": [{"type": "A", "count": n}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [{"name": "decay",
                   "reactants": [{"label": "x", "type": "A"}],
                   "rate": "k",
                   "products": [{"type": "B"}]}]
    }))
    .unwrap()
}

/// `A -> B` at rate 1 and `A -> C` at rate `tiny`: the second channel is far
/// below 2^-53 of the total.
fn unresolvable(tiny: f64, copies: usize) -> Vec<u8> {
    let mut rules = vec![json!({"name": "main",
        "reactants": [{"label": "x", "type": "A"}], "rate": "one",
        "products": [{"type": "B"}]})];
    for i in 0..copies {
        rules.push(json!({"name": format!("rare{i}"),
            "reactants": [{"label": "x", "type": "A"}], "rate": "tiny",
            "products": [{"type": "C"}]}));
    }
    serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "unresolvable",
        "parameters": {"one": 1.0, "tiny": tiny},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C"}],
        "initial": [{"type": "A", "count": 1000}],
        "observables": {"A": {"type": "A"}},
        "rules": rules
    }))
    .unwrap()
}

fn request(model: &[u8], method: Value, end: f64, max_events: u64) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "method": method,
        "ensemble": {"replicates": 2, "root_seed": 5},
        "time": {"end": end, "report_times": [0.0, end]},
        "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"max_events_per_replicate": max_events,
                   "max_sampling_steps_per_replicate": 100_000_000,
                   "max_instances": 8, "threads": 1, "max_output_bytes": 1u64 << 24}
    }))
    .unwrap()
}

fn run(model: &[u8], request: &[u8], control: &Control) -> RulesResult {
    let prepared = prepare(model, request).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, control).unwrap_or_else(|d| panic!("{d:?}"))
}

fn methods() -> Vec<(&'static str, Value)> {
    let static_direct = |selector: &str, waiting_time: &str| json!({"kind": "static_direct", "selector": selector, "waiting_time": waiting_time});
    vec![
        ("default", json!({"kind": "static_direct"})),
        ("linear/inverse", static_direct("linear", "inverse")),
        ("linear/ziggurat", static_direct("linear", "ziggurat")),
        (
            "composition_rejection",
            static_direct("composition_rejection", "inverse"),
        ),
        ("sorting", static_direct("sorting", "inverse")),
        ("rssa_cr", static_direct("rssa_cr", "inverse")),
        ("next_reaction", json!({"kind": "next_reaction"})),
        ("tree", json!({"kind": "tree"})),
    ]
}

#[test]
fn an_event_cap_equal_to_the_events_needed_completes_and_one_less_does_not() {
    let model = decay(5, 1.0);
    for (name, method) in methods() {
        let exact = run(
            &model,
            &request(&model, method.clone(), 1e6, 5),
            &Control::default(),
        );
        assert_eq!(exact.termination, Termination::Completed, "{name}");
        assert!(exact.replicates.iter().all(|r| r.events == 5), "{name}");
        let short = run(
            &model,
            &request(&model, method, 1e6, 4),
            &Control::default(),
        );
        assert_eq!(short.termination, Termination::EventLimit, "{name}");
        assert!(short.statistics.is_none(), "{name}");
    }
}

#[test]
fn a_tau_leap_that_would_pass_the_event_cap_stops_at_the_cap() {
    let model = decay(100_000, 1.0);
    let method = json!({"kind": "tau", "epsilon": 0.03});
    let capped = run(
        &model,
        &request(&model, method.clone(), 10.0, 10),
        &Control::default(),
    );
    assert_eq!(capped.termination, Termination::EventLimit);
    assert!(capped.replicates.iter().all(|r| r.events <= 10));
    let free = run(
        &model,
        &request(&model, method, 10.0, 10_000_000),
        &Control::default(),
    );
    assert_eq!(free.termination, Termination::Completed);
}

#[test]
fn the_sampler_identity_says_which_methods_are_approximate() {
    let model = decay(50, 1.0);
    let sampler = |method: Value| {
        run(
            &model,
            &request(&model, method, 5.0, 1000),
            &Control::default(),
        )
        .sampler
    };
    for exact in [json!({"kind": "static_direct"}), json!({"kind": "tree"})] {
        let name = sampler(exact);
        assert!(!name.contains("approximate"), "{name}");
    }
    let tau = sampler(json!({"kind": "tau", "epsilon": 0.03}));
    assert!(
        tau.contains("approximate") && tau.contains("poisson"),
        "{tau}"
    );
}

#[test]
fn a_channel_below_the_resolution_of_the_total_is_refused_not_misselected() {
    // Static channels: a rate below 2^-53 of the total cannot be selected by
    // a 53-bit uniform. The linear selector scans below 32 channels and
    // fuses the total and the smallest rate from 32 on; both refuse.
    for copies in [1, 40] {
        let model = unresolvable(1e-20, copies);
        let method = json!({"kind": "static_direct", "selector": "linear",
                            "waiting_time": "inverse"});
        let result = run(
            &model,
            &request(&model, method, 1.0, 1_000_000),
            &Control::default(),
        );
        assert_ne!(result.termination, Termination::Completed, "{copies}");
        let code = result
            .replicates
            .iter()
            .find_map(|r| r.diagnostic.as_ref())
            .map(|d| d.code.clone());
        assert_eq!(code.as_deref(), Some("CRN-R031"), "{copies} rare channels");
    }
    // The tree method resolves rules and instances separately: a rule whose
    // weight is below 2^-53 of its total is refused at the rule level.
    let model = unresolvable(1e-20, 1);
    let result = run(
        &model,
        &request(&model, json!({"kind": "tree"}), 1.0, 1_000_000),
        &Control::default(),
    );
    assert_ne!(result.termination, Termination::Completed);
    assert!(
        result
            .replicates
            .iter()
            .filter_map(|r| r.diagnostic.as_ref())
            .any(|d| d.code == "CRN-R031" && d.path.starts_with("/rules/")),
        "{:?}",
        result.replicates.first().map(|r| &r.diagnostic)
    );
}

#[test]
fn a_cancel_during_a_long_run_stops_it_within_the_poll_interval() {
    // Millions of events, so a run that never polls would finish only at the
    // event cap; a run that polls every 1,024 events stops far sooner.
    let model = decay(4_000_000, 1.0);
    for (name, method) in methods() {
        let control = Arc::new(Control::default());
        let flag = Arc::clone(&control);
        let stopper = std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(20));
            flag.cancel();
        });
        let result = run(&model, &request(&model, method, 1e9, 4_000_000), &control);
        stopper.join().unwrap();
        assert_eq!(result.termination, Termination::Cancelled, "{name}");
        assert!(
            result.replicates.iter().all(|r| r.events < 4_000_000),
            "{name}: {:?}",
            result
                .replicates
                .iter()
                .map(|r| r.events)
                .collect::<Vec<_>>()
        );
    }
}

/// M8b: the engine's stream against the independent scalar ChaCha8 of
/// `validation/reference_rng.py` (nine streams, 40 words each, across
/// blocks and stream ids). The network engine's `tests/ssa.rs` pinned its
/// own stream this way; that stream left with the network engine.
#[test]
fn the_replicate_stream_matches_independent_known_answers() {
    use rand_core::RngCore;
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/rng-vectors-rules.json")).unwrap();
    let streams = fixture["streams"].as_array().unwrap();
    assert_eq!(streams.len(), 9);
    for row in streams {
        let mut rng = crn_craft::rules::stream(
            row["root"].as_u64().unwrap(),
            row["replicate"].as_u64().unwrap(),
        );
        for word in row["words"].as_array().unwrap() {
            assert_eq!(rng.next_u64(), word.as_u64().unwrap());
        }
    }
}

/// Immigration-death (`0 -> X` at 4, `X -> 0` at 0.2 per molecule) and a
/// Static Direct request with `replicates`, `report_times` and `keep_paths`.
fn immigration(
    replicates: u64,
    report_times: &[f64],
    keep: u64,
    max_events: u64,
) -> (Vec<u8>, Vec<u8>) {
    let model = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "immigration_death",
        "parameters": {"birth": 4.0, "death": 0.2},
        "types": [{"name": "X"}], "initial": [],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {"name": "birth", "reactants": [], "rate": "birth", "products": [{"type": "X"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "death", "products": []}
        ]
    }))
    .unwrap();
    let end = *report_times.last().unwrap();
    let request = serde_json::to_vec(&json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "m.json", "sha256": crn_craft::sha256(&model)},
        "method": {"kind": "static_direct"},
        "ensemble": {"replicates": replicates, "root_seed": 77},
        "time": {"end": end, "report_times": report_times},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": keep, "keep_instances": 0},
        "budget": {"max_events_per_replicate": max_events,
                   "max_sampling_steps_per_replicate": 100_000_000,
                   "max_instances": 1000, "threads": 2, "max_output_bytes": 1u64 << 24}
    }))
    .unwrap();
    (model, request)
}

/// M8b, from the network engine's `tests/ssa.rs`: report times draw no
/// random words, so a denser grid leaves every event and every state at the
/// common times unchanged.
#[test]
fn a_denser_report_grid_changes_no_event_or_common_state() {
    let sparse_times = [0.0, 0.5, 1.0, 2.0];
    let dense_times = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0];
    let (model, sparse) = immigration(32, &sparse_times, 32, 100_000);
    let (_, dense) = immigration(32, &dense_times, 32, 100_000);
    let sparse = run(&model, &sparse, &Control::default());
    let dense = run(&model, &dense, &Control::default());
    let events = |r: &RulesResult| r.replicates.iter().map(|x| x.events).collect::<Vec<_>>();
    assert_eq!(events(&sparse), events(&dense));
    for (a, b) in sparse.paths.iter().zip(&dense.paths) {
        assert_eq!(a.id, b.id);
        for (i, j) in [0, 3, 5, 7].into_iter().enumerate() {
            assert_eq!(a.values[i], b.values[j], "replicate {}", a.id);
        }
    }
}

/// M8b, from `tests/ssa.rs`: a rule without reactants or products changes
/// no count, yet every firing is an event and counts against the cap.
#[test]
fn a_rule_without_reactants_or_products_counts_events_and_meets_the_cap() {
    let model = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "null",
        "parameters": {"k": 1e6}, "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": 3}], "observables": {"A": {"type": "A"}},
        "rules": [{"name": "null", "reactants": [], "rate": "k", "products": []}]
    }))
    .unwrap();
    let request = request(&model, json!({"kind": "static_direct"}), 1.0, 2);
    let result = run(&model, &request, &Control::default());
    assert_eq!(result.termination, Termination::EventLimit);
    assert!(result.replicates.iter().all(|r| r.events == 2));
    assert!(result.statistics.is_none());
}

/// M8b, from `tests/ssa.rs`: the kept paths are the first replicates, and
/// the statistics are the plain moments of the completed replicates (here
/// all kept); one replicate has no variance and no standard error.
#[test]
fn kept_paths_are_the_first_replicates_and_carry_the_statistics() {
    let times = [0.0, 1.0, 5.0];
    let (model, first) = immigration(32, &times, 2, 100_000);
    let result = run(&model, &first, &Control::default());
    assert_eq!(result.replicates.len(), 32);
    assert_eq!(
        result.paths.iter().map(|p| p.id).collect::<Vec<_>>(),
        [0, 1]
    );
    let (_, all) = immigration(32, &times, 32, 100_000);
    let result = run(&model, &all, &Control::default());
    let statistics = result.statistics.as_ref().unwrap();
    let (mean, variance, error) = (
        statistics.mean.as_ref().unwrap(),
        statistics.sample_variance.as_ref().unwrap(),
        statistics.mean_standard_error.as_ref().unwrap(),
    );
    for t in 0..times.len() {
        let x: Vec<f64> = result.paths.iter().map(|p| p.values[t][0]).collect();
        let m = x.iter().sum::<f64>() / 32.0;
        let v = x.iter().map(|y| (y - m) * (y - m)).sum::<f64>() / 31.0;
        assert!((mean[t][0] - m).abs() <= 1e-12 * m.max(1.0));
        assert!((variance[t][0] - v).abs() <= 1e-12 * v.max(1.0));
        assert!((error[t][0] - (v / 32.0).sqrt()).abs() <= 1e-12);
    }
    let (_, single) = immigration(1, &times, 1, 100_000);
    let result = run(&model, &single, &Control::default());
    let statistics = result.statistics.unwrap();
    assert!(statistics.mean.is_some());
    assert!(statistics.sample_variance.is_none() && statistics.mean_standard_error.is_none());
}

/// M8b, from `tests/ssa.rs`: when some replicates complete and others stop
/// at the event cap, no conditional statistics of the completed ones are
/// reported.
#[test]
fn a_partly_completed_ensemble_reports_no_statistics() {
    let (model, request) = immigration(64, &[0.0, 0.3], 0, 1);
    let result = run(&model, &request, &Control::default());
    assert!(
        result.completed_replicates > 0 && result.completed_replicates < 64,
        "{}",
        result.completed_replicates
    );
    assert_eq!(result.termination, Termination::EventLimit);
    assert!(result.statistics.is_none());
}

/// M8b, from `tests/ssa.rs`: a run cancelled before it starts, or past its
/// wall-time limit, accounts for every replicate, and none of them fires an
/// event or draws a word. Unlike the network engine, which started none, the
/// unified engine marks the first replicate of each worker as started (it
/// stops at its first poll), so at most `threads` replicates are started.
#[test]
fn unstarted_replicates_are_accounted_for() {
    let (model, request) = immigration(32, &[0.0, 1.0], 0, 100_000);
    let idle = |result: &RulesResult| {
        result.replicates.len() == 32
            && result
                .replicates
                .iter()
                .all(|r| r.events == 0 && r.rng_words == 0 && r.termination == result.termination)
            && result.replicates.iter().filter(|r| r.started).count() <= 2
            && result.completed_replicates == 0
            && result.statistics.is_none()
    };
    let control = Control::default();
    control.cancel();
    let cancelled = run(&model, &request, &control);
    assert_eq!(cancelled.termination, Termination::Cancelled);
    assert!(idle(&cancelled), "{:?}", cancelled.replicates);
    let mut expired: Value = serde_json::from_slice(&request).unwrap();
    expired["budget"]["wall_time_s"] = json!(1e-300);
    let expired = run(
        &model,
        &serde_json::to_vec(&expired).unwrap(),
        &Control::default(),
    );
    assert_eq!(expired.termination, Termination::WallTimeLimit);
    assert!(idle(&expired), "{:?}", expired.replicates);
}

/// M8b, from the network engine's `tests/ssa_alternatives.rs`: an absorbing
/// state (here a pair rule with one molecule) ends the path at the horizon
/// without a random word, also where an RSSA-CR box has enabled neighbors.
#[test]
fn absorbing_states_draw_no_random_words() {
    let model = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "absorbed",
        "parameters": {"k": 1.0}, "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 1}], "observables": {"X": {"type": "X"}},
        "rules": [{"name": "pair", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "X"}],
                   "rate": "k", "products": []}]
    }))
    .unwrap();
    for (name, method) in methods() {
        let mut document: Value =
            serde_json::from_slice(&request(&model, method, 10.0, 100)).unwrap();
        document["budget"]["max_sampling_steps_per_replicate"] = json!(1);
        let result = run(
            &model,
            &serde_json::to_vec(&document).unwrap(),
            &Control::default(),
        );
        assert_eq!(result.termination, Termination::Completed, "{name}");
        assert!(
            result
                .replicates
                .iter()
                .all(|r| r.events == 0 && r.rng_words == 0),
            "{name}: {:?}",
            result.replicates
        );
    }
}

/// M8b, from `tests/ssa_alternatives.rs`: every exact method replays the
/// same result whatever the thread count; only the request hash (it covers
/// `threads`) differs.
#[test]
fn every_method_replays_across_thread_counts() {
    let model = decay(40, 0.5);
    for (name, method) in methods() {
        let results: Vec<Value> = [1, 3]
            .into_iter()
            .map(|threads| {
                let mut document: Value =
                    serde_json::from_slice(&request(&model, method.clone(), 5.0, 1000)).unwrap();
                document["ensemble"]["replicates"] = json!(16);
                document["outputs"]["keep_paths"] = json!(16);
                document["budget"]["threads"] = json!(threads);
                let result = run(
                    &model,
                    &serde_json::to_vec(&document).unwrap(),
                    &Control::default(),
                );
                let mut value = serde_json::to_value(result).unwrap();
                value.as_object_mut().unwrap().remove("run_sha256");
                value
            })
            .collect();
        assert_eq!(results[0], results[1], "{name}");
    }
}

/// M8b, from the network engine's `tests/tau.rs`: tau leaping never passes
/// the last reachable count of a critical homodimer (`2X -> X` from three
/// molecules ends at one), stays non-negative, and replays across thread
/// counts.
#[test]
fn tau_never_passes_the_last_reachable_count_of_a_critical_homodimer() {
    let model = serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0", "name": "dimer",
        "parameters": {"k": 100.0}, "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 3}], "observables": {"X": {"type": "X"}},
        "rules": [{"name": "critical", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "X"}],
                   "rate": "k", "products": [{"type": "X"}]}]
    }))
    .unwrap();
    let results: Vec<Value> = [1, 3]
        .into_iter()
        .map(|threads| {
            let mut document: Value = serde_json::from_slice(&request(
                &model,
                json!({"kind": "tau", "epsilon": 0.03}),
                10.0,
                1000,
            ))
            .unwrap();
            document["ensemble"]["replicates"] = json!(64);
            document["outputs"]["keep_paths"] = json!(64);
            document["budget"]["threads"] = json!(threads);
            let result = run(
                &model,
                &serde_json::to_vec(&document).unwrap(),
                &Control::default(),
            );
            assert_eq!(result.termination, Termination::Completed);
            assert!(result.paths.iter().all(|p| p.values[1] == [1.0]));
            assert!(
                result
                    .paths
                    .iter()
                    .all(|p| p.values.iter().all(|row| row[0] >= 1.0))
            );
            let mut value = serde_json::to_value(result).unwrap();
            value.as_object_mut().unwrap().remove("run_sha256");
            value
        })
        .collect();
    assert_eq!(results[0], results[1]);
}
