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
