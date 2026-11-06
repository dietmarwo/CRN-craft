//! Generated event loops for ten fixture models are checked in, so the
//! generated-loop path runs under `cargo test` without a separate build.
//! Regenerate a fixture with
//! `cargo run --release -p crn-craft-rules --example generate_rules_loop -- MODEL REQUEST`
//! when the generator changes; the first test names the stale fixture.
use crn_craft_rules::{
    Control, RulesLoop, execute, execute_with_rules_loop, generate_rules_loop, prepare,
};
use serde_json::{Value, json};

mod stop {
    include!("generated/mtdna_three_class_stop_loop.rs");
}
mod pairs {
    include!("generated/loop_pairs_loop.rs");
}
mod resolution {
    include!("generated/loop_resolution_loop.rs");
}
mod invalid {
    include!("generated/loop_invalid_rate_loop.rs");
}
mod edges_nan {
    include!("generated/loop_edges_nan_loop.rs");
}
mod edges_fraction {
    include!("generated/loop_edges_fraction_loop.rs");
}
mod count_range {
    include!("generated/loop_count_range_loop.rs");
}
mod overflow {
    include!("generated/loop_overflow_loop.rs");
}
mod empty {
    include!("generated/loop_empty_loop.rs");
}
mod operators {
    include!("generated/loop_stop_operators_loop.rs");
}

const STOP: &[u8] = include_bytes!("generated/mtdna_three_class_stop.rules.json");
const PAIRS: &[u8] = include_bytes!("generated/loop_pairs.rules.json");
const RESOLUTION: &[u8] = include_bytes!("generated/loop_resolution.rules.json");
const INVALID: &[u8] = include_bytes!("generated/loop_invalid_rate.rules.json");
const EDGES_NAN: &[u8] = include_bytes!("generated/loop_edges_nan.rules.json");
const EDGES_FRACTION: &[u8] = include_bytes!("generated/loop_edges_fraction.rules.json");
const COUNT_RANGE: &[u8] = include_bytes!("generated/loop_count_range.rules.json");
const OVERFLOW: &[u8] = include_bytes!("generated/loop_overflow.rules.json");
const EMPTY: &[u8] = include_bytes!("generated/loop_empty.rules.json");
const OPERATORS: &[u8] = include_bytes!("generated/loop_stop_operators.rules.json");

const PLAIN: &[u8] = include_bytes!("../../../models/rules/mtdna-three-class.rules.json");

const DEFAULT: (&str, &str) = ("crn-craft.rules-run/0.2", r#"{"kind": "static_direct"}"#);
/// Every waiting-time and tracker policy, including the rules-run/0.6
/// defaults (Ziggurat waits, sparse trackers where eligible) and the explicit
/// options that restore the earlier default.
const METHODS: [(&str, &str); 6] = [
    DEFAULT,
    (
        "crn-craft.rules-run/0.4",
        r#"{"kind": "static_direct", "waiting_time": "ziggurat"}"#,
    ),
    (
        "crn-craft.rules-run/0.4",
        r#"{"kind": "static_direct", "sparse_trackers": true}"#,
    ),
    (
        "crn-craft.rules-run/0.4",
        r#"{"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true}"#,
    ),
    ("crn-craft.rules-run/0.6", r#"{"kind": "static_direct"}"#),
    (
        "crn-craft.rules-run/0.6",
        r#"{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": false}"#,
    ),
];

struct Run {
    method: (&'static str, &'static str),
    end: f64,
    reports: Vec<f64>,
    replicates: u64,
    threads: u32,
    events: u64,
    words: u64,
    instances: u64,
    keep_instances: u64,
    trackers: Option<Value>,
    event_log: Option<Value>,
}

impl Run {
    fn new(end: f64, reports: Vec<f64>) -> Self {
        Self {
            method: DEFAULT,
            end,
            reports,
            replicates: 8,
            threads: 1,
            events: 1_000_000,
            words: 10_000_000,
            instances: 1000,
            keep_instances: 2,
            trackers: None,
            event_log: None,
        }
    }

    fn request(&self, model: &[u8]) -> Vec<u8> {
        let (schema, method) = self.method;
        let mut request = json!({
            "schema": schema,
            "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(model)},
            "time": {"end": self.end, "report_times": self.reports},
            "ensemble": {"replicates": self.replicates, "root_seed": 23},
            "outputs": {"statistics": ["mean", "variance"], "keep_paths": self.replicates,
                        "keep_instances": self.keep_instances},
            "budget": {"threads": self.threads, "max_events_per_replicate": self.events,
                       "max_sampling_steps_per_replicate": self.words,
                       "max_instances": self.instances, "max_output_bytes": 1u64 << 28},
            "method": serde_json::from_str::<Value>(method).unwrap(),
        });
        if let Some(trackers) = &self.trackers {
            request["trackers"] = trackers.clone();
        }
        if let Some(log) = &self.event_log {
            request["event_log"] = log.clone();
        }
        serde_json::to_vec(&request).unwrap()
    }
}

/// Runs the request with the interpreter and the generated loop, asserts that
/// the results are byte-identical and returns them.
fn replay(model: &[u8], run: &Run, kernel: &impl RulesLoop) -> String {
    let prepared = prepare(model, &run.request(model)).unwrap();
    let interpreted = serde_json::to_string(&execute(&prepared, &Control::default()).unwrap());
    let generated = serde_json::to_string(
        &execute_with_rules_loop(&prepared, &Control::default(), kernel).unwrap(),
    );
    let (interpreted, generated) = (interpreted.unwrap(), generated.unwrap());
    assert!(
        interpreted == generated,
        "{} end {} differs from the interpreter",
        run.method.1,
        run.end
    );
    interpreted
}

fn trackers() -> Option<Value> {
    Some(json!([{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]))
}

#[test]
fn generator_output_matches_the_checked_in_loops() {
    for (model, golden, name) in [
        (
            STOP,
            include_str!("generated/mtdna_three_class_stop_loop.rs"),
            "mtdna_three_class_stop",
        ),
        (
            PAIRS,
            include_str!("generated/loop_pairs_loop.rs"),
            "loop_pairs",
        ),
        (
            RESOLUTION,
            include_str!("generated/loop_resolution_loop.rs"),
            "loop_resolution",
        ),
        (
            INVALID,
            include_str!("generated/loop_invalid_rate_loop.rs"),
            "loop_invalid_rate",
        ),
        (
            EDGES_NAN,
            include_str!("generated/loop_edges_nan_loop.rs"),
            "loop_edges_nan",
        ),
        (
            EDGES_FRACTION,
            include_str!("generated/loop_edges_fraction_loop.rs"),
            "loop_edges_fraction",
        ),
        (
            COUNT_RANGE,
            include_str!("generated/loop_count_range_loop.rs"),
            "loop_count_range",
        ),
        (
            OVERFLOW,
            include_str!("generated/loop_overflow_loop.rs"),
            "loop_overflow",
        ),
        (
            EMPTY,
            include_str!("generated/loop_empty_loop.rs"),
            "loop_empty",
        ),
        (
            OPERATORS,
            include_str!("generated/loop_stop_operators_loop.rs"),
            "loop_stop_operators",
        ),
    ] {
        let prepared = prepare(model, &Run::new(10.0, vec![0.0, 10.0]).request(model)).unwrap();
        let source = generate_rules_loop(&prepared).unwrap();
        assert!(
            source == golden,
            "generator output changed; regenerate {name}_loop.rs"
        );
    }
}

#[test]
fn the_mtdna_loop_replays_every_static_direct_policy() {
    let grid = vec![0.0, 250.0, 500.0, 1000.0, 2000.0];
    for method in METHODS {
        for (tracked, threads) in [(None, 1), (trackers(), 1), (trackers(), 3)] {
            let run = Run {
                method,
                threads,
                trackers: tracked,
                ..Run::new(2000.0, grid.clone())
            };
            replay(STOP, &run, &stop::GeneratedRulesLoop);
        }
    }
}

#[test]
fn the_mtdna_loop_stops_like_the_interpreter() {
    // At this horizon the 90% takeover stop fires in some cells.
    let run = Run {
        method: METHODS[3],
        replicates: 24,
        trackers: trackers(),
        ..Run::new(25_920.0, vec![0.0, 8640.0, 17_280.0, 25_920.0])
    };
    let result = replay(STOP, &run, &stop::GeneratedRulesLoop);
    assert!(result.contains("\"stopped\""), "no replicate stopped");
}

#[test]
fn budgets_and_short_horizons_end_paths_identically() {
    let grid = vec![0.0, 1.0, 50.0, 100.0];
    for method in METHODS {
        let runs = [
            // The first wait lies beyond the horizon.
            Run {
                method,
                ..Run::new(1e-9, vec![0.0, 1e-9])
            },
            Run {
                method,
                events: 25,
                ..Run::new(100.0, grid.clone())
            },
            Run {
                method,
                words: 41,
                ..Run::new(100.0, grid.clone())
            },
            // The first mutation creates a second class.
            Run {
                method,
                instances: 1,
                ..Run::new(5000.0, vec![0.0, 5000.0])
            },
        ];
        for run in &runs {
            replay(STOP, run, &stop::GeneratedRulesLoop);
        }
    }
    let events = Run {
        events: 25,
        ..Run::new(100.0, grid.clone())
    };
    assert!(replay(STOP, &events, &stop::GeneratedRulesLoop).contains("CRN-R040"));
    let words = Run {
        words: 41,
        ..Run::new(100.0, grid)
    };
    assert!(replay(STOP, &words, &stop::GeneratedRulesLoop).contains("CRN-R041"));
    let instances = Run {
        instances: 1,
        ..Run::new(5000.0, vec![0.0, 5000.0])
    };
    assert!(replay(STOP, &instances, &stop::GeneratedRulesLoop).contains("CRN-R042"));
}

#[test]
fn pairs_extinction_and_a_general_stop_replay() {
    for method in METHODS {
        for (end, instances) in [(40.0, 1000), (40.0, 1)] {
            let run = Run {
                method,
                replicates: 32,
                instances,
                ..Run::new(end, vec![0.0, 1.0, 10.0, 40.0])
            };
            replay(PAIRS, &run, &pairs::GeneratedRulesLoop);
        }
    }
    let result = replay(
        PAIRS,
        &Run {
            replicates: 32,
            ..Run::new(40.0, vec![0.0, 40.0])
        },
        &pairs::GeneratedRulesLoop,
    );
    assert!(
        result.contains("\"stopped\""),
        "the stop predicate never fired"
    );
    assert!(
        result.contains("\"completed\""),
        "no path reached the horizon"
    );
}

#[test]
fn an_unresolvable_channel_fails_identically() {
    for method in METHODS {
        let run = Run {
            method,
            replicates: 16,
            ..Run::new(200.0, vec![0.0, 100.0, 200.0])
        };
        let result = replay(RESOLUTION, &run, &resolution::GeneratedRulesLoop);
        assert!(result.contains("CRN-R031"), "no replicate lost resolution");
    }
}

#[test]
fn an_invalid_rate_fails_identically() {
    for method in METHODS {
        let run = Run {
            method,
            replicates: 8,
            ..Run::new(100.0, vec![0.0, 100.0])
        };
        let result = replay(INVALID, &run, &invalid::GeneratedRulesLoop);
        assert!(result.contains("invalid fixed-channel propensity"));
    }
}

#[test]
fn a_run_cancelled_before_it_starts_is_reported_identically() {
    let prepared = prepare(STOP, &Run::new(100.0, vec![0.0, 100.0]).request(STOP)).unwrap();
    let control = Control::default();
    control.cancel();
    let interpreted = execute(&prepared, &control);
    let generated = execute_with_rules_loop(&prepared, &control, &stop::GeneratedRulesLoop);
    assert_eq!(
        serde_json::to_string(&interpreted.unwrap()).unwrap(),
        serde_json::to_string(&generated.unwrap()).unwrap()
    );
}

#[test]
fn a_loop_for_other_model_bytes_is_rejected() {
    // Parameters and the stop predicate are compiled into the loop.
    let prepared = prepare(PLAIN, &Run::new(10.0, vec![0.0, 10.0]).request(PLAIN)).unwrap();
    let Err(error) =
        execute_with_rules_loop(&prepared, &Control::default(), &stop::GeneratedRulesLoop)
    else {
        panic!("stale loop accepted");
    };
    assert_eq!(error.code, "CRN-R020");
    assert_eq!(error.path, "/model");
}

#[test]
fn other_methods_and_rates_are_refused() {
    let run = Run {
        method: ("crn-craft.rules-run/0.2", r#"{"kind": "next_reaction"}"#),
        ..Run::new(10.0, vec![0.0, 10.0])
    };
    let prepared = prepare(STOP, &run.request(STOP)).unwrap();
    let error = generate_rules_loop(&prepared).unwrap_err();
    assert_eq!(error.path, "/method");
    let Err(error) =
        execute_with_rules_loop(&prepared, &Control::default(), &stop::GeneratedRulesLoop)
    else {
        panic!("next reaction accepted");
    };
    assert_eq!(error.path, "/method");

    let mut model: Value = serde_json::from_slice(PAIRS).unwrap();
    model["rules"][0]["rate"] = json!("d * (1 + nB)");
    let model = serde_json::to_vec(&model).unwrap();
    let prepared = prepare(&model, &Run::new(10.0, vec![0.0, 10.0]).request(&model)).unwrap();
    assert!(generate_rules_loop(&prepared).is_err());
}

/// Trackers, first-passage trackers and event logs, over the edge fixtures'
/// observables.
fn observations() -> [(Option<Value>, Option<Value>); 4] {
    [
        (None, None),
        (
            Some(json!([{"observable": "nC"}, {"observable": "a2", "keep_episodes": 2}])),
            None,
        ),
        (
            Some(json!([{"observable": "neg", "first_passage": "nB < 5", "keep_episodes": 1}])),
            None,
        ),
        (
            Some(json!([{"observable": "nC"}])),
            Some(json!({"max_entries": 5, "rules": ["cross", "pairB"]})),
        ),
    ]
}

/// Every policy, observation set and budget, with reports close together.
fn edge_matrix(model: &[u8], kernel: &impl RulesLoop, ends: [f64; 2]) -> Vec<String> {
    let mut results = Vec::new();
    for method in METHODS {
        for (trackers, event_log) in observations() {
            for end in ends {
                let reports = vec![
                    0.0,
                    end * 1e-6,
                    end * 0.01,
                    end * 0.010_000_1,
                    end * 0.5,
                    end,
                ];
                for (events, words, instances) in [
                    (1_000_000, 10_000_000, 1000),
                    (37, 10_000_000, 1000),
                    (1_000_000, 53, 1000),
                    (1_000_000, 10_000_000, 3),
                ] {
                    let run = Run {
                        method,
                        replicates: 12,
                        threads: 2,
                        events,
                        words,
                        instances,
                        keep_instances: 12,
                        trackers: trackers.clone(),
                        event_log: event_log.clone(),
                        ..Run::new(end, reports.clone())
                    };
                    results.push(replay(model, &run, kernel));
                }
            }
        }
    }
    results
}

#[test]
fn shared_vacancy_rates_logs_and_an_undefined_stop_replay() {
    let results = edge_matrix(EDGES_NAN, &edges_nan::GeneratedRulesLoop, [5.0, 400.0]);
    assert!(results.iter().any(|r| r.contains("CRN-R040")));
    assert!(results.iter().any(|r| r.contains("undefined (NaN)")));
}

#[test]
fn shared_vacancy_rates_and_an_emitted_fraction_stop_replay() {
    let results = edge_matrix(
        EDGES_FRACTION,
        &edges_fraction::GeneratedRulesLoop,
        [5.0, 400.0],
    );
    assert!(results.iter().any(|r| r.contains("\"stopped\"")));
}

#[test]
fn count_range_and_projection_overflow_fail_identically() {
    let run = Run::new(10.0, vec![0.0, 1.0, 10.0]);
    let result = replay(COUNT_RANGE, &run, &count_range::GeneratedRulesLoop);
    assert!(result.contains("outside 0..2^53"));
    for method in METHODS {
        let run = Run {
            method,
            ..Run::new(100.0, vec![0.0, 1.0, 100.0])
        };
        let result = replay(OVERFLOW, &run, &overflow::GeneratedRulesLoop);
        assert!(result.contains("observable overflows"));
    }
}

#[test]
fn a_network_without_channels_replays() {
    for method in METHODS {
        let run = Run {
            method,
            ..Run::new(10.0, vec![0.0, 10.0])
        };
        replay(EMPTY, &run, &empty::GeneratedRulesLoop);
    }
}

#[test]
fn cancellation_is_polled_inside_the_loop() {
    // The engine polls every 1024 events, so a run cancelled while a
    // replicate is in the generated loop ends that replicate at a multiple
    // of 1024 events.
    let run = Run {
        replicates: 512,
        ..Run::new(25_920.0, vec![0.0, 25_920.0])
    };
    let prepared = prepare(STOP, &run.request(STOP)).unwrap();
    let control = Control::default();
    let result = std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(200));
            control.cancel();
        });
        execute_with_rules_loop(&prepared, &control, &stop::GeneratedRulesLoop).unwrap()
    });
    let result = serde_json::to_value(&result).unwrap();
    let interrupted = result["replicates"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|r| r["termination"] == "cancelled" && r["events"].as_u64().unwrap() > 0)
        .map(|r| r["events"].as_u64().unwrap())
        .collect::<Vec<_>>();
    assert!(
        !interrupted.is_empty(),
        "no replicate was cancelled mid-run"
    );
    assert!(interrupted.iter().all(|events| events % 1024 == 0));
}

#[test]
fn densely_coupled_networks_are_refused_before_they_reach_the_compiler() {
    // 128 instances share one vacancy base, so every event updates every
    // replication rate and the emitted source grows with the square.
    let model = json!({
        "schema": "crn-craft.rules-model/0.1", "name": "dense",
        "parameters": {"K": 100_000, "r": 1.0, "d": 0.1},
        "types": [{"name": "A", "attributes": ["s"]}],
        "initial": (0..128)
            .map(|i| json!({"type": "A", "attributes": {"s": i}, "count": 10}))
            .collect::<Vec<_>>(),
        "observables": {"N": {"type": "A"}, "D": {"type": "A", "weight": "1 + x.s"}},
        "rules": [
            {"name": "rep", "reactants": [{"label": "x", "type": "A"}],
             "rate": "max(0, K - N) / D * r",
             "products": [{"type": "A", "copy": "x"}, {"type": "A", "copy": "x"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "A"}], "rate": "d",
             "products": []}
        ]
    });
    let model = serde_json::to_vec(&model).unwrap();
    let prepared = prepare(&model, &Run::new(1.0, vec![0.0, 1.0]).request(&model)).unwrap();
    let error = generate_rules_loop(&prepared).unwrap_err();
    assert_eq!(error.path, "/method");
    assert!(error.message.contains("1 MiB"), "{}", error.message);
}

#[test]
fn every_stop_operator_is_emitted_like_the_bytecode() {
    for method in METHODS {
        let run = Run {
            method,
            replicates: 32,
            ..Run::new(40.0, vec![0.0, 1.0, 10.0, 40.0])
        };
        let result = replay(OPERATORS, &run, &operators::GeneratedRulesLoop);
        assert!(result.contains("\"stopped\""));
    }
}
