//! M4f: Static Direct accepts a positive wait that leaves the rounded clock
//! unchanged under `crn-craft.run/1.0`.
//!
//! At a late time t, a positive, finite waiting time below half an ulp of t
//! gives t + wait == t. Every `rules-run` version fails the replicate there
//! (`CRN-R030`, "clock failed to advance"), as 0.6.0 did. Under
//! `crn-craft.run/1.0` every selector processes the event (or the RSSA-CR
//! trial) at the same representable time, and fails with `CRN-R030` only
//! when the clock has not advanced for `CLOCK_STALL_LIMIT` (1,000,000)
//! consecutive steps, which keeps runaway rates detected.
use crn_craft::rules::{Control, RulesResult, Termination, execute, prepare};
use serde_json::{Value, json};

const SELECTORS: [&str; 4] = ["linear", "composition_rejection", "sorting", "rssa_cr"];
const LIMIT: u64 = 1_000_000;

/// G fires once (rate 1) and leaves `burst` molecules of A (at most 64
/// products), which then react at 1e20 per molecule: every wait after the
/// burst is below half an ulp of the burst time, so the clock stays at that
/// time. No other channel is positive then, so no rate is lost in the total.
/// Every rate reads two lowered species (with weight 0; without `flip` one B
/// molecule is present, so B is lowered), so RSSA-CR pins them and its
/// bounds are the exact values: no stale positive bound of G remains.
fn model(burst: usize, flip: bool) -> Vec<u8> {
    let products: Vec<Value> = (0..burst).map(|_| json!({"type": "A"})).collect();
    let mut rules = vec![
        json!({"name": "burst", "reactants": [{"label": "g", "type": "G"}],
                                "rate": "1 + 0 * (G + B)", "products": products}),
    ];
    if flip {
        // A and B turn into each other forever: a runaway rate.
        rules.push(
            json!({"name": "flip", "reactants": [{"label": "a", "type": "A"}],
                          "rate": "1e20 + 0 * (A + B)", "products": [{"type": "B"}]}),
        );
        rules.push(
            json!({"name": "flop", "reactants": [{"label": "b", "type": "B"}],
                          "rate": "1e20 + 0 * (A + B)", "products": [{"type": "A"}]}),
        );
    } else {
        rules.push(
            json!({"name": "decay", "reactants": [{"label": "a", "type": "A"}],
                          "rate": "1e20 + 0 * (A + B)", "products": []}),
        );
    }
    let mut initial = vec![json!({"type": "G", "count": 1})];
    if !flip {
        initial.push(json!({"type": "B", "count": 1}));
    }
    serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "clock",
        "types": [{"name": "G"}, {"name": "A"}, {"name": "B"}],
        "initial": initial,
        "observables": {"G": {"type": "G"}, "A": {"type": "A"}, "B": {"type": "B"}},
        "rules": rules
    }))
    .unwrap()
}

fn request(model: &[u8], schema: &str, method: &Value, replicates: u64) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": schema,
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 1000.0, "report_times": [0.0, 1000.0]},
        "ensemble": {"replicates": replicates, "root_seed": 5},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates,
                    "keep_instances": 0},
        "budget": {"max_events_per_replicate": 10_000_000,
                   "max_sampling_steps_per_replicate": 1_000_000_000,
                   "max_instances": 3, "threads": 1, "max_output_bytes": 1u64 << 26},
        "method": method
    }))
    .unwrap()
}

fn run(model: &[u8], request: &[u8]) -> RulesResult {
    let prepared = prepare(model, request).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap()
}

fn final_value(result: &RulesResult, replicate: usize, observable: &str) -> f64 {
    let column = result
        .observables
        .iter()
        .position(|o| o == observable)
        .unwrap();
    *result.paths[replicate]
        .values
        .last()
        .unwrap()
        .get(column)
        .unwrap()
}

#[test]
fn run_1_0_processes_events_at_an_unchanged_clock() {
    let model = model(50, false);
    for selector in SELECTORS {
        for wait in ["inverse", "ziggurat"] {
            let method = json!({"kind": "static_direct", "selector": selector,
                                "waiting_time": wait});
            let result = run(&model, &request(&model, "crn-craft.run/1.0", &method, 4));
            assert_eq!(
                result.termination,
                Termination::Completed,
                "{selector}/{wait}: {:?}",
                result.replicates[0].diagnostic
            );
            for (k, record) in result.replicates.iter().enumerate() {
                // The burst and the decay of all 50 molecules.
                assert_eq!(record.events, 51, "{selector}/{wait} replicate {k}");
                assert_eq!(record.rule_events, vec![1, 50]);
                assert_eq!(final_value(&result, k, "A"), 0.0);
            }
        }
    }
}

#[test]
fn rules_run_versions_keep_the_strict_clock() {
    let model = model(50, false);
    // The same model as a rules-model, so rules-run versions accept it.
    let mut rules_model: Value = serde_json::from_slice(&model).unwrap();
    rules_model["schema"] = json!("crn-craft.rules-model/0.4");
    let rules_model = serde_json::to_vec(&rules_model).unwrap();
    for version in ["0.2", "0.4", "0.6"] {
        for wait in ["inverse", "ziggurat"] {
            if version == "0.2" && wait == "ziggurat" {
                continue;
            }
            let method = json!({"kind": "static_direct", "waiting_time": wait});
            let method = if version == "0.2" {
                json!({"kind": "static_direct"})
            } else {
                method
            };
            let schema = format!("crn-craft.rules-run/{version}");
            let result = run(&rules_model, &request(&rules_model, &schema, &method, 2));
            for record in &result.replicates {
                assert_eq!(record.termination, Termination::Failed, "{version}/{wait}");
                let diagnostic = record.diagnostic.as_ref().unwrap();
                assert_eq!(diagnostic.code, "CRN-R030");
                assert_eq!(diagnostic.path, "/time");
                assert_eq!(diagnostic.message, "clock failed to advance");
                // The burst fired; the first decay could not advance.
                assert_eq!(record.events, 1, "{version}/{wait}");
            }
        }
    }
}

#[test]
fn run_1_0_still_fails_a_runaway_rate() {
    let model = model(1, true);
    for selector in SELECTORS {
        let method = json!({"kind": "static_direct", "selector": selector});
        let result = run(&model, &request(&model, "crn-craft.run/1.0", &method, 1));
        let record = &result.replicates[0];
        assert_eq!(record.termination, Termination::Failed, "{selector}");
        let diagnostic = record.diagnostic.as_ref().unwrap();
        assert_eq!(diagnostic.code, "CRN-R030", "{selector}");
        assert_eq!(diagnostic.path, "/time");
        assert_eq!(
            diagnostic.message,
            format!("clock failed to advance in {LIMIT} consecutive steps")
        );
        // Every step after the burst stalls; the limit's steps were accepted
        // (RSSA-CR's exact bounds accept every trial).
        assert_eq!(record.events, LIMIT + 1, "{selector}");
    }
}
