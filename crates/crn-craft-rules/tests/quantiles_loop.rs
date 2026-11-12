//! Quantile statistics (`crn-craft.run/1.0`) with a generated event loop:
//! the loop's result is byte-identical to the interpreter's, quantiles
//! included, for any thread count.
use crn_craft_rules::{Control, execute, execute_with_rules_loop, prepare};
use serde_json::json;

mod sir {
    include!("generated/loop_sir_loop.rs");
}

const SIR: &[u8] = include_bytes!("generated/loop_sir.rules.json");

#[test]
fn a_generated_loop_reports_the_interpreters_quantiles() {
    for threads in [1, 2, 4] {
        for waiting_time in ["inverse", "ziggurat"] {
            let request = serde_json::to_vec(&json!({
                "schema": "crn-craft.run/1.0",
                "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(SIR)},
                "time": {"end": 200.0, "report_times": [0.0, 20.0, 50.0, 100.0, 200.0]},
                "ensemble": {"replicates": 40, "root_seed": 5},
                "outputs": {"statistics": ["mean"], "keep_paths": 3, "keep_instances": 0,
                            "quantiles": [0.025, 0.5, 0.975]},
                "budget": {"threads": threads, "max_events_per_replicate": 1_000_000,
                           "max_sampling_steps_per_replicate": 10_000_000,
                           "max_instances": 16, "max_output_bytes": 1 << 24},
                "method": {"kind": "static_direct", "waiting_time": waiting_time}
            }))
            .unwrap();
            let prepared = prepare(SIR, &request).unwrap();
            let interpreted =
                serde_json::to_string(&execute(&prepared, &Control::default()).unwrap()).unwrap();
            let generated = serde_json::to_string(
                &execute_with_rules_loop(&prepared, &Control::default(), &sir::GeneratedRulesLoop)
                    .unwrap(),
            )
            .unwrap();
            assert!(
                interpreted == generated,
                "{threads} threads, {waiting_time}"
            );
            let value: serde_json::Value = serde_json::from_str(&generated).unwrap();
            assert_eq!(value["schema"], "crn-craft.result/1.0");
            let quantiles = value["statistics"]["quantiles"].as_array().unwrap();
            assert_eq!(quantiles.len(), 3);
            assert_eq!(quantiles[1]["p"], 0.5);
        }
    }
}
