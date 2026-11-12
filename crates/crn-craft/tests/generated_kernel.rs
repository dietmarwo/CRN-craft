//! The generator's output for the three-class fixture is checked in, so the
//! compiled-kernel path runs under `cargo test` without a separate build.
//! Regenerate it with
//! `cargo run --release -p crn-craft-rules --example generate_rules_kernel -- MODEL REQUEST`
//! when the generator changes; the first test reports the difference.
use crn_craft::rules::{
    Control, execute, execute_with_rules_kernel, generate_rules_kernel, prepare,
};
use serde_json::{Value, json};

include!("generated/mtdna_three_class_kernel.rs");

const MODEL: &[u8] = include_bytes!("../../../models/rules/mtdna-three-class.rules.json");
const GOLDEN: &str = include_str!("generated/mtdna_three_class_kernel.rs");

fn request(model: &[u8], method: Value, schema: &str, threads: u32) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": schema,
        "model": {"path": "mtdna-three-class.rules.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 2000.0, "report_times": [0.0, 500.0, 1000.0, 2000.0]},
        "ensemble": {"replicates": 12, "root_seed": 41},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 12, "keep_instances": 0},
        "budget": {"threads": threads, "max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 10_000_000u64, "max_instances": 1000u64,
                   "max_output_bytes": 1u64 << 28},
        "method": method,
        "trackers": [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]
    }))
    .unwrap()
}

fn with_mu(mu: f64) -> Vec<u8> {
    let mut model: Value = serde_json::from_slice(MODEL).unwrap();
    model["parameters"]["mu"] = json!(mu);
    serde_json::to_vec(&model).unwrap()
}

#[test]
fn generator_output_matches_the_checked_in_kernel() {
    let req = request(
        MODEL,
        json!({"kind": "static_direct"}),
        "crn-craft.rules-run/0.2",
        1,
    );
    let prepared = prepare(MODEL, &req).unwrap();
    let source = generate_rules_kernel(&prepared).unwrap();
    assert!(
        source == GOLDEN,
        "generator output changed; regenerate the fixture"
    );
}

#[test]
fn generated_kernel_replays_the_interpreter_byte_for_byte() {
    for (method, schema) in [
        (json!({"kind": "static_direct"}), "crn-craft.rules-run/0.2"),
        (
            json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true}),
            "crn-craft.rules-run/0.4",
        ),
        // The rules-run/0.6 default: Ziggurat waits and sparse trackers.
        (json!({"kind": "static_direct"}), "crn-craft.rules-run/0.6"),
    ] {
        for threads in [1, 3] {
            let prepared =
                prepare(MODEL, &request(MODEL, method.clone(), schema, threads)).unwrap();
            let interpreted = execute(&prepared, &Control::default()).unwrap();
            let compiled =
                execute_with_rules_kernel(&prepared, &Control::default(), &GeneratedRulesKernel)
                    .unwrap();
            assert!(
                serde_json::to_vec(&interpreted).unwrap() == serde_json::to_vec(&compiled).unwrap(),
                "{method} with {threads} workers"
            );
        }
    }
}

#[test]
fn a_kernel_for_other_model_bytes_is_rejected() {
    // Parameters are compiled into the kernel, so a rebinding is a new model.
    let other = with_mu(4e-5);
    let req = request(
        &other,
        json!({"kind": "static_direct"}),
        "crn-craft.rules-run/0.2",
        1,
    );
    let prepared = prepare(&other, &req).unwrap();
    let Err(error) =
        execute_with_rules_kernel(&prepared, &Control::default(), &GeneratedRulesKernel)
    else {
        panic!("stale kernel accepted");
    };
    assert_eq!(error.code, "CRN-R020");
    assert_eq!(error.path, "/model");
}
