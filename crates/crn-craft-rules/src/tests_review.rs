//! Regression tests from the review of the mtDNA A0–A6 patch.
use super::*;
use serde_json::{Value, json};

const THREE_CLASS: &str = include_str!("../../../models/rules/mtdna-three-class.rules.json");

fn three_class() -> Value {
    serde_json::from_str(THREE_CLASS).unwrap()
}

fn request_v(schema: &str, model: &[u8], replicates: u64, end: f64, times: &[f64]) -> Value {
    json!({
        "schema": schema,
        "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": times},
        "ensemble": {"replicates": replicates, "root_seed": 23},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 10_000_000u64,
                   "max_sampling_steps_per_replicate": 100_000_000u64,
                   "max_instances": 10_000u64, "threads": 2, "max_output_bytes": 1u64 << 30}
    })
}

fn prepare_values(model: &Value, request: &Value) -> Result<Prepared, Vec<Diagnostic>> {
    prepare(
        &serde_json::to_vec(model).unwrap(),
        &serde_json::to_vec(request).unwrap(),
    )
}

/// Prepare `model` with a request edited by `edit`; the request hash is kept.
fn prepare_with(
    model: &Value,
    schema: &str,
    edit: impl FnOnce(&mut Value),
) -> Result<Prepared, Vec<Diagnostic>> {
    let bytes = serde_json::to_vec(model).unwrap();
    let mut req = request_v(schema, &bytes, 8, 50.0, &[0.0, 25.0, 50.0]);
    edit(&mut req);
    prepare_values(model, &req)
}

fn first_code(result: Result<Prepared, Vec<Diagnostic>>) -> Option<(String, String)> {
    result.err().map(|d| (d[0].code.clone(), d[0].path.clone()))
}

const V01: &str = "crn-craft.rules-run/0.1";
const V02: &str = "crn-craft.rules-run/0.2";

#[test]
fn method_objects_reject_unknown_fields() {
    let model = three_class();
    for method in [
        json!({"kind": "static_direct", "epsilon": 0.5}),
        json!({"kind": "tree", "bogus": 1}),
        json!({"kind": "next_reaction", "step": 1.0}),
        json!({"kind": "thinning", "x": 0}),
    ] {
        let result = prepare_with(&model, V02, |r| r["method"] = method.clone());
        assert!(result.is_err(), "accepted {method}");
    }
    assert!(
        prepare_with(&model, V02, |r| r["method"] =
            json!({"kind": "static_direct"}))
        .is_ok()
    );
}

#[test]
fn version_01_contracts_reject_new_keys_even_when_empty() {
    let model = three_class();
    for (key, value) in [
        ("trackers", json!([])),
        ("method", Value::Null),
        ("event_log", Value::Null),
    ] {
        let result = prepare_with(&model, V01, |r| r[key] = value.clone());
        assert!(result.is_err(), "rules-run/0.1 accepted {key}");
    }
    // A present value must not be null in later versions either.
    assert!(prepare_with(&model, V02, |r| r["method"] = Value::Null).is_err());
    let mut with_atp = three_class();
    with_atp["atp"] = Value::Null;
    assert!(!check_model(&serde_json::to_vec(&with_atp).unwrap()).valid);
}

#[test]
fn method_and_log_diagnostics_name_the_right_field() {
    let model = three_class();
    let thinning = prepare_with(&model, V02, |r| r["method"] = json!({"kind": "thinning"}));
    let diagnostics = thinning.err().unwrap();
    assert_eq!(diagnostics[0].path, "/method");
    assert!(
        diagnostics[0].message.contains("thinning requires an ATP"),
        "{}",
        diagnostics[0].message
    );
    let log = prepare_with(&model, V02, |r| {
        r["event_log"] = json!({"max_entries": 4, "rules": ["no_such_rule"]});
    });
    assert_eq!(
        first_code(log),
        Some(("CRN-R020".into(), "/event_log".into()))
    );
}
