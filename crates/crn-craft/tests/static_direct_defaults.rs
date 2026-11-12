//! rules-run/0.6: Static Direct defaults to Ziggurat waits and to sparse
//! integration of every eligible tracker. Explicit options restore the
//! rules-run/0.5 behavior, and every other method runs as under 0.5.
use crn_craft::rules::{
    Control, Prepared, RulesResult, Termination, execute, generate_rules_kernel,
    generate_rules_loop, prepare,
};
use serde_json::{Value, json};
use std::collections::BTreeMap;

const THREE_CLASS: &[u8] = include_bytes!("../../../models/rules/mtdna-three-class.rules.json");
const CONTINUOUS: &[u8] =
    include_bytes!("../../../models/rules/mtdna-three-class-atp-continuous.rules.json");

const ZIGGURAT_SPARSE: &str = "chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1";
const ZIGGURAT: &str = "chacha8-u53-ziggurat/rules-static-direct/1";
const INVERSE: &str = "chacha8-u53/rules-static-direct/1";

fn counts() -> Value {
    json!([{"observable": "neutral_mutants", "keep_episodes": 2},
           {"observable": "advantaged_mutants", "keep_episodes": 2}])
}

fn passage() -> Value {
    json!([{"observable": "mutants", "first_passage": "N > 0 && mutants >= 0.002 * N"}])
}

fn request(model: &[u8], version: &str, method: Option<Value>, trackers: Option<Value>) -> Value {
    let mut request = json!({
        "schema": format!("crn-craft.rules-run/{version}"),
        "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 2000.0, "report_times": [0.0, 1000.0, 2000.0]},
        "ensemble": {"replicates": 6, "root_seed": 11},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 3, "keep_instances": 2},
        "budget": {"max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 10_000_000u64,
                   "max_instances": 10u64, "threads": 2, "max_output_bytes": 1u64 << 24},
    });
    if let Some(method) = method {
        request["method"] = method;
    }
    if let Some(trackers) = trackers {
        request["trackers"] = trackers;
    }
    request
}

fn prepared(model: &[u8], request: &Value) -> Prepared {
    prepare(model, &serde_json::to_vec(request).unwrap()).unwrap()
}

fn run(model: &[u8], request: &Value) -> RulesResult {
    let result = execute(&prepared(model, request), &Control::default()).unwrap();
    assert_eq!(result.termination, Termination::Completed);
    result
}

/// The numbers of a result: without its schema, request hash and the echoed
/// method, which name the request rather than what ran.
fn numbers(result: &RulesResult) -> Value {
    let mut value = serde_json::to_value(result).unwrap();
    let object = value.as_object_mut().unwrap();
    for key in ["schema", "run_sha256", "method"] {
        object.remove(key);
    }
    value
}

fn r051(result: &RulesResult) -> bool {
    result.diagnostics.iter().any(|d| d.code == "CRN-R051")
}

#[test]
fn static_direct_defaults_to_ziggurat_waits_and_sparse_trackers() {
    let default = json!({"kind": "static_direct"});
    let new = run(
        THREE_CLASS,
        &request(THREE_CLASS, "0.6", Some(default.clone()), Some(counts())),
    );
    assert_eq!(new.sampler, ZIGGURAT_SPARSE);
    assert!(new.diagnostics.is_empty(), "{:?}", new.diagnostics);
    // The result is rules-result/0.3 and echoes the request's method.
    let value = serde_json::to_value(&new).unwrap();
    assert_eq!(value["schema"], "crn-craft.rules-result/0.3");
    assert_eq!(value["method"], default);
    // The same numbers as the explicit opt-ins of rules-run/0.4 and 0.5.
    let both =
        json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true});
    for version in ["0.4", "0.5"] {
        let old = run(
            THREE_CLASS,
            &request(THREE_CLASS, version, Some(both.clone()), Some(counts())),
        );
        assert_eq!(old.sampler, ZIGGURAT_SPARSE);
        assert!(
            numbers(&new) == numbers(&old),
            "0.6 default differs from {version}"
        );
    }
}

#[test]
fn without_an_eligible_tracker_the_default_is_silent() {
    let both = json!({"kind": "static_direct", "waiting_time": "ziggurat"});
    for trackers in [None, Some(passage())] {
        let new = run(
            THREE_CLASS,
            &request(
                THREE_CLASS,
                "0.6",
                Some(json!({"kind": "static_direct"})),
                trackers.clone(),
            ),
        );
        // Sparse integration is the default, not a request: no CRN-R051.
        assert_eq!(new.sampler, ZIGGURAT);
        assert!(new.diagnostics.is_empty(), "{:?}", new.diagnostics);
        let old = run(
            THREE_CLASS,
            &request(THREE_CLASS, "0.4", Some(both.clone()), trackers),
        );
        assert!(numbers(&new) == numbers(&old));
    }
}

#[test]
fn explicit_options_restore_the_rules_run_05_behavior() {
    let off = json!({"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": false});
    for trackers in [Some(counts()), Some(passage()), None] {
        let new = run(
            THREE_CLASS,
            &request(THREE_CLASS, "0.6", Some(off.clone()), trackers.clone()),
        );
        let old = run(
            THREE_CLASS,
            &request(
                THREE_CLASS,
                "0.5",
                Some(json!({"kind": "static_direct"})),
                trackers,
            ),
        );
        assert_eq!(new.sampler, INVERSE);
        assert!(numbers(&new) == numbers(&old));
    }
    // Each option alone keeps the 0.6 default of the other one.
    let cases = [
        (
            json!({"kind": "static_direct", "waiting_time": "inverse"}),
            json!({"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": true}),
            "chacha8-u53/rules-static-direct-sparse-trackers/1",
        ),
        (
            json!({"kind": "static_direct", "sparse_trackers": false}),
            json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": false}),
            ZIGGURAT,
        ),
        (
            json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true}),
            json!({"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": true}),
            ZIGGURAT_SPARSE,
        ),
    ];
    for (method, explicit, sampler) in cases {
        let new = run(
            THREE_CLASS,
            &request(THREE_CLASS, "0.6", Some(method.clone()), Some(counts())),
        );
        let old = run(
            THREE_CLASS,
            &request(THREE_CLASS, "0.5", Some(explicit), Some(counts())),
        );
        assert_eq!(new.sampler, sampler, "{method}");
        assert!(numbers(&new) == numbers(&old), "{method}");
    }
}

#[test]
fn an_explicit_true_without_an_eligible_tracker_still_warns() {
    let method = json!({"kind": "static_direct", "sparse_trackers": true});
    let new = run(
        THREE_CLASS,
        &request(THREE_CLASS, "0.6", Some(method), Some(passage())),
    );
    assert_eq!(new.sampler, ZIGGURAT);
    let warning = new
        .diagnostics
        .iter()
        .find(|d| d.code == "CRN-R051")
        .expect("ineligible explicit request not reported");
    assert_eq!(warning.path, "/method/sparse_trackers");
}

#[test]
fn other_methods_run_as_under_rules_run_05() {
    let methods = [
        None,
        Some(json!({"kind": "tree"})),
        Some(json!({"kind": "next_reaction"})),
        Some(json!({"kind": "tau", "epsilon": 0.03})),
        Some(json!({"kind": "binomial", "step": 1.0})),
        Some(json!({"kind": "mtdna_direct"})),
        Some(json!({"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9})),
    ];
    for method in methods {
        for trackers in [Some(counts()), None] {
            let new = run(
                THREE_CLASS,
                &request(THREE_CLASS, "0.6", method.clone(), trackers.clone()),
            );
            let old = run(
                THREE_CLASS,
                &request(THREE_CLASS, "0.5", method.clone(), trackers),
            );
            assert!(numbers(&new) == numbers(&old), "{method:?}");
            assert!(
                serde_json::to_value(&new).unwrap()["method"]
                    == serde_json::to_value(&old).unwrap()["method"]
            );
        }
    }
    // hybrid with continuous variables is a rules-run/0.5 feature and stays.
    let hybrid = json!({"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9});
    let mut short = request(CONTINUOUS, "0.6", Some(hybrid.clone()), None);
    short["time"] = json!({"end": 200.0, "report_times": [0.0, 100.0, 200.0]});
    let new = run(CONTINUOUS, &short);
    short["schema"] = json!("crn-craft.rules-run/0.5");
    let old = run(CONTINUOUS, &short);
    assert_eq!(new.continuous, ["ATP"]);
    assert!(numbers(&new) == numbers(&old));
}

#[test]
fn mtdna_direct_falls_back_to_the_rules_run_05_static_direct() {
    // A stop outside the certified form makes mtdna_direct run Static Direct.
    let mut model: Value = serde_json::from_slice(THREE_CLASS).unwrap();
    model["stop"] = json!("mutants >= 2");
    let bytes = serde_json::to_vec(&model).unwrap();
    let fallback = run(
        &bytes,
        &request(
            &bytes,
            "0.6",
            Some(json!({"kind": "mtdna_direct"})),
            Some(counts()),
        ),
    );
    assert_eq!(fallback.sampler, INVERSE);
    assert!(r051(&fallback));
    let old = run(
        &bytes,
        &request(
            &bytes,
            "0.6",
            Some(
                json!({"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": false}),
            ),
            Some(counts()),
        ),
    );
    let mut a = numbers(&fallback);
    let b = numbers(&old);
    a["diagnostics"] = b["diagnostics"].clone();
    assert!(a == b);
}

#[test]
fn bound_candidates_keep_the_defaults_of_their_version() {
    let template = prepared(
        THREE_CLASS,
        &request(
            THREE_CLASS,
            "0.6",
            Some(json!({"kind": "static_direct"})),
            Some(counts()),
        ),
    );
    let candidate = template
        .bind_parameters(&BTreeMap::from([("mu".to_string(), 1e-4)]))
        .unwrap();
    let result = execute(&candidate, &Control::default()).unwrap();
    assert_eq!(result.sampler, ZIGGURAT_SPARSE);
}

#[test]
fn generated_code_serves_the_new_defaults() {
    // One generated loop and kernel serve every Static Direct request of the
    // model bytes, so the 0.6 default uses the 0.5 source.
    let old = prepared(
        THREE_CLASS,
        &request(
            THREE_CLASS,
            "0.5",
            Some(json!({"kind": "static_direct"})),
            Some(counts()),
        ),
    );
    let new = prepared(
        THREE_CLASS,
        &request(
            THREE_CLASS,
            "0.6",
            Some(json!({"kind": "static_direct"})),
            Some(counts()),
        ),
    );
    assert_eq!(
        generate_rules_loop(&new).unwrap(),
        generate_rules_loop(&old).unwrap()
    );
    assert_eq!(
        generate_rules_kernel(&new).unwrap(),
        generate_rules_kernel(&old).unwrap()
    );
}

#[test]
fn discovery_lists_the_new_version() {
    let entry = crn_craft::rules::describe();
    assert_eq!(
        entry["static_direct_defaults_request"],
        "crn-craft.rules-run/0.6"
    );
    let requests = entry["accepted_requests"].as_array().unwrap();
    assert_eq!(requests.len(), 6);
    assert_eq!(requests[5], "crn-craft.rules-run/0.6");
    let schema = crn_craft::rules::schema("rules-run").unwrap();
    assert!(schema.to_string().contains("crn-craft.rules-run/0.6"));
}
