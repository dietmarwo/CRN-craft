//! Stored requests of every earlier version keep their meaning: each request
//! under `stored/` replays to the result bytes recorded in
//! `stored/expected.json`, recorded with the 0.3.0 engine before
//! rules-run/0.6 changed the Static Direct defaults. The requests name their
//! models by file name in `models/rules`.
use crn_craft_rules::{Control, execute, prepare};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[test]
fn stored_requests_replay_byte_identically() {
    let stored = root().join("crates/crn-craft-rules/tests/stored");
    let expected: BTreeMap<String, String> =
        serde_json::from_slice(&std::fs::read(stored.join("expected.json")).unwrap()).unwrap();
    let mut names: Vec<String> = std::fs::read_dir(&stored)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.ends_with(".run.json"))
        .collect();
    names.sort();
    // Every stored request has a recorded result, and every record a request.
    assert_eq!(names, expected.keys().cloned().collect::<Vec<_>>());
    let mut differ = Vec::new();
    for name in &names {
        let request = std::fs::read(stored.join(name)).unwrap();
        let parsed: Value = serde_json::from_slice(&request).unwrap();
        let model_name = parsed["model"]["path"].as_str().unwrap();
        let model = std::fs::read(root().join("models/rules").join(model_name)).unwrap();
        let prepared = prepare(&model, &request).unwrap();
        let result = execute(&prepared, &Control::default()).unwrap();
        let bytes = serde_json::to_vec(&result).unwrap();
        let actual = crn_craft::sha256(&bytes);
        if actual != expected[name] {
            differ.push(format!("{name}: {actual}"));
        }
    }
    assert!(differ.is_empty(), "results changed:\n{}", differ.join("\n"));
}
