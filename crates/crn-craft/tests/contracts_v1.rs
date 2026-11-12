//! The 1.0 contract names: `crn-craft.model/1.0`, `crn-craft.run/1.0`,
//! `crn-craft.result/1.0` and `crn-craft.check/1.0`.
//!
//! A `model/1.0` document means what it means as `rules-model/0.4`, and a
//! `run/1.0` request what it means as `rules-run/0.6`; a `run/1.0` request
//! returns `result/1.0`. Network commands refuse 1.0 documents with
//! `CRN-E032`, and `rules-run` requests refuse `model/1.0` models.
use crn_craft::diagnostic::Diagnostic;
use crn_craft::rules::{Control, RulesResult, check_model, execute, prepare};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn read(path: &str) -> Vec<u8> {
    std::fs::read(root().join(path)).unwrap()
}

/// `bytes` with its one `"schema": "..."` value replaced; every other byte
/// is kept.
fn with_schema(bytes: &[u8], schema: &str) -> Vec<u8> {
    let text = String::from_utf8(bytes.to_vec()).unwrap();
    let old = serde_json::from_slice::<Value>(bytes).unwrap()["schema"]
        .as_str()
        .unwrap()
        .to_string();
    assert_eq!(text.matches(&format!("\"{old}\"")).count(), 1);
    text.replacen(&format!("\"{old}\""), &format!("\"{schema}\""), 1)
        .into_bytes()
}

fn rehash(request: &mut Value, model: &[u8]) {
    request["model"]["sha256"] = crn_craft::sha256(model).into();
}

fn run(model: &[u8], request: &Value) -> RulesResult {
    let prepared =
        prepare(model, &serde_json::to_vec(request).unwrap()).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap()
}

/// The result without the fields that name the documents rather than what
/// ran: the schema string and the hashes that cover the rewritten bytes.
fn numbers(result: &RulesResult, hashes: &[&str]) -> Vec<u8> {
    let mut value = serde_json::to_value(result).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("schema").unwrap();
    for key in hashes {
        object.remove(*key).unwrap();
    }
    serde_json::to_vec(&value).unwrap()
}

fn immigration_request(model: &[u8], schema: &str) -> Value {
    json!({
        "schema": schema,
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 5.0, "report_times": [0.0, 5.0]},
        "ensemble": {"replicates": 16, "root_seed": 7},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 2, "keep_instances": 1},
        "budget": {"max_events_per_replicate": 100_000, "max_sampling_steps_per_replicate": 1_000_000,
                   "max_instances": 4, "threads": 2, "max_output_bytes": 1 << 22},
        "method": {"kind": "static_direct"}
    })
}

#[test]
fn every_stored_request_as_run_1_0_gives_the_rules_run_0_6_result() {
    let stored = root().join("crates/crn-craft-rules/tests/stored");
    let mut names: Vec<String> = std::fs::read_dir(&stored)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.ends_with(".run.json"))
        .collect();
    names.sort();
    // Every stored request of `stored_requests.rs` (I1).
    let expected: BTreeMap<String, String> =
        serde_json::from_slice(&std::fs::read(stored.join("expected.json")).unwrap()).unwrap();
    assert_eq!(names, expected.keys().cloned().collect::<Vec<_>>());
    for name in &names {
        let request: Value =
            serde_json::from_slice(&std::fs::read(stored.join(name)).unwrap()).unwrap();
        let model = read(&format!(
            "models/rules/{}",
            request["model"]["path"].as_str().unwrap()
        ));
        let mut old = request.clone();
        old["schema"] = "crn-craft.rules-run/0.6".into();
        let mut new = request.clone();
        new["schema"] = "crn-craft.run/1.0".into();
        let reference = run(&model, &old);
        let result = run(&model, &new);
        assert_eq!(
            serde_json::to_value(reference.schema).unwrap(),
            "crn-craft.rules-result/0.3"
        );
        assert_eq!(
            serde_json::to_value(result.schema).unwrap(),
            "crn-craft.result/1.0",
            "{name}"
        );
        assert_ne!(result.run_sha256, reference.run_sha256);
        assert_eq!(result.model_sha256, reference.model_sha256);
        assert!(
            numbers(&result, &["run_sha256"]) == numbers(&reference, &["run_sha256"]),
            "{name}: run/1.0 differs from rules-run/0.6"
        );
        // The same model bytes as model/1.0.
        let model_1_0 = with_schema(&model, "crn-craft.model/1.0");
        rehash(&mut new, &model_1_0);
        let result = run(&model_1_0, &new);
        assert!(
            numbers(&result, &["run_sha256", "model_sha256"])
                == numbers(&reference, &["run_sha256", "model_sha256"]),
            "{name}: model/1.0 differs from its rules-model version"
        );
    }
}

#[test]
fn run_1_0_accepts_every_rules_model_and_rules_run_refuses_model_1_0() {
    let model = read("models/rules/immigration-death.rules.json");
    for version in ["0.1", "0.2", "0.3", "0.4"] {
        let model = with_schema(&model, &format!("crn-craft.rules-model/{version}"));
        let result = run(&model, &immigration_request(&model, "crn-craft.run/1.0"));
        assert_eq!(
            serde_json::to_value(result.schema).unwrap(),
            "crn-craft.result/1.0"
        );
        // A run/1.0 result has the shape of rules-result/0.3.
        assert_eq!(
            result.method,
            Some(crn_craft::rules::Method::StaticDirect {
                waiting_time: None,
                sparse_trackers: None
            })
        );
    }
    let model = with_schema(&model, "crn-craft.model/1.0");
    for version in ["0.2", "0.4", "0.5", "0.6"] {
        let request = immigration_request(&model, &format!("crn-craft.rules-run/{version}"));
        let diagnostics = prepare(&model, &serde_json::to_vec(&request).unwrap())
            .err()
            .unwrap();
        assert!(
            diagnostics
                .iter()
                .any(|d| d.code == "CRN-R020" && d.path == "/schema"),
            "{version}: {diagnostics:?}"
        );
    }
    // Parameter candidates keep the 1.0 versions.
    let request = immigration_request(&model, "crn-craft.run/1.0");
    let template = prepare(&model, &serde_json::to_vec(&request).unwrap()).unwrap();
    let candidate = template
        .bind_parameters(&BTreeMap::from([("k_in".to_string(), 3.0)]))
        .unwrap();
    let result = execute(&candidate, &Control::default()).unwrap();
    assert_eq!(
        serde_json::to_value(result.schema).unwrap(),
        "crn-craft.result/1.0"
    );
}

#[test]
fn checks_of_model_1_0_report_check_1_0() {
    let model = read("models/rules/immigration-death.rules.json");
    let report = check_model(&model);
    assert_eq!(
        serde_json::to_value(report.schema).unwrap(),
        "crn-craft.rules-check/0.1"
    );
    let model = with_schema(&model, "crn-craft.model/1.0");
    let report = check_model(&model);
    assert!(report.valid, "{:?}", report.diagnostics);
    assert_eq!(
        serde_json::to_value(report.schema).unwrap(),
        "crn-craft.check/1.0"
    );
    // A model/1.0 document that does not parse still gets a check/1.0 report.
    let broken = String::from_utf8(model)
        .unwrap()
        .replacen("\"rules\"", "\"rulez\"", 1);
    let report = check_model(broken.as_bytes());
    assert!(!report.valid);
    assert_eq!(
        serde_json::to_value(report.schema).unwrap(),
        "crn-craft.check/1.0"
    );
}

#[test]
fn the_converter_writes_model_1_0() {
    let converted =
        crn_craft::rules::convert::from_network(&read("models/examples/repressilator.model.json"))
            .unwrap();
    let document: Value = serde_json::from_slice(&converted).unwrap();
    assert_eq!(document["schema"], "crn-craft.model/1.0");
    let report = check_model(&converted);
    assert!(report.valid);
    assert_eq!(
        serde_json::to_value(report.schema).unwrap(),
        "crn-craft.check/1.0"
    );
}

fn has_e032(diagnostics: &[Diagnostic]) -> bool {
    diagnostics
        .iter()
        .any(|d| d.code == "CRN-E032" && d.path == "/schema")
}

#[test]
fn network_entry_points_refuse_1_0_and_rules_documents() {
    let network = read("models/examples/immigration-death.model.json");
    let run = read("models/examples/immigration-death.run.json");
    let rules = read("models/rules/immigration-death.rules.json");
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    for model in [converted.clone(), rules.clone()] {
        let report = crn_craft::check_bytes(&model, None);
        assert!(!report.valid);
        assert!(has_e032(&report.diagnostics), "{:?}", report.diagnostics);
        assert!(has_e032(&[crn_craft::parse_model(&model).unwrap_err()]));
        // Analysis reads every model version since M5b.
        let analysis = crn_craft::analyze_bytes(&model, None);
        assert!(analysis.success, "{:?}", analysis.diagnostics);
        // So does the GillesPy2 export since M5f.
        let export = crn_craft::export::export_gillespy2(&model);
        assert!(export.success, "{:?}", export.diagnostics);
        // A request bound to the 1.0 model still fails on the model.
        let mut request: Value = serde_json::from_slice(&run).unwrap();
        rehash(&mut request, &model);
        let request = serde_json::to_vec(&request).unwrap();
        let control = crn_craft::simulation::Control::default();
        let errors = crn_craft::simulation::simulate_bytes(&model, &request, &control)
            .err()
            .unwrap();
        assert!(has_e032(&errors), "{errors:?}");
    }
    // Requests of the rules family are refused as requests.
    let request = immigration_request(&network, "crn-craft.run/1.0");
    let request = serde_json::to_vec(&request).unwrap();
    let report = crn_craft::check_bytes(&network, Some(&request));
    assert!(has_e032(&report.diagnostics), "{:?}", report.diagnostics);
    assert!(has_e032(&[crn_craft::parse_run(&request).unwrap_err()]));
    let errors = crn_craft::ssa::prepare(&network, &request).err().unwrap();
    assert!(has_e032(&errors), "{errors:?}");
    let control = crn_craft::simulation::Control::default();
    let errors = crn_craft::tau::solve(&network, &request, &control)
        .err()
        .unwrap();
    assert!(has_e032(&errors), "{errors:?}");
    let errors = crn_craft::fsp::solve(&network, &request, &control)
        .err()
        .unwrap();
    assert!(has_e032(&errors), "{errors:?}");
    // Other shape errors keep their codes.
    let report = crn_craft::check_bytes(br#"{"schema": "crn-craft.model/0.2"}"#, None);
    assert_eq!(report.diagnostics[0].code, "CRN-E003");
    let report = crn_craft::check_bytes(br#"{"schema": "crn-craft.model/1.0""#, None);
    assert_eq!(report.diagnostics[0].code, "CRN-E001");
}

#[test]
fn schemas_match_their_snapshots_and_list_only_their_versions() {
    for (kind, versions) in [
        (
            "rules-model",
            vec![
                "crn-craft.rules-model/0.1",
                "crn-craft.rules-model/0.2",
                "crn-craft.rules-model/0.3",
                "crn-craft.rules-model/0.4",
            ],
        ),
        (
            "rules-run",
            vec![
                "crn-craft.rules-run/0.1",
                "crn-craft.rules-run/0.2",
                "crn-craft.rules-run/0.3",
                "crn-craft.rules-run/0.4",
                "crn-craft.rules-run/0.5",
                "crn-craft.rules-run/0.6",
            ],
        ),
        (
            "rules-result",
            vec![
                "crn-craft.rules-result/0.1",
                "crn-craft.rules-result/0.2",
                "crn-craft.rules-result/0.3",
            ],
        ),
        ("rules-check", vec!["crn-craft.rules-check/0.1"]),
        ("model-1.0", vec!["crn-craft.model/1.0"]),
        ("run-1.0", vec!["crn-craft.run/1.0"]),
        ("result-1.0", vec!["crn-craft.result/1.0"]),
        ("check-1.0", vec!["crn-craft.check/1.0"]),
    ] {
        let schema = crn_craft::rules::schema(kind).unwrap();
        let snapshot: Value =
            serde_json::from_slice(&read(&format!("schemas/{kind}.schema.json"))).unwrap();
        assert!(
            schema == snapshot,
            "{kind}: regenerate schemas/{kind}.schema.json"
        );
        // run-1.0 lists the simulation and the FSP request under `oneOf`
        // (M5a); both reference the same version definition.
        let documents: Vec<&Value> = match schema.get("oneOf") {
            Some(Value::Array(shapes)) => shapes
                .iter()
                .map(|shape| {
                    let name = shape["$ref"].as_str().unwrap().rsplit('/').next().unwrap();
                    &schema["definitions"][name]
                })
                .collect(),
            _ => vec![&schema],
        };
        assert_eq!(documents.len(), if kind == "run-1.0" { 2 } else { 1 });
        let reference = documents[0]["properties"]["schema"]["$ref"]
            .as_str()
            .unwrap();
        for document in &documents {
            assert_eq!(
                document["properties"]["schema"]["$ref"], reference,
                "{kind}"
            );
        }
        let definition = reference.rsplit('/').next().unwrap();
        assert_eq!(
            schema["definitions"][definition]["enum"],
            json!(versions),
            "{kind}"
        );
    }
    // The FSP result of run/1.0 comes from the crate-level schemas.
    let fsp = crn_craft::schema("fsp-result-1.0").unwrap();
    let snapshot: Value =
        serde_json::from_slice(&read("schemas/fsp-result-1.0.schema.json")).unwrap();
    assert!(
        fsp == snapshot,
        "regenerate schemas/fsp-result-1.0.schema.json"
    );
    assert_eq!(
        fsp["definitions"]["FspVersion"]["enum"],
        json!(["crn-craft.fsp-result/1.0"])
    );
    // The 1.0 model schema carries the corrected `continuous` description.
    let model = crn_craft::rules::schema("model-1.0").unwrap();
    let text = model["properties"]["continuous"]["description"]
        .as_str()
        .unwrap();
    assert!(!text.contains("crate::"), "{text}");
    assert!(text.contains("crn_craft::rules::continuous"), "{text}");
}

#[test]
fn discovery_names_the_1_0_contracts_and_keeps_the_rules_lists() {
    let entry = crn_craft::rules::describe();
    assert_eq!(entry["contracts_1_0"]["model"], "crn-craft.model/1.0");
    assert_eq!(entry["contracts_1_0"]["request"], "crn-craft.run/1.0");
    assert_eq!(entry["contracts_1_0"]["result"], "crn-craft.result/1.0");
    assert_eq!(entry["contracts_1_0"]["check"], "crn-craft.check/1.0");
    assert_eq!(
        entry["contracts_1_0"]["schemas"],
        json!(["model-1.0", "run-1.0", "result-1.0", "check-1.0"])
    );
    let accepted = entry["contracts_1_0"]["request_accepts_models"]
        .as_array()
        .unwrap();
    assert_eq!(accepted.last().unwrap(), "crn-craft.model/1.0");
    assert_eq!(accepted.len(), 5);
    assert_eq!(entry["convert"]["to"], "crn-craft.model/1.0");
    // The rules-family lists are unchanged.
    assert_eq!(entry["accepted_requests"].as_array().unwrap().len(), 6);
    assert_eq!(entry["accepted_models"].as_array().unwrap().len(), 4);
}
