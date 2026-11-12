use serde_json::Value;
use std::{
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn invoke(args: &[&str]) -> (i32, Value) {
    let output = Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .current_dir(root())
        .args(args)
        .output()
        .unwrap();
    assert!(output.stderr.is_empty(), "{output:?}");
    (
        output.status.code().unwrap(),
        serde_json::from_slice(&output.stdout).unwrap(),
    )
}

#[test]
fn cli_description_and_example_check() {
    let (code, value) = invoke(&["describe"]);
    assert_eq!(code, 0);
    // The unified discovery document (M7): crn-craft.describe/1.0 replaces
    // the native describe/0.15 and the core describe/0.11.
    assert_eq!(value, crn_craft::describe());
    assert_eq!(value["schema"], "crn-craft.describe/1.0");
    assert_eq!(value["milestone"], "unified-engine");
    assert_eq!(value["crate"], "crn-craft");
    assert_eq!(
        value["legacy"]["rules"]["convert"]["from"],
        "crn-craft.model/0.1"
    );
    assert_eq!(
        value["legacy"]["rules"]["static_direct_defaults_request"],
        "crn-craft.rules-run/0.6"
    );
    let requests = value["legacy"]["requests"].as_array().unwrap();
    assert_eq!(requests.first().unwrap(), "crn-craft.rules-run/0.1");
    assert_eq!(requests.last().unwrap(), "crn-craft.rules-run/0.6");
    // A network model with a network run/0.1 request: no request names a
    // network model since the network engine's check/0.1 was removed (M8b).
    let (code, value) = invoke(&[
        "check",
        "models/legacy/immigration-death.model.json",
        "models/legacy/immigration-death.run.json",
    ]);
    assert_eq!(code, 1);
    assert_eq!(value["schema"], "crn-craft.check/1.0");
    assert_eq!(value["valid"], false);
    assert_eq!(value["diagnostics"][0]["code"], "CRN-E033");
    // The converted example with a run/1.0 request.
    let temp = staged();
    let (code, value) = invoke(&[
        "check",
        temp.0
            .join("immigration-death.model.json")
            .to_str()
            .unwrap(),
        temp.0.join("immigration-death.run.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{value}");
    assert_eq!(value["schema"], "crn-craft.check/1.0");
    assert_eq!(value["valid"], true);
}

#[test]
fn cli_errors_are_json_and_exit_one() {
    for args in [
        vec![],
        vec!["run"],
        vec!["check"],
        vec!["describe", "missing"],
        vec!["describe", "model", "extra"],
    ] {
        let (code, value) = invoke(&args);
        assert_eq!(code, 1);
        assert_eq!(value["diagnostics"][0]["code"], "CRN-E031");
    }
    for path in ["not-present.json", "models"] {
        let (code, value) = invoke(&["check", path]);
        assert_eq!(code, 1);
        assert_eq!(value["diagnostics"][0]["code"], "CRN-E030");
    }
}

#[test]
fn cli_negative_fixtures_have_expected_exit_code_and_fields() {
    let base = root().join("models/invalid");
    let entries: Vec<Value> =
        serde_json::from_slice(&std::fs::read(base.join("expected.json")).unwrap()).unwrap();
    for entry in entries {
        let model = base.join(entry["model"].as_str().unwrap());
        let run = entry["run"].as_str().map(|name| base.join(name));
        let mut args = vec!["check", model.to_str().unwrap()];
        if let Some(run) = &run {
            args.push(run.to_str().unwrap());
        }
        let (code, value) = invoke(&args);
        assert_eq!(code, 1, "{entry}");
        assert!(
            value["diagnostics"]
                .as_array()
                .unwrap()
                .iter()
                .any(|d| d["code"] == entry["code"] && d["path"] == entry["path"]),
            "{entry}: {value}"
        );
    }
}

#[test]
fn analyze_supports_named_state_and_reports_real_capabilities() {
    let temp = staged();
    let state = temp.0.join("state.json");
    std::fs::write(
        &state,
        br#"{"schema":"crn-craft.state/0.1","counts":{"A":5}}"#,
    )
    .unwrap();
    let model = temp.0.join("immigration-death.model.json");
    let (code, report) = invoke(&["analyze", model.to_str().unwrap(), state.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert_eq!(report["schema"], "crn-craft.analysis/1.0");
    assert_eq!(report["success"], true);
    assert_eq!(report["simulation_available"], false);
    assert_eq!(
        report["analysis"]["propensities"],
        serde_json::json!([4.0, 1.0])
    );
}

#[test]
fn analyze_state_errors_remain_coded_json() {
    let temp = staged();
    let state = temp.0.join("state.json");
    let model = temp.0.join("immigration-death.model.json");
    for (text, expected) in [
        (
            r#"{"schema":"crn-craft.state/0.1","counts":{"A":-1}}"#,
            "CRN-E041",
        ),
        (
            r#"{"schema":"crn-craft.state/0.1","counts":{"A":1,"A":2}}"#,
            "CRN-E002",
        ),
        (
            r#"{"schema":"crn-craft.state/0.1","counts":{"A":1},"extra":0}"#,
            "CRN-E003",
        ),
        (
            r#"{"schema":"crn-craft.state/99","counts":{"A":1}}"#,
            "CRN-E003",
        ),
    ] {
        std::fs::write(&state, text).unwrap();
        let (code, report) = invoke(&["analyze", model.to_str().unwrap(), state.to_str().unwrap()]);
        assert_eq!(code, 1);
        assert_eq!(report["success"], false);
        assert!(report["analysis"].is_null());
        assert_eq!(report["diagnostics"][0]["code"], expected);
    }
    let (code, report) = invoke(&["analyze", model.to_str().unwrap(), "missing-state.json"]);
    assert_eq!(code, 1);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E030");
}

struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        static SERIAL: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "crn-craft-test-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The converted immigration-death example and a `crn-craft.run/1.0`
/// Static Direct request for it, side by side in a new directory.
fn staged() -> Temp {
    let temp = Temp::new();
    let model = std::fs::read(root().join("models/examples/immigration-death.model.json")).unwrap();
    std::fs::write(temp.0.join("immigration-death.model.json"), &model).unwrap();
    let request = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "immigration-death.model.json", "sha256": crn_craft::sha256(&model)},
        "method": {"kind": "static_direct"},
        "time": {"end": 10.0, "report_times": [0.0, 1.0, 2.0, 5.0, 10.0]},
        "ensemble": {"replicates": 64, "root_seed": 42},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 2, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1000000, "max_sampling_steps_per_replicate": 8000000,
                   "max_instances": 100000, "threads": 1, "max_output_bytes": 1048576}
    });
    std::fs::write(
        temp.0.join("immigration-death.run.json"),
        serde_json::to_vec_pretty(&request).unwrap(),
    )
    .unwrap();
    temp
}

#[test]
fn reference_is_relative_to_run_directory_not_working_directory() {
    let temp = staged();
    let (code, value) = invoke(&[
        "check",
        temp.0
            .join("immigration-death.model.json")
            .to_str()
            .unwrap(),
        temp.0.join("immigration-death.run.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{value}");
}

#[test]
fn different_file_with_identical_bytes_is_not_the_referenced_path() {
    let temp = staged();
    let copy = temp.0.join("other.json");
    std::fs::copy(temp.0.join("immigration-death.model.json"), &copy).unwrap();
    let (code, value) = invoke(&[
        "check",
        copy.to_str().unwrap(),
        temp.0.join("immigration-death.run.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 1);
    assert!(
        value["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "CRN-R020" && d["path"] == "/model/path")
    );
}

#[cfg(unix)]
#[test]
fn symlink_outside_run_directory_rejected_even_when_hash_matches() {
    let temp = staged();
    let inside = temp.0.join("immigration-death.model.json");
    std::fs::remove_file(&inside).unwrap();
    std::os::unix::fs::symlink(
        root()
            .join("models/examples/immigration-death.model.json")
            .canonicalize()
            .unwrap(),
        &inside,
    )
    .unwrap();
    let (code, value) = invoke(&[
        "check",
        inside.to_str().unwrap(),
        temp.0.join("immigration-death.run.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 1);
    assert!(
        value["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "CRN-R020" && d["path"] == "/model/path")
    );
}

#[test]
fn analyze_reads_network_rules_and_1_0_models_and_analyze_v1_is_gone() {
    let temp = Temp::new();
    let network = "models/legacy/enzyme-mechanism.model.json";
    let (code, report) = invoke(&["analyze", network]);
    assert_eq!(code, 0, "{report}");
    assert_eq!(report["schema"], "crn-craft.analysis/1.0");
    let expected = crn_craft::analyze_bytes(&std::fs::read(root().join(network)).unwrap(), None);
    assert_eq!(report, serde_json::to_value(expected).unwrap());
    // The converted model gives the same analysis.
    let converted =
        crn_craft::rules::convert::from_network(&std::fs::read(root().join(network)).unwrap())
            .unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let path = temp.0.join("model.json");
    let (code, other) = invoke(&["analyze", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{other}");
    assert_eq!(other["analysis"], report["analysis"]);
    assert_eq!(other["model_sha256"], crn_craft::sha256(&converted));
    // A rules model of an older version.
    let (code, rules) = invoke(&["analyze", "models/rules/immigration-death.rules.json"]);
    assert_eq!(code, 0, "{rules}");
    assert_eq!(rules["schema"], "crn-craft.analysis/1.0");
    // The legacy command is gone; its schema is no longer described.
    let (code, report) = invoke(&["analyze-v1", network]);
    assert_eq!(code, 1);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E031");
    let (code, _) = invoke(&["describe", "analysis-v1"]);
    assert_eq!(code, 1);
    let (code, schema) = invoke(&["describe", "analysis-1.0"]);
    assert_eq!(code, 0);
    let snapshot: Value = serde_json::from_slice(
        &std::fs::read(root().join("schemas/analysis-1.0.schema.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(schema, snapshot);
    let (_, description) = invoke(&["describe"]);
    assert!(
        !description["commands"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c == "analyze-v1")
    );
    assert_eq!(
        description["contracts"]["analysis"],
        "crn-craft.analysis/1.0"
    );
}

#[test]
fn export_cli_is_json_and_rejects_unsupported_or_missing_inputs() {
    let (code, report) = invoke(&[
        "export",
        "models/legacy/immigration-death.model.json",
        "gillespy2",
    ]);
    assert_eq!(code, 0);
    assert_eq!(report["schema"], "crn-craft.export/0.1");
    assert_eq!(report["simulation_available"], false);
    assert!(
        report["artifact"]["source"]
            .as_str()
            .unwrap()
            .contains("def build_model():")
    );
    // Michaelis-Menten, refused by 0.6.0, exports a custom propensity
    // (M5f); a model with ATP has no flat view.
    let (code, report) = invoke(&[
        "export",
        "models/legacy/michaelis-menten.model.json",
        "gillespy2",
    ]);
    assert_eq!(code, 0, "{report}");
    assert!(
        report["artifact"]["source"]
            .as_str()
            .unwrap()
            .contains("propensity_function=")
    );
    for (args, expected) in [
        (
            vec![
                "export",
                "models/rules/mtdna-three-class-atp.rules.json",
                "gillespy2",
            ],
            "CRN-R017",
        ),
        (
            vec!["export", "missing-model.json", "gillespy2"],
            "CRN-E030",
        ),
        (vec!["export"], "CRN-E031"),
        (vec!["analyze-v1"], "CRN-E031"),
        (
            vec![
                "export",
                "models/legacy/immigration-death.model.json",
                "unknown",
            ],
            "CRN-E031",
        ),
    ] {
        let (code, report) = invoke(&args);
        assert_eq!(code, 1);
        assert!(
            report["diagnostics"]
                .as_array()
                .unwrap()
                .iter()
                .any(|d| d["code"] == expected)
        );
    }
}

#[test]
fn tau_accepts_a_cancel_file_and_stops_before_any_replicate() {
    let temp = Temp::new();
    let model = br#"{"schema":"crn-craft.model/1.0","name":"T","parameters":{"k":20.0},"types":[{"name":"X"}],"initial":[{"type":"X","count":30}],"observables":{"X":{"type":"X"}},"rules":[{"name":"birth","reactants":[],"rate":"k","products":[{"type":"X"}]}]}"#;
    std::fs::write(temp.0.join("model.json"), model).unwrap();
    let request = serde_json::json!({"schema":"crn-craft.run/1.0","model":{"path":"model.json","sha256":crn_craft::sha256(model)},
        "method":{"kind":"tau","epsilon":0.05},"ensemble":{"replicates":4,"root_seed":1},
        "time":{"end":1.0,"report_times":[0.0,1.0]},"outputs":{"statistics":["mean"],"keep_paths":0,"keep_instances":0},
        "budget":{"max_events_per_replicate":100000,"max_sampling_steps_per_replicate":1000000,"max_instances":1000,"threads":1,"max_output_bytes":100000}});
    std::fs::write(
        temp.0.join("tau.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    std::fs::write(temp.0.join("cancel"), b"").unwrap();
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .arg("run")
        .arg(temp.0.join("model.json"))
        .arg(temp.0.join("tau.json"))
        .arg("--cancel-file")
        .arg(temp.0.join("cancel"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    let value: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(value["termination"], "cancelled");
    assert_eq!(value["completed_replicates"], 0);
    assert!(value["statistics"].is_null());
}

#[test]
fn convert_prints_the_converted_model_bytes_or_coded_refusals() {
    let path = "models/legacy/repressilator.model.json";
    let output = Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .current_dir(root())
        .args(["convert", path])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert!(output.stderr.is_empty());
    let network = std::fs::read(root().join(path)).unwrap();
    assert_eq!(
        output.stdout,
        crn_craft::rules::convert::from_network(&network).unwrap()
    );
    assert!(crn_craft::rules::check_model(&output.stdout).valid);
    // A population limit has no rules equivalent.
    let temp = Temp::new();
    let limited = temp.0.join("limited.model.json");
    let mut model: Value = serde_json::from_slice(&network).unwrap();
    model["species"][0]["population_limit"] = 1000.into();
    std::fs::write(&limited, serde_json::to_vec(&model).unwrap()).unwrap();
    let (code, value) = invoke(&["convert", limited.to_str().unwrap()]);
    assert_eq!(code, 1);
    assert_eq!(value["success"], false);
    assert_eq!(value["diagnostics"][0]["code"], "CRN-R016");
    assert_eq!(
        value["diagnostics"][0]["path"],
        "/species/0/population_limit"
    );
    let (code, value) = invoke(&["convert", "not-present.json"]);
    assert_eq!(code, 1);
    assert_eq!(value["diagnostics"][0]["code"], "CRN-E030");
    let (code, value) = invoke(&["convert"]);
    assert_eq!(code, 1);
    assert_eq!(value["diagnostics"][0]["code"], "CRN-E031");
}

#[test]
fn removed_network_commands_are_usage_errors_and_rules_commands_run_1_0_documents() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let network = std::fs::read(root().join("models/legacy/immigration-death.model.json")).unwrap();
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let mut run: Value = serde_json::from_slice(
        &std::fs::read(root().join("models/legacy/immigration-death.run.json")).unwrap(),
    )
    .unwrap();
    run["model"] =
        serde_json::json!({"path": "model.json", "sha256": crn_craft::sha256(&converted)});
    std::fs::write(temp.0.join("run.json"), serde_json::to_vec(&run).unwrap()).unwrap();
    let request = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(&converted)},
        "time": {"end": 5.0, "report_times": [0.0, 5.0]},
        "ensemble": {"replicates": 8, "root_seed": 1},
        "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 1000000,
                   "max_instances": 1, "threads": 1, "max_output_bytes": 1048576},
        "method": {"kind": "static_direct"}
    });
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let names = [
        path("model.json"),
        path("run.json"),
        path("request.json"),
        path("out"),
    ];
    let [model, run, request, out] = names.each_ref().map(String::as_str);
    // The network engine's commands were removed with it (M8b).
    for args in [
        vec!["simulate", model, run, out],
        vec!["simulate-cached", model, run, out],
        vec!["fsp", model, run],
        vec!["ssa", model, request],
        vec!["tau", model, request],
    ] {
        let (code, report) = invoke(&args);
        assert_eq!(code, 1, "{args:?}");
        assert_eq!(report["diagnostics"][0]["code"], "CRN-E031", "{args:?}");
        assert_eq!(report["schema"], "crn-craft.check/1.0", "{args:?}");
    }
    assert!(!temp.0.join("out").exists());
    // `check` reads every model version (M7): model/1.0 alone is valid; a
    // network run/0.1 request is no request of the rules engine.
    let (code, report) = invoke(&["check", model]);
    assert_eq!(
        (code, &report["schema"]),
        (0, &Value::from("crn-craft.check/1.0"))
    );
    let (code, report) = invoke(&["check", model, run]);
    assert_eq!(
        (code, &report["valid"]),
        (1, &Value::Bool(false)),
        "{report}"
    );
    let (code, report) = invoke(&["rules-check", model]);
    assert_eq!(
        (code, &report["schema"]),
        (0, &Value::from("crn-craft.check/1.0"))
    );
    let (code, result) = invoke(&["rules", model, request]);
    assert_eq!(code, 0, "{result}");
    assert_eq!(result["schema"], "crn-craft.result/1.0");
    for kind in ["model-1.0", "run-1.0", "result-1.0", "check-1.0"] {
        let (code, schema) = invoke(&["describe", kind]);
        assert_eq!(code, 0);
        let snapshot: Value = serde_json::from_slice(
            &std::fs::read(root().join(format!("schemas/{kind}.schema.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(schema, snapshot, "{kind}");
    }
}

#[test]
fn rules_solves_fsp_requests_and_describes_their_contracts() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let network = std::fs::read(root().join("models/legacy/immigration-death.model.json")).unwrap();
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let mut request = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(&converted)},
        "time": {"end": 5.0, "report_times": [0.0, 1.0, 5.0]},
        "method": {"kind": "fsp", "bounds": {"A": 100}, "max_states": 200, "tolerance": 1e-12},
        "budget": {"max_output_bytes": 10000000}
    });
    let bytes = serde_json::to_vec(&request).unwrap();
    std::fs::write(temp.0.join("request.json"), &bytes).unwrap();
    let (model, file) = (path("model.json"), path("request.json"));
    let output = Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .args(["rules", &model, &file])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert!(output.stderr.is_empty());
    // Compact JSON with sorted keys (as every CLI output) and a newline,
    // with the engine's numbers.
    let expected =
        crn_craft::fsp::v1::solve(&converted, &bytes, &crn_craft::rules::Control::default())
            .unwrap();
    let mut text = serde_json::to_vec(&serde_json::to_value(&expected).unwrap()).unwrap();
    text.push(b'\n');
    assert_eq!(output.stdout, text);
    let result: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["schema"], "crn-craft.fsp-result/1.0");
    assert_eq!(result["species"], serde_json::json!(["A"]));
    // A refused request exits 1 with its diagnostics.
    request["method"]["bounds"] = serde_json::json!({"B": 3});
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let (code, report) = invoke(&["rules", &model, &file]);
    assert_eq!(code, 1);
    assert_eq!(report["success"], false);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-R020");
    assert_eq!(report["diagnostics"][0]["path"], "/method/bounds/B");
    // An existing cancel file stops the solve.
    std::fs::write(temp.0.join("request.json"), &bytes).unwrap();
    std::fs::write(temp.0.join("cancel"), b"").unwrap();
    let cancel = path("cancel");
    let (code, report) = invoke(&["rules", &model, &file, "--cancel-file", &cancel]);
    assert_eq!(code, 1);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E070");
    for kind in ["fsp-result-1.0", "run-1.0"] {
        let (code, schema) = invoke(&["describe", kind]);
        assert_eq!(code, 0);
        let snapshot: Value = serde_json::from_slice(
            &std::fs::read(root().join(format!("schemas/{kind}.schema.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(schema, snapshot, "{kind}");
    }
    let (_, description) = invoke(&["describe"]);
    assert_eq!(
        description["contracts"]["fsp_result"],
        "crn-craft.fsp-result/1.0"
    );
}

#[test]
fn rules_archives_runs_and_verify_reads_1_0_archives() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let network = std::fs::read(root().join("models/legacy/immigration-death.model.json")).unwrap();
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let mut request = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(&converted)},
        "method": {"kind": "static_direct"},
        "time": {"end": 5.0, "report_times": [0.0, 1.0, 5.0]},
        "ensemble": {"replicates": 16, "root_seed": 9},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 16, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 10000000,
                   "max_instances": 1000, "threads": 2, "max_output_bytes": 10000000}
    });
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let (model, file) = (path("model.json"), path("request.json"));
    let plain = Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .args(["rules", &model, &file])
        .output()
        .unwrap();
    assert!(plain.status.success(), "{plain:?}");
    // Both option orders; the archived result is the printed result.
    for (name, args) in [
        ("first", vec!["--archive", "OUT"]),
        (
            "second",
            vec!["--cancel-file", "CANCEL", "--archive", "OUT"],
        ),
    ] {
        let out = path(name);
        let cancel = path("absent-cancel-file");
        let mut argv = vec!["rules", model.as_str(), file.as_str()];
        argv.extend(args.into_iter().map(|a| match a {
            "OUT" => out.as_str(),
            "CANCEL" => cancel.as_str(),
            other => other,
        }));
        let (code, report) = invoke(&argv);
        assert_eq!(code, 0, "{report}");
        assert_eq!(report["schema"], "crn-craft.archive-check/1.0");
        assert_eq!(report["valid"], true);
        assert_eq!(report["simulation_complete"], true);
        assert_eq!(report["result"], "crn-craft.result/1.0");
        assert_eq!(
            std::fs::read(temp.0.join(name).join("result.json")).unwrap(),
            plain.stdout
        );
        let (code, verified) = invoke(&["verify", &out]);
        assert_eq!(code, 0, "{verified}");
        assert_eq!(verified["schema"], "crn-craft.archive-check/1.0");
        // An existing directory is never reused.
        let (code, again) = invoke(&argv);
        assert_eq!((code, &again["valid"]), (1, &Value::Bool(false)));
    }
    // An incomplete run: a valid archive, exit 1; verify exits 0.
    request["budget"]["max_events_per_replicate"] = 1.into();
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let partial = path("partial");
    let (code, report) = invoke(&["rules", &model, &file, "--archive", &partial]);
    assert_eq!(code, 1, "{report}");
    assert_eq!(report["valid"], true);
    assert_eq!(report["simulation_complete"], false);
    assert_eq!(report["termination"], "event_limit");
    let (code, report) = invoke(&["verify", &partial]);
    assert_eq!(
        (code, &report["simulation_complete"]),
        (0, &Value::Bool(false))
    );
    std::fs::remove_file(temp.0.join("partial").join("COMPLETE")).unwrap();
    let (code, report) = invoke(&["verify", &partial]);
    assert_eq!(
        (code, &report["schema"]),
        (1, &"crn-craft.archive-check/1.0".into())
    );
    // Malformed options are usage errors.
    for args in [
        vec!["rules", &model, &file, "--archive"],
        vec!["rules", &model, &file, "--archive", "a", "--archive", "b"],
        vec!["rules", &model, &file, "--output", "a"],
    ] {
        let (code, report) = invoke(&args);
        assert_eq!(code, 1);
        assert_eq!(report["diagnostics"][0]["code"], "CRN-E031", "{args:?}");
    }
    // The file set of a network-engine archive (crn-craft.archive/0.1),
    // which verify no longer reads (M8b).
    let old = temp.0.join("old");
    std::fs::create_dir(&old).unwrap();
    for name in [
        "model.json",
        "run.json",
        "result.json",
        "provenance.json",
        "manifest.json",
        "COMPLETE",
    ] {
        std::fs::write(old.join(name), b"{}\n").unwrap();
    }
    let (code, report) = invoke(&["verify", old.to_str().unwrap()]);
    assert_eq!(code, 1);
    assert_eq!(report["schema"], "crn-craft.archive-check/1.0");
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E065");
    for kind in ["archive-1.0", "archive-check-1.0"] {
        let (code, schema) = invoke(&["describe", kind]);
        assert_eq!(code, 0);
        let snapshot: Value = serde_json::from_slice(
            &std::fs::read(root().join(format!("schemas/{kind}.schema.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(schema, snapshot, "{kind}");
    }
}

/// The visible commands, in `--help` order.
const COMMANDS: [&str; 7] = [
    "describe", "check", "run", "verify", "analyze", "convert", "export",
];
/// The network engine's commands, removed with it in M8b.
const REMOVED: [&str; 5] = ["simulate", "simulate-cached", "ssa", "tau", "fsp"];
/// The `describe` kinds of the network engine's contracts, removed with it.
const REMOVED_SCHEMAS: [&str; 13] = [
    "run",
    "result",
    "result-v1",
    "check",
    "fsp-result",
    "tau-request",
    "tau-result",
    "ssa-run",
    "ssa-result",
    "archive",
    "archive-check",
    "analysis",
    "analysis-v1",
];

fn text(args: &[&str]) -> (i32, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .current_dir(root())
        .args(args)
        .output()
        .unwrap();
    assert!(output.stderr.is_empty(), "{output:?}");
    (
        output.status.code().unwrap(),
        String::from_utf8(output.stdout).unwrap(),
    )
}

#[test]
fn help_usage_and_describe_list_the_unified_commands_only() {
    let (code, help) = text(&["--help"]);
    assert_eq!(code, 0);
    let listed: Vec<&str> = help
        .lines()
        .filter_map(|line| line.strip_prefix("  crn-craft "))
        .map(|rest| rest.split_whitespace().next().unwrap())
        .collect();
    assert_eq!(listed, [&COMMANDS[..], &["--version"]].concat());
    assert!(help.contains("Aliases: rules = run, rules-check = check."));
    for removed in REMOVED {
        assert!(
            !help.contains(&format!("crn-craft {removed} ")),
            "{removed} in --help"
        );
    }
    let (code, version) = invoke(&["--version"]);
    assert_eq!(code, 0);
    assert_eq!(
        version,
        serde_json::json!({"version": env!("CARGO_PKG_VERSION"), "milestone": "unified-engine"})
    );
    let (code, usage) = invoke(&["no-such-command"]);
    assert_eq!(code, 1);
    let hint = usage["diagnostics"][0]["hint"].as_str().unwrap();
    for command in COMMANDS {
        assert!(hint.contains(command), "{command} missing from {hint}");
    }
    for removed in ["simulate", "ssa", "tau"] {
        assert!(!hint.contains(removed), "{removed} in the usage hint");
    }
    let (code, description) = invoke(&["describe"]);
    assert_eq!(code, 0);
    assert_eq!(description["commands"], serde_json::json!(COMMANDS));
    assert_eq!(
        description["aliases"],
        serde_json::json!({"rules": "run", "rules-check": "check"})
    );
    assert_eq!(description["contracts"]["model"], "crn-craft.model/1.0");
    assert_eq!(description["contracts"]["request"], "crn-craft.run/1.0");
    assert_eq!(description["contracts"]["result"], "crn-craft.result/1.0");
    assert_eq!(description["contracts"]["check"], "crn-craft.check/1.0");
    assert_eq!(description["legacy"]["rules"]["crate"], "crn-craft");
    let legacy_models = description["legacy"]["models"].as_array().unwrap();
    assert!(legacy_models.contains(&"crn-craft.model/0.1".into()));
    assert!(legacy_models.contains(&"crn-craft.rules-model/0.4".into()));
    // The 1.0 contracts come first; every listed kind is served and equals
    // its committed schema.
    let keys: Vec<&String> = description.as_object().unwrap().keys().collect();
    assert!(keys.contains(&&"legacy".to_string()));
    let mut kinds = description["schemas"].as_array().unwrap().clone();
    assert_eq!(kinds[0], "model-1.0");
    kinds.extend(description["legacy"]["schemas"].as_array().unwrap().clone());
    for kind in kinds {
        let kind = kind.as_str().unwrap();
        let (code, schema) = invoke(&["describe", kind]);
        assert_eq!(code, 0, "{kind}");
        let snapshot: Value = serde_json::from_slice(
            &std::fs::read(root().join(format!("schemas/{kind}.schema.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(schema, snapshot, "{kind}");
    }
    let text = description.to_string();
    for removed in REMOVED {
        assert!(
            !description["commands"]
                .as_array()
                .unwrap()
                .contains(&removed.into()),
            "{removed}"
        );
        assert!(!text.contains(&format!("\"{removed}\"")) || removed == "fsp" || removed == "tau");
    }
    assert!(!text.contains("crn-craft.archive/0.1"));
}

/// A converted immigration-death model with a run/1.0 Static Direct
/// request and an FSP request, written next to each other.
fn unified_inputs(temp: &Temp) -> (Vec<u8>, Value, Value) {
    let network = std::fs::read(root().join("models/legacy/immigration-death.model.json")).unwrap();
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let request = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(&converted)},
        "method": {"kind": "static_direct"},
        "time": {"end": 5.0, "report_times": [0.0, 1.0, 5.0]},
        "ensemble": {"replicates": 16, "root_seed": 9},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 2, "keep_instances": 0, "quantiles": [0.5]},
        "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 10000000,
                   "max_instances": 1000, "threads": 2, "max_output_bytes": 10000000}
    });
    let fsp = serde_json::json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crn_craft::sha256(&converted)},
        "time": {"end": 5.0, "report_times": [0.0, 5.0]},
        "method": {"kind": "fsp", "bounds": {"A": 100}, "max_states": 200, "tolerance": 1e-12},
        "budget": {"max_output_bytes": 10000000}
    });
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    std::fs::write(temp.0.join("fsp.json"), serde_json::to_vec(&fsp).unwrap()).unwrap();
    (converted, request, fsp)
}

#[test]
fn run_and_check_are_rules_and_rules_check() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let (converted, mut request, mut fsp) = unified_inputs(&temp);
    let (model, file, fsp_file) = (path("model.json"), path("request.json"), path("fsp.json"));
    // run and its alias print the same bytes, for simulations and FSP.
    for request_file in [&file, &fsp_file] {
        let (code, by_run) = text(&["run", &model, request_file]);
        let (alias_code, by_alias) = text(&["rules", &model, request_file]);
        assert_eq!((code, alias_code), (0, 0), "{by_run}");
        assert_eq!(by_run, by_alias);
        assert!(by_run.ends_with('\n') && !by_run.trim_end().contains('\n'));
    }
    let (_, result) = invoke(&["run", &model, &file]);
    assert_eq!(result["schema"], "crn-craft.result/1.0");
    assert_eq!(result["statistics"]["quantiles"][0]["p"], 0.5);
    // run --archive is rules --archive.
    let (code, report) = invoke(&["run", &model, &file, "--archive", &path("archive")]);
    assert_eq!(code, 0, "{report}");
    assert_eq!(report["schema"], "crn-craft.archive-check/1.0");
    let (code, verified) = invoke(&["verify", &path("archive")]);
    assert_eq!((code, &verified["valid"]), (0, &Value::Bool(true)));
    // check and its alias: the model alone, then with each request.
    let (code, by_check) = text(&["check", &model]);
    let (_, by_alias) = text(&["rules-check", &model]);
    assert_eq!((code, &by_check), (0, &by_alias));
    for (request_file, bytes) in [
        (&file, serde_json::to_vec(&request).unwrap()),
        (&fsp_file, serde_json::to_vec(&fsp).unwrap()),
    ] {
        let (code, report) = invoke(&["check", &model, request_file]);
        assert_eq!(code, 0, "{report}");
        assert_eq!(
            report,
            serde_json::json!({
                "schema": "crn-craft.check/1.0", "valid": true, "diagnostics": [],
                "model_sha256": crn_craft::sha256(&converted),
                "run_sha256": crn_craft::sha256(&bytes)
            })
        );
        let (_, alias) = invoke(&["rules-check", &model, request_file]);
        assert_eq!(alias, report);
    }
    // check reports what run refuses, before running anything.
    request["budget"]["threads"] = 0.into();
    std::fs::write(
        temp.0.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let (code, checked) = invoke(&["check", &model, &file]);
    let (run_code, refused) = invoke(&["run", &model, &file]);
    assert_eq!((code, run_code), (1, 1));
    assert_eq!(checked["valid"], false);
    assert_eq!(checked["diagnostics"], refused["diagnostics"], "{checked}");
    fsp["method"]["bounds"] = serde_json::json!({});
    std::fs::write(temp.0.join("fsp.json"), serde_json::to_vec(&fsp).unwrap()).unwrap();
    let (code, checked) = invoke(&["check", &model, &fsp_file]);
    let (_, refused) = invoke(&["run", &model, &fsp_file]);
    assert_eq!(code, 1);
    assert_eq!(checked["diagnostics"][0]["path"], "/method/bounds/A");
    assert_eq!(checked["diagnostics"], refused["diagnostics"]);
    // A request that names another file is refused as run refuses it.
    let other = path("other.json");
    std::fs::write(&other, &converted).unwrap();
    let (code, checked) = invoke(&["check", &other, &fsp_file]);
    assert_eq!(code, 1);
    assert!(
        checked["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "CRN-R020" && d["path"] == "/model/path"),
        "{checked}"
    );
    // A request file that cannot be read.
    let (code, checked) = invoke(&["check", &model, &path("missing.json")]);
    assert_eq!(
        (code, &checked["schema"]),
        (1, &"crn-craft.check/1.0".into())
    );
    assert_eq!(checked["diagnostics"][0]["code"], "CRN-E030");
}

#[test]
fn check_reads_every_model_version() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let (_, request, _) = unified_inputs(&temp);
    // A rules model alone: the rules-check/0.1 report of release 0.6.0.
    let rules = "models/rules/immigration-death.rules.json";
    let (code, report) = invoke(&["check", rules]);
    assert_eq!(code, 0);
    let expected = crn_craft::rules::check_model(&std::fs::read(root().join(rules)).unwrap());
    assert_eq!(report, serde_json::to_value(expected).unwrap());
    assert_eq!(report["schema"], "crn-craft.rules-check/0.1");
    assert!(report.get("run_sha256").is_none());
    // A network model alone: check/1.0 of its conversion.
    let network = "models/legacy/repressilator.model.json";
    let bytes = std::fs::read(root().join(network)).unwrap();
    let (code, report) = invoke(&["check", network]);
    assert_eq!(code, 0, "{report}");
    assert_eq!(
        report,
        serde_json::json!({
            "schema": "crn-craft.check/1.0", "valid": true, "diagnostics": [],
            "model_sha256": crn_craft::sha256(&bytes)
        })
    );
    // An invalid network model keeps its validation diagnostics; one the
    // converter cannot hold gets CRN-R016.
    let (code, report) = invoke(&["check", "models/invalid/unknown-species.model.json"]);
    assert_eq!(
        (code, &report["schema"]),
        (1, &"crn-craft.check/1.0".into())
    );
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E012");
    let mut limited: Value = serde_json::from_slice(&bytes).unwrap();
    limited["species"][0]["population_limit"] = 1000.into();
    std::fs::write(temp.0.join("limited.json"), limited.to_string()).unwrap();
    let (code, report) = invoke(&["check", &path("limited.json")]);
    assert_eq!(code, 1);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-R016", "{report}");
    // A network model with a run/1.0 request, and run with a network
    // model: CRN-E033, convert it first.
    std::fs::copy(root().join(network), temp.0.join("model.json")).unwrap();
    let (code, report) = invoke(&["check", &path("model.json"), &path("request.json")]);
    assert_eq!(code, 1);
    assert_eq!(report["schema"], "crn-craft.check/1.0");
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E033");
    assert_eq!(report["diagnostics"][0]["path"], "/schema");
    let (code, report) = invoke(&["run", &path("model.json"), &path("request.json")]);
    assert_eq!(code, 1);
    assert_eq!(report["diagnostics"][0]["code"], "CRN-E033");
    assert_eq!(request["schema"], "crn-craft.run/1.0");
}

#[test]
fn removed_network_commands_and_schemas_are_gone() {
    let temp = staged();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let (model, run) = (
        path("immigration-death.model.json"),
        path("immigration-death.run.json"),
    );
    let network = "models/legacy/immigration-death.model.json";
    let legacy_run = "models/legacy/immigration-death.run.json";
    for command in REMOVED {
        for args in [
            vec![command, network, legacy_run],
            vec![command, network, legacy_run, &path("out")],
            vec![command, &model, &run],
            vec![command, &model, &run, "--cancel-file", &path("stop")],
        ] {
            let (code, report) = invoke(&args);
            assert_eq!(code, 1, "{args:?}");
            assert_eq!(report["diagnostics"][0]["code"], "CRN-E031", "{args:?}");
        }
    }
    assert!(!temp.0.join("out").exists());
    for kind in REMOVED_SCHEMAS {
        let (code, report) = invoke(&["describe", kind]);
        assert_eq!(code, 1, "{kind}");
        assert_eq!(report["diagnostics"][0]["code"], "CRN-E031", "{kind}");
        assert!(
            !root().join(format!("schemas/{kind}.schema.json")).exists(),
            "{kind}"
        );
    }
}
