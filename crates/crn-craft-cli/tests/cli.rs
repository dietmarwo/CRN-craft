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
    assert_eq!(value["milestone"], "SSA-acceleration");
    // The rules entry gained `hybrid_request` and `hybrid_result` in 0.3.0,
    // so the native discovery moved from 0.11 to 0.12; it gained
    // `static_direct_defaults_request` and `accepted_requests` with
    // rules-run/0.6, so it moved to 0.13. The `convert` command and the
    // rules entry's `convert` key moved it to 0.14; the analysis/1.0 report
    // and the removal of `analyze-v1` (M5b) moved it to 0.15.
    assert_eq!(value["schema"], "crn-craft.describe/0.15");
    assert_eq!(
        value["commands"].as_array().unwrap().last().unwrap(),
        "convert"
    );
    assert_eq!(
        value["dynamic_species"]["convert"]["from"],
        "crn-craft.model/0.1"
    );
    assert_eq!(
        value["dynamic_species"]["hybrid_request"],
        "crn-craft.rules-run/0.5"
    );
    assert_eq!(
        value["dynamic_species"]["hybrid_result"],
        "crn-craft.rules-result/0.3"
    );
    assert_eq!(
        value["dynamic_species"]["static_direct_defaults_request"],
        "crn-craft.rules-run/0.6"
    );
    let requests = value["dynamic_species"]["accepted_requests"]
        .as_array()
        .unwrap();
    assert_eq!(requests.first().unwrap(), "crn-craft.rules-run/0.1");
    assert_eq!(requests.last().unwrap(), "crn-craft.rules-run/0.6");
    let (code, value) = invoke(&[
        "check",
        "models/examples/immigration-death.model.json",
        "models/examples/immigration-death.run.json",
    ]);
    assert_eq!(code, 0);
    assert_eq!(value["valid"], true);
    assert_eq!(value["execution_available"], false);
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

fn staged() -> Temp {
    let temp = Temp::new();
    for name in ["immigration-death.model.json", "immigration-death.run.json"] {
        std::fs::copy(root().join("models/examples").join(name), temp.0.join(name)).unwrap();
    }
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
            .any(|d| d["code"] == "CRN-E022")
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
            .any(|d| d["code"] == "CRN-E022")
    );
}

#[test]
fn analyze_reads_network_rules_and_1_0_models_and_analyze_v1_is_gone() {
    let temp = Temp::new();
    let network = "models/examples/enzyme-mechanism.model.json";
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
        description["dynamic_species"]["contracts_1_0"]["analysis"]["result"],
        "crn-craft.analysis/1.0"
    );
}

#[test]
fn export_cli_is_json_and_rejects_unsupported_or_missing_inputs() {
    let (code, report) = invoke(&[
        "export",
        "models/examples/immigration-death.model.json",
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
    for (args, expected) in [
        (
            vec![
                "export",
                "models/examples/michaelis-menten.model.json",
                "gillespy2",
            ],
            "CRN-E050",
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
                "models/examples/immigration-death.model.json",
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
fn simulate_exits_zero_only_for_a_complete_ensemble_and_verify_for_valid_archives() {
    let temp = staged();
    let model = temp.0.join("immigration-death.model.json");
    let run = temp.0.join("immigration-death.run.json");
    let out = temp.0.join("complete");
    let paths = |p: &PathBuf| p.to_str().unwrap().to_owned();
    let (code, report) = invoke(&["simulate", &paths(&model), &paths(&run), &paths(&out)]);
    assert_eq!(code, 0, "{report}");
    assert_eq!(report["simulation_complete"], true);
    assert_eq!(invoke(&["verify", &paths(&out)]).0, 0);
    // Existing output directories are never reused.
    let (code, report) = invoke(&["simulate", &paths(&model), &paths(&run), &paths(&out)]);
    assert_eq!((code, &report["valid"]), (1, &Value::Bool(false)));

    let mut limited: Value = serde_json::from_slice(&std::fs::read(&run).unwrap()).unwrap();
    limited["budget"]["max_events_per_replicate"] = 1.into();
    let limited_run = temp.0.join("limited.run.json");
    std::fs::write(&limited_run, serde_json::to_vec_pretty(&limited).unwrap()).unwrap();
    let partial = temp.0.join("partial");
    let (code, report) = invoke(&[
        "simulate",
        &paths(&model),
        &paths(&limited_run),
        &paths(&partial),
    ]);
    assert_eq!(code, 1, "{report}");
    assert_eq!(report["valid"], true);
    assert_eq!(report["simulation_complete"], false);
    assert_eq!(report["termination"], "event_limit");
    let (code, report) = invoke(&["verify", &paths(&partial)]);
    assert_eq!(
        (code, &report["simulation_complete"]),
        (0, &Value::Bool(false))
    );
    std::fs::remove_file(partial.join("COMPLETE")).unwrap();
    assert_eq!(invoke(&["verify", &paths(&partial)]).0, 1);
}

#[test]
fn tau_accepts_a_cancel_file_and_stops_before_any_replicate() {
    let dir = std::env::temp_dir().join(format!("crn-cli-tau-cancel-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir(&dir).unwrap();
    let model = br#"{"schema":"crn-craft.model/0.1","name":"T","convention":"stochastic_counts","omega":1.0,"parameters":{},"species":[{"name":"X","initial_count":30}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":20.0}}]}"#;
    std::fs::write(dir.join("model.json"), model).unwrap();
    let run = serde_json::json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(model)},"method":{"kind":"ssa","aggregator":"direct"},"ensemble":{"replicates":4,"root_seed":1},"time":{"end":1.0,"save_times":[0.0,1.0]},"outputs":{"statistics":["mean"],"keep_trajectories":0},"budget":{"threads":1,"max_events_per_replicate":100000,"max_output_bytes":100000}});
    let request = serde_json::json!({"schema":"crn-craft.tau-request/0.1","epsilon":0.05,"max_leaps_per_replicate":1000,"reference_run":run});
    std::fs::write(dir.join("tau.json"), serde_json::to_vec(&request).unwrap()).unwrap();
    std::fs::write(dir.join("cancel"), b"").unwrap();
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_crn-craft"))
        .arg("tau")
        .arg(dir.join("model.json"))
        .arg(dir.join("tau.json"))
        .arg("--cancel-file")
        .arg(dir.join("cancel"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    let value: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(value["termination"], "cancelled");
    assert!(
        value["replicates"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r["started"] == false)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn convert_prints_the_converted_model_bytes_or_coded_refusals() {
    let path = "models/examples/repressilator.model.json";
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
fn network_commands_refuse_1_0_documents_and_rules_commands_run_them() {
    let temp = Temp::new();
    let path = |name: &str| temp.0.join(name).to_str().unwrap().to_string();
    let network =
        std::fs::read(root().join("models/examples/immigration-death.model.json")).unwrap();
    let converted = crn_craft::rules::convert::from_network(&network).unwrap();
    std::fs::write(temp.0.join("model.json"), &converted).unwrap();
    let mut run: Value = serde_json::from_slice(
        &std::fs::read(root().join("models/examples/immigration-death.run.json")).unwrap(),
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
    for args in [
        vec!["check", model],
        vec!["check", model, run],
        vec!["export", model, "gillespy2"],
        vec!["simulate", model, run, out],
        vec!["fsp", model, run],
        vec!["ssa", model, request],
        vec!["tau", model, request],
    ] {
        let (code, report) = invoke(&args);
        assert_eq!(code, 1, "{args:?}");
        let diagnostics = report["diagnostics"].as_array().unwrap();
        assert!(
            diagnostics
                .iter()
                .any(|d| d["code"] == "CRN-E032" && d["path"] == "/schema"),
            "{args:?}: {report}"
        );
    }
    assert!(!temp.0.join("out").exists());
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
    let network =
        std::fs::read(root().join("models/examples/immigration-death.model.json")).unwrap();
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
        description["dynamic_species"]["contracts_1_0"]["fsp"]["result"],
        "crn-craft.fsp-result/1.0"
    );
}
