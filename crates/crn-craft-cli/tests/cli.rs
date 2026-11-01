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
    assert_eq!(value["milestone"], "M10");
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
    assert_eq!(report["schema"], "crn-craft.analysis/0.2");
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
fn legacy_analysis_retains_version_and_common_payload() {
    let path = "models/examples/enzyme-mechanism.model.json";
    let (code, mut current) = invoke(&["analyze", path]);
    assert_eq!(code, 0);
    let (code, legacy) = invoke(&["analyze-v1", path]);
    assert_eq!(code, 0);
    assert_eq!(current["schema"], "crn-craft.analysis/0.2");
    assert_eq!(legacy["schema"], "crn-craft.analysis/0.1");
    let data = current["analysis"].as_object_mut().unwrap();
    for key in [
        "nonnegative_conservation",
        "linkage",
        "initial_class_upper_bounds",
    ] {
        assert!(data.remove(key).is_some());
    }
    assert_eq!(current["analysis"], legacy["analysis"]);
    assert_eq!(current["model_sha256"], legacy["model_sha256"]);
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
