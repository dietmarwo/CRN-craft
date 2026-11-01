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
    assert_eq!(value["milestone"], "M0");
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
