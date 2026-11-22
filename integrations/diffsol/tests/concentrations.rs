use serde_json::{Value, json};
use std::io::Write;
use std::process::{Command, Stdio};

const BINARY: &str = env!("CARGO_BIN_EXE_crn-diffsol");

#[test]
fn concentration_batches_reuse_a_model_and_preserve_amount() {
    let path = std::env::temp_dir().join(format!("crn-concentration-{}.json", std::process::id()));
    let model = json!({
        "schema": "crn-craft.concentration-model/0.1",
        "states": ["A", "B"],
        "parameters": ["k"],
        "initial": ["0.375", "0.125"],
        "assignments": [["flux", "k*2*A"]],
        "derivatives": ["-flux/2", "flux/3"],
        "observations": ["2*A+3*B", "B+time"]
    });
    std::fs::write(&path, model.to_string()).unwrap();
    let mut child = Command::new(BINARY)
        .arg("concentrations")
        .arg(&path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    {
        let mut input = child.stdin.take().unwrap();
        writeln!(input, "{{\"schema\":\"x\",\"schema\":\"y\",\"cases\":[]}}").unwrap();
        let request = |rate| {
            json!({"parameters":[rate], "times":[0.,1.,2.],
                                   "rtol":1e-9, "atol":1e-11})
        };
        writeln!(
            input,
            "{}",
            json!({"schema":"crn-craft.concentration-request/0.1",
                                   "cases":[request(0.4),request(0.9)]})
        )
        .unwrap();
    }
    let output = child.wait_with_output().unwrap();
    std::fs::remove_file(path).unwrap();
    assert!(output.status.success());
    let rows: Vec<Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|row| serde_json::from_str(row).unwrap())
        .collect();
    assert!(rows[0]["error"].is_string());
    assert_eq!(rows[1]["schema"], "crn-craft.concentration-result/0.1");
    for (i, rate) in [0.4_f64, 0.9].into_iter().enumerate() {
        for (j, t) in [0., 1., 2.].into_iter().enumerate() {
            let case = &rows[1]["cases"][i]["values"];
            let a = case["states"][j][0].as_f64().unwrap();
            let b = case["states"][j][1].as_f64().unwrap();
            assert!((a - 0.375 * (-rate * t).exp()).abs() < 1e-8);
            assert!((2. * a + 3. * b - 1.125).abs() < 1e-11);
            assert!((case["observations"][j][1].as_f64().unwrap() - b - t).abs() < 1e-11);
        }
    }
}

#[test]
fn concentration_discovery_matches_checked_in_schemas() {
    let result = Command::new(BINARY)
        .args(["concentrations", "describe"])
        .output()
        .unwrap();
    assert!(result.status.success());
    let description: Value = serde_json::from_slice(&result.stdout).unwrap();
    for (key, source) in [
        (
            "model_schema",
            include_str!("../schemas/concentration-model.schema.json"),
        ),
        (
            "request_schema",
            include_str!("../schemas/concentration-request.schema.json"),
        ),
        (
            "result_schema",
            include_str!("../schemas/concentration-result.schema.json"),
        ),
    ] {
        assert_eq!(
            description[key],
            serde_json::from_str::<Value>(source).unwrap()
        );
    }
}
