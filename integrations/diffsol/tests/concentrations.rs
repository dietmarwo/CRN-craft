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
    assert_eq!(rows[0]["error"]["code"], "CRN-E002");
    assert_eq!(rows[0]["schema"], "crn-craft.concentration-result/0.1");
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

/// Run the worker on `model` with `lines` on stdin; parsed result lines.
fn exchange(model: &Value, lines: &[String]) -> Vec<Value> {
    let path = std::env::temp_dir().join(format!(
        "crn-concentration-{}-{}.json",
        std::process::id(),
        lines.len()
    ));
    std::fs::write(&path, model.to_string()).unwrap();
    let mut child = Command::new(BINARY)
        .arg("concentrations")
        .arg(&path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    let text = lines.join("\n") + "\n";
    // Write from a thread so that a large output cannot block the input.
    let writer = std::thread::spawn(move || input.write_all(text.as_bytes()).unwrap());
    let output = child.wait_with_output().unwrap();
    writer.join().unwrap();
    std::fs::remove_file(path).unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|row| serde_json::from_str(row).unwrap())
        .collect()
}

#[test]
fn failures_are_coded_data_and_the_worker_continues() {
    let model = json!({
        "schema": "crn-craft.concentration-model/0.1",
        "states": ["A"], "parameters": ["k"], "initial": ["1"],
        "derivatives": ["-k*A"], "observations": ["A/k"]
    });
    let case = |k: f64, times: Value, rtol: f64| json!({"parameters": [k], "times": times, "rtol": rtol, "atol": 1e-12});
    let request = |cases: Vec<Value>| {
        json!({"schema": "crn-craft.concentration-request/0.1", "cases": cases}).to_string()
    };
    let good = case(0.5, json!([0.0, 1.0]), 1e-9);
    let rows = exchange(
        &model,
        &[
            "not json".into(),
            request(vec![]),
            request(vec![good.clone(); 33]),
            "x".repeat((1 << 20) + 10),
            json!({"schema": "crn-craft.concentration-request/0.1",
                   "cases": [good.clone()], "extra": 1})
            .to_string(),
            request(vec![
                good.clone(),
                case(0.5, json!([0.0, 1.0, 1.0]), 1e-9),
                case(0.5, json!([0.0, 1.0]), 0.5),
                case(0.0, json!([0.0, 1.0]), 1e-9),
            ]),
            request(vec![good.clone()]),
        ],
    );
    let codes: Vec<_> = rows[..5]
        .iter()
        .map(|r| r["error"]["code"].clone())
        .collect();
    assert_eq!(
        codes,
        ["CRN-E001", "CRN-E004", "CRN-E004", "CRN-E004", "CRN-E003"]
    );
    let cases = &rows[5]["cases"];
    assert!((cases[0]["values"]["states"][1][0].as_f64().unwrap() - (-0.5_f64).exp()).abs() < 1e-8);
    assert_eq!(cases[1]["error"]["code"], "CRN-E015");
    assert_eq!(cases[2]["error"]["code"], "CRN-E015");
    assert_eq!(cases[3]["error"]["code"], "CRN-E080");
    // The same request gives the same bytes.
    assert_eq!(rows[6]["cases"][0], cases[0]);
}

#[test]
fn time_dependent_input_and_equilibrium() {
    let model = json!({
        "schema": "crn-craft.concentration-model/0.1",
        "states": ["X", "A", "B"], "parameters": ["k", "kf", "kr"],
        "initial": ["0", "1", "0"],
        "derivatives": ["1 + time - k*X", "-kf*A + kr*B", "kf*A - kr*B"]
    });
    let rows = exchange(
        &model,
        &[json!({"schema": "crn-craft.concentration-request/0.1",
                 "cases": [{"parameters": [0.5, 2.0, 0.5], "times": [0.0, 2.0, 20.0],
                            "rtol": 1e-11, "atol": 1e-13}]})
        .to_string()],
    );
    let states = &rows[0]["cases"][0]["values"]["states"];
    for (j, t) in [(1, 2.0_f64), (2, 20.0)] {
        // X(t) = (1 + t)/k - 1/k^2 + (1/k^2 - 1/k) exp(-k t), k = 0.5.
        let exact = (1.0 + t) / 0.5 - 4.0 + (4.0 - 2.0) * (-0.5 * t).exp();
        assert!((states[j][0].as_f64().unwrap() - exact).abs() < 1e-8);
        let a = 0.2 + 0.8 * (-2.5 * t).exp();
        assert!((states[j][1].as_f64().unwrap() - a).abs() < 1e-8);
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
