//! The M5e reference gate, declared in `validation/merge/RESULTS.md`: the
//! kinetics on the flat view against the 0.6.0 kinetics references of M0c
//! (`validation/merge/reference/<name>.kinetics.json`), at the initial and
//! the probe state of each of the twelve merge models, for the network file
//! and for the model converted by this engine.
//!
//! Nonzero references: relative difference at most 1e-12. Zero references:
//! absolute value at most 1e-12 times the largest absolute value of the same
//! array at the same point (exactly 0 when that array is all zeros). With
//! `CRN_KINETICS_REPORT=PATH` the test also writes its report there.
use crn_craft::kinetics::{Derivatives, Kinetics};
use crn_craft::{rules::convert, sha256};
use serde_json::{Value, json};
use std::path::PathBuf;

const RELATIVE: f64 = 1e-12;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[derive(Default)]
struct Run {
    values: usize,
    identical: usize,
    largest: f64,
    zeros: usize,
    zero_exact: usize,
    zero_inexact: Vec<Value>,
    failures: Vec<Value>,
}

impl Run {
    fn fail(&mut self, what: Value) {
        self.failures.push(what);
    }

    fn compare(&mut self, point: &str, array: &str, reference: &[f64], actual: &[f64]) {
        if reference.len() != actual.len() {
            self.fail(
                json!({"point": point, "array": array, "lengths": [reference.len(), actual.len()]}),
            );
            return;
        }
        let scale = reference.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
        for (index, (&r, &x)) in reference.iter().zip(actual).enumerate() {
            self.values += 1;
            if r.to_bits() == x.to_bits() {
                self.identical += 1;
            }
            if r != 0.0 {
                let relative = (x - r).abs() / r.abs();
                self.largest = self.largest.max(relative);
                if !(relative <= RELATIVE) {
                    self.fail(json!({"point": point, "array": array, "index": index,
                                     "reference": r, "value": x, "relative": relative}));
                }
            } else {
                self.zeros += 1;
                if x == 0.0 {
                    self.zero_exact += 1;
                } else {
                    let entry = json!({"point": point, "array": array, "index": index,
                                       "value": x, "array_scale": scale});
                    self.zero_inexact.push(entry.clone());
                    if !(x.abs() <= RELATIVE * scale) {
                        self.fail(entry);
                    }
                }
            }
        }
    }

    fn report(&self, name: &str) -> Value {
        json!({
            "run": name,
            "values": self.values,
            "bit_identical": self.identical,
            "largest_relative_difference": self.largest,
            "zero_references": self.zeros,
            "zero_references_exact": self.zero_exact,
            "zero_references_not_exact": self.zero_inexact,
            "failures": self.failures,
        })
    }
}

fn floats(value: &Value) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn matrix(value: &Value) -> Vec<f64> {
    value.as_array().unwrap().iter().flat_map(floats).collect()
}

fn flat(rows: &[Vec<f64>]) -> Vec<f64> {
    rows.iter().flatten().copied().collect()
}

/// Compares one kinetics object with the stored reference of a model.
fn check(kinetics: &Kinetics, reference: &Value, run: &mut Run) {
    let strings = |v: &Value| -> Vec<String> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect()
    };
    if kinetics.species() != strings(&reference["species"]).as_slice() {
        run.fail(json!({"species": kinetics.species()}));
    }
    if kinetics.parameters != strings(&reference["parameters"]) {
        run.fail(json!({"parameters": kinetics.parameters}));
    }
    let p = floats(&reference["parameter_values"]);
    if kinetics.parameter_values() != p.as_slice() {
        run.fail(json!({"parameter_values": kinetics.parameter_values()}));
    }
    for point in reference["points"].as_array().unwrap() {
        let label = point["label"].as_str().unwrap();
        let y = floats(&point["y"]);
        let result: Result<Derivatives, _> = kinetics.evaluate(&y, &p);
        if point.get("diagnostic").is_some() {
            if !matches!(&result, Err(d) if d.code == "CRN-E080") {
                run.fail(json!({"point": label, "expected": "CRN-E080"}));
            }
            continue;
        }
        let d = match result {
            Ok(d) => d,
            Err(e) => {
                run.fail(json!({"point": label, "diagnostic": e}));
                continue;
            }
        };
        run.compare(label, "rhs", &floats(&point["rhs"]), &d.rhs);
        run.compare(
            label,
            "state_jacobian",
            &matrix(&point["state_jacobian"]),
            &flat(&d.state_jacobian),
        );
        run.compare(
            label,
            "parameter_jacobian",
            &matrix(&point["parameter_jacobian"]),
            &flat(&d.parameter_jacobian),
        );
    }
}

#[test]
fn kinetics_agree_with_the_m0c_references() {
    let root = root();
    let manifest: Value =
        serde_json::from_slice(&std::fs::read(root.join("validation/merge/models.json")).unwrap())
            .unwrap();
    let mut models = Vec::new();
    let (mut values, mut failures, mut acceptance) = (0, 0, 0);
    for row in manifest["rows"].as_array().unwrap() {
        let name = row["name"].as_str().unwrap();
        let path = row["path"].as_str().unwrap();
        let stored: Value = serde_json::from_slice(
            &std::fs::read(root.join(format!("validation/merge/reference/{name}.kinetics.json")))
                .unwrap(),
        )
        .unwrap();
        let reference = &stored["model"];
        let bytes = std::fs::read(root.join(path)).unwrap();
        let mut problems = Vec::new();
        if reference["model"] != path || reference["model_sha256"] != sha256(&bytes) {
            problems.push(json!({"reference": "model path or SHA-256 differs"}));
        }
        let expected = reference["accepted"].as_bool().unwrap();
        let built = Kinetics::from_bytes(&bytes);
        acceptance += 1;
        let codes: Vec<String> = match &built {
            Ok(_) => Vec::new(),
            Err(diagnostics) => diagnostics.iter().map(|d| d.code.clone()).collect(),
        };
        let accepted = built.is_ok();
        if accepted != expected || codes.iter().any(|c| c != "CRN-E080") {
            problems.push(json!({"acceptance": accepted, "codes": codes}));
        }
        let mut runs = Vec::new();
        if let (Ok(kinetics), true) = (&built, expected) {
            let mut network = Run::default();
            check(kinetics, reference, &mut network);
            let converted = convert::from_network(&bytes).unwrap();
            let mut run = Run::default();
            match Kinetics::from_bytes(&converted) {
                Ok(k) => check(&k, reference, &mut run),
                Err(e) => run.fail(json!({"converted": e})),
            }
            for (label, run) in [("network", network), ("converted", run)] {
                values += run.values;
                failures += run.failures.len();
                runs.push(run.report(label));
            }
        }
        failures += problems.len();
        models.push(json!({
            "name": name,
            "reference_accepted": expected,
            "accepted": accepted,
            "refusal_codes": codes,
            "problems": problems,
            "runs": runs,
        }));
    }
    let report = json!({
        "schema": "crn-craft.merge-kinetics-check/0.1",
        "declared_in": "validation/merge/RESULTS.md (M5e)",
        "relative_tolerance": RELATIVE,
        "zero_reference_rule": "|x| <= 1e-12 * largest |reference| of the same array at the same point",
        "values": values,
        "acceptance_checks": acceptance,
        "failures": failures,
        "passed": failures == 0,
        "models": models,
    });
    if let Ok(path) = std::env::var("CRN_KINETICS_REPORT") {
        let mut text = serde_json::to_string_pretty(&report).unwrap();
        text.push('\n');
        std::fs::write(path, text).unwrap();
    }
    println!(
        "kinetics references: {values} values, {acceptance} acceptance checks, {failures} failures"
    );
    assert_eq!(failures, 0, "{report:#}");
}
