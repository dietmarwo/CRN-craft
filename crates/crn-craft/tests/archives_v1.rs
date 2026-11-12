//! `crn-craft.archive/1.0`: the archive tests of `tests/archives.rs`, ported
//! to rules runs, plus the result bytes, finite state projections and the
//! dispatch of `verify_any`.
use crn_craft::{
    archive::{self, AnyReport, v1},
    rules::{Control, Termination},
    sha256,
};
use serde_json::{Value, json};
use std::{
    fs,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

static ID: AtomicU64 = AtomicU64::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "crn-archive-v1-{}-{}",
            std::process::id(),
            ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

const MODEL: &[u8] = b"{\n \"schema\":\"crn-craft.model/1.0\",\"name\":\"archive\",\"parameters\":{\"k\":10.0},\"types\":[{\"name\":\"X\"}],\"initial\":[{\"type\":\"X\",\"count\":3}],\"observables\":{\"X\":{\"type\":\"X\"}},\"rules\":[{\"name\":\"birth\",\"reactants\":[],\"rate\":\"k\",\"products\":[{\"type\":\"X\"}]}]}\n";

fn request(model: &[u8]) -> Value {
    json!({"schema":"crn-craft.run/1.0","model":{"path":"model.json","sha256":sha256(model)},
        "method":{"kind":"static_direct"},
        "ensemble":{"replicates":8,"root_seed":42},"time":{"end":2.0,"report_times":[0.0,1.0,2.0]},
        "outputs":{"statistics":["mean","variance"],"keep_paths":8,"keep_instances":2,"quantiles":[0.25,0.5]},
        "budget":{"max_events_per_replicate":10000,"max_sampling_steps_per_replicate":1000000,
                  "max_instances":100,"threads":2,"max_output_bytes":1000000}})
}

fn inputs(dir: &Path, edit: impl FnOnce(&mut Value)) -> (Vec<u8>, Vec<u8>) {
    let mut run = request(MODEL);
    edit(&mut run);
    let run = serde_json::to_vec_pretty(&run).unwrap();
    fs::write(dir.join("model.json"), MODEL).unwrap();
    fs::write(dir.join("request.json"), &run).unwrap();
    (MODEL.to_vec(), run)
}

fn create(dir: &Path, control: &Control) -> v1::ArchiveReport {
    v1::run_files(
        &dir.join("model.json"),
        &dir.join("request.json"),
        &dir.join("out"),
        control,
    )
}

/// Rewrites the manifest and the marker for the current payload bytes, as a
/// forger would.
fn reseal(dir: &Path) {
    let files: Vec<_> = ["model.json", "request.json", "result.json"]
        .into_iter()
        .map(|name| {
            let b = fs::read(dir.join(name)).unwrap();
            json!({"path":name,"bytes":b.len(),"sha256":sha256(&b)})
        })
        .collect();
    let m = serde_json::to_vec(&json!({"schema":"crn-craft.archive/1.0","files":files})).unwrap();
    fs::write(dir.join("COMPLETE"), format!("{}\n", sha256(&m))).unwrap();
    fs::write(dir.join("manifest.json"), m).unwrap();
}

/// The result as `crn-craft rules` prints it: compact JSON, sorted keys.
fn printed(value: &Value) -> Vec<u8> {
    let mut bytes = serde_json::to_vec(value).unwrap();
    bytes.push(b'\n');
    bytes
}

fn message(report: &v1::ArchiveReport) -> &str {
    report
        .diagnostics
        .first()
        .map_or("", |d| d.message.as_str())
}

#[test]
fn completed_archive_preserves_inputs_and_byte_budget_and_refuses_overwrite() {
    let tmp = Temp::new();
    let (m, r) = inputs(&tmp.0, |_| {});
    let report = create(&tmp.0, &Control::default());
    assert!(
        report.valid && report.archive_complete && report.simulation_complete,
        "{report:?}"
    );
    assert_eq!(report.result.as_deref(), Some("crn-craft.result/1.0"));
    assert_eq!(report.termination, Some(Termination::Completed));
    assert_eq!(report.model_sha256, Some(sha256(&m)));
    assert_eq!(report.run_sha256, Some(sha256(&r)));
    let out = tmp.0.join("out");
    assert_eq!(fs::read(out.join("model.json")).unwrap(), m);
    assert_eq!(fs::read(out.join("request.json")).unwrap(), r);
    let names: Vec<_> = fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    assert_eq!(
        names,
        [
            "COMPLETE",
            "manifest.json",
            "model.json",
            "request.json",
            "result.json"
        ]
    );
    // The result is the engine's result for the same bytes, as printed by
    // `crn-craft rules`.
    let expected = crn_craft::rules::solve(&m, &r, &Control::default()).unwrap();
    let bytes = fs::read(out.join("result.json")).unwrap();
    assert_eq!(bytes, printed(&serde_json::to_value(&expected).unwrap()));
    assert!(bytes.len() <= 1_000_000);
    let verify = v1::verify(&out);
    assert!(verify.valid && verify.simulation_complete, "{verify:?}");
    assert_eq!(verify.run_sha256, Some(sha256(&r)));
    let before = fs::read(out.join("COMPLETE")).unwrap();
    let again = create(&tmp.0, &Control::default());
    assert!(!again.valid);
    assert_eq!(again.diagnostics[0].code, "CRN-E065");
    assert_eq!(before, fs::read(out.join("COMPLETE")).unwrap());
}

#[test]
fn the_byte_api_writes_the_same_archive() {
    let tmp = Temp::new();
    let (m, r) = inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let written = v1::write(&m, &r, &tmp.0.join("bytes"), &Control::default());
    assert!(written.valid && written.simulation_complete, "{written:?}");
    for name in [
        "model.json",
        "request.json",
        "result.json",
        "manifest.json",
        "COMPLETE",
    ] {
        assert_eq!(
            fs::read(tmp.0.join("out").join(name)).unwrap(),
            fs::read(tmp.0.join("bytes").join(name)).unwrap(),
            "{name}"
        );
    }
    // A missing parent is refused without creating anything.
    let nested = tmp.0.join("missing").join("out");
    assert!(!v1::write(&m, &r, &nested, &Control::default()).valid);
    assert!(!tmp.0.join("missing").exists());
}

#[test]
fn interrupted_or_tampered_archives_are_never_accepted() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let marker = fs::read(out.join("COMPLETE")).unwrap();
    fs::remove_file(out.join("COMPLETE")).unwrap();
    assert!(!v1::verify(&out).valid);
    fs::write(out.join("COMPLETE"), marker).unwrap();
    fs::write(out.join("extra"), b"x").unwrap();
    assert!(!v1::verify(&out).valid);
    fs::remove_file(out.join("extra")).unwrap();
    fs::create_dir(out.join("nested")).unwrap();
    assert!(!v1::verify(&out).valid);
    fs::remove_dir(out.join("nested")).unwrap();
    let original = fs::read(out.join("result.json")).unwrap();
    fs::write(out.join("result.json"), b"{}").unwrap();
    assert!(!v1::verify(&out).valid);
    fs::write(out.join("result.json"), original).unwrap();
    assert!(v1::verify(&out).valid);
    let m = fs::read(out.join("manifest.json")).unwrap();
    let mut manifest: Value = serde_json::from_slice(&m).unwrap();
    manifest["files"][0]["path"] = json!("../model.json");
    let bytes = serde_json::to_vec(&manifest).unwrap();
    fs::write(out.join("COMPLETE"), format!("{}\n", sha256(&bytes))).unwrap();
    fs::write(out.join("manifest.json"), bytes).unwrap();
    let report = v1::verify(&out);
    assert!(!report.valid);
    assert_eq!(report.diagnostics[0].code, "CRN-E065");
}

#[test]
fn checksum_and_completion_hash_are_checked_independently_of_json_semantics() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let marker = fs::read(out.join("COMPLETE")).unwrap();
    fs::write(out.join("COMPLETE"), format!("{}\n", "0".repeat(64))).unwrap();
    assert!(!v1::verify(&out).valid);
    fs::write(out.join("COMPLETE"), marker).unwrap();
    // The same JSON value with one string of equal length changed, without
    // resealing: the digest catches it before any JSON is read.
    let result = fs::read_to_string(out.join("result.json")).unwrap();
    let changed = result.replacen("\"birth\"", "\"bxrth\"", 1);
    assert_ne!(result, changed);
    assert_eq!(result.len(), changed.len());
    fs::write(out.join("result.json"), changed).unwrap();
    let report = v1::verify(&out);
    assert!(!report.valid);
    assert_eq!(message(&report), "payload checksum or length mismatch");
}

#[test]
fn refreshed_hashes_do_not_hide_false_accounting_or_statistics() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let original: Value =
        serde_json::from_slice(&fs::read(out.join("result.json")).unwrap()).unwrap();
    let identity = "result identity or dimensions differ from inputs";
    let accounting = "invalid replicate accounting";
    let aggregate = "aggregate accounting mismatch";
    let recomputed = "statistics disagree with the complete kept ensemble";
    type Mutation = Box<dyn Fn(&mut Value)>;
    let mutations: Vec<(&str, &str, Mutation)> = vec![
        (
            "completed count",
            aggregate,
            Box::new(|r| r["completed_replicates"] = json!(7)),
        ),
        (
            "events",
            accounting,
            Box::new(|r| r["replicates"][0]["events"] = json!(u64::MAX)),
        ),
        (
            "saved points",
            "incomplete replicate marked complete",
            Box::new(|r| r["replicates"][0]["saved_points"] = json!(2)),
        ),
        (
            "mean",
            recomputed,
            Box::new(|r| r["statistics"]["mean"][1][0] = json!(999.0)),
        ),
        (
            "variance alone",
            "inconsistent statistics",
            Box::new(|r| r["statistics"]["sample_variance"][2][0] = json!(0.5)),
        ),
        (
            "variance with its standard error",
            recomputed,
            Box::new(|r| {
                r["statistics"]["sample_variance"][2][0] = json!(8.0);
                r["statistics"]["mean_standard_error"][2][0] = json!(1.0);
            }),
        ),
        (
            "standard error",
            "inconsistent statistics",
            Box::new(|r| r["statistics"]["mean_standard_error"][2][0] = json!(0.5)),
        ),
        (
            "decreasing quantile",
            "inconsistent statistics",
            Box::new(|r| r["statistics"]["quantiles"][1]["values"][2][0] = json!(1.0)),
        ),
        (
            "quantile",
            recomputed,
            Box::new(|r| r["statistics"]["quantiles"][1]["values"][2][0] = json!(1e9)),
        ),
        (
            "kept path",
            recomputed,
            Box::new(|r| r["paths"][0]["values"][1][0] = json!(1.0)),
        ),
        (
            "kept path start",
            "kept path starts at the wrong initial state",
            Box::new(|r| r["paths"][3]["values"][0][0] = json!(4.0)),
        ),
        (
            "record id",
            accounting,
            Box::new(|r| r["replicates"][0]["id"] = json!(1)),
        ),
        (
            "sampler",
            identity,
            Box::new(|r| r["sampler"] = json!("unknown")),
        ),
        ("seed", identity, Box::new(|r| r["root_seed"] = json!(43))),
        (
            "run hash",
            identity,
            Box::new(|r| r["run_sha256"] = json!("0".repeat(64))),
        ),
        (
            "rule events",
            accounting,
            Box::new(|r| r["replicates"][0]["rule_events"][0] = json!(0)),
        ),
        (
            "total events",
            aggregate,
            Box::new(|r| {
                let e = r["events"].as_u64().unwrap();
                r["events"] = json!(e + 1)
            }),
        ),
        (
            "termination",
            "termination evidence mismatch",
            Box::new(|r| r["termination"] = json!("failed")),
        ),
        (
            "diagnostics",
            "termination evidence mismatch",
            Box::new(|r| {
                r["diagnostics"] = json!([{"severity": "warning", "code": "CRN-R050",
                    "path": "/replicates", "message": "x", "hint": "y"}])
            }),
        ),
        (
            "stopped count",
            aggregate,
            Box::new(|r| r["stopped_replicates"] = json!(1)),
        ),
        (
            "kept path dropped",
            "kept path identity or length mismatch",
            Box::new(|r| {
                r["paths"].as_array_mut().unwrap().pop();
            }),
        ),
        (
            "final instance count",
            "final instance identity or content mismatch",
            Box::new(|r| r["final_instances"][0]["instances"][0]["count"] = json!(0)),
        ),
        (
            "statistics withheld",
            "statistics present for incomplete work or absent for completed work",
            Box::new(|r| r["statistics"] = Value::Null),
        ),
        (
            "estimator dropped",
            "statistic estimator/sample count mismatch",
            Box::new(|r| r["statistics"]["sample_variance"] = Value::Null),
        ),
        (
            "report time",
            identity,
            Box::new(|r| r["report_times"][1] = json!(1.5)),
        ),
    ];
    for (name, expected, mutation) in mutations {
        let mut result = original.clone();
        mutation(&mut result);
        assert_ne!(result, original, "{name} changed nothing");
        fs::write(out.join("result.json"), printed(&result)).unwrap();
        reseal(&out);
        let report = v1::verify(&out);
        assert!(!report.valid, "mutation {name} survived");
        assert_eq!(report.diagnostics[0].code, "CRN-E065", "{name}");
        assert_eq!(message(&report), expected, "{name}");
    }
    fs::write(out.join("result.json"), printed(&original)).unwrap();
    reseal(&out);
    assert!(v1::verify(&out).valid);
}

#[test]
fn a_result_that_is_not_the_engines_serialization_is_refused() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let original: Value =
        serde_json::from_slice(&fs::read(out.join("result.json")).unwrap()).unwrap();
    // Pretty-printed: the same value, other bytes.
    fs::write(
        out.join("result.json"),
        serde_json::to_vec_pretty(&original).unwrap(),
    )
    .unwrap();
    reseal(&out);
    let report = v1::verify(&out);
    assert!(!report.valid);
    assert_eq!(
        message(&report),
        "result bytes are not the engine's serialization"
    );
    // An unknown field fails the strict parse.
    let mut extra = original.clone();
    extra["note"] = json!("x");
    fs::write(out.join("result.json"), printed(&extra)).unwrap();
    reseal(&out);
    assert_eq!(v1::verify(&out).diagnostics[0].code, "CRN-E003");
}

#[test]
fn duplicate_keys_and_manifest_entries_fail_even_with_matching_marker() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let original = fs::read(out.join("manifest.json")).unwrap();
    let mut m: Value = serde_json::from_slice(&original).unwrap();
    m["files"][1] = m["files"][0].clone();
    let bytes = serde_json::to_vec(&m).unwrap();
    fs::write(out.join("manifest.json"), &bytes).unwrap();
    fs::write(out.join("COMPLETE"), format!("{}\n", sha256(&bytes))).unwrap();
    assert!(!v1::verify(&out).valid);
    let text = String::from_utf8(original.clone()).unwrap().replacen(
        '{',
        "{\"schema\":\"crn-craft.archive/1.0\",",
        1,
    );
    fs::write(out.join("manifest.json"), text.as_bytes()).unwrap();
    fs::write(
        out.join("COMPLETE"),
        format!("{}\n", sha256(text.as_bytes())),
    )
    .unwrap();
    let bad = v1::verify(&out);
    assert!(!bad.valid);
    assert_eq!(bad.diagnostics[0].code, "CRN-E002");
    // The 0.1 manifest version is not a 1.0 archive.
    let old = String::from_utf8(original)
        .unwrap()
        .replace("crn-craft.archive/1.0", "crn-craft.archive/0.1");
    fs::write(out.join("manifest.json"), old.as_bytes()).unwrap();
    fs::write(
        out.join("COMPLETE"),
        format!("{}\n", sha256(old.as_bytes())),
    )
    .unwrap();
    assert_eq!(v1::verify(&out).diagnostics[0].code, "CRN-E003");
}

#[test]
fn cancellation_and_event_limits_have_complete_archives_but_incomplete_ensembles() {
    for cancel in [true, false] {
        let tmp = Temp::new();
        inputs(&tmp.0, |r| {
            r["budget"]["max_events_per_replicate"] = json!(1)
        });
        let control = Control::default();
        if cancel {
            control.cancel();
        }
        let r = create(&tmp.0, &control);
        assert!(
            r.valid && r.archive_complete && !r.simulation_complete,
            "{r:?}"
        );
        assert_eq!(
            r.termination,
            Some(if cancel {
                Termination::Cancelled
            } else {
                Termination::EventLimit
            })
        );
        let verified = v1::verify(&tmp.0.join("out"));
        assert!(
            verified.valid && !verified.simulation_complete,
            "{verified:?}"
        );
        assert_eq!(verified.termination, r.termination);
    }
}

#[cfg(unix)]
#[test]
fn archive_symlinks_are_rejected() {
    use std::os::unix::fs::symlink;
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    fs::remove_file(out.join("model.json")).unwrap();
    symlink(tmp.0.join("model.json"), out.join("model.json")).unwrap();
    assert!(!v1::verify(&out).valid);
    symlink(&out, tmp.0.join("alias")).unwrap();
    assert!(!v1::verify(&tmp.0.join("alias")).valid);
}

#[test]
fn incorrect_path_or_output_reservation_leaves_no_archive() {
    for bad in ["path", "budget", "hash", "syntax"] {
        let tmp = Temp::new();
        inputs(&tmp.0, |r| match bad {
            "path" => r["model"]["path"] = json!("absent.json"),
            "budget" => r["budget"]["max_output_bytes"] = json!(1000),
            "hash" => r["model"]["sha256"] = json!("0".repeat(64)),
            _ => r["unknown"] = json!(1),
        });
        let report = create(&tmp.0, &Control::default());
        assert!(!report.valid, "{bad}");
        assert!(!tmp.0.join("out").exists(), "{bad}");
        let code = report.diagnostics[0].code.as_str();
        assert_eq!(
            code,
            if bad == "syntax" {
                "CRN-E003"
            } else {
                "CRN-R020"
            },
            "{bad}"
        );
    }
}

#[test]
fn oversized_archive_payload_is_rejected_before_reading() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    fs::OpenOptions::new()
        .write(true)
        .open(out.join("result.json"))
        .unwrap()
        .set_len(crn_craft::rules::MAX_RESULT_BYTES + 1)
        .unwrap();
    let report = v1::verify(&out);
    assert!(!report.valid);
    assert_eq!(report.diagnostics[0].code, "CRN-E065");
    assert!(message(&report).contains("oversized"), "{report:?}");
}

#[test]
fn refreshed_hashes_do_not_hide_inconsistent_partial_evidence_or_budget() {
    // Event-limited paths, two of eight kept: no statistics can be
    // recomputed, so each record and path rule must hold on its own.
    let tmp = Temp::new();
    inputs(&tmp.0, |r| {
        r["budget"]["max_events_per_replicate"] = json!(15);
        r["outputs"]["keep_paths"] = json!(2);
    });
    let report = create(&tmp.0, &Control::default());
    assert!(report.valid && !report.simulation_complete, "{report:?}");
    let out = tmp.0.join("out");
    let original = fs::read(out.join("result.json")).unwrap();
    let parsed: Value = serde_json::from_slice(&original).unwrap();
    let records = parsed["replicates"].as_array().unwrap();
    let limited = records
        .iter()
        .position(|r| r["termination"] == "event_limit")
        .expect("an event-limited replicate");
    let budget = 1_000_000_usize;
    let expected = [
        "invalid replicate accounting",
        "failed or halted replicate has no diagnostic",
        "kept path identity or length mismatch",
        "kept path starts at the wrong initial state",
        "result exceeds the request's output budget",
    ];
    for (mutation, text) in expected.into_iter().enumerate() {
        let mut result = parsed.clone();
        let mut bytes = match mutation {
            0 => {
                result["replicates"][limited]["started"] = json!(false);
                printed(&result)
            }
            1 => {
                result["replicates"][limited]["diagnostic"] = Value::Null;
                printed(&result)
            }
            2 => {
                result["paths"][1]["values"]
                    .as_array_mut()
                    .unwrap()
                    .push(json!([1.0]));
                printed(&result)
            }
            3 => {
                result["paths"][0]["values"][0][0] = json!(1.0);
                printed(&result)
            }
            _ => original.clone(),
        };
        if mutation == 4 {
            // The same JSON value, but longer than the request's budget.
            bytes.pop();
            bytes.resize(budget, b' ');
            bytes.push(b'\n');
        }
        fs::write(out.join("result.json"), &bytes).unwrap();
        reseal(&out);
        let verified = v1::verify(&out);
        assert!(
            !verified.valid && message(&verified) == text,
            "mutation {mutation}: {verified:?}"
        );
    }
    fs::write(out.join("result.json"), &original).unwrap();
    reseal(&out);
    assert!(v1::verify(&out).valid);
}

#[test]
fn large_report_grids_with_many_threads_remain_verifiable_when_cancelled() {
    let tmp = Temp::new();
    inputs(&tmp.0, |run| {
        run["time"] = json!({"end": 99999.0,
            "report_times": (0..100_000).map(f64::from).collect::<Vec<_>>()});
        run["outputs"]["keep_paths"] = json!(0);
        run["outputs"]["keep_instances"] = json!(0);
        run["budget"]["threads"] = json!(8);
        run["budget"]["max_output_bytes"] = json!(20_000_000);
    });
    let control = Control::default();
    control.cancel();
    let created = create(&tmp.0, &control);
    assert!(created.valid, "{:?}", created.diagnostics);
    assert!(!created.simulation_complete);
    assert!(v1::verify(&tmp.0.join("out")).valid);
}

#[test]
fn every_method_and_rules_run_version_archives_and_verifies() {
    let methods = [
        json!({"kind": "tree"}),
        json!({"kind": "static_direct", "waiting_time": "inverse"}),
        json!({"kind": "next_reaction"}),
        json!({"kind": "tau", "epsilon": 0.05}),
        json!({"kind": "binomial", "step": 0.1}),
        json!({"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}),
    ];
    for method in methods {
        for version in ["crn-craft.run/1.0", "crn-craft.rules-run/0.6"] {
            let tmp = Temp::new();
            // Single-parent unit births and deaths, which binomial accepts.
            let model = r#"{"schema":"crn-craft.rules-model/0.4","name":"branching","parameters":{},
                "types":[{"name":"X"}],"initial":[{"type":"X","count":5}],"observables":{"X":{"type":"X"}},
                "rules":[{"name":"birth","reactants":[{"label":"x","type":"X"}],"rate":"0.6","products":[{"type":"X"},{"type":"X"}]},
                         {"name":"death","reactants":[{"label":"x","type":"X"}],"rate":"0.5","products":[]}]}"#;
            let mut run = request(model.as_bytes());
            run["schema"] = json!(version);
            run["method"] = method.clone();
            if version != "crn-craft.run/1.0" {
                run["outputs"].as_object_mut().unwrap().remove("quantiles");
            }
            fs::write(tmp.0.join("model.json"), model).unwrap();
            fs::write(
                tmp.0.join("request.json"),
                serde_json::to_vec(&run).unwrap(),
            )
            .unwrap();
            let report = create(&tmp.0, &Control::default());
            assert!(
                report.valid && report.simulation_complete,
                "{method} {version}: {report:?}"
            );
            let verified = v1::verify(&tmp.0.join("out"));
            assert!(verified.valid, "{method} {version}: {verified:?}");
            assert_eq!(verified.result, report.result);
        }
    }
}

#[test]
fn finite_state_projections_archive_and_verify() {
    let tmp = Temp::new();
    let model = String::from_utf8(MODEL.to_vec()).unwrap().replace(
        "\"rules\":[",
        "\"rules\":[{\"name\":\"death\",\"reactants\":[{\"label\":\"x\",\"type\":\"X\"}],\"rate\":\"1.0\",\"products\":[]},",
    );
    let run = json!({"schema":"crn-craft.run/1.0","model":{"path":"model.json","sha256":sha256(model.as_bytes())},
        "time":{"end":1.0,"report_times":[0.0,0.5,1.0]},
        "method":{"kind":"fsp","bounds":{"X":60},"max_states":1000,"tolerance":1e-12},
        "budget":{"max_output_bytes":1000000}});
    fs::write(tmp.0.join("model.json"), &model).unwrap();
    fs::write(
        tmp.0.join("request.json"),
        serde_json::to_vec(&run).unwrap(),
    )
    .unwrap();
    let report = create(&tmp.0, &Control::default());
    assert!(report.valid && report.simulation_complete, "{report:?}");
    assert_eq!(report.result.as_deref(), Some("crn-craft.fsp-result/1.0"));
    assert_eq!(report.termination, None);
    let out = tmp.0.join("out");
    let expected = crn_craft::fsp::v1::solve(
        model.as_bytes(),
        &fs::read(tmp.0.join("request.json")).unwrap(),
        &Control::default(),
    )
    .unwrap();
    assert_eq!(
        fs::read(out.join("result.json")).unwrap(),
        printed(&serde_json::to_value(&expected).unwrap())
    );
    assert!(v1::verify(&out).valid);
    let original: Value =
        serde_json::from_slice(&fs::read(out.join("result.json")).unwrap()).unwrap();
    type Mutation = Box<dyn Fn(&mut Value)>;
    let mutations: Vec<Mutation> = vec![
        Box::new(|r| r["species"][0] = json!("Y")),
        Box::new(|r| r["states"][0][0] = json!(4)),
        Box::new(|r| r["snapshots"][1]["time"] = json!(0.25)),
        Box::new(|r| {
            r["snapshots"].as_array_mut().unwrap().pop();
        }),
        Box::new(|r| {
            r["snapshots"][0]["probabilities"]
                .as_array_mut()
                .unwrap()
                .pop();
        }),
        Box::new(|r| r["model_sha256"] = json!("0".repeat(64))),
    ];
    for (i, mutation) in mutations.into_iter().enumerate() {
        let mut result = original.clone();
        mutation(&mut result);
        fs::write(out.join("result.json"), printed(&result)).unwrap();
        reseal(&out);
        assert!(!v1::verify(&out).valid, "fsp mutation {i} survived");
    }
    // A failed projection (bounds below the initial count) leaves no
    // directory.
    let mut bad = run.clone();
    bad["method"]["bounds"]["X"] = json!(1);
    fs::write(
        tmp.0.join("request.json"),
        serde_json::to_vec(&bad).unwrap(),
    )
    .unwrap();
    let refused = v1::run_files(
        &tmp.0.join("model.json"),
        &tmp.0.join("request.json"),
        &tmp.0.join("refused"),
        &Control::default(),
    );
    assert!(!refused.valid);
    assert!(!tmp.0.join("refused").exists());
}

#[test]
fn verify_any_dispatches_on_the_archive_version() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    assert!(v1::is_archive(&out));
    let report = archive::verify_any(&out);
    assert!(report.valid());
    assert!(matches!(report, AnyReport::V10(_)));
    assert_eq!(
        serde_json::to_value(&report).unwrap()["schema"],
        "crn-craft.archive-check/1.0"
    );
    // Without a readable manifest, the 1.0 file names still select 1.0.
    fs::remove_file(out.join("manifest.json")).unwrap();
    assert!(v1::is_archive(&out));
    assert!(matches!(archive::verify_any(&out), AnyReport::V10(r) if !r.valid));

    // A 0.1 archive of the network engine keeps its verifier and report.
    let network = br#"{"schema": "crn-craft.model/0.1", "name": "decay",
      "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
      "species": [{"name": "A", "initial_count": 5}],
      "reactions": [{"name": "death", "reactants": {"A": 1}, "products": {},
                     "rate": {"kind": "mass_action", "coefficient": 1.0}}]}"#;
    let run = json!({"schema": "crn-craft.run/0.1",
      "model": {"path": "network.json", "sha256": sha256(network)},
      "method": {"kind": "ssa", "aggregator": "direct"},
      "time": {"end": 1.0, "save_times": [0.0, 1.0]},
      "ensemble": {"replicates": 4, "root_seed": 3},
      "outputs": {"statistics": ["mean"], "keep_trajectories": 0},
      "budget": {"max_events_per_replicate": 100, "threads": 1, "max_output_bytes": 100000}});
    fs::write(tmp.0.join("network.json"), network).unwrap();
    fs::write(tmp.0.join("run.json"), serde_json::to_vec(&run).unwrap()).unwrap();
    let old = tmp.0.join("old");
    let written = archive::simulate_files(
        &tmp.0.join("network.json"),
        &tmp.0.join("run.json"),
        &old,
        &crn_craft::simulation::Control::default(),
    );
    assert!(written.valid);
    assert!(!v1::is_archive(&old));
    let report = archive::verify_any(&old);
    assert!(report.valid());
    assert!(matches!(report, AnyReport::V01(_)));
    assert_eq!(
        serde_json::to_value(&report).unwrap()["schema"],
        "crn-craft.archive-check/0.1"
    );
    // Neither version: the 0.1 verifier reports it.
    assert!(!archive::verify_any(&tmp.0.join("absent")).valid());
}

#[test]
fn the_archive_schemas_describe_the_written_documents() {
    let manifest = crn_craft::schema("archive-1.0").unwrap();
    assert_eq!(
        manifest["definitions"]["ArchiveVersion"]["enum"][0],
        "crn-craft.archive/1.0"
    );
    let report = crn_craft::schema("archive-check-1.0").unwrap();
    assert_eq!(
        report["definitions"]["CheckVersion"]["enum"][0],
        "crn-craft.archive-check/1.0"
    );
    let entry = crn_craft::rules::describe();
    assert_eq!(
        entry["contracts_1_0"]["archive"]["manifest"],
        "crn-craft.archive/1.0"
    );
    assert_eq!(
        entry["contracts_1_0"]["archive"]["report"],
        "crn-craft.archive-check/1.0"
    );
    for (kind, file) in [
        ("archive-1.0", "archive-1.0.schema.json"),
        ("archive-check-1.0", "archive-check-1.0.schema.json"),
    ] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../schemas")
            .join(file);
        let committed: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        assert_eq!(committed, crn_craft::schema(kind).unwrap(), "{kind}");
    }
}
