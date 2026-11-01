use crn_craft::{
    archive, sha256,
    simulation::{Control, Stop},
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
            "crn-m2-{}-{}",
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

fn inputs(dir: &Path, edit: impl FnOnce(&mut Value)) -> (Vec<u8>, Vec<u8>) {
    let model=b"{\n \"schema\":\"crn-craft.model/0.1\",\"name\":\"archive\",\"convention\":\"stochastic_counts\",\"omega\":1.0,\"parameters\":{},\"species\":[{\"name\":\"X\",\"initial_count\":0}],\"reactions\":[{\"name\":\"birth\",\"reactants\":{},\"products\":{\"X\":1},\"rate\":{\"kind\":\"mass_action\",\"coefficient\":10.0}}]}\n".to_vec();
    let mut run = json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
        "ensemble":{"replicates":8,"root_seed":42},"time":{"end":2.0,"save_times":[0.0,1.0,2.0]},"outputs":{"statistics":["mean","variance"],"keep_trajectories":8},
        "budget":{"max_events_per_replicate":10000,"threads":2,"max_output_bytes":1000000}});
    edit(&mut run);
    let run = serde_json::to_vec_pretty(&run).unwrap();
    fs::write(dir.join("model.json"), &model).unwrap();
    fs::write(dir.join("run.json"), &run).unwrap();
    (model, run)
}
fn create(dir: &Path, control: &Control) -> archive::ArchiveReport {
    archive::simulate_files(
        &dir.join("model.json"),
        &dir.join("run.json"),
        &dir.join("out"),
        control,
    )
}
fn reseal(dir: &Path) {
    let files: Vec<_> = ["model.json", "run.json", "result.json", "provenance.json"]
        .into_iter()
        .map(|name| {
            let b = fs::read(dir.join(name)).unwrap();
            json!({"path":name,"bytes":b.len(),"sha256":sha256(&b)})
        })
        .collect();
    let m = serde_json::to_vec(&json!({"schema":"crn-craft.archive/0.1","files":files})).unwrap();
    fs::write(dir.join("COMPLETE"), format!("{}\n", sha256(&m))).unwrap();
    fs::write(dir.join("manifest.json"), m).unwrap();
}

#[test]
fn completed_archive_preserves_inputs_and_byte_budget_and_refuses_overwrite() {
    let tmp = Temp::new();
    let (m, r) = inputs(&tmp.0, |_| {});
    let report = create(&tmp.0, &Control::default());
    assert!(report.valid && report.simulation_complete, "{report:?}");
    let out = tmp.0.join("out");
    assert_eq!(fs::read(out.join("model.json")).unwrap(), m);
    assert_eq!(fs::read(out.join("run.json")).unwrap(), r);
    let verify = archive::verify_archive(&out);
    assert!(verify.valid && verify.simulation_complete, "{verify:?}");
    let size: u64 = fs::read_dir(&out)
        .unwrap()
        .map(|p| p.unwrap().metadata().unwrap().len())
        .sum();
    assert!(size <= 1_000_000);
    let before = fs::read(out.join("COMPLETE")).unwrap();
    assert!(!create(&tmp.0, &Control::default()).valid);
    assert_eq!(before, fs::read(out.join("COMPLETE")).unwrap());
}
#[test]
fn interrupted_or_tampered_archives_are_never_accepted() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let marker = fs::read(out.join("COMPLETE")).unwrap();
    fs::remove_file(out.join("COMPLETE")).unwrap();
    assert!(!archive::verify_archive(&out).valid);
    fs::write(out.join("COMPLETE"), marker).unwrap();
    fs::write(out.join("extra"), b"x").unwrap();
    assert!(!archive::verify_archive(&out).valid);
    fs::remove_file(out.join("extra")).unwrap();
    let original = fs::read(out.join("result.json")).unwrap();
    fs::write(out.join("result.json"), b"{}").unwrap();
    assert!(!archive::verify_archive(&out).valid);
    fs::write(out.join("result.json"), original).unwrap();
    let m = fs::read(out.join("manifest.json")).unwrap();
    let mut manifest: Value = serde_json::from_slice(&m).unwrap();
    manifest["files"][0]["path"] = json!("../model.json");
    let bytes = serde_json::to_vec(&manifest).unwrap();
    fs::write(out.join("COMPLETE"), format!("{}\n", sha256(&bytes))).unwrap();
    fs::write(out.join("manifest.json"), bytes).unwrap();
    assert!(!archive::verify_archive(&out).valid);
}

#[test]
fn checksum_and_completion_hash_are_checked_independently_of_json_semantics() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let marker = fs::read(out.join("COMPLETE")).unwrap();
    fs::write(out.join("COMPLETE"), format!("{}\n", "0".repeat(64))).unwrap();
    assert!(!archive::verify_archive(&out).valid);
    fs::write(out.join("COMPLETE"), marker).unwrap();
    let provenance = fs::read_to_string(out.join("provenance.json")).unwrap();
    let changed = provenance.replace("0.1.0-dev.4", "0.1.0-dev.9");
    assert_ne!(provenance, changed);
    assert_eq!(provenance.len(), changed.len());
    fs::write(out.join("provenance.json"), changed).unwrap();
    assert!(!archive::verify_archive(&out).valid);
}
#[test]
fn refreshed_hashes_do_not_hide_false_accounting_or_statistics() {
    let tmp = Temp::new();
    inputs(&tmp.0, |_| {});
    assert!(create(&tmp.0, &Control::default()).valid);
    let out = tmp.0.join("out");
    let original: Value =
        serde_json::from_slice(&fs::read(out.join("result.json")).unwrap()).unwrap();
    for mutation in 0..7 {
        let mut result = original.clone();
        match mutation {
            0 => result["completed_replicates"] = json!(7),
            1 => result["replicates"][0]["events"] = json!(u64::MAX),
            2 => result["replicates"][0]["saved_points"] = json!(2),
            3 => result["statistics"]["mean"][1][0] = json!(999.0),
            4 => result["trajectories"][0]["counts"][0][0] = json!(1),
            5 => result["replicates"][0]["id"] = json!(1),
            _ => result["sampler"] = json!("unknown"),
        }
        fs::write(
            out.join("result.json"),
            serde_json::to_vec(&result).unwrap(),
        )
        .unwrap();
        reseal(&out);
        assert!(
            !archive::verify_archive(&out).valid,
            "mutation {mutation} survived"
        );
    }
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
    assert!(!archive::verify_archive(&out).valid);
    let text = String::from_utf8(original).unwrap().replacen(
        '{',
        "{\"schema\":\"crn-craft.archive/0.1\",",
        1,
    );
    fs::write(out.join("manifest.json"), text.as_bytes()).unwrap();
    fs::write(
        out.join("COMPLETE"),
        format!("{}\n", sha256(text.as_bytes())),
    )
    .unwrap();
    let bad = archive::verify_archive(&out);
    assert!(!bad.valid);
    assert_eq!(bad.diagnostics[0].code, "CRN-E002");
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
                Stop::Cancelled
            } else {
                Stop::EventLimit
            })
        );
        let verified = archive::verify_archive(&tmp.0.join("out"));
        assert!(
            verified.valid && !verified.simulation_complete,
            "{verified:?}"
        );
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
    assert!(!archive::verify_archive(&out).valid);
    symlink(&out, tmp.0.join("alias")).unwrap();
    assert!(!archive::verify_archive(&tmp.0.join("alias")).valid);
}
#[test]
fn incorrect_path_or_output_reservation_leaves_no_archive() {
    for bad_path in [true, false] {
        let tmp = Temp::new();
        inputs(&tmp.0, |r| {
            if bad_path {
                r["model"]["path"] = json!("absent.json");
            } else {
                r["budget"]["max_output_bytes"] = json!(1000);
            }
        });
        assert!(!create(&tmp.0, &Control::default()).valid);
        assert!(!tmp.0.join("out").exists());
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
        .set_len(crn_craft::simulation::MAX_RESULT_BYTES + 1)
        .unwrap();
    let report = archive::verify_archive(&out);
    assert!(!report.valid);
    assert_eq!(report.diagnostics[0].code, "CRN-E065");
    assert!(report.diagnostics[0].message.contains("oversized"));
}
