//! The network model `crn-craft.model/0.1`, read only by the converter since
//! 0.7.0: strict parsing, validation codes and paths, and the fixtures of
//! `models/invalid`. The cases of the network engine's `tests/contracts.rs`
//! that concern the model; its `run/0.1` cases went with that contract.
use crn_craft::diagnostic::Diagnostic;
use crn_craft::rules::convert::from_network;
use crn_craft::{check_any, sha256};
use serde_json::{Value, json};
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn bytes() -> Vec<u8> {
    std::fs::read(root().join("models/legacy/immigration-death.model.json")).unwrap()
}
fn model_value() -> Value {
    serde_json::from_slice(&bytes()).unwrap()
}
fn has(diagnostics: &[Diagnostic], code: &str, path: &str) -> bool {
    diagnostics.iter().any(|x| x.code == code && x.path == path)
}
/// The converter's diagnostics for model bytes it refuses.
fn refused(bytes: &[u8]) -> Vec<Diagnostic> {
    from_network(bytes).expect_err("the converter must refuse the model")
}
/// The first diagnostic of the network reader (strict parsing).
fn parse_error(bytes: &[u8]) -> Diagnostic {
    refused(bytes).remove(0)
}
fn mutated_model(edit: impl FnOnce(&mut Value), code: &str, path: &str) {
    let mut value = model_value();
    edit(&mut value);
    let model = serde_json::to_vec(&value).unwrap();
    let diagnostics = refused(&model);
    assert!(has(&diagnostics, code, path), "{diagnostics:?}");
    // `crn-craft check MODEL` reports the same diagnostics (check/1.0) for a
    // network model; it reads another `schema` as a model of that version.
    if value["schema"] != "crn-craft.model/0.1" {
        return;
    }
    let report = check_any(&model, None);
    assert!(!report.valid);
    assert_eq!(
        serde_json::to_value(&report.diagnostics).unwrap(),
        serde_json::to_value(&diagnostics).unwrap()
    );
}

#[test]
fn negative_fixture_manifest_pins_codes_and_paths() {
    let base = root().join("models/invalid");
    let entries: Vec<Value> =
        serde_json::from_slice(&std::fs::read(base.join("expected.json")).unwrap()).unwrap();
    assert!(entries.len() >= 8);
    for entry in entries {
        let model = std::fs::read(base.join(entry["model"].as_str().unwrap())).unwrap();
        let diagnostics = refused(&model);
        assert!(
            has(
                &diagnostics,
                entry["code"].as_str().unwrap(),
                entry["path"].as_str().unwrap()
            ),
            "{entry}: {diagnostics:?}"
        );
    }
}

#[test]
fn byte_hash_known_answer() {
    assert_eq!(
        sha256(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

/// A request binds the exact model bytes: the same converted model,
/// formatted compactly, is another model.
#[test]
fn formatting_changes_model_identity() {
    let model = from_network(&bytes()).unwrap();
    let request = json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": sha256(&model)},
        "method": {"kind": "static_direct"},
        "time": {"end": 1.0, "report_times": [0.0, 1.0]},
        "ensemble": {"replicates": 2, "root_seed": 1},
        "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1000, "max_sampling_steps_per_replicate": 100000,
                   "max_instances": 1000, "threads": 1, "max_output_bytes": 100000}
    });
    let request = serde_json::to_vec(&request).unwrap();
    assert!(crn_craft::rules::prepare(&model, &request).is_ok());
    let value: Value = serde_json::from_slice(&model).unwrap();
    let compact = serde_json::to_vec(&value).unwrap();
    assert_ne!(compact, model);
    let errors = crn_craft::rules::prepare(&compact, &request).err().unwrap();
    assert!(has(&errors, "CRN-R020", "/model/sha256"), "{errors:?}");
}

/// The shortest round-trip decimal of a double must parse back to the same
/// double, as Rust's `str::parse` and Python's `json` do. serde_json's default
/// best-effort parser misreads about a quarter of such literals by one ulp.
/// The converter reads them and writes them again as parameters.
#[test]
fn decimal_literals_parse_to_correctly_rounded_doubles() {
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut values = vec![3.983_236_103_458_967e-14, 3.167_810_103_450_541e-16];
    while values.len() < 10_000 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let value = f64::from_bits(state & 0x7fff_ffff_ffff_ffff);
        if value.is_finite() {
            values.push(value);
        }
    }
    let entries: Vec<String> = values
        .iter()
        .enumerate()
        .map(|(i, value)| format!("\"p{i}\": {value:e}"))
        .collect();
    let original = String::from_utf8(bytes()).unwrap();
    let text = original.replacen(
        "\"birth\": 4.0,",
        &format!("\"birth\": 4.0, {},", entries.join(", ")),
        1,
    );
    assert_ne!(text, original);
    let converted = from_network(text.as_bytes()).unwrap();
    let parsed: Value = crn_craft::strict_json::parse(&converted).unwrap();
    let parameter = |i: usize| parsed["parameters"][format!("p{i}")].as_f64().unwrap();
    let wrong: Vec<String> = values
        .iter()
        .enumerate()
        .filter(|(i, value)| parameter(*i).to_bits() != value.to_bits())
        .map(|(i, value)| format!("{value:e} -> {:e}", parameter(i)))
        .collect();
    assert!(
        wrong.is_empty(),
        "{} of {} literals misparsed, e.g. {:?}",
        wrong.len(),
        values.len(),
        &wrong[..wrong.len().min(3)]
    );
}

#[test]
fn i64_overflow_and_fractional_counts_rejected() {
    for number in [json!(9223372036854775808_u64), json!(1.5)] {
        mutated_model(
            |m| m["species"][0]["initial_count"] = number,
            "CRN-E003",
            "/species/0/initial_count",
        );
    }
}

#[test]
fn unknown_fields_rejected_in_every_model_object() {
    for path in ["", "/species/0", "/reactions/0", "/reactions/0/rate"] {
        let mut value = model_value();
        value
            .pointer_mut(path)
            .unwrap()
            .as_object_mut()
            .unwrap()
            .insert("typo".into(), json!(1));
        let diagnostics = refused(&serde_json::to_vec(&value).unwrap());
        assert!(
            diagnostics.iter().any(|d| d.code == "CRN-E003"),
            "{path}: {diagnostics:?}"
        );
    }
}

#[test]
fn duplicate_keys_at_all_depths_and_escaped_spellings_rejected() {
    for (text, path) in [
        (r#"{"schema":0,"schema":1}"#, "/schema"),
        (r#"{"parameters":{"a":1,"a":2}}"#, "/parameters/a"),
        (
            r#"{"reactions":[{"products":{"A":1,"A":2}}]}"#,
            "/reactions/0/products/A",
        ),
        (r#"{"a/b~c":1,"a/b~c":2}"#, "/a~1b~0c"),
    ] {
        let error = parse_error(text.as_bytes());
        assert_eq!(error.code, "CRN-E002");
        assert_eq!(error.path, path);
    }
}

#[test]
fn excessive_depth_is_rejected_without_stack_overflow() {
    let raw = format!("{}0{}", "[".repeat(1000), "]".repeat(1000));
    assert_eq!(parse_error(raw.as_bytes()).code, "CRN-E004");
}

#[test]
fn nesting_limit_boundary_matches_its_documentation() {
    // docs/contracts.md: nesting depth at most 48, root at depth 0.
    let nested =
        |levels: usize, inner: &str| format!("{}{inner}{}", "[".repeat(levels), "]".repeat(levels));
    // Within the limit, the document fails only as a model (an array).
    for raw in [nested(48, "0"), nested(49, "")] {
        assert_eq!(parse_error(raw.as_bytes()).code, "CRN-E003");
    }
    let error = parse_error(nested(49, "0").as_bytes());
    assert_eq!(error.code, "CRN-E004");
    assert_eq!(error.path, "/0".repeat(49));
    assert_eq!(
        error.hint,
        "Keep every value within 48 levels below the root, which is level 0."
    );
}

#[test]
fn oversized_document_rejected_before_parsing() {
    let bytes = vec![b' '; crn_craft::strict_json::MAX_DOCUMENT_BYTES + 1];
    assert_eq!(parse_error(&bytes).code, "CRN-E004");
}

#[test]
fn malformed_trailing_non_json_and_non_utf8_rejected() {
    for raw in [
        b"{} {}".as_slice(),
        b"NaN",
        b"Infinity",
        b"\xff",
        b"{\"a\":",
    ] {
        assert_eq!(parse_error(raw).code, "CRN-E001");
    }
}

#[test]
fn wrong_version_is_not_silently_upgraded() {
    mutated_model(
        |m| m["schema"] = json!("crn-craft.model/99"),
        "CRN-E003",
        "/schema",
    );
}

#[test]
fn invalid_identifier_and_namespace_collision() {
    mutated_model(
        |m| m["species"][0]["name"] = json!("bad name"),
        "CRN-E010",
        "/species/0/name",
    );
    mutated_model(
        |m| m["parameters"]["A"] = json!(1),
        "CRN-E011",
        "/parameters/A",
    );
}

#[test]
fn zero_and_oversized_stoichiometry_rejected() {
    for number in [0, 2147483648_i64] {
        mutated_model(
            |m| m["reactions"][0]["products"]["A"] = json!(number),
            "CRN-E013",
            "/reactions/0/products/A",
        );
    }
}

#[test]
fn duplicate_reaction_names_rejected() {
    mutated_model(
        |m| m["reactions"][1]["name"] = json!("birth"),
        "CRN-E011",
        "/reactions/1/name",
    );
}

#[test]
fn population_limit_is_checked_separately() {
    mutated_model(
        |m| m["species"][0]["population_limit"] = json!(-1),
        "CRN-E016",
        "/species/0/population_limit",
    );
}

/// A limit below the initial count is a validation error; a valid limit is
/// refused by the converter (`CRN-R016`), which has no population limits.
#[test]
fn initial_state_must_lie_within_population_limits() {
    for (limit, invalid) in [(9, true), (10, false)] {
        let mut value = model_value();
        value["species"][0]["initial_count"] = json!(10);
        value["species"][0]["population_limit"] = json!(limit);
        let diagnostics = refused(&serde_json::to_vec(&value).unwrap());
        assert_eq!(
            has(&diagnostics, "CRN-E016", "/species/0/population_limit"),
            invalid,
            "limit {limit}"
        );
        assert_eq!(
            has(&diagnostics, "CRN-R016", "/species/0/population_limit"),
            !invalid,
            "limit {limit}: {diagnostics:?}"
        );
    }
}

#[test]
fn finite_and_positive_domains_checked() {
    for value in [json!(0.0), json!(-1.0)] {
        mutated_model(|m| m["omega"] = value, "CRN-E015", "/omega");
    }
    // A literal outside the f64 range is a parse error (CRN-E015).
    let text =
        String::from_utf8(bytes())
            .unwrap()
            .replacen("\"birth\": 4.0", "\"birth\": 1e400", 1);
    assert_eq!(parse_error(text.as_bytes()).code, "CRN-E015");
}

#[test]
fn all_rate_law_domains_and_species_reads_checked() {
    for kind in ["hill_activation", "hill_repression"] {
        let rate = json!({"kind":kind,"regulator":"A","vmax":1,"basal":0,"k":2,"n":2});
        let mut value = model_value();
        value["reactions"][0]["rate"] = rate;
        assert!(from_network(&serde_json::to_vec(&value).unwrap()).is_ok());
        value["reactions"][0]["rate"]["k"] = json!(0);
        let diagnostics = refused(&serde_json::to_vec(&value).unwrap());
        assert!(has(&diagnostics, "CRN-E015", "/reactions/0/rate/k"));
    }
    mutated_model(
        |m| {
            m["reactions"][0]["rate"] =
                json!({"kind":"michaelis_menten","substrate":"ghost","vmax":1,"km":1})
        },
        "CRN-E012",
        "/reactions/0/rate/substrate",
    );
    mutated_model(
        |m| {
            m["reactions"][0]["rate"] =
                json!({"kind":"michaelis_menten","substrate":"A","vmax":1,"km":0})
        },
        "CRN-E015",
        "/reactions/0/rate/km",
    );
    mutated_model(
        |m| m["reactions"][0]["rate"]["coefficient"] = json!(-1),
        "CRN-E015",
        "/reactions/0/rate/coefficient",
    );
}

#[test]
fn check_report_is_deterministic() {
    let one = serde_json::to_vec(&check_any(&bytes(), None)).unwrap();
    let two = serde_json::to_vec(&check_any(&bytes(), None)).unwrap();
    assert_eq!(one, two);
}

#[test]
fn full_repressilator_has_three_modules() {
    let model: Value = serde_json::from_slice(
        &std::fs::read(root().join("models/legacy/repressilator.model.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(model["species"].as_array().unwrap().len(), 6);
    let reactions = model["reactions"].as_array().unwrap();
    assert_eq!(reactions.len(), 12);
    let regulators: Vec<_> = reactions
        .iter()
        .filter(|r| r["rate"]["kind"] == "hill_repression")
        .map(|r| r["rate"]["regulator"].as_str().unwrap())
        .collect();
    assert_eq!(regulators, ["p3", "p1", "p2"]);
    // The converted model keeps the modules: one rule per reaction.
    let converted: Value =
        serde_json::from_slice(&from_network(&serde_json::to_vec(&model).unwrap()).unwrap())
            .unwrap();
    assert_eq!(converted["rules"].as_array().unwrap().len(), 12);
}
