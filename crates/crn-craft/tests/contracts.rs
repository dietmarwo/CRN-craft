use crn_craft::{
    check_bytes, contracts::*, describe, parse_model, parse_run, schema, sha256, validate,
};
use serde_json::{Value, json};
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn bytes() -> Vec<u8> {
    std::fs::read(root().join("models/examples/immigration-death.model.json")).unwrap()
}
fn model_value() -> Value {
    serde_json::from_slice(&bytes()).unwrap()
}
fn model() -> Model {
    parse_model(&bytes()).unwrap()
}
fn run() -> Run {
    parse_run(&std::fs::read(root().join("models/examples/immigration-death.run.json")).unwrap())
        .unwrap()
}
fn has(diagnostics: &[crn_craft::diagnostic::Diagnostic], code: &str, path: &str) -> bool {
    diagnostics.iter().any(|x| x.code == code && x.path == path)
}
fn mutated_model(edit: impl FnOnce(&mut Value), code: &str, path: &str) {
    let mut value = model_value();
    edit(&mut value);
    let report = check_bytes(&serde_json::to_vec(&value).unwrap(), None);
    assert!(!report.valid, "{report:?}");
    assert!(has(&report.diagnostics, code, path), "{report:?}");
}
fn mutated_run(edit: impl FnOnce(&mut Run), code: &str, path: &str) {
    let mut run = run();
    edit(&mut run);
    let diagnostics = validate::validate_run(&model(), &run, &sha256(&bytes()));
    assert!(has(&diagnostics, code, path), "{diagnostics:?}");
}

#[test]
fn all_example_models_and_run_references_pass() {
    for entry in std::fs::read_dir(root().join("models/examples")).unwrap() {
        let path = entry.unwrap().path();
        if path.to_string_lossy().ends_with(".run.json") {
            let run = parse_run(&std::fs::read(&path).unwrap()).unwrap();
            let report =
                crn_craft::check_files(&path.parent().unwrap().join(run.model.path), Some(&path));
            assert!(report.valid, "{}: {report:?}", path.display());
            assert!(!report.execution_available);
        }
    }
}

#[test]
fn negative_fixture_manifest_pins_codes_and_paths() {
    let base = root().join("models/invalid");
    let entries: Vec<Value> =
        serde_json::from_slice(&std::fs::read(base.join("expected.json")).unwrap()).unwrap();
    assert!(entries.len() >= 8);
    for entry in entries {
        let model = std::fs::read(base.join(entry["model"].as_str().unwrap())).unwrap();
        let run = entry["run"]
            .as_str()
            .map(|p| std::fs::read(base.join(p)).unwrap());
        let report = check_bytes(&model, run.as_deref());
        assert!(!report.valid, "{entry}: {report:?}");
        assert!(
            has(
                &report.diagnostics,
                entry["code"].as_str().unwrap(),
                entry["path"].as_str().unwrap()
            ),
            "{entry}: {report:?}"
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

#[test]
fn formatting_changes_model_identity() {
    let compact = serde_json::to_vec(&model_value()).unwrap();
    let report = check_bytes(
        &compact,
        Some(&std::fs::read(root().join("models/examples/immigration-death.run.json")).unwrap()),
    );
    assert!(has(&report.diagnostics, "CRN-E021", "/model/sha256"));
}

#[test]
fn roundtrip_preserves_integer_counts_above_double_precision() {
    let mut value = model_value();
    value["species"][0]["initial_count"] = json!(9_007_199_254_740_993_i64);
    let parsed = parse_model(&serde_json::to_vec(&value).unwrap()).unwrap();
    assert_eq!(parsed.species[0].initial_count, 9_007_199_254_740_993);
    let again = parse_model(&serde_json::to_vec(&parsed).unwrap()).unwrap();
    assert_eq!(
        again.species[0].initial_count,
        parsed.species[0].initial_count
    );
}

/// The shortest round-trip decimal of a double must parse back to the same
/// double, as Rust's `str::parse` and Python's `json` do. serde_json's default
/// best-effort parser misreads about a quarter of such literals by one ulp.
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
    let parsed = parse_model(text.as_bytes()).unwrap();
    let wrong: Vec<String> = values
        .iter()
        .enumerate()
        .filter(|(i, value)| parsed.parameters[&format!("p{i}")].to_bits() != value.to_bits())
        .map(|(i, value)| format!("{value:e} -> {:e}", parsed.parameters[&format!("p{i}")]))
        .collect();
    assert!(
        wrong.is_empty(),
        "{} of {} literals misparsed, e.g. {:?}",
        wrong.len(),
        values.len(),
        &wrong[..wrong.len().min(3)]
    );
    let again = parse_model(&serde_json::to_vec(&parsed).unwrap()).unwrap();
    assert!(
        parsed
            .parameters
            .iter()
            .all(|(name, value)| again.parameters[name].to_bits() == value.to_bits())
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
        let report = check_bytes(&serde_json::to_vec(&value).unwrap(), None);
        assert!(!report.valid, "{path}: {report:?}");
        assert!(report.diagnostics.iter().any(|d| d.code == "CRN-E003"));
    }
}

#[test]
fn unknown_fields_rejected_in_every_run_object() {
    let original = serde_json::to_value(run()).unwrap();
    for path in [
        "",
        "/model",
        "/method",
        "/time",
        "/ensemble",
        "/outputs",
        "/budget",
    ] {
        let mut value = original.clone();
        value
            .pointer_mut(path)
            .unwrap()
            .as_object_mut()
            .unwrap()
            .insert("typo".into(), json!(1));
        assert!(
            parse_run(&serde_json::to_vec(&value).unwrap()).is_err(),
            "{path}"
        );
    }
}

#[test]
fn unsupported_run_versions_methods_and_aggregators_rejected() {
    type Edit = fn(&mut Value);
    let original = serde_json::to_value(run()).unwrap();
    let cases: [(Edit, &str); 4] = [
        (|r| r["schema"] = json!("crn-craft.run/0.2"), "/schema"),
        (|r| r["method"] = json!({"kind": "ode"}), "/method/kind"),
        (
            |r| r["method"] = json!({"kind": "tau_leap", "epsilon": 0.03}),
            "/method/kind",
        ),
        // Tagged variants report an invalid field value at the variant object.
        (|r| r["method"]["aggregator"] = json!("nrm"), "/method"),
    ];
    for (edit, path) in cases {
        let mut value = original.clone();
        edit(&mut value);
        let error = parse_run(&serde_json::to_vec(&value).unwrap()).unwrap_err();
        assert_eq!(
            (error.code.as_str(), error.path.as_str()),
            ("CRN-E003", path)
        );
    }
}

#[test]
fn duplicate_keys_at_all_depths_and_escaped_spellings_rejected() {
    for (text, path) in [
        (r#"{"schema":0,"schema":1}"#, "/schema"),
        (r#"{"parameters":{"a":1,"\u0061":2}}"#, "/parameters/a"),
        (
            r#"{"reactions":[{"products":{"A":1,"A":2}}]}"#,
            "/reactions/0/products/A",
        ),
        (r#"{"a/b~c":1,"a/b~c":2}"#, "/a~1b~0c"),
    ] {
        let error = parse_model(text.as_bytes()).unwrap_err();
        assert_eq!(error.code, "CRN-E002");
        assert_eq!(error.path, path);
    }
}

#[test]
fn excessive_depth_is_rejected_without_stack_overflow() {
    let raw = format!("{}0{}", "[".repeat(1000), "]".repeat(1000));
    assert_eq!(parse_model(raw.as_bytes()).unwrap_err().code, "CRN-E004");
}

#[test]
fn nesting_limit_boundary_matches_its_documentation() {
    // docs/contracts.md: nesting depth at most 48, root at depth 0.
    let nested =
        |levels: usize, inner: &str| format!("{}{inner}{}", "[".repeat(levels), "]".repeat(levels));
    // Within the limit, the document fails only as a model (an array).
    for raw in [nested(48, "0"), nested(49, "")] {
        assert_eq!(parse_model(raw.as_bytes()).unwrap_err().code, "CRN-E003");
    }
    let error = parse_model(nested(49, "0").as_bytes()).unwrap_err();
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
    assert_eq!(parse_model(&bytes).unwrap_err().code, "CRN-E004");
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
        assert_eq!(parse_model(raw).unwrap_err().code, "CRN-E001");
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

#[test]
fn initial_state_must_lie_within_population_limits_and_fsp_bounds() {
    let mut model = model();
    model.species[0].initial_count = 10;
    for (limit, rejected) in [(9, true), (10, false)] {
        model.species[0].population_limit = Some(limit);
        assert_eq!(
            has(
                &validate::validate_model(&model),
                "CRN-E016",
                "/species/0/population_limit"
            ),
            rejected,
            "limit {limit}"
        );
    }
    let bytes =
        std::fs::read(root().join("models/examples/immigration-death.fsp.run.json")).unwrap();
    let mut run = parse_run(&bytes).unwrap();
    for (bound, rejected) in [(9, true), (10, false)] {
        if let Method::Fsp { bounds, .. } = &mut run.method {
            bounds.insert("A".into(), bound);
        }
        let diagnostics = validate::validate_run(&model, &run, &run.model.sha256);
        assert_eq!(
            has(&diagnostics, "CRN-E023", "/method/bounds/A"),
            rejected,
            "bound {bound}"
        );
    }
}

#[test]
fn output_budget_minimum_matches_the_documented_formula() {
    // docs/contracts.md: 8 * save_times * (1 + species * (statistics + kept_trajectories)).
    let run = run();
    let species = model().species.len() as u64;
    let minimum = 8
        * run.time.save_times.len() as u64
        * (1 + species * (run.outputs.statistics.len() as u64 + run.outputs.keep_trajectories));
    assert!(run.outputs.keep_trajectories > 0);
    for (budget, rejected) in [(minimum - 1, true), (minimum, false)] {
        mutated_run_status(
            |r| r.budget.max_output_bytes = budget,
            "CRN-E020",
            "/budget/max_output_bytes",
            rejected,
        );
    }
}

fn mutated_run_status(edit: impl FnOnce(&mut Run), code: &str, path: &str, expected: bool) {
    let mut run = run();
    edit(&mut run);
    let diagnostics = validate::validate_run(&model(), &run, &sha256(&bytes()));
    assert_eq!(has(&diagnostics, code, path), expected, "{diagnostics:?}");
}

#[test]
fn finite_and_positive_domains_checked_for_typed_callers() {
    for value in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let mut model = model();
        model.omega = value;
        assert!(has(&validate::validate_model(&model), "CRN-E015", "/omega"));
    }
    let mut model = model();
    model.parameters.insert("birth".into(), f64::NAN);
    assert!(has(
        &validate::validate_model(&model),
        "CRN-E015",
        "/parameters/birth"
    ));
}

#[test]
fn all_rate_law_domains_and_species_reads_checked() {
    for kind in ["hill_activation", "hill_repression"] {
        let rate = json!({"kind":kind,"regulator":"A","vmax":1,"basal":0,"k":2,"n":2});
        let mut value = model_value();
        value["reactions"][0]["rate"] = rate;
        assert!(check_bytes(&serde_json::to_vec(&value).unwrap(), None).valid);
        value["reactions"][0]["rate"]["k"] = json!(0);
        let report = check_bytes(&serde_json::to_vec(&value).unwrap(), None);
        assert!(has(&report.diagnostics, "CRN-E015", "/reactions/0/rate/k"));
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
fn warnings_do_not_make_valid_models_invalid_or_prove_growth() {
    let mut value = model_value();
    value["reactions"].as_array_mut().unwrap().pop();
    let report = check_bytes(&serde_json::to_vec(&value).unwrap(), None);
    assert!(report.valid);
    assert!(has(&report.diagnostics, "CRN-W001", "/species/0"));
    value["reactions"][0]["rate"]["coefficient"] = json!(0);
    let report = check_bytes(&serde_json::to_vec(&value).unwrap(), None);
    assert!(report.valid);
    assert!(has(
        &report.diagnostics,
        "CRN-W002",
        "/reactions/0/rate/coefficient"
    ));
}

#[test]
fn absorbing_model_with_no_reactions_is_valid() {
    let mut value = model_value();
    value["reactions"] = json!([]);
    assert!(check_bytes(&serde_json::to_vec(&value).unwrap(), None).valid);
}

#[test]
fn time_grid_errors_are_pinned() {
    mutated_run(|r| r.time.end = f64::INFINITY, "CRN-E020", "/time/end");
    mutated_run(
        |r| r.time.save_times[1] = f64::NAN,
        "CRN-E020",
        "/time/save_times/1",
    );
    mutated_run(
        |r| r.time.save_times.pop().map(|_| ()).unwrap(),
        "CRN-E020",
        "/time/save_times",
    );
    mutated_run(
        |r| r.time.save_times[0] = -1.0,
        "CRN-E020",
        "/time/save_times/0",
    );
}

#[test]
fn budget_and_replication_errors_are_pinned() {
    mutated_run(|r| r.budget.threads = 0, "CRN-E020", "/budget/threads");
    mutated_run(|r| r.budget.threads = 1025, "CRN-E020", "/budget/threads");
    mutated_run(
        |r| r.budget.wall_time_s = Some(f64::NAN),
        "CRN-E020",
        "/budget/wall_time_s",
    );
    mutated_run(
        |r| r.budget.max_events_per_replicate = None,
        "CRN-E020",
        "/budget/max_events_per_replicate",
    );
    mutated_run(
        |r| r.budget.max_output_bytes = 1,
        "CRN-E020",
        "/budget/max_output_bytes",
    );
    mutated_run(
        |r| r.ensemble.as_mut().unwrap().replicates = 1,
        "CRN-E020",
        "/ensemble/replicates",
    );
    mutated_run(|r| r.ensemble = None, "CRN-E020", "/ensemble");
    mutated_run(
        |r| r.outputs.keep_trajectories = u64::MAX,
        "CRN-E020",
        "/outputs/keep_trajectories",
    );
}

#[test]
fn output_selection_errors_are_pinned() {
    mutated_run(
        |r| r.outputs.statistics.push(Statistic::Mean),
        "CRN-E020",
        "/outputs/statistics/2",
    );
    mutated_run(
        |r| {
            r.outputs.statistics.clear();
            r.outputs.keep_trajectories = 0;
        },
        "CRN-E020",
        "/outputs",
    );
}

#[test]
fn invalid_relative_paths_and_hashes_rejected() {
    for path in [
        "", "/tmp/m", "../m", "a/../m", "./m", "a//m", "a/", "C:/m", "a\\m", "a\0m",
    ] {
        mutated_run(|r| r.model.path = path.into(), "CRN-E022", "/model/path");
    }
    mutated_run(
        |r| r.model.sha256 = "A".repeat(64),
        "CRN-E021",
        "/model/sha256",
    );
    mutated_run(
        |r| r.model.sha256 = "f".repeat(63),
        "CRN-E021",
        "/model/sha256",
    );
}

#[test]
fn fsp_allows_physically_unbounded_model_with_finite_truncation() {
    let bytes =
        std::fs::read(root().join("models/examples/immigration-death.fsp.run.json")).unwrap();
    let mut run = parse_run(&bytes).unwrap();
    let mut model = model();
    model.reactions.pop();
    assert!(validate::validate_run(&model, &run, &run.model.sha256).is_empty());
    if let Method::Fsp { bounds, .. } = &mut run.method {
        bounds.insert("ghost".into(), 5);
    }
    assert!(has(
        &validate::validate_run(&model, &run, &run.model.sha256),
        "CRN-E012",
        "/method/bounds/ghost"
    ));
}

#[test]
fn fsp_options_initial_state_and_method_incompatibility() {
    let bytes =
        std::fs::read(root().join("models/examples/immigration-death.fsp.run.json")).unwrap();
    let mut run = parse_run(&bytes).unwrap();
    if let Method::Fsp {
        bounds,
        max_states,
        tolerance,
    } = &mut run.method
    {
        bounds.insert("A".into(), -1);
        *max_states = 0;
        *tolerance = f64::NAN;
    }
    run.ensemble = Some(Ensemble {
        replicates: 2,
        root_seed: 0,
    });
    let diagnostics = validate::validate_run(&model(), &run, &run.model.sha256);
    assert!(has(&diagnostics, "CRN-E023", "/method/bounds/A"));
    assert!(has(&diagnostics, "CRN-E023", "/method"));
    assert!(has(&diagnostics, "CRN-E020", "/method"));
}

#[test]
fn description_distinguishes_exact_implemented_and_experimental_methods() {
    let value = describe();
    assert_eq!(value["simulation_available"], true);
    assert_eq!(
        value["commands"],
        json!([
            "describe",
            "check",
            "analyze",
            "export",
            "simulate",
            "verify",
            "simulate-cached",
            "fsp",
            "tau",
            "ssa"
        ])
    );
    assert_eq!(value["aggregators"][0]["status"], "implemented");
    assert_eq!(value["methods"][1]["status"], "implemented");
    assert_eq!(value["methods"][2]["status"], "experimental_approximation");
    for kind in [
        "model",
        "run",
        "result",
        "check",
        "fsp-result",
        "tau-request",
        "tau-result",
        "ssa-run",
        "ssa-result",
    ] {
        assert!(schema(kind).unwrap().is_object());
    }
    assert!(schema("unknown").is_none());
}

#[test]
fn check_report_is_deterministic() {
    let one = serde_json::to_vec(&check_bytes(&bytes(), None)).unwrap();
    let two = serde_json::to_vec(&check_bytes(&bytes(), None)).unwrap();
    assert_eq!(one, two);
}

#[test]
fn full_repressilator_has_three_modules() {
    let model = parse_model(
        &std::fs::read(root().join("models/examples/repressilator.model.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(model.species.len(), 6);
    assert_eq!(model.reactions.len(), 12);
    let regulators: Vec<_> = model
        .reactions
        .iter()
        .filter_map(|r| match &r.rate {
            Rate::HillRepression { regulator, .. } => Some(regulator.as_str()),
            _ => None,
        })
        .collect();
    assert_eq!(regulators, ["p3", "p1", "p2"]);
}
