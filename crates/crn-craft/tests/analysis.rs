//! Network analysis on the flat view (M5b): `crn-craft.analysis/1.0` of
//! network models (through the converter), rules models and
//! `crn-craft.model/1.0`, checked against the network engine's structure,
//! with labels of attributed models, warnings, states and refusals.
use crn_craft::rules::{FlatView, convert};
use crn_craft::{CompiledNetwork, analyze_bytes, analyze_files, parse_model, schema, sha256};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The examples, the workloads W1 to W7 and the two higher-order models.
fn network_models() -> Vec<PathBuf> {
    let mut out = Vec::new();
    for folder in [
        "models/examples",
        "validation/ssa-acceleration/models",
        "validation/merge/higher-order",
    ] {
        for entry in std::fs::read_dir(root().join(folder)).unwrap() {
            let path = entry.unwrap().path();
            if path.to_string_lossy().ends_with(".model.json") {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

fn name(path: &Path) -> String {
    path.file_name().unwrap().to_string_lossy().into_owned()
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

/// The propensity tolerance M2 and M4b declared for converted rate laws
/// against the network engine (see tests/flat_view.rs).
fn allowed(reaction: &Value, parameters: &Value, old: f64) -> f64 {
    let scalar = |v: &Value| {
        v.as_f64()
            .unwrap_or_else(|| parameters[v.as_str().unwrap()].as_f64().unwrap())
    };
    let rate = &reaction["rate"];
    let order: i64 = reaction["reactants"]
        .as_object()
        .unwrap()
        .values()
        .map(|v| v.as_i64().unwrap())
        .sum();
    match rate["kind"].as_str().unwrap() {
        "mass_action" if order == 2 => 4.0 * f64::EPSILON * old.abs(),
        "mass_action" | "michaelis_menten" => 0.0,
        "hill_activation" | "hill_repression" => {
            let n = scalar(&rate["n"]);
            (1e-14 * n.max(1.0) * old.abs()).max(scalar(&rate["vmax"]) * 2f64.powi(-1020))
        }
        other => panic!("unknown rate law {other}"),
    }
}

fn codes(report: &crn_craft::AnalysisReport) -> Vec<(String, String)> {
    report
        .diagnostics
        .iter()
        .map(|d| (d.code.clone(), d.path.clone()))
        .collect()
}

fn refused(model: &[u8], state: Option<&[u8]>) -> (String, String) {
    let report = analyze_bytes(model, state);
    assert!(!report.success);
    assert!(report.analysis.is_none());
    assert_eq!(report.model_sha256, Some(sha256(model)));
    let first = &report.diagnostics[0];
    (first.code.clone(), first.path.clone())
}

#[test]
fn network_models_keep_the_network_structure_through_the_converter() {
    let mut models = 0;
    for path in network_models() {
        let label = name(&path);
        let network_bytes = std::fs::read(&path).unwrap();
        let report = analyze_bytes(&network_bytes, None);
        assert!(report.success, "{label}: {:?}", report.diagnostics);
        let value = serde_json::to_value(&report).unwrap();
        assert_eq!(value["schema"], "crn-craft.analysis/1.0", "{label}");
        assert_eq!(report.model_sha256, Some(sha256(&network_bytes)));
        assert_eq!(report.state_sha256, None);
        // The converted model gives the same report, apart from its hash.
        let converted = convert::from_network(&network_bytes).unwrap();
        let mut other = serde_json::to_value(analyze_bytes(&converted, None)).unwrap();
        assert_eq!(other["model_sha256"], json!(sha256(&converted)));
        other["model_sha256"] = value["model_sha256"].clone();
        assert_eq!(other, value, "{label}");
        // Species, reactions, reads and dependencies of the network engine.
        let network = CompiledNetwork::compile(&parse_model(&network_bytes).unwrap()).unwrap();
        let analysis = report.analysis.unwrap();
        assert_eq!(analysis.species, network.species(), "{label}");
        assert_eq!(analysis.initial_counts, network.initial_counts(), "{label}");
        assert_eq!(analysis.evaluation_counts, network.initial_counts());
        assert_eq!(
            serde_json::to_value(&analysis.reactions).unwrap(),
            serde_json::to_value(network.reactions()).unwrap(),
            "{label}"
        );
        assert_eq!(analysis.dependencies, network.dependencies(), "{label}");
        let source: Value = serde_json::from_slice(&network_bytes).unwrap();
        let old = network.propensities(network.initial_counts()).unwrap();
        for (r, (new, old)) in analysis.propensities.iter().zip(&old).enumerate() {
            let bound = allowed(&source["reactions"][r], &source["parameters"], *old);
            assert!(
                (new - old).abs() <= bound,
                "{label} reaction {r}: {new} against {old}"
            );
        }
        models += 1;
    }
    assert_eq!(models, 14);
}

/// Two copies of `M` that switch between `wt = 0` and `wt = 1`.
fn switching() -> Value {
    json!({
        "schema": "crn-craft.model/1.0", "name": "switching",
        "parameters": {"k": 0.5, "r": 0.25},
        "types": [{"name": "M", "attributes": ["wt"]}],
        "initial": [
            {"type": "M", "attributes": {"wt": 0}, "count": 10},
            {"type": "M", "attributes": {"wt": 1}, "count": 2}
        ],
        "observables": {"mutants": {"type": "M", "where": "x.wt == 1"}},
        "rules": [
            {"name": "mutate", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 0"}],
             "rate": "k", "products": [{"type": "M", "attributes": {"wt": 1}}]},
            {"name": "revert", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 1"}],
             "rate": "r", "products": [{"type": "M", "attributes": {"wt": 0}}]}
        ]
    })
}

#[test]
fn attributed_models_are_analyzed_on_their_flat_labels() {
    let model = bytes(&switching());
    let report = analyze_bytes(&model, None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert!(report.diagnostics.is_empty(), "{:?}", report.diagnostics);
    let analysis = report.analysis.unwrap();
    assert_eq!(analysis.species, ["M{wt=0}", "M{wt=1}"]);
    assert_eq!(analysis.initial_counts, [10, 2]);
    assert_eq!(analysis.propensities, [5.0, 0.5]);
    let reactions = serde_json::to_value(&analysis.reactions).unwrap();
    assert_eq!(
        reactions,
        json!([
            {"name": "mutate", "reactants": [{"species": 0, "coefficient": 1}],
             "products": [{"species": 1, "coefficient": 1}],
             "net_changes": [{"species": 0, "coefficient": -1}, {"species": 1, "coefficient": 1}],
             "reads": [0]},
            {"name": "revert", "reactants": [{"species": 1, "coefficient": 1}],
             "products": [{"species": 0, "coefficient": 1}],
             "net_changes": [{"species": 0, "coefficient": 1}, {"species": 1, "coefficient": -1}],
             "reads": [1]}
        ])
    );
    assert_eq!(analysis.dependencies, [vec![0, 1], vec![0, 1]]);
    assert_eq!(analysis.conservation.basis[0].coefficients, ["1", "1"]);
    assert_eq!(analysis.conservation.basis[0].initial_total, "12");
    assert_eq!(
        analysis.initial_class_upper_bounds,
        [Some("12".into()), Some("12".into())]
    );
    assert!(analysis.linkage.weakly_reversible);
    assert_eq!(analysis.linkage.deficiency, Some(0));
    // A state names the species by label.
    let state = br#"{"schema": "crn-craft.state/0.1", "counts": {"M{wt=0}": 4, "M{wt=1}": 8}}"#;
    let report = analyze_bytes(&model, Some(state));
    assert!(report.success, "{:?}", report.diagnostics);
    assert_eq!(report.state_sha256, Some(sha256(state)));
    let analysis = report.analysis.unwrap();
    assert_eq!(analysis.evaluation_counts, [4, 8]);
    assert_eq!(analysis.propensities, [2.0, 2.0]);
    // A rule over several tuples has one reaction per tuple, by rule name.
    let mut loss = switching();
    loss["rules"][1] = json!({"name": "loss", "reactants": [{"label": "x", "type": "M"}],
                              "rate": "r", "products": []});
    let analysis = analyze_bytes(&bytes(&loss), None).analysis.unwrap();
    let names: Vec<&str> = analysis.reactions.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["mutate", "loss", "loss"]);
    assert_eq!(analysis.propensities, [5.0, 2.5, 0.5]);
}

#[test]
fn higher_order_and_whole_rule_reactions_keep_their_coefficients() {
    let model = json!({
        "schema": "crn-craft.model/1.0", "name": "orders",
        "parameters": {"c": 0.25},
        "types": [{"name": "A"}, {"name": "B"}],
        "initial": [{"type": "A", "count": 7}, {"type": "B", "count": 2}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [
            {"name": "triple", "reactants": [{"label": "x", "type": "A", "count": 3}], "rate": "c", "products": []},
            {"name": "whole", "reactants": [{"label": "x", "type": "A", "count": 3}, {"label": "y", "type": "B"}],
             "propensity": "c * B", "products": [{"type": "B"}, {"type": "B"}]}
        ]
    });
    let analysis = analyze_bytes(&bytes(&model), None).analysis.unwrap();
    let reactions = serde_json::to_value(&analysis.reactions).unwrap();
    assert_eq!(
        reactions[0]["reactants"],
        json!([{"species": 0, "coefficient": 3}])
    );
    assert_eq!(
        reactions[1]["reactants"],
        json!([{"species": 0, "coefficient": 3}, {"species": 1, "coefficient": 1}])
    );
    assert_eq!(
        reactions[1]["products"],
        json!([{"species": 1, "coefficient": 2}])
    );
    assert_eq!(
        reactions[1]["net_changes"],
        json!([{"species": 0, "coefficient": -3}, {"species": 1, "coefficient": 1}])
    );
    assert_eq!(reactions[1]["reads"], json!([0, 1]));
    // c C(7, 3) and c B while A >= 3 and B >= 1.
    assert_eq!(analysis.propensities, [8.75, 0.5]);
}

#[test]
fn types_that_never_hold_molecules_keep_their_species_and_reactions() {
    // B starts at 0 and nothing makes it: the network listed it, and the
    // reaction that needs it, and so does the analysis.
    let network = json!({
        "schema": "crn-craft.model/0.1", "name": "unused", "convention": "stochastic_counts",
        "omega": 1.0, "parameters": {},
        "species": [{"name": "A", "initial_count": 1}, {"name": "B", "initial_count": 0}],
        "reactions": [
            {"name": "idle", "reactants": {}, "products": {}, "rate": {"kind": "mass_action", "coefficient": 1.0}},
            {"name": "convert", "reactants": {"B": 1}, "products": {"A": 1},
             "rate": {"kind": "mass_action", "coefficient": 0}}
        ]
    });
    let report = analyze_bytes(&bytes(&network), None);
    assert!(report.success, "{:?}", report.diagnostics);
    // W003 for the reaction without effect, W002 for the zero rate; A + B
    // is conserved, so no growth warning for A.
    assert_eq!(
        codes(&report),
        [
            ("CRN-W003".to_string(), "/reactions/0".to_string()),
            ("CRN-W002".to_string(), "/reactions/1".to_string())
        ]
    );
    let analysis = report.analysis.unwrap();
    assert_eq!(analysis.species, ["A", "B"]);
    assert_eq!(analysis.reactions.len(), 2);
    assert_eq!(analysis.propensities, [1.0, 0.0]);
    assert_eq!(
        analysis.initial_class_upper_bounds,
        [Some("1".into()), Some("1".into())]
    );
}

#[test]
fn growth_warnings_name_the_species_without_a_proven_bound() {
    let network = json!({
        "schema": "crn-craft.model/0.1", "name": "growth", "convention": "stochastic_counts",
        "omega": 1.0, "parameters": {},
        "species": [{"name": "A", "initial_count": 0}, {"name": "B", "initial_count": 3}],
        "reactions": [
            {"name": "make", "reactants": {}, "products": {"A": 1}, "rate": {"kind": "mass_action", "coefficient": 1.0}},
            {"name": "swap", "reactants": {"B": 1}, "products": {"B": 1}, "rate": {"kind": "mass_action", "coefficient": 1.0}}
        ]
    });
    let report = analyze_bytes(&bytes(&network), None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert_eq!(
        codes(&report),
        [
            ("CRN-W003".to_string(), "/reactions/1".to_string()),
            ("CRN-W001".to_string(), "/species/0".to_string())
        ]
    );
}

#[test]
fn models_without_an_analysis_are_refused_with_coded_diagnostics() {
    let network = |species: Value, reactions: Value| {
        bytes(&json!({
            "schema": "crn-craft.model/0.1", "name": "m", "convention": "stochastic_counts",
            "omega": 1.0, "parameters": {}, "species": species, "reactions": reactions
        }))
    };
    let decay = json!([{"name": "decay", "reactants": {"A": 1}, "products": {},
                        "rate": {"kind": "mass_action", "coefficient": 1.0}}]);
    // A rules model needs a rule.
    assert_eq!(
        refused(
            &network(json!([{"name": "A", "initial_count": 1}]), json!([])),
            None
        ),
        ("CRN-R016".into(), "/reactions".into())
    );
    assert_eq!(
        refused(
            &network(
                json!([{"name": "A", "initial_count": 1, "population_limit": 5}]),
                decay.clone()
            ),
            None
        ),
        ("CRN-R016".into(), "/species/0/population_limit".into())
    );
    // Network validation errors are kept.
    let unknown = json!([{"name": "decay", "reactants": {"A": 1}, "products": {"B": 1},
                          "rate": {"kind": "mass_action", "coefficient": 1.0}}]);
    assert_eq!(
        refused(
            &network(json!([{"name": "A", "initial_count": 1}]), unknown),
            None
        ),
        ("CRN-E012".into(), "/reactions/0/products/B".into())
    );
    let mut extra: Value =
        serde_json::from_slice(&network(json!([{"name": "A", "initial_count": 1}]), decay))
            .unwrap();
    extra["unexpected"] = json!(true);
    assert_eq!(refused(&bytes(&extra), None).0, "CRN-E003");
    assert_eq!(refused(b"{", None).0, "CRN-E001");
    // Models without a flat view.
    // `revert` acting on wt = 1 makes wt = 2, a tuple not known before the run.
    let mut growing = switching();
    growing["rules"][1]["products"][0] =
        json!({"type": "M", "copy": "x", "attributes": {"wt": "x.wt + 1"}});
    assert_eq!(
        refused(&bytes(&growing), None),
        ("CRN-R017".into(), "/rules".into())
    );
    let atp = std::fs::read(root().join("models/rules/mtdna-three-class-atp.rules.json")).unwrap();
    assert_eq!(refused(&atp, None), ("CRN-R017".into(), "/atp".into()));
    let mut invalid = switching();
    invalid["rules"][0]["rate"] = json!("k +");
    assert_eq!(refused(&bytes(&invalid), None).0, "CRN-R012");
}

#[test]
fn states_are_checked_by_label() {
    let model = bytes(&switching());
    for (state, path) in [
        (json!({"M{wt=0}": 1}), "/counts/M{wt=1}"),
        (json!({"M{wt=0}": 1, "M{wt=1}": 1, "M": 1}), "/counts/M"),
        (json!({"M{wt=0}": -1, "M{wt=1}": 1}), "/counts/M{wt=0}"),
        (
            json!({"M{wt=0}": 9_007_199_254_740_993_i64, "M{wt=1}": 1}),
            "/counts/M{wt=0}",
        ),
    ] {
        let document = bytes(&json!({"schema": "crn-craft.state/0.1", "counts": state}));
        assert_eq!(
            refused(&model, Some(&document)),
            ("CRN-E041".into(), path.into()),
            "{state}"
        );
    }
    let at_limit = bytes(&json!({"schema": "crn-craft.state/0.1",
                                 "counts": {"M{wt=0}": 9_007_199_254_740_992_i64, "M{wt=1}": 0}}));
    assert!(analyze_bytes(&model, Some(&at_limit)).success);
    for (text, code) in [
        (
            r#"{"schema":"crn-craft.state/0.1","counts":{"A":1,"A":2}}"#,
            "CRN-E002",
        ),
        (
            r#"{"schema":"crn-craft.state/0.1","counts":{},"extra":0}"#,
            "CRN-E003",
        ),
        (r#"{"schema":"crn-craft.state/99","counts":{}}"#, "CRN-E003"),
    ] {
        assert_eq!(refused(&model, Some(text.as_bytes())).0, code, "{text}");
    }
}

#[test]
fn propensity_failures_are_engine_failures_and_underflow_is_zero() {
    let one = |reactants: Value, rate: Value, count: i64| {
        bytes(&json!({
            "schema": "crn-craft.model/0.1", "name": "m", "convention": "stochastic_counts",
            "omega": 1.0, "parameters": {},
            "species": [{"name": "A", "initial_count": count}],
            "reactions": [{"name": "r", "reactants": reactants, "products": {}, "rate": rate}]
        }))
    };
    // C(2000, 1000) overflows f64: the Static Direct replicate failure.
    let overflow = one(
        json!({"A": 1000}),
        json!({"kind": "mass_action", "coefficient": 1.0}),
        2000,
    );
    assert_eq!(refused(&overflow, None).0, "CRN-R030");
    // 0.6.0 refused a positive propensity that underflows (CRN-E042); the
    // rules engine evaluates it to 0.
    let underflow = one(
        json!({"A": 1}),
        json!({"kind": "michaelis_menten", "substrate": "A", "vmax": 5e-324, "km": 1e308}),
        1,
    );
    let report = analyze_bytes(&underflow, None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert_eq!(report.analysis.unwrap().propensities, [0.0]);
}

#[test]
fn files_and_the_generated_schema() {
    let temp = std::env::temp_dir().join(format!("crn-analysis-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    let model = temp.join("model.json");
    let state = temp.join("state.json");
    std::fs::write(&model, bytes(&switching())).unwrap();
    std::fs::write(
        &state,
        br#"{"schema": "crn-craft.state/0.1", "counts": {"M{wt=0}": 0, "M{wt=1}": 4}}"#,
    )
    .unwrap();
    let report = analyze_files(&model, Some(&state));
    assert!(report.success, "{:?}", report.diagnostics);
    assert_eq!(report.analysis.unwrap().propensities, [0.0, 1.0]);
    let missing = analyze_files(&temp.join("missing.json"), None);
    assert_eq!(missing.diagnostics[0].code, "CRN-E030");
    std::fs::remove_dir_all(&temp).unwrap();
    // The 1.0 report has a generated schema; the 0.2 and 0.1 reports are no
    // longer generated (their files stay until M8 removes them).
    let generated = schema("analysis-1.0").unwrap();
    let snapshot: Value = serde_json::from_slice(
        &std::fs::read(root().join("schemas/analysis-1.0.schema.json")).unwrap(),
    )
    .unwrap();
    assert!(
        generated == snapshot,
        "regenerate schemas/analysis-1.0.schema.json"
    );
    assert_eq!(
        generated["definitions"]["AnalysisVersion"]["enum"],
        json!(["crn-craft.analysis/1.0"])
    );
    assert!(schema("analysis").is_none() && schema("analysis-v1").is_none());
    let description = crn_craft::describe_network();
    assert_eq!(description["schema"], "crn-craft.describe/0.11");
    assert_eq!(description["schemas"]["analysis-1.0"], generated);
    assert!(
        !description["commands"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c == "analyze-v1")
    );
}

#[test]
fn the_analysis_view_of_a_converted_model_equals_the_converter_output() {
    for path in network_models() {
        let network_bytes = std::fs::read(&path).unwrap();
        let direct = FlatView::of_any_model(&network_bytes).unwrap();
        let converted = convert::from_network(&network_bytes).unwrap();
        let model = serde_json::from_slice(&converted).unwrap();
        let view = FlatView::with_every_type(&model).unwrap();
        assert_eq!(direct.species(), view.species());
        assert_eq!(direct.channels(), view.channels());
    }
}
