//! Quantile statistics of `crn-craft.run/1.0` (`outputs.quantiles`) and
//! `crn-craft.result/1.0` (`statistics.quantiles`).
//!
//! The quantile of probability p over the n completed replicates is NumPy's
//! default `linear` quantile (Hyndman and Fan type 7): with the values
//! sorted, h = (n - 1) p, j = floor(h), g = h - j, the value is
//! x[j] + g (x[j+1] - x[j]) for g < 0.5 and x[j+1] - (1 - g) (x[j+1] - x[j])
//! otherwise, and x[n-1] when j = n - 1. These tests apply that definition to
//! the kept paths and require equal bits.
use crn_craft::rules::{Control, RulesResult, execute, prepare};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn immigration_death(schema: &str) -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": schema,
        "name": "immigration_death",
        "parameters": {"k_in": 10.0, "k_out": 1.0},
        "types": [{"name": "X"}, {"name": "Y"}],
        "initial": [{"type": "Y", "count": 3}],
        "observables": {"X": {"type": "X"}, "Y": {"type": "Y"}},
        "rules": [
            {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out",
             "products": []},
            {"name": "y_birth", "reactants": [{"label": "y", "type": "Y"}], "rate": "0.3",
             "products": [{"type": "Y"}, {"type": "Y"}]},
            {"name": "y_death", "reactants": [{"label": "y", "type": "Y"}], "rate": "0.4",
             "products": []}
        ]
    }))
    .unwrap()
}

/// Unit births and deaths of one type, the scope of `binomial`.
fn birth_death() -> Vec<u8> {
    serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "birth_death",
        "parameters": {"b": 0.9, "d": 1.0},
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 30}],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {"name": "birth", "reactants": [{"label": "x", "type": "X"}], "rate": "b",
             "products": [{"type": "X"}, {"type": "X"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "d",
             "products": []}
        ]
    }))
    .unwrap()
}

fn request(model: &[u8], quantiles: Option<Value>) -> Value {
    let mut request = json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 4.0, "report_times": [0.0, 0.5, 1.0, 2.0, 4.0]},
        "ensemble": {"replicates": 200, "root_seed": 17},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 200, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1_000_000,
                   "max_sampling_steps_per_replicate": 10_000_000,
                   "max_instances": 16, "threads": 2, "max_output_bytes": 1 << 26},
        "method": {"kind": "static_direct"}
    });
    if let Some(quantiles) = quantiles {
        request["outputs"]["quantiles"] = quantiles;
    }
    request
}

fn try_run(
    model: &[u8],
    request: &Value,
) -> Result<RulesResult, Vec<crn_craft::diagnostic::Diagnostic>> {
    let prepared = prepare(model, &serde_json::to_vec(request).unwrap())?;
    Ok(execute(&prepared, &Control::default()).unwrap())
}

fn run(model: &[u8], request: &Value) -> Value {
    serde_json::to_value(try_run(model, request).unwrap_or_else(|d| panic!("{d:?}"))).unwrap()
}

fn refusal(model: &[u8], request: &Value) -> (String, String, String) {
    let diagnostics = match prepare(model, &serde_json::to_vec(request).unwrap()) {
        Ok(_) => panic!("accepted: {request}"),
        Err(d) => d,
    };
    let first = &diagnostics[0];
    (first.code.clone(), first.path.clone(), first.hint.clone())
}

/// The documented definition, applied to sorted values.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn linear(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    let h = (n - 1) as f64 * p;
    let j = h.floor() as usize;
    if j + 1 >= n {
        return sorted[n - 1];
    }
    let g = h - j as f64;
    let (a, b) = (sorted[j], sorted[j + 1]);
    if g < 0.5 {
        a + g * (b - a)
    } else {
        b - (1.0 - g) * (b - a)
    }
}

/// The quantiles of the result's kept paths, report time x observable, for
/// each probability. Every replicate must be kept.
fn from_paths(result: &Value, probabilities: &[f64]) -> Vec<Vec<Vec<f64>>> {
    let paths = result["paths"].as_array().unwrap();
    assert_eq!(
        paths.len() as u64,
        result["requested_replicates"].as_u64().unwrap()
    );
    let rows = result["report_times"].as_array().unwrap().len();
    let columns = result["observables"].as_array().unwrap().len();
    let mut out = vec![vec![vec![0.0; columns]; rows]; probabilities.len()];
    for row in 0..rows {
        for column in 0..columns {
            let mut cell: Vec<f64> = paths
                .iter()
                .map(|p| p["values"][row][column].as_f64().unwrap())
                .collect();
            cell.sort_by(f64::total_cmp);
            for (k, &p) in probabilities.iter().enumerate() {
                out[k][row][column] = linear(&cell, p);
            }
        }
    }
    out
}

/// `statistics.quantiles` as (p, values).
fn quantiles_of(result: &Value) -> Vec<(f64, Vec<Vec<f64>>)> {
    result["statistics"]["quantiles"]
        .as_array()
        .unwrap_or_else(|| panic!("no quantiles in {}", result["statistics"]))
        .iter()
        .map(|q| {
            assert_eq!(q.as_object().unwrap().len(), 2, "{q}");
            (
                q["p"].as_f64().unwrap(),
                serde_json::from_value(q["values"].clone()).unwrap(),
            )
        })
        .collect()
}

fn assert_bits(left: &[Vec<f64>], right: &[Vec<f64>], what: &str) {
    assert_eq!(left.len(), right.len(), "{what}");
    for (a, b) in left.iter().zip(right) {
        assert_eq!(a.len(), b.len(), "{what}");
        for (x, y) in a.iter().zip(b) {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: {x} != {y}");
        }
    }
}

fn check_against_paths(result: &Value, probabilities: &[f64], what: &str) {
    let reported = quantiles_of(result);
    let expected = from_paths(result, probabilities);
    assert_eq!(
        reported.iter().map(|q| q.0).collect::<Vec<_>>(),
        probabilities,
        "{what}: probabilities in request order"
    );
    for ((_, values), expected) in reported.iter().zip(&expected) {
        assert_bits(values, expected, what);
    }
}

/// The result without `run_sha256` and the quantiles, as JSON text.
fn without_quantiles(result: &Value) -> String {
    let mut value = result.clone();
    value.as_object_mut().unwrap().remove("run_sha256");
    if let Some(stats) = value["statistics"].as_object_mut() {
        stats.remove("quantiles");
    }
    value.to_string()
}

const PROBABILITIES: [f64; 12] = [
    0.001,
    0.01,
    0.025,
    0.1,
    0.25,
    1.0 / 3.0,
    0.5,
    0.6,
    0.75,
    0.9,
    0.975,
    0.999,
];

#[test]
fn quantiles_equal_the_linear_definition_on_the_kept_paths_bit_for_bit() {
    let model = immigration_death("crn-craft.model/1.0");
    let probabilities = PROBABILITIES.to_vec();
    for method in [
        json!({"kind": "static_direct"}),
        json!({"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": false}),
        json!({"kind": "tree"}),
        json!({"kind": "next_reaction"}),
        json!({"kind": "tau", "epsilon": 0.03}),
    ] {
        let mut with = request(&model, Some(json!(probabilities)));
        with["method"] = method.clone();
        let result = run(&model, &with);
        assert_eq!(result["schema"], "crn-craft.result/1.0");
        check_against_paths(&result, &probabilities, &method.to_string());
        // Quantiles change nothing else.
        let mut plain = with.clone();
        plain["outputs"]
            .as_object_mut()
            .unwrap()
            .remove("quantiles");
        assert_eq!(
            without_quantiles(&result),
            without_quantiles(&run(&model, &plain)),
            "{method}"
        );
    }
    // binomial, on its unit birth-death scope.
    let model = birth_death();
    let mut with = request(&model, Some(json!(probabilities)));
    with["method"] = json!({"kind": "binomial", "step": 0.01});
    check_against_paths(&run(&model, &with), &probabilities, "binomial");
}

#[test]
fn a_single_replicate_and_a_stop_give_the_held_values() {
    let model = immigration_death("crn-craft.model/1.0");
    let mut one = request(&model, Some(json!([0.025, 0.5, 0.975])));
    one["ensemble"]["replicates"] = 1.into();
    one["outputs"]["keep_paths"] = 1.into();
    let result = run(&model, &one);
    let path: Vec<Vec<f64>> = serde_json::from_value(result["paths"][0]["values"].clone()).unwrap();
    for (_, values) in quantiles_of(&result) {
        assert_bits(&values, &path, "one replicate");
    }
    // Stopped replicates hold their state and count as complete.
    let mut stopped: Value = serde_json::from_slice(&model).unwrap();
    stopped["stop"] = "X >= 8".into();
    let stopped = serde_json::to_vec(&stopped).unwrap();
    let result = run(&stopped, &request(&stopped, Some(json!([0.1, 0.5, 0.9]))));
    assert!(result["stopped_replicates"].as_u64().unwrap() > 0);
    check_against_paths(&result, &[0.1, 0.5, 0.9], "stop");
}

#[test]
fn every_stored_request_as_run_1_0_with_quantiles_matches_its_paths() {
    // Every method the stored requests use (tree, static_direct with each
    // option, next_reaction, tau, mtdna_direct, hybrid), on their models.
    let stored = root().join("crates/crn-craft-rules/tests/stored");
    let expected: BTreeMap<String, String> =
        serde_json::from_slice(&std::fs::read(stored.join("expected.json")).unwrap()).unwrap();
    let probabilities = [0.025, 0.25, 0.5, 0.75, 0.975];
    let mut checked = 0;
    for name in expected.keys() {
        let mut request: Value =
            serde_json::from_slice(&std::fs::read(stored.join(name)).unwrap()).unwrap();
        let model = std::fs::read(
            root()
                .join("models/rules")
                .join(request["model"]["path"].as_str().unwrap()),
        )
        .unwrap();
        request["schema"] = "crn-craft.run/1.0".into();
        request["outputs"]["keep_paths"] = request["ensemble"]["replicates"].clone();
        let plain = run(&model, &request);
        request["outputs"]["quantiles"] = json!(probabilities);
        let result = run(&model, &request);
        assert_eq!(
            without_quantiles(&result),
            without_quantiles(&plain),
            "{name}"
        );
        if result["statistics"].is_null() {
            assert!(plain["statistics"].is_null(), "{name}");
            continue;
        }
        check_against_paths(&result, &probabilities, name);
        checked += 1;
    }
    assert!(checked >= 15, "{checked} stored requests with statistics");
}

#[test]
fn thinning_with_atp_reports_quantiles_of_every_column() {
    let model =
        std::fs::read(root().join("models/rules/mtdna-three-class-atp.rules.json")).unwrap();
    let mut with = request(&model, Some(json!([0.05, 0.5, 0.95])));
    with["method"] = json!({"kind": "thinning"});
    with["time"] = json!({"end": 200.0, "report_times": [0.0, 100.0, 200.0]});
    with["ensemble"]["replicates"] = 12.into();
    with["outputs"]["keep_paths"] = 12.into();
    with["budget"]["max_instances"] = 64.into();
    let result = run(&model, &with);
    assert!(
        result["observables"]
            .as_array()
            .unwrap()
            .contains(&"ATP".into())
    );
    check_against_paths(&result, &[0.05, 0.5, 0.95], "thinning");
}

#[test]
fn quantiles_do_not_depend_on_the_thread_count() {
    let model = immigration_death("crn-craft.model/1.0");
    let mut reference = None;
    for threads in [1, 2, 4] {
        let mut with = request(&model, Some(json!(PROBABILITIES)));
        with["ensemble"]["replicates"] = 500.into();
        with["outputs"]["keep_paths"] = 0.into();
        with["budget"]["threads"] = threads.into();
        let mut result = run(&model, &with);
        assert!(result["statistics"]["quantiles"].is_array());
        result.as_object_mut().unwrap().remove("run_sha256");
        let text = result.to_string();
        match &reference {
            None => reference = Some(text),
            Some(first) => assert!(*first == text, "{threads} threads differ"),
        }
    }
}

#[test]
fn quantiles_without_other_statistics_and_without_quantiles() {
    let model = immigration_death("crn-craft.model/1.0");
    // Only quantiles: statistics are reported, mean and variance are null.
    let mut only = request(&model, Some(json!([0.5])));
    only["outputs"]["statistics"] = json!([]);
    let result = run(&model, &only);
    let stats = &result["statistics"];
    assert_eq!(stats["samples"], 200);
    assert!(stats["mean"].is_null() && stats["sample_variance"].is_null());
    check_against_paths(&result, &[0.5], "only quantiles");
    // No statistics and no quantiles: withheld as before.
    only["outputs"]["quantiles"] = json!([]);
    assert!(run(&model, &only)["statistics"].is_null());
    // An empty list is the absent field: no quantiles key.
    let empty = run(&model, &request(&model, Some(json!([]))));
    let absent = run(&model, &request(&model, None));
    assert!(empty["statistics"].get("quantiles").is_none());
    assert!(absent["statistics"].get("quantiles").is_none());
    assert_eq!(without_quantiles(&empty), without_quantiles(&absent));
    // A rules-model with a run/1.0 request also has them.
    let old = immigration_death("crn-craft.rules-model/0.1");
    check_against_paths(
        &run(&old, &request(&old, Some(json!([0.2, 0.8])))),
        &[0.2, 0.8],
        "rules-model/0.1",
    );
}

#[test]
fn quantiles_are_withheld_with_the_other_statistics() {
    let model = immigration_death("crn-craft.model/1.0");
    // An event limit: incomplete ensemble.
    let mut limited = request(&model, Some(json!([0.5])));
    limited["budget"]["max_events_per_replicate"] = 5.into();
    let result = run(&model, &limited);
    assert!(result["statistics"].is_null());
    assert_eq!(result["diagnostics"][0]["code"], "CRN-R050");
    // An output limit drops the statistics: a budget that holds the result
    // without quantiles but not with them.
    let plain = try_run(&model, &request(&model, Some(json!([])))).unwrap();
    let size = serde_json::to_vec(&plain).unwrap().len() + 1;
    let mut small = request(&model, Some(json!(PROBABILITIES)));
    small["budget"]["max_output_bytes"] = (size + 100).into();
    let result = run(&model, &small);
    assert_eq!(result["termination"], "output_limit");
    assert!(result["statistics"].is_null());
    small["outputs"]["quantiles"] = json!([]);
    assert_eq!(run(&model, &small)["termination"], "completed");
}

#[test]
fn bound_candidates_keep_the_quantiles() {
    let model = immigration_death("crn-craft.model/1.0");
    let request = serde_json::to_vec(&request(&model, Some(json!([0.025, 0.5, 0.975])))).unwrap();
    let template = prepare(&model, &request).unwrap();
    let candidate = template
        .bind_parameters(&BTreeMap::from([("k_in".to_string(), 25.0)]))
        .unwrap();
    let result = serde_json::to_value(execute(&candidate, &Control::default()).unwrap()).unwrap();
    assert_eq!(result["schema"], "crn-craft.result/1.0");
    check_against_paths(&result, &[0.025, 0.5, 0.975], "candidate");
    let median = &quantiles_of(&result)[1].1;
    // X at t = 4 has a Poisson law with mean 25 (1 - e^-4), about 24.5.
    assert!((median[4][0] - 24.5).abs() <= 3.0, "{median:?}");
}

#[test]
fn old_request_versions_refuse_the_field() {
    let model = immigration_death("crn-craft.rules-model/0.1");
    for version in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6"] {
        for quantiles in [json!([0.5]), json!([])] {
            let mut old = request(&model, Some(quantiles));
            old["schema"] = format!("crn-craft.rules-run/{version}").into();
            if version == "0.1" {
                old.as_object_mut().unwrap().remove("method");
            }
            let (code, path, _) = refusal(&model, &old);
            assert_eq!(
                (code.as_str(), path.as_str()),
                ("CRN-R020", "/schema"),
                "{version}"
            );
        }
    }
    // null is not an absent list.
    let model = immigration_death("crn-craft.model/1.0");
    let (code, _, _) = refusal(&model, &request(&model, Some(Value::Null)));
    assert_eq!(code, "CRN-E003");
}

#[test]
fn invalid_probability_lists_are_refused() {
    let model = immigration_death("crn-craft.model/1.0");
    let seventeen: Vec<f64> = (1..=17).map(|i| f64::from(i) / 18.0).collect();
    let sixteen: Vec<f64> = (1..=16).map(|i| f64::from(i) / 17.0).collect();
    assert!(try_run(&model, &request(&model, Some(json!(sixteen)))).is_ok());
    let cases = [
        (json!(seventeen), "/outputs/quantiles"),
        (json!([0.0]), "/outputs/quantiles/0"),
        (json!([1.0]), "/outputs/quantiles/0"),
        (json!([-0.5]), "/outputs/quantiles/0"),
        (json!([0.1, 1.5]), "/outputs/quantiles/1"),
        (json!([0.5, 0.5]), "/outputs/quantiles/1"),
        (json!([0.1, 0.6, 0.4]), "/outputs/quantiles/2"),
        (json!([0.1, 2e-320, 0.4]), "/outputs/quantiles/1"),
    ];
    for (quantiles, pointer) in cases {
        let (code, path, _) = refusal(&model, &request(&model, Some(quantiles.clone())));
        assert_eq!(
            (code.as_str(), path.as_str()),
            ("CRN-R020", pointer),
            "{quantiles}"
        );
    }
}

#[test]
fn the_quantile_buffer_counts_against_the_workspace_limit() {
    let model = immigration_death("crn-craft.model/1.0");
    // 1,000,000 replicates x 150 report times x 2 observables x 8 bytes =
    // 2.4 GB, above the 2 GiB workspace limit; without quantiles the request
    // needs about 0.3 GB.
    let mut large = request(&model, Some(json!([0.5])));
    large["ensemble"]["replicates"] = 1_000_000.into();
    large["outputs"]["keep_paths"] = 0.into();
    large["time"]["report_times"] =
        json!((0..150).map(|i| f64::from(i) * 0.02).collect::<Vec<_>>());
    large["budget"]["max_output_bytes"] = (1u64 << 30).into();
    let (code, path, hint) = refusal(&model, &large);
    assert_eq!(
        (code.as_str(), path.as_str()),
        ("CRN-R020", "/outputs/quantiles")
    );
    for words in ["fewer report times", "fewer replicates", "no quantiles"] {
        assert!(hint.contains(words), "{hint}");
    }
    large["outputs"]["quantiles"] = json!([]);
    assert!(prepare(&model, &serde_json::to_vec(&large).unwrap()).is_ok());
    // A smaller grid fits.
    large["outputs"]["quantiles"] = json!([0.5]);
    large["time"]["report_times"] =
        json!((0..100).map(|i| f64::from(i) * 0.02).collect::<Vec<_>>());
    assert!(prepare(&model, &serde_json::to_vec(&large).unwrap()).is_ok());
}

#[test]
fn the_1_0_schemas_describe_the_quantiles_and_the_rules_schemas_do_not() {
    let run = crn_craft::rules::schema("run-1.0").unwrap();
    let field = &run["definitions"]["Outputs"]["properties"]["quantiles"];
    assert_eq!(field["type"], "array");
    assert_eq!(field["maxItems"], 16);
    assert!(field["description"].as_str().unwrap().contains("type 7"));
    let result = crn_craft::rules::schema("result-1.0").unwrap();
    let field = &result["definitions"]["Statistics"]["properties"]["quantiles"];
    assert_eq!(field["items"]["$ref"], "#/definitions/Quantile");
    assert!(result["definitions"]["Quantile"]["properties"]["p"].is_object());
    let old_run = crn_craft::rules::schema("rules-run").unwrap();
    assert!(
        old_run["definitions"]["Outputs"]["properties"]
            .get("quantiles")
            .is_none()
    );
    let old_result = crn_craft::rules::schema("rules-result").unwrap();
    assert!(
        old_result["definitions"]["Statistics"]["properties"]
            .get("quantiles")
            .is_none()
    );
    assert!(old_result["definitions"].get("Quantile").is_none());
}
