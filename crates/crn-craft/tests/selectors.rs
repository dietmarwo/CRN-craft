//! M4c: event selectors of Static Direct in `crn-craft.run/1.0`.
//!
//! `method.static_direct.selector` chooses how Static Direct picks the next
//! channel: `linear` (the scan of every rules-run version),
//! `composition_rejection`, `sorting`, `rssa_cr`, or `auto` (the default in
//! run/1.0), which picks one of them from the lowered network. Every
//! `rules-run` version runs the linear scan and refuses the field. The
//! sampler identity names the selector that ran; `linear` keeps the 0.6.0
//! identities, so a run/1.0 request with `"selector": "linear"` gives the
//! rules-run/0.6 result.
use crn_craft::diagnostic::Diagnostic;
use crn_craft::rules::convert::from_network;
use crn_craft::rules::{
    Control, LoopHalt, LoopHost, RulesLoop, RulesResult, execute, execute_with_rules_loop, prepare,
};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::PathBuf;

const SELECTORS: [&str; 4] = ["linear", "composition_rejection", "sorting", "rssa_cr"];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn read(path: &str) -> Vec<u8> {
    std::fs::read(root().join(path)).unwrap()
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

/// `bytes` with its one `"schema": "..."` value replaced.
fn with_schema(bytes: &[u8], schema: &str) -> Vec<u8> {
    let text = String::from_utf8(bytes.to_vec()).unwrap();
    let old = serde_json::from_slice::<Value>(bytes).unwrap()["schema"]
        .as_str()
        .unwrap()
        .to_string();
    text.replacen(&format!("\"{old}\""), &format!("\"{schema}\""), 1)
        .into_bytes()
}

fn unified(path: &str) -> Vec<u8> {
    with_schema(&read(path), "crn-craft.model/1.0")
}

fn converted(path: &str) -> Vec<u8> {
    from_network(&read(path)).unwrap_or_else(|d| panic!("{d:?}"))
}

#[allow(clippy::too_many_arguments)]
fn request(
    model: &[u8],
    schema: &str,
    method: Value,
    end: f64,
    report_times: &[f64],
    replicates: u64,
    keep: u64,
    threads: u32,
) -> Value {
    json!({
        "schema": schema,
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": report_times},
        "ensemble": {"replicates": replicates, "root_seed": 11},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": keep, "keep_instances": 1},
        "budget": {"max_events_per_replicate": 10_000_000,
                   "max_sampling_steps_per_replicate": 1_000_000_000,
                   "max_instances": 4096, "threads": threads, "max_output_bytes": 1u64 << 30},
        "method": method
    })
}

fn static_direct(selector: &str) -> Value {
    json!({"kind": "static_direct", "selector": selector})
}

fn run(model: &[u8], request: &Value) -> RulesResult {
    let prepared = prepare(model, &bytes(request)).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap()
}

fn refused(model: &[u8], request: &Value) -> Vec<Diagnostic> {
    match prepare(model, &bytes(request)) {
        Ok(_) => panic!("prepared: {request}"),
        Err(d) => d,
    }
}

/// The result as JSON without `run_sha256`, for comparisons of runs that
/// differ only in settings that must not change the numbers.
fn numbers(result: &RulesResult) -> Value {
    let mut value = serde_json::to_value(result).unwrap();
    value.as_object_mut().unwrap().remove("run_sha256");
    value
}

#[test]
fn selector_is_a_run_1_0_field() {
    let model = unified("models/rules/immigration-death.rules.json");
    for selector in [
        "linear",
        "composition_rejection",
        "sorting",
        "rssa_cr",
        "auto",
    ] {
        let request = request(
            &model,
            "crn-craft.run/1.0",
            static_direct(selector),
            1.0,
            &[0.0, 1.0],
            4,
            0,
            1,
        );
        prepare(&model, &bytes(&request)).unwrap_or_else(|d| panic!("{selector}: {d:?}"));
    }
    // Every rules-run version runs the linear scan and refuses the field.
    let rules_model = read("models/rules/immigration-death.rules.json");
    for version in ["0.2", "0.3", "0.4", "0.5", "0.6"] {
        for selector in ["linear", "auto", "composition_rejection"] {
            let request = request(
                &rules_model,
                &format!("crn-craft.rules-run/{version}"),
                static_direct(selector),
                1.0,
                &[0.0, 1.0],
                4,
                0,
                1,
            );
            let diagnostics = refused(&rules_model, &request);
            assert!(
                diagnostics.iter().any(|d| d.code == "CRN-R020"
                    && d.path == "/schema"
                    && d.message == "static_direct.selector requires crn-craft.run/1.0"),
                "{version}/{selector}: {diagnostics:?}"
            );
        }
    }
    // A present field holds a value of the enumeration.
    for selector in [json!(null), json!("tree"), json!(1), json!("Linear")] {
        let request = request(
            &model,
            "crn-craft.run/1.0",
            json!({"kind": "static_direct", "selector": selector}),
            1.0,
            &[0.0, 1.0],
            4,
            0,
            1,
        );
        assert!(
            refused(&model, &request)
                .iter()
                .any(|d| d.code == "CRN-E003"),
            "{selector}"
        );
    }
    // Only Static Direct has a selector.
    for kind in ["next_reaction", "tree", "mtdna_direct"] {
        let request = request(
            &model,
            "crn-craft.run/1.0",
            json!({"kind": kind, "selector": "linear"}),
            1.0,
            &[0.0, 1.0],
            4,
            0,
            1,
        );
        assert!(
            refused(&model, &request)
                .iter()
                .any(|d| d.code == "CRN-E003"),
            "{kind}"
        );
    }
}

#[test]
fn the_sampler_names_the_selector() {
    let model = unified("models/rules/immigration-death.rules.json");
    let names = [
        ("linear", "rules-static-direct"),
        ("composition_rejection", "rules-static-direct-cr"),
        ("sorting", "rules-static-direct-sorting"),
        ("rssa_cr", "rules-static-direct-rssa-cr"),
    ];
    for (selector, name) in names {
        for (wait, prefix) in [
            ("inverse", "chacha8-u53"),
            ("ziggurat", "chacha8-u53-ziggurat"),
        ] {
            for sparse in [false, true] {
                let method = json!({"kind": "static_direct", "selector": selector,
                                    "waiting_time": wait, "sparse_trackers": sparse});
                let mut request = request(
                    &model,
                    "crn-craft.run/1.0",
                    method,
                    2.0,
                    &[0.0, 2.0],
                    8,
                    2,
                    1,
                );
                request["trackers"] = json!([{"observable": "X"}]);
                let result = run(&model, &request);
                let suffix = if sparse { "-sparse-trackers" } else { "" };
                assert_eq!(result.sampler, format!("{prefix}/{name}{suffix}/1"));
                assert_eq!(result.completed_replicates, 8);
            }
        }
    }
}

/// The run/1.0 result of a request with `"selector": "linear"` is the
/// rules-run/0.6 result of the same request without the field, apart from
/// the schema string, `run_sha256` and the echoed method's `selector`.
#[test]
fn linear_under_run_1_0_equals_rules_run_0_6() {
    let stored = root().join("crates/crn-craft/tests/stored");
    let mut names: Vec<String> = std::fs::read_dir(&stored)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.ends_with(".run.json"))
        .collect();
    names.sort();
    let mut compared = 0;
    let mut cases: Vec<(String, Vec<u8>, Value)> = Vec::new();
    for name in &names {
        let request: Value =
            serde_json::from_slice(&std::fs::read(stored.join(name)).unwrap()).unwrap();
        if request["method"]["kind"] != "static_direct" {
            continue;
        }
        let model = read(&format!(
            "models/rules/{}",
            request["model"]["path"].as_str().unwrap()
        ));
        cases.push((name.clone(), model, request));
    }
    // Converted network models with several channels per event.
    for path in [
        "validation/ssa-acceleration/models/W2.model.json",
        "validation/ssa-acceleration/models/W3.model.json",
        "models/legacy/enzyme-mechanism.model.json",
    ] {
        let model = with_schema(&converted(path), "crn-craft.rules-model/0.4");
        let request = request(
            &model,
            "crn-craft.rules-run/0.6",
            json!({"kind": "static_direct"}),
            5.0,
            &[0.0, 1.0, 5.0],
            16,
            4,
            2,
        );
        cases.push((path.to_string(), model, request));
    }
    for (name, model, request) in cases {
        let mut old = request.clone();
        old["schema"] = "crn-craft.rules-run/0.6".into();
        let mut new = request.clone();
        new["schema"] = "crn-craft.run/1.0".into();
        new["method"]["selector"] = "linear".into();
        let reference = numbers(&run(&model, &old));
        let mut result = numbers(&run(&model, &new));
        assert_eq!(result["schema"], "crn-craft.result/1.0", "{name}");
        assert_eq!(result["method"]["selector"], "linear", "{name}");
        result["method"].as_object_mut().unwrap().remove("selector");
        result["schema"] = reference["schema"].clone();
        assert!(
            result == reference,
            "{name}: linear run/1.0 differs from rules-run/0.6"
        );
        compared += 1;
    }
    // Ten stored Static Direct requests and three converted models.
    assert_eq!(compared, 13);
}

/// Mean and variance of a Poisson law with mean `mu` from `n` replicates:
/// 5 standard errors (the variance's from the fourth moment).
fn assert_poisson(name: &str, mean: f64, variance: f64, mu: f64, n: f64) {
    let se_mean = (mu / n).sqrt();
    let se_variance = ((mu + 3.0 * mu * mu - mu * mu * (n - 3.0) / (n - 1.0)) / n).sqrt();
    assert!(
        (mean - mu).abs() <= 5.0 * se_mean,
        "{name}: mean {mean}, law {mu}"
    );
    assert!(
        (variance - mu).abs() <= 5.0 * se_variance,
        "{name}: variance {variance}, law {mu}"
    );
}

/// A quick law check for every selector and both waits; the declared
/// exact-law check is `validation/merge/check_selectors.py`.
#[test]
fn every_selector_samples_the_immigration_death_law() {
    let model = unified("models/rules/immigration-death.rules.json");
    let n = 20_000.0;
    for selector in SELECTORS {
        for wait in ["inverse", "ziggurat"] {
            let method =
                json!({"kind": "static_direct", "selector": selector, "waiting_time": wait});
            let result = run(
                &model,
                &request(
                    &model,
                    "crn-craft.run/1.0",
                    method,
                    3.0,
                    &[0.0, 0.5, 3.0],
                    20_000,
                    0,
                    2,
                ),
            );
            let statistics = result.statistics.unwrap();
            let (mean, variance) = (
                statistics.mean.unwrap(),
                statistics.sample_variance.unwrap(),
            );
            for (row, t) in [(1, 0.5), (2, 3.0)] {
                let mu = 10.0 * (1.0 - f64::exp(-t));
                assert_poisson(
                    &format!("{selector}/{wait}/t={t}"),
                    mean[row][0],
                    variance[row][0],
                    mu,
                    n,
                );
            }
        }
    }
}

#[test]
fn selectors_do_not_depend_on_the_thread_count() {
    for (path, end) in [
        ("validation/ssa-acceleration/models/W5.model.json", 0.5),
        ("validation/ssa-acceleration/models/W4.model.json", 2.0),
        ("models/legacy/repressilator.model.json", 5.0),
    ] {
        let model = converted(path);
        for selector in SELECTORS {
            let results: Vec<Value> = [1, 2, 3]
                .iter()
                .map(|&threads| {
                    numbers(&run(
                        &model,
                        &request(
                            &model,
                            "crn-craft.run/1.0",
                            static_direct(selector),
                            end,
                            &[0.0, end],
                            12,
                            12,
                            threads,
                        ),
                    ))
                })
                .collect();
            assert_eq!(results[0], results[1], "{path}/{selector}");
            assert_eq!(results[0], results[2], "{path}/{selector}");
        }
    }
}

/// A model of every channel shape and rate form the selectors handle:
/// sources, single, pair and cross channels, a higher-order rule with a
/// combination count, whole-rule propensities that read one species (a
/// Michaelis–Menten law) and two species (pinned in RSSA-CR), and attributed
/// species with factors. Invariants hold at every report time.
fn mixed_model() -> Value {
    json!({
        "schema": "crn-craft.model/1.0",
        "name": "mixed",
        "parameters": {"k": 0.5, "kd": 0.002, "vmax": 4.0, "km": 30.0, "c": 0.01},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C"}, {"name": "D"},
                  {"name": "P", "attributes": ["size"]}],
        "initial": [{"type": "A", "count": 60}, {"type": "B", "count": 25},
                    {"type": "C", "count": 0}, {"type": "D", "count": 3},
                    {"type": "P", "attributes": {"size": 1}, "count": 10},
                    {"type": "P", "attributes": {"size": 2}, "count": 4}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}, "C": {"type": "C"},
                        "D": {"type": "D"}, "P": {"type": "P"},
                        "mass": {"type": "P", "weight": "x.size"}},
        "rules": [
            {"name": "make_a", "reactants": [], "rate": "k * 10", "products": [{"type": "A"}]},
            {"name": "a_to_b", "reactants": [{"label": "x", "type": "A"}], "rate": "0.05",
             "products": [{"type": "B"}]},
            {"name": "dimer", "reactants": [{"label": "x", "type": "B"}, {"label": "y", "type": "B"}],
             "rate": "kd", "products": [{"type": "C"}]},
            {"name": "split", "reactants": [{"label": "x", "type": "C"}], "rate": "0.3",
             "products": [{"type": "B"}, {"type": "B"}]},
            {"name": "trimer", "reactants": [{"label": "x", "type": "A", "count": 3}],
             "rate": "0.00001", "products": [{"type": "D"}]},
            {"name": "mm", "reactants": [{"label": "x", "type": "A"}],
             "propensity": "vmax * A / (km + A)", "products": [{"type": "C"}]},
            {"name": "ratio", "reactants": [{"label": "x", "type": "D"}],
             "propensity": "c * D * B / (B + C + 1)", "products": [{"type": "A"}]},
            {"name": "cross", "reactants": [{"label": "x", "type": "B"},
                                            {"label": "p", "type": "P", "factor": "p.size"}],
             "rate": "0.001",
             "products": [{"type": "A"}, {"type": "P", "copy": "p"}]},
            {"name": "decay", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "B"},
                                            {"label": "z", "type": "D"}],
             "rate": "0.0001", "products": [{"type": "D"}]}
        ]
    })
}

#[test]
fn every_selector_runs_every_channel_shape_with_trackers_and_logs() {
    let model = bytes(&mixed_model());
    let check = crn_craft::rules::check_model(&model);
    assert!(check.diagnostics.is_empty(), "{:?}", check.diagnostics);
    let mut means = BTreeMap::new();
    for selector in SELECTORS {
        for (wait, sparse) in [("inverse", false), ("ziggurat", true)] {
            let method = json!({"kind": "static_direct", "selector": selector,
                                "waiting_time": wait, "sparse_trackers": sparse});
            let mut request = request(
                &model,
                "crn-craft.run/1.0",
                method,
                20.0,
                &[0.0, 5.0, 10.0, 20.0],
                400,
                400,
                2,
            );
            request["trackers"] = json!([{"observable": "B", "first_passage": "B > 30"}]);
            request["event_log"] = json!({"max_entries": 8});
            let result = run(&model, &request);
            assert_eq!(result.completed_replicates, 400, "{selector}");
            let column = |name: &str| result.observables.iter().position(|o| o == name).unwrap();
            for path in &result.paths {
                for row in &path.values {
                    // The P instances only change identity, never count.
                    assert_eq!(row[column("P")], 14.0, "{selector}");
                    assert_eq!(row[column("mass")], 18.0, "{selector}");
                }
            }
            assert!(result.events > 0);
            means.insert(
                format!("{selector}/{wait}"),
                result.statistics.unwrap().mean.unwrap()[3].clone(),
            );
        }
    }
    // Every selector samples one law: the means at t = 20 agree within a
    // generous tolerance (each from 400 replicates).
    let reference = &means["linear/inverse"];
    for (name, mean) in &means {
        for (a, b) in mean.iter().zip(reference) {
            assert!(
                (a - b).abs() <= 0.15 * b.abs().max(5.0),
                "{name}: {mean:?} vs {reference:?}"
            );
        }
    }
}

/// A run whose propensities all vanish while some RSSA-CR upper bound is
/// positive ends at the horizon, like the linear scan.
#[test]
fn absorbing_states_end_every_selector() {
    let model = bytes(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "absorb",
        "parameters": {},
        "types": [{"name": "A"}, {"name": "B"}],
        "initial": [{"type": "A", "count": 6}, {"type": "B", "count": 0}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [
            {"name": "pair", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "A"}],
             "rate": "1", "products": [{"type": "B"}]}
        ]
    }));
    for selector in SELECTORS {
        let result = run(
            &model,
            &request(
                &model,
                "crn-craft.run/1.0",
                static_direct(selector),
                1e6,
                &[0.0, 1e6],
                20,
                20,
                1,
            ),
        );
        assert_eq!(result.completed_replicates, 20);
        for path in &result.paths {
            assert_eq!(path.values[1], vec![0.0, 3.0], "{selector}");
        }
        assert_eq!(result.events, 60, "{selector}");
    }
}

/// A positive channel below 2^-53 of the total: the linear scan refuses it
/// with `CRN-R031`, and so do the other selectors (their totals check the
/// smallest positive rate the same way).
#[test]
fn unresolvable_channels_fail_with_crn_r031() {
    let model = bytes(&json!({
        "schema": "crn-craft.model/1.0",
        "name": "tiny",
        "parameters": {},
        "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": 1}],
        "observables": {"A": {"type": "A"}},
        "rules": [
            {"name": "big", "reactants": [], "rate": "1e20", "products": []},
            {"name": "small", "reactants": [{"label": "x", "type": "A"}], "rate": "1e-3", "products": []}
        ]
    }));
    for selector in SELECTORS {
        let result = run(
            &model,
            &request(
                &model,
                "crn-craft.run/1.0",
                static_direct(selector),
                1.0,
                &[0.0, 1.0],
                1,
                0,
                1,
            ),
        );
        let diagnostic = result.replicates[0].diagnostic.as_ref().unwrap();
        assert_eq!(diagnostic.code, "CRN-R031", "{selector}: {diagnostic:?}");
    }
}

struct Stale;

impl RulesLoop for Stale {
    fn fingerprint(&self) -> &str {
        "stale"
    }
    fn run(
        &self,
        _host: &mut LoopHost<'_, '_>,
        _counts: &mut [u64],
        _values: &mut [f64],
    ) -> Result<Option<f64>, LoopHalt> {
        unreachable!("a refused loop never runs")
    }
}

/// Generated event loops implement the linear scan: a run/1.0 request whose
/// selector is another one, or `auto` resolving to another one, is refused
/// with `CRN-R020` at `/method/selector` before the fingerprint is checked.
#[test]
fn generated_loops_refuse_other_selectors() {
    let model = unified("models/rules/immigration-death.rules.json");
    for selector in ["composition_rejection", "sorting", "rssa_cr"] {
        let prepared = prepare(
            &model,
            &bytes(&request(
                &model,
                "crn-craft.run/1.0",
                static_direct(selector),
                1.0,
                &[0.0, 1.0],
                2,
                0,
                1,
            )),
        )
        .unwrap();
        let error = execute_with_rules_loop(&prepared, &Control::default(), &Stale).unwrap_err();
        assert_eq!(
            (error.code.as_str(), error.path.as_str()),
            ("CRN-R020", "/method/selector")
        );
        // Generation does not depend on the selector.
        let error = crn_craft::rules::generate_rules_loop(&prepared).err();
        assert!(error.is_none(), "{selector}: {error:?}");
    }
    // `linear` and an `auto` that picks the linear scan pass the selector
    // check and fail at the fingerprint.
    for selector in [Some("linear"), Some("auto"), None] {
        let method = selector.map_or(json!({"kind": "static_direct"}), static_direct);
        let prepared = prepare(
            &model,
            &bytes(&request(
                &model,
                "crn-craft.run/1.0",
                method,
                1.0,
                &[0.0, 1.0],
                2,
                0,
                1,
            )),
        )
        .unwrap();
        let error = execute_with_rules_loop(&prepared, &Control::default(), &Stale).unwrap_err();
        assert_eq!(
            (error.code.as_str(), error.path.as_str()),
            ("CRN-R020", "/model")
        );
    }
    // `auto` on a large network picks another selector: refused.
    let model = converted("validation/ssa-acceleration/models/W6.model.json");
    let prepared = prepare(
        &model,
        &bytes(&request(
            &model,
            "crn-craft.run/1.0",
            json!({"kind": "static_direct"}),
            0.1,
            &[0.0, 0.1],
            2,
            0,
            1,
        )),
    )
    .unwrap();
    let error = execute_with_rules_loop(&prepared, &Control::default(), &Stale).unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-R020", "/method/selector")
    );
}

#[test]
fn the_1_0_schemas_describe_the_selector_and_the_rules_schemas_do_not() {
    let selector = |kind: &str| {
        let schema = crn_craft::rules::schema(kind).unwrap();
        let methods = schema["definitions"]["Method"]["oneOf"]
            .as_array()
            .unwrap()
            .clone();
        let direct = methods
            .iter()
            .find(|m| m["properties"]["kind"]["enum"][0] == "static_direct")
            .unwrap()
            .clone();
        (schema, direct["properties"].get("selector").cloned())
    };
    for kind in ["run-1.0", "result-1.0"] {
        let (schema, property) = selector(kind);
        assert!(property.is_some(), "{kind}");
        let values = &schema["definitions"]["Selector"]["enum"];
        assert_eq!(
            values,
            &json!([
                "linear",
                "composition_rejection",
                "sorting",
                "rssa_cr",
                "auto"
            ]),
            "{kind}"
        );
    }
    for kind in ["rules-run", "rules-result"] {
        let (schema, property) = selector(kind);
        assert!(property.is_none(), "{kind}");
        assert!(schema["definitions"].get("Selector").is_none(), "{kind}");
    }
    let discovery = crn_craft::rules::describe();
    assert_eq!(
        discovery["contracts_1_0"]["selectors"]["values"],
        json!([
            "linear",
            "composition_rejection",
            "sorting",
            "rssa_cr",
            "auto"
        ])
    );
}

/// A ring of `n` species and channels X_i -> X_(i+1) with `count` molecules
/// each; with `two_reads` each rate reads X_i and X_(i+1).
fn ring(n: usize, two_reads: bool, count: u64) -> Vec<u8> {
    let names: Vec<String> = (0..n).map(|i| format!("X{i}")).collect();
    let rules: Vec<Value> = (0..n)
        .map(|i| {
            let (x, y) = (&names[i], &names[(i + 1) % n]);
            let mut rule = json!({"name": format!("r{i}"), "reactants": [{"label": "x", "type": x}],
                                  "products": [{"type": y}]});
            if two_reads {
                rule["propensity"] = format!("{x} * {y} / (1 + {y})").into();
            } else {
                rule["rate"] = "1".into();
            }
            rule
        })
        .collect();
    bytes(&json!({
        "schema": "crn-craft.model/1.0", "name": "ring", "parameters": {},
        "types": names.iter().map(|x| json!({"name": x})).collect::<Vec<_>>(),
        "initial": names.iter().map(|x| json!({"type": x, "count": count})).collect::<Vec<_>>(),
        "observables": names.iter().map(|x| (x.clone(), json!({"type": x}))).collect::<serde_json::Map<_, _>>(),
        "rules": rules
    }))
}

/// The rule fixed in M4f from `validation/merge/selector-rule.json`: RSSA-CR
/// from 128 channels on when the mean initial count of the species that
/// start with molecules is below 8, from 8 channels on otherwise (also when
/// every species starts at zero), unless a rule reads two or more species;
/// else composition–rejection from 128 channels on; else the linear scan.
#[test]
fn auto_follows_the_declared_threshold() {
    let sampler = |model: &[u8], method: Value| {
        run(
            model,
            &request(
                model,
                "crn-craft.run/1.0",
                method,
                0.01,
                &[0.0, 0.01],
                1,
                0,
                1,
            ),
        )
        .sampler
    };
    let linear = "chacha8-u53-ziggurat/rules-static-direct/1";
    let cr = "chacha8-u53-ziggurat/rules-static-direct-cr/1";
    let rssa = "chacha8-u53-ziggurat/rules-static-direct-rssa-cr/1";
    for (n, two_reads, count, expected) in [
        (7, false, 8, linear),
        (8, false, 8, rssa),
        (8, false, 7, linear),
        (127, false, 7, linear),
        (128, false, 7, rssa),
        (127, true, 100, linear),
        (128, true, 100, cr),
        (2048, true, 3, cr),
        (2048, false, 3, rssa),
    ] {
        let model = ring(n, two_reads, count);
        assert_eq!(
            sampler(&model, json!({"kind": "static_direct"})),
            expected,
            "{n} {two_reads} {count}"
        );
        assert_eq!(sampler(&model, static_direct("auto")), expected, "{n}");
    }
    for (name, expected) in [
        ("W1", linear),
        ("W2", rssa),
        ("W3", linear),
        ("W4", rssa),
        ("W5", rssa),
        ("W6", rssa),
        ("W7", rssa),
    ] {
        let model = converted(&format!(
            "validation/ssa-acceleration/models/{name}.model.json"
        ));
        assert_eq!(
            sampler(&model, json!({"kind": "static_direct"})),
            expected,
            "{name}"
        );
    }
    for (name, expected) in [
        ("enzyme-mechanism", linear),
        ("immigration-death", linear),
        ("michaelis-menten", linear),
        ("repressilator", linear),
        ("vilar-noise", linear),
    ] {
        let model = converted(&format!("models/legacy/{name}.model.json"));
        assert_eq!(
            sampler(&model, json!({"kind": "static_direct"})),
            expected,
            "{name}"
        );
    }
    let discovery = crn_craft::rules::describe();
    assert_eq!(
        discovery["contracts_1_0"]["selectors"]["auto"],
        json!({
            "rssa_cr_from_channels": {"population_below_cut": 128, "otherwise": 8},
            "population_cut": 8.0,
            "composition_rejection_from_channels": 128
        })
    );
}
