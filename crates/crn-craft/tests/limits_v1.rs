//! M4a: the limits of `crn-craft.model/1.0`.
//!
//! A `model/1.0` document may have up to 4,096 types, 65,536 rules, 4,096
//! observables and 16,384 parameters, and its static lowering up to 4,096
//! instance tuples and 65,536 channels. Older model versions keep their
//! limits and diagnostics (64 types, 1,024 rules, 256 observables, 1,024
//! parameters; 256 tuples and 4,096 channels).
use crn_craft::diagnostic::Diagnostic;
use crn_craft::rules::convert::from_network;
use crn_craft::rules::{Control, check_model, execute, prepare};
use serde_json::{Value, json};
use std::path::PathBuf;
use std::time::{Duration, Instant};

const UNIFIED: &str = "crn-craft.model/1.0";
const LEGACY: &str = "crn-craft.rules-model/0.4";

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// A reversible chain X0 <-> X1 <-> ... with one observable per type: `types`
/// types, 2 (types - 1) rules.
fn chain(schema: &str, types: usize) -> Value {
    let mut rules = Vec::new();
    for i in 0..types - 1 {
        for (name, from, to) in [("f", i, i + 1), ("b", i + 1, i)] {
            rules.push(json!({
                "name": format!("{name}{i}"),
                "reactants": [{"label": "x", "type": format!("X{from}")}],
                "rate": "k",
                "products": [{"type": format!("X{to}")}]
            }));
        }
    }
    json!({
        "schema": schema,
        "name": "chain",
        "parameters": {"k": 1.0},
        "types": (0..types).map(|i| json!({"name": format!("X{i}")})).collect::<Vec<_>>(),
        "initial": [{"type": "X0", "count": 1000}],
        "observables": (0..types)
            .map(|i| (format!("X{i}"), json!({"type": format!("X{i}")})))
            .collect::<serde_json::Map<_, _>>(),
        "rules": rules
    })
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

fn request(model: &[u8], method: &str, end: f64) -> Vec<u8> {
    let schema = serde_json::from_slice::<Value>(model).unwrap()["schema"] == UNIFIED;
    bytes(&json!({
        "schema": if schema { "crn-craft.run/1.0" } else { "crn-craft.rules-run/0.6" },
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": [0.0, end]},
        "ensemble": {"replicates": 2, "root_seed": 11},
        "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1_000_000, "max_sampling_steps_per_replicate": 10_000_000,
                   "max_instances": 8192, "threads": 1, "max_output_bytes": 1u64 << 28},
        "method": {"kind": method}
    }))
}

fn has(diagnostics: &[Diagnostic], code: &str, path: &str, message: &str) -> bool {
    diagnostics
        .iter()
        .any(|d| d.code == code && d.path == path && d.message == message)
}

fn check(model: &Value) -> Vec<Diagnostic> {
    check_model(&bytes(model)).diagnostics
}

/// Declared in docs/campaign/merge-0.7/M4a-log.md before the first
/// measurement: `prepare` of the 4,096-species, 8,190-reaction chain with
/// `static_direct` takes under 2 s on the development host (release build).
const CHAIN_PREPARE: Duration = Duration::from_secs(2);

#[test]
fn a_4096_species_8190_reaction_chain_compiles_in_under_two_seconds() {
    let model = bytes(&chain(UNIFIED, 4096));
    assert_eq!(
        serde_json::from_slice::<Value>(&model).unwrap()["rules"]
            .as_array()
            .unwrap()
            .len(),
        8190
    );
    for method in ["static_direct", "next_reaction"] {
        let request = request(&model, method, 0.05);
        let start = Instant::now();
        let prepared = prepare(&model, &request).unwrap_or_else(|d| panic!("{d:?}"));
        let elapsed = start.elapsed();
        eprintln!("prepare {method}: {:.3} s", elapsed.as_secs_f64());
        // Timing is meaningful only for optimized builds.
        if !cfg!(debug_assertions) {
            assert!(
                elapsed < CHAIN_PREPARE,
                "{method}: prepare took {elapsed:?}"
            );
        }
        let result = execute(&prepared, &Control::default()).unwrap();
        assert_eq!(
            serde_json::to_value(result.schema).unwrap(),
            "crn-craft.result/1.0"
        );
        assert_eq!(result.rules.len(), 8190);
        assert_eq!(result.observables.len(), 4096);
        assert_eq!(result.completed_replicates, 2);
        for record in &result.replicates {
            // About 1000 molecules leave X0 at rate 1 for 0.05 time units.
            assert!(record.events > 10, "{record:?}");
            assert_eq!(record.rule_events.len(), 8190);
        }
        // Molecules are conserved along the chain.
        let total: f64 = result.statistics.as_ref().unwrap().mean.as_ref().unwrap()[1]
            .iter()
            .sum();
        assert_eq!(total, 1000.0);
    }
}

#[test]
fn model_1_0_accepts_its_limits_and_refuses_one_more() {
    // 4,096 types and observables, 65,536 rules and 16,384 parameters.
    let mut model = chain(UNIFIED, 4096);
    let rules = model["rules"].as_array_mut().unwrap();
    for i in rules.len()..65_536 {
        rules.push(
            json!({"name": format!("s{i}"), "reactants": [], "rate": "k",
                          "products": [{"type": "X0"}]}),
        );
    }
    for i in 0..16_383 {
        model["parameters"][format!("p{i}")] = json!(1.0);
    }
    let report = check_model(&bytes(&model));
    assert!(report.valid, "{:?}", &report.diagnostics[..1]);
    assert_eq!(
        serde_json::to_value(report.schema).unwrap(),
        "crn-craft.check/1.0"
    );

    let mut over = model.clone();
    over["types"]
        .as_array_mut()
        .unwrap()
        .push(json!({"name": "Y"}));
    assert!(has(
        &check(&over),
        "CRN-R015",
        "/types",
        "4097 types exceed the limit of 4096"
    ));
    let mut over = model.clone();
    over["rules"]
        .as_array_mut()
        .unwrap()
        .push(json!({"name": "extra", "reactants": [], "rate": "k", "products": []}));
    assert!(has(
        &check(&over),
        "CRN-R015",
        "/rules",
        "65537 rules exceed the limit of 65536"
    ));
    let mut over = model.clone();
    over["parameters"]["extra"] = json!(1.0);
    assert!(has(
        &check(&over),
        "CRN-R015",
        "/parameters",
        "16385 parameters exceed the limit of 16384"
    ));
    let mut over = model;
    over["observables"]["extra"] = json!({"type": "X0"});
    assert!(has(
        &check(&over),
        "CRN-R015",
        "/observables",
        "4097 observables exceed the limit of 4096"
    ));
}

#[test]
fn older_model_versions_keep_their_limits_and_diagnostics() {
    let model = chain(LEGACY, 64);
    assert!(check_model(&bytes(&model)).valid);
    let over = chain(LEGACY, 65);
    let diagnostics = check(&over);
    assert!(has(
        &diagnostics,
        "CRN-R015",
        "/types",
        "65 types exceed the limit of 64"
    ));
    assert_eq!(diagnostics.len(), 1, "{diagnostics:?}");
    let mut over = chain(LEGACY, 2);
    for i in 0..1024 {
        over["parameters"][format!("p{i}")] = json!(1.0);
    }
    for i in 0..1023 {
        over["rules"]
            .as_array_mut()
            .unwrap()
            .push(json!({"name": format!("s{i}"), "reactants": [], "rate": "k", "products": []}));
    }
    for i in 0..255 {
        over["observables"][format!("o{i}")] = json!({"type": "X0"});
    }
    let diagnostics = check(&over);
    assert!(has(
        &diagnostics,
        "CRN-R015",
        "/rules",
        "1025 rules exceed the limit of 1024"
    ));
    assert!(has(
        &diagnostics,
        "CRN-R015",
        "/parameters",
        "1025 parameters exceed the limit of 1024"
    ));
    assert!(has(
        &diagnostics,
        "CRN-R015",
        "/observables",
        "257 observables exceed the limit of 256"
    ));
    assert_eq!(diagnostics.len(), 3, "{diagnostics:?}");
}

/// A model whose static lowering has `tuples` instance tuples of type `A`
/// (attribute `i`) and `rules` decay rules over `A`: `tuples` × `rules`
/// channels, each affecting only the channels of its own tuple.
fn tuples(schema: &str, tuples: usize, rules: usize) -> Value {
    json!({
        "schema": schema,
        "name": "tuples",
        "parameters": {"k": 1.0},
        "types": [{"name": "A", "attributes": ["i"]}],
        "initial": (0..tuples).map(|i| json!({"type": "A", "attributes": {"i": i}, "count": 1}))
            .collect::<Vec<_>>(),
        "observables": {"A": {"type": "A"}},
        "rules": (0..rules).map(|r| json!({"name": format!("decay{r}"),
            "reactants": [{"label": "a", "type": "A"}], "rate": "k", "products": []}))
            .collect::<Vec<_>>()
    })
}

fn lowering(model: &Value) -> Result<(), Vec<Diagnostic>> {
    let model = bytes(model);
    prepare(&model, &request(&model, "static_direct", 1.0)).map(|_| ())
}

#[test]
fn static_lowering_limits_follow_the_model_version() {
    // Instance tuples: 256 before 1.0, 4,096 in 1.0.
    assert!(lowering(&tuples(LEGACY, 256, 1)).is_ok());
    let refused = lowering(&tuples(LEGACY, 257, 1)).unwrap_err();
    assert!(has(
        &refused,
        "CRN-R020",
        "/method",
        "static backend supports up to 256 instances"
    ));
    assert!(lowering(&tuples(UNIFIED, 4096, 1)).is_ok());
    let refused = lowering(&tuples(UNIFIED, 4097, 1)).unwrap_err();
    assert!(has(
        &refused,
        "CRN-R020",
        "/method",
        "static backend supports up to 4096 instances"
    ));
    // Channels: 4,096 before 1.0, 65,536 in 1.0.
    assert!(lowering(&tuples(LEGACY, 256, 16)).is_ok());
    let refused = lowering(&tuples(LEGACY, 256, 17)).unwrap_err();
    assert!(has(
        &refused,
        "CRN-R020",
        "/method",
        "static channel limit exceeded"
    ));
    assert!(lowering(&tuples(UNIFIED, 4096, 16)).is_ok());
    let refused = lowering(&tuples(UNIFIED, 4096, 17)).unwrap_err();
    assert!(has(
        &refused,
        "CRN-R020",
        "/method",
        "static channel limit exceeded"
    ));
}

#[test]
fn the_converter_writes_models_up_to_the_1_0_limits() {
    for workload in ["W5", "W7"] {
        let network = std::fs::read(root().join(format!(
            "validation/ssa-acceleration/models/{workload}.model.json"
        )))
        .unwrap();
        let converted = from_network(&network).unwrap_or_else(|d| panic!("{workload}: {d:?}"));
        let report = check_model(&converted);
        assert!(report.valid, "{workload}: {:?}", report.diagnostics);
        let value: Value = serde_json::from_slice(&converted).unwrap();
        assert_eq!(value["schema"], UNIFIED);
    }
    // The network contract's own limits: 4,096 species.
    let species: Vec<Value> = (0..4096)
        .map(|i| json!({"name": format!("S{i}"), "initial_count": 1}))
        .collect();
    let network = json!({
        "schema": "crn-craft.model/0.1", "name": "wide", "convention": "stochastic_counts",
        "omega": 1.0, "parameters": {"k": 1.0}, "species": species,
        "reactions": [{"name": "decay", "reactants": {"S4095": 1}, "products": {},
                       "rate": {"kind": "mass_action", "coefficient": "k"}}]
    });
    let converted = from_network(&bytes(&network)).unwrap_or_else(|d| panic!("{d:?}"));
    assert!(check_model(&converted).valid);
}
