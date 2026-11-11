//! The converter from network models (`crn-craft.model/0.1`) to rules models.
//!
//! Every example model converts to a rules model that passes
//! `rules::check_model`, conversion is deterministic, and every model the
//! rules language cannot hold yet is refused with a coded diagnostic.
use crn_craft::diagnostic::{Diagnostic, Severity};
use crn_craft::rules::convert::from_network;
use crn_craft::rules::{Control, Termination, check_model, execute, prepare};
use serde_json::{Value, json};
use std::path::PathBuf;

/// The schema string the converter writes.
const CONVERTED: &str = "crn-craft.model/1.0";
/// The converter's refusal code.
const REFUSED: &str = "CRN-R016";

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn read(path: &str) -> Vec<u8> {
    std::fs::read(root().join(path)).unwrap()
}

fn convert(model: &Value) -> Result<Value, Vec<Diagnostic>> {
    from_network(&serde_json::to_vec(model).unwrap())
        .map(|bytes| serde_json::from_slice(&bytes).unwrap())
}

fn refused(model: &Value) -> Vec<Diagnostic> {
    convert(model).expect_err("the converter must refuse this model")
}

fn has(diagnostics: &[Diagnostic], code: &str, path: &str) -> bool {
    diagnostics
        .iter()
        .any(|d| d.code == code && d.path == path && d.severity == Severity::Error)
}

/// A small valid network model: `A -> B` at rate `k`, `A + B -> C` at 0.5.
fn base() -> Value {
    json!({
        "schema": "crn-craft.model/0.1", "name": "base",
        "convention": "stochastic_counts", "omega": 1.0,
        "parameters": {"k": 0.25},
        "species": [{"name": "A", "initial_count": 10},
                    {"name": "B", "initial_count": 0},
                    {"name": "C", "initial_count": 3}],
        "reactions": [
            {"name": "convert", "reactants": {"A": 1}, "products": {"B": 1},
             "rate": {"kind": "mass_action", "coefficient": "k"}},
            {"name": "bind", "reactants": {"A": 1, "B": 1}, "products": {"C": 1},
             "rate": {"kind": "mass_action", "coefficient": 0.5}}]
    })
}

fn species(n: usize) -> Value {
    (0..n)
        .map(|i| json!({"name": format!("S{i}"), "initial_count": 1}))
        .collect()
}

fn decay(name: &str, species: &str) -> Value {
    json!({"name": name, "reactants": {species: 1}, "products": {},
           "rate": {"kind": "mass_action", "coefficient": 1.0}})
}

#[test]
fn every_example_and_m2_workload_converts_to_a_valid_rules_model() {
    let mut paths: Vec<String> = std::fs::read_dir(root().join("models/examples"))
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.ends_with(".model.json"))
        .map(|name| format!("models/examples/{name}"))
        .collect();
    paths.sort();
    assert_eq!(paths.len(), 5, "{paths:?}");
    for workload in ["W1", "W2", "W3", "W6"] {
        paths.push(format!(
            "validation/ssa-acceleration/models/{workload}.model.json"
        ));
    }
    for path in &paths {
        let network = read(path);
        let first = from_network(&network).unwrap_or_else(|d| panic!("{path}: {d:?}"));
        let second = from_network(&network).unwrap();
        assert_eq!(first, second, "{path}: conversion is not deterministic");
        let report = check_model(&first);
        assert!(report.valid, "{path}: {:?}", report.diagnostics);
        let converted: Value = serde_json::from_slice(&first).unwrap();
        let source: Value = serde_json::from_slice(&network).unwrap();
        assert_eq!(converted["schema"], CONVERTED, "{path}");
        assert_eq!(converted["name"], source["name"], "{path}");
        assert_eq!(converted["parameters"], source["parameters"], "{path}");
        // One attribute-free type and one observable per species, in order.
        let names: Vec<&Value> = source["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| &s["name"])
            .collect();
        let types: Vec<&Value> = converted["types"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| &t["name"])
            .collect();
        assert_eq!(types, names, "{path}");
        for name in &names {
            let name = name.as_str().unwrap();
            assert_eq!(
                converted["observables"][name],
                json!({"type": name, "reduce": "sum"}),
                "{path}"
            );
        }
        let rules = converted["rules"].as_array().unwrap();
        assert_eq!(rules.len(), source["reactions"].as_array().unwrap().len());
        // Canonical bytes: the rules model type's pretty JSON and a newline.
        let typed: crn_craft::rules::Model = serde_json::from_slice(&first).unwrap();
        assert!(
            first == [serde_json::to_vec_pretty(&typed).unwrap(), b"\n".to_vec()].concat(),
            "{path}: not canonical"
        );
    }
}

#[test]
fn stoichiometry_rates_and_initial_counts_map_as_documented() {
    let mut model = base();
    model["reactions"].as_array_mut().unwrap().extend([
        json!({"name": "dimerize", "reactants": {"A": 2}, "products": {"C": 2},
               "rate": {"kind": "mass_action", "coefficient": 1e-7}}),
        json!({"name": "make", "reactants": {}, "products": {"A": 1},
               "rate": {"kind": "hill_activation", "regulator": "C", "vmax": "k",
                        "basal": 0.0, "k": 40.0, "n": 2.0}}),
        json!({"name": "repress", "reactants": {"C": 1}, "products": {"A": 1, "C": 1},
               "rate": {"kind": "hill_repression", "regulator": "B", "vmax": 2.5,
                        "basal": "k", "k": "k", "n": 3.5}}),
        json!({"name": "enzyme", "reactants": {"A": 1}, "products": {"B": 1},
               "rate": {"kind": "michaelis_menten", "substrate": "A", "vmax": 3.0,
                        "km": "k"}}),
    ]);
    let converted = convert(&model).unwrap();
    assert_eq!(
        converted["initial"],
        json!([{"type": "A", "attributes": {}, "count": 10.0},
               {"type": "C", "attributes": {}, "count": 3.0}]),
        "zero counts are omitted"
    );
    let rules = converted["rules"].as_array().unwrap();
    let pattern = |label: &str, ty: &str| json!({"label": label, "type": ty});
    let product = |ty: &str| json!({"type": ty, "attributes": {}});
    assert_eq!(rules[0]["reactants"], json!([pattern("x", "A")]));
    assert_eq!(rules[0]["rate"], "k");
    assert_eq!(rules[0]["products"], json!([product("B")]));
    assert_eq!(
        rules[1]["reactants"],
        json!([pattern("x", "A"), pattern("y", "B")])
    );
    assert_eq!(rules[1]["rate"], "0.5");
    assert_eq!(
        rules[2]["reactants"],
        json!([pattern("x", "A"), pattern("y", "A")])
    );
    assert_eq!(rules[2]["rate"], "1e-7");
    assert_eq!(rules[2]["products"], json!([product("C"), product("C")]));
    assert!(rules[2].get("propensity").is_none());
    assert_eq!(rules[3]["reactants"], json!([]));
    assert_eq!(rules[3]["propensity"], "0.0 + k / (1 + (40.0 / C)^2.0)");
    assert!(rules[3].get("rate").is_none());
    assert_eq!(rules[4]["propensity"], "k + 2.5 / (1 + (B / k)^3.5)");
    assert_eq!(rules[4]["products"], json!([product("A"), product("C")]));
    assert_eq!(rules[5]["propensity"], "3.0 * (A / (k + A))");
    assert!(check_model(&from_network(&serde_json::to_vec(&model).unwrap()).unwrap()).valid);
}

#[test]
fn literal_numbers_round_trip_exactly() {
    for value in [
        0.0,
        5e-324,
        f64::MIN_POSITIVE,
        1e-300,
        0.1,
        1.0 / 3.0,
        1.0,
        216.0,
        123_456_789.123,
        9_007_199_254_740_993.0,
        1e300,
        f64::MAX,
    ] {
        let mut model = base();
        model["reactions"][1]["rate"]["coefficient"] = json!(value);
        let converted = convert(&model).unwrap();
        let text = converted["rules"][1]["rate"].as_str().unwrap().to_string();
        let parsed: f64 = text.parse().unwrap();
        assert_eq!(parsed.to_bits(), value.to_bits(), "{value:e} became {text}");
        assert!(
            check_model(&from_network(&serde_json::to_vec(&model).unwrap()).unwrap()).valid,
            "{text}"
        );
    }
}

#[test]
fn reserved_function_names_are_renamed_deterministically() {
    let model = json!({
        "schema": "crn-craft.model/0.1", "name": "min",
        "convention": "stochastic_counts", "omega": 1.0,
        "parameters": {"log": 0.5, "exp_": 2.0},
        "species": [{"name": "exp", "initial_count": 4}, {"name": "if", "initial_count": 0}],
        "reactions": [
            {"name": "max", "reactants": {"exp": 1}, "products": {"if": 1},
             "rate": {"kind": "mass_action", "coefficient": "log"}},
            {"name": "sqrt", "reactants": {}, "products": {"exp": 1},
             "rate": {"kind": "hill_repression", "regulator": "if", "vmax": "exp_",
                      "basal": 0.0, "k": "log", "n": 1.0}}]
    });
    let bytes = from_network(&serde_json::to_vec(&model).unwrap()).unwrap();
    assert!(check_model(&bytes).valid);
    let converted: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(converted["name"], "min_");
    // `exp_` is a parameter, so the species `exp` becomes `exp_2`.
    assert_eq!(
        converted["types"],
        json!([
        {"name": "exp_2", "attributes": []}, {"name": "if_", "attributes": []}])
    );
    assert_eq!(converted["parameters"], json!({"exp_": 2.0, "log_": 0.5}));
    assert_eq!(converted["rules"][0]["name"], "max_");
    assert_eq!(converted["rules"][0]["rate"], "log_");
    assert_eq!(converted["rules"][1]["name"], "sqrt_");
    assert_eq!(
        converted["rules"][1]["propensity"],
        "0.0 + exp_ / (1 + (if_ / log_)^1.0)"
    );
    let description = converted["description"].as_str().unwrap();
    for entry in [
        "model `min` is `min_`",
        "species `exp` is `exp_2`",
        "species `if` is `if_`",
        "parameter `log` is `log_`",
        "reaction `max` is `max_`",
        "reaction `sqrt` is `sqrt_`",
    ] {
        assert!(description.contains(entry), "{entry}: {description}");
    }
    assert_eq!(
        from_network(&serde_json::to_vec(&model).unwrap()).unwrap(),
        bytes
    );
    // Without reserved names the description has no rename list.
    let plain = convert(&base()).unwrap();
    assert!(!plain["description"].as_str().unwrap().contains("Renamed"));
}

#[test]
fn omega_does_not_change_the_converted_rules() {
    let mut scaled = base();
    scaled["omega"] = json!(602.5);
    let a = convert(&base()).unwrap();
    let b = convert(&scaled).unwrap();
    assert_eq!(a["rules"], b["rules"]);
    assert_eq!(a["parameters"], b["parameters"]);
    assert!(b["description"].as_str().unwrap().contains("602.5"));
}

#[test]
fn models_beyond_the_rules_limits_are_refused() {
    // The model/1.0 limits (M4a) equal the network contract's: a network
    // model above them keeps its own diagnostic, and one at them converts.
    let mut model = base();
    model["species"] = species(4097);
    model["reactions"] = json!([decay("d", "S0")]);
    let diagnostics = refused(&model);
    assert!(
        diagnostics.iter().all(|d| d.code != REFUSED),
        "{diagnostics:?}"
    );
    let mut model = base();
    model["species"] = species(4096);
    model["reactions"] = json!([decay("d", "S0")]);
    assert!(convert(&model).is_ok());
    let mut model = base();
    model["reactions"] = (0..1025).map(|i| decay(&format!("d{i}"), "A")).collect();
    assert!(convert(&model).is_ok(), "1,025 reactions fit model/1.0");
    let mut model = base();
    for i in 0..1024 {
        model["parameters"][format!("p{i}")] = json!(1.0);
    }
    assert!(convert(&model).is_ok(), "1,025 parameters fit model/1.0");
    // No reaction: a rules model needs one rule.
    let mut model = base();
    model["reactions"] = json!([]);
    assert!(has(&refused(&model), REFUSED, "/reactions"));
}

#[test]
fn higher_order_reactions_convert_with_counts() {
    // Since M4b any stoichiometry converts; tests/higher_order.rs checks
    // the mapping. W4 has order 8.
    for reactants in [
        json!({"A": 3}),
        json!({"A": 2, "B": 1}),
        json!({"A": 1, "B": 1, "C": 1}),
    ] {
        let mut model = base();
        model["reactions"][0]["reactants"] = reactants.clone();
        let converted = convert(&model).unwrap_or_else(|d| panic!("{reactants}: {d:?}"));
        assert_eq!(converted["schema"], "crn-craft.model/1.0");
    }
    let w4 = from_network(&read("validation/ssa-acceleration/models/W4.model.json")).unwrap();
    assert!(crn_craft::rules::check_model(&w4).valid);
}

#[test]
fn products_counts_and_population_limits_are_refused() {
    let mut model = base();
    model["reactions"][0]["products"] = json!({"B": 64, "C": 1});
    assert!(has(&refused(&model), REFUSED, "/reactions/0/products"));
    let mut model = base();
    model["reactions"][0]["products"] = json!({"B": 64});
    assert!(convert(&model).is_ok(), "64 products fit");
    let mut model = base();
    model["species"][0]["initial_count"] = json!(9_007_199_254_740_993_i64);
    assert!(has(&refused(&model), REFUSED, "/species/0/initial_count"));
    let mut model = base();
    model["species"][0]["initial_count"] = json!(9_007_199_254_740_992_i64);
    assert!(convert(&model).is_ok(), "2^53 fits");
    let mut model = base();
    model["species"][2]["population_limit"] = json!(100);
    assert!(has(
        &refused(&model),
        REFUSED,
        "/species/2/population_limit"
    ));
}

#[test]
fn invalid_network_models_keep_their_network_diagnostics() {
    let mut model = base();
    model["reactions"][0]["products"] = json!({"D": 1});
    let diagnostics = refused(&model);
    assert!(
        has(&diagnostics, "CRN-E012", "/reactions/0/products/D"),
        "{diagnostics:?}"
    );
    let diagnostics = from_network(b"{\"schema\": ").unwrap_err();
    assert_eq!(diagnostics[0].code, "CRN-E001");
    // A rules-engine model is not a network model.
    let diagnostics = from_network(&read("models/rules/immigration-death.rules.json")).unwrap_err();
    assert!(has(&diagnostics, "CRN-E032", "/schema"), "{diagnostics:?}");
    let converted = from_network(&read("models/examples/repressilator.model.json")).unwrap();
    assert!(has(
        &from_network(&converted).unwrap_err(),
        "CRN-E032",
        "/schema"
    ));
    // Warnings do not block conversion.
    let mut model = base();
    model["reactions"][0]["rate"]["coefficient"] = json!(0.0);
    assert!(convert(&model).is_ok());
}

#[test]
fn a_converted_model_runs_with_static_direct() {
    let converted = from_network(&read("models/examples/immigration-death.model.json")).unwrap();
    let request = serde_json::to_vec(&json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "m.json", "sha256": crn_craft::sha256(&converted)},
        "time": {"end": 10.0, "report_times": [0.0, 10.0]},
        "ensemble": {"replicates": 400, "root_seed": 3},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 100_000, "max_sampling_steps_per_replicate": 1_000_000,
                   "max_instances": 4, "threads": 1, "max_output_bytes": 1 << 22},
        "method": {"kind": "static_direct"}
    }))
    .unwrap();
    let result = execute(&prepare(&converted, &request).unwrap(), &Control::default()).unwrap();
    assert_eq!(result.termination, Termination::Completed);
    // Immigration at 4 and death at 0.2 per molecule from 0: Poisson with
    // mean 20 (1 - e^-2), about 17.29, at t = 10; the standard error of
    // 400 replicates is about 0.21.
    let mean = result.statistics.unwrap().mean.unwrap()[1][0];
    assert!((mean - 17.29).abs() < 1.5, "{mean}");
}
