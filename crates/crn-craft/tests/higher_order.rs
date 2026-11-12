//! M4b: higher-order reactions in `crn-craft.model/1.0`.
//!
//! A reactant pattern may take `count` molecules (1 to 2,147,483,647) and a
//! rule may have any number of patterns. Every pattern with `count > 1`, and
//! every pattern of a rule with more than two patterns, is over an
//! attribute-free type without `where` or `factor`. Patterns over one type
//! add their counts. A `rate` gives rate × Π C(nᵢ, kᵢ); a `propensity`
//! applies while every type has its count. Rules of at most two molecules
//! keep the shapes and arithmetic of 0.6.0; `{"count": 2}` is the pair of two
//! identical patterns.
use crn_craft::diagnostic::Diagnostic;
use crn_craft::rules::convert::from_network;
use crn_craft::rules::{
    Control, RulesResult, check_model, execute, generate_rules_kernel, generate_rules_loop, prepare,
};
use serde_json::{Value, json};
use std::path::PathBuf;

const UNIFIED: &str = "crn-craft.model/1.0";

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

/// A model with types A (attribute-free), B (attribute-free), C (attribute
/// `size`) and one rule with the given reactants, rate and products.
fn model(schema: &str, reactants: Value, products: Value) -> Value {
    json!({
        "schema": schema,
        "name": "order",
        "parameters": {"k": 0.002},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C", "attributes": ["size"]}],
        "initial": [{"type": "A", "count": 40}, {"type": "B", "count": 20},
                    {"type": "C", "attributes": {"size": 1}, "count": 5}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}, "C": {"type": "C"}},
        "rules": [{"name": "react", "reactants": reactants, "rate": "k", "products": products}]
    })
}

fn request(model: &[u8], method: Value, replicates: u64, end: f64) -> Vec<u8> {
    let unified = serde_json::from_slice::<Value>(model).unwrap()["schema"] == UNIFIED;
    bytes(&json!({
        "schema": if unified { "crn-craft.run/1.0" } else { "crn-craft.rules-run/0.6" },
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": end, "report_times": [0.0, end / 4.0, end / 2.0, end]},
        "ensemble": {"replicates": replicates, "root_seed": 5},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": replicates,
                    "keep_instances": 1},
        "budget": {"max_events_per_replicate": 1_000_000, "max_sampling_steps_per_replicate": 100_000_000,
                   "max_instances": 64, "threads": 2, "max_output_bytes": 1u64 << 28},
        "method": method,
        "event_log": {"max_entries": 4}
    }))
}

fn run(model: &Value, method: Value, replicates: u64, end: f64) -> RulesResult {
    let model = bytes(model);
    let prepared = prepare(&model, &request(&model, method, replicates, end))
        .unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap()
}

fn refused(model: &Value, method: Value) -> Vec<Diagnostic> {
    let model = bytes(model);
    match prepare(&model, &request(&model, method, 2, 1.0)) {
        Ok(_) => panic!("prepared"),
        Err(d) => d,
    }
}

fn has(diagnostics: &[Diagnostic], code: &str, path: &str) -> bool {
    diagnostics.iter().any(|d| d.code == code && d.path == path)
}

fn check(model: &Value) -> Vec<Diagnostic> {
    check_model(&bytes(model)).diagnostics
}

fn pattern(label: &str, ty: &str, count: Option<u64>) -> Value {
    match count {
        Some(n) => json!({"label": label, "type": ty, "count": n}),
        None => json!({"label": label, "type": ty}),
    }
}

/// The result without the fields that name the model bytes.
fn numbers(result: &RulesResult) -> Value {
    let mut value = serde_json::to_value(result).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("model_sha256").unwrap();
    object.remove("run_sha256").unwrap();
    value
}

#[test]
fn count_is_a_model_1_0_field() {
    // Before 1.0 the field is unknown, with the 0.6.0 diagnostic.
    for version in ["0.1", "0.2", "0.3", "0.4"] {
        let schema = format!("crn-craft.rules-model/{version}");
        let old = model(&schema, json!([pattern("a", "A", Some(3))]), json!([]));
        let diagnostics = check(&old);
        assert_eq!(diagnostics.len(), 1, "{diagnostics:?}");
        let d = &diagnostics[0];
        assert_eq!(
            (d.code.as_str(), d.path.as_str(), d.message.as_str()),
            (
                "CRN-E003",
                "/rules/0/reactants/0/count",
                "unknown field `count`, expected one of `label`, `type`, `where`, `factor`"
            )
        );
        assert_eq!(
            d.hint,
            "Compare this field with the schema emitted by describe."
        );
        let request = request(&bytes(&old), json!({"kind": "static_direct"}), 2, 1.0);
        let prepared = prepare(&bytes(&old), &request).err().unwrap();
        assert_eq!(
            serde_json::to_value(&prepared).unwrap(),
            serde_json::to_value(&diagnostics).unwrap()
        );
        // More than two patterns keep their 0.x diagnostic.
        let three = model(
            &schema,
            json!([
                pattern("a", "A", None),
                pattern("b", "A", None),
                pattern("c", "A", None)
            ]),
            json!([]),
        );
        let diagnostics = check(&three);
        assert!(
            has(&diagnostics, "CRN-R013", "/rules/0/reactants"),
            "{diagnostics:?}"
        );
        assert_eq!(diagnostics[0].message, "more than two reactant patterns");
    }
    // In 1.0: 1 to 2,147,483,647, an integer, not null.
    for count in [1u64, 2, 3, 2_147_483_647] {
        let m = model(UNIFIED, json!([pattern("a", "A", Some(count))]), json!([]));
        assert!(check(&m).is_empty(), "count {count}");
    }
    for count in [0u64, 2_147_483_648, u64::MAX] {
        let m = model(UNIFIED, json!([pattern("a", "A", Some(count))]), json!([]));
        assert!(
            has(&check(&m), "CRN-R013", "/rules/0/reactants/0/count"),
            "count {count}"
        );
    }
    for count in [json!(-1), json!(2.5), json!(null), json!("3")] {
        let m = model(
            UNIFIED,
            json!([{"label": "a", "type": "A", "count": count}]),
            json!([]),
        );
        let diagnostics = check(&m);
        assert_eq!(diagnostics.len(), 1, "{count}");
        assert_eq!(diagnostics[0].code, "CRN-E003", "{count}");
    }
}

#[test]
fn higher_order_patterns_must_be_attribute_free_without_where_or_factor() {
    let cases = [
        // count > 1 over an attributed type
        (
            json!([pattern("c", "C", Some(2))]),
            "/rules/0/reactants/0/type",
        ),
        // count > 1 with where or factor
        (
            json!([{"label": "a", "type": "A", "count": 2, "where": "k > 0"}]),
            "/rules/0/reactants/0/where",
        ),
        (
            json!([{"label": "a", "type": "A", "count": 3, "factor": "2"}]),
            "/rules/0/reactants/0/factor",
        ),
        // more than two patterns: every one plain and attribute-free
        (
            json!([
                pattern("a", "A", None),
                pattern("b", "B", None),
                pattern("c", "C", None)
            ]),
            "/rules/0/reactants/2/type",
        ),
        (
            json!([pattern("a", "A", None), {"label": "b", "type": "B", "factor": "2"},
                   pattern("x", "A", None)]),
            "/rules/0/reactants/1/factor",
        ),
        // different patterns over one type
        (
            json!([pattern("a", "A", Some(2)), {"label": "b", "type": "A", "where": "k > 0"}]),
            "/rules/0/reactants",
        ),
    ];
    for (reactants, path) in cases {
        let m = model(UNIFIED, reactants.clone(), json!([]));
        let diagnostics = check(&m);
        assert!(
            has(&diagnostics, "CRN-R013", path),
            "{reactants}: {diagnostics:?}"
        );
    }
    // Labels are unique over every pattern.
    let m = model(
        UNIFIED,
        json!([
            pattern("a", "A", None),
            pattern("b", "B", None),
            pattern("a", "B", None)
        ]),
        json!([]),
    );
    assert!(has(&check(&m), "CRN-R010", "/rules/0/reactants/2/label"));
    // A count-1 pattern beside a count > 1 pattern may have attributes,
    // where and factor, and its attributes reach the products.
    let m = model(
        UNIFIED,
        json!([pattern("a", "A", Some(2)),
               {"label": "c", "type": "C", "where": "c.size > 0", "factor": "c.size"}]),
        json!([{"type": "C", "copy": "c"}, {"type": "C", "attributes": {"size": "c.size"}}]),
    );
    assert!(check(&m).is_empty(), "{:?}", check(&m));
    let result = run(&m, json!({"kind": "static_direct"}), 4, 5.0);
    assert_eq!(result.completed_replicates, 4);
    // 2A + C -> 2C: A falls by two per event and C rises by one.
    for path in &result.paths {
        let last = path.values.last().unwrap();
        assert_eq!(40.0 - last[0], 2.0 * (last[2] - 5.0), "{last:?}");
        assert!(last[2] > 5.0);
    }
}

#[test]
fn count_two_is_the_pair_of_two_identical_patterns() {
    // The same results, event logs included, on the tree and on Static
    // Direct: {"count": 2} keeps the pair shape and its arithmetic.
    let products = json!([{"type": "B"}]);
    let pair = model(
        UNIFIED,
        json!([pattern("a", "A", None), pattern("b", "A", None)]),
        products.clone(),
    );
    let counted = model(UNIFIED, json!([pattern("a", "A", Some(2))]), products);
    for method in [json!({"kind": "tree"}), json!({"kind": "static_direct"})] {
        let a = run(&pair, method.clone(), 8, 20.0);
        let b = run(&counted, method, 8, 20.0);
        assert_eq!(numbers(&a), numbers(&b));
        assert!(a.replicates.iter().all(|r| r.events > 0));
    }
    // Three identical patterns are the same rule as one pattern with count 3.
    let three = model(
        UNIFIED,
        json!([
            pattern("a", "A", None),
            pattern("b", "A", None),
            pattern("c", "A", None)
        ]),
        json!([]),
    );
    let counted = model(UNIFIED, json!([pattern("a", "A", Some(3))]), json!([]));
    let a = run(&three, json!({"kind": "static_direct"}), 8, 20.0);
    let b = run(&counted, json!({"kind": "static_direct"}), 8, 20.0);
    let (mut a, mut b) = (numbers(&a), numbers(&b));
    // Event logs list one entry per pattern.
    for value in [&mut a, &mut b] {
        value.as_object_mut().unwrap().remove("observations");
    }
    assert_eq!(a, b);
}

#[test]
fn every_static_method_runs_higher_order_rules_and_tree_and_mtdna_refuse() {
    // 3A -> 0 from A = 40 keeps A = 1 (mod 3); 2A + B -> C keeps A - 2B.
    let decay = model(UNIFIED, json!([pattern("a", "A", Some(3))]), json!([]));
    let bind = model(
        UNIFIED,
        json!([pattern("a", "A", Some(2)), pattern("b", "B", None)]),
        json!([{"type": "C", "attributes": {"size": 1}}]),
    );
    for method in [
        json!({"kind": "static_direct"}),
        json!({"kind": "static_direct", "waiting_time": "inverse"}),
        json!({"kind": "next_reaction"}),
        json!({"kind": "tau", "epsilon": 0.03}),
        json!({"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}),
    ] {
        for (m, invariant) in [(&decay, 0), (&bind, 1)] {
            let result = run(m, method.clone(), 16, 10.0);
            assert_eq!(result.completed_replicates, 16, "{method}");
            let a = result.observables.iter().position(|o| o == "A").unwrap();
            let b = result.observables.iter().position(|o| o == "B").unwrap();
            let mut fired = 0;
            for path in &result.paths {
                for row in &path.values {
                    if invariant == 0 {
                        assert_eq!(row[a] as u64 % 3, 1, "{method}: {row:?}");
                    } else {
                        assert_eq!(row[a] - 2.0 * row[b], 0.0, "{method}: {row:?}");
                    }
                }
                fired += u64::from(path.values.last().unwrap()[a] < 40.0);
            }
            assert!(fired > 0, "{method}: no event");
        }
    }
    for (method, path) in [
        (json!({"kind": "tree"}), "/method"),
        (json!({"kind": "mtdna_direct"}), "/method"),
        (json!({"kind": "binomial", "step": 0.1}), "/method"),
    ] {
        let d = refused(&decay, method.clone());
        assert!(has(&d, "CRN-R020", path), "{method}: {d:?}");
    }
    // The default method is the tree.
    let m = bytes(&decay);
    let mut no_method: Value =
        serde_json::from_slice(&request(&m, json!({"kind": "tree"}), 2, 1.0)).unwrap();
    no_method.as_object_mut().unwrap().remove("method");
    let d = prepare(&m, &bytes(&no_method)).err().unwrap();
    assert!(has(&d, "CRN-R020", "/method"), "{d:?}");
    // Generated rates kernels refuse; generated event loops emit the rule.
    let prepared = prepare(&m, &request(&m, json!({"kind": "static_direct"}), 2, 1.0)).unwrap();
    let d = generate_rules_kernel(&prepared).unwrap_err();
    assert_eq!((d.code.as_str(), d.path.as_str()), ("CRN-R020", "/method"));
    let source = generate_rules_loop(&prepared).unwrap();
    assert!(source.contains("crn_craft::rules::combination_rate("));
}

#[test]
fn propensities_count_combinations_and_whole_propensities_need_every_count() {
    // One replicate up to its first event: the first wait is exponential
    // with the total propensity, so compare event counts over short horizons
    // through the ensemble mean instead, against the exact first-order
    // expansion at t -> 0: E[events] ~ a0 * t.
    let decay = model(UNIFIED, json!([pattern("a", "A", Some(3))]), json!([]));
    let a0 = 0.002 * 9880.0; // c * C(40, 3)
    let t = 1e-3;
    let result = run(&decay, json!({"kind": "static_direct"}), 20_000, t);
    let events: u64 = result.replicates.iter().map(|r| r.events).sum();
    let mean = events as f64 / 20_000.0;
    let expected = a0 * t;
    let se = (expected / 20_000.0).sqrt();
    assert!((mean - expected).abs() < 5.0 * se, "{mean} vs {expected}");
    // A whole-rule propensity fires while A >= 3 and B >= 1, at its value.
    let mut whole = model(
        UNIFIED,
        json!([pattern("a", "A", Some(3)), pattern("b", "B", None)]),
        json!([]),
    );
    whole["rules"][0].as_object_mut().unwrap().remove("rate");
    whole["rules"][0]["propensity"] = json!("5");
    whole["initial"] = json!([{"type": "A", "count": 7}, {"type": "B", "count": 9}]);
    let result = run(&whole, json!({"kind": "static_direct"}), 64, 1000.0);
    for path in &result.paths {
        let last = path.values.last().unwrap();
        // Two events: A 7 -> 4 -> 1, B 9 -> 8 -> 7; then A < 3 stops it.
        assert_eq!((last[0], last[1]), (1.0, 7.0));
    }
}

#[test]
fn overflowing_combination_counts_fail_the_replicate_with_a_coded_diagnostic() {
    // Ported from the network engine: binomial(n, k) may overflow while
    // rate * binomial(n, k) is representable, and a true overflow is coded.
    let mut m = model(UNIFIED, json!([pattern("a", "A", Some(1000))]), json!([]));
    m["parameters"]["k"] = json!(1e-300);
    m["initial"] = json!([{"type": "A", "count": 2000}]);
    let result = run(&m, json!({"kind": "static_direct"}), 2, 1e-300);
    assert_eq!(result.completed_replicates, 2, "{:?}", result.replicates[0]);
    let mut m = model(
        UNIFIED,
        json!([pattern("a", "A", Some(2_147_483_647))]),
        json!([]),
    );
    m["parameters"]["k"] = json!(5e-324);
    m["initial"] = json!([{"type": "A", "count": 9_007_199_254_740_992u64}]);
    let result = run(&m, json!({"kind": "static_direct"}), 2, 1.0);
    for record in &result.replicates {
        let d = record.diagnostic.as_ref().unwrap();
        assert_eq!(
            (d.code.as_str(), d.path.as_str()),
            ("CRN-R030", "/rules/0/rate")
        );
    }
    // n < count: propensity 0, the rule never fires and its rate is never read.
    let mut m = model(UNIFIED, json!([pattern("a", "A", Some(41))]), json!([]));
    m["rules"][0]["rate"] = json!("log(0 - k)");
    let result = run(&m, json!({"kind": "static_direct"}), 2, 1.0);
    assert!(
        result
            .replicates
            .iter()
            .all(|r| r.events == 0 && r.diagnostic.is_none())
    );
}

#[test]
fn the_converter_maps_any_stoichiometry_with_count() {
    let network = |reactants: Value| {
        json!({
            "schema": "crn-craft.model/0.1", "name": "law", "convention": "stochastic_counts",
            "omega": 1.0, "parameters": {"c": 0.5},
            "species": [{"name": "A", "initial_count": 10}, {"name": "B", "initial_count": 10},
                        {"name": "C", "initial_count": 10}],
            "reactions": [{"name": "r", "reactants": reactants, "products": {},
                           "rate": {"kind": "mass_action", "coefficient": "c"}}]
        })
    };
    for (reactants, patterns) in [
        // Up to two molecules: unit patterns x and y, as in 0.6.0.
        (
            json!({"A": 2}),
            json!([{"label": "x", "type": "A"}, {"label": "y", "type": "A"}]),
        ),
        (
            json!({"A": 1, "B": 1}),
            json!([{"label": "x", "type": "A"}, {"label": "y", "type": "B"}]),
        ),
        // Higher order: one pattern per species, with its count.
        (
            json!({"A": 3}),
            json!([{"label": "x", "type": "A", "count": 3}]),
        ),
        (
            json!({"A": 2, "B": 1}),
            json!([{"label": "x", "type": "A", "count": 2}, {"label": "y", "type": "B"}]),
        ),
        (
            json!({"A": 1, "B": 1, "C": 1}),
            json!([{"label": "x1", "type": "A"}, {"label": "x2", "type": "B"},
                   {"label": "x3", "type": "C"}]),
        ),
        (
            json!({"A": 2147483647}),
            json!([{"label": "x", "type": "A", "count": 2147483647}]),
        ),
    ] {
        let converted = from_network(&bytes(&network(reactants.clone()))).unwrap();
        assert!(check_model(&converted).valid, "{reactants}");
        let value: Value = serde_json::from_slice(&converted).unwrap();
        assert_eq!(value["rules"][0]["reactants"], patterns, "{reactants}");
    }
    let w4 =
        std::fs::read(root().join("validation/ssa-acceleration/models/W4.model.json")).unwrap();
    let converted = from_network(&w4).unwrap_or_else(|d| panic!("{d:?}"));
    assert!(check_model(&converted).valid);
}

#[test]
fn the_rules_model_schema_has_no_count_and_model_1_0_has_it() {
    let old = crn_craft::rules::schema("rules-model").unwrap();
    assert!(
        old["definitions"]["Pattern"]["properties"]
            .get("count")
            .is_none()
    );
    let new = crn_craft::rules::schema("model-1.0").unwrap();
    let count = &new["definitions"]["Pattern"]["properties"]["count"];
    assert_eq!(count["minimum"], 1.0);
    assert_eq!(count["maximum"], 2_147_483_647.0);
}
