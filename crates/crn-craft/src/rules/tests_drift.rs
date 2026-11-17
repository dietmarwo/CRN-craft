//! Regression tests from user feedback on 0.7.1 (a published mtDNA drift
//! study): random labels in pair rules, the `max_count` reduction, and kept
//! final instances counted at run time instead of reserved in advance.
use super::*;
use serde_json::{Value, json};

fn request(model: &[u8], method: &str, end: f64, replicates: u64, keep: u64) -> Value {
    json!({
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": crate::sha256(model)},
        "time": {"end": end, "report_times": [0.0, end]},
        "ensemble": {"replicates": replicates, "root_seed": 11},
        "outputs": {"statistics": ["mean", "variance"],
                    "keep_paths": keep, "keep_instances": keep},
        "budget": {"threads": 2, "max_events_per_replicate": 10_000_000u64,
                   "max_sampling_steps_per_replicate": 100_000_000u64,
                   "max_instances": 100_000u64, "max_output_bytes": 64u64 << 20},
        "method": {"kind": method}
    })
}

fn run(model: &Value, request: &Value) -> RulesResult {
    let model = serde_json::to_vec(model).unwrap();
    let prepared =
        prepare(&model, &serde_json::to_vec(request).unwrap()).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap()
}

/// Two clones of `n` molecules each and pair rules over two identical
/// patterns with the given products.
fn moran(n: u64, rules: &Value) -> Value {
    json!({
        "schema": "crn-craft.model/1.0",
        "name": "moran",
        "parameters": {"g": 0.05},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [
            {"type": "m", "attributes": {"clone": 0}, "count": n},
            {"type": "m", "attributes": {"clone": 1}, "count": n}
        ],
        "observables": {
            "N": {"type": "m"},
            "c1": {"type": "m", "where": "x.clone == 1"}
        },
        "rules": rules
    })
}

fn pair() -> Value {
    json!([{"label": "x", "type": "m"}, {"label": "y", "type": "m"}])
}

fn compiled(model: &Value) -> model::Compiled {
    let model: Model = serde_json::from_value(model.clone()).unwrap();
    model::compile(&model).unwrap()
}

/// A neutral Moran step written as one pair rule ("x replaces y") keeps the
/// mean of either clone. Version 2 of the sampler gave `x` to the instance
/// earlier in its internal order, so clone 0 always replaced clone 1.
#[test]
fn one_pair_rule_is_neutral() {
    let model = moran(
        20,
        &json!([{"name": "x_replaces_y", "reactants": pair(), "rate": "g",
                 "products": [{"type": "m", "copy": "x"}, {"type": "m", "copy": "x"}]}]),
    );
    let result = run(
        &model,
        &request(&serde_json::to_vec(&model).unwrap(), "tree", 5.0, 4_000, 0),
    );
    assert_eq!(result.sampler, SAMPLER_RANDOM_LABELS);
    let stats = result.statistics.unwrap();
    let (mean, variance) = (stats.mean.unwrap(), stats.sample_variance.unwrap());
    // Observables N, c1 at the end time.
    assert_eq!(mean[1][0], 40.0);
    let (c1, error) = (mean[1][1], (variance[1][1] / 4_000.0).sqrt());
    assert!(error > 0.1, "the clones drift: standard error {error}");
    assert!(
        (c1 - 20.0).abs() < 5.0 * error,
        "clone 1 mean {c1} +- {error}"
    );
}

/// The two orientations as separate rules give the same law.
#[test]
fn both_orientations_are_neutral_too() {
    let model = moran(
        20,
        &json!([
            {"name": "x_replaces_y", "reactants": pair(), "rate": "g / 2",
             "products": [{"type": "m", "copy": "x"}, {"type": "m", "copy": "x"}]},
            {"name": "y_replaces_x", "reactants": pair(), "rate": "g / 2",
             "products": [{"type": "m", "copy": "y"}, {"type": "m", "copy": "y"}]}
        ]),
    );
    let result = run(
        &model,
        &request(&serde_json::to_vec(&model).unwrap(), "tree", 5.0, 4_000, 0),
    );
    let stats = result.statistics.unwrap();
    let c1 = stats.mean.unwrap()[1][1];
    let error = (stats.sample_variance.unwrap()[1][1] / 4_000.0).sqrt();
    assert!(
        (c1 - 20.0).abs() < 5.0 * error,
        "clone 1 mean {c1} +- {error}"
    );
}

/// Only rules whose products can depend on the labels draw the extra
/// uniform; symmetric products keep the sampler identity and stream of
/// version 2.
#[test]
fn label_dependence_is_detected_conservatively() {
    let cluster = |products: Value, reactants: Value| {
        json!({
            "schema": "crn-craft.model/1.0",
            "name": "clusters",
            "parameters": {"k": 1.0},
            "types": [{"name": "c", "attributes": ["size", "id"]}],
            "initial": [
                {"type": "c", "attributes": {"size": 1, "id": 1}, "count": 3},
                {"type": "c", "attributes": {"size": 2, "id": 2}, "count": 3}
            ],
            "observables": {"n": {"type": "c"}},
            "rules": [{"name": "merge", "reactants": reactants, "rate": "k",
                       "products": products}]
        })
    };
    let two = json!([{"label": "x", "type": "c"}, {"label": "y", "type": "c"}]);
    let cases = [
        // Merging: x.size + y.size equals y.size + x.size.
        (
            json!([{"type": "c", "attributes": {"size": "x.size + y.size", "id": "max(x.id, y.id)"}}]),
            two.clone(),
            false,
        ),
        // Both kept: the product list maps onto itself.
        (
            json!([{"type": "c", "copy": "x"}, {"type": "c", "copy": "y"}]),
            two.clone(),
            false,
        ),
        // Annihilation: no product reads a reactant.
        (json!([]), two.clone(), false),
        // x survives, y is lost.
        (json!([{"type": "c", "copy": "x"}]), two.clone(), true),
        (
            json!([{"type": "c", "attributes": {"size": "x.size - y.size", "id": "x.id"}}]),
            two.clone(),
            true,
        ),
        // A random product attribute: always drawn.
        (
            json!([{"type": "c", "attributes": {"size": "x.size + y.size", "id": "new_id()"}}]),
            two.clone(),
            true,
        ),
    ];
    for (i, (products, reactants, dependent)) in cases.into_iter().enumerate() {
        let model = cluster(products, reactants);
        assert_eq!(
            compiled(&model).rules[0].random_labels,
            dependent,
            "case {i}"
        );
        let result = run(
            &model,
            &request(&serde_json::to_vec(&model).unwrap(), "tree", 0.5, 2, 0),
        );
        let expected = if dependent {
            SAMPLER_RANDOM_LABELS
        } else {
            SAMPLER
        };
        assert_eq!(result.sampler, expected, "case {i}");
    }
}

fn clones(reduce_model: &str) -> Value {
    json!({
        "schema": reduce_model,
        "name": "clones",
        "parameters": {"d": 0.0},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [
            {"type": "m", "attributes": {"clone": 0}, "count": 5},
            {"type": "m", "attributes": {"clone": 1}, "count": 3},
            {"type": "m", "attributes": {"clone": 2}, "count": 7}
        ],
        "observables": {
            "largest": {"type": "m", "reduce": "max_count"},
            "largest_mutant": {"type": "m", "where": "x.clone > 0", "reduce": "max_count"},
            "weighted": {"type": "m", "where": "x.clone < 2", "weight": "2", "reduce": "max_count"},
            "none": {"type": "m", "where": "x.clone > 9", "reduce": "max_count"}
        },
        "rules": [{"name": "decay", "reactants": [{"label": "x", "type": "m"}],
                   "rate": "d", "products": []}]
    })
}

/// `max_count` is the largest count times weight over the matching present
/// instances, with the same value on the tree and the static methods.
#[test]
fn max_count_on_tree_and_static_methods() {
    let model = clones("crn-craft.model/1.0");
    let bytes = serde_json::to_vec(&model).unwrap();
    for method in ["tree", "static_direct", "next_reaction"] {
        let result = run(&model, &request(&bytes, method, 1.0, 1, 1));
        // Observables in name order: largest, largest_mutant, none, weighted.
        assert_eq!(
            result.paths[0].values[0],
            vec![7.0, 7.0, 0.0, 10.0],
            "{method}"
        );
    }
}

/// Counts change: the largest clone follows births and deaths on the tree.
#[test]
fn max_count_follows_the_counts() {
    let model = json!({
        "schema": "crn-craft.model/1.0",
        "name": "births",
        "parameters": {"b": 1.0},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [
            {"type": "m", "attributes": {"clone": 0}, "count": 2},
            {"type": "m", "attributes": {"clone": 1}, "count": 1}
        ],
        "observables": {
            "largest": {"type": "m", "reduce": "max_count"},
            "c0": {"type": "m", "where": "x.clone == 0"},
            "c1": {"type": "m", "where": "x.clone == 1"}
        },
        "rules": [{"name": "birth", "reactants": [{"label": "x", "type": "m"}], "rate": "b",
                   "products": [{"type": "m", "copy": "x"}, {"type": "m", "copy": "x"}]}]
    });
    let result = run(
        &model,
        &request(&serde_json::to_vec(&model).unwrap(), "tree", 2.0, 50, 50),
    );
    for path in &result.paths {
        for row in &path.values {
            // Observables in name order: c0, c1, largest.
            assert_eq!(row[2], row[0].max(row[1]));
        }
    }
}

/// Requests before `crn-craft.run/1.0` replay the label order of sampler
/// version 2, so stored requests stay byte-identical.
#[test]
fn legacy_requests_keep_version_2_labels() {
    let mut model = moran(
        20,
        &json!([{"name": "x_replaces_y", "reactants": pair(), "rate": "g",
                 "products": [{"type": "m", "copy": "x"}, {"type": "m", "copy": "x"}]}]),
    );
    model["schema"] = json!("crn-craft.rules-model/0.3");
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut legacy = request(&bytes, "tree", 5.0, 200, 0);
    legacy["schema"] = json!("crn-craft.rules-run/0.2");
    legacy["outputs"] = json!({"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0});
    let result = run(&model, &legacy);
    assert_eq!(result.sampler, SAMPLER);
    // Version 2 lets the earlier instance (clone 0) take over.
    assert!(result.statistics.unwrap().mean.unwrap()[1][1] < 2.0);
}

/// Before `crn-craft.model/1.0`, `max_count` is the strict parser's
/// unknown-variant error.
#[test]
fn max_count_needs_model_1_0() {
    let model = clones("crn-craft.rules-model/0.3");
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut legacy = request(&bytes, "tree", 1.0, 1, 0);
    legacy["schema"] = json!("crn-craft.rules-run/0.2");
    legacy["outputs"] = json!({"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0});
    let Err(errors) = prepare(&bytes, &serde_json::to_vec(&legacy).unwrap()) else {
        panic!("max_count was accepted before model 1.0");
    };
    assert_eq!(errors[0].code, "CRN-E003");
    assert_eq!(errors[0].path, "/observables/largest/reduce");
}

/// A model that creates a new instance per event, for kept-instance tests.
fn founders() -> Value {
    json!({
        "schema": "crn-craft.model/1.0",
        "name": "founders",
        "parameters": {"k": 1.0},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [{"type": "m", "attributes": {"clone": 0}, "count": 1}],
        "observables": {"n": {"type": "m"}},
        "rules": [{"name": "found", "reactants": [], "rate": "k",
                   "products": [{"type": "m", "attributes": {"clone": "new_id()"}}]}]
    })
}

/// Kept instances are no longer reserved as keep_instances × max_instances:
/// 3000 kept sets with the default max_instances of 100,000 are admitted
/// (16.8 GB under the old reservation).
#[test]
fn kept_instances_are_not_reserved_in_advance() {
    let model = founders();
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request(&bytes, "tree", 3.0, 3_000, 3_000);
    req["outputs"]["keep_paths"] = json!(0);
    req["budget"]["max_output_bytes"] = json!(1u64 << 30);
    let result = run(&model, &req);
    assert_eq!(result.termination, Termination::Completed);
    assert_eq!(result.final_instances.len(), 3_000);
}

/// Kept instances past max_output_bytes are dropped during the run and the
/// result ends as output_limit with CRN-R045, as a too-large result did.
#[test]
fn kept_instances_past_the_output_budget_end_as_output_limit() {
    let model = founders();
    let bytes = serde_json::to_vec(&model).unwrap();
    // About 2000 instances x 56 bytes per replicate, 4 replicates: 450 kB.
    let mut req = request(&bytes, "tree", 2_000.0, 4, 4);
    req["outputs"]["keep_paths"] = json!(0);
    req["budget"]["max_output_bytes"] = json!(100_000u64);
    let result = run(&model, &req);
    assert_eq!(result.termination, Termination::OutputLimit);
    assert!(result.final_instances.is_empty());
    assert!(result.diagnostics.iter().any(|d| d.code == "CRN-R045"));
}

/// The remaining reservation error itemizes what it reserved.
#[test]
fn reservation_error_is_itemized() {
    let model = founders();
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request(&bytes, "tree", 1.0, 1_000_000, 1_000_000);
    req["outputs"]["keep_instances"] = json!(0);
    req["time"]["report_times"] = json!((0..400).map(f64::from).collect::<Vec<_>>());
    req["time"]["end"] = json!(400.0);
    req["budget"]["max_output_bytes"] = json!(1u64 << 30);
    let Err(errors) = prepare(&bytes, &serde_json::to_vec(&req).unwrap()) else {
        panic!("the reservation was admitted");
    };
    let Some(error) = errors.iter().find(|d| d.path == "/budget") else {
        panic!("{errors:?}");
    };
    let message = &error.message;
    assert!(message.contains("kept paths"), "{message}");
    assert!(message.contains("400 report times"), "{message}");
}

/// Review finding: a result whose instance estimate exceeds the budget but
/// whose JSON fits is returned complete (the estimate only guards memory).
#[test]
fn kept_instances_that_fit_are_kept() {
    let model = founders();
    let bytes = serde_json::to_vec(&model).unwrap();
    // About 20,000 instances in one kept set: estimate ~1.1 MB, JSON ~0.93 MB.
    let mut req = request(&bytes, "tree", 20_000.0, 1, 1);
    req["outputs"]["keep_paths"] = json!(0);
    req["budget"]["max_output_bytes"] = json!(1_000_000u64);
    let result = run(&model, &req);
    assert_eq!(result.termination, Termination::Completed);
    let instances = result.final_instances[0].instances.len();
    assert!(
        instances > 1_000_000 / 56,
        "the estimate exceeds the budget: {instances}"
    );
}

/// Requests before crn-craft.run/1.0 keep the 0.7 reservation of kept
/// instances, so they are admitted and refused as before.
#[test]
fn legacy_requests_keep_the_instance_reservation() {
    let mut model = founders();
    model["schema"] = json!("crn-craft.rules-model/0.3");
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut legacy = request(&bytes, "tree", 1.0, 1_000, 1_000);
    legacy["schema"] = json!("crn-craft.rules-run/0.2");
    legacy["outputs"] = json!({"statistics": ["mean"], "keep_paths": 0, "keep_instances": 1_000});
    legacy["budget"]["max_output_bytes"] = json!(1u64 << 30);
    let Err(errors) = prepare(&bytes, &serde_json::to_vec(&legacy).unwrap()) else {
        panic!("the 0.7 reservation was not applied");
    };
    let error = errors.iter().find(|d| d.path == "/budget").unwrap();
    assert_eq!(error.code, "CRN-R020");
    assert!(
        error.message.contains("kept instance sets"),
        "{}",
        error.message
    );
    assert!(error.hint.contains("max_instances"), "{}", error.hint);
}

/// The hint names the largest part of the reservation.
#[test]
fn reservation_hint_names_the_largest_part() {
    let mut model = founders();
    for i in 0..100 {
        model["observables"][format!("n{i}")] = json!({"type": "m"});
    }
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut req = request(&bytes, "tree", 1.0, 64, 0);
    req["budget"]["threads"] = json!(64);
    req["time"]["report_times"] = json!(
        (0..100_000)
            .map(|i| f64::from(i) * 1e-5)
            .collect::<Vec<_>>()
    );
    let Err(errors) = prepare(&bytes, &serde_json::to_vec(&req).unwrap()) else {
        panic!("the reservation was admitted");
    };
    let error = errors.iter().find(|d| d.path == "/budget").unwrap();
    assert!(
        error.message.contains("in-flight paths"),
        "{}",
        error.message
    );
    assert!(error.hint.contains("fewer threads"), "{}", error.hint);
}

/// A pair over a type without attributes never spans two instances, so it
/// keeps sampler version 2 even with random product attributes.
#[test]
fn attribute_free_pairs_never_draw_labels() {
    let model = json!({
        "schema": "crn-craft.model/1.0",
        "name": "dimers",
        "parameters": {"k": 1.0},
        "types": [{"name": "A"}, {"name": "D", "attributes": ["u"]}],
        "initial": [{"type": "A", "count": 10}],
        "observables": {"A": {"type": "A"}, "D": {"type": "D"}},
        "rules": [{"name": "dimerize", "reactants": [{"label": "x", "type": "A", "count": 2}],
                   "rate": "k", "products": [{"type": "D", "attributes": {"u": "uniform(0, 1)"}}]}]
    });
    assert!(!compiled(&model).rules[0].random_labels);
    let result = run(
        &model,
        &request(&serde_json::to_vec(&model).unwrap(), "tree", 1.0, 2, 0),
    );
    assert_eq!(result.sampler, SAMPLER);
}

/// With both a pattern `count` and a `max_count` observable, an old model
/// gets the error 0.7.1 gave: the reduction it could not parse.
#[test]
fn legacy_errors_keep_the_0_7_order() {
    let mut model = clones("crn-craft.rules-model/0.3");
    model["rules"][0]["reactants"][0]["count"] = json!(1);
    let bytes = serde_json::to_vec(&model).unwrap();
    let mut legacy = request(&bytes, "tree", 1.0, 1, 0);
    legacy["schema"] = json!("crn-craft.rules-run/0.2");
    legacy["outputs"] = json!({"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0});
    let Err(errors) = prepare(&bytes, &serde_json::to_vec(&legacy).unwrap()) else {
        panic!("accepted");
    };
    assert_eq!(errors[0].path, "/observables/largest/reduce");
}
