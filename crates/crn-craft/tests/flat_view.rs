//! The flat view of a closed rules model (M5a): structure and propensities
//! of converted network models against the network engine of release 0.6.0
//! (`law_oracle`), labels of attributed models, and the refusal of models
//! without a flat view.
use crn_craft::rules::{FlatView, convert};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};

#[path = "../src/rules/law_oracle.rs"]
mod law_oracle;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// Every network model of the repository with an FSP or law reference:
/// the examples, the workloads W1 to W7 and the two higher-order models.
fn network_models() -> Vec<PathBuf> {
    let mut out = Vec::new();
    for folder in [
        "models/legacy",
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

/// Deterministic pseudo-random counts.
fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *seed >> 33
}

/// The difference allowed between the converted model's propensity and the
/// network engine's, as measured and declared in M2 (rate laws) and M4b
/// (higher-order mass action): bit-identical for mass action of order 0, 1
/// and 3 or more and for Michaelis-Menten; relative 4 epsilon for mass
/// action over pairs; relative 1e-14 max(1, n) for Hill laws, with an
/// absolute floor of vmax 2^-1020 where the induced term underflows.
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

#[test]
fn converted_models_have_the_network_structure_and_propensities() {
    let mut models = 0;
    let mut compared = 0;
    let mut identical = 0;
    for path in network_models() {
        let bytes = std::fs::read(&path).unwrap();
        let source: Value = serde_json::from_slice(&bytes).unwrap();
        let network = law_oracle::Network::new(&source);
        let view = FlatView::from_bytes(&convert::from_network(&bytes).unwrap())
            .unwrap_or_else(|d| panic!("{}: {d:?}", name(&path)));
        let label = name(&path);
        // One species per network species in the same order, one channel per
        // reaction.
        assert_eq!(view.species(), network.species, "{label}");
        let initial: Vec<u64> = network.initial.iter().map(|&n| n as u64).collect();
        assert_eq!(view.initial(), initial, "{label}");
        assert_eq!(view.channels(), network.reactions.len(), "{label}");
        assert_eq!(view.rules(), network.reactions.len(), "{label}");
        for (c, reaction) in network.reactions.iter().enumerate() {
            assert_eq!(view.rule(c), c, "{label}");
            assert_eq!(view.rule_name(c), reaction.name, "{label}");
            let counted = |list: &[(usize, u64)]| -> Vec<(usize, i64)> {
                list.iter().map(|&(s, k)| (s, k as i64)).collect()
            };
            assert_eq!(view.changes(c), reaction.net_changes, "{label} {c}");
            assert_eq!(
                counted(view.reactants(c)),
                reaction.reactants,
                "{label} {c}"
            );
            assert_eq!(counted(view.products(c)), reaction.products, "{label} {c}");
            assert_eq!(view.reads(c), reaction.reads, "{label} {c}");
            assert_eq!(view.dependencies(c), network.dependencies[c], "{label} {c}");
            assert_eq!(view.factor(c), 1.0, "{label} {c}");
            let whole = source["reactions"][c]["rate"]["kind"] != "mass_action";
            assert_eq!(view.whole(c), whole, "{label} {c}");
        }
        // Propensities at the initial state, a probe state, the empty state
        // and pseudo-random states.
        let n = initial.len();
        let top = initial.iter().copied().max().unwrap_or(0) * 2 + 20;
        let mut states = vec![
            initial.clone(),
            initial
                .iter()
                .enumerate()
                .map(|(i, &x)| x + i as u64 + 1)
                .collect(),
            vec![0; n],
        ];
        let mut seed = 20_261_111;
        for _ in 0..12 {
            states.push((0..n).map(|_| lcg(&mut seed) % (top + 1)).collect());
        }
        for counts in &states {
            let old_counts: Vec<i64> = counts.iter().map(|&x| x as i64).collect();
            let all = view.propensities(counts).unwrap();
            for (c, reaction) in source["reactions"].as_array().unwrap().iter().enumerate() {
                let Ok(old) = network.propensity(c, &old_counts) else {
                    continue; // the network engine refused an underflow
                };
                let new = view.propensity(c, counts).unwrap();
                assert_eq!(new.to_bits(), all[c].to_bits(), "{label} {c}");
                let limit = allowed(reaction, &source["parameters"], old);
                assert!(
                    (new - old).abs() <= limit,
                    "{label} reaction {c} at {counts:?}: {new} against {old}"
                );
                compared += 1;
                identical += usize::from(new.to_bits() == old.to_bits());
            }
        }
        models += 1;
    }
    assert_eq!(models, 14);
    assert!(compared > 10_000, "{compared}");
    eprintln!("{compared} propensities compared, {identical} bit-identical");
}

fn rules_model(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

/// A closed attributed model: copies of `M` with attributes `wt` and
/// `size`, a conversion into a constant tuple, and an attribute-free type
/// `E` made by the conversion.
fn attributed() -> Value {
    json!({
        "schema": "crn-craft.model/1.0", "name": "attributed",
        "parameters": {"k": 0.5, "big": 1e-7},
        "types": [{"name": "E"}, {"name": "M", "attributes": ["wt", "size"]}],
        "initial": [
            {"type": "M", "attributes": {"wt": 0, "size": 0.5}, "count": 3},
            {"type": "M", "attributes": {"wt": -2, "size": "big"}, "count": 1}
        ],
        "observables": {
            "E": {"type": "E"},
            "mutants": {"type": "M", "where": "x.wt == 1"},
            "total": {"type": "M", "weight": "x.size"}
        },
        "rules": [
            {"name": "mutate", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 0", "factor": "2 * x.size"}],
             "rate": "k", "products": [{"type": "M", "attributes": {"wt": 1, "size": 0.5}}, {"type": "E"}]},
            {"name": "decay", "reactants": [{"label": "x", "type": "E"}],
             "propensity": "k * E / (1 + mutants)", "products": []}
        ]
    })
}

#[test]
fn attributed_species_are_labeled_by_type_and_attribute_values() {
    let view = FlatView::from_bytes(&rules_model(&attributed())).unwrap();
    // By type in declaration order (E first), then in lowering order: the
    // initial tuples, then the tuples products create.
    assert_eq!(
        view.species(),
        [
            "E",
            "M{wt=0,size=0.5}",
            "M{wt=-2,size=1e-7}",
            "M{wt=1,size=0.5}"
        ]
    );
    assert_eq!(view.initial(), [0, 3, 1, 0]);
    assert_eq!(view.index("M{wt=1,size=0.5}"), Some(3));
    assert_eq!(view.index("M"), None);
    assert_eq!(view.species_type(3), 1);
    assert_eq!(view.attributes(3), [1.0, 0.5]);
    // `mutate` matches the wt = 0 tuple only; `decay` has one channel.
    assert_eq!(view.channels(), 2);
    assert_eq!(view.rule(0), 0);
    assert_eq!(view.changes(0), [(0, 1), (1, -1), (3, 1)]);
    assert_eq!(view.reactants(0), [(1, 1)]);
    assert_eq!(view.products(0), [(0, 1), (3, 1)]);
    assert_eq!(view.factor(0), 1.0); // 2 * 0.5
    assert!(!view.whole(0) && view.whole(1));
    assert_eq!(view.reads(1), [0, 3]); // E, and the mutants observable
    assert_eq!(view.dependencies(0), [0, 1]);
    assert_eq!(view.rate_source(1), "k * E / (1 + mutants)");
    assert_eq!(view.rate_constant(0), Some(0.5));
    assert_eq!(view.rate_constant(1), None);
    assert_eq!(view.parameter_names(), ["big", "k"]);
    assert_eq!(view.parameter_values(), [1e-7, 0.5]);
    assert_eq!(view.observable_names(), ["E", "mutants", "total"]);
    assert_eq!(view.rate_observables(1), [0, 1]);
    let (_, terms) = view.observable(2);
    assert_eq!(terms, [(1, 0.5), (2, 1e-7), (3, 0.5)]);
    // k * 3 copies * factor 1, and k * E / (1 + mutants) while E > 0.
    assert_eq!(view.propensity(0, &[0, 3, 1, 0]).unwrap(), 1.5);
    assert_eq!(view.propensity(1, &[0, 3, 1, 0]).unwrap(), 0.0);
    assert_eq!(view.propensity(1, &[4, 1, 1, 3]).unwrap(), 0.5 * 4.0 / 4.0);
    assert_eq!(view.propensities(&[4, 1, 1, 3]).unwrap(), [0.5, 0.5]);
    assert!(!view.has_stop());
    assert!(!view.stops(&[4, 1, 1, 3]).unwrap());
}

#[test]
fn higher_order_and_whole_rule_propensities_follow_static_direct() {
    let model = json!({
        "schema": "crn-craft.model/1.0", "name": "orders",
        "parameters": {"c": 0.25},
        "types": [{"name": "A"}, {"name": "B"}],
        "initial": [{"type": "A", "count": 7}, {"type": "B", "count": 2}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [
            {"name": "triple", "reactants": [{"label": "x", "type": "A", "count": 3}], "rate": "c", "products": []},
            {"name": "mixed", "reactants": [{"label": "x", "type": "A", "count": 2}, {"label": "y", "type": "B"}],
             "rate": "c", "products": [{"type": "B"}, {"type": "B"}]},
            {"name": "whole", "reactants": [{"label": "x", "type": "A", "count": 3}, {"label": "y", "type": "B"}],
             "propensity": "c * B", "products": []}
        ],
        "stop": "B >= 5"
    });
    let view = FlatView::from_bytes(&rules_model(&model)).unwrap();
    assert_eq!(view.reactants(0), [(0, 3)]);
    assert_eq!(view.reactants(1), [(0, 2), (1, 1)]);
    assert_eq!(view.changes(1), [(0, -2), (1, 1)]);
    assert_eq!(view.reactants(2), [(0, 3), (1, 1)]);
    // c C(7,3), c C(7,2) 2, and c B while A >= 3 and B >= 1.
    assert_eq!(view.propensities(&[7, 2]).unwrap(), [8.75, 10.5, 0.5]);
    assert_eq!(view.propensities(&[2, 2]).unwrap(), [0.0, 0.5, 0.0]);
    assert_eq!(view.propensities(&[7, 0]).unwrap(), [8.75, 0.0, 0.0]);
    assert!(view.has_stop());
    assert!(!view.stops(&[7, 4]).unwrap());
    assert!(view.stops(&[7, 5]).unwrap());
}

#[test]
fn models_without_a_flat_view_are_refused_with_r017() {
    let code = |value: &Value| {
        let d = FlatView::from_bytes(&rules_model(value)).unwrap_err();
        (d[0].code.clone(), d[0].path.clone())
    };
    // Attribute arithmetic creates new tuples at run time: not closed.
    let mut growing = attributed();
    growing["rules"][0]["products"][0] =
        json!({"type": "M", "copy": "x", "attributes": {"wt": "x.wt + 1"}});
    assert_eq!(code(&growing), ("CRN-R017".into(), "/rules".into()));
    let mut random = attributed();
    random["rules"][0]["products"][0]["attributes"] = json!({"wt": "bernoulli(0.5)", "size": 0.5});
    assert_eq!(code(&random), ("CRN-R017".into(), "/rules".into()));
    // A whole-rule propensity over a type with several tuples.
    let mut whole = attributed();
    whole["rules"][0] = json!({"name": "mutate", "reactants": [{"label": "x", "type": "M"}],
                               "propensity": "k", "products": []});
    assert_eq!(code(&whole), ("CRN-R017".into(), "/rules".into()));
    let continuous = json!({
        "schema": "crn-craft.rules-model/0.3", "name": "growth",
        "types": [{"name": "X"}], "initial": [{"type": "X", "count": 1}],
        "observables": {"X": {"type": "X"}},
        "continuous": [{"name": "v", "initial": 1.0, "derivative": "0.1 * v"}],
        "rules": [{"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "v", "products": []}]
    });
    assert_eq!(code(&continuous), ("CRN-R017".into(), "/continuous".into()));
    // Model errors keep their own codes.
    let mut invalid = attributed();
    invalid["rules"][0]["rate"] = json!("k +");
    assert_eq!(code(&invalid).0, "CRN-R012");
    let error = FlatView::from_bytes(b"{").unwrap_err();
    assert_eq!(error[0].code, "CRN-E001");
}

#[test]
fn atp_models_have_no_flat_view() {
    let path = root().join("models/rules");
    let mut refused = 0;
    for entry in std::fs::read_dir(path).unwrap() {
        let path = entry.unwrap().path();
        let bytes = std::fs::read(&path).unwrap();
        let Ok(model) = serde_json::from_slice::<Value>(&bytes) else {
            continue;
        };
        if model.get("atp").is_some() {
            let d = FlatView::from_bytes(&bytes).unwrap_err();
            assert_eq!(
                (d[0].code.as_str(), d[0].path.as_str()),
                ("CRN-R017", "/atp")
            );
            refused += 1;
        }
    }
    assert!(refused >= 1);
}

#[test]
fn propensity_errors_are_coded() {
    let view = FlatView::from_bytes(&rules_model(&attributed())).unwrap();
    let d = view.propensity(0, &[1, 2]).unwrap_err();
    assert_eq!((d.code.as_str(), d.path.as_str()), ("CRN-R030", "/counts"));
    let d = view.propensity(0, &[0, (1 << 53) + 1, 0, 0]).unwrap_err();
    assert_eq!(
        (d.code.as_str(), d.path.as_str()),
        ("CRN-R030", "/counts/M{wt=0,size=0.5}")
    );
    // A negative rate fails as it fails a Static Direct replicate.
    let mut negative = attributed();
    negative["rules"][1]["propensity"] = json!("k - E");
    let view = FlatView::from_bytes(&rules_model(&negative)).unwrap();
    // Without an E molecule the propensity is not evaluated.
    assert_eq!(view.propensity(1, &[0, 0, 0, 0]).unwrap(), 0.0);
    let d = view.propensity(1, &[1, 0, 0, 0]).unwrap_err();
    assert_eq!(
        (d.code.as_str(), d.path.as_str()),
        ("CRN-R030", "/rules/1/propensity")
    );
}
