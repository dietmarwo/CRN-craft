//! Deterministic kinetics on the flat view (M5e): the 0.6.0 derivative test
//! ported, models of every version, finite-difference Jacobians of
//! expression rates, higher-order (`count`) channels, whole-rule
//! propensities, the engine's own propensities and the refusals.
use crn_craft::kinetics::{Derivatives, Kinetics};
use crn_craft::parse_model;
use crn_craft::rules::{FlatView, convert};
use serde_json::{Value, json};

fn bytes(model: &Value) -> Vec<u8> {
    serde_json::to_vec(model).unwrap()
}

const EXAMPLES: [&[u8]; 5] = [
    include_bytes!("../../../models/examples/repressilator.model.json"),
    include_bytes!("../../../models/examples/enzyme-mechanism.model.json"),
    include_bytes!("../../../models/examples/michaelis-menten.model.json"),
    include_bytes!("../../../models/examples/immigration-death.model.json"),
    include_bytes!("../../../models/examples/vilar-noise.model.json"),
];

/// Central differences of the right-hand side in every state and parameter
/// component, with relative steps, against the analytic Jacobians.
fn check_jacobians(k: &Kinetics, y: &[f64], p: &[f64], label: &str) {
    let exact = k.evaluate(y, p).unwrap();
    for j in 0..y.len() {
        let h = 1e-5 * y[j].abs().max(1e-3);
        let (mut plus, mut minus) = (y.to_vec(), y.to_vec());
        plus[j] += h;
        minus[j] -= h;
        let a = k.evaluate(&plus, p).unwrap().rhs;
        let b = k.evaluate(&minus, p).unwrap().rhs;
        for i in 0..y.len() {
            let numeric = (a[i] - b[i]) / (2. * h);
            let analytic = exact.state_jacobian[i][j];
            assert!(
                (numeric - analytic).abs() <= 1e-6 * numeric.abs().max(1.),
                "{label}: d rhs[{i}] / d y[{j}] at {y:?}: {analytic} vs {numeric}"
            );
        }
    }
    for j in 0..p.len() {
        let h = 1e-5 * p[j].abs().max(1e-3);
        let (mut plus, mut minus) = (p.to_vec(), p.to_vec());
        plus[j] += h;
        minus[j] -= h;
        let a = k.evaluate(y, &plus).unwrap().rhs;
        let b = k.evaluate(y, &minus).unwrap().rhs;
        for i in 0..y.len() {
            let numeric = (a[i] - b[i]) / (2. * h);
            let analytic = exact.parameter_jacobian[i][j];
            assert!(
                (numeric - analytic).abs() <= 1e-6 * numeric.abs().max(1.),
                "{label}: d rhs[{i}] / d {} at {y:?}: {analytic} vs {numeric}",
                k.parameters[j]
            );
        }
    }
}

/// Three distinct, well-separated positive states.
fn states(n: usize) -> Vec<Vec<f64>> {
    vec![
        (0..n).map(|i| 3.7 + 2.9 * i as f64).collect(),
        (0..n).map(|i| 41.3 - 3.1 * i as f64).collect(),
        (0..n).map(|i| 9.6 + 0.7 * (i * i) as f64).collect(),
    ]
}

fn same(a: &Derivatives, b: &Derivatives) -> bool {
    a.rhs == b.rhs
        && a.state_jacobian == b.state_jacobian
        && a.parameter_jacobian == b.parameter_jacobian
}

/// The 0.6.0 test, unchanged apart from comments: `Kinetics::compile` of a
/// parsed network model, central differences at two distinct states.
#[test]
fn analytic_state_and_parameter_derivatives_match_independent_central_differences() {
    for bytes in &EXAMPLES[..3] {
        let model = parse_model(bytes).unwrap();
        let k = Kinetics::compile(&model).unwrap();
        let p: Vec<_> = k
            .parameters
            .iter()
            .map(|name| model.parameters[name])
            .collect();
        // Distinct components: a derivative taken with respect to the wrong species must not
        // coincide with the right one.
        let n = model.species.len();
        let distinct = [
            (0..n).map(|i| 3.7 + 2.9 * i as f64).collect::<Vec<_>>(),
            (0..n).map(|i| 41.3 - 3.1 * i as f64).collect(),
        ];
        for y in distinct {
            let result = k.evaluate(&y, &p).unwrap();
            for j in 0..y.len() {
                let h = 1e-4;
                let mut plus = y.clone();
                let mut minus = y.clone();
                plus[j] += h;
                minus[j] -= h;
                let a = k.evaluate(&plus, &p).unwrap().rhs;
                let b = k.evaluate(&minus, &p).unwrap().rhs;
                for i in 0..y.len() {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    assert!(
                        (numeric - result.state_jacobian[i][j]).abs()
                            < 1e-7 * numeric.abs().max(1.)
                    );
                }
            }
            for j in 0..p.len() {
                let h = p[j] * 1e-5;
                let mut plus = p.clone();
                let mut minus = p.clone();
                plus[j] += h;
                minus[j] -= h;
                let a = k.evaluate(&y, &plus).unwrap().rhs;
                let b = k.evaluate(&y, &minus).unwrap().rhs;
                for i in 0..y.len() {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    assert!(
                        (numeric - result.parameter_jacobian[i][j]).abs()
                            < 1e-6 * numeric.abs().max(1.)
                    );
                }
            }
        }
        assert!(k.evaluate(&vec![-1.; model.species.len()], &p).is_err());
    }
}

/// A network model, its parsed form and its converted `model/1.0` give the
/// same kinetics bit for bit, with the network's species and parameters.
#[test]
fn network_and_converted_models_give_the_same_kinetics() {
    for bytes in EXAMPLES {
        let network = parse_model(bytes).unwrap();
        let converted = convert::from_network(bytes).unwrap();
        let kinetics = [
            Kinetics::compile(&network).unwrap(),
            Kinetics::from_bytes(bytes).unwrap(),
            Kinetics::from_bytes(&converted).unwrap(),
            Kinetics::new(&FlatView::of_any_model(&converted).unwrap()).unwrap(),
        ];
        let names: Vec<&str> = network.species.iter().map(|s| s.name.as_str()).collect();
        let initial: Vec<u64> = network
            .species
            .iter()
            .map(|s| s.initial_count as u64)
            .collect();
        let parameters: Vec<&String> = network.parameters.keys().collect();
        let values: Vec<f64> = network.parameters.values().copied().collect();
        for k in &kinetics {
            assert_eq!(k.species(), names);
            assert_eq!(k.initial(), initial);
            assert_eq!(k.parameters.iter().collect::<Vec<_>>(), parameters);
            assert_eq!(k.parameter_values(), values);
        }
        for y in states(names.len()) {
            let first = kinetics[0].evaluate(&y, &values).unwrap();
            for k in &kinetics[1..] {
                assert!(
                    same(&first, &k.evaluate(&y, &values).unwrap()),
                    "{}",
                    network.name
                );
            }
            check_jacobians(&kinetics[0], &y, &values, &network.name);
        }
    }
}

/// Expression rates with `^`, `/`, `exp` and `max`, an observable-dependent
/// rate over a pair, a count-2 pattern and a whole-rule propensity that
/// consumes its substrate: analytic Jacobians against central differences.
fn smooth() -> Value {
    json!({
        "schema": "crn-craft.model/1.0", "name": "smooth",
        "parameters": {"K": 8.0, "k1": 0.7, "k2": 0.05, "n": 2.5, "vmax": 12.0},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C"}],
        "initial": [{"type": "A", "count": 20}, {"type": "B", "count": 5}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}, "C": {"type": "C"}},
        "rules": [
            {"name": "make_b", "reactants": [], "propensity": "vmax * A^n / (K^n + A^n)",
             "products": [{"type": "B"}]},
            {"name": "a_to_c", "reactants": [{"label": "x", "type": "A"}],
             "rate": "k1 * exp(-B / K)", "products": [{"type": "C"}]},
            {"name": "bind", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "B"}],
             "rate": "max(k2, k1 * C / (K + C))", "products": [{"type": "C"}]},
            {"name": "pair", "reactants": [{"label": "x", "type": "C", "count": 2}],
             "rate": "k2 / (1 + (A / K)^2)", "products": [{"type": "A"}]},
            {"name": "mm", "reactants": [{"label": "x", "type": "C"}],
             "propensity": "vmax * (C / (K + C))", "products": [{"type": "A"}, {"type": "B"}]}
        ]
    })
}

#[test]
fn expression_rates_have_exact_jacobians() {
    let k = Kinetics::from_bytes(&bytes(&smooth())).unwrap();
    assert_eq!(k.species(), ["A", "B", "C"]);
    assert_eq!(k.parameters, ["K", "k1", "k2", "n", "vmax"]);
    let p = k.parameter_values().to_vec();
    for y in states(3) {
        check_jacobians(&k, &y, &p, "smooth");
    }
    // Both branches of max: k1 C / (K + C) < k2 at C = 0.5.
    check_jacobians(&k, &[6.0, 2.0, 0.5], &p, "smooth, first branch of max");
    // The right-hand side by hand at one state.
    let (a, b, c) = (6.0_f64, 2.0_f64, 3.0_f64);
    let (big_k, k1, k2, n, vmax) = (8.0_f64, 0.7_f64, 0.05_f64, 2.5_f64, 12.0_f64);
    let make_b = vmax * a.powf(n) / (big_k.powf(n) + a.powf(n));
    let a_to_c = k1 * (-b / big_k).exp() * a;
    let bind = f64::max(k2, k1 * c / (big_k + c)) * a * b;
    let pair = k2 / (1.0 + (a / big_k).powi(2)) * c * c / 2.0;
    let mm = vmax * (c / (big_k + c));
    let rhs = k.evaluate(&[a, b, c], &p).unwrap().rhs;
    let expected = [
        -a_to_c - bind + pair + mm,
        make_b - bind + mm,
        a_to_c + bind - 2.0 * pair - mm,
    ];
    for (x, e) in rhs.iter().zip(expected) {
        assert!((x - e).abs() <= 1e-13 * e.abs(), "{rhs:?} vs {expected:?}");
    }
}

/// Higher-order (`Multi`) channels: rate × Π y^k / k!, against closed forms
/// and central differences; the converted higher-order network models.
#[test]
fn higher_order_channels_use_the_macroscopic_form() {
    let model = json!({
        "schema": "crn-craft.model/1.0", "name": "orders",
        "parameters": {"c": 0.25},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C"}],
        "initial": [{"type": "A", "count": 7}, {"type": "B", "count": 2}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [
            {"name": "triple", "reactants": [{"label": "x", "type": "A", "count": 3}], "rate": "c", "products": []},
            {"name": "mixed", "reactants": [{"label": "x", "type": "A", "count": 2}, {"label": "y", "type": "B"}],
             "rate": "c * A", "products": [{"type": "C"}]},
            {"name": "three", "reactants": [{"label": "x", "type": "A"}, {"label": "y", "type": "B"},
                                             {"label": "z", "type": "C"}],
             "rate": "c", "products": [{"type": "C"}, {"type": "C"}]}
        ]
    });
    let k = Kinetics::from_bytes(&bytes(&model)).unwrap();
    let (a, b, c, rate) = (5.5_f64, 3.25_f64, 1.5_f64, 0.25_f64);
    let d = k.evaluate(&[a, b, c], &[rate]).unwrap();
    let triple = rate * a * a * a / 6.0;
    let mixed = rate * a * (a * a / 2.0) * b;
    let three = rate * a * b * c;
    let expected = [
        -3.0 * triple - 2.0 * mixed - three,
        -mixed - three,
        mixed + three,
    ];
    for (x, e) in d.rhs.iter().zip(expected) {
        assert!(
            (x - e).abs() <= 1e-14 * e.abs(),
            "{:?} vs {expected:?}",
            d.rhs
        );
    }
    // d triple / dA = c A^2 / 2; d mixed / dA = c * (3 A^2 / 2) * B.
    let da = -3.0 * rate * a * a / 2.0 - 2.0 * rate * 1.5 * a * a * b - rate * b * c;
    assert!((d.state_jacobian[0][0] - da).abs() <= 1e-13 * da.abs());
    for y in states(3) {
        check_jacobians(&k, &y, &[rate], "orders");
    }
    // Zero counts: the propensity vanishes, the first-order derivative does not.
    let d = k.evaluate(&[0.0, 2.0, 1.0], &[rate]).unwrap();
    assert_eq!(d.rhs, [0.0, 0.0, 0.0]);
    assert_eq!(d.state_jacobian[1][0], -rate * 2.0 * 1.0);
    // The converted network models of M4b: 3A -> 0 and A + 2B -> C.
    for (file, check) in [("third-order-decay", 0usize), ("trimolecular", 1usize)] {
        let path = format!(
            "{}/../../validation/merge/higher-order/{file}.model.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let network = std::fs::read(path).unwrap();
        let k = Kinetics::from_bytes(&network).unwrap();
        let p = k.parameter_values().to_vec();
        for y in states(k.species().len()) {
            check_jacobians(&k, &y, &p, file);
            let rhs = k.evaluate(&y, &p).unwrap().rhs;
            if check == 0 {
                let e = -3.0 * p[0] * y[0].powi(3) / 6.0;
                assert!((rhs[0] - e).abs() <= 1e-14 * e.abs(), "{file}");
            }
        }
    }
    // A rate channel may take at most eight molecules of a species.
    let mut nine = model.clone();
    nine["rules"][0]["reactants"][0]["count"] = json!(9);
    let error = Kinetics::from_bytes(&bytes(&nine)).err().unwrap();
    assert_eq!(error[0].code, "CRN-E080");
    assert!(error[0].message.contains("order above eight"), "{error:?}");
    nine["rules"][0]["reactants"][0]["count"] = json!(8);
    assert!(Kinetics::from_bytes(&bytes(&nine)).is_ok());
}

fn network(rate: Value, reactants: Value, species: Value) -> Vec<u8> {
    bytes(&json!({
        "schema": "crn-craft.model/0.1", "name": "law", "convention": "stochastic_counts",
        "omega": 1.0, "parameters": {"v": 30.0, "K": 8.0, "b": 1.0},
        "species": species,
        "reactions": [{"name": "law", "reactants": reactants, "products": {"Y": 1}, "rate": rate}]
    }))
}

/// Whole-rule propensities: the converted Michaelis–Menten law as in
/// 0.6.0, Hill laws at a regulator count of 0, the availability gate
/// dropped, and the refusal of propensities that cannot be extended.
#[test]
fn whole_rule_propensities_extend_their_gate_only_where_they_vanish() {
    // Michaelis–Menten: v S / (km + S), consumed substrate, bit for bit.
    let k = Kinetics::from_bytes(EXAMPLES[2]).unwrap();
    let (vmax, km) = (5.0, 20.0);
    for s in [0.0, 0.5, 50.0] {
        let d = k.evaluate(&[s, 3.0], &[km, vmax]).unwrap();
        assert_eq!(d.rhs, [-(vmax * (s / (km + s))), vmax * (s / (km + s))]);
        let ds = vmax * km / ((km + s) * (km + s));
        assert!((d.state_jacobian[1][0] - ds).abs() <= 1e-14 * ds);
        assert_eq!(d.state_jacobian[1][1], 0.0);
    }
    // Hill activation at R = 0: one-sided derivatives vmax / K (n = 1) and
    // 0 (n > 1); n < 1 is singular at zero.
    let species = json!([{"name": "R", "initial_count": 0}, {"name": "Y", "initial_count": 0}]);
    for (n, expected) in [(1.0, 30.0 / 8.0), (2.0, 0.0), (2.5, 0.0)] {
        let rate = json!({"kind": "hill_activation", "regulator": "R", "vmax": "v",
                          "basal": "b", "k": "K", "n": n});
        let k = Kinetics::from_bytes(&network(rate, json!({}), species.clone())).unwrap();
        let p = k.parameter_values().to_vec();
        let d = k.evaluate(&[0.0, 0.0], &p).unwrap();
        assert_eq!(d.rhs, [0.0, 1.0]);
        assert_eq!(d.state_jacobian[1][0], expected, "n = {n}");
        // Parameters K, b, v: only the basal rate moves the rate at R = 0.
        assert_eq!(d.parameter_jacobian[1], [0.0, 1.0, 0.0], "n = {n}");
        check_jacobians(&k, &[3.0, 1.0], &p, "hill");
    }
    let rate = json!({"kind": "hill_activation", "regulator": "R", "vmax": "v",
                      "basal": "b", "k": "K", "n": 0.5});
    let k = Kinetics::from_bytes(&network(rate.clone(), json!({}), species.clone())).unwrap();
    let error = k.evaluate(&[0.0, 0.0], k.parameter_values()).err().unwrap();
    assert_eq!(error.code, "CRN-E080");
    assert!(error.message.contains("singular at zero"), "{error:?}");
    assert!(k.evaluate(&[0.5, 0.0], k.parameter_values()).is_ok());
    // A Hill law with reactants cannot drop its gate: refused, as in 0.6.0.
    let error = Kinetics::from_bytes(&network(rate, json!({"R": 1}), species.clone()))
        .err()
        .unwrap();
    assert_eq!(error[0].code, "CRN-E080");
    assert_eq!(error[0].path, "/rules/0/propensity");
    assert!(error[0].message.contains("availability gate"), "{error:?}");
    // Michaelis–Menten consuming more than its substrate: refused (W4).
    let mm = json!({"kind": "michaelis_menten", "substrate": "R", "vmax": "v", "km": "K"});
    let error = Kinetics::from_bytes(&network(
        mm.clone(),
        json!({"R": 1, "Y": 1}),
        species.clone(),
    ))
    .err()
    .unwrap();
    assert!(error[0].message.contains("molecules of Y"), "{error:?}");
    // ... and consuming two of its substrate: accepted, without the gate:
    // at R = 1 the stochastic propensity is 0, the ODE's v R / (K + R).
    let text = network(mm, json!({"R": 2}), species);
    let k = Kinetics::from_bytes(&text).unwrap();
    let rhs = k.evaluate(&[1.0, 0.0], k.parameter_values()).unwrap().rhs;
    assert_eq!(rhs, [-2.0 * (30.0 * (1.0 / 9.0)), 30.0 * (1.0 / 9.0)]);
    let view = FlatView::of_any_model(&text).unwrap();
    assert_eq!(view.propensity(0, &[1, 0]).unwrap(), 0.0);
    let w4 = std::fs::read(format!(
        "{}/../../validation/ssa-acceleration/models/W4.model.json",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let error = Kinetics::from_bytes(&w4).err().unwrap();
    assert!(error.iter().all(|d| d.code == "CRN-E080"), "{error:?}");
}

/// Channels that take at most one molecule of each species have the
/// engine's own propensities: Σ changes × `FlatView::propensity` at integer
/// states, for converted models and an attributed model with pattern
/// factors, weighted observables and a whole-rule propensity.
#[test]
fn first_order_channels_have_the_engine_propensities() {
    let attributed = json!({
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
            {"name": "mutate", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 0", "factor": "3 * x.size"}],
             "rate": "k * (1 + total)", "products": [{"type": "M", "attributes": {"wt": 1, "size": 0.5}}, {"type": "E"}]},
            {"name": "decay", "reactants": [{"label": "x", "type": "E"}],
             "propensity": "k * E / (1 + mutants)", "products": []}
        ]
    });
    let mut models = vec![bytes(&attributed)];
    models.extend(EXAMPLES[..2].iter().map(|b| b.to_vec()));
    models.push(EXAMPLES[3].to_vec());
    models.push(EXAMPLES[4].to_vec());
    for model in models {
        let k = Kinetics::from_bytes(&model).unwrap();
        let view = FlatView::of_any_model(&model).unwrap();
        let n = k.species().len();
        assert_eq!(k.species(), view.species());
        for seed in 0..6u64 {
            let counts: Vec<u64> = (0..n as u64).map(|i| (seed * 7 + i * 3) % 11 + 1).collect();
            let y: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
            let rhs = k.evaluate(&y, k.parameter_values()).unwrap().rhs;
            let mut expected = vec![0.0; n];
            for c in 0..view.channels() {
                let a = view.propensity(c, &counts).unwrap();
                for &(s, d) in view.changes(c) {
                    expected[s] += d as f64 * a;
                }
            }
            for (x, e) in rhs.iter().zip(&expected) {
                assert!(
                    (x - e).abs() <= 4.0 * f64::EPSILON * e.abs(),
                    "{rhs:?} vs {expected:?}"
                );
            }
            check_jacobians(&k, &y, k.parameter_values(), "first order");
        }
    }
}

#[test]
fn invalid_inputs_and_models_are_refused_with_coded_diagnostics() {
    let k = Kinetics::from_bytes(&bytes(&smooth())).unwrap();
    let p = k.parameter_values().to_vec();
    for (y, p) in [
        (vec![1.0, 2.0], p.clone()),
        (vec![1.0, 2.0, -1.0], p.clone()),
        (vec![1.0, f64::NAN, 1.0], p.clone()),
        (vec![1.0, 2.0, 3.0], p[..4].to_vec()),
        (vec![1.0, 2.0, 3.0], vec![1.0, 1.0, f64::INFINITY, 1.0, 1.0]),
    ] {
        let error = k.evaluate(&y, &p).err().unwrap();
        assert_eq!(
            (error.code.as_str(), error.path.as_str()),
            ("CRN-E080", "/ode")
        );
        assert!(error.message.contains("invalid ODE state"), "{error:?}");
    }
    // Rules parameters may be negative; a negative propensity is refused.
    let mut model = smooth();
    model["rules"][0]["propensity"] = json!("vmax * (5 - A)");
    let k = Kinetics::from_bytes(&bytes(&model)).unwrap();
    assert!(k.evaluate(&[4.0, 1.0, 1.0], k.parameter_values()).is_ok());
    let error = k
        .evaluate(&[6.0, 1.0, 1.0], k.parameter_values())
        .err()
        .unwrap();
    assert!(
        error.message.contains("rule `make_b` is negative"),
        "{error:?}"
    );
    let mut negative = k.parameter_values().to_vec();
    negative[4] = -1.0;
    let error = k.evaluate(&[4.0, 1.0, 1.0], &negative).err().unwrap();
    assert!(error.message.contains("is negative"), "{error:?}");
    // A propensity that is not finite, and an observable that overflows.
    model["rules"][0]["propensity"] = json!("vmax / (A - 4)");
    let k = Kinetics::from_bytes(&bytes(&model)).unwrap();
    let error = k
        .evaluate(&[4.0, 1.0, 1.0], k.parameter_values())
        .err()
        .unwrap();
    assert!(error.message.contains("not finite"), "{error:?}");
    model["observables"]["A"] = json!({"type": "A", "weight": "1e300"});
    model["rules"][0]["propensity"] = json!("vmax");
    let k = Kinetics::from_bytes(&bytes(&model)).unwrap();
    let error = k
        .evaluate(&[1e10, 1.0, 1.0], k.parameter_values())
        .err()
        .unwrap();
    assert!(
        error.message.contains("observable `A` overflows"),
        "{error:?}"
    );
    // The size cap: 65 species.
    let types: Vec<Value> = (0..65).map(|i| json!({"name": format!("S{i}")})).collect();
    let large = json!({
        "schema": "crn-craft.model/1.0", "name": "large", "parameters": {},
        "types": types, "initial": [], "observables": {},
        "rules": [{"name": "make", "reactants": [], "rate": "1", "products": [{"type": "S0"}]}]
    });
    let error = Kinetics::from_bytes(&bytes(&large)).err().unwrap();
    assert!(error[0].message.contains("ODE size cap"), "{error:?}");
    // Models without a flat view, and documents that do not parse.
    let mut hybrid = smooth();
    hybrid["continuous"] = json!([{"name": "z", "initial": 0.0, "derivative": "1"}]);
    let error = Kinetics::from_bytes(&bytes(&hybrid)).err().unwrap();
    assert_eq!(error[0].code, "CRN-R017", "{error:?}");
    assert_eq!(
        Kinetics::from_bytes(b"{").err().unwrap()[0].code,
        "CRN-E001"
    );
    // A network model's population limits are ignored, as in 0.6.0.
    let mut limited: Value = serde_json::from_slice(EXAMPLES[3]).unwrap();
    limited["species"][0]["population_limit"] = json!(5);
    let k = Kinetics::from_bytes(&bytes(&limited)).unwrap();
    let plain = Kinetics::from_bytes(EXAMPLES[3]).unwrap();
    let p = plain.parameter_values();
    assert!(same(
        &k.evaluate(&[50.0], p).unwrap(),
        &plain.evaluate(&[50.0], p).unwrap()
    ));
}
