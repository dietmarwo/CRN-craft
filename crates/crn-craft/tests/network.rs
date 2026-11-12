//! The network engine's compiled networks and rate laws. The conservation
//! tests moved to tests/structure.rs with the analysis (M5b).
use crn_craft::{CompiledNetwork, contracts::Model, parse_model};
use serde_json::{Value, json};

fn make(species: &[(&str, i64)], reactions: Vec<Value>) -> Model {
    parse_model(&serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"test","convention":"stochastic_counts","omega":1.0,"parameters":{},"species":species.iter().map(|(name,count)|json!({"name":name,"initial_count":count})).collect::<Vec<_>>(),"reactions":reactions})).unwrap()).unwrap()
}
fn reaction(name: &str, left: Value, right: Value, rate: Value) -> Value {
    json!({"name":name,"reactants":left,"products":right,"rate":rate})
}
fn ma(c: f64) -> Value {
    json!({"kind":"mass_action","coefficient":c})
}
fn hill(repression: bool, regulator: &str, vmax: f64, basal: f64, k: f64, n: f64) -> Value {
    json!({"kind":if repression{"hill_repression"}else{"hill_activation"},"regulator":regulator,"vmax":vmax,"basal":basal,"k":k,"n":n})
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs() + 8.0 * f64::from_bits(1),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn homodimer_all_counts_exact() {
    for c in [1.0, 0.5] {
        let network = CompiledNetwork::compile(&make(
            &[("A", 0)],
            vec![reaction("dimer", json!({"A":2}), json!({}), ma(c))],
        ))
        .unwrap();
        for a in 0_i64..=50 {
            assert_eq!(
                network.propensity(0, &[a]).unwrap(),
                c * (a * (a - 1) / 2) as f64
            );
        }
    }
}

#[test]
fn heterodimer_all_count_pairs_exact() {
    for c in [1.0, 0.5] {
        let network = CompiledNetwork::compile(&make(
            &[("A", 0), ("B", 0)],
            vec![reaction("bind", json!({"A":1,"B":1}), json!({}), ma(c))],
        ))
        .unwrap();
        for a in 0..=50 {
            for b in 0..=50 {
                assert_eq!(network.propensity(0, &[a, b]).unwrap(), c * (a * b) as f64);
            }
        }
    }
}

#[test]
fn high_order_and_small_coefficient_do_not_overflow_intermediates() {
    let network = CompiledNetwork::compile(&make(
        &[("A", 0)],
        vec![reaction("many", json!({"A":1000}), json!({}), ma(1e-300))],
    ))
    .unwrap();
    close(
        network.propensity(0, &[2000]).unwrap(),
        2.048151626989489e300,
    );
    assert_eq!(network.propensity(0, &[1000]).unwrap(), 1e-300);
    assert_eq!(network.propensity(0, &[999]).unwrap(), 0.0);
}

#[test]
fn representable_subnormal_source_rate_is_kept() {
    let c = f64::from_bits(1);
    let network = CompiledNetwork::compile(&make(
        &[("A", 0)],
        vec![reaction("source", json!({}), json!({"A":1}), ma(c))],
    ))
    .unwrap();
    assert_eq!(network.propensity(0, &[0]).unwrap().to_bits(), 1);
}

#[test]
fn mass_action_true_overflow_is_coded_and_bounded_in_work() {
    let network = CompiledNetwork::compile(&make(
        &[("A", 0)],
        vec![reaction(
            "many",
            json!({"A":2147483647}),
            json!({}),
            ma(f64::from_bits(1)),
        )],
    ))
    .unwrap();
    assert_eq!(
        network.propensity(0, &[i64::MAX]).unwrap_err().code,
        "CRN-E042"
    );
}

#[test]
fn hill_complement_and_zero_regulator() {
    let network = CompiledNetwork::compile(&make(
        &[("X", 0)],
        vec![
            reaction(
                "a",
                json!({}),
                json!({"X":1}),
                hill(false, "X", 7.0, 0.25, 10.0, 2.0),
            ),
            reaction(
                "r",
                json!({}),
                json!({"X":1}),
                hill(true, "X", 7.0, 0.25, 10.0, 2.0),
            ),
        ],
    ))
    .unwrap();
    for x in 0..=50 {
        let p = network.propensities(&[x]).unwrap();
        close(p[0], 0.25 + 7.0 * (x * x) as f64 / (100 + x * x) as f64);
        close(p[0] + p[1], 7.5);
    }
    assert_eq!(network.propensities(&[0]).unwrap(), [0.25, 7.25]);
}

#[test]
fn hill_extreme_parameters_and_tails_are_checked() {
    let network = CompiledNetwork::compile(&make(
        &[("X", 1)],
        vec![reaction(
            "a",
            json!({}),
            json!({"X":1}),
            hill(false, "X", 1e300, 0.0, 1e300, 2.0),
        )],
    ))
    .unwrap();
    close(network.propensity(0, &[1]).unwrap(), 1e-300);
    let network = CompiledNetwork::compile(&make(
        &[("X", 1)],
        vec![reaction(
            "r",
            json!({}),
            json!({"X":1}),
            hill(true, "X", 1.0, 0.0, f64::from_bits(1), 2.0),
        )],
    ))
    .unwrap();
    assert_eq!(network.propensity(0, &[1]).unwrap_err().code, "CRN-E042");
    let network = CompiledNetwork::compile(&make(
        &[("X", 0)],
        vec![reaction(
            "r",
            json!({}),
            json!({"X":1}),
            hill(true, "X", f64::MAX, f64::MAX, 1.0, 2.0),
        )],
    ))
    .unwrap();
    assert_eq!(network.propensity(0, &[0]).unwrap_err().code, "CRN-E042");
}

#[test]
fn steep_hill_preserves_near_equal_arguments() {
    let k = f64::from_bits(1.0_f64.to_bits() + 1);
    let network = CompiledNetwork::compile(&make(
        &[("X", 1)],
        vec![reaction(
            "a",
            json!({}),
            json!({"X":1}),
            hill(false, "X", 1.0, 0.0, k, 1e15),
        )],
    ))
    .unwrap();
    close(network.propensity(0, &[1]).unwrap(), 0.44471580577189446);
}

#[test]
fn mm_avoids_multiplication_overflow_and_reports_underflow() {
    let rate = json!({"kind":"michaelis_menten","substrate":"S","vmax":1e308,"km":2.0});
    let network = CompiledNetwork::compile(&make(
        &[("S", 4)],
        vec![reaction("convert", json!({"S":1}), json!({}), rate)],
    ))
    .unwrap();
    close(network.propensity(0, &[4]).unwrap(), 1e308 * (2.0 / 3.0));
    assert_eq!(network.propensity(0, &[0]).unwrap(), 0.0);
    let rate =
        json!({"kind":"michaelis_menten","substrate":"S","vmax":f64::from_bits(1),"km":1e308});
    let network = CompiledNetwork::compile(&make(
        &[("S", 1)],
        vec![reaction("convert", json!({"S":1}), json!({}), rate)],
    ))
    .unwrap();
    assert_eq!(network.propensity(0, &[1]).unwrap_err().code, "CRN-E042");
}

#[test]
fn custom_constant_rate_still_requires_all_reactants() {
    let network = CompiledNetwork::compile(&make(
        &[("E", 0), ("P", 0)],
        vec![reaction(
            "custom",
            json!({"E":1}),
            json!({"E":1,"P":1}),
            hill(true, "P", 0.0, 2.0, 1.0, 1.0),
        )],
    ))
    .unwrap();
    assert_eq!(network.propensities(&[0, 0]).unwrap(), [0.0]);
    assert_eq!(network.propensities(&[1, 0]).unwrap(), [2.0]);
    let mut state = [0, 0];
    assert_eq!(
        network.fire_reaction(0, &mut state).unwrap_err().code,
        "CRN-E043"
    );
    assert_eq!(state, [0, 0]);
}

#[test]
fn unavailable_gate_prevents_evaluating_overflowing_custom_rate() {
    let network = CompiledNetwork::compile(&make(
        &[("E", 0), ("X", 0)],
        vec![reaction(
            "custom",
            json!({"E":1}),
            json!({}),
            hill(true, "X", f64::MAX, f64::MAX, 1.0, 1.0),
        )],
    ))
    .unwrap();
    assert_eq!(network.propensities(&[0, 0]).unwrap(), [0.0]);
}

#[test]
fn indexed_columns_follow_declaration_order_and_have_correct_signs() {
    let network = CompiledNetwork::compile(&make(
        &[("Z", 2), ("A", 3), ("M", 0)],
        vec![reaction(
            "r",
            json!({"A":2,"Z":1}),
            json!({"A":1,"M":3}),
            ma(1.0),
        )],
    ))
    .unwrap();
    let r = &network.reactions()[0];
    assert_eq!(
        r.reactants
            .iter()
            .map(|t| (t.species, t.coefficient))
            .collect::<Vec<_>>(),
        [(0, 1), (1, 2)]
    );
    assert_eq!(
        r.net_changes
            .iter()
            .map(|t| (t.species, t.coefficient))
            .collect::<Vec<_>>(),
        [(0, -1), (1, -1), (2, 3)]
    );
    let mut counts = [2, 3, 0];
    network.fire_reaction(0, &mut counts).unwrap();
    assert_eq!(counts, [1, 2, 3]);
}

#[test]
fn dependency_graph_includes_custom_reads_and_gate_reads() {
    let network = CompiledNetwork::compile(&make(
        &[("E", 1), ("R", 1), ("P", 0)],
        vec![
            reaction("make_e", json!({}), json!({"E":1}), ma(1.0)),
            reaction("make_r", json!({}), json!({"R":1}), ma(1.0)),
            reaction(
                "custom",
                json!({"E":1}),
                json!({"E":1,"P":1}),
                hill(false, "R", 2.0, 1.0, 1.0, 2.0),
            ),
            reaction("consume", json!({"P":1}), json!({}), ma(1.0)),
        ],
    ))
    .unwrap();
    assert_eq!(network.reactions()[2].reads, [0, 1]);
    assert_eq!(network.dependencies(), [vec![2], vec![2], vec![3], vec![3]]);
}

#[test]
fn no_change_reaction_invalidates_no_propensity() {
    let network = CompiledNetwork::compile(&make(
        &[("A", 1)],
        vec![reaction("same", json!({"A":1}), json!({"A":1}), ma(1.0))],
    ))
    .unwrap();
    assert!(network.reactions()[0].net_changes.is_empty());
    assert!(network.dependencies()[0].is_empty());
}

#[test]
fn compiled_parameters_do_not_follow_later_model_edits() {
    let mut model = make(
        &[("A", 1)],
        vec![reaction(
            "r",
            json!({"A":1}),
            json!({}),
            json!({"kind":"mass_action","coefficient":"c"}),
        )],
    );
    model.parameters.insert("c".into(), 2.0);
    let network = CompiledNetwork::compile(&model).unwrap();
    model.parameters.insert("c".into(), 3.0);
    assert_eq!(network.propensity(0, &[1]).unwrap(), 2.0);
}

#[test]
fn compile_rejects_invalid_typed_model_before_indexing() {
    let mut model = make(
        &[("A", 1)],
        vec![reaction("r", json!({"ghost":1}), json!({}), ma(1.0))],
    );
    assert!(
        CompiledNetwork::compile(&model)
            .unwrap_err()
            .iter()
            .any(|d| d.code == "CRN-E012")
    );
    model.omega = f64::NAN;
    assert!(CompiledNetwork::compile(&model).is_err());
}

#[test]
fn zero_rate_and_invalid_indices_cannot_fire() {
    let network = CompiledNetwork::compile(&make(
        &[("A", 1)],
        vec![reaction("r", json!({"A":1}), json!({}), ma(0.0))],
    ))
    .unwrap();
    for reaction in [0, 1, usize::MAX] {
        let mut state = [1];
        assert_eq!(
            network
                .fire_reaction(reaction, &mut state)
                .unwrap_err()
                .code,
            "CRN-E043"
        );
        assert_eq!(state, [1]);
    }
}

#[test]
fn state_domains_checked_before_evaluation() {
    let mut model = make(&[("A", 1)], vec![]);
    model.species[0].population_limit = Some(3);
    let network = CompiledNetwork::compile(&model).unwrap();
    for state in [vec![], vec![1, 2], vec![-1], vec![4]] {
        assert_eq!(network.propensities(&state).unwrap_err().code, "CRN-E041");
    }
    assert!(network.propensities(&[3]).unwrap().is_empty());
}

#[test]
fn named_state_requires_exact_species_set() {
    let network = CompiledNetwork::compile(&make(&[("A", 1)], vec![])).unwrap();
    for value in [json!({}), json!({"A":1,"B":2}), json!({"A":-1})] {
        assert_eq!(
            network
                .named_state(&serde_json::from_value(value).unwrap())
                .unwrap_err()
                .code,
            "CRN-E041"
        );
    }
}

#[test]
fn count_overflow_is_transactional() {
    let network = CompiledNetwork::compile(&make(
        &[("A", 1), ("B", i64::MAX)],
        vec![reaction("convert", json!({"A":1}), json!({"B":1}), ma(1.0))],
    ))
    .unwrap();
    let mut state = [1, i64::MAX];
    assert_eq!(
        network.fire_reaction(0, &mut state).unwrap_err().code,
        "CRN-E044"
    );
    assert_eq!(state, [1, i64::MAX]);
}

#[test]
fn count_updates_remain_exact_above_two_to_the_53() {
    let count = 9_007_199_254_740_993_i64;
    let network = CompiledNetwork::compile(&make(
        &[("A", count)],
        vec![reaction("source", json!({}), json!({"A":1}), ma(1.0))],
    ))
    .unwrap();
    let mut state = [count];
    network.fire_reaction(0, &mut state).unwrap();
    assert_eq!(state, [count + 1]);
}

#[test]
fn population_limit_is_transactional() {
    let mut model = make(
        &[("A", 1), ("B", 2)],
        vec![reaction("convert", json!({"A":1}), json!({"B":1}), ma(1.0))],
    );
    model.species[1].population_limit = Some(2);
    let network = CompiledNetwork::compile(&model).unwrap();
    let mut state = [1, 2];
    assert_eq!(
        network.fire_reaction(0, &mut state).unwrap_err().code,
        "CRN-E044"
    );
    assert_eq!(state, [1, 2]);
}

#[test]
fn dependency_edge_limit_is_explicit() {
    let reactions: Vec<_> = (0..1001)
        .map(|i| reaction(&format!("r{i}"), json!({"A":1}), json!({"A":2}), ma(1.0)))
        .collect();
    let at_limit = CompiledNetwork::compile(&make(&[("A", 1)], reactions[..1000].to_vec()))
        .expect("exactly one million dependency edges must be supported");
    assert_eq!(
        at_limit.dependencies().iter().map(Vec::len).sum::<usize>(),
        1_000_000
    );
    assert!(
        CompiledNetwork::compile(&make(&[("A", 1)], reactions))
            .unwrap_err()
            .iter()
            .any(|d| d.code == "CRN-E040")
    );
}

#[test]
fn hill_near_threshold_uses_log1p_not_a_difference_of_logs() {
    // ln(x) - ln(k) loses about 2e-7 relative accuracy here; with n = 1e8 that
    // moves the fraction far beyond 2e-12. References: mpmath at 60 digits.
    let (x, k, n) = (1_000_000_007_i64, 1_000_000_000.5, 1e8);
    for (repression, expected) in [(false, 0.6570104621242138), (true, 0.3429895378757863)] {
        let network = CompiledNetwork::compile(&make(
            &[("X", 0)],
            vec![reaction(
                "r",
                json!({}),
                json!({}),
                hill(repression, "X", 1.0, 0.0, k, n),
            )],
        ))
        .unwrap();
        close(network.propensity(0, &[x]).unwrap(), expected);
    }
}

#[test]
fn high_order_reactant_uses_binomial_symmetry() {
    // binomial(k + 1, k) = k + 1 must come from one exact factor. A product of
    // k rounded ratios is inexact (relative error about 4e-13 for k = 10^7) and
    // takes k steps; the first reaction fails fast if symmetry is lost.
    let orders = [10_000_000_i64, i64::from(i32::MAX)];
    let network = CompiledNetwork::compile(&make(
        &[("A", 0), ("B", 0)],
        vec![
            reaction("r", json!({"A": orders[0]}), json!({}), ma(1.0)),
            reaction("s", json!({"B": orders[1]}), json!({}), ma(1.0)),
        ],
    ))
    .unwrap();
    for (reaction, order) in orders.into_iter().enumerate() {
        let state = |a: i64| if reaction == 0 { [a, 0] } else { [0, a] };
        let rate = |a: i64| network.propensity(reaction, &state(a)).unwrap();
        assert_eq!(rate(order + 1), (order + 1) as f64);
        assert_eq!(rate(order), 1.0);
        assert_eq!(rate(order - 1), 0.0);
    }
}
