//! Conservation laws, the non-negative conservation cone and the complex
//! graph of flat views (M5b), ported from the 0.6.0 network tests: the
//! networks are `crn-craft.model/1.0` models of attribute-free types whose
//! rule r has the net changes `columns[r]`.
//!
//! Two kinds of 0.6.0 fixture cannot be written as rules models and are
//! adapted: a model without reactions gets one rule without reactants and
//! products (a rules model needs a rule; a zero column leaves the
//! conservation space and the cone unchanged), and coefficients of
//! 2147483647 sit on the reactant side (a rule has at most 64 products).
use crn_craft::{
    analyze_bytes,
    cone::{
        ConeLimits, ConeStatus, MAX_CONE_RAYS, MAX_CONE_WORK, analyze_cone,
        analyze_cone_with_limits,
    },
    conservation::{ConservationStatus, analyze_conservation},
    linkage::analyze_linkage,
    rules::{FlatView, model::Model},
};
use serde_json::{Value, json};

const LARGE: i64 = 2_147_483_647;

/// A pattern taking `count` molecules of `name`.
fn pattern(label: usize, name: &str, count: i64) -> Value {
    if count == 1 {
        json!({"label": format!("x{label}"), "type": name})
    } else {
        json!({"label": format!("x{label}"), "type": name, "count": count})
    }
}

/// A rule from reactant and product coefficients by species name.
fn rule(name: &str, left: &[(&str, i64)], right: &[(&str, i64)], rate: &str) -> Value {
    let reactants: Vec<Value> = left
        .iter()
        .enumerate()
        .map(|(i, (s, c))| pattern(i, s, *c))
        .collect();
    let products: Vec<Value> = right
        .iter()
        .flat_map(|(s, c)| (0..*c).map(move |_| json!({"type": s})))
        .collect();
    json!({"name": name, "reactants": reactants, "rate": rate, "products": products})
}

/// The model with types `names`, initial `counts` and one rule `r{r}` with
/// the net changes of each column (negative entries are reactants).
fn model(names: &[&str], counts: &[u64], columns: &[Vec<i64>]) -> Value {
    let rules: Vec<Value> = columns
        .iter()
        .enumerate()
        .map(|(r, column)| {
            let left: Vec<(&str, i64)> = names
                .iter()
                .zip(column)
                .filter(|(_, c)| **c < 0)
                .map(|(n, c)| (*n, -c))
                .collect();
            let right: Vec<(&str, i64)> = names
                .iter()
                .zip(column)
                .filter(|(_, c)| **c > 0)
                .map(|(n, c)| (*n, *c))
                .collect();
            rule(&format!("r{r}"), &left, &right, "1")
        })
        .collect();
    with_rules(names, counts, rules)
}

fn with_rules(names: &[&str], counts: &[u64], rules: Vec<Value>) -> Value {
    json!({
        "schema": "crn-craft.model/1.0", "name": "structure",
        "types": names.iter().map(|n| json!({"name": n})).collect::<Vec<_>>(),
        "initial": names.iter().zip(counts).filter(|(_, c)| **c > 0)
            .map(|(n, c)| json!({"type": n, "count": c})).collect::<Vec<_>>(),
        "rules": rules
    })
}

/// A column without net change: the one rule of a model "without reactions".
fn idle(n: usize) -> Vec<Vec<i64>> {
    vec![vec![0; n]]
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

fn view(value: &Value) -> FlatView {
    let model: Model = serde_json::from_value(value.clone()).unwrap();
    FlatView::with_every_type(&model).unwrap()
}

fn counterexample() -> Value {
    model(
        &["X0", "X1", "X2", "X3"],
        &[0, 0, 0, 1],
        &[vec![1, 2, -1, 0], vec![1, 0, 2, -1]],
    )
}
fn adjacency_fixture() -> Value {
    model(
        &["A", "B", "C", "D"],
        &[1, 2, 3, 4],
        &[vec![1, 1, -1, -1], vec![1, -1, 1, -1]],
    )
}

/// The chain X0 -> X1 -> ... in which each step consumes `LARGE` molecules,
/// so the conservation weights grow as powers of `LARGE`.
fn large_chain(n: usize) -> Vec<Vec<i64>> {
    (0..n - 1)
        .map(|i| {
            let mut column = vec![0; n];
            column[i] = -LARGE;
            column[i + 1] = 1;
            column
        })
        .collect()
}

#[test]
fn every_type_has_its_species_in_the_analysis_view() {
    // X0, X1 and X2 never hold molecules, yet they are the species of the
    // network; the default view (for FSP) leaves them out.
    let value = counterexample();
    let model: Model = serde_json::from_value(value.clone()).unwrap();
    let complete = FlatView::with_every_type(&model).unwrap();
    assert_eq!(complete.species(), ["X0", "X1", "X2", "X3"]);
    assert_eq!(complete.initial(), [0, 0, 0, 1]);
    assert_eq!(complete.channels(), 2);
    assert_eq!(complete.changes(0), [(0, 1), (1, 2), (2, -1)]);
    let lowered = FlatView::new(&model).unwrap();
    assert_eq!(lowered.species(), ["X0", "X1", "X2", "X3"]);
    let unused = with_rules(
        &["A", "Z"],
        &[2, 0],
        vec![
            rule("death", &[("A", 1)], &[], "1"),
            rule("unused", &[("Z", 1)], &[("A", 1)], "1"),
        ],
    );
    let model: Model = serde_json::from_value(unused).unwrap();
    let lowered = FlatView::new(&model).unwrap();
    assert_eq!(
        (lowered.species(), lowered.channels()),
        (&["A".to_string()][..], 1)
    );
    let complete = FlatView::with_every_type(&model).unwrap();
    assert_eq!(complete.species(), ["A", "Z"]);
    assert_eq!(complete.channels(), 2);
    assert_eq!(complete.changes(1), [(0, 1), (1, -1)]);
    assert_eq!(complete.propensities(&[2, 0]).unwrap(), [2.0, 0.0]);
    assert_eq!(complete.propensities(&[2, 3]).unwrap(), [2.0, 3.0]);
}

#[test]
fn review_counterexample_gets_complete_bounds_and_suppresses_growth_warning() {
    let value = counterexample();
    let network = view(&value);
    assert!(
        analyze_conservation(&network)
            .initial_class_upper_bounds
            .iter()
            .all(Option::is_none)
    );
    let cone = analyze_cone(&network);
    assert_eq!(cone.status, ConeStatus::Exact);
    assert_eq!(
        cone.extreme_rays
            .iter()
            .map(|r| r.coefficients.clone())
            .collect::<Vec<_>>(),
        vec![vec!["0", "1", "2", "4"], vec!["1", "0", "1", "3"]]
    );
    assert_eq!(
        cone.initial_class_upper_bounds,
        [
            Some("3".into()),
            Some("4".into()),
            Some("2".into()),
            Some("1".into())
        ]
    );
    let report = analyze_bytes(&bytes(&value), None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert!(!report.diagnostics.iter().any(|d| d.code == "CRN-W001"));
}
#[test]
fn adjacency_excludes_nonextreme_combination() {
    let cone = analyze_cone(&view(&adjacency_fixture()));
    assert_eq!(cone.status, ConeStatus::Exact);
    assert_eq!(
        cone.extreme_rays
            .iter()
            .map(|r| r.coefficients.clone())
            .collect::<Vec<_>>(),
        vec![vec!["0", "1", "1", "0"], vec!["1", "0", "0", "1"]]
    );
}
#[test]
fn weighted_cone_bounds_use_integer_division() {
    let cone = analyze_cone(&view(&model(
        &["A", "B"],
        &[3, 2],
        &[vec![-2, 1], vec![2, -1]],
    )));
    assert_eq!(cone.extreme_rays[0].coefficients, ["1", "2"]);
    assert_eq!(cone.extreme_rays[0].initial_total, "7");
    assert_eq!(
        cone.initial_class_upper_bounds,
        [Some("7".into()), Some("3".into())]
    );
}
#[test]
fn duplicate_reversed_and_null_constraints_preserve_rays() {
    let mut value = counterexample();
    let expected = analyze_cone(&view(&value));
    let rules = value["rules"].as_array_mut().unwrap();
    let mut duplicate = rules[0].clone();
    duplicate["name"] = json!("duplicate");
    rules.push(duplicate);
    rules.push(rule("reverse", &[("X0", 1), ("X1", 2)], &[("X2", 1)], "1"));
    let actual = analyze_cone(&view(&value));
    assert_eq!(
        serde_json::to_value(expected.extreme_rays).unwrap(),
        serde_json::to_value(actual.extreme_rays).unwrap()
    );
}
#[test]
fn zero_cone_and_empty_network_are_distinguished() {
    for col in [vec![1, 1], vec![1, 0]] {
        let cone = analyze_cone(&view(&model(
            &["A", "B"],
            &[2, 3],
            std::slice::from_ref(&col),
        )));
        assert_eq!(cone.status, ConeStatus::Exact);
        assert_eq!(cone.extreme_rays.len(), usize::from(col[1] == 0));
        assert_eq!(cone.initial_class_upper_bounds[0], None);
    }
    let cone = analyze_cone(&view(&model(&["A", "B"], &[2, 3], &idle(2))));
    assert_eq!(cone.extreme_rays.len(), 2);
    assert_eq!(
        cone.initial_class_upper_bounds,
        [Some("2".into()), Some("3".into())]
    );
}
#[test]
fn cone_ray_cap_counts_intermediate_rays_without_partial_results() {
    let network = view(&model(
        &["A", "B", "C", "D", "E", "F"],
        &[1; 6],
        &[vec![1, 1, 1, -1, -1, -1]],
    ));
    let exact = analyze_cone_with_limits(
        &network,
        ConeLimits {
            max_rays: 9,
            max_work: MAX_CONE_WORK,
        },
    );
    assert_eq!(exact.status, ConeStatus::Exact);
    assert_eq!(exact.extreme_rays.len(), 9);
    let limited = analyze_cone_with_limits(
        &network,
        ConeLimits {
            max_rays: 8,
            max_work: MAX_CONE_WORK,
        },
    );
    assert_eq!(limited.status, ConeStatus::RayLimit);
    assert!(limited.extreme_rays.is_empty());
    assert!(
        limited
            .initial_class_upper_bounds
            .iter()
            .all(Option::is_none)
    );
    assert_eq!(limited.diagnostics[0].code, "CRN-W013");
}
#[test]
fn cone_work_cap_has_exact_boundary_and_is_clamped() {
    let network = view(&adjacency_fixture());
    let exact = analyze_cone(&network);
    let at = analyze_cone_with_limits(
        &network,
        ConeLimits {
            max_rays: MAX_CONE_RAYS,
            max_work: exact.work_units,
        },
    );
    assert_eq!(at.status, ConeStatus::Exact);
    let below = analyze_cone_with_limits(
        &network,
        ConeLimits {
            max_rays: MAX_CONE_RAYS,
            max_work: exact.work_units - 1,
        },
    );
    assert_eq!(below.status, ConeStatus::WorkLimit);
    assert!(below.extreme_rays.is_empty());
    assert_eq!(below.work_units, exact.work_units - 1);
    let clamped = analyze_cone_with_limits(
        &network,
        ConeLimits {
            max_rays: usize::MAX,
            max_work: u64::MAX,
        },
    );
    assert_eq!(clamped.limits.max_rays, MAX_CONE_RAYS);
    assert_eq!(clamped.limits.max_work, MAX_CONE_WORK);
}
#[test]
fn cone_size_caps_include_both_boundaries() {
    let names: Vec<_> = (0..65).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    for n in [64, 65] {
        let cone = analyze_cone(&view(&model(&refs[..n], &vec![0; n], &idle(n))));
        assert_eq!(
            cone.status,
            if n == 64 {
                ConeStatus::Exact
            } else {
                ConeStatus::SizeLimit
            }
        );
    }
    for n in [256, 257] {
        let cone = analyze_cone(&view(&model(&["A"], &[0], &vec![vec![0]; n])));
        assert_eq!(
            cone.status,
            if n == 256 {
                ConeStatus::Exact
            } else {
                ConeStatus::SizeLimit
            }
        );
    }
}
#[test]
fn cone_arithmetic_and_total_overflow_are_explicit() {
    // Seven species overflow the rays; five overflow only the conserved
    // totals, with counts of 2^53 (0.6.0 used i64::MAX, above the rules
    // count limit).
    for (n, count) in [(7, 0), (5, 1_u64 << 53)] {
        let names: Vec<_> = (0..n).map(|i| format!("X{i}")).collect();
        let refs: Vec<_> = names.iter().map(String::as_str).collect();
        let cone = analyze_cone(&view(&model(&refs, &vec![count; n], &large_chain(n))));
        assert_eq!(cone.status, ConeStatus::ArithmeticOverflow);
        assert!(cone.extreme_rays.is_empty());
        assert_eq!(cone.diagnostics[0].code, "CRN-W014");
    }
}
#[test]
fn supplied_state_does_not_change_initial_class_bounds() {
    let report = analyze_bytes(
        &bytes(&counterexample()),
        Some(br#"{"schema":"crn-craft.state/0.1","counts":{"X0":99,"X1":99,"X2":99,"X3":99}}"#),
    );
    assert!(report.success, "{:?}", report.diagnostics);
    let data = report.analysis.unwrap();
    assert_eq!(data.evaluation_counts, [99; 4]);
    assert_eq!(data.initial_class_upper_bounds[0], Some("3".into()));
}
#[test]
fn linkage_components_direction_and_terminal_classes() {
    let network = view(&model(
        &["A", "B", "C", "D", "E"],
        &[0; 5],
        &[
            vec![-1, 1, 0, 0, 0],
            vec![1, -1, 0, 0, 0],
            vec![0, -1, 1, 0, 0],
            vec![0, 0, -1, 1, 0],
            vec![0, 0, 1, -1, 0],
        ],
    ));
    let graph = analyze_linkage(&network, Some(3));
    assert_eq!(graph.linkage_classes, vec![vec![0, 1, 2, 3]]);
    assert_eq!(graph.strong_linkage_classes, vec![vec![0, 1], vec![2, 3]]);
    assert_eq!(graph.terminal_strong_linkage_classes, vec![vec![2, 3]]);
    assert!(!graph.weakly_reversible);
    assert_eq!(graph.deficiency, Some(0));
}
#[test]
fn linkage_keeps_catalysts_stoichiometry_empty_complex_and_zero_channels() {
    let value = with_rules(
        &["A", "B", "C"],
        &[2, 0, 1],
        vec![
            rule("r0", &[("A", 1)], &[("B", 1)], "1"),
            rule("r1", &[("A", 1), ("C", 1)], &[("B", 1), ("C", 1)], "0"),
            rule("r2", &[], &[("A", 1)], "1"),
            rule("r3", &[("A", 1)], &[], "1"),
        ],
    );
    let graph = analyze_linkage(&view(&value), Some(2));
    assert_eq!(graph.complexes.len(), 5);
    assert_eq!(graph.linkage_classes, vec![vec![0, 1, 4], vec![2, 3]]);
    assert!(graph.complexes[4].is_empty());
    assert_eq!(graph.complexes[2].len(), 2);
}
#[test]
fn weak_reversibility_null_and_positive_deficiency() {
    // 0.6.0 also checked a network without reactions (no complexes); a
    // rules model has at least one rule.
    let mut value = with_rules(
        &["A"],
        &[1],
        vec![
            rule("r0", &[("A", 1)], &[("A", 2)], "1"),
            rule("r1", &[("A", 2)], &[("A", 3)], "1"),
        ],
    );
    assert_eq!(analyze_linkage(&view(&value), Some(1)).deficiency, Some(1));
    assert_eq!(analyze_linkage(&view(&value), None).deficiency, None);
    value["rules"][1] = rule("r1", &[("A", 2)], &[("A", 1)], "1");
    assert!(analyze_linkage(&view(&value), Some(1)).weakly_reversible);
    value["rules"] = json!([rule("null", &[("A", 1)], &[("A", 1)], "1")]);
    let null = analyze_linkage(&view(&value), Some(0));
    assert_eq!(null.complexes.len(), 1);
    assert!(null.weakly_reversible);
    assert_eq!(null.deficiency, Some(0));
}
#[test]
fn linkage_walk_is_iterative_for_long_chain() {
    let names: Vec<_> = (0..4096).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let rules: Vec<Value> = (0..4095)
        .map(|i| rule(&format!("r{i}"), &[(refs[i], 1)], &[(refs[i + 1], 1)], "1"))
        .collect();
    let graph = analyze_linkage(&view(&with_rules(&refs, &vec![0; 4096], rules)), None);
    assert_eq!(graph.linkage_classes.len(), 1);
    assert_eq!(graph.strong_linkage_classes.len(), 4096);
    assert_eq!(graph.terminal_strong_linkage_classes, vec![vec![4095]]);
}
#[test]
fn cone_failure_preserves_exact_linear_analysis_and_fallback_bound() {
    // The final source equation kills the large chain's conservation space.
    // Linear elimination stays small, while DD exhausts i128 before that equation.
    let names: Vec<_> = (0..8).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let mut columns: Vec<_> = large_chain(7)
        .into_iter()
        .map(|mut column| {
            column.push(0);
            column
        })
        .collect();
    let mut source = vec![0; 8];
    source[6] = 1;
    columns.push(source);
    let mut initial = vec![0; 8];
    initial[7] = 7;
    let report = analyze_bytes(&bytes(&model(&refs, &initial, &columns)), None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert!(report.diagnostics.iter().any(|d| d.code == "CRN-W014"));
    let data = report.analysis.unwrap();
    assert_eq!(data.conservation.status, ConservationStatus::Exact);
    assert_eq!(data.conservation.rank, Some(7));
    assert_eq!(
        data.nonnegative_conservation.status,
        ConeStatus::ArithmeticOverflow
    );
    assert!(data.nonnegative_conservation.extreme_rays.is_empty());
    assert_eq!(data.initial_class_upper_bounds[7], Some("7".into()));
    assert!(
        data.initial_class_upper_bounds[..7]
            .iter()
            .all(Option::is_none)
    );
}
#[test]
fn combined_rays_are_reduced_to_primitive_coefficients() {
    // Combining the two rays of the last constraint gives 2 * (9, 23, 9, 6);
    // the report must carry the primitive law, as its contract states.
    let value = model(
        &["X0", "X1", "X2", "X3"],
        &[0, 6, 6, 9],
        &[vec![-2, 3, -3, -4], vec![-4, 0, 2, 3], vec![3, 0, -3, 0]],
    );
    let cone = analyze_cone(&view(&value));
    assert_eq!(cone.status, ConeStatus::Exact);
    assert_eq!(cone.extreme_rays.len(), 1);
    assert_eq!(cone.extreme_rays[0].coefficients, ["9", "23", "9", "6"]);
    assert_eq!(cone.extreme_rays[0].initial_total, "246");
    assert_eq!(
        cone.initial_class_upper_bounds,
        [
            Some("27".into()),
            Some("10".into()),
            Some("27".into()),
            Some("41".into())
        ]
    );
}
#[test]
fn initial_orthant_rays_respect_a_lower_ray_cap() {
    // Without net changes the n unit rays are the whole answer; a cap below
    // n must still report ray_limit rather than an exact cone with n rays.
    let network = view(&model(&["A", "B", "C"], &[1, 2, 3], &idle(3)));
    let limits = |max_rays| ConeLimits {
        max_rays,
        max_work: MAX_CONE_WORK,
    };
    let at = analyze_cone_with_limits(&network, limits(3));
    assert_eq!(at.status, ConeStatus::Exact);
    assert_eq!(at.extreme_rays.len(), 3);
    let below = analyze_cone_with_limits(&network, limits(2));
    assert_eq!(below.status, ConeStatus::RayLimit);
    assert!(below.extreme_rays.is_empty());
    // A view without channels (0.6.0: a network without reactions): the
    // only rule matches no tuple, so the enumeration loop never runs and
    // the unit rays alone must respect the cap.
    let model: Model = serde_json::from_value(json!({
        "schema": "crn-craft.model/1.0", "name": "idle",
        "types": [{"name": "M", "attributes": ["wt"]}],
        "initial": [
            {"type": "M", "attributes": {"wt": 0}, "count": 1},
            {"type": "M", "attributes": {"wt": 1}, "count": 2},
            {"type": "M", "attributes": {"wt": 2}, "count": 3}
        ],
        "rules": [{"name": "never", "reactants": [{"label": "x", "type": "M", "where": "x.wt == 5"}],
                   "rate": "1", "products": []}]
    }))
    .unwrap();
    let network = FlatView::with_every_type(&model).unwrap();
    assert_eq!((network.species().len(), network.channels()), (3, 0));
    let at = analyze_cone_with_limits(&network, limits(3));
    assert_eq!(at.status, ConeStatus::Exact);
    assert_eq!(
        at.initial_class_upper_bounds,
        [Some("1".into()), Some("2".into()), Some("3".into())]
    );
    let below = analyze_cone_with_limits(&network, limits(2));
    assert_eq!(below.status, ConeStatus::RayLimit);
    assert!(below.extreme_rays.is_empty());
}
#[test]
fn linkage_classes_join_converging_reactions() {
    // A -> B and C -> B: complex C reaches the others only through an edge
    // into B, so a walk over outgoing edges alone would split the class.
    let graph = analyze_linkage(
        &view(&model(
            &["A", "B", "C"],
            &[0; 3],
            &[vec![-1, 1, 0], vec![0, 1, -1]],
        )),
        Some(2),
    );
    assert_eq!(graph.linkage_classes, vec![vec![0, 1, 2]]);
    assert_eq!(
        graph.strong_linkage_classes,
        vec![vec![0], vec![1], vec![2]]
    );
    assert_eq!(graph.terminal_strong_linkage_classes, vec![vec![1]]);
    assert!(!graph.weakly_reversible);
    assert_eq!(graph.deficiency, Some(0));
}

// Conservation tests ported from tests/network.rs.

#[test]
fn known_enzyme_laws_and_bounds_are_found() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../models/legacy/enzyme-mechanism.model.json");
    let bytes = std::fs::read(path).unwrap();
    let c = analyze_conservation(&FlatView::of_any_model(&bytes).unwrap());
    assert_eq!(c.status, ConservationStatus::Exact);
    assert_eq!(c.rank, Some(2));
    let laws: Vec<Vec<i128>> = c
        .basis
        .iter()
        .map(|l| l.coefficients.iter().map(|s| s.parse().unwrap()).collect())
        .collect();
    assert!(laws.contains(&vec![0, 1, 1, 0]));
    assert!(laws.contains(&vec![1, 0, 1, 1]));
    assert_eq!(
        c.initial_class_upper_bounds,
        [
            Some("50".into()),
            Some("10".into()),
            Some("10".into()),
            Some("50".into())
        ]
    );
    let report = analyze_bytes(&bytes, None);
    assert!(report.success, "{:?}", report.diagnostics);
    assert!(!report.diagnostics.iter().any(|d| d.code == "CRN-W001"));
}

#[test]
fn signed_law_does_not_prove_individual_count_bounds() {
    let c = analyze_conservation(&view(&model(&["A", "B"], &[0, 0], &[vec![1, 1]])));
    assert_eq!(c.status, ConservationStatus::Exact);
    assert_eq!(c.basis[0].coefficients, ["1", "-1"]);
    assert_eq!(c.initial_class_upper_bounds, [None, None]);
}

#[test]
fn empty_network_has_complete_identity_conservation_basis() {
    let c = analyze_conservation(&view(&model(&["A", "B"], &[7, 9], &idle(2))));
    assert_eq!(c.rank, Some(0));
    assert_eq!(c.basis.len(), 2);
    assert_eq!(
        c.initial_class_upper_bounds,
        [Some("7".into()), Some("9".into())]
    );
}

#[test]
fn full_rank_network_has_no_conservation_basis() {
    let c = analyze_conservation(&view(&model(&["A"], &[0], &[vec![1]])));
    assert_eq!(c.rank, Some(1));
    assert!(c.basis.is_empty());
}

#[test]
fn conservation_size_limit_is_honest() {
    let names: Vec<_> = (0..65).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let at_limit = analyze_conservation(&view(&model(&refs[..64], &[0; 64], &idle(64))));
    assert_eq!(at_limit.status, ConservationStatus::Exact);
    assert_eq!(at_limit.basis.len(), 64);
    let c = analyze_conservation(&view(&model(&refs, &[0; 65], &idle(65))));
    assert_eq!(c.status, ConservationStatus::SizeLimit);
    assert_eq!(c.rank, None);
    assert!(c.basis.is_empty());
    assert_eq!(c.diagnostics[0].code, "CRN-W010");
    for count in [256, 257] {
        let report = analyze_bytes(&bytes(&model(&["A"], &[1], &vec![vec![0]; count])), None);
        assert!(report.success, "{:?}", report.diagnostics);
        let conservation = report.analysis.unwrap().conservation;
        if count == 256 {
            assert_eq!(conservation.status, ConservationStatus::Exact);
            assert_eq!(conservation.rank, Some(0));
        } else {
            assert_eq!(conservation.status, ConservationStatus::SizeLimit);
            assert_eq!(conservation.rank, None);
            assert!(conservation.basis.is_empty());
            assert!(report.diagnostics.iter().any(|d| d.code == "CRN-W010"));
        }
    }
}

#[test]
fn conservation_overflow_never_reports_partial_basis() {
    let names: Vec<_> = (0..7).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let c = analyze_conservation(&view(&model(&refs, &[0; 7], &large_chain(7))));
    assert_eq!(c.status, ConservationStatus::ArithmeticOverflow);
    assert!(c.basis.is_empty());
    assert_eq!(c.rank, None);
    assert_eq!(c.diagnostics[0].code, "CRN-W011");
}

#[test]
fn analysis_failure_and_state_identity_are_honest() {
    let bytes = bytes(&model(&["A"], &[1], &idle(1)));
    let state = br#"{"schema":"crn-craft.state/0.1","counts":{"A":9}}"#;
    let report = analyze_bytes(&bytes, Some(state));
    assert!(report.success, "{:?}", report.diagnostics);
    assert!(!report.simulation_available);
    assert_eq!(report.state_sha256, Some(crn_craft::sha256(state)));
    let data = report.analysis.unwrap();
    assert_eq!(data.evaluation_counts, [9]);
    assert_eq!(data.initial_counts, [1]);
    assert_eq!(
        data.conservation.initial_class_upper_bounds,
        [Some("1".into())]
    );
    let bad = analyze_bytes(&bytes, Some(b"{}"));
    assert!(!bad.success);
    assert!(bad.analysis.is_none());
}

#[test]
fn weighted_nonnegative_law_divides_its_bound_by_the_weight() {
    // 2A <-> B conserves A + 2B = 7, so B <= floor(7 / 2).
    let c = analyze_conservation(&view(&model(
        &["A", "B"],
        &[5, 1],
        &[vec![-2, 1], vec![2, -1]],
    )));
    assert_eq!(c.status, ConservationStatus::Exact);
    assert_eq!(c.basis.len(), 1);
    assert_eq!(c.basis[0].coefficients, ["1", "2"]);
    assert_eq!(c.basis[0].initial_total, "7");
    assert_eq!(
        c.initial_class_upper_bounds,
        [Some("7".into()), Some("3".into())]
    );
}
