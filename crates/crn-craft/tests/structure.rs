use crn_craft::{
    CompiledNetwork, analysis_v1, analyze_bytes,
    cone::{
        ConeLimits, ConeStatus, MAX_CONE_RAYS, MAX_CONE_WORK, analyze_cone,
        analyze_cone_with_limits,
    },
    conservation::analyze_conservation,
    export::export_gillespy2,
    linkage::analyze_linkage,
    parse_model,
};
use serde_json::{Value, json};
fn model(names: &[&str], counts: &[i64], columns: &[Vec<i64>]) -> Value {
    let reactions: Vec<_> = columns.iter().enumerate().map(|(r, col)| {
        let mut left = serde_json::Map::new(); let mut right = serde_json::Map::new();
        for (name, c) in names.iter().zip(col) {
            if *c < 0 { left.insert((*name).into(), json!(-c)); }
            if *c > 0 { right.insert((*name).into(), json!(c)); }
        }
        json!({"name":format!("r{r}"),"reactants":left,"products":right,"rate":{"kind":"mass_action","coefficient":1.0}})
    }).collect();
    json!({"schema":"crn-craft.model/0.1","name":"structure","convention":"stochastic_counts","omega":1.0,"parameters":{},
        "species":names.iter().zip(counts).map(|(name,count)|json!({"name":name,"initial_count":count})).collect::<Vec<_>>(),"reactions":reactions})
}
fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}
fn compile(value: &Value) -> CompiledNetwork {
    CompiledNetwork::compile(&parse_model(&bytes(value)).unwrap()).unwrap()
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

#[test]
fn review_counterexample_gets_complete_bounds_and_suppresses_growth_warning() {
    let value = counterexample();
    let network = compile(&value);
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
    assert!(report.success);
    assert!(!report.diagnostics.iter().any(|d| d.code == "CRN-W001"));
    assert!(
        analysis_v1::analyze_bytes(&bytes(&value), None)
            .diagnostics
            .iter()
            .any(|d| d.code == "CRN-W001")
    );
}
#[test]
fn adjacency_excludes_nonextreme_combination() {
    let cone = analyze_cone(&compile(&adjacency_fixture()));
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
    let cone = analyze_cone(&compile(&model(
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
    let expected = analyze_cone(&compile(&value));
    let reactions = value["reactions"].as_array_mut().unwrap();
    let mut duplicate = reactions[0].clone();
    duplicate["name"] = json!("duplicate");
    reactions.push(duplicate);
    let mut reverse = reactions[0].clone();
    reverse["name"] = json!("reverse");
    let temp = reverse["products"].clone();
    reverse["products"] = reverse["reactants"].clone();
    reverse["reactants"] = temp;
    reactions.push(reverse);
    let actual = analyze_cone(&compile(&value));
    assert_eq!(
        serde_json::to_value(expected.extreme_rays).unwrap(),
        serde_json::to_value(actual.extreme_rays).unwrap()
    );
}
#[test]
fn zero_cone_and_empty_network_are_distinguished() {
    for col in [vec![1, 1], vec![1, 0]] {
        let cone = analyze_cone(&compile(&model(
            &["A", "B"],
            &[2, 3],
            std::slice::from_ref(&col),
        )));
        assert_eq!(cone.status, ConeStatus::Exact);
        assert_eq!(cone.extreme_rays.len(), usize::from(col[1] == 0));
        assert_eq!(cone.initial_class_upper_bounds[0], None);
    }
    let cone = analyze_cone(&compile(&model(&["A", "B"], &[2, 3], &[])));
    assert_eq!(cone.extreme_rays.len(), 2);
    assert_eq!(
        cone.initial_class_upper_bounds,
        [Some("2".into()), Some("3".into())]
    );
}
#[test]
fn cone_ray_cap_counts_intermediate_rays_without_partial_results() {
    let network = compile(&model(
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
    let network = compile(&adjacency_fixture());
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
        let cone = analyze_cone(&compile(&model(&refs[..n], &vec![0; n], &[])));
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
        let cone = analyze_cone(&compile(&model(&["A"], &[0], &vec![vec![0]; n])));
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
    for (n, count) in [(7, 0), (5, i64::MAX)] {
        let names: Vec<_> = (0..n).map(|i| format!("X{i}")).collect();
        let refs: Vec<_> = names.iter().map(String::as_str).collect();
        let columns: Vec<_> = (0..n - 1)
            .map(|i| {
                let mut row = vec![0; n];
                row[i] = -1;
                row[i + 1] = 2147483647;
                row
            })
            .collect();
        let cone = analyze_cone(&compile(&model(&refs, &vec![count; n], &columns)));
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
    assert!(report.success);
    let data = report.analysis.unwrap();
    assert_eq!(data.evaluation_counts, [99; 4]);
    assert_eq!(data.initial_class_upper_bounds[0], Some("3".into()));
}
#[test]
fn linkage_components_direction_and_terminal_classes() {
    let network = compile(&model(
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
    let mut value = model(
        &["A", "B", "C"],
        &[2, 0, 1],
        &[
            vec![-1, 1, 0],
            vec![-1, 1, 0],
            vec![1, 0, 0],
            vec![-1, 0, 0],
        ],
    );
    value["reactions"][1]["reactants"]["C"] = json!(1);
    value["reactions"][1]["products"]["C"] = json!(1);
    value["reactions"][1]["rate"]["coefficient"] = json!(0);
    let graph = analyze_linkage(&compile(&value), Some(2));
    assert_eq!(graph.complexes.len(), 5);
    assert_eq!(graph.linkage_classes, vec![vec![0, 1, 4], vec![2, 3]]);
    assert!(graph.complexes[4].is_empty());
    assert_eq!(graph.complexes[2].len(), 2);
}
#[test]
fn weak_reversibility_empty_null_and_positive_deficiency() {
    let empty = analyze_linkage(&compile(&model(&["A"], &[1], &[])), Some(0));
    assert!(empty.weakly_reversible);
    assert_eq!(empty.deficiency, Some(0));
    assert!(empty.complexes.is_empty());
    let mut value = model(&["A"], &[1], &[vec![1], vec![1]]);
    value["reactions"][0]["reactants"] = json!({"A":1});
    value["reactions"][0]["products"] = json!({"A":2});
    value["reactions"][1]["reactants"] = json!({"A":2});
    value["reactions"][1]["products"] = json!({"A":3});
    assert_eq!(
        analyze_linkage(&compile(&value), Some(1)).deficiency,
        Some(1)
    );
    assert_eq!(analyze_linkage(&compile(&value), None).deficiency, None);
    value["reactions"][1]["products"] = json!({"A":1});
    assert!(analyze_linkage(&compile(&value), Some(1)).weakly_reversible);
    value["reactions"] = json!([{"name":"null","reactants":{"A":1},"products":{"A":1},"rate":{"kind":"mass_action","coefficient":1}}]);
    let null = analyze_linkage(&compile(&value), Some(0));
    assert_eq!(null.complexes.len(), 1);
    assert!(null.weakly_reversible);
    assert_eq!(null.deficiency, Some(0));
}
#[test]
fn linkage_walk_is_iterative_for_long_chain() {
    let names: Vec<_> = (0..4096).map(|i| format!("X{i}")).collect();
    let species: Vec<_> = names
        .iter()
        .map(|n| json!({"name":n,"initial_count":0}))
        .collect();
    let reactions:Vec<_>=(0..4095).map(|i|json!({"name":format!("r{i}"),"reactants":{names[i].clone():1},"products":{names[i+1].clone():1},"rate":{"kind":"mass_action","coefficient":1}})).collect();
    let mut value = model(&["A"], &[0], &[]);
    value["species"] = json!(species);
    value["reactions"] = json!(reactions);
    let graph = analyze_linkage(&compile(&value), None);
    assert_eq!(graph.linkage_classes.len(), 1);
    assert_eq!(graph.strong_linkage_classes.len(), 4096);
    assert_eq!(graph.terminal_strong_linkage_classes, vec![vec![4095]]);
}
#[test]
fn export_halves_homodimer_and_uses_safe_aliases_unit_volume() {
    let mut value = model(&["vol", "K0"], &[3, 0], &[vec![-2, 1]]);
    value["omega"] = json!(1000);
    let report = export_gillespy2(&bytes(&value));
    assert!(report.success);
    assert_eq!(report.model_sha256, Some(crn_craft::sha256(&bytes(&value))));
    let artifact = report.artifact.unwrap();
    assert_eq!(artifact.target_version, "1.8.3");
    assert!(artifact.source.contains("expression=0.5"));
    assert!(artifact.source.contains("volume=1.0"));
    assert_eq!(artifact.species_aliases["vol"], "S0");
    assert_eq!(artifact.species_aliases["K0"], "S1");
}
#[test]
fn export_rejects_unsupported_semantics_before_emitting_code() {
    let base = model(&["A"], &[2], &[vec![-1]]);
    for mutate in [
        |v: &mut Value| v["species"][0]["population_limit"] = json!(10),
        |v: &mut Value| v["species"][0]["initial_count"] = json!(9_007_199_254_740_993_i64),
        |v: &mut Value| v["reactions"][0]["reactants"]["A"] = json!(3),
        |v: &mut Value| {
            v["reactions"][0]["rate"] =
                json!({"kind":"michaelis_menten","substrate":"A","vmax":1,"km":1})
        },
    ] {
        let mut value = base.clone();
        mutate(&mut value);
        let report = export_gillespy2(&bytes(&value));
        assert!(!report.success);
        assert!(report.artifact.is_none());
        assert!(report.diagnostics.iter().any(|d| d.code == "CRN-E050"));
    }
}
#[test]
fn export_rejects_inexact_subnormal_halving_but_accepts_exact_boundary() {
    let mut value = model(&["A"], &[2], &[vec![-2]]);
    for bits in [1, 2, 3, 4] {
        value["reactions"][0]["rate"]["coefficient"] = json!(f64::from_bits(bits));
        let report = export_gillespy2(&bytes(&value));
        assert_eq!(report.success, bits % 2 == 0);
        if bits % 2 != 0 {
            assert_eq!(report.diagnostics[0].code, "CRN-E051");
            assert!(report.artifact.is_none());
        }
    }
}
#[test]
fn export_validates_typed_contract_and_exact_float_count_boundary() {
    let mut value = model(&["A"], &[9_007_199_254_740_992], &[]);
    assert!(export_gillespy2(&bytes(&value)).success);
    value["unexpected"] = json!(true);
    assert_eq!(
        export_gillespy2(&bytes(&value)).diagnostics[0].code,
        "CRN-E003"
    );
    assert!(!export_gillespy2(b"{}").success);
}

#[test]
fn cone_failure_preserves_exact_linear_analysis_and_fallback_bound() {
    // The final source equation kills the large chain's conservation space.
    // Linear elimination stays small, while DD exhausts i128 before that equation.
    let names: Vec<_> = (0..8).map(|i| format!("X{i}")).collect();
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let mut columns: Vec<_> = (0..6)
        .map(|i| {
            let mut row = vec![0; 8];
            row[i] = -1;
            row[i + 1] = 2147483647;
            row
        })
        .collect();
    let mut source = vec![0; 8];
    source[6] = 1;
    columns.push(source);
    let mut initial = vec![0; 8];
    initial[7] = 7;
    let report = analyze_bytes(&bytes(&model(&refs, &initial, &columns)), None);
    assert!(report.success);
    assert!(report.diagnostics.iter().any(|d| d.code == "CRN-W014"));
    let data = report.analysis.unwrap();
    assert_eq!(
        data.conservation.status,
        crn_craft::conservation::ConservationStatus::Exact
    );
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
    let cone = analyze_cone(&compile(&value));
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
    // Without reactions the n unit rays are the whole answer; a cap below n
    // must still report ray_limit rather than an exact cone with n rays.
    let network = compile(&model(&["A", "B", "C"], &[1, 2, 3], &[]));
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
}
#[test]
fn linkage_classes_join_converging_reactions() {
    // A -> B and C -> B: complex C reaches the others only through an edge
    // into B, so a walk over outgoing edges alone would split the class.
    let graph = analyze_linkage(
        &compile(&model(
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
