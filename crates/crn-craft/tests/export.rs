//! The GillesPy2 reference export (moved from tests/structure.rs in M5b).
use crn_craft::export::export_gillespy2;
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
