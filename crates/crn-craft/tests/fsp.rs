use crn_craft::{fsp, simulation::Control};
use serde_json::{Value, json};
fn fixture(bound: i64, k: f64, end: f64) -> (Vec<u8>, Vec<u8>) {
    let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"Decay","convention":"stochastic_counts","omega":1.0,"parameters":{},"species":[{"name":"A","initial_count":4}],"reactions":[{"name":"death","reactants":{"A":1},"products":{},"rate":{"kind":"mass_action","coefficient":k}}]})).unwrap();
    let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(&model)},"method":{"kind":"fsp","bounds":{"A":bound},"max_states":100,"tolerance":1e-14},"time":{"end":end,"save_times":[0.0,end]},"outputs":{"statistics":["mean"],"keep_trajectories":0},"budget":{"threads":1,"max_output_bytes":1000000}})).unwrap();
    (model, run)
}
#[test]
fn binomial_decay_and_large_lambda_t() {
    for end in [0.5, 1000.0] {
        let (m, r) = fixture(4, 2.0, end);
        let f = fsp::solve(&m, &r, &Control::default()).unwrap();
        let q = (-2.0 * end).exp();
        for (state, &p) in f.states.iter().zip(&f.snapshots[1].probabilities) {
            let k = state[0] as i32;
            let binomial = [1.0, 4.0, 6.0, 4.0, 1.0][k as usize];
            let expected = binomial * q.powi(k) * (1.0 - q).powi(4 - k);
            assert!((p - expected).abs() < 1e-12, "{state:?} {p} {expected}");
        }
        assert_eq!(f.snapshots[1].physical_leakage_estimate, 0.0);
        assert!(f.snapshots[1].series_tail_upper_bound <= 1.01e-14);
        assert!(f.snapshots[1].series_tail_upper_bound > 0.0);
        assert!(!f.numerical_error_certified);
    }
}
#[test]
fn leakage_is_not_reflected_or_double_counted() {
    let (m, r) = fixture(4, 2.0, 1.0);
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["reactions"][0]["reactants"] = json!({});
    model["reactions"][0]["products"] = json!({"A":1});
    let m = serde_json::to_vec(&model).unwrap();
    let mut run: Value = serde_json::from_slice(&r).unwrap();
    run["model"]["sha256"] = json!(crn_craft::sha256(&m));
    let f = fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).unwrap();
    assert_eq!(f.states.len(), 1);
    let s = &f.snapshots[1];
    assert!((s.probabilities[0] - (-2.0_f64).exp()).abs() < 1e-12);
    assert!((s.physical_leakage_estimate - (1.0 - (-2.0_f64).exp())).abs() < 1e-12);
    assert!(s.mass_defect.abs() < 1e-12);
}
#[test]
fn capacity_cancellation_bounds_and_null_reaction() {
    let (m, r) = fixture(4, 1.0, 1.0);
    let mut run: Value = serde_json::from_slice(&r).unwrap();
    run["method"]["max_states"] = json!(2);
    assert!(fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).is_err());
    let control = Control::default();
    control.cancel();
    assert!(
        fsp::solve(&m, &r, &control).unwrap_err()[0]
            .message
            .contains("interrupted")
    );
    run["method"]["bounds"]["A"] = json!(3);
    assert!(fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).is_err());
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["reactions"][0]["products"] = json!({"A":1});
    let m = serde_json::to_vec(&model).unwrap();
    let mut run: Value = serde_json::from_slice(&r).unwrap();
    run["model"]["sha256"] = json!(crn_craft::sha256(&m));
    let f = fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).unwrap();
    assert_eq!(f.uniformization_rate, 0.0);
    assert_eq!(f.states.len(), 1);
    assert_eq!(f.snapshots[1].probabilities, vec![1.0]);
}
#[test]
fn omitted_weight_is_summed_not_cancelled_and_stays_below_the_tail_bound() {
    // About 1000 split steps: 1 - (kept weight) would accumulate roundoff.
    let (m, r) = fixture(4, 2.0, 1000.0);
    let f = fsp::solve(&m, &r, &Control::default()).unwrap();
    for s in &f.snapshots {
        assert!(s.discarded_weight_estimate >= 0.0);
        assert!(
            s.discarded_weight_estimate <= s.series_tail_upper_bound,
            "{} > {}",
            s.discarded_weight_estimate,
            s.series_tail_upper_bound
        );
        assert!(
            s.roundoff_mass_residual < 1e-13,
            "{}",
            s.roundoff_mass_residual
        );
    }
}
#[test]
fn the_tail_bound_bounds_the_missing_mass_at_a_loose_tolerance() {
    let (m, r) = fixture(4, 2.0, 3.0);
    let mut run: Value = serde_json::from_slice(&r).unwrap();
    run["method"]["tolerance"] = json!(0.5);
    let f = fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).unwrap();
    let last = f.snapshots.last().unwrap();
    assert!(
        last.mass_defect > 1e-6,
        "tolerance 0.5 should lose visible mass"
    );
    assert!(
        last.mass_defect <= last.series_tail_upper_bound + 1e-12,
        "{} > {}",
        last.mass_defect,
        last.series_tail_upper_bound
    );
    assert!(last.discarded_weight_estimate <= last.series_tail_upper_bound);
    assert!((last.discarded_weight_estimate - last.mass_defect).abs() < 1e-12);
}
#[test]
fn many_snapshots_respect_the_output_budget() {
    // One reachable state: the early estimate counted 32 bytes per snapshot,
    // while each snapshot record serializes to about 240 bytes.
    let (m, _) = fixture(4, 2.0, 1.0);
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["species"][0]["initial_count"] = json!(0);
    let m = serde_json::to_vec(&model).unwrap();
    let (_, r) = fixture(4, 2.0, 1.0);
    let mut run: Value = serde_json::from_slice(&r).unwrap();
    run["model"]["sha256"] = json!(crn_craft::sha256(&m));
    let saves: Vec<f64> = (0..20000).map(|i| i as f64 * 1e-4).collect();
    run["time"] = json!({"end": saves[saves.len() - 1], "save_times": saves});
    run["method"]["tolerance"] = json!(1e-10);
    run["budget"]["max_output_bytes"] = json!(2_000_000);
    let errors =
        fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).unwrap_err();
    assert_eq!(errors[0].code, "CRN-E070");
    assert!(
        errors[0].message.contains("output"),
        "{}",
        errors[0].message
    );
    run["budget"]["max_output_bytes"] = json!(100_000_000);
    let f = fsp::solve(&m, &serde_json::to_vec(&run).unwrap(), &Control::default()).unwrap();
    let bytes = serde_json::to_vec(&f).unwrap().len();
    assert!(bytes > 2_000_000, "fixture output is only {bytes} bytes");
}
