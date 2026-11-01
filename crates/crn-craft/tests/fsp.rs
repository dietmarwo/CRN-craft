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
