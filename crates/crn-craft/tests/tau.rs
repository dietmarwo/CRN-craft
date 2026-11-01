use crn_craft::{
    simulation::{Control, Stop},
    tau,
};
use serde_json::{Value, json};
fn fixture() -> (Vec<u8>, Value) {
    let m=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"Leap","convention":"stochastic_counts","omega":1.,"parameters":{},"species":[{"name":"X","initial_count":30}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":20.}},{"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":0.5}}]})).unwrap();
    let run = json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(&m)},"method":{"kind":"ssa","aggregator":"direct"},"ensemble":{"replicates":64,"root_seed":91},"time":{"end":2.,"save_times":[0.,1.,2.]},"outputs":{"statistics":["mean","variance"],"keep_trajectories":64},"budget":{"threads":1,"max_events_per_replicate":10000,"max_output_bytes":1000000}});
    (
        m,
        json!({"schema":"crn-craft.tau-request/0.1","epsilon":0.05,"max_leaps_per_replicate":10000,"reference_run":run}),
    )
}
#[test]
fn approximate_assessment_nonnegative_repeatability_and_limits() {
    let (m, mut r) = fixture();
    let bytes = serde_json::to_vec(&r).unwrap();
    let a = tau::solve(&m, &bytes, &Control::default()).unwrap();
    let b = tau::solve(&m, &bytes, &Control::default()).unwrap();
    assert_eq!(
        serde_json::to_value(&a).unwrap(),
        serde_json::to_value(&b).unwrap()
    );
    assert_eq!(a.assessment, "approximate_ensemble");
    assert_eq!(a.termination, Stop::Completed);
    assert!(
        a.trajectories
            .iter()
            .flat_map(|t| &t.counts)
            .flatten()
            .all(|x| *x >= 0)
    );
    assert!(a.statistics.is_some());
    let final_mean = a.statistics.as_ref().unwrap().mean.as_ref().unwrap()[2][0];
    assert!((32.0..40.0).contains(&final_mean));
    r["reference_run"]["budget"]["max_events_per_replicate"] = json!(1);
    let a = tau::solve(&m, &serde_json::to_vec(&r).unwrap(), &Control::default()).unwrap();
    assert_eq!(a.termination, Stop::EventLimit);
    assert!(a.replicates.iter().all(|record| record.events <= 1));
    assert!(a.statistics.is_none());
    let c = Control::default();
    c.cancel();
    let a = tau::solve(&m, &bytes, &c).unwrap();
    assert_eq!(a.termination, Stop::Cancelled);
    assert!(a.statistics.is_none());
    r["epsilon"] = json!(0.);
    assert!(tau::solve(&m, &serde_json::to_vec(&r).unwrap(), &Control::default()).is_err());
    assert!(tau::solve(&m, b"{\"schema\":0,\"schema\":1}", &Control::default()).is_err());
}
#[test]
fn critical_homodimer_never_passes_the_last_reachable_count() {
    let (mut m, r) = fixture();
    let mut model: Value = serde_json::from_slice(&m).unwrap();
    model["species"][0]["initial_count"] = json!(3);
    model["reactions"] = json!([{"name":"critical","reactants":{"X":2},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":100.}}]);
    m = serde_json::to_vec(&model).unwrap();
    let mut r = r;
    r["reference_run"]["model"]["sha256"] = json!(crn_craft::sha256(&m));
    let a = tau::solve(&m, &serde_json::to_vec(&r).unwrap(), &Control::default()).unwrap();
    assert_eq!(a.termination, Stop::Completed);
    assert!(
        a.trajectories
            .iter()
            .all(|t| t.counts.last().unwrap()[0] == 1)
    );
}
