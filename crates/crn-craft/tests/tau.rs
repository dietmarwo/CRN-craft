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
fn request(model: &Value, replicates: u64, end: f64, epsilon: f64) -> (Vec<u8>, Vec<u8>) {
    let m = serde_json::to_vec(model).unwrap();
    let run = json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":crn_craft::sha256(&m)},"method":{"kind":"ssa","aggregator":"direct"},"ensemble":{"replicates":replicates,"root_seed":5},"time":{"end":end,"save_times":[0.,end]},"outputs":{"statistics":["mean"],"keep_trajectories":replicates},"budget":{"threads":1,"max_events_per_replicate":100000000u64,"max_output_bytes":100000000u64}});
    let r = json!({"schema":"crn-craft.tau-request/0.1","epsilon":epsilon,"max_leaps_per_replicate":1000000,"reference_run":run});
    (m, serde_json::to_vec(&r).unwrap())
}
fn birth_death(initial: i64, births: f64, deaths: usize) -> Value {
    let mut reactions = vec![
        json!({"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":births}}),
    ];
    for i in 0..deaths {
        reactions.push(json!({"name":format!("death{i}"),"reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":1.0}}));
    }
    json!({"schema":"crn-craft.model/0.1","name":"Leap","convention":"stochastic_counts","omega":1.,"parameters":{},"species":[{"name":"X","initial_count":initial}],"reactions":reactions})
}
#[test]
fn critical_counts_take_direct_steps_and_rejected_leaps_are_halved() {
    // Five molecules are critical at every step: tau runs Direct events only.
    let (m, r) = request(&birth_death(5, 0.0, 1), 50, 10.0, 0.05);
    let a = tau::solve(&m, &r, &Control::default()).unwrap();
    assert_eq!(a.termination, Stop::Completed);
    assert_eq!(a.events, 250);
    assert_eq!(a.leaps, a.events);
    assert_eq!(a.rejected_leaps, 0);
    // Just above the critical count and with three consumers, a few Poisson
    // proposals overdraw and must be rejected and halved, never clipped.
    let (m, r) = request(&birth_death(15, 45.0, 3), 500, 20.0, 0.2);
    let a = tau::solve(&m, &r, &Control::default()).unwrap();
    assert_eq!(a.termination, Stop::Completed);
    assert!(a.rejected_leaps > 0, "no rejected leap");
    assert!(
        a.leaps < a.events / 2,
        "{} leaps for {} events",
        a.leaps,
        a.events
    );
    assert!(
        a.trajectories
            .iter()
            .all(|t| t.counts.iter().flatten().all(|x| *x >= 0))
    );
}
#[test]
fn large_populations_complete_with_capped_leap_means() {
    // At equilibrium the drift vanishes and the variance bound allows dt = 5000,
    // a leap mean of 5e9: the step is capped so each mean stays within 1e7.
    let (m, r) = request(&birth_death(1_000_000, 1e6, 1), 1, 20.0, 0.1);
    let a = tau::solve(&m, &r, &Control::default()).unwrap();
    assert_eq!(
        a.termination,
        Stop::Completed,
        "{:?}",
        a.replicates[0].diagnostic
    );
    assert!(a.leaps >= 2);
    assert!(a.events > 10_000_000);
    let x = a.trajectories[0].counts[1][0] as f64;
    assert!((x - 1e6).abs() < 6.0 * 1e3, "{x}");
}
#[test]
fn tau_and_direct_with_one_seed_use_different_random_words() {
    let (m, r) = request(&birth_death(30, 20.0, 1), 400, 1.0, 0.05);
    let tau = tau::solve(&m, &r, &Control::default()).unwrap();
    let request: Value = serde_json::from_slice(&r).unwrap();
    let run = serde_json::to_vec(&request["reference_run"]).unwrap();
    let exact = crn_craft::simulation::simulate_bytes(&m, &run, &Control::default())
        .unwrap()
        .result;
    let same = tau
        .trajectories
        .iter()
        .zip(&exact.trajectories)
        .filter(|(a, b)| a.counts[1] == b.counts[1])
        .count();
    // Independent final counts agree by chance about 6% of the time.
    assert!(same < 80, "{same} of 400 final counts coincide");
    assert_eq!(tau.algorithm, tau::ALGORITHM);
}
