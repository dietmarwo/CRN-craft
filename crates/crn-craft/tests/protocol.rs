//! Scheduled interventions and fixed-delay completions (`crn-craft.run/1.1`).
//!
//! Exact laws with predeclared bounds, the boundary rules of the protocol
//! scheduler and its request validation. Every statistical bound is a
//! simultaneous 99% Dvoretzky–Kiefer–Wolfowitz bound or a four standard
//! error screen, fixed before the tests were run.
use crn_craft::rules::{Control, RulesResult, Termination, execute, prepare};
use serde_json::{Value, json};

/// The three protocol sampler identities: Static Direct with Ziggurat
/// waits (the default), with inverse-transform waits, and next reaction.
const METHODS: [&str; 3] = ["static_direct", "static_direct_inverse", "next_reaction"];

fn model(name: &str, params: &Value, types: &[&str], initial: &Value, rules: &Value) -> Vec<u8> {
    let types: Vec<_> = types.iter().map(|t| json!({"name": t})).collect();
    let observables: serde_json::Map<String, Value> = types
        .iter()
        .map(|t| {
            let name = t["name"].as_str().unwrap().to_owned();
            (name.clone(), json!({"type": name}))
        })
        .collect();
    serde_json::to_vec(&json!({
        "schema": "crn-craft.model/1.0",
        "name": name,
        "parameters": params,
        "types": types,
        "initial": initial,
        "observables": observables,
        "rules": rules,
    }))
    .unwrap()
}

/// A run/1.1 request; `method` is a method kind.
fn request(
    model: &[u8],
    method: &str,
    replicates: u64,
    end: f64,
    reports: &[f64],
    protocol: &Value,
    threads: u64,
) -> Value {
    let method = if method == "static_direct" {
        json!({"kind": "static_direct", "selector": "linear"})
    } else if method == "static_direct_inverse" {
        json!({"kind": "static_direct", "selector": "linear", "waiting_time": "inverse"})
    } else {
        json!({"kind": method})
    };
    json!({
        "schema": "crn-craft.run/1.1",
        "model": {"path": "m.json", "sha256": crn_craft::sha256(model)},
        "method": method,
        "ensemble": {"replicates": replicates, "root_seed": 2027},
        "time": {"end": end, "report_times": reports},
        "outputs": {"statistics": ["mean"], "keep_paths": replicates, "keep_instances": 0},
        "budget": {"max_events_per_replicate": 1_000_000,
                   "max_sampling_steps_per_replicate": 100_000_000,
                   "max_instances": 64, "threads": threads,
                   "max_output_bytes": 1u64 << 28},
        "protocol": protocol,
    })
}

fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

fn run(model: &[u8], request: &Value) -> RulesResult {
    let prepared = prepare(model, &bytes(request)).unwrap_or_else(|d| panic!("{d:?}"));
    execute(&prepared, &Control::default()).unwrap_or_else(|d| panic!("{d:?}"))
}

fn refusal(model: &[u8], request: &Value) -> String {
    match prepare(model, &bytes(request)) {
        Ok(_) => panic!("request was accepted"),
        Err(d) => format!("{} {} {}", d[0].code, d[0].path, d[0].message),
    }
}

/// Values of observable `k` at report `i` on every kept path.
fn column(result: &RulesResult, k: usize, i: usize) -> Vec<u64> {
    result
        .paths
        .iter()
        .map(|p| {
            let x = p.values[i][k];
            assert!(x >= 0.0 && x.fract() == 0.0, "count {x}");
            x as u64
        })
        .collect()
}

fn ln_choose(n: u64, k: u64) -> f64 {
    (1..=k).map(|j| ((n - k + j) as f64 / j as f64).ln()).sum()
}

fn binomial_cdf(n: u64, p: f64) -> Vec<f64> {
    let mut total = 0.0;
    (0..=n)
        .map(|k| {
            total += (ln_choose(n, k) + k as f64 * p.ln() + (n - k) as f64 * (-p).ln_1p()).exp();
            total
        })
        .collect()
}

fn poisson_cdf(mean: f64, top: usize) -> Vec<f64> {
    let mut term = (-mean).exp();
    let mut total = term;
    let mut out = vec![total];
    for k in 1..=top {
        term *= mean / k as f64;
        total += term;
        out.push(total);
    }
    out
}

/// Kolmogorov distance between a sample and a CDF tabulated on 0, 1, ….
fn cdf_gap(sample: &[u64], cdf: &[f64]) -> f64 {
    let top = (*sample.iter().max().unwrap() as usize).max(cdf.len() - 1);
    let mut counts = vec![0u64; top + 1];
    for &x in sample {
        counts[x as usize] += 1;
    }
    let mut seen = 0;
    let mut gap: f64 = 0.0;
    for (k, &c) in counts.iter().enumerate() {
        seen += c;
        let law = cdf.get(k).copied().unwrap_or(1.0);
        gap = gap.max((seen as f64 / sample.len() as f64 - law).abs());
    }
    gap
}

/// Simultaneous 99% DKW bound for `curves` CDFs of `n` paths each.
fn dkw(n: usize, curves: usize) -> f64 {
    ((2.0 * curves as f64 / 0.01).ln() / (2.0 * n as f64)).sqrt()
}

fn death_model() -> Vec<u8> {
    model(
        "timed_death",
        &json!({"d": 0.0}),
        &["X"],
        &json!([{"type": "X", "count": 20}]),
        &json!([{"name": "death", "reactants": [{"label": "x", "type": "X"}],
                 "rate": "d", "products": []}]),
    )
}

/// E1: a switched-off channel becomes active, then its rate changes. A wait
/// drawn from the old hazard (or an old clock kept) changes both laws.
#[test]
fn rate_changes_use_the_new_hazard() {
    let m = death_model();
    let protocol = json!({"actions": [
        {"time": 1.0, "parameters": {"d": 0.2}},
        {"time": 3.0, "parameters": {"d": 1.0}}]});
    let n = 20_000;
    let bound = dkw(n, 2);
    for method in METHODS {
        let r = run(
            &m,
            &request(&m, method, n as u64, 4.0, &[1.0, 3.0, 4.0], &protocol, 2),
        );
        assert_eq!(r.termination, Termination::Completed);
        assert!(column(&r, 0, 0).iter().all(|&x| x == 20), "{method}: off");
        let g3 = cdf_gap(&column(&r, 0, 1), &binomial_cdf(20, (-0.4f64).exp()));
        let g4 = cdf_gap(&column(&r, 0, 2), &binomial_cdf(20, (-1.4f64).exp()));
        eprintln!("E1 {method}: gaps {g3:.5} {g4:.5} bound {bound:.5}");
        assert!(
            g3 < bound && g4 < bound,
            "{method}: gaps {g3} {g4} > {bound}"
        );
    }
}

/// E2: E[X+] = p E[X-] and Var(X+) = p² Var(X-) + p(1-p) E[X-] across a
/// binomial division of a non-Poisson (linear birth–death) state.
#[test]
fn division_moment_recurrences() {
    let m = model(
        "birth_death",
        &json!({"b": 1.0, "d": 0.5}),
        &["X"],
        &json!([{"type": "X", "count": 5}]),
        &json!([
            {"name": "birth", "reactants": [{"label": "x", "type": "X"}], "rate": "b",
             "products": [{"type": "X"}, {"type": "X"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "d",
             "products": []}]),
    );
    let p = 0.3;
    let protocol = json!({"actions": [{"time": 2.0, "partition": {"X": p}}]});
    let growth = 1.0f64.exp();
    let mean = 5.0 * growth;
    let var = 5.0 * 3.0 * growth * (growth - 1.0);
    let (want_mean, want_var) = (p * mean, p * p * var + p * (1.0 - p) * mean);
    let n = 20_000;
    for method in METHODS {
        let r = run(&m, &request(&m, method, n, 2.0, &[2.0], &protocol, 2));
        let x: Vec<f64> = column(&r, 0, 0).into_iter().map(|x| x as f64).collect();
        let k = x.len() as f64;
        let mu = x.iter().sum::<f64>() / k;
        let s2 = x.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / (k - 1.0);
        let m4 = x.iter().map(|v| (v - mu).powi(4)).sum::<f64>() / k;
        let z_mean = (mu - want_mean) / (s2 / k).sqrt();
        let z_var = (s2 - want_var) / ((m4 - s2 * s2) / k).sqrt();
        eprintln!("E2 {method}: z_mean {z_mean:.3} z_var {z_var:.3}");
        assert!(
            z_mean.abs() < 4.0 && z_var.abs() < 4.0,
            "{method}: mean {mu} vs {want_mean} (z {z_mean}), var {s2} vs {want_var} (z {z_var})"
        );
    }
}

/// Right-continuous reports at time zero, at an interior time with two
/// groups, and at the horizon; p = 1 and p = 0 partitions.
#[test]
fn boundary_rules() {
    let m = death_model();
    let protocol = json!({"actions": [
        {"time": 0.0, "set": {"X": 4}},
        {"time": 1.0, "set": {"X": 10}},
        {"time": 1.0, "change": {"X": -3}},
        {"time": 1.5, "partition": {"X": 1.0}},
        {"time": 2.0, "partition": {"X": 0.0}}]});
    for method in METHODS {
        let r = run(
            &m,
            &request(&m, method, 3, 2.0, &[0.0, 1.0, 1.5, 2.0], &protocol, 1),
        );
        for p in &r.paths {
            assert_eq!(p.values, vec![vec![4.0], vec![7.0], vec![7.0], vec![0.0]]);
        }
        // Groups at one time run in the declared order.
        let mut swapped = protocol.clone();
        swapped["actions"][1] = json!({"time": 1.0, "change": {"X": -3}});
        swapped["actions"][2] = json!({"time": 1.0, "set": {"X": 10}});
        let r = run(&m, &request(&m, method, 1, 2.0, &[1.0], &swapped, 1));
        assert_eq!(r.paths[0].values, vec![vec![10.0]]);
    }
}

/// Reports between interventions consume no random words, and the worker
/// count does not change a replicate.
#[test]
fn report_grid_and_workers_do_not_change_paths() {
    let m = death_model();
    let protocol = json!({"actions": [
        {"time": 0.5, "parameters": {"d": 0.7}, "partition": {"X": 0.6}},
        {"time": 1.5, "change": {"X": 5}}]});
    for method in METHODS {
        let a = run(&m, &request(&m, method, 64, 3.0, &[1.0, 3.0], &protocol, 1));
        let b = run(
            &m,
            &request(
                &m,
                method,
                64,
                3.0,
                &[0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0],
                &protocol,
                2,
            ),
        );
        for (pa, pb) in a.paths.iter().zip(&b.paths) {
            assert_eq!(pa.values[0], pb.values[3], "{method}");
            assert_eq!(pa.values[1], pb.values[6], "{method}");
        }
        let words = |r: &RulesResult| r.replicates.iter().map(|x| x.rng_words).collect::<Vec<_>>();
        assert_eq!(words(&a), words(&b), "{method}: reports drew random words");
    }
}

/// A removal beyond the available copies fails the replicate with a
/// diagnostic that names the action group, the species and its count.
#[test]
fn removal_beyond_copies_names_the_action() {
    let m = death_model();
    let protocol = json!({"actions": [{"time": 1.0, "change": {"X": -25}}]});
    let r = run(
        &m,
        &request(&m, "static_direct", 1, 2.0, &[2.0], &protocol, 1),
    );
    let record = &r.replicates[0];
    assert_eq!(record.termination, Termination::Failed);
    let d = record.diagnostic.as_ref().unwrap();
    assert_eq!(
        (d.code.as_str(), d.path.as_str()),
        ("CRN-R030", "/protocol/actions/0")
    );
    assert!(d.message.contains("\"X\" has 20 copies"), "{}", d.message);
}

/// The stop predicate is tested after each committed group, and a stopped
/// trajectory holds its counts.
#[test]
fn stop_after_a_group_holds_the_state() {
    let mut m: Value = serde_json::from_slice(&death_model()).unwrap();
    m["stop"] = json!("X >= 30");
    let m = bytes(&m);
    let protocol = json!({"actions": [
        {"time": 1.0, "change": {"X": 15}},
        {"time": 1.0, "set": {"X": 0}}]});
    let r = run(
        &m,
        &request(&m, "static_direct", 1, 2.0, &[0.5, 1.0, 2.0], &protocol, 1),
    );
    assert_eq!(r.replicates[0].termination, Termination::Stopped);
    assert_eq!(r.replicates[0].stop_time, 1.0);
    assert_eq!(r.paths[0].values, vec![vec![20.0], vec![35.0], vec![35.0]]);
}

#[test]
fn request_validation() {
    let m = death_model();
    let ok = request(
        &m,
        "static_direct",
        1,
        2.0,
        &[2.0],
        &json!({"actions": []}),
        1,
    );
    let edit = |f: &dyn Fn(&mut Value)| {
        let mut r = ok.clone();
        f(&mut r);
        refusal(&m, &r)
    };
    assert!(
        edit(&|r| {
            r.as_object_mut().unwrap().remove("protocol");
        })
        .contains("run/1.1 requires a protocol")
    );
    assert!(
        edit(&|r| r["schema"] = json!("crn-craft.run/1.0"))
            .contains("protocol requires crn-craft.run/1.1")
    );
    for kind in ["tree", "fsp", "mtdna_direct"] {
        assert!(edit(&|r| r["method"] = json!({"kind": kind})).starts_with("CRN-"));
    }
    assert!(
        edit(&|r| r["method"] = json!({"kind": "static_direct", "selector": "auto"}))
            .contains("CRN-R020 /method")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] =
            json!([{"time": 1.0, "set": {"X": 1}, "change": {"X": 1}}]))
        .contains("conflicting operations")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] = json!([{"time": 1.0, "parameters": {"q": 1.0}}]))
            .contains("unknown parameter")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] = json!([{"time": 2.5, "change": {"X": 1}}]))
            .contains("/protocol/actions/0/time")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] = json!([{"time": 1.0, "partition": {"X": 1.5}}]))
            .contains("partition probability")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] = json!([{"time": 1.0, "change": {"Y": 1}}]))
            .contains("unknown species")
    );
    assert!(
        edit(&|r| r["protocol"]["actions"] = json!([{"time": 1.0, "colour": 1}]))
            .starts_with("CRN-E")
    );
    // A parameter that sets an initial count is not a kinetic parameter.
    let mut counted: Value = serde_json::from_slice(&m).unwrap();
    counted["parameters"]["n0"] = json!(20);
    counted["initial"][0]["count"] = json!("n0");
    let counted = bytes(&counted);
    let usable = prepare(
        &counted,
        &bytes(&request(
            &counted,
            "static_direct",
            1,
            2.0,
            &[2.0],
            &json!({}),
            1,
        )),
    );
    eprintln!("parameter initial count accepted: {}", usable.is_ok());
    if usable.is_ok() {
        let r = request(
            &counted,
            "static_direct",
            1,
            2.0,
            &[2.0],
            &json!({"actions": [{"time": 1.0, "parameters": {"n0": 3.0}}]}),
            1,
        );
        assert!(refusal(&counted, &r).contains("initial counts or channel topology"));
    }
    // Duplicate keys are rejected by the strict reader.
    let text = String::from_utf8(bytes(&ok)).unwrap().replace(
        "\"actions\":[]",
        "\"actions\":[{\"time\":1.0,\"change\":{\"X\":1,\"X\":2}}]",
    );
    assert!(prepare(&m, text.as_bytes()).is_err());
}

fn copy_model() -> Vec<u8> {
    model(
        "busy_copy",
        &json!({"r": 2.5, "loss": 0.7}),
        &["Free", "Busy", "Product"],
        &json!([]),
        &json!([
            {"name": "start", "reactants": [], "rate": "r", "products": [{"type": "Busy"}]},
            {"name": "lose", "reactants": [{"label": "x", "type": "Busy"}], "rate": "loss",
             "products": []}]),
    )
}

/// E4: Poisson starts at rate r with an irrevocable delay τ give
/// Poisson(r max(0, t - τ)) completed products by t.
#[test]
fn irrevocable_completions_are_poisson() {
    let m = model(
        "delayed_arrivals",
        &json!({"r": 2.5}),
        &["Product"],
        &json!([]),
        &json!([{"name": "start", "reactants": [], "rate": "r", "products": []}]),
    );
    let protocol = json!({"delays": [{"reaction": "start", "duration": 1.0,
                                      "products": {"Product": 1}}]});
    let n = 20_000;
    for method in METHODS {
        let r = run(
            &m,
            &request(&m, method, n as u64, 3.0, &[0.5, 3.0], &protocol, 2),
        );
        assert!(column(&r, 0, 0).iter().all(|&x| x == 0), "{method}: early");
        let gap = cdf_gap(&column(&r, 0, 1), &poisson_cdf(5.0, 40));
        eprintln!("E4 {method}: gap {gap:.5} bound {:.5}", dkw(n, 1));
        assert!(gap < dkw(n, 1), "{method}: gap {gap}");
    }
}

/// Owner death at rate δ multiplies the mean by exp(-δτ); a single start
/// yields at most one completion.
#[test]
fn owned_completions_are_cancelled_by_owner_loss() {
    let m = copy_model();
    let protocol = json!({"delays": [{"reaction": "start", "duration": 1.0, "owner": "Busy",
                                      "products": {"Product": 1}}]});
    let n = 20_000;
    for method in METHODS {
        let r = run(
            &m,
            &request(&m, method, n as u64, 3.0, &[3.0], &protocol, 2),
        );
        let gap = cdf_gap(&column(&r, 2, 0), &poisson_cdf(5.0 * (-0.7f64).exp(), 40));
        eprintln!("E4-owner {method}: gap {gap:.5} bound {:.5}", dkw(n, 1));
        assert!(gap < dkw(n, 1), "{method}: gap {gap}");
    }
    let one = model(
        "single",
        &json!({}),
        &["Free", "Busy", "Product"],
        &json!([{"type": "Free", "count": 1}]),
        &json!([{"name": "start", "reactants": [{"label": "f", "type": "Free"}], "rate": "1",
                 "products": [{"type": "Busy"}]}]),
    );
    let r = run(
        &one,
        &request(&one, "next_reaction", 200, 5.0, &[5.0], &protocol, 1),
    );
    assert!(column(&r, 2, 0).iter().all(|&x| x <= 1));
}

/// E3: a partition that keeps one of two busy owners keeps a uniformly
/// chosen one: the first owner's product appears in half of such paths.
#[test]
fn partition_cancels_a_uniform_subset_of_owners() {
    let m = model(
        "two_owners",
        &json!({"fast": 1.0e6}),
        &["Free", "Busy", "Product"],
        &json!([{"type": "Free", "count": 1}]),
        &json!([{"name": "start", "reactants": [{"label": "f", "type": "Free"}],
                 "rate": "fast", "products": [{"type": "Busy"}]}]),
    );
    let protocol = json!({
        "delays": [{"reaction": "start", "duration": 1.0, "owner": "Busy",
                    "products": {"Product": 1}}],
        "actions": [{"time": 0.25, "change": {"Free": 1}},
                    {"time": 0.5, "partition": {"Busy": 0.5}}]});
    let n = 20_000;
    for method in METHODS {
        let r = run(&m, &request(&m, method, n, 2.0, &[1.1, 1.4], &protocol, 2));
        let (early, late) = (column(&r, 2, 0), column(&r, 2, 1));
        let kept_one: Vec<_> = early.iter().zip(&late).filter(|&(_, &l)| l == 1).collect();
        let first = kept_one.iter().filter(|&&(&e, _)| e == 1).count() as f64;
        let k = kept_one.len() as f64;
        let z = (first / k - 0.5) / (0.25 / k).sqrt();
        eprintln!("E3 {method}: {first} of {k} first kept, z {z:.3}");
        assert!(
            k > 8_000.0 && z.abs() < 4.0,
            "{method}: {first} of {k} (z {z})"
        );
    }
}

/// The pending queue is bounded by `max_pending`.
#[test]
fn pending_budget() {
    let m = copy_model();
    let protocol = json!({"max_pending": 2, "delays": [{"reaction": "start", "duration": 10.0,
        "owner": "Busy", "products": {"Product": 1}}]});
    let r = run(
        &m,
        &request(&m, "static_direct", 4, 5.0, &[5.0], &protocol, 1),
    );
    assert!(
        r.replicates
            .iter()
            .all(|x| x.termination == Termination::Failed)
    );
    let mut over = protocol.clone();
    over["max_pending"] = json!(0);
    assert!(
        refusal(&m, &request(&m, "static_direct", 1, 5.0, &[5.0], &over, 1)).contains("capacity")
    );
}
