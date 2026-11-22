use super::*;
use crate::rules::FlatView;
use serde_json::{Value, json};

fn view(model: &Value) -> FlatView {
    FlatView::from_bytes(&serde_json::to_vec(model).unwrap()).unwrap()
}

fn birth_death(birth: f64, death: f64, start: u64) -> FlatView {
    view(&json!({
        "schema": "crn-craft.model/1.0", "name": "bd",
        "parameters": {"b": birth, "d": death},
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": start}],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {"name": "birth", "reactants": [], "rate": "b", "products": [{"type": "X"}]},
            {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "d", "products": []}
        ]
    }))
}

fn spec(view: &FlatView, bounds: Vec<u64>, initial: Vec<Vec<i64>>) -> Spec {
    Spec {
        species: view.species().to_vec(),
        bounds,
        initial,
        channels: ChannelShape::of_view(view),
        interventions: Vec::new(),
        max_states: 1000,
        max_edges: 10_000,
    }
}

fn never() -> bool {
    false
}

fn poisson(k: u64, mean: f64) -> f64 {
    special::poisson_pmf(k as f64, mean)
}

#[test]
fn enumeration_is_breadth_first_with_the_sink_last() {
    let v = birth_death(1.0, 1.0, 1);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![3], vec![vec![1]]), &interrupt).unwrap();
    assert_eq!(t.states(), [1, 2, 0, 3]);
    // From 1: birth -> 2, death -> 0; from 2: 3, 1; from 0: 1; from 3: sink, 2.
    assert_eq!(t.sources(), [0, 0, 1, 1, 2, 3, 3]);
    assert_eq!(t.targets(), [1, 2, 3, 0, 0, 4, 1]);
    assert_eq!(t.channels(), [0, 1, 0, 1, 0, 0, 1]);
    assert_eq!(t.multiplicities(), [1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0]);
    assert_eq!(t.index_of(&[3]), Some(3));
    assert_eq!(t.index_of(&[4]), None);
}

#[test]
fn second_order_multiplicity_and_interventions_expand_the_states() {
    let v = view(&json!({
        "schema": "crn-craft.model/1.0", "name": "dimer",
        "parameters": {"k": 0.5},
        "types": [{"name": "A"}, {"name": "B"}],
        "initial": [{"type": "A", "count": 4}],
        "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
        "rules": [{"name": "bind", "reactants": [{"label": "x", "type": "A", "count": 2}],
                   "rate": "k", "products": [{"type": "B"}]}]
    }));
    let interrupt = Interrupt::new(&never);
    let mut s = spec(&v, vec![6, 3], vec![vec![4, 0]]);
    s.interventions = vec![vec![(0, Operation::Change(2))]];
    let t = Topology::enumerate(&s, &interrupt).unwrap();
    // (4,0) -> (2,1) -> (0,2); then +2 A on each: (6,0) (4,1) (2,2) -> ...
    assert_eq!(&t.states()[..6], [4, 0, 2, 1, 0, 2]);
    assert_eq!(&t.states()[6..12], [6, 0, 4, 1, 2, 2]);
    assert_eq!(t.multiplicities()[0], 6.0); // C(4, 2)
    let generator = Generator::build(&t, &v, &[], &interrupt).unwrap();
    let q = generator.constant();
    // Columns of a generator sum to zero, the sink column is empty.
    for j in 0..q.dim {
        let sum: f64 = q.values[q.col_ptr[j]..q.col_ptr[j + 1]].iter().sum();
        assert!(sum.abs() < 1e-12);
    }
    assert_eq!(q.diagonal(0), -3.0); // 0.5 * C(4, 2)
}

#[test]
fn budgets_are_refused_with_their_messages() {
    let v = birth_death(1.0, 1.0, 1);
    let interrupt = Interrupt::new(&never);
    let mut s = spec(&v, vec![100], vec![vec![1]]);
    s.max_states = 5;
    let error = Topology::enumerate(&s, &interrupt).unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.message.as_str()),
        ("CRN-E070", "reachable-state budget exceeded")
    );
    s.max_states = 1000;
    s.max_edges = 3;
    let error = Topology::enumerate(&s, &interrupt).unwrap_err();
    assert_eq!(error.message, "projection edge budget exceeded");
    let s = spec(&v, vec![3], vec![vec![5]]);
    let error = Topology::enumerate(&s, &interrupt).unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.message.as_str()),
        ("CRN-E071", "initial law lies outside the projection")
    );
    let s = spec(&v, vec![3], vec![vec![-1]]);
    let error = Topology::enumerate(&s, &interrupt).unwrap_err();
    assert_eq!(error.message, "invalid initial or reachable counts");
    let cancelled = || true;
    let mut s = spec(&v, vec![10_000], vec![vec![1]]);
    s.max_states = 100_000;
    let error = Topology::enumerate(&s, &Interrupt::new(&cancelled)).unwrap_err();
    assert_eq!(error.message, "projection interrupted");
}

#[test]
fn uniformization_reproduces_immigration_death() {
    // From 0 copies, birth 10, death 1: Poisson(10 (1 - e^-t)).
    let v = birth_death(10.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![80], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[], &interrupt).unwrap();
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 10_000_000);
    let tail = uniform(&g, &mut p, 2.0, 1e-12, &mut work, &interrupt).unwrap();
    assert!(tail <= 1e-12 && tail > 0.0);
    let mean = 10.0 * (1.0 - (-2.0f64).exp());
    for (i, &x) in t.states().iter().enumerate() {
        assert!((p[i] - poisson(x, mean)).abs() < 1e-11, "x={x}");
    }
    // Work: one unit per series term.
    assert!(work.count > 0);
    let mut small = Work::new(0, 20);
    let mut q = vec![0.0; t.len() + 1];
    q[0] = 1.0;
    let error = uniform(&g, &mut q, 2.0, 1e-12, &mut small, &interrupt).unwrap_err();
    assert_eq!(
        (error.path.as_str(), error.message.as_str()),
        (
            "/projection/work",
            "projection work limit of 20 exceeded; raise max_work"
        )
    );
    // Pieces beyond the remaining work are refused before any term.
    let mut tiny = Work::new(0, 1);
    let error = uniform(&g, &mut q, 1000.0, 1e-12, &mut tiny, &interrupt).unwrap_err();
    assert_eq!(error.message, "uniformization exceeds the work budget");
    assert_eq!(tiny.count, 0);
}

struct Linear;
impl SignalValues for Linear {
    fn evaluate(
        &mut self,
        time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        out[0] = 2.0 + time;
        Ok(())
    }
}

#[test]
fn midpoint_integrates_a_noncommuting_forced_law() {
    // Birth 2 + t, death 1, from 0: Poisson with mean 3 - e^-2 at t = 2.
    let v = birth_death(2.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![48], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let merged = Merged::new(&g);
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 1_000_000);
    let controls = MidpointControls {
        start: 0.0,
        end: 2.0,
        tolerance: 1e-7,
        max_step: f64::INFINITY,
    };
    midpoint(
        &g,
        &merged,
        &mut p,
        controls,
        &mut Linear,
        &mut work,
        &interrupt,
    )
    .unwrap();
    let mean = 3.0 - (-2.0f64).exp();
    let gap: f64 = t
        .states()
        .iter()
        .enumerate()
        .map(|(i, &x)| (p[i] - poisson(x, mean)).abs())
        .sum();
    assert!(gap < 1e-6, "{gap}");
    assert_eq!(work.count % 3, 0);
    // A goal below the rounding level of the estimate is refused.
    let mut q = vec![0.0; t.len() + 1];
    q[0] = 1.0;
    let tight = MidpointControls {
        tolerance: 1e-14,
        ..controls
    };
    let error = midpoint(
        &g,
        &merged,
        &mut q,
        tight,
        &mut Linear,
        &mut work,
        &interrupt,
    )
    .unwrap_err();
    assert!(error.message.contains("rounding level"));
}

struct Step {
    seen: Vec<f64>,
}
impl SignalValues for Step {
    fn evaluate(
        &mut self,
        time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        self.seen.push(time);
        out[0] = if time < 1.0 { 1.0 } else { -1.0 };
        Ok(())
    }
}

#[test]
fn the_end_of_an_interval_sees_the_incoming_signal() {
    // max_step 2 over [0, 1]: the first step's midpoint is inside, and
    // the end value is never requested at exactly 1.
    let v = birth_death(1.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![5], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let merged = Merged::new(&g);
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 1000);
    let mut signal = Step { seen: Vec::new() };
    let controls = MidpointControls {
        start: 0.0,
        end: 1.0,
        tolerance: 1e-2,
        max_step: 2.0,
    };
    midpoint(
        &g,
        &merged,
        &mut p,
        controls,
        &mut signal,
        &mut work,
        &interrupt,
    )
    .unwrap();
    assert!(signal.seen.iter().all(|&s| s < 1.0));
    // A negative value is refused.
    let late = MidpointControls {
        start: 1.0,
        end: 2.0,
        ..controls
    };
    let error = midpoint(
        &g,
        &merged,
        &mut p,
        late,
        &mut signal,
        &mut work,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.message.as_str()),
        (
            "CRN-E071",
            "time-dependent reaction rate must be finite and nonnegative"
        )
    );
}

#[test]
fn kicks_partition_set_and_refuse() {
    let v = birth_death(1.0, 1.0, 4);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![6], vec![vec![4]]), &interrupt).unwrap();
    let mut p = vec![0.0; t.len() + 1];
    p[t.index_of(&[4]).unwrap()] = 0.75;
    p[t.len()] = 0.25;
    let mut work = Work::new(0, 1000);
    let half = kick(
        &t,
        &p,
        &[(0, Operation::Partition(0.5))],
        &mut work,
        &interrupt,
    )
    .unwrap();
    assert_eq!(work.count, 5);
    for k in 0..=4u64 {
        let exact = [1.0, 4.0, 6.0, 4.0, 1.0][k as usize] / 16.0 * 0.75;
        assert!((half[t.index_of(&[k]).unwrap()] - exact).abs() < 1e-16);
    }
    assert_eq!(half[t.len()], 0.25);
    let all = kick(
        &t,
        &p,
        &[(0, Operation::Partition(1.0))],
        &mut work,
        &interrupt,
    )
    .unwrap();
    assert_eq!(all, p);
    let high = kick(&t, &p, &[(0, Operation::Set(7))], &mut work, &interrupt).unwrap();
    assert_eq!(high[t.len()], 1.0);
    let error = kick(&t, &p, &[(0, Operation::Change(-5))], &mut work, &interrupt).unwrap_err();
    assert_eq!(error.message, "intervention removes more X than available");
    let mut small = Work::new(0, 2);
    let error = kick(
        &t,
        &p,
        &[(0, Operation::Partition(0.5))],
        &mut small,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(error.path, "/projection/work");
}

#[test]
fn signals_need_count_independent_rules_and_stop_states_absorb() {
    let v = view(&json!({
        "schema": "crn-craft.model/1.0", "name": "stopper",
        "parameters": {"k": 1.0},
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 0}],
        "observables": {"X": {"type": "X"}},
        "stop": "X >= 2",
        "rules": [
            {"name": "grow", "reactants": [], "rate": "k", "products": [{"type": "X"}]},
            {"name": "feedback", "reactants": [], "rate": "k / (1 + X)", "products": [{"type": "X"}]}
        ]
    }));
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![3], vec![vec![0]]), &interrupt).unwrap();
    let error = Generator::build(&t, &v, &[vec![1]], &interrupt).unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-E071", "/signals/0")
    );
    let g = Generator::build(&t, &v, &[], &interrupt).unwrap();
    let stop = t.index_of(&[2]).unwrap();
    assert!(g.stopped()[stop]);
    let q = g.constant();
    assert_eq!(q.col_ptr[stop], q.col_ptr[stop + 1]);
    // State-dependent rate: k / (1 + X) at X = 1 is 0.5, plus k.
    assert_eq!(q.diagonal(t.index_of(&[1]).unwrap()), -1.5);
}
