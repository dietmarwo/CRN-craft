use super::implicit::implicit_alpha;
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

struct Constant;
impl SignalValues for Constant {
    fn evaluate(
        &mut self,
        _time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        out[0] = 3.0;
        Ok(())
    }
}

#[test]
fn midpoint_steps_end_exactly_at_the_interval_end() {
    // Birth 3, death 1, from 0: Poisson with mean 3 (1 - e^-t). Twenty
    // steps of 0.05 do not sum to 1 in floating point; the last one ends
    // at 1 anyway and no shorter step follows.
    let v = birth_death(2.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![48], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let merged = Merged::new(&g);
    for (end, max_step, actions) in [
        (1.0, 0.05, 60),
        (0.7, 0.1, 21),
        (1.0, 0.3, 12),
        (1.0, 0.6, 6),
        (1.0, f64::INFINITY, 3),
    ] {
        let mut p = vec![0.0; t.len() + 1];
        p[0] = 1.0;
        let mut work = Work::new(0, 1000);
        let controls = MidpointControls {
            start: 0.0,
            end,
            tolerance: 1e-6,
            max_step,
        };
        midpoint(
            &g,
            &merged,
            &mut p,
            controls,
            &mut Constant,
            &mut work,
            &interrupt,
        )
        .unwrap();
        let mean = 3.0 * (1.0 - (-end).exp());
        let gap: f64 = t
            .states()
            .iter()
            .enumerate()
            .map(|(i, &x)| (p[i] - poisson(x, mean)).abs())
            .sum();
        assert!(gap < 1e-6, "{end} {max_step} {gap}");
        assert_eq!(work.count, actions, "{end} {max_step}");
    }
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

/// A dense LU with partial pivoting: the test double of the sparse solver.
#[derive(Default)]
struct DenseLu {
    n: usize,
    lu: Vec<f64>,
    pivots: Vec<usize>,
    analyzed: Option<(Vec<usize>, Vec<usize>)>,
    factorizations: usize,
}

impl SparseLinearSolver for DenseLu {
    fn analyze(&mut self, pattern: &Csc) -> Result<u64, crate::diagnostic::Diagnostic> {
        self.analyzed = Some((pattern.col_ptr.clone(), pattern.row_idx.clone()));
        Ok((pattern.dim * pattern.dim * 8) as u64)
    }

    fn factorize(&mut self, matrix: &Csc) -> Result<(), crate::diagnostic::Diagnostic> {
        let (col_ptr, row_idx) = self.analyzed.as_ref().expect("analyzed");
        assert_eq!((col_ptr, row_idx), (&matrix.col_ptr, &matrix.row_idx));
        let n = matrix.dim;
        self.n = n;
        self.lu = vec![0.0; n * n];
        for j in 0..n {
            for k in matrix.col_ptr[j]..matrix.col_ptr[j + 1] {
                self.lu[matrix.row_idx[k] * n + j] = matrix.values[k];
            }
        }
        self.pivots = (0..n).collect();
        for k in 0..n {
            let p = (k..n)
                .max_by(|&a, &b| {
                    self.lu[a * n + k]
                        .abs()
                        .total_cmp(&self.lu[b * n + k].abs())
                })
                .unwrap();
            if self.lu[p * n + k] == 0.0 {
                return Err(crate::diagnostic::Diagnostic::error(
                    "CRN-E070", "/test", "singular", "",
                ));
            }
            if p != k {
                for j in 0..n {
                    self.lu.swap(p * n + j, k * n + j);
                }
                self.pivots.swap(p, k);
            }
            for i in k + 1..n {
                let f = self.lu[i * n + k] / self.lu[k * n + k];
                self.lu[i * n + k] = f;
                for j in k + 1..n {
                    self.lu[i * n + j] -= f * self.lu[k * n + j];
                }
            }
        }
        self.factorizations += 1;
        Ok(())
    }

    fn solve_in_place(&mut self, rhs: &mut [f64]) -> Result<(), crate::diagnostic::Diagnostic> {
        let n = self.n;
        let mut x: Vec<f64> = self.pivots.iter().map(|&p| rhs[p]).collect();
        for i in 0..n {
            for j in 0..i {
                x[i] -= self.lu[i * n + j] * x[j];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= self.lu[i * n + j] * x[j];
            }
            x[i] /= self.lu[i * n + i];
        }
        rhs.copy_from_slice(&x);
        Ok(())
    }
}

/// Birth `a (1 + sin(w t))`, death `k`: a Poisson law with the mean of
/// `m' = b(t) - k m`, `m(0) = 0`.
struct Wave;
const WAVE: (f64, f64, f64) = (4.0, 1.3, 0.7);
impl SignalValues for Wave {
    fn evaluate(
        &mut self,
        time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        out[0] = WAVE.0 * (1.0 + (WAVE.1 * time).sin());
        Ok(())
    }
}

fn wave_mean(t: f64) -> f64 {
    let (a, w, k) = WAVE;
    let decay = (-k * t).exp();
    let swing = k * (w * t).sin() - w * (w * t).cos() + w * decay;
    a / k * (1.0 - decay) + a * swing / (k * k + w * w)
}

fn bdf_controls(end: f64, tolerance: f64, max_step: f64) -> BdfControls {
    BdfControls {
        start: 0.0,
        end,
        tolerance,
        max_step,
        max_factor_bytes: u64::MAX,
    }
}

#[test]
fn bdf_follows_a_forced_poisson_law_at_reports_and_the_end() {
    let v = birth_death(4.0, WAVE.2, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![60], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let mut implicit = Implicit::new(&g, DenseLu::default()).unwrap();
    for (tolerance, bound) in [(1e-6, 2e-6), (1e-9, 5e-9)] {
        let before = implicit.counts();
        let mut p = vec![0.0; t.len() + 1];
        p[0] = 1.0;
        let mut work = Work::new(0, 10_000_000);
        let reports = [0.7, 2.0, 3.3];
        let out = bdf(
            &mut implicit,
            &mut p,
            bdf_controls(5.0, tolerance, 0.5),
            &reports,
            &mut Wave,
            &mut work,
            &interrupt,
        )
        .unwrap();
        let gap = |law: &[f64], at: f64| -> f64 {
            t.states()
                .iter()
                .enumerate()
                .map(|(i, &x)| (law[i] - poisson(x, wave_mean(at))).abs())
                .fold(0.0, f64::max)
        };
        for (law, &at) in out.iter().zip(&reports) {
            assert!(gap(law, at) < bound, "{tolerance} {at} {}", gap(law, at));
        }
        assert!(gap(&p, 5.0) < bound, "{tolerance} {}", gap(&p, 5.0));
        let c = implicit.counts();
        assert_eq!(
            work.count as u64,
            c.products + c.solves + c.numeric - before.products - before.solves - before.numeric,
            "work unit"
        );
    }
    let c = implicit.counts();
    assert_eq!(c.symbolic, 1);
    // Equal steps reuse a factorization: fewer factorizations than steps.
    assert!(c.numeric < c.steps, "{c:?}");
    assert_eq!(c.numeric as usize, implicit.solver().factorizations);
}

#[test]
fn bdf_with_a_constant_signal_matches_uniformization() {
    let v = birth_death(2.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![40], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let mut implicit = Implicit::new(&g, DenseLu::default()).unwrap();
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 1_000_000);
    bdf(
        &mut implicit,
        &mut p,
        bdf_controls(2.0, 1e-10, f64::INFINITY),
        &[],
        &mut Constant,
        &mut work,
        &interrupt,
    )
    .unwrap();
    let mean = 3.0 * (1.0 - (-2.0f64).exp());
    let gap: f64 = t
        .states()
        .iter()
        .enumerate()
        .map(|(i, &x)| (p[i] - poisson(x, mean)).abs())
        .sum();
    assert!(gap < 1e-8, "{gap}");
    // Exact factorizations: no iteration, one solve per step attempt.
    let c = implicit.counts();
    assert_eq!(c.solves, c.steps + c.rejected, "{c:?}");
}

#[test]
fn bdf_refuses_memory_work_and_negative_signals() {
    let v = birth_death(1.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![5], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let mut implicit = Implicit::new(&g, DenseLu::default()).unwrap();
    // The dense double's 7 x 7 values, the 18 entries of the iteration
    // matrix and 24 vectors of 7 entries.
    let bound = 7 * 7 * 8 + 24 * 18 + 8 * 24 * 7;
    assert_eq!(implicit.factor_bytes(), bound);
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 1000);
    let mut controls = bdf_controls(1.0, 1e-6, f64::INFINITY);
    controls.max_factor_bytes = bound - 1;
    let error = bdf(
        &mut implicit,
        &mut p,
        controls,
        &[],
        &mut Constant,
        &mut work,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-E070", "/projection/memory")
    );
    assert!(error.message.contains(&format!(
        "{bound} bytes, above the memory budget max_factor_bytes={}",
        bound - 1
    )));
    assert_eq!(implicit.counts().numeric, 0);
    let mut small = Work::new(0, 5);
    let error = bdf(
        &mut implicit,
        &mut p,
        bdf_controls(1.0, 1e-6, f64::INFINITY),
        &[],
        &mut Constant,
        &mut small,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(error.path, "/projection/work");
    let mut signal = Step { seen: Vec::new() };
    let late = BdfControls {
        start: 1.0,
        end: 2.0,
        ..bdf_controls(1.0, 1e-6, f64::INFINITY)
    };
    let error = bdf(
        &mut implicit,
        &mut p,
        late,
        &[],
        &mut signal,
        &mut work,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(error.code, "CRN-E071");
    // The interval [0, 1] never sees the value at 1 itself.
    let mut signal = Step { seen: Vec::new() };
    let mut work = Work::new(0, 100_000);
    bdf(
        &mut implicit,
        &mut p,
        bdf_controls(1.0, 1e-6, f64::INFINITY),
        &[],
        &mut signal,
        &mut work,
        &interrupt,
    )
    .unwrap();
    assert!(signal.seen.iter().all(|&s| s < 1.0));
    let stop = || true;
    let error = bdf(
        &mut implicit,
        &mut p,
        bdf_controls(1.0, 1e-6, f64::INFINITY),
        &[],
        &mut Constant,
        &mut work,
        &Interrupt::new(&stop),
    )
    .unwrap_err();
    assert_eq!(error.message, "projection interrupted");
}

fn chain(rates: &[(usize, usize, f64)], n: usize) -> Generator {
    // States 0..n and the sink n; rate (s, d, r) from s to d.
    Generator::from_rates(n, rates)
}

#[test]
fn stationary_finds_closed_classes_and_anchors_one() {
    // 0 <-> 1 closed, 2 -> 0 and 2 -> 3, 3 <-> 4 closed: two closed classes
    // and one transient state.
    let rates = [
        (0, 1, 2.0),
        (1, 0, 1.0),
        (2, 0, 1.0),
        (2, 3, 1.0),
        (3, 4, 1.0),
        (4, 3, 3.0),
    ];
    let g = chain(&rates, 5);
    let interrupt = Interrupt::new(&never);
    let mut solver = DenseLu::default();
    let error = stationary(&g, None, None, &mut solver, u64::MAX, &interrupt).unwrap_err();
    assert_eq!(
        (error.code.as_str(), error.message.as_str()),
        ("CRN-E071", "projection has 2 closed communicating classes")
    );
    let law = stationary(&g, None, Some(4), &mut solver, u64::MAX, &interrupt).unwrap();
    let expected = [0.0, 0.0, 0.0, 0.75, 0.25];
    for (a, b) in law.probabilities.iter().zip(expected) {
        assert!((a - b).abs() < 1e-15, "{law:?}");
    }
    let law = stationary(&g, None, Some(1), &mut solver, u64::MAX, &interrupt).unwrap();
    assert!((law.probabilities[0] - 1.0 / 3.0).abs() < 1e-15);
    let error = stationary(&g, None, Some(2), &mut solver, u64::MAX, &interrupt).unwrap_err();
    assert_eq!(
        error.message,
        "selected state is outside every closed class"
    );
    let error = stationary(&g, None, Some(4), &mut solver, 3, &interrupt).unwrap_err();
    assert_eq!(error.path, "/projection/memory");
}

#[test]
fn stationary_redirects_escapes_to_the_return_state() {
    // 0 -> 1 -> 2 -> sink, returning to 0: a cycle with rates 1, 2, 4.
    let g = chain(&[(0, 1, 1.0), (1, 2, 2.0), (2, 3, 4.0)], 3);
    assert_eq!(g.escapes(), vec![0.0, 0.0, 4.0]);
    let interrupt = Interrupt::new(&never);
    let mut solver = DenseLu::default();
    let law = stationary(&g, Some(0), None, &mut solver, u64::MAX, &interrupt).unwrap();
    let expected = [4.0 / 7.0, 2.0 / 7.0, 1.0 / 7.0];
    for (a, b) in law.probabilities.iter().zip(expected) {
        assert!((a - b).abs() < 1e-15);
    }
    assert!((law.flux - 4.0 / 7.0).abs() < 1e-15);
    // Return to the absorbing-free middle: 1 is in the closed class too.
    let error = stationary(
        &chain(&[(0, 1, 1.0), (1, 0, 1.0), (2, 3, 1.0)], 3),
        Some(2),
        None,
        &mut solver,
        u64::MAX,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(
        error.message,
        "the redirected projection has multiple closed classes"
    );
    let error = stationary(
        &chain(&[(0, 1, 1.0), (1, 2, 1.0)], 3),
        Some(0),
        None,
        &mut solver,
        u64::MAX,
        &interrupt,
    )
    .unwrap_err();
    assert_eq!(
        error.message,
        "the designated state is outside every closed class"
    );
    // An absorbing state is its own closed class.
    let law = stationary(
        &chain(&[(0, 1, 1.0), (1, 2, 1.0)], 3),
        None,
        None,
        &mut solver,
        u64::MAX,
        &interrupt,
    )
    .unwrap();
    assert_eq!(law.probabilities, vec![0.0, 0.0, 1.0]);
}

/// Birth 1 until t = 0.5, then 40 (a jump the integrator is not told
/// about), death 1.
struct Jump;
impl SignalValues for Jump {
    fn evaluate(
        &mut self,
        time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        out[0] = if time < 0.5 { 1.0 } else { 40.0 };
        Ok(())
    }
}

#[test]
fn bdf_converges_through_an_undeclared_jump() {
    // Each corrector must solve its own equation: an unconverged corrector
    // after the jump (a factorization from before it) shows in the law.
    let v = birth_death(1.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![120], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let mut implicit = Implicit::new(&g, DenseLu::default()).unwrap();
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 10_000_000);
    bdf(
        &mut implicit,
        &mut p,
        bdf_controls(1.5, 1e-9, 0.05),
        &[],
        &mut Jump,
        &mut work,
        &interrupt,
    )
    .unwrap();
    let at_jump = 1.0 - (-0.5f64).exp();
    let mean = at_jump * (-1.0f64).exp() + 40.0 * (1.0 - (-1.0f64).exp());
    let gap: f64 = t
        .states()
        .iter()
        .enumerate()
        .map(|(i, &x)| (p[i] - poisson(x, mean)).abs())
        .sum();
    assert!(gap < 1e-6, "{gap}");
}

#[test]
fn ndf_coefficients_are_the_published_ones() {
    // Shampine and Reichelt (1997), table 1: kappa and the error constants.
    let kappa = [-0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0];
    for (k, &expected) in kappa.iter().enumerate() {
        let order = k + 1;
        let gamma: f64 = (1..=order).map(|j| 1.0 / j as f64).sum();
        // alpha_k = (1 - kappa_k) gamma_k, so kappa_k = 1 - alpha_k / gamma_k.
        let alpha = implicit_alpha(order);
        assert!((1.0 - alpha / gamma - expected).abs() < 1e-15, "{order}");
    }
}

/// A stiff forced law: birth `2000 (1 + sin(t) / 2)`, death 1000 per
/// molecule.
struct StiffWave;
impl SignalValues for StiffWave {
    fn evaluate(
        &mut self,
        time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        out[0] = 2000.0 * (1.0 + 0.5 * time.sin());
        Ok(())
    }
}

#[test]
fn bdf_iterates_stale_factorizations_of_a_stiff_forced_law() {
    // Equal steps reuse the factorization of an earlier time; with stiff
    // rates its error is large, and only converged iterations keep the law.
    let v = birth_death(2000.0, 1000.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![30], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    let mut implicit = Implicit::new(&g, DenseLu::default()).unwrap();
    let mut p = vec![0.0; t.len() + 1];
    p[0] = 1.0;
    let mut work = Work::new(0, 10_000_000);
    let end = 3.0;
    bdf(
        &mut implicit,
        &mut p,
        bdf_controls(end, 1e-8, f64::INFINITY),
        &[],
        &mut StiffWave,
        &mut work,
        &interrupt,
    )
    .unwrap();
    let (a, k) = (2000.0, 1000.0);
    let decay = (-k * end).exp();
    let mean =
        a / k * (1.0 - decay) + 0.5 * a * (k * end.sin() - end.cos() + decay) / (k * k + 1.0);
    let gap: f64 = t
        .states()
        .iter()
        .enumerate()
        .map(|(i, &x)| (p[i] - poisson(x, mean)).abs())
        .sum();
    let c = implicit.counts();
    assert!(gap < 1e-6, "{gap} {c:?}");
    assert!(c.solves > c.steps + c.rejected, "{c:?}");
}

/// No forced channels.
struct Silent;
impl SignalValues for Silent {
    fn evaluate(
        &mut self,
        _time: f64,
        out: &mut [f64],
    ) -> Result<(), crate::diagnostic::Diagnostic> {
        assert!(out.is_empty());
        Ok(())
    }
}

/// Forward sensitivity of an immigration-death law to the birth rate: the
/// Poisson law with mean `m = b (1 - e^{-k t}) / k` has derivative
/// `p(x) (x / m - 1) dm/db`. On a domain cut low enough that the sink takes
/// visible mass the augmented solve still conserves probability, keeps the
/// sensitivity sum (sink included) at zero, and equals the centered
/// difference of two truncated solves.
#[test]
fn tangents_match_the_poisson_derivative_and_conserve_mass() {
    let (b, k, h) = (3.0, 0.8, 1e-4);
    let interrupt = Interrupt::new(&never);
    let reports = [0.5, 1.5];
    let solve = |generator: &Generator, tangent: Option<usize>, y: &mut Vec<f64>| {
        let mut work = Work::new(0, 1_000_000);
        let controls = bdf_controls(3.0, 1e-10, f64::INFINITY);
        let mut out = match tangent {
            Some(n) => {
                let solver = TangentSolver::new(n, DenseLu::default());
                let mut implicit = Implicit::new(generator, solver).unwrap();
                bdf(
                    &mut implicit,
                    y,
                    controls,
                    &reports,
                    &mut Silent,
                    &mut work,
                    &interrupt,
                )
            }
            None => {
                let mut implicit = Implicit::new(generator, DenseLu::default()).unwrap();
                bdf(
                    &mut implicit,
                    y,
                    controls,
                    &reports,
                    &mut Silent,
                    &mut work,
                    &interrupt,
                )
            }
        }
        .unwrap();
        out.push(y.clone());
        out
    };
    for bound in [40, 9] {
        let base = birth_death(b, k, 0);
        let t = Topology::enumerate(&spec(&base, vec![bound], vec![vec![0]]), &interrupt).unwrap();
        let build = |rate: f64| Generator::build(&t, &birth_death(rate, k, 0), &[], &interrupt);
        let (g, plus, minus) = (
            build(b).unwrap(),
            build(b + h).unwrap(),
            build(b - h).unwrap(),
        );
        let augmented = g.tangents(&[(&plus, &minus, 2.0 * h)]).unwrap();
        let n = g.dim();
        assert_eq!(augmented.dim(), 2 * n);
        let mut y = vec![0.0; 2 * n];
        y[0] = 1.0;
        let states = solve(&augmented, Some(n), &mut y);
        let mut start = vec![0.0; n];
        start[0] = 1.0;
        let (mut up, mut down) = (start.clone(), start);
        let ups = solve(&plus, None, &mut up);
        let downs = solve(&minus, None, &mut down);
        let mut sink = 0.0_f64;
        for (j, (state, at)) in states.iter().zip([0.5, 1.5, 3.0]).enumerate() {
            let (law, tangent) = state.split_at(n);
            let total: f64 = law.iter().sum();
            let drift: f64 = tangent.iter().sum();
            assert!((total - 1.0).abs() < 1e-12, "{bound} {at} {total}");
            assert!(drift.abs() < 1e-10, "{bound} {at} {drift}");
            sink = sink.max(law[n - 1]);
            let pairs = ups[j].iter().zip(&downs[j]);
            for (i, (&slope, (&above, &below))) in tangent.iter().zip(pairs).enumerate() {
                let difference = (above - below) / (2.0 * h);
                assert!((slope - difference).abs() < 1e-6, "{bound} {at} {i}");
            }
            if bound == 40 {
                let mean = b * (1.0 - (-k * at).exp()) / k;
                let slope = (1.0 - (-k * at).exp()) / k;
                for (i, &x) in t.states().iter().enumerate() {
                    let exact = poisson(x, mean) * (x as f64 / mean - 1.0) * slope;
                    assert!((tangent[i] - exact).abs() < 1e-8, "{at} {x}");
                }
            }
        }
        if bound == 9 {
            assert!(sink > 1e-3, "{sink}");
        }
    }
}

#[test]
fn tangents_and_their_solver_refuse_invalid_systems() {
    let v = birth_death(1.0, 1.0, 0);
    let interrupt = Interrupt::new(&never);
    let t = Topology::enumerate(&spec(&v, vec![4], vec![vec![0]]), &interrupt).unwrap();
    let g = Generator::build(&t, &v, &[], &interrupt).unwrap();
    let small = Topology::enumerate(&spec(&v, vec![2], vec![vec![0]]), &interrupt).unwrap();
    let other = Generator::build(&small, &v, &[], &interrupt).unwrap();
    let forced = Generator::build(&t, &v, &[vec![0]], &interrupt).unwrap();
    for pairs in [
        vec![],
        vec![(&g, &g, 0.0)],
        vec![(&g, &g, f64::NAN)],
        vec![(&g, &other, 1.0)],
        vec![(&forced, &g, 1.0)],
    ] {
        let error = g.tangents(&pairs).unwrap_err();
        assert_eq!(error.code, "CRN-E071", "{pairs:?}", pairs = pairs.len());
    }
    let seventeen = vec![(&g, &g, 1.0); 17];
    assert_eq!(g.tangents(&seventeen).unwrap_err().code, "CRN-E071");
    // The tangent solver accepts only the block-triangular shape.
    let mut solver = TangentSolver::new(g.dim(), DenseLu::default());
    assert!(solver.analyze(g.constant()).is_err());
    let mut wrong = TangentSolver::new(g.dim() + 1, DenseLu::default());
    let augmented = g.tangents(&[(&g, &g, 1.0)]).unwrap();
    assert!(wrong.analyze(augmented.constant()).is_err());
    let mut fresh = TangentSolver::new(g.dim(), DenseLu::default());
    assert!(fresh.factorize(augmented.constant()).is_err());
    assert!(fresh.solve_in_place(&mut vec![0.0; 2 * g.dim()]).is_err());
}

#[test]
fn a_retained_set_sends_omitted_states_to_the_sink() {
    let v = birth_death(1.0, 1.0, 1);
    let interrupt = Interrupt::new(&never);
    let s = spec(&v, vec![5], vec![vec![1]]);
    let full = Topology::enumerate(&s, &interrupt).unwrap();
    let all: std::collections::HashSet<Vec<u64>> = (0..=5).map(|x| vec![x]).collect();
    let same = Topology::enumerate_retained(&s, Some(&all), &interrupt).unwrap();
    assert_eq!(same.states(), full.states());
    assert_eq!(same.targets(), full.targets());
    // Omit 3: the walk stops at 2, and 2 -> 3 goes to the sink.
    let some: std::collections::HashSet<Vec<u64>> =
        [0, 1, 2, 4, 5].iter().map(|&x| vec![x]).collect();
    let cut = Topology::enumerate_retained(&s, Some(&some), &interrupt).unwrap();
    assert_eq!(cut.states(), [1, 2, 0]);
    let g = Generator::build(&cut, &v, &[], &interrupt).unwrap();
    let escapes = g.escapes();
    assert_eq!(escapes[1], 1.0, "{escapes:?}");
    assert_eq!(escapes[0] + escapes[2], 0.0, "{escapes:?}");
    for bad in [
        std::collections::HashSet::new(),
        [vec![6]].into_iter().collect(),
        [vec![1, 0]].into_iter().collect(),
    ] {
        let error = Topology::enumerate_retained(&s, Some(&bad), &interrupt).unwrap_err();
        assert_eq!(error.code, "CRN-E071");
    }
    let mut capped = s.clone();
    capped.max_states = 2;
    let error = Topology::enumerate_retained(&capped, Some(&all), &interrupt).unwrap_err();
    assert_eq!(error.code, "CRN-E071");
}
