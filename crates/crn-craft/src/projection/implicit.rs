//! A variable-step, variable-order implicit integrator for `p' = Q(t) p`:
//! the numerical differentiation formulas (NDF) of orders 1 to 5 of
//! Shampine and Reichelt, "The MATLAB ODE Suite", SIAM J. Sci. Comput. 18
//! (1997), in backward-difference form with a quasi-constant step size
//! (Byrne and Hindmarsh, ACM TOMS 1 (1975); Hairer and Wanner, "Solving
//! Ordinary Differential Equations II", section III.5 and IV.8).
//!
//! The problem is linear, so the corrector equation of a step is linear:
//! with the iteration matrix `I - c Q(t)` factorized at the step's own
//! time and `c`, one solve gives the corrector exactly (up to rounding).
use super::generator::{Csc, Generator, Merged};
use super::propagate::SignalValues;
use super::{Interrupt, SolverCounts, SparseLinearSolver, Work, budget, input, memory};
use crate::diagnostic::Diagnostic;

/// The highest order of the formulas.
pub const MAX_ORDER: usize = 5;

/// The NDF coefficients `kappa_k` of Shampine and Reichelt (1997), table 1;
/// `kappa_5 = 0` makes order 5 the BDF.
const KAPPA: [f64; MAX_ORDER + 1] = [0.0, -0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0];

/// `gamma_k = sum_{j=1}^{k} 1/j`.
fn gamma(k: usize) -> f64 {
    (1..=k).map(|j| 1.0 / j as f64).sum()
}

/// `alpha_k = (1 - kappa_k) gamma_k`: the corrector reads
/// `alpha_k (y - y0) + sum_{j=1}^{k} gamma_j D_j = h f(t, y)`.
fn alpha(k: usize) -> f64 {
    (1.0 - KAPPA[k]) * gamma(k)
}

/// The local error constant of order `k`: the error is about
/// `(kappa_k gamma_k + 1/(k+1)) (y - y0)`.
fn error_constant(k: usize) -> f64 {
    KAPPA[k] * gamma(k) + 1.0 / (k as f64 + 1.0)
}

/// The most simplified-Newton iterations with a factorization from an
/// earlier time before refactorizing at the step's own time.
const STALE_ITERATIONS: usize = 4;
/// Safety factor of a step-size proposal.
const SAFETY: f64 = 0.9;
/// Largest growth and smallest shrink of the step size in one change.
const MAX_GROWTH: f64 = 10.0;
const MIN_SHRINK: f64 = 0.2;
/// A proposed growth below this keeps the step (and its factorization).
const MIN_GROWTH: f64 = 1.2;

/// `C_j(s) = prod_{m=0}^{j-1} (s + m) / (m + 1)`: the weight of the `j`-th
/// backward difference in the interpolant `y(t_n + s h) = sum_j C_j(s) D_j`.
fn newton_weight(j: usize, s: f64) -> f64 {
    (0..j).fold(1.0, |w, m| w * (s + m as f64) / (m as f64 + 1.0))
}

/// The controls of one implicit integration.
#[derive(Debug, Clone, Copy)]
pub struct BdfControls {
    /// Integration start.
    pub start: f64,
    /// Integration end, `> start`.
    pub end: f64,
    /// Relative tolerance; the absolute tolerance is `tolerance * 1e-4`.
    pub tolerance: f64,
    /// Largest step, positive (may be infinite).
    pub max_step: f64,
    /// Largest bound, in bytes, of one numeric factorization.
    pub max_factor_bytes: u64,
}

/// The iteration matrix `I - c Q(t)` of a generator on a fixed pattern
/// (the union pattern of [`Merged`] with every diagonal entry), its sparse
/// solver after the symbolic analysis, and the linear-algebra counts. One
/// value serves every implicit integration with the same generator.
pub struct Implicit<S> {
    merged: Merged,
    solver: S,
    matrix: Csc,
    /// The [`Merged`] entry of each matrix entry, `usize::MAX` for an added
    /// diagonal.
    source: Vec<usize>,
    /// The position of each diagonal entry in `matrix`.
    diagonal: Vec<usize>,
    factor_bytes: u64,
    /// `c` and the signal values of the current factorization.
    factored: Option<(f64, Vec<f64>)>,
    counts: SolverCounts,
}

impl<S: SparseLinearSolver> Implicit<S> {
    /// Builds the pattern of `I - c Q` for `generator` and analyzes it.
    ///
    /// # Errors
    ///
    /// Any error of the analysis.
    pub fn new(generator: &Generator, mut solver: S) -> Result<Self, Diagnostic> {
        let merged = Merged::new(generator);
        let dim = merged.dim();
        let (row_ptr, col_idx) = (merged.row_ptr(), merged.col_idx());
        let mut entries: Vec<(usize, usize, usize)> = Vec::with_capacity(col_idx.len() + dim);
        for r in 0..dim {
            for (k, &c) in col_idx
                .iter()
                .enumerate()
                .take(row_ptr[r + 1])
                .skip(row_ptr[r])
            {
                entries.push((c, r, k));
            }
        }
        for (j, &k) in merged.diagonal_positions().iter().enumerate() {
            if k == usize::MAX {
                entries.push((j, j, usize::MAX));
            }
        }
        entries.sort_unstable_by_key(|&(c, r, _)| (c, r));
        let mut col_ptr = vec![0; dim + 1];
        let mut row_idx = Vec::with_capacity(entries.len());
        let mut source = Vec::with_capacity(entries.len());
        let mut diagonal = vec![0; dim];
        for &(c, r, k) in &entries {
            if r == c {
                diagonal[c] = row_idx.len();
            }
            col_ptr[c + 1] += 1;
            row_idx.push(r);
            source.push(k);
        }
        for j in 0..dim {
            col_ptr[j + 1] += col_ptr[j];
        }
        let values = vec![0.0; row_idx.len()];
        let matrix = Csc {
            dim,
            col_ptr,
            row_idx,
            values,
        };
        let factor_bytes = solver.analyze(&matrix)?;
        let counts = SolverCounts {
            symbolic: 1,
            ..SolverCounts::default()
        };
        Ok(Self {
            merged,
            solver,
            matrix,
            source,
            diagonal,
            factor_bytes,
            factored: None,
            counts,
        })
    }

    /// The linear-algebra work so far.
    #[must_use]
    pub fn counts(&self) -> SolverCounts {
        self.counts
    }

    /// The analysis bound of one numeric factorization, in bytes.
    #[must_use]
    pub fn factor_bytes(&self) -> u64 {
        self.factor_bytes
    }

    /// The solver.
    #[must_use]
    pub fn solver(&self) -> &S {
        &self.solver
    }
}

/// The state of one integration.
struct Run<'a, S> {
    implicit: &'a mut Implicit<S>,
    signals: &'a mut dyn SignalValues,
    work: &'a mut Work,
    interrupt: &'a Interrupt<'a>,
    controls: BdfControls,
    /// Signal values at the last evaluated time.
    signal_values: Vec<f64>,
    /// `Q` at the last evaluated time, on the [`Merged`] pattern.
    q: Vec<f64>,
    rtol: f64,
    atol: f64,
}

impl<S: SparseLinearSolver> Run<'_, S> {
    /// Evaluates the signals and `Q` at `t`; at the end of the interval the
    /// incoming value (the next double toward the start), so a
    /// discontinuity at the end is approached from the left.
    fn at(&mut self, t: f64) -> Result<(), Diagnostic> {
        let BdfControls { start, end, .. } = self.controls;
        let at = if t == end {
            libm::nextafter(end, start)
        } else {
            t
        };
        self.signals.evaluate(at, &mut self.signal_values)?;
        if self
            .signal_values
            .iter()
            .any(|&v| v < 0.0 || !v.is_finite())
        {
            return Err(input(
                "time-dependent reaction rate must be finite and nonnegative",
            ));
        }
        self.implicit
            .merged
            .values_at(&self.signal_values, &mut self.q);
        Ok(())
    }

    /// `out = Q x` with the last evaluated `Q`, each row summed over
    /// increasing columns. One work unit.
    fn product(&mut self, x: &[f64], out: &mut [f64]) -> Result<(), Diagnostic> {
        self.work.charge(1)?;
        self.implicit.counts.products += 1;
        let merged = &self.implicit.merged;
        let (row_ptr, col_idx) = (merged.row_ptr(), merged.col_idx());
        for (r, value) in out.iter_mut().enumerate() {
            let range = row_ptr[r]..row_ptr[r + 1];
            let mut sum = 0.0;
            for (&q, &j) in self.q[range.clone()].iter().zip(&col_idx[range]) {
                sum += q * x[j];
            }
            *value = sum;
        }
        Ok(())
    }

    /// Factorizes `I - c Q` with the last evaluated `Q`. One work unit.
    fn refactor(&mut self, c: f64) -> Result<(), Diagnostic> {
        self.work.charge(1)?;
        self.interrupt.check()?;
        let implicit = &mut *self.implicit;
        for (value, &k) in implicit.matrix.values.iter_mut().zip(&implicit.source) {
            *value = if k == usize::MAX { 0.0 } else { -c * self.q[k] };
        }
        for &k in &implicit.diagonal {
            implicit.matrix.values[k] += 1.0;
        }
        implicit.factored = None;
        implicit.solver.factorize(&implicit.matrix)?;
        implicit.counts.numeric += 1;
        implicit.factored = Some((c, self.signal_values.clone()));
        Ok(())
    }

    /// Solves with the current factorization, in place. One work unit.
    fn solve(&mut self, rhs: &mut [f64]) -> Result<(), Diagnostic> {
        self.work.charge(1)?;
        self.implicit.counts.solves += 1;
        self.implicit.solver.solve_in_place(rhs)
    }

    /// Solves the corrector `(I - c Q(t)) d = c Q(t) y0 - psi` for `d`,
    /// `Q(t)` evaluated by the caller. With a factorization at this `c`
    /// and these signal values, one solve is exact. With one at this `c`
    /// but other signal values (an earlier time of a run of equal steps),
    /// simplified Newton iterations `d += M^{-1} (c Q (y0 + d) - psi - d)`
    /// run until the estimated iteration error `rate / (1 - rate) |delta|`
    /// is at most `newton_tol` in the error norm, at most
    /// [`STALE_ITERATIONS`] times; a rate of 1 or more, or no convergence,
    /// refactorizes at the step's own time. Otherwise it refactorizes.
    #[allow(clippy::too_many_arguments)]
    fn correct(
        &mut self,
        c: f64,
        y0: &[f64],
        psi: &[f64],
        scale: &[f64],
        newton_tol: f64,
        d: &mut [f64],
        scratch: &mut [Vec<f64>; 2],
    ) -> Result<(), Diagnostic> {
        let reuse = match &self.implicit.factored {
            Some((fc, values)) if *fc == c => Some(*values == self.signal_values),
            _ => None,
        };
        if reuse == Some(false) {
            d.fill(0.0);
            let mut previous: Option<f64> = None;
            for _ in 0..STALE_ITERATIONS {
                let [sum, residual] = scratch;
                for ((s, &a), &b) in sum.iter_mut().zip(y0).zip(d.iter()) {
                    *s = a + b;
                }
                self.product(sum, residual)?;
                for ((r, &p), &x) in residual.iter_mut().zip(psi).zip(d.iter()) {
                    *r = c * *r - p - x;
                }
                self.solve(residual)?;
                for (x, &delta) in d.iter_mut().zip(residual.iter()) {
                    *x += delta;
                }
                let norm = rms(residual, scale);
                if norm == 0.0 {
                    return Ok(());
                }
                if let Some(before) = previous {
                    let rate = norm / before;
                    if rate >= 1.0 || !rate.is_finite() {
                        break;
                    }
                    if rate / (1.0 - rate) * norm <= newton_tol {
                        return Ok(());
                    }
                }
                previous = Some(norm);
            }
        }
        if reuse != Some(true) {
            self.refactor(c)?;
        }
        self.product(y0, d)?;
        for (x, &p) in d.iter_mut().zip(psi) {
            *x = c * *x - p;
        }
        self.solve(d)
    }
}

/// The root mean square of `x / scale`.
fn rms(x: &[f64], scale: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    let sum: f64 = x
        .iter()
        .zip(scale)
        .map(|(&v, &s)| {
            let r = v / s;
            r * r
        })
        .sum();
    (sum / x.len() as f64).sqrt()
}

/// Changes the step of the differences `d[0..=order]` by `ratio`: the
/// interpolant is evaluated at `t_n - m * ratio * h`, `m = 0..=order`, and
/// the backward differences of these values with the new step replace
/// `d[1..=order]`.
fn rescale(d: &mut [Vec<f64>], order: usize, ratio: f64, scratch: &mut [Vec<f64>]) {
    // T[i][j] = sum_m (-1)^m binom(i, m) C_j(-m ratio).
    let mut t = [[0.0; MAX_ORDER + 1]; MAX_ORDER + 1];
    for (i, row) in t.iter_mut().enumerate().take(order + 1).skip(1) {
        for (j, entry) in row.iter_mut().enumerate().take(order + 1).skip(1) {
            let mut binomial = 1.0;
            let mut sum = 0.0;
            for m in 0..=i {
                let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
                sum += sign * binomial * newton_weight(j, -(m as f64) * ratio);
                binomial = binomial * (i - m) as f64 / (m as f64 + 1.0);
            }
            *entry = sum;
        }
    }
    for i in 1..=order {
        let out = &mut scratch[i];
        out.fill(0.0);
        for (j, column) in d.iter().enumerate().take(order + 1).skip(1) {
            let w = t[i][j];
            if w != 0.0 {
                for (o, &x) in out.iter_mut().zip(column) {
                    *o += w * x;
                }
            }
        }
    }
    for i in 1..=order {
        std::mem::swap(&mut d[i], &mut scratch[i]);
    }
}

/// Integrates `p` from `controls.start` to `controls.end` with the
/// variable-order NDF of orders 1 to 5 and returns the law at each time of
/// `reports` (increasing, inside `(start, end)`); `p` becomes the law at
/// `end`. The error per step is the local error estimate
/// `(kappa_k gamma_k + 1/(k+1)) (y - y0)` in the root mean square of
/// `error / (atol + rtol max(|y_n|, |y_{n+1}|))`, with `rtol = tolerance`
/// and `atol = 1e-4 tolerance`, at most 1. The step size and the order
/// change only after `k + 1` steps of equal size: the order (`k - 1`, `k`
/// or `k + 1`) whose error estimate allows the largest step is taken, and
/// the step grows by at most 10 and only by 1.2 or more (a smaller gain
/// does not pay for the factorization of the new iteration matrix); a
/// rejected step shrinks by at most 5. Steps are at most `max_step` and end
/// exactly at `end`; reports between steps come from the interpolating
/// polynomial of the last step. Signal values at the end are taken just
/// before it.
///
/// The iteration matrix `I - c Q(t)` is factorized when `c = h / alpha_k`
/// changes (a new step size or order) and, during a run of equal steps,
/// only when simplified Newton iterations with the factorization of an
/// earlier time do not converge (see `Run::correct`).
///
/// Work: one unit for each product of `Q(t)` with a vector, each solve and
/// each numeric factorization.
///
/// # Errors
///
/// `CRN-E070` when the factorization bound exceeds
/// `controls.max_factor_bytes` (at `/projection/memory`), when the clock
/// cannot advance, for the work limit or an interruption; `CRN-E071` for a
/// negative or non-finite signal value; any error of `signals` or the
/// solver.
#[allow(clippy::too_many_lines, clippy::cast_precision_loss)]
pub fn bdf<S: SparseLinearSolver>(
    implicit: &mut Implicit<S>,
    p: &mut Vec<f64>,
    controls: BdfControls,
    reports: &[f64],
    signals: &mut dyn SignalValues,
    work: &mut Work,
    interrupt: &Interrupt<'_>,
) -> Result<Vec<Vec<f64>>, Diagnostic> {
    memory(implicit.factor_bytes, controls.max_factor_bytes)?;
    let BdfControls {
        start,
        end,
        tolerance,
        max_step,
        ..
    } = controls;
    let n = implicit.merged.dim();
    let signal_count = implicit.merged.signals();
    let mut run = Run {
        implicit,
        signals,
        work,
        interrupt,
        controls,
        signal_values: vec![0.0; signal_count],
        q: Vec::new(),
        rtol: tolerance,
        atol: tolerance * 1e-4,
    };
    let (rtol, atol) = (run.rtol, run.atol);
    let span = end - start;
    let zero = || vec![0.0; n];
    // Backward differences D_0 = y_n, D_j = nabla^j y_n, two more for the
    // order selection.
    let mut d: Vec<Vec<f64>> = (0..MAX_ORDER + 3).map(|_| zero()).collect();
    let mut scratch: Vec<Vec<f64>> = (0..=MAX_ORDER).map(|_| zero()).collect();
    let mut pair = [zero(), zero()];
    let (mut y0, mut psi, mut corr, mut scale, mut y_new) =
        (zero(), zero(), zero(), zero(), zero());
    d[0].clone_from(p);

    // Initial step (Hairer, Norsett and Wanner, vol. I, section II.4).
    run.at(start)?;
    let mut f0 = zero();
    run.product(&d[0], &mut f0)?;
    for (s, &y) in scale.iter_mut().zip(&d[0]) {
        *s = atol + rtol * y.abs();
    }
    let (d0, d1) = (rms(&d[0], &scale), rms(&f0, &scale));
    let mut h0 = if d0 < 1e-5 || d1 < 1e-5 {
        1e-6
    } else {
        0.01 * d0 / d1
    };
    h0 = h0.min(span);
    for ((y, &a), &f) in y_new.iter_mut().zip(&d[0]).zip(&f0) {
        *y = a + h0 * f;
    }
    run.at(start + h0)?;
    let mut f1 = zero();
    run.product(&y_new, &mut f1)?;
    for (a, &b) in f1.iter_mut().zip(&f0) {
        *a -= b;
    }
    let d2 = rms(&f1, &scale) / h0;
    let h1 = if d1.max(d2) <= 1e-15 {
        (h0 * 1e-3).max(1e-6)
    } else {
        (0.01 / d1.max(d2)).sqrt()
    };
    let mut h = (100.0 * h0).min(h1).min(max_step).min(span);
    drop(f1);
    for (x, &f) in d[1].iter_mut().zip(&f0) {
        *x = h * f;
    }
    drop(f0);

    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mut order = 1;
    let mut equal = 0;
    let mut t = start;
    let mut out = Vec::with_capacity(reports.len());
    let mut next_report = 0;
    while t < end {
        run.interrupt.check()?;
        // End exactly at `end`.
        let mut t_new = t + h;
        if t_new >= end || end - t <= h * (1.0 + 1e-9) {
            let last = end - t;
            if last != h {
                rescale(&mut d, order, last / h, &mut scratch);
                h = last;
                equal = 0;
            }
            t_new = end;
        }
        if t_new <= t || h <= 0.0 {
            return Err(budget("BDF step cannot advance the clock"));
        }
        let c = h / alpha(order);
        for (j, column) in d.iter().enumerate().take(order + 1) {
            if j == 0 {
                y0.clone_from(column);
            } else {
                for (a, &b) in y0.iter_mut().zip(column) {
                    *a += b;
                }
            }
        }
        psi.fill(0.0);
        for (j, column) in d.iter().enumerate().take(order + 1).skip(1) {
            let w = gamma(j) / alpha(order);
            for (a, &b) in psi.iter_mut().zip(column) {
                *a += w * b;
            }
        }
        for (s, &y) in scale.iter_mut().zip(&y0) {
            *s = atol + rtol * y.abs();
        }
        run.at(t_new)?;
        run.correct(c, &y0, &psi, &scale, newton_tol, &mut corr, &mut pair)?;
        for ((y, &a), &b) in y_new.iter_mut().zip(&y0).zip(&corr) {
            *y = a + b;
        }
        for ((s, &a), &b) in scale.iter_mut().zip(&d[0]).zip(&y_new) {
            *s = atol + rtol * a.abs().max(b.abs());
        }
        let e = error_constant(order);
        let error = e.abs() * rms(&corr, &scale);
        if error.is_nan() || error > 1.0 {
            run.implicit.counts.rejected += 1;
            let factor = if error.is_finite() {
                (SAFETY * error.powf(-1.0 / (order as f64 + 1.0))).max(MIN_SHRINK)
            } else {
                MIN_SHRINK
            };
            rescale(&mut d, order, factor, &mut scratch);
            h *= factor;
            equal = 0;
            continue;
        }
        // Accept: update the differences with the correction.
        run.implicit.counts.steps += 1;
        {
            let (low, high) = d.split_at_mut(order + 2);
            for ((top, &x), &below) in high[0].iter_mut().zip(&corr).zip(&low[order + 1]) {
                *top = x - below;
            }
        }
        d[order + 1].clone_from(&corr);
        for j in (0..=order).rev() {
            let (low, high) = d.split_at_mut(j + 1);
            for (a, &b) in low[j].iter_mut().zip(&high[0]) {
                *a += b;
            }
        }
        t = t_new;
        while next_report < reports.len() && reports[next_report] <= t {
            let at = reports[next_report];
            if at == t {
                out.push(d[0].clone());
            } else {
                let s = (at - t) / h;
                let mut value = d[0].clone();
                for (j, column) in d.iter().enumerate().take(order + 1).skip(1) {
                    let w = newton_weight(j, s);
                    for (a, &b) in value.iter_mut().zip(column) {
                        *a += w * b;
                    }
                }
                out.push(value);
            }
            next_report += 1;
        }
        equal += 1;
        if t >= end || equal < order + 1 {
            continue;
        }
        // Order and step selection after k + 1 equal steps.
        let ratio = |error: f64, k: usize| {
            if error == 0.0 {
                f64::INFINITY
            } else {
                error.powf(-1.0 / (k as f64 + 1.0))
            }
        };
        let mut best = (ratio(error, order), order);
        if order > 1 {
            let lower = error_constant(order - 1).abs() * rms(&d[order], &scale);
            let candidate = ratio(lower, order - 1);
            if candidate > best.0 {
                best = (candidate, order - 1);
            }
        }
        if order < MAX_ORDER {
            let higher = error_constant(order + 1).abs() * rms(&d[order + 2], &scale);
            let candidate = ratio(higher, order + 1);
            if candidate > best.0 {
                best = (candidate, order + 1);
            }
        }
        let (gain, new_order) = best;
        let factor = (SAFETY * gain).min(MAX_GROWTH).min(max_step / h);
        if new_order == order && factor < MIN_GROWTH {
            continue;
        }
        order = new_order;
        let factor = factor.max(MIN_SHRINK);
        if factor != 1.0 {
            rescale(&mut d, order, factor, &mut scratch);
            h *= factor;
        }
        equal = 0;
    }
    std::mem::swap(p, &mut d[0]);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A polynomial of degree at most `order` keeps its backward
    /// differences under a change of step.
    #[test]
    fn rescaling_differences_is_exact_for_polynomials() {
        let y = |t: f64| 1.0 + 2.0 * t - 0.5 * t * t + 0.25 * t * t * t;
        let order = 3;
        let differences = |h: f64| {
            let values: Vec<f64> = (0..=order).map(|m| y(-(m as f64) * h)).collect();
            let mut d = vec![values.clone()];
            let mut current = values;
            for _ in 0..order {
                current = current.windows(2).map(|w| w[0] - w[1]).collect();
                d.push(current.clone());
            }
            d.iter().map(|v| vec![v[0]]).collect::<Vec<_>>()
        };
        let mut d = differences(0.3);
        d.resize(MAX_ORDER + 3, vec![0.0]);
        let mut scratch = vec![vec![0.0]; MAX_ORDER + 1];
        rescale(&mut d, order, 0.5 / 0.3, &mut scratch);
        let expected = differences(0.5);
        for j in 0..=order {
            assert!((d[j][0] - expected[j][0]).abs() < 1e-12, "{j}");
        }
        // The interpolant reproduces the polynomial between the points.
        let s = -0.37;
        let value: f64 = (0..=order).map(|j| newton_weight(j, s) * d[j][0]).sum();
        assert!((value - y(s * 0.5)).abs() < 1e-12);
    }

    #[test]
    fn coefficients_match_the_published_table() {
        // Order 5 is the BDF: error constant 1/6 - ... = gamma_5 kappa_5 + 1/6.
        assert!((error_constant(5) - 1.0 / 6.0).abs() < 1e-15);
        assert!((alpha(1) - 1.185).abs() < 1e-15);
        assert!((gamma(3) - 11.0 / 6.0).abs() < 1e-15);
    }
}
