//! Propagation of a projected law: uniformization of a constant generator,
//! the adaptive exponential midpoint rule for time-dependent signals, and
//! count interventions on a probability vector.
use super::generator::{Generator, Merged};
use super::special::{binomial_pmf, poisson_tail};
use super::topology::Topology;
use super::{Interrupt, Operation, Work, budget, input};
use crate::diagnostic::Diagnostic;
use std::collections::HashMap;

/// A uniformized transition matrix `P = I + Q / lambda` in row-major form,
/// every diagonal entry present.
struct Transition {
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
}

impl Transition {
    /// `next = P term`, each row summed over increasing columns from 0.
    fn apply(&self, term: &[f64], next: &mut [f64]) {
        for (i, out) in next.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[k] * term[self.col_idx[k]];
            }
            *out = sum;
        }
    }
}

/// The largest exit rate of the constant part.
fn exit_rate(generator: &Generator) -> f64 {
    let q = generator.constant();
    (0..q.dim)
        .map(|j| -q.diagonal(j))
        .fold(f64::NEG_INFINITY, f64::max)
}

/// The number of uniformization pieces of `rate * duration`, each with a
/// Poisson mean of at most 16, checked against the remaining work.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn pieces(rate: f64, duration: f64, work: &Work) -> Result<u64, Diagnostic> {
    let pieces = (rate * duration / 16.0).ceil().max(1.0);
    let remaining = work.limit - work.count;
    // NaN and values beyond i64 fail the comparison below as too large.
    let representable = pieces < 9.0e18;
    if !representable || i128::from(pieces as i64) > remaining {
        return Err(budget("uniformization exceeds the work budget"));
    }
    Ok(pieces as u64)
}

/// Propagates `p` (states and sink) over `duration` with the constant part
/// of `generator` by positive uniformization: the interval splits into
/// `ceil(rate * duration / 16)` pieces of Poisson mean `mu`, and each piece
/// sums the series until the omitted Poisson weight `P(K > k)` is at most
/// `tolerance / pieces`. Each series term costs one work unit. Returns the
/// sum of the omitted weights, an upper bound of the mass the truncation
/// drops.
///
/// # Errors
///
/// `CRN-E070` when the pieces exceed the remaining work or the work limit
/// is reached (at `/projection/work`), and on interruption.
#[allow(clippy::cast_precision_loss)]
pub fn uniform(
    generator: &Generator,
    p: &mut [f64],
    duration: f64,
    tolerance: f64,
    work: &mut Work,
    interrupt: &Interrupt<'_>,
) -> Result<f64, Diagnostic> {
    let rate = exit_rate(generator);
    if rate == 0.0 {
        return Ok(0.0);
    }
    let pieces = pieces(rate, duration, work)?;
    let transition = uniformized(generator, rate);
    let mu = rate * duration / pieces as f64;
    let threshold = tolerance / pieces as f64;
    let mut term = vec![0.0; p.len()];
    let mut next = vec![0.0; p.len()];
    let mut tail = 0.0;
    for _ in 0..pieces {
        term.copy_from_slice(p);
        let mut weight = libm::exp(-mu);
        for (out, &x) in p.iter_mut().zip(&term) {
            *out = weight * x;
        }
        let mut k: u64 = 0;
        let remaining = loop {
            let remaining = poisson_tail(k + 1, mu);
            // NaN stops the series like a small tail.
            if remaining.partial_cmp(&threshold) != Some(std::cmp::Ordering::Greater) {
                break remaining;
            }
            work.charge(1)?;
            interrupt.check()?;
            k += 1;
            transition.apply(&term, &mut next);
            std::mem::swap(&mut term, &mut next);
            weight *= mu / k as f64;
            for (out, &x) in p.iter_mut().zip(&term) {
                *out += weight * x;
            }
        };
        tail += remaining;
    }
    Ok(tail)
}

/// `I + Q / rate` of the constant part, by rows.
fn uniformized(generator: &Generator, rate: f64) -> Transition {
    let q = generator.constant();
    let dim = q.dim;
    let mut counts = vec![0usize; dim + 1];
    for j in 0..dim {
        let mut has_diagonal = false;
        for k in q.col_ptr[j]..q.col_ptr[j + 1] {
            counts[q.row_idx[k] + 1] += 1;
            has_diagonal |= q.row_idx[k] == j;
        }
        if !has_diagonal {
            counts[j + 1] += 1;
        }
    }
    for i in 0..dim {
        counts[i + 1] += counts[i];
    }
    let row_ptr = counts.clone();
    let mut fill = counts;
    let nnz = row_ptr[dim];
    let mut col_idx = vec![0; nnz];
    let mut values = vec![0.0; nnz];
    // Columns in increasing order, so each row's columns are increasing.
    for j in 0..dim {
        let mut has_diagonal = false;
        let mut place = |row: usize, value: f64| {
            col_idx[fill[row]] = j;
            values[fill[row]] = value;
            fill[row] += 1;
        };
        for k in q.col_ptr[j]..q.col_ptr[j + 1] {
            let (row, scaled) = (q.row_idx[k], q.values[k] / rate);
            if row == j {
                has_diagonal = true;
                place(row, 1.0 + scaled);
            } else {
                place(row, scaled);
            }
        }
        if !has_diagonal {
            place(j, 1.0);
        }
    }
    Transition {
        row_ptr,
        col_idx,
        values,
    }
}

/// Signal values at a time, one per signal of the generator.
pub trait SignalValues {
    /// Writes the value of every signal at `time` into `out`.
    ///
    /// # Errors
    ///
    /// Any failure of the signal source; it stops the integration.
    fn evaluate(&mut self, time: f64, out: &mut [f64]) -> Result<(), Diagnostic>;
}

/// The controls of one exponential-midpoint integration.
#[derive(Debug, Clone, Copy)]
pub struct MidpointControls {
    /// Integration start.
    pub start: f64,
    /// Integration end, `> start`.
    pub end: f64,
    /// Total step-error goal over `[start, end]`.
    pub tolerance: f64,
    /// Largest step, positive (may be infinite).
    pub max_step: f64,
}

/// The share of a step's error goal that the Poisson truncation of each
/// matrix-exponential action may use.
pub const MIDPOINT_TRUNCATION_SHARE: f64 = 1.0 / 1024.0;

struct Midpoint<'a> {
    merged: &'a Merged,
    controls: MidpointControls,
    signals: &'a mut dyn SignalValues,
    work: &'a mut Work,
    interrupt: &'a Interrupt<'a>,
    signal_values: Vec<f64>,
    q: Vec<f64>,
    transition: Vec<f64>,
    term: Vec<f64>,
    next: Vec<f64>,
}

impl Midpoint<'_> {
    /// Signal values at `t`; at the end of the interval the incoming value
    /// (the next double toward the start), so a discontinuity at the end
    /// is approached from the left.
    fn weights(&mut self, t: f64) -> Result<(), Diagnostic> {
        let MidpointControls { start, end, .. } = self.controls;
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
        Ok(())
    }

    /// `exp(h Q(at + h/2)) v` by uniformization of the frozen generator,
    /// with the omitted Poisson weight at most `share * goal` per action.
    #[allow(clippy::cast_precision_loss)]
    fn advance(&mut self, at: f64, h: f64, v: &[f64], goal: f64) -> Result<Vec<f64>, Diagnostic> {
        self.work.charge(1)?;
        self.weights(at + h / 2.0)?;
        let merged = self.merged;
        merged.values_at(&self.signal_values, &mut self.q);
        let diagonal = merged.diagonal_positions();
        let rate = diagonal
            .iter()
            .map(|&k| if k == usize::MAX { 0.0 } else { -self.q[k] })
            .fold(f64::NEG_INFINITY, f64::max);
        let mut p = v.to_vec();
        if rate == 0.0 {
            return Ok(p);
        }
        let pieces = pieces(rate, h, self.work)?;
        let (row_ptr, col_idx) = (merged.row_ptr(), merged.col_idx());
        self.transition.clear();
        self.transition.extend(self.q.iter().map(|&x| x / rate));
        for &k in diagonal {
            if k != usize::MAX {
                self.transition[k] += 1.0;
            }
        }
        let mu = rate * h / pieces as f64;
        let threshold = goal * MIDPOINT_TRUNCATION_SHARE / pieces as f64;
        for _ in 0..pieces {
            self.term.clone_from(&p);
            let mut weight = libm::exp(-mu);
            for (out, &x) in p.iter_mut().zip(&self.term) {
                *out = weight * x;
            }
            let mut k: u64 = 0;
            while poisson_tail(k + 1, mu) > threshold {
                self.interrupt.check()?;
                k += 1;
                for (i, out) in self.next.iter_mut().enumerate() {
                    let mut sum = if diagonal[i] == usize::MAX {
                        self.term[i]
                    } else {
                        0.0
                    };
                    let range = row_ptr[i]..row_ptr[i + 1];
                    for (&value, &j) in self.transition[range.clone()].iter().zip(&col_idx[range]) {
                        sum += value * self.term[j];
                    }
                    *out = sum;
                }
                std::mem::swap(&mut self.term, &mut self.next);
                weight *= mu / k as f64;
                for (out, &x) in p.iter_mut().zip(&self.term) {
                    *out += weight * x;
                }
            }
        }
        Ok(p)
    }
}

/// Integrates `p` from `controls.start` to `controls.end` with the adaptive
/// exponential midpoint rule: a step `h` from `t` is `exp(h Q(t + h/2))`;
/// the L1 gap between one full step and two half steps must be at most
/// `tolerance * h / (end - start)`; the step then grows or shrinks by
/// `clip(0.9 (goal / gap)^(1/3), 0.2, 2)`, at most `max_step`. Each
/// matrix-exponential action costs one work unit and is computed by
/// uniformization of the frozen generator, truncated where the omitted
/// Poisson weight is at most [`MIDPOINT_TRUNCATION_SHARE`] of the step goal
/// (so the accepted steps drop at most `tolerance / 512` of mass in all).
/// Signal values at the end of the interval are taken just before it.
///
/// # Errors
///
/// `CRN-E070` when the clock cannot advance, when the goal is below the
/// rounding level of the error estimate, for the work limit or an
/// interruption; `CRN-E071` for a negative or non-finite signal value;
/// any error of `signals`.
pub fn midpoint(
    generator: &Generator,
    merged: &Merged,
    p: &mut Vec<f64>,
    controls: MidpointControls,
    signals: &mut dyn SignalValues,
    work: &mut Work,
    interrupt: &Interrupt<'_>,
) -> Result<(), Diagnostic> {
    let MidpointControls {
        start,
        end,
        tolerance,
        max_step,
    } = controls;
    let dim = merged.dim();
    let mut state = Midpoint {
        merged,
        controls,
        signals,
        work,
        interrupt,
        signal_values: vec![0.0; generator.signals().len()],
        q: Vec::new(),
        transition: Vec::new(),
        term: vec![0.0; dim],
        next: vec![0.0; dim],
    };
    let mut t = start;
    let mut h = (end - start).min(max_step);
    while t < end {
        h = h.min(end - t);
        if t + h <= t {
            return Err(budget("midpoint step cannot advance the clock"));
        }
        let goal = tolerance * h / (end - start);
        let full = state.advance(t, h, p, goal)?;
        let first = state.advance(t, h / 2.0, p, goal)?;
        let half = state.advance(t + h / 2.0, h / 2.0, &first, goal)?;
        let gap: f64 = full.iter().zip(&half).map(|(a, b)| (a - b).abs()).sum();
        if gap <= goal {
            *p = half;
            t += h;
        } else if goal <= 64.0 * f64::EPSILON {
            // The step error estimate cannot resolve this goal in floating point.
            return Err(budget(
                "midpoint tolerance is below the rounding level of its step \
                 error estimate; use a looser tolerance or solver='bdf'",
            ));
        }
        let factor = if gap == 0.0 {
            2.0
        } else {
            (0.9 * libm::pow(goal / gap, 1.0 / 3.0)).clamp(0.2, 2.0)
        };
        h = max_step.min(h * factor);
    }
    Ok(())
}

/// Applies the count operations of one intervention to a law over the
/// states and the sink: `set` and `change` move each state's mass, a
/// partition with `0 < p < 1` spreads it with exact binomial weights; mass
/// above a bound moves to the sink. Operations apply one species at a time
/// in the given order. Each state costs one work unit per outcome, charged
/// before the outcomes are formed.
///
/// # Errors
///
/// `CRN-E071` for a decrement below zero or an outcome outside the prepared
/// states; `CRN-E070` when the intermediate law exceeds the state budget,
/// for the work limit or an interruption.
#[allow(clippy::cast_precision_loss)]
pub fn kick(
    topology: &Topology,
    p: &[f64],
    changes: &[(usize, Operation)],
    work: &mut Work,
    interrupt: &Interrupt<'_>,
) -> Result<Vec<f64>, Diagnostic> {
    let n = topology.len();
    let mut mass: Vec<(Vec<u64>, f64)> = (0..n)
        .filter(|&j| p[j] != 0.0)
        .map(|j| (topology.state(j).to_vec(), p[j]))
        .collect();
    let mut escape = p[n];
    for &(i, operation) in changes {
        let mut next: Vec<(Vec<u64>, f64)> = Vec::new();
        let mut positions: HashMap<Vec<u64>, usize> = HashMap::new();
        let bound = i128::from(topology.bounds()[i]);
        for (counts, amount) in &mass {
            interrupt.check()?;
            let count = counts[i];
            let (first, last, binomial) = match operation {
                Operation::Partition(q) if q > 0.0 && q < 1.0 => (0, i128::from(count), true),
                Operation::Partition(q) => {
                    let n = if q == 0.0 { 0 } else { i128::from(count) };
                    (n, n, false)
                }
                Operation::Set(n) => (i128::from(n), i128::from(n), false),
                Operation::Change(d) => {
                    let n = i128::from(count) + i128::from(d);
                    (n, n, false)
                }
            };
            work.charge((last - first + 1).unsigned_abs())?;
            if first < 0 {
                return Err(input(format!(
                    "intervention removes more {} than available",
                    topology.species()[i]
                )));
            }
            let mut place = |outcome: i128, weight: f64| -> Result<(), Diagnostic> {
                if outcome > bound {
                    escape += amount * weight;
                } else if weight != 0.0 {
                    let mut target = counts.clone();
                    // Inside the bound the outcome fits in u64.
                    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                    {
                        target[i] = outcome as u64;
                    }
                    if let Some(&at) = positions.get(&target) {
                        next[at].1 += amount * weight;
                    } else {
                        positions.insert(target.clone(), next.len());
                        next.push((target, 0.0 + amount * weight));
                        if next.len() > topology.max_states() {
                            return Err(budget("intervention workspace budget exceeded"));
                        }
                    }
                }
                Ok(())
            };
            if binomial && count > 0 {
                let Operation::Partition(q) = operation else {
                    unreachable!("binomial outcomes come from a partition")
                };
                let trials = count as f64;
                for x in 0..=count {
                    place(i128::from(x), binomial_pmf(x as f64, trials, q, 1.0 - q))?;
                }
            } else {
                place(first, 1.0)?;
            }
        }
        mass = next;
    }
    let mut out = vec![0.0; n + 1];
    out[n] = escape;
    for (counts, value) in mass {
        if value != 0.0 {
            match topology.index_of(&counts) {
                Some(j) => out[j] = value,
                None => return Err(input("intervention reaches an unprepared state")),
            }
        }
    }
    Ok(out)
}
