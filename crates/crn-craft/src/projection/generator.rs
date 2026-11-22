//! Generator assembly: the constant part, one matrix per time-dependent
//! signal, and the stop states.
use super::topology::Topology;
use super::{Interrupt, budget, input};
use crate::diagnostic::Diagnostic;
use crate::rules::FlatView;

/// Most entries of the state-dependent rate table: states times
/// state-dependent channels (at least one).
pub const MAX_RATE_TABLE: usize = 10_000_000;

/// A square sparse matrix in compressed sparse column form, rows sorted
/// within each column. This is the form a sparse direct solver consumes.
#[derive(Debug, Clone, PartialEq)]
pub struct Csc {
    /// Rows and columns.
    pub dim: usize,
    /// Column `j` holds entries `col_ptr[j]..col_ptr[j + 1]`.
    pub col_ptr: Vec<usize>,
    /// Row of each entry.
    pub row_idx: Vec<usize>,
    /// Value of each entry.
    pub values: Vec<f64>,
}

impl Csc {
    /// The generator part of the edges `selected` with the given rates: for
    /// edge `e` from `s` to `d`, `+rate` at `(d, s)` and `-rate` at
    /// `(s, s)`, duplicates summed in edge order. The sink column is empty.
    fn from_edges(topology: &Topology, rates: &[f64], selected: &dyn Fn(usize) -> bool) -> Self {
        let dim = topology.len() + 1;
        let (src, dst) = (topology.sources(), topology.targets());
        let mut col_ptr = Vec::with_capacity(dim + 1);
        let mut row_idx = Vec::new();
        let mut values = Vec::new();
        col_ptr.push(0);
        let mut e = 0;
        let mut column: Vec<(usize, usize, f64)> = Vec::new();
        for s in 0..dim {
            column.clear();
            let mut diagonal: Option<f64> = None;
            // Edges are grouped by source in increasing order.
            while e < src.len() && src[e] == s {
                if selected(e) {
                    column.push((dst[e], e, rates[e]));
                    diagonal = Some(diagonal.map_or(-rates[e], |d| d - rates[e]));
                }
                e += 1;
            }
            if let Some(d) = diagonal {
                column.push((s, usize::MAX, d));
            }
            // Sort by row; equal rows keep edge order, so sums follow it.
            column.sort_by_key(|&(row, order, _)| (row, order));
            for &(row, _, value) in &column {
                if row_idx.len() > col_ptr[col_ptr.len() - 1] && row_idx.last() == Some(&row) {
                    *values.last_mut().expect("entry") += value;
                } else {
                    row_idx.push(row);
                    values.push(value);
                }
            }
            col_ptr.push(row_idx.len());
        }
        Self {
            dim,
            col_ptr,
            row_idx,
            values,
        }
    }

    /// The number of stored entries.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// The diagonal entry of column `j` (0 when absent).
    #[must_use]
    pub fn diagonal(&self, j: usize) -> f64 {
        let range = self.col_ptr[j]..self.col_ptr[j + 1];
        match self.row_idx[range.clone()].binary_search(&j) {
            Ok(k) => self.values[range.start + k],
            Err(_) => 0.0,
        }
    }
}

/// The generator of a prepared projection under one parameter binding:
/// `dp/dt = (constant + sum_s v_s(t) signal_s) p` over the states and the
/// sink (index `n`, last). Edges out of a stop state are absent.
#[derive(Debug, Clone)]
pub struct Generator {
    constant: Csc,
    signals: Vec<Csc>,
    stopped: Vec<bool>,
}

impl Generator {
    /// Assembles the generator of `topology` with the rates of `view` (the
    /// model lowered under the binding; the caller checks that its channels
    /// match the prepared ones). `signals[k]` lists the channels driven by
    /// signal `k`; their rate constant is replaced by the signal value.
    ///
    /// # Errors
    ///
    /// `CRN-E070` when the state-dependent rate table exceeds
    /// [`MAX_RATE_TABLE`]; the rate evaluator's diagnostics (`CRN-R030`);
    /// `CRN-E071` at `/signals/k` for a signal that names no channel or one
    /// with a count-dependent rate, and for a negative or non-finite rate.
    pub fn build(
        topology: &Topology,
        view: &FlatView,
        signals: &[Vec<usize>],
        interrupt: &Interrupt<'_>,
    ) -> Result<Self, Diagnostic> {
        let n = topology.len();
        let channels = view.channels();
        let constants: Vec<Option<f64>> = (0..channels)
            .map(|c| view.rate_constant(view.rule(c)))
            .collect();
        let factors: Vec<f64> = (0..channels).map(|c| view.factor(c)).collect();
        let (src, edge_channel) = (topology.sources(), topology.channels());
        let basis: Vec<f64> = topology
            .multiplicities()
            .iter()
            .zip(edge_channel)
            .map(|(&m, &c)| m * factors[c])
            .collect();
        let mut rates: Vec<f64> = basis
            .iter()
            .zip(edge_channel)
            .map(|(&b, &c)| constants[c].unwrap_or(0.0) * b)
            .collect();
        let variable: Vec<usize> = (0..channels).filter(|&c| constants[c].is_none()).collect();
        let mut stopped = vec![false; n];
        if !variable.is_empty() || view.has_stop() {
            if n.saturating_mul(variable.len().max(1)) > MAX_RATE_TABLE {
                return Err(budget(
                    "state-dependent rate table exceeds the workspace budget",
                ));
            }
            let mut slot = vec![usize::MAX; channels];
            for (i, &c) in variable.iter().enumerate() {
                slot[c] = i;
            }
            let mut evaluator = view.evaluator();
            let mut values = vec![0.0; variable.len()];
            let mut e = 0;
            for (s, stop) in stopped.iter_mut().enumerate() {
                if s % 1024 == 0 {
                    interrupt.check()?;
                }
                evaluator.load(topology.state(s))?;
                *stop = evaluator.stops()?;
                for (value, &c) in values.iter_mut().zip(&variable) {
                    *value = if *stop { 0.0 } else { evaluator.channel(c)? };
                }
                while e < src.len() && src[e] == s {
                    let i = slot[edge_channel[e]];
                    if i != usize::MAX {
                        rates[e] = values[i];
                    }
                    e += 1;
                }
            }
        }
        let mut matrices = Vec::with_capacity(signals.len());
        for (k, ids) in signals.iter().enumerate() {
            if ids.is_empty() || ids.iter().any(|&c| c >= channels || constants[c].is_none()) {
                return Err(Diagnostic::error(
                    "CRN-E071",
                    &format!("/signals/{k}"),
                    "the signal names no rule with a count-independent rate",
                    "Drive only rules whose rate depends on parameters alone.",
                ));
            }
            let mut driven = vec![false; channels];
            for &c in ids {
                driven[c] = true;
            }
            matrices.push(Csc::from_edges(topology, &basis, &|e| {
                driven[edge_channel[e]] && !stopped[src[e]]
            }));
            for (rate, &c) in rates.iter_mut().zip(edge_channel) {
                if driven[c] {
                    *rate = 0.0;
                }
            }
        }
        for (rate, &s) in rates.iter_mut().zip(src) {
            if stopped[s] {
                *rate = 0.0;
            }
        }
        if rates.iter().any(|&r| r < 0.0 || !r.is_finite()) {
            return Err(input("invalid rate in the projection"));
        }
        let constant = Csc::from_edges(topology, &rates, &|e| rates[e] > 0.0);
        Ok(Self {
            constant,
            signals: matrices,
            stopped,
        })
    }

    /// A generator without signals from `(source, target, rate)` triples
    /// on states `0..n` and the sink `n` (test fixtures).
    #[cfg(test)]
    pub(crate) fn from_rates(n: usize, rates: &[(usize, usize, f64)]) -> Self {
        let mut columns: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n + 1];
        for &(s, d, r) in rates {
            columns[s].push((d, r));
            columns[s].push((s, -r));
        }
        let mut col_ptr = vec![0];
        let (mut row_idx, mut values) = (Vec::new(), Vec::new());
        for mut column in columns {
            column.sort_by_key(|&(r, _)| r);
            for (r, v) in column {
                if row_idx.len() > col_ptr[col_ptr.len() - 1] && row_idx.last() == Some(&r) {
                    *values.last_mut().expect("entry") += v;
                } else {
                    row_idx.push(r);
                    values.push(v);
                }
            }
            col_ptr.push(row_idx.len());
        }
        Self {
            constant: Csc {
                dim: n + 1,
                col_ptr,
                row_idx,
                values,
            },
            signals: Vec::new(),
            stopped: vec![false; n],
        }
    }

    /// The count-independent and state-dependent rates, signals excluded.
    #[must_use]
    pub fn constant(&self) -> &Csc {
        &self.constant
    }

    /// The unit-signal matrix of signal `k`.
    #[must_use]
    pub fn signals(&self) -> &[Csc] {
        &self.signals
    }

    /// Whether the stop condition holds in each state.
    #[must_use]
    pub fn stopped(&self) -> &[bool] {
        &self.stopped
    }

    /// The dimension: states plus the sink.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.constant.dim
    }

    /// The rate of each state into the sink (row `n` of the constant
    /// part), states only.
    #[must_use]
    pub fn escapes(&self) -> Vec<f64> {
        let q = &self.constant;
        let n = q.dim - 1;
        (0..n)
            .map(|j| {
                let end = q.col_ptr[j + 1];
                if end > q.col_ptr[j] && q.row_idx[end - 1] == n {
                    q.values[end - 1]
                } else {
                    0.0
                }
            })
            .collect()
    }
}

/// The generator at one time on a fixed pattern: the union of the entries
/// of the constant part and every signal matrix, in row-major form for the
/// gathering product. The pattern never changes between times, so a sparse
/// factorization of `I - h Q(t)` can reuse one symbolic analysis.
#[derive(Debug, Clone)]
pub struct Merged {
    dim: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    /// Per entry: constant value (if present), then per signal (if present).
    parts: Vec<Option<f64>>,
    signals: usize,
    diagonal: Vec<usize>,
}

impl Merged {
    /// The union pattern of a generator.
    #[must_use]
    pub fn new(generator: &Generator) -> Self {
        let dim = generator.dim();
        let signals = generator.signals.len();
        let stride = signals + 1;
        let mut triplets: Vec<(usize, usize, usize, f64)> = Vec::new();
        for (m, matrix) in std::iter::once(&generator.constant)
            .chain(&generator.signals)
            .enumerate()
        {
            for j in 0..dim {
                for k in matrix.col_ptr[j]..matrix.col_ptr[j + 1] {
                    triplets.push((matrix.row_idx[k], j, m, matrix.values[k]));
                }
            }
        }
        triplets.sort_by_key(|&(r, c, m, _)| (r, c, m));
        let mut row_ptr = vec![0; dim + 1];
        let mut col_idx = Vec::new();
        let mut parts: Vec<Option<f64>> = Vec::new();
        let mut last = None;
        for &(r, c, m, v) in &triplets {
            if last != Some((r, c)) {
                last = Some((r, c));
                col_idx.push(c);
                parts.extend(std::iter::repeat_n(None, stride));
                row_ptr[r + 1] += 1;
            }
            let at = (col_idx.len() - 1) * stride + m;
            parts[at] = Some(v);
        }
        for r in 0..dim {
            row_ptr[r + 1] += row_ptr[r];
        }
        let mut diagonal = vec![usize::MAX; dim];
        for (r, slot) in diagonal.iter_mut().enumerate() {
            let range = row_ptr[r]..row_ptr[r + 1];
            if let Some(k) = col_idx[range.clone()].iter().position(|&c| c == r) {
                *slot = range.start + k;
            }
        }
        Self {
            dim,
            row_ptr,
            col_idx,
            parts,
            signals,
            diagonal,
        }
    }

    /// The dimension.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// The number of signals.
    #[must_use]
    pub fn signals(&self) -> usize {
        self.signals
    }

    /// Row pointers of the union pattern.
    #[must_use]
    pub fn row_ptr(&self) -> &[usize] {
        &self.row_ptr
    }

    /// Column of each entry of the union pattern.
    #[must_use]
    pub fn col_idx(&self) -> &[usize] {
        &self.col_idx
    }

    /// Writes `Q(t) = constant + sum_s v_s signal_s` on the pattern, summed
    /// in that order (a part absent from an entry adds nothing).
    pub fn values_at(&self, signal_values: &[f64], out: &mut Vec<f64>) {
        let stride = self.signals + 1;
        out.clear();
        out.extend(self.parts.chunks_exact(stride).map(|parts| {
            let mut value: Option<f64> = parts[0];
            for (part, &v) in parts[1..].iter().zip(signal_values) {
                if let Some(m) = part {
                    let term = v * m;
                    value = Some(value.map_or(term, |x| x + term));
                }
            }
            value.unwrap_or(0.0)
        }));
    }

    /// The position of the diagonal entry of each row, `usize::MAX` if absent.
    #[must_use]
    pub fn diagonal_positions(&self) -> &[usize] {
        &self.diagonal
    }
}
