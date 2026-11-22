//! Constant-segment tangent equations and a shared sparse diagonal solve.
use super::{Csc, Generator, SparseLinearSolver, budget, input};
use crate::diagnostic::Diagnostic;
use std::collections::BTreeMap;

impl Generator {
    /// Build a triangular generator for [p, dp/dtheta_1, ...] from base,
    /// plus and minus generator bindings on the same prepared graph. Each
    /// third entry is the distance between those bindings in the requested
    /// parameter coordinate.
    ///
    /// # Errors
    /// Dimensions or finite-difference steps that do not match, reachable
    /// stop states, forced signals or an augmented problem above its cap.
    pub fn tangents(&self, pairs: &[(&Self, &Self, f64)]) -> Result<Self, Diagnostic> {
        let base = self.constant();
        let n = base.dim;
        let compatible =
            |g: &Self| g.dim() == n && g.signals().is_empty() && !g.stopped().iter().any(|&x| x);
        if pairs.is_empty()
            || pairs.len() > 16
            || !compatible(self)
            || pairs.iter().any(|(plus, minus, width)| {
                !compatible(plus) || !compatible(minus) || !width.is_finite() || *width <= 0.0
            })
        {
            return Err(input(
                "invalid generator bindings for projection sensitivities",
            ));
        }
        let dim = n * (1 + pairs.len());
        let entries = (1 + pairs.len()).saturating_mul(base.nnz()).saturating_add(
            pairs
                .iter()
                .map(|(a, b, _)| a.constant().nnz() + b.constant().nnz())
                .sum(),
        );
        if dim > 250_000 || entries > 2_000_000 {
            return Err(budget(
                "augmented sensitivity generator exceeds the workspace limit",
            ));
        }
        let mut col_ptr = Vec::with_capacity(dim + 1);
        let (mut row_idx, mut values) = (Vec::with_capacity(entries), Vec::with_capacity(entries));
        col_ptr.push(0);
        for block in 0..=pairs.len() {
            for c in 0..n {
                let range = base.col_ptr[c]..base.col_ptr[c + 1];
                row_idx.extend(base.row_idx[range.clone()].iter().map(|r| r + block * n));
                values.extend_from_slice(&base.values[range]);
                if block == 0 {
                    for (j, &(a, b, width)) in pairs.iter().enumerate() {
                        let mut derivative = BTreeMap::new();
                        for (matrix, sign) in [(a.constant(), 1.), (b.constant(), -1.)] {
                            for k in matrix.col_ptr[c]..matrix.col_ptr[c + 1] {
                                *derivative.entry(matrix.row_idx[k]).or_insert(0.) +=
                                    sign * matrix.values[k];
                            }
                        }
                        for (row, numerator) in derivative.into_iter().filter(|(_, v)| *v != 0.) {
                            let value = numerator / width;
                            if !value.is_finite() {
                                return Err(input("nonfinite generator derivative"));
                            }
                            row_idx.push((j + 1) * n + row);
                            values.push(value);
                        }
                    }
                }
                col_ptr.push(row_idx.len());
            }
        }
        Ok(Self::from_augmented(Csc {
            dim,
            col_ptr,
            row_idx,
            values,
        }))
    }
}

/// Back-substitution for the triangular iteration matrix. One factorization
/// of the diagonal block serves the law and every sensitivity direction.
pub struct TangentSolver<S> {
    size: usize,
    shape: (usize, usize),
    solver: S,
    diagonal: Option<Csc>,
    entries: Vec<usize>,
    coupling: Vec<(usize, usize, usize, f64)>,
    factored: bool,
}

impl<S: SparseLinearSolver> TangentSolver<S> {
    /// Wrap the sparse solver for one probability block of the given size.
    pub fn new(size: usize, solver: S) -> Self {
        Self {
            size,
            shape: (0, 0),
            solver,
            diagonal: None,
            entries: Vec::new(),
            coupling: Vec::new(),
            factored: false,
        }
    }
}

impl<S: SparseLinearSolver> SparseLinearSolver for TangentSolver<S> {
    fn analyze(&mut self, pattern: &Csc) -> Result<u64, Diagnostic> {
        self.factored = false;
        self.entries.clear();
        self.coupling.clear();
        let n = self.size;
        if n == 0 || pattern.dim <= n || pattern.dim % n != 0 {
            return Err(input("invalid tangent system size"));
        }
        let mut col_ptr = vec![0];
        let mut row_idx = Vec::new();
        for c in 0..n {
            for (index, &row) in pattern.row_idx[pattern.col_ptr[c]..pattern.col_ptr[c + 1]]
                .iter()
                .enumerate()
            {
                let k = pattern.col_ptr[c] + index;
                if row < n {
                    self.entries.push(k);
                    row_idx.push(row);
                } else {
                    self.coupling.push((c, row - n, k, 0.));
                }
            }
            col_ptr.push(row_idx.len());
        }
        // Accept only the triangular structure with identical diagonal
        // blocks assembled by Generator::tangents.
        for c in n..pattern.dim {
            let column = c % n;
            let offset = c - column;
            let actual = pattern.row_idx[pattern.col_ptr[c]..pattern.col_ptr[c + 1]]
                .iter()
                .copied();
            let expected = row_idx[col_ptr[column]..col_ptr[column + 1]]
                .iter()
                .map(|r| r + offset);
            if !actual.eq(expected) {
                return Err(input("invalid diagonal sparsity in tangent system"));
            }
        }
        let values = vec![0.; row_idx.len()];
        let base = Csc {
            dim: n,
            col_ptr,
            row_idx,
            values,
        };
        let bound = self.solver.analyze(&base)?;
        self.diagonal = Some(base);
        self.shape = (pattern.dim, pattern.nnz());
        Ok(bound.saturating_add(32 * self.coupling.len() as u64 + 24 * self.entries.len() as u64))
    }

    fn factorize(&mut self, matrix: &Csc) -> Result<(), Diagnostic> {
        self.factored = false;
        if (matrix.dim, matrix.values.len()) != self.shape {
            return Err(input("tangent factorization size differs from analysis"));
        }
        let Some(base) = self.diagonal.as_mut() else {
            return Err(input("tangent factorization requires symbolic analysis"));
        };
        for (value, &k) in base.values.iter_mut().zip(&self.entries) {
            *value = matrix.values[k];
        }
        for (_, _, k, value) in &mut self.coupling {
            *value = matrix.values[*k];
        }
        self.solver.factorize(base)?;
        self.factored = true;
        Ok(())
    }

    fn solve_in_place(&mut self, rhs: &mut [f64]) -> Result<(), Diagnostic> {
        if !self.factored || rhs.len() != self.shape.0 {
            return Err(input("invalid factorization or tangent right-hand side"));
        }
        let (law, sensitivities) = rhs.split_at_mut(self.size);
        self.solver.solve_in_place(law)?;
        for &(column, row, _, value) in &self.coupling {
            sensitivities[row] -= value * law[column];
        }
        for derivative in sensitivities.chunks_exact_mut(self.size) {
            self.solver.solve_in_place(derivative)?;
        }
        Ok(())
    }
}
