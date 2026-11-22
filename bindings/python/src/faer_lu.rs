//! The sparse LU of faer behind `crn_craft::projection::SparseLinearSolver`.
//!
//! One symbolic analysis (COLAMD column ordering, the column elimination
//! tree of `A^T A`) per pattern serves every numeric factorization with
//! that pattern; workspaces are kept between factorizations and solves.
//! Everything runs single-threaded (`Par::Seq`; faer is built without
//! rayon).
//!
//! The memory bound of a numeric factorization comes from the same
//! analysis: with partial pivoting, the pattern of `U` lies in that of the
//! Cholesky factor `R` of `A^T A` and the pattern of `L` in that of the
//! Householder vectors `H` of a QR factorization of `A` (George and Ng,
//! SIAM J. Sci. Stat. Comput. 8 (1987)), whose column counts the analysis
//! gives. The bound is `16 (nnz(H) + nnz(R))` bytes (a value and an index
//! per entry) plus the numeric workspace and the permutations.
use crn_craft::diagnostic::Diagnostic;
use crn_craft::projection::{Csc, SparseLinearSolver};
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::sparse::linalg::{lu, qr};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Conj, MatMut, Par};

fn failure(message: String) -> Diagnostic {
    Diagnostic::error(
        "CRN-E070",
        "/projection/solver",
        message,
        "Lower the bounds or raise max_factor_bytes within the memory of the machine.",
    )
}

fn workspace(req: StackReq) -> Result<MemBuffer, Diagnostic> {
    MemBuffer::try_new(req).map_err(|_| failure("sparse LU workspace allocation failed".into()))
}

/// faer's sparse LU with a reused symbolic analysis.
#[derive(Default)]
pub(crate) struct FaerLu {
    pattern: Option<(Vec<usize>, Vec<usize>)>,
    symbolic: Option<lu::SymbolicLu<usize>>,
    numeric: lu::NumericLu<usize, f64>,
    factored: bool,
    factor_workspace: Option<MemBuffer>,
    solve_workspace: Option<MemBuffer>,
}

/// `nnz(H) + nnz(R)` of `A` (square, `n` columns) under the column
/// permutation `perm` of the symbolic analysis.
fn fill_bound(
    a: SymbolicSparseColMatRef<'_, usize>,
    perm: faer::perm::PermRef<'_, usize>,
) -> Result<u128, Diagnostic> {
    let n = a.ncols();
    // A^T in compressed column form.
    let mut counts = vec![0usize; n + 1];
    for &r in a.row_idx() {
        counts[r + 1] += 1;
    }
    for i in 0..n {
        counts[i + 1] += counts[i];
    }
    let t_ptr = counts.clone();
    let mut fill = counts;
    let mut t_idx = vec![0usize; a.row_idx().len()];
    for j in 0..n {
        for &r in a.row_idx_of_col_raw(j) {
            t_idx[fill[r]] = j;
            fill[r] += 1;
        }
    }
    let at = SymbolicSparseColMatRef::new_checked(n, n, &t_ptr, None, &t_idx);
    let req = StackReq::any_of(&[
        qr::col_etree_scratch::<usize>(n, n),
        qr::postorder_scratch::<usize>(n),
        qr::column_counts_aat_scratch::<usize>(n, n),
    ]);
    let mut buffer = workspace(req)?;
    let stack = MemStack::new(&mut buffer);
    let mut etree = vec![0usize; n];
    let mut post = vec![0usize; n];
    let mut col_counts = vec![0usize; n];
    let mut min_col = vec![0usize; n];
    let tree = qr::col_etree(a, Some(perm), &mut etree, stack);
    qr::postorder(&mut post, tree, stack);
    qr::column_counts_ata(
        &mut col_counts,
        &mut min_col,
        at,
        Some(perm),
        tree,
        &post,
        stack,
    );
    let parents = tree.into_inner();
    // Column counts of H: each row starts at its first column, and a column
    // passes all but one of its rows to its parent.
    let mut h = vec![0usize; n];
    for &first in &min_col {
        if (first as isize) >= 0 {
            h[first] += 1;
        }
    }
    for j in 0..n {
        let parent = parents[j];
        if parent >= 0 {
            h[parent as usize] += h[j].saturating_sub(1);
        }
    }
    Ok(h.iter().chain(&col_counts).map(|&c| c as u128).sum())
}

impl SparseLinearSolver for FaerLu {
    fn analyze(&mut self, pattern: &Csc) -> Result<u64, Diagnostic> {
        let n = pattern.dim;
        self.factored = false;
        self.symbolic = None;
        self.factor_workspace = None;
        self.solve_workspace = None;
        let a =
            SymbolicSparseColMatRef::new_checked(n, n, &pattern.col_ptr, None, &pattern.row_idx);
        let symbolic = lu::factorize_symbolic_lu(a, Default::default())
            .map_err(|e| failure(format!("sparse LU analysis failed: {e:?}")))?;
        let entries = fill_bound(a, symbolic.col_perm())?;
        let scratch = symbolic
            .factorize_numeric_lu_scratch::<f64>(Par::Seq, Default::default())
            .size_bytes() as u128;
        let bytes = 16 * entries + scratch + 32 * n as u128;
        self.pattern = Some((pattern.col_ptr.clone(), pattern.row_idx.clone()));
        self.symbolic = Some(symbolic);
        Ok(u64::try_from(bytes).unwrap_or(u64::MAX))
    }

    fn factorize(&mut self, matrix: &Csc) -> Result<(), Diagnostic> {
        self.factored = false;
        let (Some(symbolic), Some((col_ptr, row_idx))) = (&self.symbolic, &self.pattern) else {
            return Err(failure("sparse LU factorization before an analysis".into()));
        };
        if *col_ptr != matrix.col_ptr || *row_idx != matrix.row_idx {
            return Err(failure("sparse LU factorization of another pattern".into()));
        }
        let n = matrix.dim;
        if self.factor_workspace.is_none() {
            self.factor_workspace = Some(workspace(
                symbolic.factorize_numeric_lu_scratch::<f64>(Par::Seq, Default::default()),
            )?);
        }
        let buffer = self.factor_workspace.as_mut().expect("workspace");
        let structure = SymbolicSparseColMatRef::new_checked(n, n, col_ptr, None, row_idx);
        let a = SparseColMatRef::new(structure, &matrix.values);
        symbolic
            .factorize_numeric_lu(
                &mut self.numeric,
                a,
                Par::Seq,
                MemStack::new(buffer),
                Default::default(),
            )
            .map_err(|e| failure(format!("sparse LU factorization failed: {e:?}")))?;
        self.factored = true;
        Ok(())
    }

    fn solve_in_place(&mut self, rhs: &mut [f64]) -> Result<(), Diagnostic> {
        let Some(symbolic) = &self.symbolic else {
            return Err(failure("sparse LU solve without a factorization".into()));
        };
        if !self.factored || rhs.len() != symbolic.nrows() {
            return Err(failure("sparse LU solve without a factorization".into()));
        }
        if self.solve_workspace.is_none() {
            self.solve_workspace = Some(workspace(
                symbolic.solve_in_place_scratch::<f64>(1, Par::Seq),
            )?);
        }
        let buffer = self.solve_workspace.as_mut().expect("workspace");
        let n = rhs.len();
        lu::LuRef::new_unchecked(symbolic, &self.numeric).solve_in_place_with_conj(
            Conj::No,
            MatMut::from_column_major_slice_mut(rhs, n, 1),
            Par::Seq,
            MemStack::new(buffer),
        );
        Ok(())
    }
}
