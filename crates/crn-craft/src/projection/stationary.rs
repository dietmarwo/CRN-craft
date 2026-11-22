//! Stationary laws of a projection generator: closed communicating
//! classes and one anchored sparse solve per class.
use super::generator::{Csc, Generator};
use super::{BUDGET, Interrupt, SolverCounts, SparseLinearSolver, input, memory};
use crate::diagnostic::Diagnostic;

/// A stationary law on the states of a projection (the sink excluded).
#[derive(Debug, Clone)]
pub struct StationaryLaw {
    /// Probability of each state, zero outside the solved class.
    pub probabilities: Vec<f64>,
    /// `|Q p|_1` of the solved generator.
    pub residual: f64,
    /// `sum_j escape_j p_j`: the rate at which the law leaves the bounds.
    pub flux: f64,
    /// The linear-algebra work.
    pub counts: SolverCounts,
}

fn failed(message: &str) -> Diagnostic {
    Diagnostic::error(
        BUDGET,
        "/projection",
        message,
        "Check that the generator has one closed class with a unique stationary law.",
    )
}

/// The generator restricted to the states (the sink row and column
/// dropped), with the escapes of each column added at `redirect`.
fn restricted(generator: &Generator, redirect: Option<usize>) -> Csc {
    let q = generator.constant();
    let n = q.dim - 1;
    let escapes = generator.escapes();
    let mut col_ptr = Vec::with_capacity(n + 1);
    let mut row_idx = Vec::new();
    let mut values = Vec::new();
    col_ptr.push(0);
    for (j, &escape) in escapes.iter().enumerate() {
        let mut column: Vec<(usize, f64)> = (q.col_ptr[j]..q.col_ptr[j + 1])
            .filter(|&k| q.row_idx[k] < n)
            .map(|k| (q.row_idx[k], q.values[k]))
            .collect();
        if let Some(anchor) = redirect {
            if escape != 0.0 {
                match column.binary_search_by_key(&anchor, |&(r, _)| r) {
                    Ok(at) => column[at].1 += escape,
                    Err(at) => column.insert(at, (anchor, escape)),
                }
            }
        }
        for (r, v) in column {
            row_idx.push(r);
            values.push(v);
        }
        col_ptr.push(row_idx.len());
    }
    Csc {
        dim: n,
        col_ptr,
        row_idx,
        values,
    }
}

/// The strongly connected components of the graph with an edge from `j`
/// to `i` for every stored off-diagonal entry `(i, j)` of `q`, by an
/// iterative Tarjan search. Returns the component of each node.
fn components(q: &Csc, interrupt: &Interrupt<'_>) -> Result<(Vec<usize>, usize), Diagnostic> {
    let n = q.dim;
    const NONE: usize = usize::MAX;
    let mut index = vec![NONE; n];
    let mut low = vec![0; n];
    let mut on_stack = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut label = vec![NONE; n];
    let mut labels = 0;
    let mut counter = 0;
    // Call stack: (node, next entry of its column).
    let mut calls: Vec<(usize, usize)> = Vec::new();
    for root in 0..n {
        if index[root] != NONE {
            continue;
        }
        if root % 1024 == 0 {
            interrupt.check()?;
        }
        index[root] = counter;
        low[root] = counter;
        counter += 1;
        stack.push(root);
        on_stack[root] = true;
        calls.push((root, q.col_ptr[root]));
        while let Some(&mut (v, ref mut next)) = calls.last_mut() {
            let end = q.col_ptr[v + 1];
            let mut descended = false;
            while *next < end {
                let w = q.row_idx[*next];
                *next += 1;
                if w == v {
                    continue;
                }
                if index[w] == NONE {
                    index[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push(w);
                    on_stack[w] = true;
                    calls.push((w, q.col_ptr[w]));
                    descended = true;
                    break;
                } else if on_stack[w] {
                    low[v] = low[v].min(index[w]);
                }
            }
            if descended {
                continue;
            }
            calls.pop();
            if let Some(&(parent, _)) = calls.last() {
                low[parent] = low[parent].min(low[v]);
            }
            if low[v] == index[v] {
                loop {
                    let w = stack.pop().expect("Tarjan stack holds the component");
                    on_stack[w] = false;
                    label[w] = labels;
                    if w == v {
                        break;
                    }
                }
                labels += 1;
            }
        }
    }
    Ok((label, labels))
}

/// The closed communicating classes of `q`, each as increasing state
/// indices, in the order of their smallest state. A class is closed when
/// no positive entry leads from it to another class.
fn closed_classes(q: &Csc, interrupt: &Interrupt<'_>) -> Result<Vec<Vec<usize>>, Diagnostic> {
    let (label, count) = components(q, interrupt)?;
    let mut open = vec![false; count];
    for j in 0..q.dim {
        for k in q.col_ptr[j]..q.col_ptr[j + 1] {
            if q.values[k] > 0.0 && label[q.row_idx[k]] != label[j] {
                open[label[j]] = true;
            }
        }
    }
    let mut position = vec![usize::MAX; count];
    let mut classes: Vec<Vec<usize>> = Vec::new();
    for (state, &l) in label.iter().enumerate() {
        if open[l] {
            continue;
        }
        if position[l] == usize::MAX {
            position[l] = classes.len();
            classes.push(Vec::new());
        }
        classes[position[l]].push(state);
    }
    Ok(classes)
}

/// The submatrix of `q` on the increasing indices `ids`.
fn submatrix(q: &Csc, ids: &[usize]) -> Csc {
    let mut local = vec![usize::MAX; q.dim];
    for (i, &g) in ids.iter().enumerate() {
        local[g] = i;
    }
    let mut col_ptr = Vec::with_capacity(ids.len() + 1);
    let mut row_idx = Vec::new();
    let mut values = Vec::new();
    col_ptr.push(0);
    for &j in ids {
        for k in q.col_ptr[j]..q.col_ptr[j + 1] {
            let r = local[q.row_idx[k]];
            if r != usize::MAX {
                row_idx.push(r);
                values.push(q.values[k]);
            }
        }
        col_ptr.push(row_idx.len());
    }
    Csc {
        dim: ids.len(),
        col_ptr,
        row_idx,
        values,
    }
}

/// The stationary law of `q` on its single closed class (or, with
/// `anchor`, on the class holding it): the probability of the anchor (the
/// largest state of the class without one) is fixed to 1, the balance
/// equations of the other states of the class are solved, and the result
/// is clipped at zero and normalized. Returns the law and `|q p|_1`.
fn anchored<S: SparseLinearSolver>(
    q: &Csc,
    anchor: Option<usize>,
    solver: &mut S,
    max_factor_bytes: u64,
    counts: &mut SolverCounts,
    interrupt: &Interrupt<'_>,
) -> Result<(Vec<f64>, f64), Diagnostic> {
    let classes = closed_classes(q, interrupt)?;
    let (ids, anchor) = match anchor {
        None => {
            if classes.len() != 1 {
                return Err(input(format!(
                    "projection has {} closed communicating classes",
                    classes.len()
                )));
            }
            let ids = &classes[0];
            (ids, ids[ids.len() - 1])
        }
        Some(anchor) => {
            let Some(ids) = classes
                .iter()
                .find(|ids| ids.binary_search(&anchor).is_ok())
            else {
                return Err(input("the designated state is outside every closed class"));
            };
            if classes.len() != 1 {
                return Err(input(
                    "the redirected projection has multiple closed classes",
                ));
            }
            (ids, anchor)
        }
    };
    let remaining: Vec<usize> = ids.iter().copied().filter(|&s| s != anchor).collect();
    let mut probability = vec![0.0; q.dim];
    probability[anchor] = 1.0;
    if !remaining.is_empty() {
        let matrix = submatrix(q, &remaining);
        let mut rhs = vec![0.0; remaining.len()];
        let mut local = vec![usize::MAX; q.dim];
        for (i, &g) in remaining.iter().enumerate() {
            local[g] = i;
        }
        for k in q.col_ptr[anchor]..q.col_ptr[anchor + 1] {
            let r = local[q.row_idx[k]];
            if r != usize::MAX {
                rhs[r] = -q.values[k];
            }
        }
        interrupt.check()?;
        // The factorization bound and the restricted matrices and vectors.
        let bytes = solver
            .analyze(&matrix)?
            .saturating_add(32 * (matrix.nnz() + q.nnz()) as u64)
            .saturating_add(48 * q.dim as u64);
        counts.symbolic += 1;
        memory(bytes, max_factor_bytes)?;
        interrupt.check()?;
        if solver.factorize(&matrix).is_err() {
            return Err(failed(
                "stationary sparse solve returned invalid probabilities",
            ));
        }
        counts.numeric += 1;
        solver.solve_in_place(&mut rhs)?;
        counts.solves += 1;
        for (&g, &x) in remaining.iter().zip(&rhs) {
            probability[g] = x;
        }
    }
    if probability.iter().any(|x| !x.is_finite())
        || probability.iter().copied().fold(f64::INFINITY, f64::min) < -1e-10
    {
        return Err(failed(
            "stationary sparse solve returned invalid probabilities",
        ));
    }
    for x in &mut probability {
        *x = x.max(0.0);
    }
    let total: f64 = probability.iter().sum();
    for x in &mut probability {
        *x /= total;
    }
    let mut product = vec![0.0; q.dim];
    let mut scale: f64 = f64::NEG_INFINITY;
    for (j, &x) in probability.iter().enumerate() {
        for k in q.col_ptr[j]..q.col_ptr[j + 1] {
            product[q.row_idx[k]] += q.values[k] * x;
        }
        scale = scale.max(-q.diagonal(j));
    }
    counts.products += 1;
    let residual: f64 = product.iter().map(|x| x.abs()).sum();
    if residual > 1e-9 * scale.max(1.0) {
        return Err(failed("stationary generator residual exceeded tolerance"));
    }
    Ok((probability, residual))
}

/// The stationary law of the constant part of `generator` on its states.
///
/// - Without `redirect`, the generator restricted to the states must have
///   one closed communicating class, or `selected` names a state of the
///   class to solve.
/// - With `redirect` (a state index), the rate of each state into the sink
///   returns to that state, and the redirected generator must have one
///   closed class, which holds `redirect`.
///
/// The caller checks that the indices are states and refuses escapes
/// without a redirect. The closed classes come from an iterative
/// strongly-connected-components pass over the positive off-diagonal
/// entries; a class is closed when no entry leaves it. The solve anchors
/// the largest state of the class (or `redirect`) at probability 1 and
/// solves the balance equations of the others with one sparse
/// factorization, refused above `max_factor_bytes`; the law is clipped at
/// zero, normalized and accepted when `|Q p|_1 <= 1e-9 max(1, max exit
/// rate)`.
///
/// # Errors
///
/// `CRN-E071` for a class count other than one or a state outside every
/// closed class; `CRN-E070` for a factorization above the budget, a law
/// with a non-finite or negative entry below `-1e-10`, a residual above
/// the tolerance, or an interruption; any error of the solver's analysis.
pub fn stationary<S: SparseLinearSolver>(
    generator: &Generator,
    redirect: Option<usize>,
    selected: Option<usize>,
    solver: &mut S,
    max_factor_bytes: u64,
    interrupt: &Interrupt<'_>,
) -> Result<StationaryLaw, Diagnostic> {
    let q = restricted(generator, redirect);
    let mut counts = SolverCounts::default();
    let probabilities = if let Some(target) = selected {
        let classes = closed_classes(&q, interrupt)?;
        let Some(ids) = classes
            .iter()
            .find(|ids| ids.binary_search(&target).is_ok())
        else {
            return Err(input("selected state is outside every closed class"));
        };
        let sub = submatrix(&q, ids);
        let (p, residual) = anchored(&sub, None, solver, max_factor_bytes, &mut counts, interrupt)?;
        let mut full = vec![0.0; q.dim];
        for (&g, &x) in ids.iter().zip(&p) {
            full[g] = x;
        }
        (full, residual)
    } else {
        anchored(
            &q,
            redirect,
            solver,
            max_factor_bytes,
            &mut counts,
            interrupt,
        )?
    };
    let (probabilities, residual) = probabilities;
    let flux = generator
        .escapes()
        .iter()
        .zip(&probabilities)
        .map(|(e, p)| e * p)
        .sum();
    Ok(StationaryLaw {
        probabilities,
        residual,
        flux,
        counts,
    })
}
