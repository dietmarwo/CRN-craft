//! Prepared finite state projections: one enumerated topology, many solves.
//!
//! This is the numerical engine of the Python API `crn_craft.prepare_fsp`
//! (`PreparedFsp`). It is independent of the JSON `fsp` method of
//! `crn-craft.run/1.0` ([`crate::fsp`]), whose enumeration, results and
//! budgets it does not share or change.
//!
//! - [`Topology::enumerate`] finds the states reachable inside inclusive
//!   bounds from an initial support, with scheduled count interventions
//!   expanding the reachable set, and the edges between them (one per state
//!   and channel that changes the state; leaving the bounds goes to one
//!   absorbing sink). The order of states, edges and multiplicities is the
//!   breadth-first discovery order, so a projection is identified by its
//!   state table.
//! - [`Generator::build`] assembles, for one parameter binding, the
//!   constant part of the generator, one matrix per time-dependent signal
//!   and the stop states, in compressed sparse column form ([`Csc`]).
//! - [`uniform`] propagates a law under the constant part by positive
//!   uniformization and returns the omitted Poisson weight.
//! - [`midpoint`] integrates time-dependent signals with the adaptive
//!   exponential midpoint rule, each matrix-exponential action by
//!   uniformization of the frozen generator.
//! - [`bdf`] integrates time-dependent signals with the variable-order
//!   numerical differentiation formulas (orders 1 to 5) on an
//!   [`Implicit`] iteration matrix, through a [`SparseLinearSolver`].
//! - [`stationary`] solves the stationary law of one closed communicating
//!   class with one anchored sparse solve.
//! - [`kick`] applies a count intervention to a law with exact binomial
//!   partition weights.
//!
//! Failures are diagnostics: `CRN-E070` for budgets, work, the
//! integrators' stalls, the factorization memory budget, a stationary
//! solve that fails and interruption; `CRN-E071` for an input or law the
//! projection cannot represent (including a class count other than one);
//! and the rate evaluator's codes (`CRN-R030`).
//! Everything runs on the calling thread; results do not depend on timing.
//!
//! ```
//! use crn_craft::projection::{ChannelShape, Generator, Interrupt, Spec, Topology, Work, uniform};
//! use crn_craft::rules::FlatView;
//! use serde_json::json;
//!
//! // Pure death from 3 copies at rate 1.
//! let model = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.model/1.0", "name": "death",
//!     "parameters": {"k": 1.0},
//!     "types": [{"name": "X"}],
//!     "initial": [{"type": "X", "count": 3}],
//!     "observables": {"X": {"type": "X"}},
//!     "rules": [{"name": "death", "reactants": [{"label": "x", "type": "X"}],
//!                "rate": "k", "products": []}]
//! }))
//! .unwrap();
//! let view = FlatView::from_bytes(&model).unwrap();
//! let spec = Spec {
//!     species: view.species().to_vec(),
//!     bounds: vec![3],
//!     initial: vec![vec![3]],
//!     channels: ChannelShape::of_view(&view),
//!     interventions: Vec::new(),
//!     max_states: 10,
//!     max_edges: 10,
//! };
//! let never = || false;
//! let interrupt = Interrupt::new(&never);
//! let topology = Topology::enumerate(&spec, &interrupt).unwrap();
//! assert_eq!(topology.states(), [3, 2, 1, 0]);
//! let generator = Generator::build(&topology, &view, &[], &interrupt).unwrap();
//! let mut p = vec![1.0, 0.0, 0.0, 0.0, 0.0]; // four states and the sink
//! let mut work = Work::new(0, 1_000_000);
//! let tail = uniform(&generator, &mut p, 1.0, 1e-12, &mut work, &interrupt).unwrap();
//! assert!(tail <= 1e-12);
//! let survive = (-1.0f64).exp();
//! assert!((p[0] - survive.powi(3)).abs() < 1e-12); // all three survive
//! ```
mod generator;
mod implicit;
mod propagate;
mod sensitivity;
pub mod special;
mod stationary;
pub use sensitivity::TangentSolver;
mod topology;

#[cfg(test)]
mod tests;

pub use generator::{Csc, Generator, MAX_RATE_TABLE, Merged};
pub use implicit::{BdfControls, Implicit, MAX_ORDER, bdf};
pub use propagate::{
    LAST_STEP_SLACK, MIDPOINT_PIECE_MEAN, MIDPOINT_TRUNCATION_SHARE, MidpointControls,
    SignalValues, kick, midpoint, uniform,
};
pub use stationary::{StationaryLaw, stationary};
pub use topology::{
    ChannelShape, MAX_COUNT, MAX_EDGES, MAX_STATES, STATE_TABLE_BYTES, Spec, Topology,
};

use crate::diagnostic::Diagnostic;

/// The diagnostic code of a projection budget, work, stall or interruption.
pub const BUDGET: &str = "CRN-E070";
/// The diagnostic code of a projection input or law it cannot represent.
pub const INPUT: &str = "CRN-E071";

fn budget(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        BUDGET,
        "/projection",
        message,
        "Lower the bounds, raise the budget within its cap, or loosen the tolerance.",
    )
}

fn input(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        INPUT,
        "/projection",
        message,
        "Revise the initial law, the interventions or the rates.",
    )
}

/// A count operation of a scheduled intervention on one species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Operation {
    /// Set the count.
    Set(u64),
    /// Add a (possibly negative) number of molecules. Increments beyond the
    /// range of `i64` behave like any increment that leaves the bounds or
    /// goes below zero.
    Change(i64),
    /// Keep each molecule independently with this probability in `[0, 1]`.
    Partition(f64),
}

/// The deterministic work counter of one solve: each uniformization series
/// term, each matrix-exponential action and each intervention outcome costs
/// one unit.
#[derive(Debug, Clone, Copy)]
pub struct Work {
    /// Units used.
    pub count: i128,
    /// Most units.
    pub limit: i128,
}

impl Work {
    /// A counter that has used `count` of `limit` units.
    #[must_use]
    pub fn new(count: i128, limit: i128) -> Self {
        Self { count, limit }
    }

    /// Charges `amount` units.
    ///
    /// # Errors
    ///
    /// `CRN-E070` at `/projection/work` when the count exceeds the limit.
    pub fn charge(&mut self, amount: u128) -> Result<(), Diagnostic> {
        self.count = self
            .count
            .saturating_add(i128::try_from(amount).unwrap_or(i128::MAX));
        if self.count > self.limit {
            return Err(Diagnostic::error(
                BUDGET,
                "/projection/work",
                format!(
                    "projection work limit of {} exceeded; raise max_work",
                    self.limit
                ),
                "Raise max_work or loosen the tolerance.",
            ));
        }
        Ok(())
    }
}

/// A cooperative cancellation poll.
pub struct Interrupt<'a> {
    cancelled: &'a dyn Fn() -> bool,
}

impl<'a> Interrupt<'a> {
    /// Polls `cancelled`.
    pub fn new(cancelled: &'a dyn Fn() -> bool) -> Self {
        Self { cancelled }
    }

    /// # Errors
    ///
    /// `CRN-E070` once the poll reports a cancellation.
    pub fn check(&self) -> Result<(), Diagnostic> {
        if (self.cancelled)() {
            return Err(budget("projection interrupted"));
        }
        Ok(())
    }
}

/// The linear-algebra boundary of an implicit integrator or a stationary
/// solve on a projection generator.
///
/// No implementation lives in this crate: a sparse direct solver is an
/// adapter outside it (the Python binding implements it on faer's sparse
/// LU), so the engine keeps its dependency set. The matrices it receives
/// are [`Csc`] with sorted rows. One [`analyze`](Self::analyze) fixes a
/// pattern; every later [`factorize`](Self::factorize) has that pattern
/// with new values (an implicit step factorizes `I - c Q(t)` on the fixed
/// pattern of [`Merged`] with the diagonal added; a stationary solve
/// analyzes once per closed class).
pub trait SparseLinearSolver {
    /// The symbolic analysis of the pattern of `pattern` (its values are
    /// ignored). Returns an upper bound of the bytes one numeric
    /// factorization on this pattern needs, so that the caller can refuse
    /// before factorizing.
    ///
    /// # Errors
    ///
    /// A diagnostic when the analysis fails.
    fn analyze(&mut self, pattern: &Csc) -> Result<u64, Diagnostic>;
    /// Factorizes `matrix`, whose pattern is the analyzed one.
    ///
    /// # Errors
    ///
    /// A diagnostic for a singular or unsupported matrix, or when no
    /// analysis matches.
    fn factorize(&mut self, matrix: &Csc) -> Result<(), Diagnostic>;
    /// Solves with the last factorization, in place.
    ///
    /// # Errors
    ///
    /// A diagnostic when no factorization is available.
    fn solve_in_place(&mut self, rhs: &mut [f64]) -> Result<(), Diagnostic>;
}

/// The linear-algebra work of an implicit integration or a stationary
/// solve.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SolverCounts {
    /// Symbolic analyses (one per pattern).
    pub symbolic: u64,
    /// Numeric factorizations.
    pub numeric: u64,
    /// Solves with a factorization.
    pub solves: u64,
    /// Products of the generator with a vector.
    pub products: u64,
    /// Accepted steps.
    pub steps: u64,
    /// Rejected step attempts (error test failed).
    pub rejected: u64,
}

/// Refuses a factorization whose bound `bytes` exceeds `budget`.
fn memory(bytes: u64, budget: u64) -> Result<(), Diagnostic> {
    if bytes > budget {
        return Err(Diagnostic::error(
            BUDGET,
            "/projection/memory",
            format!(
                "the sparse LU factorization and its workspace may need up to \
                 {bytes} bytes, above the memory budget max_factor_bytes={budget}"
            ),
            "Raise max_factor_bytes or lower the bounds.",
        ));
    }
    Ok(())
}
