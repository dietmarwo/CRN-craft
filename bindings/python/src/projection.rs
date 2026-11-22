//! Prepared projections: thin wrappers around `crn_craft::projection`.
//!
//! Vectors cross as little-endian bytes (`numpy.frombuffer` reads them
//! without a further copy); long work runs without the GIL on a helper
//! thread and Ctrl-C cancels it.
use crate::faer_lu::FaerLu;
use crate::{failure, run_interruptible};
use crn_craft::diagnostic::Diagnostic;
use crn_craft::projection::{
    BdfControls, ChannelShape, Generator, Implicit, Interrupt, Merged, MidpointControls, Operation,
    SignalValues, SolverCounts, Spec, TangentSolver, Topology, Work, bdf, kick, midpoint,
    stationary, uniform,
};
use crn_craft::rules::{Control, FlatView};
use pyo3::{
    prelude::*,
    types::{PyByteArray, PyBytes},
};
use std::sync::{Arc, Mutex};

/// Read-only rate evaluator shared with the native finite state projection.
#[pyclass(frozen)]
pub(crate) struct ProjectionView(Arc<FlatView>);

#[pymethods]
impl ProjectionView {
    #[new]
    fn new(py: Python<'_>, model: Vec<u8>) -> PyResult<Self> {
        py.detach(move || FlatView::of_any_model(&model))
            .map(|view| Self(Arc::new(view)))
            .map_err(|diagnostics| failure(&diagnostics))
    }

    fn layout(&self) -> String {
        let view = &self.0;
        let channels: Vec<_> = (0..view.channels())
            .map(|i| {
                let rule = view.rule(i);
                serde_json::json!({
                    "name": view.rule_name(rule),
                    "needs": view.reactants(i),
                    "changes": view.changes(i),
                    "whole": view.whole(i),
                    "factor": view.factor(i),
                    "rate": view.rate_constant(rule),
                })
            })
            .collect();
        serde_json::json!({
            "species": view.species(),
            "initial": view.initial(),
            "channels": channels,
            "has_stop": view.has_stop(),
        })
        .to_string()
    }

    fn evaluate(
        &self,
        py: Python<'_>,
        states: Vec<Vec<u64>>,
        channels: Vec<usize>,
    ) -> PyResult<(Vec<bool>, Vec<f64>)> {
        py.detach(move || self.0.projection_rates(&states, &channels))
            .map_err(|diagnostic| failure(&[diagnostic]))
    }
}

fn one(diagnostic: &Diagnostic) -> PyErr {
    failure(&[diagnostic])
}

fn floats(bytes: &[u8]) -> Vec<f64> {
    bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().expect("8 bytes")))
        .collect()
}

fn float_bytes<'py>(py: Python<'py>, values: &[f64]) -> PyResult<Bound<'py, PyByteArray>> {
    PyByteArray::new_with(py, values.len() * 8, |buffer| {
        for (chunk, value) in buffer.chunks_exact_mut(8).zip(values) {
            chunk.copy_from_slice(&value.to_le_bytes());
        }
        Ok(())
    })
}

fn index_bytes<'py>(py: Python<'py>, values: &[usize]) -> PyResult<Bound<'py, PyBytes>> {
    PyBytes::new_with(py, values.len() * 8, |buffer| {
        for (chunk, &value) in buffer.chunks_exact_mut(8).zip(values) {
            chunk.copy_from_slice(&(value as i64).to_le_bytes());
        }
        Ok(())
    })
}

fn work(count: i64, limit: i64) -> Work {
    Work::new(i128::from(count), i128::from(limit))
}

fn count(work: &Work) -> i64 {
    i64::try_from(work.count).unwrap_or(i64::MAX)
}

/// A count operation from Python: (species, kind, integer, probability),
/// kind 0 set, 1 change, 2 partition.
type PyOperation = (usize, u8, i64, f64);

fn operation(&(species, kind, integer, probability): &PyOperation) -> PyResult<(usize, Operation)> {
    let op = match kind {
        0 => Operation::Set(
            u64::try_from(integer)
                .map_err(|_| pyo3::exceptions::PyValueError::new_err("invalid set count"))?,
        ),
        1 => Operation::Change(integer),
        2 => Operation::Partition(probability),
        _ => return Err(pyo3::exceptions::PyValueError::new_err("unknown operation")),
    };
    Ok((species, op))
}

/// The reachable states and edges of a prepared projection.
#[pyclass(frozen)]
pub(crate) struct ProjectionTopology(Arc<Topology>);

#[pymethods]
impl ProjectionTopology {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (view, bounds, initial, interventions, max_states, max_edges, retained=None))]
    fn new(
        py: Python<'_>,
        view: &ProjectionView,
        bounds: Vec<u64>,
        initial: Vec<Vec<i64>>,
        interventions: Vec<Vec<PyOperation>>,
        max_states: usize,
        max_edges: usize,
        retained: Option<Vec<Vec<u64>>>,
    ) -> PyResult<Self> {
        let interventions = interventions
            .iter()
            .map(|ops| ops.iter().map(operation).collect::<PyResult<Vec<_>>>())
            .collect::<PyResult<Vec<_>>>()?;
        let view = Arc::clone(&view.0);
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let spec = Spec {
                species: view.species().to_vec(),
                bounds,
                initial,
                channels: ChannelShape::of_view(&view),
                interventions,
                max_states,
                max_edges,
            };
            let poll = || c.cancel_flagged();
            let retained = retained.map(|rows| rows.into_iter().collect());
            Topology::enumerate_retained(&spec, retained.as_ref(), &Interrupt::new(&poll))
        })?;
        result.map(|t| Self(Arc::new(t))).map_err(|d| one(&d))
    }

    /// The states, row-major little-endian int64.
    fn states<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let states = self.0.states();
        PyBytes::new_with(py, states.len() * 8, |buffer| {
            for (chunk, &value) in buffer.chunks_exact_mut(8).zip(states) {
                chunk.copy_from_slice(&value.to_le_bytes());
            }
            Ok(())
        })
    }

    /// The index of a state, `None` outside the projection.
    fn index(&self, counts: Vec<u64>) -> Option<usize> {
        self.0.index_of(&counts)
    }

    /// Builds the generator under a parameter binding (`view`).
    fn generator(
        &self,
        py: Python<'_>,
        view: &ProjectionView,
        signals: Vec<Vec<usize>>,
    ) -> PyResult<ProjectionGenerator> {
        let (topology, view) = (Arc::clone(&self.0), Arc::clone(&view.0));
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            Generator::build(&topology, &view, &signals, &Interrupt::new(&poll))
        })?;
        result
            .map(|g| ProjectionGenerator {
                generator: Arc::new(g),
                implicit: Arc::new(Mutex::new(None)),
            })
            .map_err(|d| one(&d))
    }

    /// Applies one intervention to a law (states and sink).
    fn kick<'py>(
        &self,
        py: Python<'py>,
        p: &[u8],
        changes: Vec<PyOperation>,
        work_count: i64,
        work_limit: i64,
    ) -> PyResult<(Bound<'py, PyByteArray>, i64)> {
        let changes = changes
            .iter()
            .map(operation)
            .collect::<PyResult<Vec<_>>>()?;
        let (topology, p) = (Arc::clone(&self.0), floats(p));
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let mut w = work(work_count, work_limit);
            kick(&topology, &p, &changes, &mut w, &Interrupt::new(&poll)).map(|p| (p, count(&w)))
        })?;
        let (p, used) = result.map_err(|d| one(&d))?;
        Ok((float_bytes(py, &p)?, used))
    }
}

/// A native signal source: a table, or a Python callable called with the
/// GIL reacquired.
enum Source {
    Table {
        step: bool,
        times: Vec<f64>,
        values: Vec<f64>,
    },
    Callable(Py<PyAny>),
}

/// `numpy.interp` (held endpoint values) or a right-continuous step table.
fn table_value(step: bool, t: &[f64], v: &[f64], at: f64) -> f64 {
    let last = t.len() - 1;
    if at.is_nan() {
        return f64::NAN;
    }
    // Index of the last knot at or before `at`, if any.
    let after = t.partition_point(|&x| x <= at);
    if step {
        return v[after.saturating_sub(1)];
    }
    if after == 0 {
        return v[0];
    }
    let j = after - 1;
    if j == last || t[j] == at {
        return v[j];
    }
    let slope = (v[j + 1] - v[j]) / (t[j + 1] - t[j]);
    let value = slope * (at - t[j]) + v[j];
    if value.is_nan() {
        let other = slope * (at - t[j + 1]) + v[j + 1];
        if other.is_nan() && v[j] == v[j + 1] {
            return v[j];
        }
        return other;
    }
    value
}

struct Sources {
    sources: Vec<Source>,
    error: Mutex<Option<PyErr>>,
}

impl SignalValues for Sources {
    fn evaluate(&mut self, time: f64, out: &mut [f64]) -> Result<(), Diagnostic> {
        let callback_failed = || {
            Diagnostic::error(
                "CRN-E070",
                "/signals",
                "a signal callable raised",
                "See the Python exception.",
            )
        };
        let mut pending: Option<PyErr> = None;
        for (value, source) in out.iter_mut().zip(&self.sources) {
            *value = match source {
                Source::Table {
                    step,
                    times,
                    values,
                } => table_value(*step, times, values, time),
                Source::Callable(function) => {
                    let result = Python::attach(|py| -> PyResult<f64> {
                        let answer = function.call1(py, (time,))?;
                        py.import("builtins")?
                            .getattr("float")?
                            .call1((answer,))?
                            .extract::<f64>()
                    });
                    match result {
                        Ok(v) => v,
                        Err(error) => {
                            pending = Some(error);
                            break;
                        }
                    }
                }
            };
        }
        if let Some(error) = pending {
            *self.error.lock().expect("signal error lock") = Some(error);
            return Err(callback_failed());
        }
        Ok(())
    }
}

/// A CSC matrix for Python: dimension, column pointers, rows, values.
type CscArrays<'py> = (
    usize,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyByteArray>,
);

/// The parsed signal sources of a call: tables or Python callables.
fn parse_sources(sources: Vec<Bound<'_, PyAny>>, expected: usize) -> PyResult<Sources> {
    let mut native = Vec::with_capacity(sources.len());
    for source in sources {
        if let Ok((step, times, values)) = source.extract::<(bool, Vec<f64>, Vec<f64>)>() {
            if times.is_empty() || times.len() != values.len() {
                return Err(pyo3::exceptions::PyValueError::new_err("invalid table"));
            }
            native.push(Source::Table {
                step,
                times,
                values,
            });
        } else {
            native.push(Source::Callable(source.unbind()));
        }
    }
    if native.len() != expected {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "one source per signal",
        ));
    }
    Ok(Sources {
        sources: native,
        error: Mutex::new(None),
    })
}

/// Linear-algebra counts for Python: symbolic analyses, numeric
/// factorizations, solves, products, accepted and rejected steps.
type PyCounts = (u64, u64, u64, u64, u64, u64);

fn py_counts(c: SolverCounts) -> PyCounts {
    (
        c.symbolic, c.numeric, c.solves, c.products, c.steps, c.rejected,
    )
}

/// The generator of a prepared projection under one binding, with the
/// analyzed iteration matrix of the implicit integrator once it is needed.
#[pyclass(frozen)]
pub(crate) struct ProjectionGenerator {
    generator: Arc<Generator>,
    implicit: Arc<Mutex<Option<Implicit<FaerLu>>>>,
}

#[pymethods]
impl ProjectionGenerator {
    /// The stop flag of each state, one byte per state.
    fn stopped<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        let flags: Vec<u8> = self
            .generator
            .stopped()
            .iter()
            .map(|&s| u8::from(s))
            .collect();
        PyBytes::new(py, &flags)
    }

    /// The number of signal matrices.
    fn signals(&self) -> usize {
        self.generator.signals().len()
    }

    /// One matrix as CSC arrays (int64 pointers, int64 rows, float64
    /// values): `None` for the constant part, `k` for signal `k`.
    #[pyo3(signature = (signal=None))]
    fn csc<'py>(&self, py: Python<'py>, signal: Option<usize>) -> PyResult<CscArrays<'py>> {
        let matrix = match signal {
            None => self.generator.constant(),
            Some(k) => self
                .generator
                .signals()
                .get(k)
                .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("no such signal"))?,
        };
        Ok((
            matrix.dim,
            index_bytes(py, &matrix.col_ptr)?,
            index_bytes(py, &matrix.row_idx)?,
            float_bytes(py, &matrix.values)?,
        ))
    }

    /// Uniformization of the constant part over `duration`: the law, the
    /// omitted Poisson weight and the work count.
    fn uniform<'py>(
        &self,
        py: Python<'py>,
        p: &[u8],
        duration: f64,
        tolerance: f64,
        work_count: i64,
        work_limit: i64,
    ) -> PyResult<(Bound<'py, PyByteArray>, f64, i64)> {
        let (generator, mut p) = (Arc::clone(&self.generator), floats(p));
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let mut w = work(work_count, work_limit);
            uniform(
                &generator,
                &mut p,
                duration,
                tolerance,
                &mut w,
                &Interrupt::new(&poll),
            )
            .map(|tail| (p, tail, count(&w)))
        })?;
        let (p, tail, used) = result.map_err(|d| one(&d))?;
        Ok((float_bytes(py, &p)?, tail, used))
    }

    /// The exponential midpoint rule from `start` to `end`. `sources` has
    /// one entry per signal: `(step, times, values)` for a table or a
    /// callable.
    #[allow(clippy::too_many_arguments)]
    fn midpoint<'py>(
        &self,
        py: Python<'py>,
        p: &[u8],
        start: f64,
        end: f64,
        tolerance: f64,
        max_step: f64,
        sources: Vec<Bound<'py, PyAny>>,
        work_count: i64,
        work_limit: i64,
    ) -> PyResult<(Bound<'py, PyByteArray>, i64)> {
        let parsed = parse_sources(sources, self.generator.signals().len())?;
        let (generator, mut p) = (Arc::clone(&self.generator), floats(p));
        let shared = Arc::new(Mutex::new(Some(parsed)));
        let inside = Arc::clone(&shared);
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let mut w = work(work_count, work_limit);
            let merged = Merged::new(&generator);
            let mut guard = inside.lock().expect("sources lock");
            let sources = guard.as_mut().expect("sources");
            let controls = MidpointControls {
                start,
                end,
                tolerance,
                max_step,
            };
            midpoint(
                &generator,
                &merged,
                &mut p,
                controls,
                sources,
                &mut w,
                &Interrupt::new(&poll),
            )
            .map(|()| (p, count(&w)))
        })?;
        let sources = shared.lock().expect("sources lock").take();
        if let Some(error) = sources.and_then(|s| s.error.into_inner().expect("signal error lock"))
        {
            return Err(error);
        }
        let (p, used) = result.map_err(|d| one(&d))?;
        Ok((float_bytes(py, &p)?, used))
    }

    /// Constant-segment forward sensitivity solve with a shared LU block.
    #[allow(clippy::too_many_arguments)]
    fn tangents<'py>(
        &self,
        py: Python<'py>,
        state: &[u8],
        pairs: Vec<(Py<Self>, Py<Self>, f64)>,
        start: f64,
        end: f64,
        reports: Vec<f64>,
        tolerance: f64,
        work_limit: i64,
        max_factor_bytes: u64,
    ) -> PyResult<(Bound<'py, PyByteArray>, Bound<'py, PyByteArray>, i64)> {
        let mut state = floats(state);
        let base = Arc::clone(&self.generator);
        let bindings: Vec<_> = pairs
            .iter()
            .map(|(a, b, width)| {
                (
                    Arc::clone(&a.borrow(py).generator),
                    Arc::clone(&b.borrow(py).generator),
                    *width,
                )
            })
            .collect();
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let interrupt = Interrupt::new(&poll);
            let n = base.dim();
            let views: Vec<_> = bindings
                .iter()
                .map(|(a, b, w)| (a.as_ref(), b.as_ref(), *w))
                .collect();
            if state.len() != n * (1 + views.len())
                || state.iter().any(|x| !x.is_finite())
                || !start.is_finite()
                || !end.is_finite()
                || start >= end
                || !(0.0 < tolerance && tolerance < 1.0)
                || reports.iter().any(|&t| !(start < t && t < end))
                || reports.windows(2).any(|w| w[0] >= w[1])
            {
                return Err(Diagnostic::error(
                    "CRN-E071",
                    "/projection/sensitivities",
                    "invalid tangent initial state, times or tolerance",
                    "",
                ));
            }
            let generator = base.tangents(&views)?;
            let solver = TangentSolver::new(n, FaerLu::default());
            let mut implicit = Implicit::new(&generator, solver)?;
            let mut sources = Sources {
                sources: Vec::new(),
                error: Mutex::new(None),
            };
            let controls = BdfControls {
                start,
                end,
                tolerance,
                max_step: f64::INFINITY,
                max_factor_bytes,
            };
            let mut w = work(0, work_limit);
            let result = bdf(
                &mut implicit,
                &mut state,
                controls,
                &reports,
                &mut sources,
                &mut w,
                &interrupt,
            )?;
            Ok((result.concat(), state, count(&w)))
        })?;
        let (reported, final_state, used) = result.map_err(|d| one(&d))?;
        Ok((
            float_bytes(py, &reported)?,
            float_bytes(py, &final_state)?,
            used,
        ))
    }

    /// The rate of each state into the sink, little-endian float64.
    fn escapes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyByteArray>> {
        float_bytes(py, &self.generator.escapes())
    }

    /// The implicit NDF integrator from `start` to `end` with reports at the
    /// increasing interior times `reports`: the reported laws (one after
    /// another), the law at `end`, the work count and the linear-algebra
    /// counts of this call. The analysis of the iteration matrix is kept
    /// with the generator for later calls.
    #[allow(clippy::too_many_arguments)]
    fn bdf<'py>(
        &self,
        py: Python<'py>,
        p: &[u8],
        start: f64,
        end: f64,
        tolerance: f64,
        max_step: f64,
        reports: Vec<f64>,
        sources: Vec<Bound<'py, PyAny>>,
        work_count: i64,
        work_limit: i64,
        max_factor_bytes: u64,
    ) -> PyResult<(
        Bound<'py, PyByteArray>,
        Bound<'py, PyByteArray>,
        i64,
        PyCounts,
    )> {
        let parsed = parse_sources(sources, self.generator.signals().len())?;
        let (generator, mut p) = (Arc::clone(&self.generator), floats(p));
        let implicit = Arc::clone(&self.implicit);
        let shared = Arc::new(Mutex::new(Some(parsed)));
        let inside = Arc::clone(&shared);
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let interrupt = Interrupt::new(&poll);
            let mut w = work(work_count, work_limit);
            let mut slot = implicit.lock().expect("implicit lock");
            let before = match slot.as_ref() {
                Some(implicit) => implicit.counts(),
                None => {
                    interrupt.check()?;
                    *slot = Some(Implicit::new(&generator, FaerLu::default())?);
                    SolverCounts::default()
                }
            };
            let implicit = slot.as_mut().expect("implicit");
            let mut guard = inside.lock().expect("sources lock");
            let sources = guard.as_mut().expect("sources");
            let controls = BdfControls {
                start,
                end,
                tolerance,
                max_step,
                max_factor_bytes,
            };
            let out = bdf(
                implicit, &mut p, controls, &reports, sources, &mut w, &interrupt,
            )?;
            let after = implicit.counts();
            let used = SolverCounts {
                symbolic: after.symbolic - before.symbolic,
                numeric: after.numeric - before.numeric,
                solves: after.solves - before.solves,
                products: after.products - before.products,
                steps: after.steps - before.steps,
                rejected: after.rejected - before.rejected,
            };
            Ok((out.concat(), p, count(&w), used))
        })?;
        let sources = shared.lock().expect("sources lock").take();
        if let Some(error) = sources.and_then(|s| s.error.into_inner().expect("signal error lock"))
        {
            return Err(error);
        }
        let (reported, p, used, counts) = result.map_err(|d: Diagnostic| one(&d))?;
        Ok((
            float_bytes(py, &reported)?,
            float_bytes(py, &p)?,
            used,
            py_counts(counts),
        ))
    }

    /// The stationary law of the constant part (see the Rust
    /// `stationary`): probabilities, residual, outward flux and the
    /// linear-algebra counts.
    #[pyo3(signature = (redirect, selected, max_factor_bytes))]
    fn stationary<'py>(
        &self,
        py: Python<'py>,
        redirect: Option<usize>,
        selected: Option<usize>,
        max_factor_bytes: u64,
    ) -> PyResult<(Bound<'py, PyByteArray>, f64, f64, PyCounts)> {
        let generator = Arc::clone(&self.generator);
        let result = run_interruptible(py, Control::default(), Control::cancel, move |c| {
            let poll = || c.cancel_flagged();
            let mut solver = FaerLu::default();
            stationary(
                &generator,
                redirect,
                selected,
                &mut solver,
                max_factor_bytes,
                &Interrupt::new(&poll),
            )
        })?;
        let law = result.map_err(|d| one(&d))?;
        Ok((
            float_bytes(py, &law.probabilities)?,
            law.residual,
            law.flux,
            py_counts(law.counts),
        ))
    }
}

/// The binomial probability the partition kernel uses (for validation).
#[pyfunction]
#[pyo3(name = "_binomial_pmf")]
fn binomial_pmf(x: f64, n: f64, p: f64) -> f64 {
    crn_craft::projection::special::binomial_pmf(x, n, p)
}

/// The Poisson upper tail `P(K >= first)` uniformization uses (for validation).
#[pyfunction]
#[pyo3(name = "_poisson_tail")]
fn poisson_tail(first: u64, mu: f64) -> f64 {
    crn_craft::projection::special::poisson_tail(first, mu)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(binomial_pmf, module)?)?;
    module.add_function(wrap_pyfunction!(poisson_tail, module)?)?;
    module.add_class::<ProjectionView>()?;
    module.add_class::<ProjectionTopology>()?;
    module.add_class::<ProjectionGenerator>()?;
    Ok(())
}
