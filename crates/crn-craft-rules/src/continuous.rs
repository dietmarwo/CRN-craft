//! Continuous variables (rules-model/0.3 and 0.4) and the `hybrid` method.
//!
//! A model may declare continuous variables ([`ContinuousSpec`]): each has an
//! initial value, a derivative written as an expression over parameters,
//! observables and continuous variables, and an optional lower bound. Rule
//! rates and propensities may read them, so the jump rates change between
//! events. Results carry each variable as an observable column after the
//! model's observables, in declaration order.
//!
//! Such a model is a piecewise-deterministic Markov process (Davis 1984,
//! DOI [10.1111/j.2517-6161.1984.tb01308.x](https://doi.org/10.1111/j.2517-6161.1984.tb01308.x)).
//! The `hybrid` method simulates it by integrated propensities: after each
//! event it draws a unit exponential threshold E and integrates the
//! variables together with the cumulative total propensity
//! Λ(t) = ∫ a₀(s) ds until Λ reaches E; the next event happens there and
//! selects a channel with probability a_j / a₀ at that time (Alfonsi,
//! Cancès, Turinici, Di Ventura and Huisinga 2005, DOI
//! [10.1051/proc:2005001](https://doi.org/10.1051/proc:2005001)). The
//! integrator is the embedded Runge–Kutta pair of Dormand and Prince, 5(4)
//! with local extrapolation (DOI
//! [10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3)),
//! with the standard step-size control of Hairer, Nørsett and Wanner
//! (*Solving Ordinary Differential Equations I*, 2nd ed., Springer 1993,
//! sections II.4 and II.6) and the fourth-order continuous extension given
//! there after Shampine (1986, DOI
//! [10.1090/S0025-5718-1986-0815836-3](https://doi.org/10.1090/S0025-5718-1986-0815836-3)).
//! Event times, report values and bound crossings are located on that dense
//! output. The method is approximate at the request's tolerances; with
//! constant rates between events it is exact up to rounding.
//!
//! A lower bound with `on_bound: "stop"` ends the replicate as `stopped` when
//! the variable reaches the bound while its derivative there is not
//! positive, as the legacy ATP block does at zero; the variable reads the
//! bound value on every later report. With `"clamp"` the variable is held at
//! the bound while its derivative there is negative and follows the ODE
//! again once it is positive.
//!
//! Rules-model/0.4 adds [`ContinuousSpec::independent`]. An independent
//! variable's derivative reads only parameters and other independent
//! variables, so its trajectory is the same in every replicate. `hybrid`
//! integrates the independent subsystem once per prepared run over
//! [0, end], with the same integrator, step control and bound handling at
//! the request's tolerances, and keeps the dense output of every accepted
//! step. Each replicate integrates only the other variables and Λ and reads
//! the independent ones from that dense output, so the results do not
//! depend on the number of threads. A stop bound of an independent variable
//! ends every replicate that reaches it at the same time. A failure of the
//! shared solve (a non-finite derivative, an underflowing step or the step
//! budget) fails every replicate that reaches its time, with its diagnostic.

use crate::expr::{Env, Expr, Scope};
use crate::model::{Value, invalid, parse_expr, pointer};
use crn_craft::diagnostic::Diagnostic;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Most continuous variables in a model (`CRN-R015`).
pub const MAX_VARIABLES: usize = 64;

/// What happens when a continuous variable reaches its lower bound.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum OnBound {
    // Variant docs use `concat!` to stay out of the schema's string enum.
    #[doc = concat!(
        "End the replicate as `stopped` when the variable reaches the bound ",
        "with a derivative that is not positive; the variable reads the ",
        "bound on every later report."
    )]
    Stop,
    #[doc = concat!(
        "Hold the variable at the bound while its derivative there is ",
        "negative; it follows the ODE again once the derivative is positive."
    )]
    Clamp,
}

/// A continuous variable of a rules-model/0.3 or 0.4 model.
///
/// Its name shares one namespace with parameters and observables
/// (`CRN-R010`). Rule rates and propensities may read it, and so may the
/// derivatives of all variables; the model's `stop`, observables, patterns
/// and products may not (`CRN-R012`), and neither may run-request trackers.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContinuousSpec {
    /// The variable's name and result column.
    pub name: String,
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!(
        "Value at time 0: a finite number or an expression over parameters ",
        "(`CRN-R014` when not finite or below `lower_bound`)."
    )]
    pub initial: Value,
    /// dx/dt between events: an expression over parameters, observables and
    /// continuous variables. A non-finite value fails the replicate
    /// (`CRN-R030`).
    pub derivative: String,
    #[doc = concat!(
        "Optional lower bound: a finite number or an expression over ",
        "parameters. It requires `on_bound`, and `on_bound` requires it ",
        "(`CRN-R013`). `null` is rejected."
    )]
    #[serde(
        default,
        deserialize_with = "crate::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Value")]
    pub lower_bound: Option<Value>,
    #[doc = concat!(
        "What reaching `lower_bound` does: `stop` or `clamp`. `null` is ",
        "rejected."
    )]
    #[serde(
        default,
        deserialize_with = "crate::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "OnBound")]
    pub on_bound: Option<OnBound>,
    /// rules-model/0.4 (`CRN-R015` before): `true` declares an independent
    /// ODE, whose `derivative` reads only parameters and other independent
    /// variables (`CRN-R012` otherwise). `hybrid` integrates the independent
    /// variables once per prepared run and shares the solution by every
    /// replicate. Default `false`; `null` is rejected.
    #[serde(
        default,
        deserialize_with = "crate::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "bool")]
    pub independent: Option<bool>,
}

/// A compiled continuous variable.
#[derive(Debug)]
pub(crate) struct Variable {
    pub name: String,
    /// Column in the observable values.
    pub index: usize,
    pub initial: f64,
    pub derivative: Expr,
    pub bound: Option<(f64, OnBound)>,
    /// Declared `independent: true` (rules-model/0.4).
    pub independent: bool,
}

fn parameter_value(
    value: &Value,
    params: &BTreeMap<String, usize>,
    values: &[f64],
    path: &str,
    out: &mut Vec<Diagnostic>,
) -> Option<f64> {
    let x = match value {
        Value::Number(x) => *x,
        Value::Expression(source) => {
            let scope = Scope {
                params,
                observables: None,
                reactants: &[],
                random: false,
            };
            parse_expr(source, &scope, path, out)?.eval(&Env {
                params: values,
                observables: &[],
                attributes: [&[], &[]],
            })
        }
    };
    if x.is_finite() {
        Some(if x == 0.0 { 0.0 } else { x })
    } else {
        out.push(invalid(
            "CRN-R014",
            path,
            format!("value is {x}"),
            "Use a finite number.",
        ));
        None
    }
}

impl Variable {
    /// Compiles the declared variables; their columns start at `first`.
    pub(crate) fn compile_all(
        specs: &[ContinuousSpec],
        version: crate::model::ModelVersion,
        params: &BTreeMap<String, usize>,
        names: &BTreeMap<String, usize>,
        values: &[f64],
        first: usize,
        out: &mut Vec<Diagnostic>,
    ) -> Vec<Self> {
        let scope = Scope {
            params,
            observables: Some(names),
            reactants: &[],
            random: false,
        };
        let mut compiled = Vec::new();
        for (k, spec) in specs.iter().enumerate() {
            let path = pointer("/continuous", k);
            let initial = parameter_value(
                &spec.initial,
                params,
                values,
                &pointer(&path, "initial"),
                out,
            );
            let derivative =
                parse_expr(&spec.derivative, &scope, &pointer(&path, "derivative"), out);
            let bound = match (&spec.lower_bound, spec.on_bound) {
                (None, None) => Some(None),
                (Some(value), Some(action)) => {
                    parameter_value(value, params, values, &pointer(&path, "lower_bound"), out)
                        .map(|b| Some((b, action)))
                }
                (Some(_), None) | (None, Some(_)) => {
                    out.push(invalid(
                        "CRN-R013",
                        &path,
                        "lower_bound and on_bound go together",
                        "Give both lower_bound and on_bound, or neither.",
                    ));
                    None
                }
            };
            if let (Some(x), Some(Some((b, _)))) = (initial, bound) {
                if x < b {
                    out.push(invalid(
                        "CRN-R014",
                        &pointer(&path, "initial"),
                        format!("initial value {x} is below the lower bound {b}"),
                        "Start at or above the lower bound.",
                    ));
                }
            }
            let independent = spec.independent == Some(true);
            if spec.independent.is_some() && version < crate::model::ModelVersion::V04 {
                out.push(invalid(
                    "CRN-R015",
                    &pointer(&path, "independent"),
                    "independent requires rules-model/0.4",
                    "Select the 0.4 model contract.",
                ));
            }
            if let (true, Some(derivative)) = (independent, &derivative) {
                // Only parameters and independent variables, so that one
                // solution serves every replicate.
                let mut flags = vec![false; first + specs.len()];
                derivative.reads(&mut flags);
                let name = |i: usize| {
                    names
                        .iter()
                        .find(|&(_, &j)| j == i)
                        .map_or("?", |(n, _)| n.as_str())
                };
                if let Some(i) = flags.iter().enumerate().position(|(i, &read)| {
                    read && !(i >= first && specs[i - first].independent == Some(true))
                }) {
                    let what = if i < first {
                        "observable"
                    } else {
                        "dependent continuous variable"
                    };
                    out.push(invalid(
                        "CRN-R012",
                        &pointer(&path, "derivative"),
                        format!(
                            "the derivative of independent variable `{}` reads {what} `{}`",
                            spec.name,
                            name(i)
                        ),
                        "An independent derivative reads only parameters and independent continuous variables; drop `independent` or the read.",
                    ));
                }
            }
            if let (Some(initial), Some(derivative), Some(bound)) = (initial, derivative, bound) {
                compiled.push(Self {
                    name: spec.name.clone(),
                    index: first + k,
                    initial,
                    derivative,
                    bound,
                    independent,
                });
            }
        }
        compiled
    }
}

// ---------------------------------------------------------------------------
// Dormand–Prince 5(4)

/// Nodes c_i of the Dormand–Prince tableau.
const C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];

/// Coefficients a_ij below the diagonal; row 7 equals the fifth-order
/// weights (first same as last).
const A: [[f64; 6]; 7] = [
    [0.0; 6],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];

/// Fifth-order solution minus the embedded fourth-order one, per stage.
const E: [f64; 7] = [
    71.0 / 57600.0,
    0.0,
    -71.0 / 16695.0,
    71.0 / 1920.0,
    -17253.0 / 339200.0,
    22.0 / 525.0,
    -1.0 / 40.0,
];

/// Coefficients of the fourth-order continuous extension (Hairer, Nørsett
/// and Wanner, section II.6, after Shampine 1986).
const D: [f64; 7] = [
    -12715105075.0 / 11282082432.0,
    0.0,
    87487479700.0 / 32700410799.0,
    -10690763975.0 / 1880347072.0,
    701980252875.0 / 199316789632.0,
    -1453857185.0 / 822651844.0,
    69997945.0 / 29380423.0,
];

const SAFETY: f64 = 0.9;
const MIN_FACTOR: f64 = 0.2;
const MAX_FACTOR: f64 = 10.0;

/// The dense output of one accepted step over [t, t + h].
#[derive(Debug, Clone, Default)]
pub(crate) struct Dense {
    pub t: f64,
    pub h: f64,
    r: [Vec<f64>; 5],
}

/// The continuous extension at θ in [0, 1] from the five coefficients of
/// one component.
#[inline]
pub(crate) fn interpolate(r: [f64; 5], theta: f64) -> f64 {
    let other = 1.0 - theta;
    r[0] + theta * (r[1] + other * (r[2] + theta * (r[3] + other * r[4])))
}

impl Dense {
    /// Component `i` at θ in [0, 1].
    pub fn component(&self, i: usize, theta: f64) -> f64 {
        interpolate(self.coefficients(i), theta)
    }

    /// The five coefficients of component `i`.
    pub fn coefficients(&self, i: usize) -> [f64; 5] {
        let [r1, r2, r3, r4, r5] = &self.r;
        [r1[i], r2[i], r3[i], r4[i], r5[i]]
    }

    /// The whole state at θ.
    pub fn state(&self, theta: f64, out: &mut [f64]) {
        for (i, x) in out.iter_mut().enumerate() {
            *x = self.component(i, theta);
        }
    }

    /// The smallest θ in (0, 1] found by the Illinois method with
    /// `component(i, θ) >= target`, given `component(i, 0) < target` and
    /// `component(i, 1) >= target`: the right end of the final bracket.
    pub fn crossing_up(&self, i: usize, target: f64) -> f64 {
        self.crossing(|theta| self.component(i, theta) - target)
    }

    /// As [`Dense::crossing_up`] for a decreasing crossing of `target`.
    pub fn crossing_down(&self, i: usize, target: f64) -> f64 {
        self.crossing(|theta| target - self.component(i, theta))
    }

    fn crossing(&self, g: impl Fn(f64) -> f64) -> f64 {
        let (mut a, mut b) = (0.0f64, 1.0f64);
        let (mut ga, mut gb) = (g(a), g(b));
        if gb < 0.0 || ga >= 0.0 {
            return if ga >= 0.0 { 0.0 } else { 1.0 };
        }
        let mut side = 0i8;
        for _ in 0..200 {
            if b - a <= 2.0 * f64::EPSILON * b {
                break;
            }
            // Regula falsi with the Illinois modification, safeguarded by
            // bisection when the secant point leaves the bracket.
            let mut c = b - gb * (b - a) / (gb - ga);
            if !(c > a && c < b) {
                c = 0.5 * (a + b);
            }
            let gc = g(c);
            if gc >= 0.0 {
                b = c;
                gb = gc;
                if side == 1 {
                    ga *= 0.5;
                }
                side = 1;
            } else {
                a = c;
                ga = gc;
                if side == -1 {
                    gb *= 0.5;
                }
                side = -1;
            }
            if gb == 0.0 {
                break;
            }
        }
        b
    }
}

/// Work arrays of the integrator for a state of fixed dimension.
pub(crate) struct Stepper {
    pub rtol: f64,
    pub atol: f64,
    k: [Vec<f64>; 7],
    stage: Vec<f64>,
    pub next: Vec<f64>,
    error: Vec<f64>,
    /// Whether `k[0]` holds f(t, y) for the current state (first same as
    /// last after an accepted step).
    pub fresh: bool,
}

/// The outcome of one attempted step.
pub(crate) enum Attempt {
    /// Accepted; `next` holds y(t + h) and the dense output is filled.
    Accepted {
        /// The suggested size of the following step.
        next_h: f64,
    },
    /// Rejected; retry with this smaller step.
    Rejected(f64),
}

impl Stepper {
    pub fn new(dimension: usize, rtol: f64, atol: f64) -> Self {
        let zeros = || vec![0.0; dimension];
        Self {
            rtol,
            atol,
            k: std::array::from_fn(|_| zeros()),
            stage: zeros(),
            next: zeros(),
            error: zeros(),
            fresh: false,
        }
    }

    fn norm(&self, y: &[f64], other: &[f64], values: &[f64]) -> f64 {
        let n = y.len() as f64;
        let sum: f64 = values
            .iter()
            .zip(y.iter().zip(other))
            .map(|(e, (a, b))| {
                let scale = self.atol + self.rtol * a.abs().max(b.abs());
                (e / scale) * (e / scale)
            })
            .sum();
        libm::sqrt(sum / n)
    }

    /// f(t, y), evaluated unless the last accepted step left it.
    pub fn derivative<F: FnMut(f64, &[f64], &mut [f64]) -> Result<(), H> + ?Sized, H>(
        &mut self,
        t: f64,
        y: &[f64],
        f: &mut F,
    ) -> Result<&[f64], H> {
        if !self.fresh {
            f(t, y, &mut self.k[0])?;
            self.fresh = true;
        }
        Ok(&self.k[0])
    }

    /// Hairer, Nørsett and Wanner's starting step (section II.4) for a
    /// fifth-order method, at most `limit`.
    pub fn initial_step<F: FnMut(f64, &[f64], &mut [f64]) -> Result<(), H> + ?Sized, H>(
        &mut self,
        t: f64,
        y: &[f64],
        limit: f64,
        f: &mut F,
    ) -> Result<f64, H> {
        if !self.fresh {
            f(t, y, &mut self.k[0])?;
            self.fresh = true;
        }
        let d0 = self.norm(y, y, y);
        let d1 = self.norm(y, y, &self.k[0]);
        let mut h0 = if d0 < 1e-5 || d1 < 1e-5 {
            1e-6
        } else {
            0.01 * d0 / d1
        };
        h0 = h0.min(limit);
        for (s, (x, k)) in self.stage.iter_mut().zip(y.iter().zip(&self.k[0])) {
            *s = x + h0 * k;
        }
        let stage = std::mem::take(&mut self.stage);
        f(t + h0, &stage, &mut self.k[1])?;
        self.stage = stage;
        for (e, (a, b)) in self.error.iter_mut().zip(self.k[1].iter().zip(&self.k[0])) {
            *e = a - b;
        }
        let d2 = self.norm(y, y, &self.error) / h0;
        let largest = d1.max(d2);
        let h1 = if largest <= 1e-15 {
            (h0 * 1e-3).max(1e-6)
        } else {
            libm::pow(0.01 / largest, 1.0 / 5.0)
        };
        Ok((100.0 * h0).min(h1).min(limit))
    }

    /// One Dormand–Prince step of size `h` from (t, y). On acceptance the
    /// dense output is written to `dense` and `next` holds the new state;
    /// `k[0]` then holds f(t + h, next).
    pub fn attempt<F: FnMut(f64, &[f64], &mut [f64]) -> Result<(), H> + ?Sized, H>(
        &mut self,
        t: f64,
        y: &[f64],
        h: f64,
        f: &mut F,
        dense: &mut Dense,
    ) -> Result<Attempt, H> {
        let n = y.len();
        if !self.fresh {
            f(t, y, &mut self.k[0])?;
            self.fresh = true;
        }
        for s in 1..7 {
            for (i, (stage, &x)) in self.stage.iter_mut().zip(y).enumerate() {
                let mut sum = 0.0;
                for (j, a) in A[s].iter().enumerate().take(s) {
                    sum += a * self.k[j][i];
                }
                *stage = x + h * sum;
            }
            let stage = std::mem::take(&mut self.stage);
            let result = f(t + C[s] * h, &stage, &mut self.k[s]);
            self.stage = stage;
            result?;
        }
        // The seventh stage is evaluated at y(t + h): the stage vector.
        self.next.copy_from_slice(&self.stage);
        for i in 0..n {
            let mut sum = 0.0;
            for (s, e) in E.iter().enumerate() {
                sum += e * self.k[s][i];
            }
            self.error[i] = h * sum;
        }
        let err = self.norm(y, &self.next, &self.error);
        if !err.is_finite() {
            return Ok(Attempt::Rejected(h * MIN_FACTOR));
        }
        if err > 1.0 {
            let factor = (SAFETY * libm::pow(err, -0.2)).max(MIN_FACTOR);
            return Ok(Attempt::Rejected(h * factor.min(1.0)));
        }
        dense.t = t;
        dense.h = h;
        for r in &mut dense.r {
            r.resize(n, 0.0);
        }
        for (i, &x) in y.iter().enumerate() {
            let difference = self.next[i] - x;
            let spline = h * self.k[0][i] - difference;
            let mut d = 0.0;
            for (s, c) in D.iter().enumerate() {
                d += c * self.k[s][i];
            }
            dense.r[0][i] = x;
            dense.r[1][i] = difference;
            dense.r[2][i] = spline;
            dense.r[3][i] = difference - h * self.k[6][i] - spline;
            dense.r[4][i] = h * d;
        }
        // First same as last: k7 is f(t + h, next).
        self.k.swap(0, 6);
        let factor = if err == 0.0 {
            MAX_FACTOR
        } else {
            (SAFETY * libm::pow(err, -0.2)).clamp(MIN_FACTOR, MAX_FACTOR)
        };
        Ok(Attempt::Accepted { next_h: h * factor })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Rhs<'a> = &'a mut dyn FnMut(f64, &[f64], &mut [f64]) -> Result<(), ()>;

    fn integrate(f: Rhs<'_>, y0: &[f64], end: f64, rtol: f64, atol: f64) -> (Vec<f64>, usize) {
        let mut stepper = Stepper::new(y0.len(), rtol, atol);
        let mut dense = Dense::default();
        let mut y = y0.to_vec();
        let mut t = 0.0;
        let mut h = stepper.initial_step(t, &y, end, f).unwrap();
        let mut steps = 0;
        while t < end {
            h = h.min(end - t);
            steps += 1;
            match stepper.attempt(t, &y, h, f, &mut dense).unwrap() {
                Attempt::Accepted { next_h } => {
                    t += h;
                    y.copy_from_slice(&stepper.next);
                    h = next_h;
                }
                Attempt::Rejected(smaller) => {
                    stepper.fresh = true;
                    h = smaller;
                }
            }
        }
        (y, steps)
    }

    /// The tableau is consistent: rows sum to the nodes, and the weights of
    /// both solutions sum to one.
    #[test]
    fn tableau_is_consistent() {
        for s in 1..7 {
            let row: f64 = A[s].iter().sum();
            assert!((row - C[s]).abs() < 1e-15, "row {s}");
        }
        let weights: f64 = A[6].iter().sum();
        assert!((weights - 1.0).abs() < 1e-15);
        let error: f64 = E.iter().sum();
        assert!(error.abs() < 1e-15);
        // The dense output weights vanish at both ends of the step.
        let d: f64 = D.iter().sum();
        assert!(d.abs() < 1e-9, "{d}");
    }

    /// Fifth-order convergence on y' = -y at fixed steps, and the error
    /// control reaches the requested tolerance.
    #[test]
    fn order_and_tolerance() {
        let mut f = |_t: f64, y: &[f64], dy: &mut [f64]| -> Result<(), ()> {
            dy[0] = -y[0];
            Ok(())
        };
        let mut fixed = |h: f64| {
            let mut stepper = Stepper::new(1, 1e300, 1e300);
            let mut dense = Dense::default();
            let mut y = vec![1.0];
            let mut t = 0.0;
            while t < 1.0 - 1e-12 {
                let Attempt::Accepted { .. } =
                    stepper.attempt(t, &y, h, &mut f, &mut dense).unwrap()
                else {
                    panic!("rejected")
                };
                y.copy_from_slice(&stepper.next);
                t += h;
            }
            (y[0] - libm::exp(-1.0)).abs()
        };
        let ratio = fixed(0.1) / fixed(0.05);
        assert!((ratio.log2() - 5.0).abs() < 0.3, "order {}", ratio.log2());
        for rtol in [1e-4, 1e-7, 1e-10] {
            let (y, _) = integrate(&mut f, &[1.0], 10.0, rtol, rtol * 1e-3);
            let exact = libm::exp(-10.0);
            assert!(
                (y[0] - exact).abs() <= 100.0 * rtol * exact,
                "{rtol}: {}",
                y[0]
            );
        }
    }

    /// The dense output matches the step's ends, has local error of order
    /// five in the step size, and crossings are located on it.
    #[test]
    fn dense_output_and_crossings() {
        let mut f = |_t: f64, y: &[f64], dy: &mut [f64]| -> Result<(), ()> {
            dy[0] = y[1];
            dy[1] = -y[0];
            Ok(())
        };
        // Tolerances that accept any step.
        let step = |h: f64, f: Rhs<'_>| {
            let mut stepper = Stepper::new(2, 1e300, 1e300);
            let mut dense = Dense::default();
            let Attempt::Accepted { .. } =
                stepper.attempt(0.0, &[0.0, 1.0], h, f, &mut dense).unwrap()
            else {
                panic!("rejected")
            };
            (dense, stepper.next.clone())
        };
        let worst = |h: f64, f: Rhs<'_>| {
            let (dense, _) = step(h, f);
            (1..10)
                .map(|k| {
                    let theta = f64::from(k) / 10.0;
                    (dense.component(0, theta) - libm::sin(theta * h)).abs()
                })
                .fold(0.0, f64::max)
        };
        let order = (worst(0.4, &mut f) / worst(0.2, &mut f)).log2();
        assert!((order - 5.0).abs() < 0.5, "dense output order {order}");
        let h = 0.4;
        let (dense, next) = step(h, &mut f);
        assert_eq!(dense.component(0, 0.0), 0.0);
        assert_eq!(dense.component(0, 1.0), next[0]);
        // sin t = 0.25 at asin(0.25) = 0.2526802551420786.
        let theta = dense.crossing_up(0, 0.25);
        assert!((theta * h - 0.252_680_255_142_078_6).abs() < 1e-5);
        assert!(dense.component(0, theta) >= 0.25);
        assert!(dense.component(0, theta) - 0.25 < 1e-15);
        // cos t = 0.95 at acos(0.95) = 0.3175604292915215.
        let theta = dense.crossing_down(1, 0.95);
        assert!((theta * h - 0.317_560_429_291_521_5).abs() < 1e-5);
        assert!(dense.component(1, theta) <= 0.95);
    }
}
