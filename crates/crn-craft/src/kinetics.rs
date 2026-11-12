//! Deterministic count kinetics, not exact nonlinear stochastic means.
//! Mass action is coefficient * product(y^order/order!), the macroscopic law.
//!
//! [`Kinetics`] evaluates the reaction-rate ODE dy/dt = Σ changes ·
//! propensity(y, p) in molecule counts, with its state Jacobian and its
//! Jacobian with respect to the named parameters, for the optional diffsol
//! BDF integration. It reads any closed model through its flat view
//! ([`FlatView`]): a network `crn-craft.model/0.1` (converted first), a rules
//! model or a `crn-craft.model/1.0`. Species are the view's species, the
//! state `y` holds their counts as real numbers, and each channel of the
//! view contributes its net change times its propensity:
//!
//! - a rule with a `rate`: rate × factor × Π yᵢ^kᵢ / kᵢ! over the species
//!   the channel takes (kᵢ molecules of species i). This is the macroscopic
//!   form of the combination count Π C(nᵢ, kᵢ) the stochastic engine uses;
//!   for kᵢ = 1 the two agree;
//! - a rule with a whole `propensity`: the expression × factor, evaluated at
//!   real counts. Its availability gate (every taken count present) is a
//!   step with no smooth extension, so it is dropped, and a whole-rule
//!   propensity that takes molecules must be 0 whenever a taken species is
//!   0, as the Michaelis–Menten form `vmax * (S / (km + S))` is; otherwise
//!   the model is refused. Hill laws (`basal + vmax / (1 + (k / R)^n)`) are
//!   accepted on production-only reactions.
//!
//! Expressions are evaluated with dual numbers ([`crate::rules`]'s
//! expression language, all operators and functions): the values are the
//! engine's, and the Jacobians are the dual parts. Where an expression has
//! no derivative, the derivative of the branch the evaluation selects is
//! used (`if`, `min`, `max`, `abs`), and `floor`, `ceil`, comparisons and
//! logic have derivative 0. Observables are sums of weighted counts (with
//! derivative weight) or presence-based reductions (`distinct`, `max`,
//! `min`: a species is present when its count is above 0; derivative 0).
//! Hill laws at a regulator count of 0 take their one-sided derivatives (0
//! for n > 1, vmax / k for n = 1); for n < 1 the derivative is singular.
//!
//! For nonlinear networks the solution differs from the mean of the
//! stochastic process. Failures are `CRN-E080` at `/ode`, or at the rule for
//! a whole-rule propensity that cannot be extended.
//!
//! ```
//! use crn_craft::kinetics::Kinetics;
//!
//! // A network model: decay at rate k per molecule.
//! let kinetics = Kinetics::from_bytes(br#"{"schema": "crn-craft.model/0.1", "name": "decay",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {"k": 0.5},
//!   "species": [{"name": "A", "initial_count": 10}],
//!   "reactions": [{"name": "death", "reactants": {"A": 1}, "products": {},
//!                  "rate": {"kind": "mass_action", "coefficient": "k"}}]}"#).unwrap();
//! let d = kinetics.evaluate(&[4.0], &[0.5]).unwrap();
//! assert_eq!(d.rhs, vec![-2.0]); // -k * y
//! assert_eq!(d.state_jacobian, vec![vec![-0.5]]); // -k
//! assert_eq!(d.parameter_jacobian, vec![vec![-4.0]]); // -y
//!
//! // A crn-craft.model/1.0: 2A -> B at rate k per pair, propensity k A^2 / 2.
//! let kinetics = Kinetics::from_bytes(br#"{"schema": "crn-craft.model/1.0", "name": "dimer",
//!   "parameters": {"k": 0.5}, "types": [{"name": "A"}, {"name": "B"}],
//!   "initial": [{"type": "A", "count": 10}],
//!   "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
//!   "rules": [{"name": "bind", "reactants": [{"label": "x", "type": "A", "count": 2}],
//!              "rate": "k", "products": [{"type": "B"}]}]}"#).unwrap();
//! assert_eq!(kinetics.species(), ["A", "B"]);
//! let d = kinetics.evaluate(&[4.0, 0.0], &[0.5]).unwrap();
//! assert_eq!(d.rhs, vec![-8.0, 4.0]); // -2 k A^2 / 2, k A^2 / 2
//! assert_eq!(d.state_jacobian[0][0], -4.0); // -2 k A
//! ```
use crate::{
    contracts::Model,
    diagnostic::Diagnostic,
    rules::{
        FlatView,
        dual::{Evaluation, Input, Program},
        model::Reduce,
    },
};

/// The largest number of molecules of one species a channel with a `rate`
/// may take (its order in that species).
const MAX_ORDER: u64 = 8;

/// A model prepared for deterministic evaluation. Supports at most 64
/// species, 32 parameters and 256 channels (reactions).
#[derive(Clone)]
pub struct Kinetics {
    /// Parameter names in sorted order: the order of `p` and of Jacobian columns.
    pub parameters: Vec<String>,
    values: Vec<f64>,
    species: Vec<String>,
    initial: Vec<u64>,
    rules: Vec<Rule>,
    observables: Vec<Observable>,
    channels: Vec<Channel>,
}

#[derive(Clone)]
struct Rule {
    name: String,
    program: Program,
}

#[derive(Clone)]
struct Observable {
    name: String,
    reduce: Reduce,
    terms: Vec<(usize, f64)>,
}

#[derive(Clone)]
struct Channel {
    rule: usize,
    whole: bool,
    factor: f64,
    reactants: Vec<(usize, u64)>,
    changes: Vec<(usize, i64)>,
}

/// Right-hand side and Jacobians at one state, all indexed by species order.
pub struct Derivatives {
    /// dy/dt in counts per second, one entry per species.
    pub rhs: Vec<f64>,
    /// ∂(dyᵢ/dt)/∂yⱼ as rows `i` and columns `j`.
    pub state_jacobian: Vec<Vec<f64>>,
    /// ∂(dyᵢ/dt)/∂pⱼ as rows `i` (species) and columns `j` in [`Kinetics::parameters`] order.
    pub parameter_jacobian: Vec<Vec<f64>>,
}

const HINT: &str = "Use supported deterministic kinetics and finite non-negative states/parameters; ODE trajectories are not exact nonlinear SSA means.";

fn error(s: impl Into<String>) -> Diagnostic {
    Diagnostic::error("CRN-E080", "/ode", s, HINT)
}

/// y^k / k! and y^(k−1) / (k−1)!, as running products of y / j.
#[allow(clippy::cast_precision_loss)]
fn macroscopic(y: f64, k: u64) -> (f64, f64) {
    let (mut previous, mut product) = (1.0, 1.0);
    for j in 1..=k {
        previous = product;
        product *= y / j as f64;
    }
    (product, previous)
}

impl Kinetics {
    /// Prepares the kinetics of a model's flat view: its species, channels,
    /// rate and propensity expressions, observables and parameters.
    ///
    /// # Errors
    ///
    /// `CRN-E080` at `/ode` when the view has more than 64 species, 32
    /// parameters or 256 channels, or a channel with a `rate` takes more
    /// than eight molecules of one species; `CRN-E080` at
    /// `/rules/i/propensity` for a whole-rule propensity that takes
    /// molecules and is not 0 whenever a taken species is 0 (judged from
    /// the expression's structure: a product with that species' observable
    /// as a factor, for example).
    pub fn new(view: &FlatView) -> Result<Self, Vec<Diagnostic>> {
        let parameters: Vec<String> = view
            .parameter_names()
            .into_iter()
            .map(String::from)
            .collect();
        if view.species().len() > 64 || parameters.len() > 32 || view.channels() > 256 {
            return Err(vec![error(
                "ODE size cap: 64 species, 32 parameters, 256 channels (reactions)",
            )]);
        }
        let rules: Vec<Rule> = (0..view.rules())
            .map(|r| Rule {
                name: view.rule_name(r).to_string(),
                program: Program::new(view.rate_expression(r)),
            })
            .collect();
        let observables: Vec<Observable> = view
            .observable_names()
            .into_iter()
            .enumerate()
            .map(|(o, name)| {
                let (reduce, terms) = view.observable(o);
                Observable {
                    name: name.to_string(),
                    reduce,
                    terms,
                }
            })
            .collect();
        let channels: Vec<Channel> = (0..view.channels())
            .map(|c| Channel {
                rule: view.rule(c),
                whole: view.whole(c),
                factor: view.factor(c),
                reactants: view.reactants(c).to_vec(),
                changes: view.changes(c).to_vec(),
            })
            .collect();
        let mut out = Vec::new();
        for channel in &channels {
            let rule = &rules[channel.rule];
            if !channel.whole {
                if let Some(&(s, k)) = channel.reactants.iter().find(|&&(_, k)| k > MAX_ORDER) {
                    out.push(error(format!(
                        "rule `{}` takes {k} molecules of {}: ODE mass-action order above eight is unsupported",
                        rule.name,
                        view.species()[s]
                    )));
                }
                continue;
            }
            for &(s, _) in &channel.reactants {
                // An observable is 0 when species s is 0 if all its terms
                // are on s (an observable without terms is always 0).
                let zero = |o: usize| observables[o].terms.iter().all(|&(t, _)| t == s);
                if !rule.program.vanishes(&zero) {
                    out.push(Diagnostic::error(
                        "CRN-E080",
                        &format!("/rules/{}/propensity", channel.rule),
                        format!(
                            "the propensity of rule `{}` takes molecules of {} but is not 0 when that count is 0; its availability gate has no smooth extension",
                            rule.name,
                            view.species()[s]
                        ),
                        "Write the propensity with the taken species' observable as a factor (as vmax * (S / (km + S))), or make the rule production-only.",
                    ));
                    break;
                }
            }
        }
        if !out.is_empty() {
            out.dedup_by(|a, b| a.path == b.path && a.message == b.message);
            return Err(out);
        }
        Ok(Self {
            parameters,
            values: view.parameter_values().to_vec(),
            species: view.species().to_vec(),
            initial: view.initial().to_vec(),
            rules,
            observables,
            channels,
        })
    }

    /// Reads any model the engine reads and prepares its kinetics: a
    /// network `crn-craft.model/0.1` as [`Self::compile`] does, otherwise
    /// the view of [`FlatView::of_any_model`] (every species of a type
    /// without attributes, also one that stays at 0).
    ///
    /// # Errors
    ///
    /// The errors of [`Self::compile`] for a network model, the errors of
    /// [`FlatView::of_any_model`] (parse and validation errors, `CRN-R017`
    /// for a model without a flat view), and those of [`Self::new`].
    pub fn from_bytes(model_bytes: &[u8]) -> Result<Self, Vec<Diagnostic>> {
        if crate::schema_field(model_bytes).as_deref() == Some("crn-craft.model/0.1") {
            return Self::compile(&crate::parse_model(model_bytes).map_err(|d| vec![d])?);
        }
        Self::new(&FlatView::of_any_model(model_bytes)?)
    }

    /// Prepares the kinetics of a parsed network model (the 0.6.0 entry
    /// point): the model is converted to `crn-craft.model/1.0` and read
    /// through its flat view, so its species and parameters keep their
    /// order and its reactions become one channel each. Population limits
    /// are safety stops of a simulation; the kinetics never read them, so
    /// they are left out of the conversion (which refuses them).
    ///
    /// # Errors
    ///
    /// The diagnostics of [`crate::rules::convert::from_network`] (the
    /// model's validation errors, `CRN-R016` for a model without reactions
    /// or with an initial count above 2^53) and those of [`Self::new`].
    pub fn compile(model: &Model) -> Result<Self, Vec<Diagnostic>> {
        let mut model = model.clone();
        for species in &mut model.species {
            species.population_limit = None;
        }
        let bytes = serde_json::to_vec(&model).expect("a parsed model serializes");
        let converted = crate::rules::convert::from_network(&bytes)?;
        Self::new(&FlatView::of_any_model(&converted)?)
    }

    /// Species labels in state order (the view's species; for a network
    /// model, its species in declaration order).
    #[must_use]
    pub fn species(&self) -> &[String] {
        &self.species
    }

    /// The model's initial counts in state order.
    #[must_use]
    pub fn initial(&self) -> &[u64] {
        &self.initial
    }

    /// The model's parameter values in [`Self::parameters`] order.
    #[must_use]
    pub fn parameter_values(&self) -> &[f64] {
        &self.values
    }

    /// Evaluates the right-hand side and both Jacobians at counts `y` (state
    /// order, finite and non-negative) and parameter values `p` (finite, in
    /// [`Self::parameters`] order).
    ///
    /// # Errors
    ///
    /// `CRN-E080` for wrong lengths or invalid inputs, an observable that
    /// overflows, a propensity that is negative or not finite, a derivative
    /// that is not finite (a Hill law with n < 1 at a regulator count of 0
    /// is singular at zero), or an accumulation that overflows.
    pub fn evaluate(&self, y: &[f64], p: &[f64]) -> Result<Derivatives, Diagnostic> {
        let n = self.species.len();
        let m = self.parameters.len();
        if y.len() != n
            || p.len() != m
            || y.iter().any(|v| !v.is_finite() || *v < 0.0)
            || p.iter().any(|v| !v.is_finite())
        {
            return Err(error("invalid ODE state or parameter vector"));
        }
        let mut values = Vec::with_capacity(self.observables.len());
        for observable in &self.observables {
            let present = observable.terms.iter().filter(|&&(s, _)| y[s] > 0.0);
            let value = match observable.reduce {
                Reduce::Sum => observable.terms.iter().map(|&(s, w)| y[s] * w).sum(),
                Reduce::Distinct => present.map(|&(_, w)| w).sum(),
                Reduce::Max => present.map(|&(_, w)| w).reduce(f64::max).unwrap_or(0.0),
                Reduce::Min => present.map(|&(_, w)| w).reduce(f64::min).unwrap_or(0.0),
            };
            if !value.is_finite() {
                return Err(error(format!("observable `{}` overflows", observable.name)));
            }
            values.push(value);
        }
        let mut evaluations: Vec<Option<Evaluation>> = vec![None; self.rules.len()];
        let mut result = Derivatives {
            rhs: vec![0.; n],
            state_jacobian: vec![vec![0.; n]; n],
            parameter_jacobian: vec![vec![0.; m]; n],
        };
        let mut dy = vec![0.; n];
        let mut dp = vec![0.; m];
        for channel in &self.channels {
            let rule = &self.rules[channel.rule];
            let e =
                evaluations[channel.rule].get_or_insert_with(|| rule.program.evaluate(p, &values));
            dy.fill(0.);
            dp.fill(0.);
            // The propensity is e.value * scale; the combination term's own
            // derivatives go straight to dy.
            let scale = if channel.whole {
                1.0 * channel.factor
            } else {
                let factors: Vec<(f64, f64)> = channel
                    .reactants
                    .iter()
                    .map(|&(s, k)| macroscopic(y[s], k))
                    .collect();
                let combination = factors.iter().fold(1.0, |x, f| x * f.0);
                for (i, &(s, _)) in channel.reactants.iter().enumerate() {
                    let others = factors
                        .iter()
                        .enumerate()
                        .filter(|&(j, _)| j != i)
                        .fold(1.0, |x, (_, f)| x * f.0);
                    dy[s] += e.value * channel.factor * (others * factors[i].1);
                }
                combination * channel.factor
            };
            let value = e.value * scale;
            for (input, g) in rule.program.inputs().iter().zip(&e.grad) {
                match *input {
                    Input::Parameter(j) => dp[j] += g * scale,
                    Input::Observable(o) => {
                        let observable = &self.observables[o];
                        if observable.reduce == Reduce::Sum {
                            for &(s, w) in &observable.terms {
                                dy[s] += g * scale * w;
                            }
                        }
                    }
                }
            }
            if !value.is_finite() {
                return Err(error(format!(
                    "the propensity of rule `{}` is not finite at this state",
                    rule.name
                )));
            }
            if value < 0.0 {
                return Err(error(format!(
                    "the propensity of rule `{}` is negative at this state",
                    rule.name
                )));
            }
            if dy.iter().chain(&dp).any(|v| !v.is_finite()) {
                return Err(error(if e.singular {
                    format!(
                        "the derivative of rule `{}` is singular at zero (a power with an exponent below 1, a square root or a logarithm at 0)",
                        rule.name
                    )
                } else {
                    format!(
                        "the derivative of rule `{}` is not finite at this state",
                        rule.name
                    )
                }));
            }
            for &(i, change) in &channel.changes {
                #[allow(clippy::cast_precision_loss)]
                let c = change as f64;
                result.rhs[i] += c * value;
                for (out, d) in result.state_jacobian[i].iter_mut().zip(&dy) {
                    *out += c * d;
                }
                for (out, d) in result.parameter_jacobian[i].iter_mut().zip(&dp) {
                    *out += c * d;
                }
            }
        }
        if result
            .rhs
            .iter()
            .chain(result.state_jacobian.iter().flatten())
            .chain(result.parameter_jacobian.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(error("ODE derivative accumulation overflows"));
        }
        Ok(result)
    }
}
