//! Deterministic count kinetics, not exact nonlinear stochastic means.
//! Mass action is coefficient * product(y^order/order!), the macroscopic law.
//!
//! [`Kinetics`] evaluates the reaction-rate ODE dy/dt = N·a(y, p) in molecule
//! counts, with its state Jacobian and its Jacobian with respect to the named
//! parameters, for the optional diffsol BDF integration. Hill and
//! Michaelis–Menten laws use their stochastic formulas with continuous
//! counts. For nonlinear networks the solution differs from the mean of the
//! stochastic process. Failures are `CRN-E080` at `/ode`.
//!
//! ```
//! use crn_craft::{kinetics::Kinetics, parse_model};
//!
//! let model = parse_model(br#"{"schema": "crn-craft.model/0.1", "name": "decay",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {"k": 0.5},
//!   "species": [{"name": "A", "initial_count": 10}],
//!   "reactions": [{"name": "death", "reactants": {"A": 1}, "products": {},
//!                  "rate": {"kind": "mass_action", "coefficient": "k"}}]}"#).unwrap();
//! let kinetics = Kinetics::compile(&model).unwrap();
//! let d = kinetics.evaluate(&[4.0], &[0.5]).unwrap();
//! assert_eq!(d.rhs, vec![-2.0]); // -k * y
//! assert_eq!(d.state_jacobian, vec![vec![-0.5]]); // -k
//! assert_eq!(d.parameter_jacobian, vec![vec![-4.0]]); // -y
//! ```
use crate::{
    CompiledNetwork,
    contracts::{Model, Rate, Scalar},
    diagnostic::Diagnostic,
};
use std::collections::BTreeMap;
/// A model prepared for deterministic evaluation. Supports at most 64
/// species, 32 parameters and 256 reactions.
#[derive(Clone)]
pub struct Kinetics {
    /// The source model; its parameter values are not used by [`Kinetics::evaluate`].
    pub model: Model,
    /// The compiled network, which supplies the net stoichiometry.
    pub network: CompiledNetwork,
    /// Parameter names in sorted order: the order of `p` and of Jacobian columns.
    pub parameters: Vec<String>,
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
fn error(s: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E080",
        "/ode",
        s,
        "Use supported deterministic kinetics and finite non-negative states/parameters; ODE trajectories are not exact nonlinear SSA means.",
    )
}
impl Kinetics {
    /// Compiles `model` for ODE evaluation.
    ///
    /// # Errors
    ///
    /// The diagnostics of [`CompiledNetwork::compile`], or `CRN-E080` when the
    /// model exceeds the size cap, a Hill law has reactants, a
    /// Michaelis–Menten law does not consume exactly one unit of its
    /// substrate and nothing else, or a mass-action reactant order exceeds 8.
    pub fn compile(model: &Model) -> Result<Self, Vec<Diagnostic>> {
        let network = CompiledNetwork::compile(model)?;
        if model.species.len() > 64 || model.parameters.len() > 32 || model.reactions.len() > 256 {
            return Err(vec![error(
                "ODE size cap: 64 species, 32 parameters, 256 reactions",
            )]);
        }
        for r in &model.reactions {
            match &r.rate {
                Rate::HillActivation { .. } | Rate::HillRepression { .. }
                    if !r.reactants.is_empty() =>
                {
                    return Err(vec![error(
                        "Hill ODE laws require production-only reactions; discrete reactant gates have no smooth extension",
                    )]);
                }
                Rate::MichaelisMenten { substrate, .. }
                    if r.reactants.len() != 1 || r.reactants.get(substrate) != Some(&1) =>
                {
                    return Err(vec![error("MM ODE laws must consume only one substrate")]);
                }
                Rate::MassAction { .. } if r.reactants.values().any(|&n| n > 8) => {
                    return Err(vec![error(
                        "ODE mass-action order above eight is unsupported",
                    )]);
                }
                _ => {}
            }
        }
        Ok(Self {
            model: model.clone(),
            network,
            parameters: model.parameters.keys().cloned().collect(),
        })
    }
    /// Evaluates the right-hand side and both Jacobians at counts `y` (species
    /// order) and parameter values `p` ([`Self::parameters`] order).
    ///
    /// # Errors
    ///
    /// `CRN-E080` for wrong lengths, negative or non-finite inputs, a Hill
    /// `k` or `n` or a Michaelis–Menten `km` that is not positive, a Hill law
    /// with n < 1 at a zero regulator, or a non-finite result.
    pub fn evaluate(&self, y: &[f64], p: &[f64]) -> Result<Derivatives, Diagnostic> {
        let n = self.model.species.len();
        let m = self.parameters.len();
        if y.len() != n || p.len() != m || y.iter().chain(p).any(|v| !v.is_finite() || *v < 0.0) {
            return Err(error("invalid ODE state or parameter vector"));
        }
        let indices: BTreeMap<_, _> = self
            .model
            .species
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.as_str(), i))
            .collect();
        let scalar = |s: &Scalar| -> (f64, Option<usize>) {
            match s {
                Scalar::Literal(x) => (*x, None),
                Scalar::Parameter(s) => {
                    let j = self.parameters.binary_search(s).unwrap();
                    (p[j], Some(j))
                }
            }
        };
        let mut result = Derivatives {
            rhs: vec![0.; n],
            state_jacobian: vec![vec![0.; n]; n],
            parameter_jacobian: vec![vec![0.; m]; n],
        };
        for (ri, r) in self.model.reactions.iter().enumerate() {
            let mut dy = vec![0.; n];
            let mut dp = vec![0.; m];
            let value = match &r.rate {
                Rate::MassAction { coefficient } => {
                    let (c, ci) = scalar(coefficient);
                    let mut product = 1.;
                    for (s, &order) in &r.reactants {
                        for k in 1..=order {
                            product *= y[indices[s.as_str()]] / k as f64;
                        }
                    }
                    if let Some(j) = ci {
                        dp[j] += product;
                    }
                    for (s, &order) in &r.reactants {
                        let index = indices[s.as_str()];
                        let mut derivative = c;
                        for (other, &degree) in &r.reactants {
                            let power = degree - i64::from(other == s);
                            for k in 1..=degree {
                                derivative /= k as f64;
                            }
                            derivative *= y[indices[other.as_str()]].powi(power as i32);
                        }
                        dy[index] = derivative * order as f64;
                    }
                    c * product
                }
                Rate::HillActivation {
                    regulator,
                    vmax,
                    basal,
                    k,
                    n: power,
                }
                | Rate::HillRepression {
                    regulator,
                    vmax,
                    basal,
                    k,
                    n: power,
                } => {
                    let i = indices[regulator.as_str()];
                    let x = y[i];
                    let (v, vi) = scalar(vmax);
                    let (b, bi) = scalar(basal);
                    let (k, ki) = scalar(k);
                    let (power, ni) = scalar(power);
                    if k <= 0. || power <= 0. {
                        return Err(error("Hill K and n must remain positive"));
                    }
                    let repression = matches!(r.rate, Rate::HillRepression { .. });
                    let (activation, dx, dk, dn) = if x == 0. {
                        if power < 1. {
                            return Err(error("Hill derivative is singular at zero for n<1"));
                        }
                        (0., if power == 1. { 1. / k } else { 0. }, 0., 0.)
                    } else {
                        let log = (x / k).ln();
                        let z = power * log;
                        let a = if z >= 0. {
                            1. / (1. + (-z).exp())
                        } else {
                            let t = z.exp();
                            t / (1. + t)
                        };
                        let factor = a * (1. - a);
                        (a, factor * power / x, -factor * power / k, factor * log)
                    };
                    let sign = if repression { -1. } else { 1. };
                    let fraction = if repression {
                        1. - activation
                    } else {
                        activation
                    };
                    dy[i] = v * sign * dx;
                    if let Some(j) = ki {
                        dp[j] += v * sign * dk;
                    }
                    if let Some(j) = ni {
                        dp[j] += v * sign * dn;
                    }
                    if let Some(j) = vi {
                        dp[j] += fraction;
                    }
                    if let Some(j) = bi {
                        dp[j] += 1.;
                    }
                    b + v * fraction
                }
                Rate::MichaelisMenten {
                    substrate,
                    vmax,
                    km,
                } => {
                    let i = indices[substrate.as_str()];
                    let x = y[i];
                    let (v, vi) = scalar(vmax);
                    let (k, ki) = scalar(km);
                    if k <= 0. {
                        return Err(error("MM Km must remain positive"));
                    }
                    let d = k + x;
                    dy[i] = v * (k / d) / d;
                    if let Some(j) = vi {
                        dp[j] += x / d;
                    }
                    if let Some(j) = ki {
                        dp[j] -= v * (x / d) / d;
                    }
                    v * (x / d)
                }
            };
            if !value.is_finite() || dy.iter().chain(&dp).any(|v| !v.is_finite()) {
                return Err(error("non-finite deterministic rate/derivative"));
            }
            for change in &self.network.reactions()[ri].net_changes {
                let i = change.species;
                let c = change.coefficient as f64;
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
