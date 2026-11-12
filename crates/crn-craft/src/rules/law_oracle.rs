//! The network engine of release 0.6.0 as a test oracle: its structure
//! (species, stoichiometry, reads and dependencies of `network.rs`) and its
//! rate laws (`rates.rs`), both deleted in 0.7.0. The converter, the flat
//! view, the analysis and the combination rate are checked against it. The
//! arithmetic is the 0.6.0 code, including its fast branches and their
//! multiplication order; only the memo table, the limits and the
//! diagnostics (here the message text) are left out. Unit tests include
//! this file as `crate::rules::law_oracle`, integration tests by path.
#![allow(clippy::cast_precision_loss, dead_code)]

use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

/// One reaction's structure: `(species, coefficient)` terms in species
/// order, as the network engine compiled them.
pub(crate) struct Reaction {
    pub(crate) name: String,
    pub(crate) reactants: Vec<(usize, i64)>,
    pub(crate) products: Vec<(usize, i64)>,
    pub(crate) net_changes: Vec<(usize, i64)>,
    pub(crate) reads: Vec<usize>,
}

/// A network model document read as the 0.6.0 engine compiled it.
pub(crate) struct Network {
    pub(crate) species: Vec<String>,
    pub(crate) initial: Vec<i64>,
    pub(crate) reactions: Vec<Reaction>,
    /// Per reaction, the reactions that read a species it changes.
    pub(crate) dependencies: Vec<Vec<usize>>,
    model: Value,
}

impl Network {
    /// Compiles a network model document (assumed valid).
    pub(crate) fn new(model: &Value) -> Self {
        let species: Vec<String> = model["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["name"].as_str().unwrap().to_string())
            .collect();
        let initial = model["species"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["initial_count"].as_i64().unwrap())
            .collect();
        let index = |name: &str| species.iter().position(|s| s == name).unwrap();
        let terms = |side: &Value| {
            let mut out: Vec<(usize, i64)> = side
                .as_object()
                .unwrap()
                .iter()
                .map(|(name, k)| (index(name), k.as_i64().unwrap()))
                .collect();
            out.sort_by_key(|&(species, _)| species);
            out
        };
        let mut readers = vec![Vec::new(); species.len()];
        let mut reactions = Vec::new();
        for (id, reaction) in model["reactions"].as_array().unwrap().iter().enumerate() {
            let reactants = terms(&reaction["reactants"]);
            let products = terms(&reaction["products"]);
            let mut changes = BTreeMap::new();
            for &(s, k) in &reactants {
                changes.insert(s, -k);
            }
            for &(s, k) in &products {
                *changes.entry(s).or_insert(0) += k;
            }
            let net_changes = changes.into_iter().filter(|&(_, k)| k != 0).collect();
            let mut reads: BTreeSet<usize> = reactants.iter().map(|&(s, _)| s).collect();
            let rate = &reaction["rate"];
            for key in ["regulator", "substrate"] {
                if let Some(name) = rate[key].as_str() {
                    reads.insert(index(name));
                }
            }
            for &s in &reads {
                readers[s].push(id);
            }
            reactions.push(Reaction {
                name: reaction["name"].as_str().unwrap().to_string(),
                reactants,
                products,
                net_changes,
                reads: reads.into_iter().collect(),
            });
        }
        let dependencies = reactions
            .iter()
            .map(|reaction| {
                let mut affected = BTreeSet::new();
                for &(s, _) in &reaction.net_changes {
                    affected.extend(readers[s].iter().copied());
                }
                affected.into_iter().collect()
            })
            .collect();
        Self {
            species,
            initial,
            reactions,
            dependencies,
            model: model.clone(),
        }
    }

    /// The 0.6.0 propensity of `reaction` at `state` (counts in species order).
    pub(crate) fn propensity(&self, reaction: usize, state: &[i64]) -> Result<f64, String> {
        let rate = &self.model["reactions"][reaction]["rate"];
        let parameters = &self.model["parameters"];
        let arg = |key: &str| match &rate[key] {
            Value::String(name) => parameters[name].as_f64().unwrap(),
            value => value.as_f64().unwrap(),
        };
        let reactants = &self.reactions[reaction].reactants;
        let index = |key: &str| {
            let name = rate[key].as_str().unwrap();
            self.species.iter().position(|s| s == name).unwrap()
        };
        let kind = rate["kind"].as_str().unwrap();
        if kind == "mass_action" {
            if let Some(value) = kernel(arg("coefficient"), reactants, state) {
                return Ok(value);
            }
        }
        // Gate before evaluation, including custom and constant propensities.
        if reactants.iter().any(|&(s, k)| state[s] < k) {
            return Ok(0.0);
        }
        match kind {
            "mass_action" => mass_action(arg("coefficient"), state, reactants),
            "hill_activation" | "hill_repression" => hill(
                state[index("regulator")] as f64,
                arg("vmax"),
                arg("basal"),
                arg("k"),
                arg("n"),
                kind == "hill_repression",
            ),
            "michaelis_menten" => {
                let x = state[index("substrate")] as f64;
                let (vmax, km) = (arg("vmax"), arg("km"));
                checked(vmax * (x / (km + x)), vmax > 0.0 && x > 0.0)
            }
            other => panic!("unknown rate law {other}"),
        }
    }
}

/// The 0.6.0 propensity of reaction 0 of a network model document at
/// `state` (counts in species order).
pub(crate) fn propensity(model: &Value, state: &[i64]) -> Result<f64, String> {
    Network::new(model).propensity(0, state)
}

/// `MassKernel::compile` and `MassKernel::value` of 0.6.0.
fn kernel(c: f64, reactants: &[(usize, i64)], state: &[i64]) -> Option<f64> {
    if reactants.is_empty() {
        return Some(c);
    }
    if !c.is_normal() {
        return None;
    }
    if c == 0.0 {
        return Some(0.0);
    }
    let value = match *reactants {
        [(s, 1)] => {
            let n = state[s];
            if n < 1 {
                return Some(0.0);
            }
            c * n as f64
        }
        [(s, 2)] => {
            let n = state[s];
            if n < 2 {
                return Some(0.0);
            }
            let first = if n == 2 {
                c
            } else if n == 3 {
                c * 3.0
            } else {
                c * (n - 1) as f64
            };
            let result = if n <= 3 {
                first
            } else {
                first * (n as f64 / 2.0)
            };
            if !first.is_normal() {
                return None;
            }
            result
        }
        [(s, 1), (t, 1)] => {
            if state[s] < 1 || state[t] < 1 {
                return Some(0.0);
            }
            let first = c * state[s] as f64;
            if !first.is_normal() {
                return None;
            }
            first * state[t] as f64
        }
        _ => return None,
    };
    value.is_normal().then_some(value)
}

fn hill(x: f64, vmax: f64, basal: f64, k: f64, n: f64, repression: bool) -> Result<f64, String> {
    let (value, positive) = if x == 0.0 {
        (
            basal + if repression { vmax } else { 0.0 },
            basal > 0.0 || (repression && vmax > 0.0),
        )
    } else if vmax == 0.0 {
        (basal, basal > 0.0)
    } else {
        // Near x=k, log1p preserves differences erased by log(x)-log(k).
        let relative = (x - k) / k;
        let log_ratio = if relative.abs() <= 0.5 {
            libm::log1p(relative)
        } else {
            libm::log(x) - libm::log(k)
        };
        let odds = n * log_ratio;
        let t = if repression { odds } else { -odds };
        let log_fraction = if t > 0.0 {
            -t - libm::log1p(libm::exp(-t))
        } else {
            -libm::log1p(libm::exp(t))
        };
        // Do not underflow the fraction before a large vmax rescales it.
        let induced = if log_fraction >= -700.0 {
            vmax * libm::exp(log_fraction)
        } else {
            libm::exp(libm::log(vmax) + log_fraction)
        };
        (basal + induced, true)
    };
    checked(value, positive)
}

fn checked(value: f64, positive: bool) -> Result<f64, String> {
    if !value.is_finite() || value < 0.0 {
        Err("non-finite or negative evaluated propensity".into())
    } else if positive && value == 0.0 {
        Err("positive propensity underflows to zero".into())
    } else {
        Ok(value)
    }
}

fn mass_action(coefficient: f64, state: &[i64], reactants: &[(usize, i64)]) -> Result<f64, String> {
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    // Common elementary rates: exact same rounded product when normal.
    // Exceptional/subnormal ranges retain the original scaled implementation.
    if reactants.is_empty() {
        return Ok(coefficient);
    }
    if coefficient.is_normal() {
        if let [(s, 1)] = *reactants {
            let value = coefficient * state[s] as f64;
            if value.is_normal() {
                return Ok(value);
            }
        }
        if let [(s, 2)] = *reactants {
            let count = state[s];
            let first = if count == 2 {
                coefficient
            } else if count == 3 {
                coefficient * 3.0
            } else {
                coefficient * (count - 1) as f64
            };
            let value = if count <= 3 {
                first
            } else {
                first * (count as f64 / 2.0)
            };
            if first.is_normal() && value.is_normal() {
                return Ok(value);
            }
        }
        if let [(s, 1), (t, 1)] = *reactants {
            let partial = coefficient * state[s] as f64;
            let value = partial * state[t] as f64;
            if partial.is_normal() && value.is_normal() {
                return Ok(value);
            }
        }
    }
    // Keep the coefficient in the scaled product: binomial(n,k) may overflow
    // independently even when coefficient*binomial(n,k) is representable.
    let (mut mantissa, mut exponent) = libm::frexp(coefficient);
    for &(species, k) in reactants {
        let count = state[species];
        let order = k.min(count - k);
        for i in 1..=order {
            let ratio = (count - order + i) as f64 / i as f64;
            let (fraction, shift) = libm::frexp(mantissa * ratio);
            mantissa = fraction;
            exponent += shift;
            // All remaining binomial factors are >=1, so this cannot recover.
            if exponent > 1024 {
                return Err("mass-action propensity overflows f64".into());
            }
        }
    }
    let value = libm::scalbn(mantissa, exponent);
    if !value.is_finite() || value <= 0.0 {
        Err("mass-action propensity is outside the positive f64 range".into())
    } else {
        Ok(value)
    }
}
