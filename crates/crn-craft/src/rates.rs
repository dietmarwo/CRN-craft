//! Resolved rates: no name lookup, parameter mutation or expression evaluation.
use crate::diagnostic::Diagnostic;
use crate::network::Term;

#[derive(Debug, Clone)]
pub(crate) enum CompiledRate {
    MassAction(f64),
    Hill {
        regulator: usize,
        vmax: f64,
        basal: f64,
        k: f64,
        n: f64,
        repression: bool,
    },
    MichaelisMenten {
        substrate: usize,
        vmax: f64,
        km: f64,
    },
}

fn range_error(reaction: usize, reason: &str) -> Diagnostic {
    Diagnostic::error(
        "CRN-E042",
        &format!("/reactions/{reaction}/rate"),
        reason,
        "Rescale or revise the model; a positive propensity must fit in f64 without becoming zero.",
    )
}

impl CompiledRate {
    pub(crate) fn extra_read(&self) -> Option<usize> {
        match self {
            Self::MassAction(_) => None,
            Self::Hill { regulator, .. } => Some(*regulator),
            Self::MichaelisMenten { substrate, .. } => Some(*substrate),
        }
    }

    pub(crate) fn evaluate(
        &self,
        state: &[i64],
        reactants: &[Term],
        reaction: usize,
    ) -> Result<f64, Diagnostic> {
        // Gate *before* evaluation, including custom and constant propensities.
        if reactants
            .iter()
            .any(|term| state[term.species] < term.coefficient)
        {
            return Ok(0.0);
        }
        let (value, positive) = match *self {
            Self::MassAction(coefficient) => {
                return mass_action(coefficient, state, reactants, reaction);
            }
            Self::Hill {
                regulator,
                vmax,
                basal,
                k,
                n,
                repression,
            } => {
                let x = state[regulator] as f64;
                if x == 0.0 {
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
                }
            }
            Self::MichaelisMenten {
                substrate,
                vmax,
                km,
            } => {
                let x = state[substrate] as f64;
                (vmax * (x / (km + x)), vmax > 0.0 && x > 0.0)
            }
        };
        if !value.is_finite() || value < 0.0 {
            Err(range_error(
                reaction,
                "non-finite or negative evaluated propensity",
            ))
        } else if positive && value == 0.0 {
            Err(range_error(
                reaction,
                "positive propensity underflows to zero",
            ))
        } else {
            Ok(value)
        }
    }
}

fn mass_action(
    coefficient: f64,
    state: &[i64],
    reactants: &[Term],
    reaction: usize,
) -> Result<f64, Diagnostic> {
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    // Keep the coefficient in the scaled product: binomial(n,k) may overflow
    // independently even when coefficient*binomial(n,k) is representable.
    let (mut mantissa, mut exponent) = libm::frexp(coefficient);
    for term in reactants {
        let count = state[term.species];
        let order = term.coefficient.min(count - term.coefficient);
        for i in 1..=order {
            let ratio = (count - order + i) as f64 / i as f64;
            let (fraction, shift) = libm::frexp(mantissa * ratio);
            mantissa = fraction;
            exponent += shift;
            // All remaining binomial factors are >=1, so this cannot recover.
            if exponent > 1024 {
                return Err(range_error(
                    reaction,
                    "mass-action propensity overflows f64",
                ));
            }
        }
    }
    let value = libm::scalbn(mantissa, exponent);
    if !value.is_finite() || value <= 0.0 {
        Err(range_error(
            reaction,
            "mass-action propensity is outside the positive f64 range",
        ))
    } else {
        Ok(value)
    }
}
