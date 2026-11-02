//! Resolved rates: no name lookup, parameter mutation or expression evaluation.
use crate::diagnostic::Diagnostic;
use crate::network::Term;
use std::sync::{
    Arc, OnceLock,
    atomic::{AtomicU64, Ordering},
};

/// Regulator counts below this bound share one lazily filled table per law.
pub(crate) const MEMO_COUNTS: usize = 8192;
/// At most this many Hill laws per network get a table (64 KiB each when used).
pub(crate) const MEMO_LAWS: usize = 128;
/// A NaN payload that a checked rate never returns.
const EMPTY: u64 = u64::MAX;

/// A one-variable Hill law is a pure function of an integer count, so a
/// computed, checked value can be returned again bit for bit. Failed
/// evaluations are not stored and are recomputed with their diagnostic.
#[derive(Default)]
pub(crate) struct Memo(OnceLock<Box<[AtomicU64]>>);
impl std::fmt::Debug for Memo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Memo")
    }
}

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
        memo: Option<Arc<Memo>>,
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
        match *self {
            Self::MassAction(coefficient) => mass_action(coefficient, state, reactants, reaction),
            Self::Hill {
                regulator,
                vmax,
                basal,
                k,
                n,
                repression,
                ref memo,
            } => {
                let count = state[regulator];
                let slot = memo.as_ref().and_then(|memo| {
                    let index = usize::try_from(count).ok().filter(|&i| i < MEMO_COUNTS)?;
                    let table = memo
                        .0
                        .get_or_init(|| (0..MEMO_COUNTS).map(|_| AtomicU64::new(EMPTY)).collect());
                    Some(&table[index])
                });
                if let Some(slot) = slot {
                    let bits = slot.load(Ordering::Relaxed);
                    if bits != EMPTY {
                        return Ok(f64::from_bits(bits));
                    }
                }
                let value = hill(count as f64, vmax, basal, k, n, repression, reaction)?;
                if let Some(slot) = slot {
                    slot.store(value.to_bits(), Ordering::Relaxed);
                }
                Ok(value)
            }
            Self::MichaelisMenten {
                substrate,
                vmax,
                km,
            } => {
                let x = state[substrate] as f64;
                checked(vmax * (x / (km + x)), vmax > 0.0 && x > 0.0, reaction)
            }
        }
    }
}

fn hill(
    x: f64,
    vmax: f64,
    basal: f64,
    k: f64,
    n: f64,
    repression: bool,
    reaction: usize,
) -> Result<f64, Diagnostic> {
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
    checked(value, positive, reaction)
}

fn checked(value: f64, positive: bool, reaction: usize) -> Result<f64, Diagnostic> {
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

fn mass_action(
    coefficient: f64,
    state: &[i64],
    reactants: &[Term],
    reaction: usize,
) -> Result<f64, Diagnostic> {
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    // Common elementary rates: exact same rounded product when normal.
    // Exceptional/subnormal ranges retain the original scaled implementation.
    if reactants.is_empty() {
        return Ok(coefficient);
    }
    if let [term] = reactants {
        if term.coefficient == 1 && coefficient.is_normal() {
            let value = coefficient * state[term.species] as f64;
            if value.is_normal() {
                return Ok(value);
            }
        }
    }
    if coefficient.is_normal() {
        if let [term] = reactants {
            if term.coefficient == 2 {
                let count = state[term.species];
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
        }
        if let [first, second] = reactants {
            if first.coefficient == 1 && second.coefficient == 1 {
                let partial = coefficient * state[first.species] as f64;
                let value = partial * state[second.species] as f64;
                if partial.is_normal() && value.is_normal() {
                    return Ok(value);
                }
            }
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn law(vmax: f64, basal: f64, k: f64, n: f64, repression: bool, memo: bool) -> CompiledRate {
        CompiledRate::Hill {
            regulator: 0,
            vmax,
            basal,
            k,
            n,
            repression,
            memo: memo.then(Default::default),
        }
    }

    #[test]
    fn hill_memo_returns_the_computed_bits_and_never_stores_failures() {
        let cases = [
            (216.0, 0.216, 40.0, 2.0, true),
            (3.0, 0.0, 20.0, 1.0, false),
            (1e300, 0.0, 7.5, 3.7, false),
            (5.0, 1e-300, 1.0, 60.0, true),
        ];
        for (vmax, basal, k, n, repression) in cases {
            let direct = law(vmax, basal, k, n, repression, false);
            let cached = law(vmax, basal, k, n, repression, true);
            // Twice: the first pass fills the table, the second reads it.
            for _ in 0..2 {
                for count in (0..MEMO_COUNTS as i64 + 40).chain([1 << 40]) {
                    let a = direct.evaluate(&[count], &[], 0);
                    let b = cached.evaluate(&[count], &[], 0);
                    match (a, b) {
                        (Ok(a), Ok(b)) => assert_eq!(a.to_bits(), b.to_bits(), "count {count}"),
                        (Err(a), Err(b)) => assert_eq!(a.message, b.message),
                        (a, b) => panic!("count {count}: {a:?} versus {b:?}"),
                    }
                }
            }
        }
        // Underflow is a diagnostic each time, from an empty slot.
        let failing = law(5e-324, 0.0, 2.0, 60.0, false, true);
        for _ in 0..2 {
            assert!(failing.evaluate(&[1], &[], 0).is_err());
        }
        let CompiledRate::Hill {
            memo: Some(memo), ..
        } = &failing
        else {
            unreachable!()
        };
        assert_eq!(memo.0.get().unwrap()[1].load(Ordering::Relaxed), EMPTY);
    }
}
