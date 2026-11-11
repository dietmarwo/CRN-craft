//! Combination counts of higher-order rules (`crn-craft.model/1.0`).
//!
//! A rule whose patterns take kᵢ molecules of types with nᵢ molecules fires
//! at rate × factor × Π C(nᵢ, kᵢ). The binomials can overflow `f64` even
//! when the product is representable (a tiny rate times a huge count), so
//! the product is accumulated as a mantissa and a binary exponent and
//! rounded once at the end. This is the scaled product of the network
//! engine's mass-action law (`crate::rates`, "binomial(n,k) may overflow"),
//! with the channel factor folded in.

/// rate × factor × Π C(nᵢ, kᵢ) over `reactants`, `(species, k)` pairs with
/// `counts[species] >= k` (checked by the caller).
///
/// Each binomial C(n, k) is the product of (n − m + i) / i for i = 1..m with
/// m = min(k, n − k), every ratio at least 1. The result is `+inf` as soon
/// as the exponent passes 1024, since no later ratio can bring it back, so
/// the work is bounded by about 3,200 ratios whatever k is. A NaN, infinite
/// or negative rate or factor carries through to the result, which the
/// caller rejects; a result below the smallest subnormal is 0.
///
/// Generated event loops call this function, so their arithmetic is the
/// interpreter's.
#[must_use]
pub fn combination_rate(rate: f64, factor: f64, counts: &[u64], reactants: &[(usize, u64)]) -> f64 {
    scaled(rate, factor, counts, reactants).0
}

/// The value and the number of ratios multiplied in.
#[allow(clippy::cast_precision_loss)]
fn scaled(rate: f64, factor: f64, counts: &[u64], reactants: &[(usize, u64)]) -> (f64, u64) {
    let (mut mantissa, mut exponent) = libm::frexp(rate);
    let (fraction, shift) = libm::frexp(factor);
    let (fraction, normal) = libm::frexp(mantissa * fraction);
    mantissa = fraction;
    exponent += shift + normal;
    let mut ratios = 0;
    for &(species, k) in reactants {
        let n = counts[species];
        let order = k.min(n - k);
        for i in 1..=order {
            // Counts are at most 2^53, so n - order + i and i are exact.
            let ratio = (n - order + i) as f64 / i as f64;
            let (fraction, shift) = libm::frexp(mantissa * ratio);
            mantissa = fraction;
            exponent += shift;
            ratios += 1;
            if exponent > 1024 {
                return (f64::INFINITY.copysign(mantissa), ratios);
            }
        }
    }
    (libm::scalbn(mantissa, exponent), ratios)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact(n: u64, k: u64) -> u128 {
        let k = k.min(n - k);
        (1..=u128::from(k)).fold(1u128, |acc, i| {
            acc.checked_mul(u128::from(n - k) + i)
                .expect("exact in u128")
                / i
        })
    }

    /// Declared in docs/campaign/merge-0.7/M4b-log.md: a binomial of k
    /// ratios is within k + 1 roundings of the exact value.
    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn binomials_are_exact_to_a_rounding_per_ratio() {
        for n in 0..=120u64 {
            for k in 0..=n {
                let value = combination_rate(1.0, 1.0, &[n], &[(0, k)]);
                let truth = exact(n, k) as f64;
                let order = k.min(n - k) as f64;
                assert!(
                    (value - truth).abs() <= (order + 1.0) * f64::EPSILON * truth,
                    "C({n}, {k}) = {value}, exact {truth}"
                );
                if k <= 2 || n - k <= 2 {
                    assert_eq!(value, truth, "C({n}, {k})");
                }
            }
        }
        // Two types: C(10, 2) * C(7, 3) * 0.5 * 3.
        let value = combination_rate(0.5, 3.0, &[10, 7], &[(0, 2), (1, 3)]);
        assert_eq!(value, 0.5 * 3.0 * 45.0 * 35.0);
    }

    #[test]
    fn high_order_and_small_coefficient_do_not_overflow_intermediates() {
        // Ported from the network engine (tests/network.rs): C(2000, 1000)
        // overflows f64, 1e-300 times it does not.
        let value = combination_rate(1e-300, 1.0, &[2000], &[(0, 1000)]);
        let close = (value - 2.048151626989489e300).abs() <= 1e-9 * 2.048151626989489e300;
        assert!(close, "{value}");
        assert_eq!(combination_rate(1e-300, 1.0, &[1000], &[(0, 1000)]), 1e-300);
        // The factor joins the scaled product: neither 1e-300 * 1e-300 nor
        // C(2000, 1000) is representable, their product is.
        let value = combination_rate(1e-300, 1e-300, &[2000], &[(0, 1000)]);
        assert!((value / 2.048151626989489e0 - 1.0).abs() < 1e-9, "{value}");
    }

    #[test]
    fn true_overflow_is_infinite_and_bounded_in_work() {
        // Ported from the network engine: count 2^31 - 1 of 2^53 molecules
        // with the smallest subnormal rate overflows after a few thousand
        // ratios, not 2^31.
        let n = 1u64 << 53;
        let (value, ratios) = scaled(f64::from_bits(1), 1.0, &[n], &[(0, 2_147_483_647)]);
        assert_eq!(value, f64::INFINITY);
        assert!(ratios < 4000, "{ratios} ratios");
        // With the smallest factor as well.
        let (value, ratios) = scaled(
            f64::from_bits(1),
            f64::from_bits(1),
            &[n],
            &[(0, 2_147_483_647)],
        );
        assert_eq!(value, f64::INFINITY);
        assert!(ratios < 4000, "{ratios} ratios");
        // Just below overflow and just above underflow.
        assert_eq!(combination_rate(f64::MAX, 1.0, &[5], &[(0, 5)]), f64::MAX);
        assert_eq!(
            combination_rate(f64::MAX, 2.0, &[5], &[(0, 5)]),
            f64::INFINITY
        );
        assert_eq!(
            combination_rate(f64::from_bits(1), 0.5, &[3], &[(0, 3)]),
            0.0
        );
        assert_eq!(
            combination_rate(f64::from_bits(1), 3.0, &[3], &[(0, 3)]).to_bits(),
            3
        );
    }

    #[test]
    fn invalid_rates_carry_through_and_zero_stays_zero() {
        assert!(combination_rate(f64::NAN, 1.0, &[5], &[(0, 3)]).is_nan());
        assert!(combination_rate(1.0, f64::NAN, &[5], &[(0, 3)]).is_nan());
        assert_eq!(
            combination_rate(f64::INFINITY, 1.0, &[5], &[(0, 3)]),
            f64::INFINITY
        );
        assert!(combination_rate(-1.0, 1.0, &[5], &[(0, 3)]) < 0.0);
        assert_eq!(combination_rate(0.0, 1.0, &[1 << 53], &[(0, 1 << 30)]), 0.0);
        assert_eq!(combination_rate(1.0, 0.0, &[5], &[(0, 3)]), 0.0);
    }

    #[test]
    fn the_network_engine_gives_the_same_bits() {
        // The network engine's general mass-action path (orders of three and
        // more) runs the same scaled product; with factor 1 the bits agree.
        for (reactants, counts) in [
            (vec![("A", 3)], vec![40_i64]),
            (vec![("A", 2), ("B", 1)], vec![30, 20]),
            (vec![("A", 1), ("B", 1), ("C", 1)], vec![7, 11, 13]),
            (vec![("A", 1000)], vec![2000]),
            (vec![("A", 5), ("B", 4)], vec![1 << 40, 1 << 30]),
        ] {
            for c in [0.002, 1.0, 1e-300, 3.7e10] {
                let species: Vec<_> = ["A", "B", "C"]
                    .iter()
                    .map(|s| serde_json::json!({"name": s, "initial_count": 0}))
                    .collect();
                let terms: serde_json::Map<_, _> = reactants
                    .iter()
                    .map(|(s, k)| ((*s).to_string(), serde_json::json!(k)))
                    .collect();
                let model = serde_json::json!({
                    "schema": "crn-craft.model/0.1", "name": "m", "convention": "stochastic_counts",
                    "omega": 1.0, "parameters": {}, "species": species,
                    "reactions": [{"name": "r", "reactants": terms, "products": {},
                                   "rate": {"kind": "mass_action", "coefficient": c}}]
                });
                let network = crate::CompiledNetwork::compile(
                    &crate::parse_model(&serde_json::to_vec(&model).unwrap()).unwrap(),
                )
                .unwrap();
                let mut state = vec![0_i64; 3];
                state[..counts.len()].copy_from_slice(&counts);
                let old = network.propensity(0, &state);
                let list: Vec<(usize, u64)> = reactants
                    .iter()
                    .enumerate()
                    .map(|(i, (_, k))| (i, *k as u64))
                    .collect();
                #[allow(clippy::cast_sign_loss)]
                let counts: Vec<u64> = state.iter().map(|&n| n as u64).collect();
                let new = combination_rate(c, 1.0, &counts, &list);
                match old {
                    Ok(old) => assert_eq!(old.to_bits(), new.to_bits(), "{reactants:?} c {c}"),
                    // The network engine refuses an overflow; here it is
                    // infinite, which the channel refuses (`CRN-R030`).
                    Err(d) => {
                        assert_eq!(d.message, "mass-action propensity overflows f64");
                        assert_eq!(new, f64::INFINITY, "{reactants:?} c {c}");
                    }
                }
            }
        }
    }
}
