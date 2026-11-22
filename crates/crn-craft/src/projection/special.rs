//! Binomial and Poisson probabilities for the projection.
//!
//! The binomial and Poisson probability mass functions use the saddle-point
//! expansion of C. Loader, "Fast and accurate computation of binomial
//! probabilities" (2000): the Stirling-series remainder `stirlerr` and the
//! deviance term `bd0`, evaluated without cancellation. The Poisson upper
//! tail is summed term by term from its first omitted index (or, below the
//! mean, as one minus the lower sum, which is then at most about one half),
//! never as `1 - kept weight` near zero. `validation/check_projection_native.py`
//! compares both with SciPy on a declared grid (relative error at most
//! 1e-12 where the value is at least 1e-280).

/// `ln(sqrt(2 pi))`.
const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;
/// `ln(2 pi)`.
const LN_2PI: f64 = 1.837_877_066_409_345_5;

/// `stirlerr(n) = ln(n!) - (n + 1/2) ln(n) + n - ln(sqrt(2 pi))` for
/// `n = 0..=15`, computed in 60-digit decimal arithmetic and rounded.
const STIRLERR: [f64; 16] = [
    0.0,
    0.081_061_466_795_327_26,
    0.041_340_695_955_409_3,
    0.027_677_925_684_998_34,
    0.020_790_672_103_765_093,
    0.016_644_691_189_821_193,
    0.013_876_128_823_070_748,
    0.011_896_709_945_891_77,
    0.010_411_265_261_972_096,
    0.009_255_462_182_712_733,
    0.008_330_563_433_362_87,
    0.007_573_675_487_951_841,
    0.006_942_840_107_209_53,
    0.006_408_994_188_004_207,
    0.005_951_370_112_758_847_5,
    0.005_554_733_551_962_801,
];

/// The remainder of Stirling's formula for `ln(n!)`, `n >= 0`.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub(crate) fn stirlerr(n: f64) -> f64 {
    const S0: f64 = 1.0 / 12.0;
    const S1: f64 = 1.0 / 360.0;
    const S2: f64 = 1.0 / 1260.0;
    const S3: f64 = 1.0 / 1680.0;
    const S4: f64 = 1.0 / 1188.0;
    if n <= 15.0 {
        // Only integers reach this function.
        return STIRLERR[n as usize];
    }
    let nn = n * n;
    if n > 500.0 {
        (S0 - S1 / nn) / n
    } else if n > 80.0 {
        (S0 - (S1 - S2 / nn) / nn) / n
    } else if n > 35.0 {
        (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n
    } else {
        (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n
    }
}

/// The deviance term `x ln(x / m) + m - x`, by its series when `x` is near
/// `m`.
pub(crate) fn bd0(x: f64, m: f64) -> f64 {
    if (x - m).abs() < 0.1 * (x + m) {
        let v = (x - m) / (x + m);
        let mut s = (x - m) * v;
        let mut ej = 2.0 * x * v;
        let v2 = v * v;
        for j in 1..1000 {
            ej *= v2;
            let next = s + ej / f64::from(2 * j + 1);
            if next == s {
                return next;
            }
            s = next;
        }
        return s;
    }
    x * libm::log(x / m) + m - x
}

/// The binomial probability of `x` successes in `n` trials of success
/// probability `p` (failure probability `q = 1 - p`), for integers
/// `0 <= x <= n`.
pub(crate) fn binomial_pmf(x: f64, n: f64, p: f64, q: f64) -> f64 {
    if p == 0.0 {
        return if x == 0.0 { 1.0 } else { 0.0 };
    }
    if q == 0.0 {
        return if x == n { 1.0 } else { 0.0 };
    }
    if x == 0.0 {
        if n == 0.0 {
            return 1.0;
        }
        let lc = if p < 0.1 {
            -bd0(n, n * q) - n * p
        } else {
            n * libm::log(q)
        };
        return libm::exp(lc);
    }
    if x == n {
        let lc = if q < 0.1 {
            -bd0(n, n * p) - n * q
        } else {
            n * libm::log(p)
        };
        return libm::exp(lc);
    }
    let lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
    let lf = LN_2PI + libm::log(x) + libm::log1p(-x / n);
    libm::exp(lc - 0.5 * lf)
}

/// The Poisson probability of `x` (an integer `>= 0`) at mean `mu >= 0`.
pub(crate) fn poisson_pmf(x: f64, mu: f64) -> f64 {
    if mu == 0.0 {
        return if x == 0.0 { 1.0 } else { 0.0 };
    }
    if x == 0.0 {
        return libm::exp(-mu);
    }
    libm::exp(-stirlerr(x) - bd0(x, mu) - LN_SQRT_2PI - 0.5 * libm::log(x))
}

/// `P(K >= first)` for `K ~ Poisson(mu)`, `first >= 1`: the regularized
/// lower incomplete gamma function `P(first, mu)`.
#[allow(clippy::cast_precision_loss)]
pub(crate) fn poisson_tail(first: u64, mu: f64) -> f64 {
    if mu <= 0.0 {
        return 0.0;
    }
    let m = first as f64;
    if m > mu {
        // Terms decrease from `first` on; stop once the geometric bound of
        // the rest is below the rounding of the sum.
        let mut term = poisson_pmf(m, mu);
        let mut sum = 0.0;
        let mut j = m;
        while term > 0.0 {
            sum += term;
            j += 1.0;
            let ratio = mu / j;
            term *= ratio;
            if ratio < 1.0 && term / (1.0 - ratio) <= sum * (f64::EPSILON / 4.0) {
                sum += term;
                break;
            }
        }
        return sum;
    }
    // Below the mean the tail is at least about one half: subtract the lower
    // sum, accumulated from its largest term down.
    let mut term = poisson_pmf(m - 1.0, mu);
    let mut lower = 0.0;
    let mut j = m - 1.0;
    loop {
        lower += term;
        if j == 0.0 || term <= lower * (f64::EPSILON / 4.0) {
            break;
        }
        term *= j / mu;
        j -= 1.0;
    }
    (1.0 - lower).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact_binomial(x: u64, n: u64, p: f64) -> f64 {
        // Small n: exact integer binomial coefficient times powers.
        let mut c = 1.0f64;
        for j in 0..x {
            c = c * (n - j) as f64 / (j + 1) as f64;
        }
        c * p.powi(x as i32) * (1.0 - p).powi((n - x) as i32)
    }

    #[test]
    fn binomial_matches_small_exact_values() {
        for n in 0..30u64 {
            for &p in &[0.1, 0.25, 0.5, 0.7, 0.95] {
                let mut total = 0.0;
                for x in 0..=n {
                    let got = binomial_pmf(x as f64, n as f64, p, 1.0 - p);
                    let want = exact_binomial(x, n, p);
                    assert!(
                        (got - want).abs() <= 1e-13 * want.max(1e-300),
                        "n={n} x={x} p={p}: {got} vs {want}"
                    );
                    total += got;
                }
                assert!((total - 1.0).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn binomial_mass_sums_to_one_for_large_n() {
        let n = 100_000.0;
        let total: f64 = (0..=100_000)
            .map(|x| binomial_pmf(f64::from(x), n, 0.3, 0.7))
            .sum();
        assert!((total - 1.0).abs() < 1e-11, "{total}");
    }

    #[test]
    fn poisson_tail_matches_complement_and_known_values() {
        // P(K >= 1) = 1 - e^-mu.
        for &mu in &[1e-12, 1e-3, 0.5, 4.0, 16.0] {
            let want = -libm::expm1(-mu);
            let got = poisson_tail(1, mu);
            assert!((got - want).abs() <= 1e-14 * want, "mu={mu}");
        }
        // Tails decrease and the pmf sums to the complement.
        let mu = 7.5;
        let mut lower = 0.0;
        for k in 0..60u64 {
            lower += poisson_pmf(k as f64, mu);
            let tail = poisson_tail(k + 1, mu);
            assert!((lower + tail - 1.0).abs() < 1e-14, "k={k}");
        }
        assert_eq!(poisson_tail(3, 0.0), 0.0);
    }
}
