//! Random-word mappings and the Poisson sampler of the unified engine.
//!
//! The engine draws 64-bit words from one ChaCha8 stream per replicate
//! ([`crate::rules::stream`]). [`uniform_closed_open`] maps a word to a
//! 53-bit uniform in [0, 1) and `poisson_with` draws a Poisson count from
//! the caller's open uniforms (tau leaping and binomial steps). Release
//! 0.7.0 removed the network engine's Direct sampler
//! (`crn-craft.chacha8-direct/0.1`), its streams and its `fast` policy.
//!
//! ```
//! use crn_craft::sampling::{UNIFORM_STEP, poisson_with, uniform_closed_open};
//!
//! assert_eq!(uniform_closed_open(0), 0.0);
//! assert_eq!(uniform_closed_open(u64::MAX), 1.0 - UNIFORM_STEP);
//! // The multiplication method: 0.5 · 0.5 is below e^-1, so one event.
//! assert_eq!(poisson_with(1.0, || Ok::<_, ()>(0.5)), Ok(1));
//! ```

/// 2⁻⁵³: the spacing of the uniform lattice of [`uniform_closed_open`].
pub const UNIFORM_STEP: f64 = 1.0 / 9_007_199_254_740_992.0;

/// Maps a random word to [0, 1): its high 53 bits times 2⁻⁵³.
pub fn uniform_closed_open(word: u64) -> f64 {
    (word >> 11) as f64 * UNIFORM_STEP
}

/// A ChaCha8 stream keyed by SHA-256 of the domain `crn-craft.chacha8/0.1`,
/// one NUL byte and the root seed, with the replicate as stream id: the
/// network engine's stream of release 0.6.0, kept as a fixed random source
/// for the unit tests of the selectors and samplers.
#[cfg(test)]
pub(crate) fn stream(root_seed: u64, replicate: u64) -> rand_chacha::ChaCha8Rng {
    keyed_stream(b"crn-craft.chacha8/0.1\0", root_seed, replicate)
}

#[cfg(test)]
fn keyed_stream(domain: &[u8], root_seed: u64, replicate: u64) -> rand_chacha::ChaCha8Rng {
    use rand_core::SeedableRng;
    use sha2::{Digest, Sha256};
    let mut digest = Sha256::new();
    digest.update(domain);
    digest.update(root_seed.to_le_bytes());
    let mut rng = rand_chacha::ChaCha8Rng::from_seed(digest.finalize().into());
    rng.set_stream(replicate);
    rng
}

/// Leap means above this use transformed rejection; at or below it the
/// multiplication method draws one Poisson variate directly.
const PTRS_MEAN: f64 = 16.0;

/// Shared Poisson kernel. Callers validate the mean and supply open uniforms
/// from their own budgeted sampler.
#[doc(hidden)]
pub fn poisson_with<E>(mean: f64, mut open: impl FnMut() -> Result<f64, E>) -> Result<u64, E> {
    if mean == 0.0 {
        return Ok(0);
    }
    if mean > PTRS_MEAN {
        return ptrs(mean, open);
    }
    let threshold = libm::exp(-mean);
    let mut count = 0;
    let mut product = open()?;
    while product > threshold {
        count += 1;
        product *= open()?;
    }
    Ok(count)
}

/// Hörmann's PTRS for means of at least 10, with his published constants.
/// Uniforms are drawn from (0, 1), so `us` is positive and every logarithm is
/// finite. Each pair of uniforms is accepted by the squeeze, rejected, or
/// accepted by the exact log-density comparison.
fn ptrs<E>(mean: f64, mut open: impl FnMut() -> Result<f64, E>) -> Result<u64, E> {
    let root = mean.sqrt();
    let b = 0.931 + 2.53 * root;
    let a = -0.059 + 0.02483 * b;
    let inverse_alpha = 1.1239 + 1.1328 / (b - 3.4);
    let v_r = 0.9277 - 3.6224 / (b - 2.0);
    let log_mean = mean.ln();
    loop {
        let u = open()? - 0.5;
        let v = open()?;
        let us = 0.5 - u.abs();
        let k = ((2.0 * a / us + b) * u + mean + 0.43).floor();
        if us >= 0.07 && v <= v_r {
            return Ok(k as u64);
        }
        if k < 0.0 || (us < 0.013 && v > us) {
            continue;
        }
        if v.ln() + inverse_alpha.ln() - (a / (us * us) + b).ln()
            <= -mean + k * log_mean - libm::lgamma(k + 1.0)
        {
            return Ok(k as u64);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::RngCore;

    /// The stream of 0.6.0's tau leaping, the source of these tests' words.
    fn tau_stream(root_seed: u64, replicate: u64) -> rand_chacha::ChaCha8Rng {
        keyed_stream(b"crn-craft.chacha8-tau/0.1\0", root_seed, replicate)
    }

    /// A uniform in (0, 1): [`uniform_closed_open`], redrawing a zero.
    fn uniform_open(rng: &mut impl RngCore) -> f64 {
        loop {
            let value = uniform_closed_open(rng.next_u64());
            if value > 0.0 {
                return value;
            }
        }
    }

    /// Open uniforms of a tau stream, never failing.
    fn open(rng: &mut impl RngCore) -> impl FnMut() -> Result<f64, ()> + '_ {
        move || Ok(uniform_open(rng))
    }

    /// The law test of `tau.rs` (P19, part a, declared in
    /// docs/campaign/performance-p11-p20.md) on the kernel itself: 200,000
    /// draws at each of five means from fixed streams, a two-sided DKW band
    /// at family alpha 0.001 (Bonferroni over the means) against the exact
    /// Poisson CDF, and every sample mean within five standard errors.
    #[test]
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    fn poisson_with_follows_the_poisson_law_on_both_sides_of_the_switch() {
        let draws = 200_000u32;
        let alpha = 0.001 / 5.0;
        let band = ((2.0f64 / alpha).ln() / (2.0 * f64::from(draws))).sqrt();
        for (index, mean) in [16.5f64, 30.0, 100.0, 1000.0, 1e5].into_iter().enumerate() {
            let mut rng = tau_stream(20261101, index as u64);
            let mut counts = std::collections::BTreeMap::<u64, u32>::new();
            let mut sum = 0.0;
            for _ in 0..draws {
                let k = poisson_with(mean, open(&mut rng)).unwrap();
                *counts.entry(k).or_default() += 1;
                sum += k as f64;
            }
            let sample_mean = sum / f64::from(draws);
            let standard_error = (mean / f64::from(draws)).sqrt();
            assert!(
                (sample_mean - mean).abs() <= 5.0 * standard_error,
                "mean {mean}: sample mean {sample_mean}"
            );
            let spread = 40.0 * mean.sqrt();
            let low = (mean - spread).max(0.0).floor() as u64;
            let high = (mean + spread).ceil() as u64;
            let (first, last) = (
                *counts.keys().next().unwrap(),
                *counts.keys().next_back().unwrap(),
            );
            assert!(
                first >= low && last <= high,
                "draws outside 40 sd at {mean}"
            );
            let (mut cdf, mut empirical, mut worst) = (0.0f64, 0u64, 0.0f64);
            for k in low..=high {
                let x = k as f64;
                cdf += (-mean + x * mean.ln() - libm::lgamma(x + 1.0)).exp();
                empirical += u64::from(counts.get(&k).copied().unwrap_or(0));
                worst = worst.max((empirical as f64 / f64::from(draws) - cdf).abs());
            }
            assert!(worst <= band, "mean {mean}: DKW distance {worst} > {band}");
        }
    }

    #[test]
    fn poisson_with_keeps_the_multiplication_method_for_small_means() {
        assert_eq!(poisson_with(0.0, || Err::<f64, ()>(())), Ok(0));
        for mean in [0.5, 3.0, 16.0] {
            let (mut a, mut b) = (tau_stream(7, 1), tau_stream(7, 1));
            for _ in 0..1000 {
                let k = poisson_with(mean, open(&mut a)).unwrap();
                let threshold = libm::exp(-mean);
                let mut expected = 0;
                let mut product = uniform_open(&mut b);
                while product > threshold {
                    expected += 1;
                    product *= uniform_open(&mut b);
                }
                assert_eq!(k, expected);
            }
            assert_eq!(a.next_u64(), b.next_u64(), "same words consumed at {mean}");
        }
    }

    #[test]
    fn poisson_with_draws_a_few_uniforms_for_large_means_and_passes_errors_on() {
        for mean in [17.0, 1e3, 1e5, 1e7] {
            let mut rng = tau_stream(3, 0);
            let mut uniforms = 0_u64;
            for _ in 0..10_000 {
                poisson_with(mean, || {
                    uniforms += 1;
                    Ok::<_, ()>(uniform_open(&mut rng))
                })
                .unwrap();
            }
            assert!(uniforms <= 30_000, "{uniforms} uniforms at {mean}");
        }
        // A failing sampler (an exhausted budget) ends the draw with its error.
        assert_eq!(
            poisson_with(3.0, || Err::<f64, &str>("budget")),
            Err("budget")
        );
        assert_eq!(
            poisson_with(300.0, || Err::<f64, &str>("budget")),
            Err("budget")
        );
    }
}
