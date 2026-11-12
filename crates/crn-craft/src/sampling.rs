//! Pinned random-word streams and explicit floating-point mappings, version 0.1.
//!
//! These are the building blocks of the checked Direct sampler
//! [`SAMPLER`] (`crn-craft.chacha8-direct/0.1`), specified in
//! [docs/direct-ssa.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/direct-ssa.md):
//! one ChaCha8 stream per replicate, 53-bit uniforms, a strict ordered
//! propensity total, an inverse-transform exponential wait and first-channel
//! selection. Numerical failures are diagnostics at `/simulation`: `CRN-E061`
//! for totals and selection, `CRN-E062` for waits and event times. The rules
//! engine reuses [`uniform_closed_open`] and, for its leaps, the Poisson
//! sampler [`poisson_with`].
//!
//! ```
//! use crn_craft::sampling::{select_channel, stream, total_propensity, uniform_closed_open,
//!     uniform_open, waiting_time};
//! use rand_core::RngCore;
//!
//! let rates = [1.0, 0.0, 3.0];
//! let total = total_propensity(&rates).unwrap();
//! let mut rng = stream(42, 0); // root seed 42, replicate 0
//! let wait = waiting_time(total, uniform_open(&mut rng)).unwrap();
//! assert!(wait > 0.0);
//! // A threshold of 0.5 * total = 2 lies in the third channel's interval.
//! assert_eq!(select_channel(&rates, total, 0.5).unwrap(), 2);
//! assert!(uniform_closed_open(rng.next_u64()) < 1.0);
//! ```
use crate::diagnostic::Diagnostic;
use rand_chacha::ChaCha8Rng;
use rand_core::{RngCore, SeedableRng};
use sha2::{Digest, Sha256};

/// Sampler identity recorded in Direct results (`simulate`, `simulate-cached`).
pub const SAMPLER: &str = "crn-craft.chacha8-direct/0.1";
/// 2⁻⁵³: the spacing of the uniform lattice, and the smallest share of the
/// total that a positive channel may have (`CRN-E061` below it).
pub const UNIFORM_STEP: f64 = 1.0 / 9_007_199_254_740_992.0;

/// SHA256(domain || root_seed little endian) is the key; replicate is the
/// ChaCha 64-bit stream id. No hashing of thread count, model or save grid.
/// The domain is the ASCII string `crn-craft.chacha8/0.1` and one NUL byte;
/// the stream starts at word position zero.
pub fn stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    keyed_stream(b"crn-craft.chacha8/0.1\0", root_seed, replicate)
}

/// Tau leaping's own key: a tau ensemble and an exact ensemble with the same
/// root seed do not share random words, so comparing them is not coupled.
/// The domain is `crn-craft.chacha8-tau/0.1` and one NUL byte.
pub fn tau_stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    keyed_stream(b"crn-craft.chacha8-tau/0.1\0", root_seed, replicate)
}

/// Xoshiro256++ (Blackman & Vigna, public domain), the generator Julia's
/// `Random.Xoshiro` uses. Only `fast` alternative-SSA runs draw from it.
#[derive(Clone, Debug)]
pub(crate) struct Xoshiro256pp {
    s: [u64; 4],
}
impl RngCore for Xoshiro256pp {
    #[inline]
    fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }
    fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        rand_core::impls::fill_bytes_via_next(self, dest);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

/// One Xoshiro256++ state per replicate, seeded from its own keyed ChaCha8
/// stream, so fast results stay independent of the worker count and share no
/// words with the exact methods at the same root seed.
pub(crate) fn fast_stream(root_seed: u64, replicate: u64) -> Xoshiro256pp {
    let mut seeder = keyed_stream(b"crn-craft.xoshiro256pp-fast/0.1\0", root_seed, replicate);
    let mut s = [0u64; 4];
    while s == [0; 4] {
        for word in &mut s {
            *word = seeder.next_u64();
        }
    }
    Xoshiro256pp { s }
}

/// Running sums with additions only, as Julia's Direct keeps them: no
/// per-channel loss or 53-bit resolution checks. A positive rate below half an
/// ulp of its running sum can then never be selected. Only `fast` runs use it.
pub(crate) fn fast_prefix(
    rates: &[f64],
    prefix: &mut [f64],
    from: usize,
) -> Result<f64, Diagnostic> {
    let mut total = if from == 0 { 0.0 } else { prefix[from - 1] };
    for (sum, &rate) in prefix[from..].iter_mut().zip(&rates[from..]) {
        total += rate;
        *sum = total;
    }
    if total.is_finite() {
        Ok(total)
    } else {
        Err(numeric("CRN-E061", "fast Direct sum overflows"))
    }
}

fn keyed_stream(domain: &[u8], root_seed: u64, replicate: u64) -> ChaCha8Rng {
    let mut digest = Sha256::new();
    digest.update(domain);
    digest.update(root_seed.to_le_bytes());
    let mut rng = ChaCha8Rng::from_seed(digest.finalize().into());
    rng.set_stream(replicate);
    rng
}

/// Maps a random word to [0, 1): its high 53 bits times 2⁻⁵³.
pub fn uniform_closed_open(word: u64) -> f64 {
    (word >> 11) as f64 * UNIFORM_STEP
}

/// Draws a uniform in (0, 1) with [`uniform_closed_open`], redrawing a zero.
pub fn uniform_open(rng: &mut impl RngCore) -> f64 {
    loop {
        let value = uniform_closed_open(rng.next_u64());
        if value > 0.0 {
            return value;
        }
    }
}

fn numeric(code: &str, message: &str) -> Diagnostic {
    Diagnostic::error(
        code,
        "/simulation",
        message,
        "Rescale rates/time or revise the model; this strict solver does not clip or silently drop probability.",
    )
}

/// Sums `rates` in order, strictly. An all-zero slice gives 0, which callers
/// treat as absorbing.
///
/// # Errors
///
/// `CRN-E061` for a negative or non-finite rate, an overflowing sum, a
/// positive rate lost in rounding, or a positive rate below
/// [`UNIFORM_STEP`] times the total.
pub fn total_propensity(rates: &[f64]) -> Result<f64, Diagnostic> {
    let mut total = 0.0;
    let mut minimum = f64::INFINITY;
    for &rate in rates {
        if !rate.is_finite() || rate < 0.0 {
            return Err(numeric("CRN-E061", "invalid propensity in Direct sum"));
        }
        let next = total + rate;
        if !next.is_finite() || (rate > 0.0 && next <= total) {
            return Err(numeric(
                "CRN-E061",
                "Direct sum overflows or loses a positive channel",
            ));
        }
        total = next;
        if rate > 0.0 {
            minimum = minimum.min(rate);
        }
    }
    resolution(minimum, total)?;
    Ok(total)
}

/// `total_propensity` that also records each running sum: `prefix[i]` is the
/// ordered sum of `rates[..=i]`, the same values `select_channel` accumulates.
pub(crate) fn total_with_prefix(rates: &[f64], prefix: &mut Vec<f64>) -> Result<f64, Diagnostic> {
    prefix.resize(rates.len(), 0.0);
    let mut total = 0.0;
    let mut minimum = f64::INFINITY;
    for (slot, &rate) in prefix.iter_mut().zip(rates) {
        if !rate.is_finite() || rate < 0.0 {
            return Err(numeric("CRN-E061", "invalid propensity in Direct sum"));
        }
        let next = total + rate;
        if !next.is_finite() || (rate > 0.0 && next <= total) {
            return Err(numeric(
                "CRN-E061",
                "Direct sum overflows or loses a positive channel",
            ));
        }
        total = next;
        *slot = total;
        if rate > 0.0 {
            minimum = minimum.min(rate);
        }
    }
    resolution(minimum, total)?;
    Ok(total)
}

/// The strict ordered total resumed at `from` after only `rates[from..]`
/// changed. `prefix[i]` and `minima[i]` hold the running sum and the running
/// minimum positive rate of `rates[..=i]`; entries before `from` come from an
/// earlier successful pass over the same values. Running sums depend only on
/// earlier rates and the minimum is exact, so the total, the recorded sums and
/// the first failure equal those of a full pass from index zero.
pub(crate) fn resume_total(
    rates: &[f64],
    prefix: &mut [f64],
    minima: &mut [f64],
    from: usize,
) -> Result<f64, Diagnostic> {
    let (mut total, mut minimum) = if from == 0 {
        (0.0, f64::INFINITY)
    } else {
        (prefix[from - 1], minima[from - 1])
    };
    for ((sum, least), &rate) in prefix[from..]
        .iter_mut()
        .zip(&mut minima[from..])
        .zip(&rates[from..])
    {
        if !rate.is_finite() || rate < 0.0 {
            return Err(numeric("CRN-E061", "invalid propensity in Direct sum"));
        }
        let next = total + rate;
        if !next.is_finite() || (rate > 0.0 && next <= total) {
            return Err(numeric(
                "CRN-E061",
                "Direct sum overflows or loses a positive channel",
            ));
        }
        total = next;
        *sum = total;
        if rate > 0.0 {
            minimum = minimum.min(rate);
        }
        *least = minimum;
    }
    resolution(minimum, total)?;
    Ok(total)
}

fn resolution(minimum: f64, total: f64) -> Result<(), Diagnostic> {
    // A 53-bit categorical uniform cannot resolve arbitrarily small channels.
    if minimum.is_finite() && minimum / total < UNIFORM_STEP {
        return Err(numeric(
            "CRN-E061",
            "positive channel is below the 53-bit selection resolution",
        ));
    }
    Ok(())
}

/// The channel `select_channel` returns for the same rates, found by binary
/// search over the running sums recorded by a successful `total_with_prefix`.
/// A zero rate leaves the running sum unchanged and the strict total rejects a
/// lost positive increment, so the first running sum above the threshold
/// always belongs to a positive channel.
pub(crate) fn select_from_prefix(
    prefix: &[f64],
    total: f64,
    uniform: f64,
) -> Result<usize, Diagnostic> {
    if !(0.0..1.0).contains(&uniform) || !total.is_finite() || total <= 0.0 {
        return Err(numeric("CRN-E061", "invalid categorical sampler arguments"));
    }
    let threshold = uniform * total;
    if threshold >= total {
        return Err(numeric(
            "CRN-E061",
            "selection threshold rounded to the total",
        ));
    }
    let index = prefix.partition_point(|&sum| sum <= threshold);
    if index < prefix.len() {
        Ok(index)
    } else {
        Err(numeric(
            "CRN-E061",
            "no channel contains the categorical threshold",
        ))
    }
}

/// Exponential wait `-ln(uniform) / total` in seconds, using pinned `libm`.
///
/// # Errors
///
/// `CRN-E062` unless `total` is finite and positive, `uniform` lies in
/// (0, 1) and the wait is finite and positive.
pub fn waiting_time(total: f64, uniform: f64) -> Result<f64, Diagnostic> {
    if !total.is_finite() || total <= 0.0 || !(0.0 < uniform && uniform < 1.0) {
        return Err(numeric("CRN-E062", "invalid exponential sampler arguments"));
    }
    let wait = -libm::log(uniform) / total;
    if !wait.is_finite() || wait <= 0.0 {
        return Err(numeric(
            "CRN-E062",
            "exponential wait is outside the positive finite range",
        ));
    }
    Ok(wait)
}

/// Advances the clock: `time + wait`.
///
/// # Errors
///
/// `CRN-E062` unless `time` is finite and non-negative, `wait` is finite and
/// positive and the sum is finite and strictly greater than `time`.
pub fn event_time(time: f64, wait: f64) -> Result<f64, Diagnostic> {
    let next = time + wait;
    if !time.is_finite()
        || time < 0.0
        || !wait.is_finite()
        || wait <= 0.0
        || !next.is_finite()
        || next <= time
    {
        return Err(numeric(
            "CRN-E062",
            "event time cannot advance in finite f64 arithmetic",
        ));
    }
    Ok(next)
}

/// Caller supplies the same ordered rate array whose total was checked.
///
/// Returns the first positive channel whose running sum exceeds
/// `uniform * total`; `uniform` must lie in [0, 1). There is no fallback to
/// the last channel.
///
/// # Errors
///
/// `CRN-E061` for invalid arguments, a threshold that rounds to the total,
/// or no channel above the threshold.
pub fn select_channel(rates: &[f64], total: f64, uniform: f64) -> Result<usize, Diagnostic> {
    if !(0.0..1.0).contains(&uniform) || !total.is_finite() || total <= 0.0 {
        return Err(numeric("CRN-E061", "invalid categorical sampler arguments"));
    }
    let threshold = uniform * total;
    if threshold >= total {
        return Err(numeric(
            "CRN-E061",
            "selection threshold rounded to the total",
        ));
    }
    let mut cumulative = 0.0;
    for (index, &rate) in rates.iter().enumerate() {
        cumulative += rate;
        if threshold < cumulative && rate > 0.0 {
            return Ok(index);
        }
    }
    Err(numeric(
        "CRN-E061",
        "no channel contains the categorical threshold",
    ))
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

    fn same(a: Result<usize, Diagnostic>, b: Result<usize, Diagnostic>) {
        match (a, b) {
            (Ok(x), Ok(y)) => assert_eq!(x, y),
            (Err(x), Err(y)) => assert_eq!((x.code, x.message), (y.code, y.message)),
            (x, y) => panic!("linear {x:?} differs from binary {y:?}"),
        }
    }

    #[test]
    fn binary_search_selects_the_linear_scan_channel_at_exact_boundaries() {
        // Dyadic rates summing to one: uniform*total lands exactly on running sums.
        let rates = [0.0, 0.25, 0.0, 0.0, 0.125, 0.5, 0.0, 0.125, 0.0];
        let mut prefix = Vec::new();
        let total = total_with_prefix(&rates, &mut prefix).unwrap();
        assert_eq!(total, total_propensity(&rates).unwrap());
        assert_eq!(total, 1.0);
        assert_eq!(prefix.len(), rates.len());
        for k in 0..1024 {
            let uniform = k as f64 / 1024.0;
            same(
                select_channel(&rates, total, uniform),
                select_from_prefix(&prefix, total, uniform),
            );
        }
        // Exactly at a running sum, the next positive channel is selected.
        assert_eq!(select_from_prefix(&prefix, total, 0.25).unwrap(), 4);
        assert_eq!(select_from_prefix(&prefix, total, 0.875).unwrap(), 7);
        same(
            select_channel(&rates, total, 1.0),
            select_from_prefix(&prefix, total, 1.0),
        );
    }

    #[test]
    fn binary_search_matches_the_linear_scan_on_random_rates() {
        let mut rng = stream(20261101, 0);
        let mut prefix = Vec::new();
        for _ in 0..4000 {
            let n = 1 + (rng.next_u64() % 40) as usize;
            let rates: Vec<f64> = (0..n)
                .map(|_| match rng.next_u64() % 4 {
                    0 => 0.0,
                    1 => 1e-30 * uniform_closed_open(rng.next_u64()),
                    _ => libm::exp(40.0 * (uniform_closed_open(rng.next_u64()) - 0.5)),
                })
                .collect();
            let linear = total_propensity(&rates);
            let binary = total_with_prefix(&rates, &mut prefix);
            match (&linear, &binary) {
                (Ok(a), Ok(b)) => assert_eq!(a.to_bits(), b.to_bits()),
                (Err(a), Err(b)) => {
                    assert_eq!((&a.code, &a.message), (&b.code, &b.message));
                    continue;
                }
                _ => panic!("total {linear:?} differs from {binary:?}"),
            }
            let total = linear.unwrap();
            if total == 0.0 {
                continue;
            }
            for _ in 0..64 {
                let uniform = uniform_closed_open(rng.next_u64());
                same(
                    select_channel(&rates, total, uniform),
                    select_from_prefix(&prefix, total, uniform),
                );
            }
        }
    }

    #[test]
    fn resumed_totals_equal_full_passes_after_suffix_changes() {
        let mut rng = stream(20261102, 0);
        let draw = |rng: &mut rand_chacha::ChaCha8Rng| match rng.next_u64() % 12 {
            0..=2 => 0.0,
            3 => 1e-300 * (1.0 + uniform_closed_open(rng.next_u64())),
            4 => 1e300 * (1.0 + uniform_closed_open(rng.next_u64())),
            _ => libm::exp(30.0 * (uniform_closed_open(rng.next_u64()) - 0.5)),
        };
        let (mut prefix, mut minima, mut full) = (Vec::new(), Vec::new(), Vec::new());
        let mut checked = 0;
        for _ in 0..1000 {
            let n = 1 + (rng.next_u64() % 30) as usize;
            let mut rates: Vec<f64> = (0..n).map(|_| draw(&mut rng)).collect();
            prefix.resize(n, 0.0);
            minima.resize(n, 0.0);
            if resume_total(&rates, &mut prefix, &mut minima, 0).is_err() {
                continue;
            }
            for _ in 0..20 {
                let from = (rng.next_u64() % n as u64) as usize;
                for rate in &mut rates[from..] {
                    if rng.next_u64() % 3 == 0 {
                        *rate = draw(&mut rng);
                    }
                }
                rates[from] = draw(&mut rng);
                let resumed = resume_total(&rates, &mut prefix, &mut minima, from);
                let expected = total_with_prefix(&rates, &mut full);
                match (&resumed, &expected) {
                    (Ok(a), Ok(b)) => {
                        assert_eq!(a.to_bits(), b.to_bits());
                        assert_eq!(prefix, full);
                        assert_eq!(a.to_bits(), total_propensity(&rates).unwrap().to_bits());
                        checked += 1;
                    }
                    (Err(a), Err(b)) => {
                        assert_eq!((&a.code, &a.message), (&b.code, &b.message));
                        break;
                    }
                    _ => panic!("resumed {resumed:?} differs from full {expected:?}"),
                }
            }
        }
        assert!(checked > 1000, "only {checked} successful resumed passes");
    }

    #[test]
    fn xoshiro256pp_matches_the_reference_sequence() {
        // Reference outputs for the state [1, 2, 3, 4] (Blackman & Vigna).
        let mut rng = Xoshiro256pp { s: [1, 2, 3, 4] };
        assert_eq!(rng.next_u64(), 41_943_041);
        assert_eq!(rng.next_u64(), 58_720_359);
        assert_eq!(rng.next_u64(), 3_588_806_011_781_223);
        assert_eq!(rng.next_u64(), 3_591_011_842_654_386);
        let (mut a, mut b) = (fast_stream(7, 3), fast_stream(7, 3));
        assert_eq!(a.next_u64(), b.next_u64());
        assert_ne!(fast_stream(7, 3).next_u64(), fast_stream(7, 4).next_u64());
        assert_ne!(fast_stream(7, 3).next_u64(), stream(7, 3).next_u64());
    }

    #[test]
    fn fast_prefix_sums_only_add() {
        let rates = [1.0, 0.0, 2.0, 1e-17, 3.0];
        let mut prefix = [0.0; 5];
        assert_eq!(fast_prefix(&rates, &mut prefix, 0).unwrap(), 6.0);
        assert_eq!(prefix, [1.0, 1.0, 3.0, 3.0, 6.0]);
        let changed = [1.0, 0.0, 5.0, 0.0, 3.0];
        assert_eq!(fast_prefix(&changed, &mut prefix, 2).unwrap(), 9.0);
        assert!(fast_prefix(&[f64::MAX, f64::MAX], &mut [0.0; 2], 0).is_err());
    }

    #[test]
    fn the_resolution_check_uses_the_smallest_positive_rate() {
        // Not lost in the sum (more than half an ulp of 1.9), yet below 2^-53
        // of the total: the first positive rate would wrongly pass the check.
        let rates = [1.9, 1.5e-16];
        assert!(rates[0] + rates[1] > rates[0]);
        let mut prefix = Vec::new();
        let (mut sums, mut minima) = (vec![0.0; 2], vec![0.0; 2]);
        for error in [
            total_propensity(&rates).unwrap_err(),
            total_with_prefix(&rates, &mut prefix).unwrap_err(),
            resume_total(&rates, &mut sums, &mut minima, 0).unwrap_err(),
        ] {
            assert_eq!(error.code, "CRN-E061");
            assert!(error.message.contains("resolution"), "{}", error.message);
        }
    }
}
