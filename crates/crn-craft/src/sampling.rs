//! Pinned random-word streams and explicit floating-point mappings, version 0.1.
use crate::diagnostic::Diagnostic;
use rand_chacha::ChaCha8Rng;
use rand_core::{RngCore, SeedableRng};
use sha2::{Digest, Sha256};

pub const SAMPLER: &str = "crn-craft.chacha8-direct/0.1";
pub const UNIFORM_STEP: f64 = 1.0 / 9_007_199_254_740_992.0;

/// SHA256(domain || root_seed little endian) is the key; replicate is the
/// ChaCha 64-bit stream id. No hashing of thread count, model or save grid.
pub fn stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    keyed_stream(b"crn-craft.chacha8/0.1\0", root_seed, replicate)
}

/// Tau leaping's own key: a tau ensemble and an exact ensemble with the same
/// root seed do not share random words, so comparing them is not coupled.
pub fn tau_stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    keyed_stream(b"crn-craft.chacha8-tau/0.1\0", root_seed, replicate)
}

fn keyed_stream(domain: &[u8], root_seed: u64, replicate: u64) -> ChaCha8Rng {
    let mut digest = Sha256::new();
    digest.update(domain);
    digest.update(root_seed.to_le_bytes());
    let mut rng = ChaCha8Rng::from_seed(digest.finalize().into());
    rng.set_stream(replicate);
    rng
}

pub fn uniform_closed_open(word: u64) -> f64 {
    (word >> 11) as f64 * UNIFORM_STEP
}

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
/// changed. `prefix[i]` holds the running sum of `rates[..=i]`; entries before
/// `from` come from an earlier successful pass over the same values, and
/// `floor` never exceeds any positive rate those passes saw. Running sums
/// depend only on earlier rates, so the total, the recorded sums and the first
/// failure equal those of a full pass from index zero.
///
/// The common pass only adds. When every rate is non-negative, the total is
/// normal and finite, and the smallest positive rate is at least 2^-52 of the
/// total, no positive rate can be lost (that would need at most half an ulp,
/// i.e. 2^-53 of a running sum, and running sums never exceed the total) and
/// the 53-bit resolution test passes. Anything else repeats the checked pass
/// from index zero, which decides and reports exactly as before.
pub(crate) fn resume_total(
    rates: &[f64],
    prefix: &mut [f64],
    floor: &mut f64,
    from: usize,
) -> Result<f64, Diagnostic> {
    let mut total = if from == 0 { 0.0 } else { prefix[from - 1] };
    let mut low = if from == 0 { f64::INFINITY } else { *floor };
    let mut valid = true;
    for (sum, &rate) in prefix[from..].iter_mut().zip(&rates[from..]) {
        valid &= rate >= 0.0;
        total += rate;
        *sum = total;
        low = low.min(if rate > 0.0 { rate } else { f64::INFINITY });
    }
    if valid && total.is_normal() && low >= total * (1.0 / (1u64 << 52) as f64) {
        *floor = low;
        return Ok(total);
    }
    checked_resume(rates, prefix, floor)
}

#[cold]
fn checked_resume(rates: &[f64], prefix: &mut [f64], floor: &mut f64) -> Result<f64, Diagnostic> {
    let mut total = 0.0;
    let mut minimum = f64::INFINITY;
    for (sum, &rate) in prefix.iter_mut().zip(rates) {
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
    }
    *floor = minimum;
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

#[cfg(test)]
mod tests {
    use super::*;

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
        let (mut prefix, mut floor, mut full) = (Vec::new(), f64::INFINITY, Vec::new());
        let mut checked = 0;
        for _ in 0..1000 {
            let n = 1 + (rng.next_u64() % 30) as usize;
            let mut rates: Vec<f64> = (0..n).map(|_| draw(&mut rng)).collect();
            prefix.resize(n, 0.0);
            if resume_total(&rates, &mut prefix, &mut floor, 0).is_err() {
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
                let resumed = resume_total(&rates, &mut prefix, &mut floor, from);
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
    fn resumed_totals_match_full_passes_at_the_fast_path_boundaries() {
        let tiny = f64::MIN_POSITIVE;
        let cases: [&[f64]; 10] = [
            &[1.0, 2.0f64.powi(-52)],
            &[1.0, 2.0f64.powi(-53)],
            &[1.0, 2.0f64.powi(-53) * (1.0 + f64::EPSILON)],
            &[2.0f64.powi(-53), 1.0],
            &[1.0, 0.0, 2.0f64.powi(-54)],
            &[tiny, tiny * 2.0f64.powi(-40)],
            &[tiny * 2.0f64.powi(52), tiny],
            &[f64::MAX, f64::MAX],
            &[0.0, 0.0],
            &[f64::from_bits(1), f64::from_bits(1)],
        ];
        for rates in cases {
            let (mut resumed, mut full) = (vec![0.0; rates.len()], Vec::new());
            let mut floor = f64::INFINITY;
            let a = resume_total(rates, &mut resumed, &mut floor, 0);
            let b = total_with_prefix(rates, &mut full);
            match (a, b) {
                (Ok(a), Ok(b)) => {
                    assert_eq!(a.to_bits(), b.to_bits(), "{rates:?}");
                    assert_eq!(resumed, full);
                }
                (Err(a), Err(b)) => assert_eq!((a.code, a.message), (b.code, b.message)),
                (a, b) => panic!("{rates:?}: {a:?} versus {b:?}"),
            }
        }
    }

    #[test]
    fn the_resolution_check_uses_the_smallest_positive_rate() {
        // Not lost in the sum (more than half an ulp of 1.9), yet below 2^-53
        // of the total: the first positive rate would wrongly pass the check.
        let rates = [1.9, 1.5e-16];
        assert!(rates[0] + rates[1] > rates[0]);
        let mut prefix = Vec::new();
        let (mut sums, mut floor) = (vec![0.0; 2], f64::INFINITY);
        for error in [
            total_propensity(&rates).unwrap_err(),
            total_with_prefix(&rates, &mut prefix).unwrap_err(),
            resume_total(&rates, &mut sums, &mut floor, 0).unwrap_err(),
        ] {
            assert_eq!(error.code, "CRN-E061");
            assert!(error.message.contains("resolution"), "{}", error.message);
        }
    }
}
