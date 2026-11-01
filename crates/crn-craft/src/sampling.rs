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
    let mut digest = Sha256::new();
    digest.update(b"crn-craft.chacha8/0.1\0");
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
    }
    // A 53-bit categorical uniform cannot resolve arbitrarily small channels.
    if rates
        .iter()
        .any(|&rate| rate > 0.0 && rate / total < UNIFORM_STEP)
    {
        return Err(numeric(
            "CRN-E061",
            "positive channel is below the 53-bit selection resolution",
        ));
    }
    Ok(total)
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
