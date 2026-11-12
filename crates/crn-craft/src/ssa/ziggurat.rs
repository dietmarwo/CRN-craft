//! Exponential Ziggurat adapted from Julia Random; see THIRD_PARTY_NOTICES.md.
use super::ziggurat_tables::{FE, KE, R, WE};

/// Draw a unit exponential. The caller accounts for each word and can poll
/// for cancellation before a redraw following rejection (the flag is true).
#[inline]
pub fn exponential_from_words<E>(
    mut next_word: impl FnMut(bool) -> Result<u64, E>,
) -> Result<f64, E> {
    let mut retry = false;
    loop {
        let word = next_word(retry)? & 0x000f_ffff_ffff_ffff;
        // Exclude exact zero, as in the legacy inverse sampler.
        if word != 0 {
            let index = (word & 255) as usize;
            let x = word as f64 * WE[index];
            if word < KE[index] {
                return Ok(x);
            }
            let u = crate::sampling::uniform_closed_open(next_word(false)?);
            if index == 0 {
                return Ok(R - libm::log1p(-u));
            }
            let y = (FE[index - 1] - FE[index]) * u + FE[index];
            if y < libm::exp(-x) {
                return Ok(x);
            }
        }
        retry = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::RngCore;

    /// A unit exponential from a ChaCha8 stream, never failing.
    fn exponential(rng: &mut impl RngCore) -> f64 {
        exponential_from_words(|_| Ok::<_, ()>(rng.next_u64())).unwrap()
    }

    #[test]
    fn exponential_million_sample_cdf_gate() {
        let mut rng = crate::sampling::stream(4242, 778);
        let mut bins = [0u64; 91];
        let mut mean = 0.0;
        let mut second = 0.0;
        const N: u64 = 1_000_000;
        for _ in 0..N {
            let x = exponential(&mut rng);
            assert!(x > 0.0 && x.is_finite());
            bins[(x * 10.0).floor().min(90.0) as usize] += 1;
            mean += x / N as f64;
            second += x * x / N as f64;
        }
        let bound = (libm::log(2.0 * 100.0 / 0.001) / (2.0 * N as f64)).sqrt();
        let mut cumulative = 0;
        for (k, count) in bins.iter().take(90).enumerate() {
            cumulative += count;
            let expected = -libm::expm1(-((k + 1) as f64) / 10.0);
            assert!((cumulative as f64 / N as f64 - expected).abs() <= bound);
        }
        eprintln!("Ziggurat N={N} mean={mean} second_moment={second}");
    }

    /// Why a draw ended early: the caller's word budget or its stop poll.
    #[derive(Debug, PartialEq)]
    enum Bail {
        SamplingLimit,
        Stopped,
    }

    /// Feeds `words` cyclically within a budget of `limit` words; before a
    /// redraw (the flag is true) it polls `stopped`, as the engine's
    /// callers do.
    fn draw(words: &[u64], limit: u64, stopped: bool, used: &mut u64) -> Result<f64, Bail> {
        let mut position = 0;
        exponential_from_words(|retry| {
            if retry && stopped {
                return Err(Bail::Stopped);
            }
            if *used == limit {
                return Err(Bail::SamplingLimit);
            }
            *used += 1;
            let word = words[position % words.len()];
            position += 1;
            Ok(word)
        })
    }

    #[test]
    fn zero_tail_and_rejection_branches_respect_budgets_and_stops() {
        // Zero words are redrawn until the budget ends the draw.
        let mut used = 0;
        assert_eq!(draw(&[0], 4, false, &mut used), Err(Bail::SamplingLimit));
        assert_eq!(used, 4);
        // Layer zero with a uniform of 1/2 takes the tail: R + ln 2.
        let mut used = 0;
        let x = draw(&[0x000f_ffff_ffff_ff00, 1u64 << 63], 4, false, &mut used);
        assert_eq!(x, Ok(R + libm::log(2.0)));
        assert_eq!(used, 2);
        // Layer one always enters the slow branch; a high ordinate rejects.
        let mut used = 0;
        let words = [0x000f_ffff_ffff_ff01, u64::MAX];
        assert_eq!(draw(&words, 4, false, &mut used), Err(Bail::SamplingLimit));
        assert_eq!(used, 4);
        // The redraw after a rejection polls first.
        let mut used = 0;
        assert_eq!(draw(&words, 4, true, &mut used), Err(Bail::Stopped));
        assert_eq!(used, 2);
    }
}
