//! Quantile statistics of `crn-craft.result/1.0`.
//!
//! The engine collects the reported values of every completed replicate in
//! one buffer, replicate × report time × observable, written at the
//! replicate's index, so the buffer does not depend on the thread count.
//! After the ensemble each (report time, observable) column is sorted with
//! [`f64::total_cmp`] and every requested probability is evaluated with
//! [`linear`].

use super::Quantile;

/// NumPy's default `linear` quantile (Hyndman and Fan type 7) of `sorted`,
/// which must be non-empty and sorted ascending.
///
/// With n values, h = (n − 1)·p, j = ⌊h⌋ and g = h − j, the quantile is
/// x\[j\] + g·(x\[j+1\] − x\[j\]). It is evaluated as `numpy.quantile` does:
/// as written for g < 0.5, as x\[j+1\] − (1 − g)·(x\[j+1\] − x\[j\]) for
/// g ≥ 0.5, and as x\[n−1\] when h ≥ n − 1, so the value equals
/// `numpy.quantile(x, p)` bit for bit. A NaN among the values makes the
/// quantile NaN, as in NumPy.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
pub(crate) fn linear(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    // total_cmp puts a NaN with the sign bit first and any other NaN last.
    if sorted[0].is_nan() || sorted[n - 1].is_nan() {
        return f64::NAN;
    }
    let last = (n - 1) as f64;
    let h = last * p;
    // NumPy takes index -1 (the last value) for both neighbors at or past
    // the end, with the weight h - (-1).
    let (below, above, g) = if h >= last {
        (n - 1, n - 1, h + 1.0)
    } else {
        let j = h.floor();
        (j as usize, j as usize + 1, h - j)
    };
    let (a, b) = (sorted[below], sorted[above]);
    let difference = b - a;
    if g >= 0.5 {
        b - difference * (1.0 - g)
    } else {
        a + difference * g
    }
}

/// The quantiles of `samples` replicates' values for each probability, in
/// the given order. `buffer` holds `samples × rows × columns` values,
/// replicate-major. The values of each quantile are shaped like the mean:
/// report time × observable, and an empty list when there are no
/// observables.
pub(crate) fn statistics(
    buffer: &[f64],
    samples: usize,
    rows: usize,
    columns: usize,
    probabilities: &[f64],
) -> Vec<Quantile> {
    let cells = rows * columns;
    debug_assert_eq!(buffer.len(), samples * cells);
    let mut values = vec![vec![0.0; cells]; probabilities.len()];
    let mut column = Vec::with_capacity(samples);
    for cell in 0..cells {
        column.clear();
        column.extend((0..samples).map(|i| buffer[i * cells + cell]));
        column.sort_unstable_by(f64::total_cmp);
        for (out, &p) in values.iter_mut().zip(probabilities) {
            out[cell] = linear(&column, p);
        }
    }
    probabilities
        .iter()
        .zip(values)
        .map(|(&p, flat)| Quantile {
            p,
            values: if columns == 0 {
                Vec::new()
            } else {
                flat.chunks(columns).map(<[f64]>::to_vec).collect()
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sorted(mut values: Vec<f64>) -> Vec<f64> {
        values.sort_by(f64::total_cmp);
        values
    }

    #[test]
    fn linear_matches_numpy_on_small_samples() {
        // numpy.quantile([1, 2, 3, 4], p) for these p, printed with repr.
        let x = [1.0, 2.0, 3.0, 4.0];
        for (p, expected) in [
            (0.0, 1.0),
            (0.1, 1.3),
            (0.25, 1.75),
            (0.5, 2.5),
            (0.75, 3.25),
            (0.9, 3.7),
            (1.0, 4.0),
        ] {
            assert_eq!(linear(&x, p), expected, "p = {p}");
        }
        // Both branches of the evaluation: numpy.quantile([0, 1], p).
        assert_eq!(linear(&[0.0, 1.0], 0.3), 0.3);
        assert_eq!(linear(&[0.0, 1.0], 0.7), 0.7);
        assert_eq!(linear(&[7.0], 0.025), 7.0);
        assert_eq!(linear(&[2.0, 2.0, 2.0], 0.5), 2.0);
    }

    #[test]
    fn the_upper_branch_gives_numpys_bits() {
        // numpy.quantile([16, 25], 0.925) == 24.325, while the formula as
        // written gives 16 + 0.925 * 9 == 24.325000000000003.
        let x = [16.0, 25.0];
        assert_eq!(linear(&x, 0.925), 24.325);
        assert_eq!(16.0 + 0.925 * 9.0, 24.325_000_000_000_003);
        // numpy.quantile([0, 0.1], 0.7) == 0.07; the formula gives
        // 0.06999999999999999.
        assert_eq!(linear(&[0.0, 0.1], 0.7), 0.07);
    }

    #[test]
    fn non_finite_values_follow_numpy() {
        assert!(linear(&sorted(vec![1.0, f64::NAN, 3.0]), 0.1).is_nan());
        assert!(linear(&sorted(vec![1.0, -f64::NAN, 3.0]), 0.9).is_nan());
        // numpy.quantile: [1, inf] at 0.25 is inf and at 0.5 nan
        // (inf - inf on the upper branch); [-inf, 1] at 0.25 nan and at
        // 0.75 -inf; [inf, inf] and [-inf, inf] at 0.5 nan.
        assert_eq!(linear(&[1.0, f64::INFINITY], 0.25), f64::INFINITY);
        assert!(linear(&[1.0, f64::INFINITY], 0.5).is_nan());
        assert!(linear(&[f64::NEG_INFINITY, 1.0], 0.25).is_nan());
        assert_eq!(linear(&[f64::NEG_INFINITY, 1.0], 0.75), f64::NEG_INFINITY);
        assert!(linear(&[f64::INFINITY, f64::INFINITY], 0.5).is_nan());
        assert!(linear(&[f64::NEG_INFINITY, f64::INFINITY], 0.5).is_nan());
    }

    #[test]
    fn statistics_reads_columns_across_replicates() {
        // Three replicates, two report times, two observables.
        let buffer = [
            1.0, 10.0, 100.0, 1000.0, //
            3.0, 30.0, 300.0, 3000.0, //
            2.0, 20.0, 200.0, 2000.0,
        ];
        let quantiles = statistics(&buffer, 3, 2, 2, &[0.5, 0.75]);
        assert_eq!(quantiles[0].p, 0.5);
        assert_eq!(quantiles[0].values, [vec![2.0, 20.0], vec![200.0, 2000.0]]);
        assert_eq!(quantiles[1].values, [vec![2.5, 25.0], vec![250.0, 2500.0]]);
        let empty = statistics(&[], 3, 2, 0, &[0.5]);
        assert!(empty[0].values.is_empty());
    }
}
