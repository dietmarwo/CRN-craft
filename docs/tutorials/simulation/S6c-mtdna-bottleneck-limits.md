# Choices, checks and limits

## Reconstruction choices

The paper leaves these points open, or the example departs from it. The
measured effect comes from the stored review files.

| Choice | Paper | Example | Measured effect |
|---|---|---|---|
| Ages of the cell-count rows at 7, 14 and 21 h | printed as 0.29, 0.58, 0.88 dpc | 7/24, 14/24, 21/24 dpc, after the division | Large: with the printed values about half of the accepted states leave the threshold (`review-age-shift.csv`) |
| Observation exactly at a division | not stated | sees the daughter | part of the row above |
| Degradation in phase 6 | replication fixed at 0; degradation under the [0, 1] prior | fixed at 0 | Small on Figure 3 and the groups (`review-final-degradation.csv`) |
| Log base in the distance | not stated | natural | not measured; base 10 would weight copy number 5.3 times less |
| Variance term in the distance | written `V(h)` in the appendix formula, normalized variance in the text and Table 1 | `V'(h)` | not measured |
| Heteroplasmy variance in the distance | first-order expansion | first-order expansion of exact count moments | below 3.4% relative error (`review-vprime.csv`) |
| Sampler | Metropolis chain, small normal steps | 4 ensembles of 16 differential evolution walkers | same target; ensembles disagree in parameter tails |
| First division | not stated | at 7 h | not measured |
| Extra rejection | none | first-order `V'` above 1 at a data age | not measured |
| Figure 3 reference | the published figure | central curve and two envelopes read by eye at 6 ages | ±3 pixels on the center only |
| Figure 4 comparison | half of the groups fitted, transformed single cells of the other half against a 95% range | no group fitted; group `V'` against a sample-size range | a new validation, not Figure 4D |

## What the check guards

`validation/check_mtdna_bottleneck.py` runs in about a second without
Numba and without a stored run. It checks:

- the binomial partition law;
- the founder-family law of the descendant sampler against the CRN-craft
  finite state projection;
- the native lineage driver;
- the mean recurrence across all six phases and 36 divisions, and the
  first-order variance added by divisions;
- the reference files: 74 Table 1 rows and their count per study, the four
  division-timed ages, 339 cells in 15 groups, and group means and `V'`
  against the spreadsheet summary;
- the structure of the stored results and their distances;
- every table, list and quoted number on these pages against the stored
  result files.

`--write` regenerates the tables after a new run.

## What was not verified

- The 74 Table 1 values, sample sizes and studies were compared with the
  first published version of the paper. Table 2 of the final version and
  the Figure 4 source spreadsheet could not be opened during the review.
  The HB cell values were checked only against the summary rows of the
  same spreadsheet.
- The Figure 3 image could not be viewed during the review. The
  digitized values and their calibration are unchecked.
- The paper does not state the width of its two Figure 3 envelopes or the
  number of accepted parameter sets behind them.
- The cause of the 100 dpc copy-number discrepancy is not found. Freeing
  the phase-6 degradation and the first-order variance are ruled out.
- The draft and full sizes were not timed without Numba.

## Data and license

The paper, its tables and its source data are published under the
Creative Commons Attribution 4.0 license. The reference CSV files are
derived from them and keep that license; cite the paper when reusing
them. `fetch.py` downloads and converts the sources again. No code of the
authors was used.

## Exercises

1. Run `study.py --size draft` and compare its effective sizes with the
   stored full run. Why does the effective size per sweep drop for the
   longer run?
2. Change the three early ages back to 0.29, 0.58 and 0.88 dpc in
   `reference-table1.csv` and rerun the draft size. Which part of Figure 3
   moves?
3. Add a copy-number constraint at 100 dpc (for example, the published
   central value with a factor-two error) and see whether the late
   `V'` values can still be fitted.
4. Replace the group-sized band with exact lineages from `review.py` and
   count the groups again.
