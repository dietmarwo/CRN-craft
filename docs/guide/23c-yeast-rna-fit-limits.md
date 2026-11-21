# 23c. Choices, checks and limits

[Chapter 23](23-yeast-rna-fit.md) ·
[Data and model](23a-yeast-rna-fit-model.md) ·
[Fits and predictions](23b-yeast-rna-fit-results.md) · Choices, checks
and limits

## Reconstruction choices

The paper and the authors' files leave these points open, or the example
decided differently from the authors' code. "Effect" is measured where a
number is given (`results/review-data-choices.csv` scores the saved
vectors with both variants); otherwise it was not measured.

| Choice | Paper / authors' code | This example | Effect |
| --- | --- | --- | --- |
| Count array | Paper: 39,943 cells at 0.4 M; code: mature array with 40,132 | mature array | 189 cells, two replicates; not measured |
| 0.4 M replicates fitted | code: 201301004 and 20160509 only | all three | archived vector 2.01804 on the authors' cells; the fit here moves toward 20131022 (see 23b) |
| Total RNA | code: nuclear + cytoplasmic | recorded total (235 cells differ) | +0.0002 nats per cell |
| Time zero | code: excluded from both objectives; paper silent | excluded | not measured |
| Mean likelihood groups | code: replicates pooled per time | one group per replicate and time | not measured |
| Mean variance | paper Eq. 23: `s²/N` | `s²/N + 1/(12N)` | not measured |
| Histogram tail | no tail bin; cap large enough | one tail bin above the cap | leakage ≤ 6.6e-6 |
| ON fraction | "more than three" and "three or more" | `X > 3` | figure only |
| Search box | authors: MCMC and several fit stages | one box around the archived histogram fit | see 23b |
| Held-out replicates | not in the paper | 20140429 and 20160509, refitted at draft size | two of five |
| Spatial fit | 15 coefficients fitted | 2 varied on a 3 × 3 grid; result is the archived vector | not a refit |
| Site prediction | elongation 63 nt/s, `L` from Table S1, active above 2 | same; `L = 1630 nt`; delay from the archived vector | agrees with the authors' curve to 0.0024 |
| Site reference | Figure 4D–E | the numerical series in the authors' figure files | differ from the raw site array by up to 0.27 |

The bioRxiv supplement (Table S3) prints the nonspatial histogram fit
with `k_3 = 1.1e2 /s`; the archived vector at the pinned revision has
564 /s. The other printed entries agree with the archive to their two
digits. The journal version of the tables was not reachable.

## What the check verifies

`validation/check_yeast_rna_fit.py` runs in about a second, offline:

- the constant-signal law of CRN-craft's FSP against an explicit generator
  for the total and the joint model, the BDF2 solver against that law and
  the BDF2 means against the moment equations (skipped with a message if
  no C++ compiler is found; everything else runs without one);
- the hand-calculated histogram score;
- the shipped rows: cells per dose (25,511 and 40,132), usable cells after
  time zero (59,557), the 235 cells whose total differs from the
  compartment sum, and the site array totals (25,511 and 39,943);
- the archived spatial delay and the site detection threshold;
- that the result files have the expected rows and that the two tables of
  page 23b equal what `render.py` makes from them.

It does not refit anything. The scores themselves depend on the
platform's floating point and are not compared bit for bit.

## Limits

- The mean-versus-histogram comparison is made near the authors' answer;
  see [what it shows](23b-yeast-rna-fit-results.md#what-the-comparison-shows).
- The held-out scores are not clean: the box center was fitted by the
  authors to data that include the held-out replicates. Three of five
  replicates were never held out.
- The spatial model and the transcription-site prediction use the
  authors' archived vector. The site prediction is a reimplementation, not
  a forecast from a model fitted here.
- The fourfold budget check ran once per objective (one 8,000-evaluation
  restart); it did not improve on the best 2,000-evaluation restart.
- The full-size and spatial runs were not repeated by the review; the
  draft runs were, and reproduce.
- The measured site traces of the figure files and the raw site array
  disagree; which one the paper's Figure 4D shows was not resolved.

## Data and license

The authors' repository is licensed GPL-3.0; it has no separate license
for its data. `source-counts.csv`, `source-sites.csv`,
`reference-site-*.csv` and `source-parameters.json` are tables derived
from its MATLAB arrays and figure files. Whether these derived tables may
be shipped in this MIT-licensed repository has not been settled;
`fetch.py` rebuilds all of them from the pinned revision. The underlying
images are deposited at the BioImage Archive (S-BIAD1) and the Image Data
Resource (idr0047).

## Exercises

1. Refit with the authors' replicate selection (drop 20131022) and the
   total as nuclear plus cytoplasmic. Does the fit still move away from
   the archived vector?
2. Hold out each of the other three replicates in turn
   (`manifest.json`, `holdouts`). Does the histogram fit win on every one?
3. Fit the mean objective in the wide box of `review.py`. Where does it
   end, and how does that fit score on the histograms?
4. Change the ON threshold to `X ≥ 3` and redraw the summary figure.
