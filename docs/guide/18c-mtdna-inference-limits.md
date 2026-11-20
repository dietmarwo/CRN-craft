# Choices, checks and limits

Part of [chapter 18](18-mtdna-inference.md).

## Reconstruction choices

The paper describes its method in words, equations and figures. Where it
leaves something open, the reconstruction had to choose. All choices are
listed here, with their effect where it was measured.

| Choice | What the paper says | Reconstruction | Effect where measured |
| --- | --- | --- | --- |
| Days to Moran steps | `Δt` days are about `N*ln(2)/halfL*Δt` steps | `floor(days*1000*ln(2)/10)`: 74,860 steps for 1,080 days | Less than one step; not varied |
| Report ages | Three years, exported every 30 days | 0, 30, …, 1,080 days (37 ages); the companion preprint's copy of Figure 6B ends at 1,080 days | — |
| Forward law in the fit | Forward master equation of the Moran process (discrete steps) | Continuous-time (Poisson-clock) law, 5-day implicit split steps | At most 1.2e-5 from the discrete law at the true values, 3.4e-5 at the box corners |
| Neutral-count cap | Not stated | 192 (baseline), 320 (strong); escaped probability recorded, never renormalized | 5.6e-5 against the complete lattice at the baseline truth; escaped mass at most 3.2e-4 at any grid node |
| Spline surrogate | The fit evaluates the forward law | Cubic spline through 192 / 216 precomputed parameter points | At most 1.1e-4 over the whole box, 1.3e-6 at the fits |
| Search box | Only the bound `f ≤ 1` is stated (the companion preprint fits `s` in 0–2) | A box around the truth, the extent of the spline grid ([table](18b-mtdna-inference-results.md#the-protocol)) | Four refits in much wider boxes without the spline: same fits within 0.2% |
| Optimizer settings | fcmaes differential evolution | Population 32, 1,800 evaluations, one seed per cohort | Ten times the budget or other seeds: parameters change by at most 2e-6 relative |
| How synthetic cohorts are made | Methods: simulations of the Moran model. Results: output of Model-1 and Model-2 | The exact discrete Moran law | Model-1 data fitted once: `s` 0.104, `f` 0.757 (inconclusive) |
| Same or different cells across ages | Not stated | The same 10,000 cells at every age | A new set of cells at every age: intervals 50% to 65% narrower, CV of `s` 2.0% instead of 5.8% (baseline) |
| Percentile definition | 2.5th and 97.5th percentiles of 40 fits | Linear interpolation (NumPy default) | Other definitions widen the intervals toward the extremes, e.g. baseline `s` 0.088–0.114 (Hazen) or 0.083–0.115 (Weibull) instead of 0.092–0.112 |
| Takeover | Reaching 90% mutant copies | At least 900 of 1,000 copies of the founder's type | — |
| Start of a takeover | A single advantageous mutant in a wild-type cell | One `A` copy, 999 `W`, no neutral mutants, mutation switched off | — |
| Figure 5 bins | Histograms of accumulation times | 30-day bins named by their first day, among all founders that ever take over | The published bars agree within 0.0085 per bin; the other bin naming would differ by up to 0.19 ([takeover times](18b-mtdna-inference-results.md#from-a-selection-advantage-to-a-takeover-time)) |
| Model-1 time stepping | Stochastic simulation in one-hour steps | Chapter 17's order: degrade, replicate survivors, mutate wild-type daughters, relax ATP for one hour with the new counts | Half-hour steps: at most 0.28 percentage points, 2.3 standard errors |
| Model-1 constants | Equations without values in the text read by the review | 2014 PNAS Table 1: `f_ATP = 0.1`, `v2 = 0.01`, `v1 = 0.2`; `c` from the 1,000-copy steady state | Not varied |
| Model-1 neutral mutants | A fraction `fracAdv` of mutants is advantageous | Neutral mutants replicate like wild type | Not varied |
| Model-1 cell whose ATP reaches zero | Takeover coincides with ATP collapse | The cell keeps its counts from that hour on (chapter 17's choice) | 394 of 100,000 cells collapse by day 1,080 |

## How the reconstruction was checked

`validation/check_kk26.py` runs in CI, normally and under `python -O`, in
a few seconds:

- **Exact laws in a seven-copy cell**: the discrete forward solver against
  matrix powers of an independently built Equation 7 chain (within 1e-11);
  the strict thresholds; a neutral cap of one copy that must bound the
  exact tails from below and, with the escaped mass, from above;
  CRN-craft's FSP against the matrix exponential (within 1e-10); 16,000
  Static Direct cells and 16,000 cells of the geometric-skip sampler
  against the exact law, with the Dvoretzky–Kiefer–Wolfowitz bound for
  error probability 0.001 (0.0154).
- **Replay and takeover**: the same cohort with one and four threads is
  identical; a neutral and a selected ten-copy takeover chain match the
  gambler's-ruin probabilities.
- **Optimizer**: fcmaes differential evolution recovers a known target of
  a small artificial problem.
- **Saved evidence**: the solver hash and build flags in the manifest match
  `moran.cpp` and `model.py`; the five review files are complete; all 80
  fits are present and the summary
  statistics replay from them; all 44 interpolation checks and the 9
  forward checks are within 0.001; the Figure 6 differences use the
  digitized values; the modal bins and the order of the medians of Figure
  5 hold; the table of [the results page](18b-mtdna-inference-results.md#results)
  equals the table generated from the saved files.

The check does not rerun the study and compares no simulated cell, curve
or fit with a saved value: those depend on the platform. The forward grids
and the 80 fits change when a single bit of the grids changes (a 6e-15
change of the grid values moved the fits by up to 3e-8 relative). The
library is therefore built with `-ffp-contract=off` and without
`-march=native`, so that compilers do not fuse multiply-adds on CPUs that
support them. The grid solver then uses only IEEE arithmetic; the one
math-library result it depends on is `10**x` for the mutation probability
of each node, computed in Python. The cohort sampler uses the C math library (`log1p`) and
the ATP model `std::binomial_distribution`, whose results can differ
between standard libraries.

The saved results were produced on an Intel x86-64 host with GCC 13.3,
Python 3.13, NumPy 2.3.5, SciPy 1.18.1 and fcmaes 2.2.0. fcmaes 2.1.0
reproduces all 80 fits byte for byte on the same grids and cohorts; other
versions may give statistically equivalent, not identical, fits.

## Limits

- Each saved cohort follows the same cells through all ages. Real
  snapshots see different cells at every age and, with the same number of
  cells, pin the parameters down two to three times more precisely.
- The cohorts are simulated with the Moran model and fitted with the
  Moran model. The recovery shows that the method is consistent when its
  model is right; it says nothing about data from cells that do not follow
  it. The paper's misspecification experiments (wrong turnover or copy
  number, its Tables 4 and 5) and its 100,000-cell Table 1 are not
  reproduced.
- The fits search a box around the true values on a spline of
  precomputed curves. The review's wider refits cover four of 80 cohorts.
- The published intervals are wider in the baseline regime and the
  published means are shifted; the reason is not known (see
  [the results page](18b-mtdna-inference-results.md#how-to-read-the-recovery)).
- The authors' scripts and data (linked from the paper) could not be
  reached, so their optimizer settings, bounds and data generation are
  unknown.
- The published values and figures were checked against the preprint of
  the paper, not the journal article: the printed values of Tables 1 to 3
  all agree, as do the labels of Figure 5 and the readings of both panels
  of Figure 6. The published takeover distributions of Figure 5 lie 1.5 to
  3 days later than the reconstruction's; the cause is not known.
- The Model-1 reconstruction inherits chapter 17's choices and constants;
  the 2026 paper's own constants were not available to the review.
- Matching the published curves and intervals does not show that real
  cells follow the Moran model; the companion preprint finds that measured
  curves do not fully fit it.

## Exercises

1. Fit the saved first cohort (`study.py fit
   examples/kk26/results/reconstruction/snapshots-baseline.csv`) and then
   the same counts with bins 1 and 2 merged into bin 1 (drop the 5%
   threshold). Which parameter loses most precision?
2. Rerun `bootstrap --size quick` with `--cells 100000`. Predict, before
   you run it, how the spread of the four fits compares with 2,000 cells.
3. Compute the takeover table for a half-life of 20 or 30 days (the
   paper's sensitivity values) with `model.passage(Model(0, 0.1, 1,
   half_days=20))`. How do the success probability and the median scale?
4. Fit the 30,000 half-hour ATP cells of `fig6-curves.csv` the way
   `review.py atp-fit` fits the hourly ones and compare the two estimates.

See also: [Reproducing a published simulation](08-reproducibility.md#reproducing-a-published-simulation),
[Choices, checks and limits of chapter 17](17c-mtdna-mechanisms-limits.md).
