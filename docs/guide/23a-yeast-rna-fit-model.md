# 23a. Data and model

[Chapter 23](23-yeast-rna-fit.md) · Data and model ·
[Fits and predictions](23b-yeast-rna-fit-results.md) ·
[Choices, checks and limits](23c-yeast-rna-fit-limits.md)

## The cells

Yeast cells were shocked with 0.2 M or 0.4 M NaCl and fixed at 0, 1, 2,
4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50 and 55 minutes. Single-molecule
FISH counted the STL1 RNA in the nucleus and in the cytoplasm of each cell.
`fetch.py` reads the authors' array `Hog_Data/Analyzed_FISH_Data_2016_06.mat`
(the file their fitting code loads) and writes one row per distinct
(dose, replicate, minute, total, nuclear, cytoplasmic) tuple with its cell
count to `source-counts.csv` (7,374 rows). It records the SHA-256 of every
file it reads.

| Dose | Replicate | Cells | Cells after time zero | In the authors' fit |
| --- | --- | ---: | ---: | --- |
| 0.2 M | 20140425 | 14,382 | 12,452 | yes |
| 0.2 M | 20140429 | 11,129 | 9,101 | yes |
| 0.4 M | 201301004 | 18,471 | 17,544 | yes |
| 0.4 M | 20131022 | 7,457 | 7,089 | no |
| 0.4 M | 20160509 (no 1-minute sample) | 14,204 | 13,371 | yes |
| All | | 65,643 | 59,557 | |

Three facts about these rows:

- **Cell totals.** The paper reports 65,454 cells, 25,511 at 0.2 M and
  39,943 at 0.4 M. The mature-RNA array holds 25,511 and 40,132. The
  authors' transcription-site array (`RNA_TS_Brian.mat`) holds the paper's
  totals. The two arrays agree cell for cell in 0.2 M and in 20160509. The
  other two 0.4 M replicates differ at almost every time point:
  20131022 has 1,051 more cells in the site array, 201301004 has 1,240
  fewer. The 189 cells are the net of these two differences, not one
  duplicated sample. The example uses the mature-RNA array, as the
  authors' fitting code does.
- **Replicates in the fit.** The authors' fitting code
  (`Process_Hog_Data_Class.m`, option `_AC`) pools both 0.2 M replicates
  but only the first and third 0.4 M replicate. Replicate 20131022 is not
  in their fit. This example fits all five replicates (see
  [the choices](23c-yeast-rna-fit-limits.md#reconstruction-choices)).
- **Totals.** In 235 cells (228 rows) the recorded total is one or two
  smaller than nuclear plus cytoplasmic. The total-RNA fits use the
  recorded total; the authors' code uses nuclear plus cytoplasmic. The
  difference changes the scores by 0.0002 nats per cell.

Time zero is shown in the figures and left out of both objectives. The
authors' likelihood code also sums from the second time point on; the
paper does not say so. All fitted models start with no RNA in promoter
state 1; 99.0% of the 6,086 time-zero cells have no RNA.

The paper calls a cell ON when it has more than three RNAs (its Methods
also say "three or more"). The example uses `X > 3`. Recounted from the
rows, the ON fraction peaks at 0.699 at 15 minutes for 0.2 M and at 0.897
at 20 minutes for 0.4 M.

## The model

The STL1 promoter switches between four states `S1` to `S4`. In state
`Si` it starts new RNAs at rate `k_i` (per second). In the nonspatial
model an RNA decays at rate `γ`:

    S1 <-> S2 <-> S3 <-> S4          rates k12, k21(t), k23, k32, k34, k43   [1/s]
    Si -> Si + RNA                   rate  k_i, i = 1..4                     [1/s]
    RNA -> 0                         rate  γ per molecule                    [1/s]

Only the switch from `S2` back to `S1` depends on the stress. It falls as
the nuclear Hog1 kinase signal rises, after a delay `t0`:

    k21(t) = max(0, k21⁰ − k21¹ · H(t − t0))
    H(s)   = A · [ u(s) / (1 + u(s)/M) ]^η,   u(s) = (1 − e^(−r1·s)) · e^(−r2·s),   s ≥ 0

`H` is dimensionless. Its constants were fitted by the authors to measured
Hog1 nuclear enrichment and are the same for both doses except `r2`:
`r1 = 6.09e-3 /s`, `r2 = α/(c + 0.0410)` with `α = 1.68e-3 M/s` and the
NaCl concentration `c` in M (so `r2 = 6.97e-3 /s` at 0.2 M and
`3.81e-3 /s` at 0.4 M), `η = 5.89`, `M = 0.0215`, `A = 9.32e9`. These
match Table S2 of the supplement to its two printed digits; the example
uses the full-precision values of the authors' parameter file. The
nonspatial model has 13 free coefficients: six switching rates and
`k21¹`, four initiation rates, `γ` and `t0`.

The spatial model splits RNA into nuclear `N` and cytoplasmic `C`.
Initiation makes `N`, export `N -> C` has rate `k_exp`, and the two
compartments decay at `γ_nuc` and `γ_cyt`. That is 15 coefficients.

`fetch.py` maps the authors' 39-entry parameter vector to these names by
the layout of their `PARS_Definition_Connected_Models.m`. Free entries 1,
4, 5, 8, 9, 12 (switching), 16 (signal on `k21`, stored with a negative
sign), 25–28 (initiation), 33 (decay) and 39 (delay); the spatial model
adds 35 and 37. `source-parameters.json` holds the four archived vectors.

## How it is solved

The time-varying rate rules out a single matrix exponential. The example
uses three local solvers; none is part of CRN-craft.

- **Total RNA** (`cme.cpp`): the master equation of the 4 × (cap + 1)
  states of promoter and count, integrated with fixed-step BDF2 (one
  backward-Euler start step) and a banded direct solve, with `k21` set at
  the end of each step. Initiation at the cap leaves the domain. Nothing
  is stored as a separate sink; the leakage is `1 − Σp`, which is the mass
  that crossed the cap because BDF2 conserves the sum of a conservative
  system. The fits use steps of 12 s (draft) or 6 s (full) and caps of
  250 or 280 RNAs. Every reported score is recomputed with a 1.5 s step
  and the cap doubled from 360 until the leakage is below 1e-5.
- **Means**: the first moments close exactly (all rates are linear in the
  counts). The fit integrates them with the same BDF2; the reported mean
  costs use SciPy's adaptive BDF with relative tolerance 1e-7.
- **Joint nuclear and cytoplasmic RNA**: SciPy's sparse BDF on the
  4 × 22 × 241 joint states (relative tolerance 2e-6), checked with
  4 × 29 × 306 states and tolerance 4e-7.

The library is built with `-std=c++17 -O3 -ffp-contract=off` and no
`-march=native`, so compilers do not fuse multiply-adds differently on
different machines.

The reusable projection in `native_fit.py` reads the same model and
prescribed rate through CRN-craft. Its checks compare it with the two
original solvers. Both share the input waveform function. Each original
solver has its own integrator: the total-RNA one fixed-step BDF2, the
joint one SciPy's sparse BDF, while the projection uses CRN-craft's
native variable-order integrator (`solver="bdf"`). [Prepared probability laws](../reference/experiments.md#prepared-probability-laws)

What the check verifies about these solvers:

- With a constant signal the model is an ordinary reaction network. The
  check builds it with `crn_craft.Model`, solves it with CRN-craft's finite
  state projection, and compares with an explicit generator and SciPy's
  `expm_multiply`, for total RNA (cap 5) and for the joint model (caps 3
  and 4), to 5e-10 in probability and leakage.
- The BDF2 solver with 0.05 s steps agrees with that law to 2e-5 after
  one minute, and the BDF2 means agree with the adaptive moment equations
  to a relative 1e-4.
- The histogram score of a three-bin example equals the hand calculation
  above to 1e-12.

## The two objectives

Both objectives use the same groups: one dose, one replicate, one time
after zero (79 − 5 = 74 groups).

- **Means.** For group `g` with `n_g` cells, sample mean `x̄_g` and sample
  variance `s_g²`, the cost is `Σ_g (x̄_g − μ_g(θ))² / (2 v_g)` with
  `v_g = s_g²/n_g + 1/(12 n_g)`. `μ_g` is the model mean. The paper's
  Equation 23 uses `s²/N` without the last term; the term keeps groups
  with zero variance finite. The joint model uses the 2 × 2 covariance.
  The paper pools the replicates of one time into one group; the example
  keeps them separate.
- **Histograms.** `−Σ_g Σ_x n_gx log p_g(x; θ)` over the observed counts
  `x`. Each count below the cap is its own bin, including zero; cells
  above the cap form one tail bin with probability `1 − Σp`. If an
  occupied bin gets probability zero, or a tail bin is occupied and no
  mass left the cap, the candidate is rejected as invalid (cost +∞) and
  counted. In the saved full fits no candidate was invalid; in the draft
  fits at most 2 of about 500 evaluations per restart were.

Scores are divided by the number of cells to give nats per cell.
