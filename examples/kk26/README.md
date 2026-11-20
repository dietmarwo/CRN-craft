# Inferring mtDNA takeover times from cell snapshots

[Guide chapter 18](../../docs/guide/18-mtdna-inference.md) reconstructs
Equation 7, Figures 5 and 6 and the parameter recovery of Tables 2 and 3 of
Kowald & Kirkwood (2026), *npj Aging* 12
([doi:10.1038/s41514-026-00431-4](https://doi.org/10.1038/s41514-026-00431-4)).
Its [last page](../../docs/guide/18c-mtdna-inference-limits.md#reconstruction-choices)
lists every reconstruction choice in one table. The paper's authors did not
write or review this example.

Run the commands from the CRN-craft root with Python 3.10 or newer, a
C++17 compiler on Linux or macOS (set `CXX` to choose it), and these
packages:

```sh
python -m pip install 'crn-craft==0.8.0' 'fcmaes==2.2.0' scipy matplotlib
python examples/kk26/study.py bootstrap --size quick --threads 2
python examples/kk26/study.py fit examples/kk26/results/reconstruction/snapshots-baseline.csv
python examples/kk26/study.py figures --images docs/guide/img/kk26
python validation/check_kk26.py --output kk26-validation.json
```

`quick` fits four baseline cohorts of 2,000 cells. `fit` fits a CSV of
counts (columns `age_days`, `bin_0` to `bin_4`) with the saved baseline
grid; the ages must be 0, 30, …, 1,080 days, otherwise a new grid is built.
`figures` draws the two guide figures from the saved results.

The full study and its evidence:

```sh
python examples/kk26/study.py bootstrap --size full --threads 8
python examples/kk26/study.py evidence --threads 8 --complete-grid --results examples/kk26/results/local/full
python examples/kk26/review.py budget
python examples/kk26/review.py wide-box
python examples/kk26/review.py atp-fit
python examples/kk26/review.py design --threads 2
```

`bootstrap --size full` fits 40 cohorts of 10,000 cells in each regime.
Forward grids are cached under `results/local/grids` by a key of the
solver source, the build flags and the numerical configuration; a matching
grid in `results/reconstruction/grids` is used instead of building one.
Cohort counts are cached by seed and solver source, so a repeated or
interrupted run reuses them; grid nodes are saved one by one while a grid
is built. Thread count does not change any result. `evidence` recomputes
the discrete reference, the solver and corner checks, the complete lattice
(`--complete-grid`) and the 100,000 hourly and 30,000 half-hour ATP cells.
`review.py` reruns the review checks of guide page 18b and writes
`review-*.csv` to `results/local/review` (`--results` changes it).
`wide-box` solves the forward law directly, 15 to 20 minutes per cohort on
one core; it memoizes the solves under `results/local/review-cache`, so an
interrupted run continues, and `--cohort strong:2` selects one cohort.

| Saved output (`results/reconstruction/`) | What it records |
| --- | --- |
| `fits.csv` | The 80 fits: cohort seed, parameters, objective at the fit and at the truth, evaluations, optimizer status and seed |
| `bootstrap.csv`, `parameter-differences.csv`, `recovery-table.md` | Mean, CV, percentiles, search box and fits near its edge; the differences from Tables 2 and 3; the guide table |
| `snapshots-*.csv` | The first cohort of each regime as five bin counts per age, an input for `fit` |
| `fraction-profile.csv` | The objective of the first cohort with `f` fixed and `p`, `s` refitted |
| `grid-checks.csv` | The spline against the solver at 3 fixed points, 16 points over the whole box and 3 fits per regime |
| `forward-checks.csv` | Solver against discrete reference at the truths and two box corners, neutral cap against the complete lattice, Figure 6 differences, ATP step comparison |
| `fig6-curves.csv`, `fig6-differences.csv` | The Figure 6 curves of every backend and their differences from the digitized points |
| `fig5-passage.csv`, `takeover-summary.csv` | Takeover CDFs, success probabilities, quantiles, modal bins and the share reached by day 1,080 |
| `fig5-differences.csv` | The published bars of Figure 5 and the reconstruction's probability of the same 30-day bins |
| `review-*.csv` | The review checks: larger optimizer budget, wider box, Model-1 data, independent cells per age |
| `manifest.json`, `grids/` | Configuration, seeds, solver hash, build flags and the two forward grids |

## Reproducibility

The C++ library is built with `-std=c++17 -O3 -ffp-contract=off` and no
`-march=native`: fused multiply-adds change the grids in the last bits,
and the fits follow such changes (by about 3e-8 relative). The saved
results were produced on Intel x86-64 with GCC 13.3, Python 3.13, NumPy
2.3.5, SciPy 1.18.1 and fcmaes 2.2.0; fcmaes 2.1.0 gives the same fits
byte for byte. Other compilers, math libraries or optimizer versions may
give statistically equivalent, not identical, cohorts and fits: the
sampler uses `log1p`, the ATP model `std::binomial_distribution`.
`validation/check_kk26.py` therefore compares no simulated value with a
saved one.

## Reference data

- `paper-figure6.csv`: 296 points of both panels of Figure 6, read from
  the 2,050 × 726 px figure image of the article PDF. Each panel is
  calibrated on its 0–1,000 day and 0.01–0.05 ticks; the curves are read at
  the 30-day ages with a reading interval of 3 to 5 pixels. Panel B was
  checked against the copy of Figure 6B in the authors' companion preprint
  (bioRxiv [10.64898/2026.01.30.702741](https://doi.org/10.64898/2026.01.30.702741),
  Figure 4B): 132 readable points agree within 1.05 pixels of that image
  (0.017 percentage points). Both panels were read again from the figure
  in the preprint of the paper (bioRxiv [10.64898/2026.01.30.702736](https://doi.org/10.64898/2026.01.30.702736), version 1; 878 × 321 px): after day 0, 127
  points per panel are readable; panel A agrees within 1.05 pixels of that
  image (0.023 percentage points), panel B within 1.7 pixels.
- `paper-tables.csv`: the means, coefficients of variation and 95%
  intervals of Tables 1 to 3, transcribed from the article PDF. All 27
  values and the nine true values agree with the tables of the preprint
  of the paper (version 1).
- The modal bins of Figure 5 (990, 540, 390 and 240 days, in `evidence.py`)
  are the four labels printed in that figure.
- `paper-figure5.csv`: the 60 bars of Figure 5, read from the 768 × 465 px
  figure image of the preprint of the paper. The ordinate is calibrated on
  the grid lines 0.00 (row 419) and 0.35 (row 12); a bar's value is the row
  of its top edge, with a reading interval of one pixel (0.00086). A bar
  is assigned to its 30-day bin by its position in the group of four.
- The authors' repository (linked from the paper) returned server errors
  when the example was written.
