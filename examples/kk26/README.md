# Cross-sectional mtDNA parameter inference

[Guide chapter 18](../../docs/guide/18-mtdna-inference.md) reproduces parts
of Kowald & Kirkwood (2026) at 1,000 mtDNA copies per cell. The Python driver
uses native CRN-craft for small-law checks, the paper's Moran model for
synthetic snapshots, and fcmaes differential evolution for fitting.
[Paper and methods](https://doi.org/10.1038/s41514-026-00431-4)

Run from the repository root with Python 3.10 or newer, a C++17 compiler
on Linux or macOS, and these Python packages. `CXX` can select the compiler.

```sh
python -m pip install 'crn-craft==0.8.0' fcmaes scipy matplotlib
python examples/kk26/study.py bootstrap --size quick --threads 2
python examples/kk26/study.py bootstrap --size full --threads 8
python examples/kk26/study.py fit examples/kk26/results/reconstruction/snapshots-baseline.csv
python examples/kk26/study.py evidence --threads 8 --complete-grid --results examples/kk26/results/local/full
python examples/kk26/study.py figures --images docs/guide/img/kk26
python validation/check_kk26.py --output kk26-validation.json
```

The quick run fits four 2,000-cell cohorts in the baseline regime. The
full run fits 40 cohorts of 10,000 cells in each regime. The example
keeps synthetic cohort counts and forward-grid nodes by configuration
and seed, and reuses completed data when a run is repeated. The saved
grids in `results/reconstruction/grids` make the default fit inexpensive.
A different set of ages requires a new grid. The evaluation budget is
passed to fcmaes, and the fit rows give actual evaluations, which can
extend to the end of a generation. The evidence command
recomputes the direct and implicit comparison curves and simulates
100,000 ATP cells at a 1 h step and 30,000 at a 0.5 h step. The complete
baseline lattice is an additional run.

## Fit and output data

Input CSVs for `fit` have columns `age_days` and `bin_0` through `bin_4`.
The bins used for fitting are exactly zero, (0%, 5%], (5%, 20%],
(20%, 50%] and (50%, 100%]. Counts are nonnegative integers.
Ages increase, and all five counts at an age share one denominator.

| Regime | Mutation probability bound | Selection advantage bound | Advantageous fraction bound | Truth (p, s, f) |
| --- | --- | --- | --- | --- |
| Baseline | 3.8e-6–1.12e-5 | 0.052–0.21 | 0.32–1 | (6.8e-6, 0.10, 0.80) |
| Strong | 1.65e-5–5.2e-5 | 0.26–0.87 | 0.065–0.47 | (3e-5, 0.50, 0.20) |

| Saved output | What to read |
| --- | --- |
| `fig6-curves.csv`, `fig6-differences.csv` | The simulated and numerically propagated threshold curves and their differences from the digitized published values |
| `fits.csv`, `bootstrap.csv`, `parameter-differences.csv` | Individual parameter fits, their seed and objective, summary intervals, and comparisons to the printed tables |
| `fraction-profile.csv` | The least-squares objective with the advantageous fraction fixed. The mutation probability and selection advantage are fitted again. |
| `grid-checks.csv`, `forward-checks.csv` | Off-grid checks, direct checks at selected fitted coefficients, numerical refinement, neutral-count leakage, and differences between clocks |
| `fig5-passage.csv`, `takeover-summary.csv` | Takeover CDFs, success probability, duration quantiles, modal bins, finite-window coverage and remaining tails |
| `manifest.json`, `grids/` | Reproducible choices and the prepared forward data. Wall time is written separately for live runs. |

## Numerical and source details

Equation 7 uses independent choices of a reproducing and a replaced
copy. Exact discrete simulation skips unchanged ticks with geometric
waiting times. The full-scale master-equation reference propagates the
six moves and the stay probability. The accelerated solver uses a
Poisson clock, symmetric splitting and Crank–Nicolson along each of
the three axes of the triangular lattice. Its default substep is five
days. Caps of 192 and 320 on the neutral count, respectively, apply
to the two regimes. Escaped probability stays excluded from the
reported tails. The stored fit uses tensor-product cubic interpolation
over 192 and 216 forward nodes. See the two numerical check files.

Both panels of Figure 6 are in a 2,050 × 726 px image extracted from
the article PDF. The digitization uses the 0–1,000 day ticks and
ordinate ticks 0.01–0.05 to calibrate each panel. Read the traces
at the published 30-day ages. Source reading intervals allow
3–5 ordinate pixels. They describe the process of reading the image,
while the bootstrap describes the synthetic cohort variation.
The bins in Figure 5 are labeled by their lower day.
[Article PDF](https://europepmc.org/articles/PMC13272806?pdf=render)

The hourly ATP reconstruction uses the ATP constants of the 2014
model, survival and birth probabilities derived from hourly hazards,
and a closed-form ATP update with the post-reaction counts held
fixed. A cell that reaches zero ATP retains its counts. Only
wild-type daughters can mutate in this reconstruction. The authors'
source repository returned a service error during preparation.
[ATP model](https://pmc.ncbi.nlm.nih.gov/articles/PMC3939916/)
· [Published data and code](https://gitlab.uni-rostock.de/IBIMA/MitoMutantAccTimes1)
