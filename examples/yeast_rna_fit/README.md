# Fit single-cell RNA count distributions

[Chapter 23](../../docs/guide/23-yeast-rna-fit.md) explains the study: the
data, the model, what the fits show and what they do not. The scripts
need NumPy, SciPy, pandas, fcmaes, Matplotlib, CRN-craft 0.8.0 and, for the
fits, a C++17 compiler (`CXX`). `cme.cpp` solves the total-RNA master
equation with fixed-step BDF2 and a banded solve; it is built with
`-std=c++17 -O3 -DNDEBUG -ffp-contract=off`. The joint nuclear and
cytoplasmic model uses SciPy's sparse BDF.

```sh
python examples/yeast_rna_fit/study.py --size quick
python examples/yeast_rna_fit/study.py --size draft --workers 4 --cross-validate --budget-check
python examples/yeast_rna_fit/study.py --size full --workers 4 --budget-check
python examples/yeast_rna_fit/study.py --size full --spatial-only --spatial-fit --workers 4
python examples/yeast_rna_fit/review.py start
python examples/yeast_rna_fit/review.py choices
python examples/yeast_rna_fit/study.py --figures docs/guide/img/yeast_rna_fit
```

Both nonspatial fits vary all 13 coefficients in one multiplicative box
centered on the authors' archived histogram fit, and each restart first
evaluates that archived vector. The comparison therefore ranks models near
the published fit; it is not a search from no prior knowledge (guide page
23b). The spatial run varies two of 15 coefficients on a 3 × 3 grid; its
best histogram point is the authors' archived spatial vector, so the
spatial model is evaluated, not refitted. The transcription-site
prediction recomputes the authors' prediction from that archived vector.

| File | Content |
| --- | --- |
| `source-counts.csv` | one row per (dose, replicate, minute, total, nuclear, cytoplasmic) with its cell count, from `Analyzed_FISH_Data_2016_06.mat` |
| `source-audit.csv` | cells per sample and cells whose total differs from nuclear plus cytoplasmic |
| `source-sites.csv` | brightest-site intensity bins per cell, from `RNA_TS_Brian.mat` |
| `reference-site-*.csv` | numerical series saved in the authors' Figure 4 `.fig` files |
| `source-parameters.json` | the four archived parameter vectors and the SHA-256 of every source file |
| `results/scores.csv`, `results/optimizer-work.csv` | scores and restarts of the full and draft runs (`preset` column) |
| `results/parameters.json` | the two full-size fits |
| `results/histograms.csv`, `results/model-moments.csv` | predicted laws and summaries for the figures |
| `results/spatial-*.csv`, `results/spatial*.json` | the spatial grid and its accounting |
| `results/sites.csv`, `results/site-activity-comparison.csv` | the site prediction |
| `results/review-optimizer-start.csv` | the search-start review (page 23b) |
| `results/review-data-choices.csv` | saved vectors scored on other data choices (page 23c) |

`study.py` writes `scores.csv` and `restarts.csv` per run; the shipped
`scores.csv` and `optimizer-work.csv` join the full and draft runs and add
the `preset` column.

`fetch.py` rebuilds the source tables from the pinned authors' archive and
records the SHA-256 of each MATLAB file. The authors' repository is
licensed GPL-3.0. The derived tables (`source-*.csv`,
`source-parameters.json`, `reference-site-*.csv`) are not covered by this
repository's MIT license and stay under the terms of their source; see
[THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md#reference-data-of-the-guide-case-studies). No upstream code is included. The checks fetch nothing.

[Paper](https://doi.org/10.1073/pnas.1804060115) · [Author archive](https://github.com/MunskyGroup/Munsky_et_al_2018/tree/1aa0b5ce25d21484c7c4e8a894477b1719209ce7)


## Prepared projection and fitting

`native_fit.py` runs on the native engine (`solve(solver="bdf")`) and
needs SciPy only for the chapter's reference solvers it compares with.
`native_fit.py check OUTPUT`
solves both readouts at both salt doses at the archived vectors with
`prepare_fsp` and compares them with the chapter's solvers. The
total-RNA reference is the fixed-step BDF2 solve of `cme.cpp`. The joint
reference is the chapter's SciPy sparse BDF with the same tolerance and
input waveform; the projection integrates with CRN-craft's native
variable-order integrator, so the check compares two integrators on the
same state space. The check requires a maximum marginal CDF difference of
0.001 and a score within 0.01 nats per cell of the recorded 2.04126
(total RNA) and 2.65027 (joint counts). Total-RNA leakage must stay
below 1e-5.

The total-RNA projection has 1,444 states and takes a few seconds. Its
CDF differences of about 1.8e-4 come mostly from the reference's 1.5 s
step. At a 0.25 s step they fall to about 4.6e-5. Tightening the
projection tolerance from 1e-8 to 1e-10 changes its law by less than
4e-8. The joint projection has 35,496 states. One salt dose (tolerance
4e-7, `max_step` 20 s) took 77 s on one core in one measurement, with 638
steps and 126 sparse LU factorizations, and peaked at about 420 MB of
resident memory.

`fit OUTPUT` (joint model, `--budget 0`) scores the chapter's 3 × 3 grid
in the transcription rate `i3` and the decay rate on the (21, 240)
projection, plus the center and one random point, and records the time
of each evaluation. With a positive `--budget` it then starts searches
in all 15 coefficients. **The 15-coefficient joint refit has not been
completed.** Attempts with this script ended with exit status 137
(the process was killed) after the grid pilot and with three workers.
A single-worker draft (`--draft-joint`, search on the (11, 120)
projection) produced a candidate that scored worse than the archived
vector on the refined projection. Each worker thread holds its own
solver workspace, and every candidate is rescored on the refined
projection, so plan at least 1 GB of memory per worker.

`--total --holdout` runs the total-count searches with the two chapter
splits. These searches start at the archived vector inside the
archive-centered box. That vector was fitted to all replicates,
including the held-out ones. A short budget therefore exercises the
plumbing; it is not a fit and not a held-out comparison. A search that
does not improve on its start returns the start vector. Add
`--training-only` to center the box on a fixed preset instead of the
archived vector. The original example's data and likelihood
conventions apply. Results include the actual number of evaluations,
which can exceed `--budget`.
