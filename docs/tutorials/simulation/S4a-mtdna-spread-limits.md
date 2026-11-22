# Choices, checks and limits

This page belongs to [tutorial S4](S4-mtdna-spread.md). It lists what the
paper and the authors' files leave open, the checks that guard the
example, and the review experiments on the Figure 2E speed.

## Reconstruction choices

| Choice | Value used | Alternative | Measured effect |
| --- | --- | --- | --- |
| Figure 1C `Nss` and `c` (not in the archive headers) | 150 and 0.001/day in every panel | equal density with `Nss = 200` | with 200, the per-region spread at day 700 is 15% below the archive's; with 150 it is within 2% |
| Figure 1C start | `h = 0.5` on the steady state | other copy numbers | not run |
| Faster mutant degradation | `0.07 * 0.001` per day extra | header "0.001 µ" read otherwise | the archived deterministic curve is not reproduced (−0.005 at day 700) |
| Hop rate | `gamma` per copy and directed edge | `gamma` in total | the archived C code and notebook use per edge |
| Birth-rate clipping | `max(0, …)` | the archived C code never clips | none: the clip needs a fluctuation of about 9 (Figure 1C) or 40 (Figure 2E) standard deviations |
| Figure 2E `Nss` | 3,500 (Figure 2 caption) | 3,000 (preprint supplement) | 0.143 vs 0.146 µm/day, within cohort noise |
| Chain length and spacing | 500 regions of 30 µm | not stated in the paper | the archived day-0 window matches this geometry |
| Left boundary | fixed reservoir at `h = 1` | dynamic region held at `h = 1` | none expected: diffusion reaches about 520 µm |
| Window alignment | both archive files start at the same region | shifted by one region | one region moves the archived speed by 0.020 µm/day; not testable from the files |
| Wave estimator | crossing of the ensemble-mean profile, three-point smoother | no smoother | at most 0.004 µm/day |
| Wave simulator | approximate steps of 0.0625 to 0.5 days | exact path simulation | eight exact fibers: 0.160 µm/day (99% interval 0.064–0.199), inconclusive; 800 fibers outstanding |

## The approximate wave scheme

`wave.py` advances all regions in steps of `dt`. Within a step it holds
the birth rate of each region at its start-of-step value. Births are
Poisson draws from the start-of-step counts. Deaths and hops compete for
the same start-of-step copies, so each copy leaves at most once per step,
and copies born in the step can neither die nor hop. Counts cannot become
negative. The scheme is first order in `dt`.

## Review experiments on the speed gap

Targets and decision rules were written down before each run. The rows are
in `results/review-experiments.csv`; `review.py` reruns them.

- **Deterministic part** (`review.py deterministic`). The noise-free chain
  moves its half-height crossing at 0.04199 µm/day for `Nss = 3,500` and
  0.04196 for 3,000. The mean of the step scheme gives 0.03932 at a
  0.5-day step and 0.04130 at 0.125 days. Step bias in the mean is
  therefore at most 0.003 µm/day.
- **Copy number** (`review.py cohort`, 96 fibers, 0.125-day step).
  `Nss = 3,000` gives 0.14293 µm/day (99% interval 0.12424–0.16245). The
  interval excludes 0.18973, so the copy number alone does not explain the
  gap. Label: unresolved discrepancy.
- **Smaller step** (32 fibers, 0.0625-day step). 0.14921 µm/day
  (0.12571–0.17367): the gap persists at half the smallest delivered step.
  Label: unresolved discrepancy.
- **Step effect on the whole profile** (32 fibers at 0.5 and 0.0625 days).
  The mean area gained by the mutant, divided by 1,480 days, is 0.12816
  and 0.14140 µm/day. The difference, −0.0132 (99% interval −0.0266 to
  +0.0002), includes zero but is not shown to be below 0.02. Label:
  inconclusive.
- **Noise of the scheme against the exact simulator** (`review.py noise`;
  16 exact fibers, 64 approximate fibers per step, 20 days, full geometry).
  The variance of the regional mutant fraction over the 274 front regions
  has a ratio approximate/exact of 0.987 (99% interval 0.949–1.189) at
  0.125 days and 0.981 (0.940–1.178) at 0.5 days. The scheme shows no
  sign of distorting the local noise that drives the spread.
- **Equal-density setting** (`review.py equal`). See the first row of the
  choices table; the shipped Figure 1C results use `Nss = 150`.

- **Exact pilot** (`study.py --source-chain --wave-sampler local --fibers 8`).
  Eight exact fibers of the full chain, fibers 0 to 7 of the 800-fiber
  design. Half-height speed 0.15951 µm/day (99% interval 0.06445–0.19865,
  no ambiguous resample). Label: inconclusive. Without the smoother, the
  eight-fiber mean profile crosses 0.5 twice (at 5,699 and 5,731 µm), and
  100 of 500 resamples are ambiguous. The fibers are in
  `results/exact-wave-fibers.csv`.

What remains: the exact simulator over 1,480 days was run for eight
fibers, not 800, and the authors' Figure 2E run itself is not documented in the
archive (no run time, no step, no window position). Before the full exact
run, the most informative check is whether the authors' two profile
windows start at the same region.

## Checks

`validation/check_mtdna_spread.py` runs in a few seconds without network
access. It checks:

- the reference reader: 48 Figure 1C rows and the three day-700 means;
  139 + 140 Figure 2E rows, the crossings 1,695.52 and 1,976.33 µm, and the
  speeds 0.18973, 0.2035 and 0.1747 µm/day at 0.5, 0.4 and 0.6;
- the exact law of one copy hopping between two regions (finite state
  projection, `(1 - exp(-0.8)) / 2` at day 2) and conservation of each
  genotype under migration alone;
- the row counts of every saved result table;
- that both tables of tutorial S4 equal the saved results;
- for the extension of pages 19b and 19c: the exact bookkeeping laws of the
  lineage sampler, a moment law, the quick calibration, the saved tables
  and the two generated tables of page S4c.

## Limits

- The paper's model has no de novo mutations: every run starts with
  mutants present. Pages [S4b](S4b-mtdna-spread-de-novo.md) and
  [S4c](S4c-mtdna-spread-de-novo-results.md) add them as an extension of
  this guide.
- The published supplement was not reachable for this review. Chain length,
  spacing, initial profile and boundary handling come from the first
  preprint's supplement; they match the archived day-0 profile exactly.
- The archive's Figure 1C files do not state `Nss` or `c`. The values used
  reproduce the archived means and per-region spread, but they are
  inferred.
- The archived speed has its own sampling error from 800 fibers, which the
  archive does not report. From the spread of the fresh fibers it is about
  0.002 µm/day, too small to close the gap.
- The upstream repository has no license file. The two reference tables
  are short numerical extracts of its data. They are not covered by this
  repository's MIT license and stay under the authors' terms
  ([notice](../../../THIRD_PARTY_NOTICES.md#reference-data-of-the-guide-case-studies));
  `fetch.py` regenerates them from the pinned revision.
- The rhesus speed of Figure 2B is an observation from another study. The
  data are not in the archive and were not refitted.

## Exercises

1. Change `--step` in `wave.py` and plot the half-height speed against the
   step. How many fibers do you need before the step effect is larger
   than the cohort noise?
2. Remove the noise: `review.py deterministic` gives the speed without
   survival of the densest. Why does the crossing at 0.4 move faster than
   the crossing at 0.6?
3. Run `study.py --source-chain --wave-sampler local --fibers 1` (about
   5 minutes) and compare the single fiber's crossing with the spread of
   the approximate fibers.
