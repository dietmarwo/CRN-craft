# Saved C5 tables

`draw(tables, folder)` in `publications/rna_counts.py` renders
`img/C5-contour.png` and `img/C5-recovery.png` from these tables. All data
are regenerated with CRN-craft; the paper's random samples are not
available.

## recovery.csv

One row per optimizer start: run, sparse-data score (`midrank`,
`published` for the individually ranked Equation 9 reading, `population`
for the 10,000-trajectory ECDF), calls, cells, dataset index and seed,
start, fitted rates, best score, the score of the generating rates
(0.0277, 0.1667, 0.4) on the same data, and whether the best score is at
or below that value (`reached_truth_score`; `false` means the search
stopped short of the optimum).

- `run=delivered-300`: recounted from the implementer's study runs on 8
  threads, not rerun here: `run.py C5 --preset study` with 300 calls per
  start, and the same with `--published-rank`. DE population 30, seed
  `data_seed + 91 + start`. Scores agree with the minima of the saved
  traces. The 10,000-trajectory fit is shared by both runs and listed once.
- `run=review-3000`: rerun in review, one thread, with the tutorial's own
  model, data generator, score and search driver and the study seeds, 3,000
  calls per start: the 100-trajectory datasets 0 to 3 (start 0) and the
  10-, 20- and 10,000-trajectory datasets. These are the arms that
  `run.py C5 --preset study --datasets 4 --starts 1` runs; the full command
  was not run here (about 21,000 calls).

## contour-m10.csv, contour-m100.csv

Midrank score on a 31 x 31 grid of k2 and k3 in [0.1, 1] with k1 = 0.0277,
the setting of the paper's Figure 2A, over all 150 times. `m10` uses the
10-trajectory dataset of seed 5203 (the tutorial's contour data), `m100`
the 100-trajectory dataset 0 of seed 5293. Rerun in review, one thread.
