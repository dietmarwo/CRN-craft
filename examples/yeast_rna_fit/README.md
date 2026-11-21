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
licensed GPL-3.0; the derived tables are kept here pending a decision on
their redistribution (page 23c). The checks fetch nothing.

[Paper](https://doi.org/10.1073/pnas.1804060115) · [Author archive](https://github.com/MunskyGroup/Munsky_et_al_2018/tree/1aa0b5ce25d21484c7c4e8a894477b1719209ce7)
