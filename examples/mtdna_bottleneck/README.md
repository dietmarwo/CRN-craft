# Birth, death and partition

The example of [tutorial S6](../../docs/tutorials/simulation/S6-mtdna-bottleneck.md):
copy number and heteroplasmy variance through mouse development, fitted
to Table 1 of Johnston et al. (2015), *eLife* 4:e07464,
[doi:10.7554/eLife.07464](https://doi.org/10.7554/eLife.07464).

```sh
python examples/mtdna_bottleneck/study.py --size quick
python examples/mtdna_bottleneck/study.py --size full
python examples/mtdna_bottleneck/study.py --size full --free-final-degradation
python examples/mtdna_bottleneck/review.py
python examples/mtdna_bottleneck/study.py --figures docs/guide/img/mtdna_bottleneck
```

Run outputs go to `runs/<size>` (with `-free-final` for the sensitivity
run). Numba is optional; without it only the quick size is practical.

| File | Content |
|---|---|
| `model.py` | CRN-craft builder and lineage driver, exact founder-family sampler, moment recursion, ABC distance and samplers |
| `study.py` | ABC fit of Table 1 at distance 40, Figure 3 bands, HB group comparison, plot |
| `review.py` | first-order variance against exact lineages, early-age timing, ensemble agreement, phase-six degradation |
| `fetch.py` | converts the paper's XML (Table 1) and the Figure 4 source spreadsheet |
| `digitize.py` | pixel readings of the Figure 3 central curve and envelopes (unchecked; see the chapter) |
| `reference-table1.csv` | 74 rows of Table 1; the cell-count rows at 7, 14 and 21 h are placed at the division times |
| `reference-table2.csv` | the 15 group summaries of the Figure 4 spreadsheet |
| `reference-hb.csv` | the 339 single-cell heteroplasmy values of the Figure 4 spreadsheet |
| `reference-figure3.csv` | digitized central values and envelopes at 6 ages |
| `results/abc.csv` | acceptance, distances and effective sizes per ensemble |
| `results/trajectories.csv` | quantiles of mean copies and `V'(h)` over accepted sets |
| `results/figure3-comparison.csv` | the 12 anchors against the reconstructed bands |
| `results/hb-comparison.csv` | the 15 groups against population and group-sized bands |
| `results/review-*.csv` | the review checks described in the chapter's results page |

The `holdout` column of `hb-comparison.csv` marks the two litters listed in
`manifest.json`. No HB group enters the fit, so all 15 groups are a new
validation. The group-sized band uses a founder law matched to two
moments; it is not an exact sample of the six-phase law.

The reference files are derived from the paper and its source data,
published under the Creative Commons Attribution 4.0 license; cite the
paper when reusing them.
