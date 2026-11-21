# Mitochondrial network and early drift

The model, the protocol and the results are explained in
[chapter 20](../../docs/guide/20-mtdna-network.md) and
[its results page](../../docs/guide/20a-mtdna-network-results.md).
Time is in days.

```sh
python examples/mtdna_network/study.py --size quick
python examples/mtdna_network/study.py --size full
python examples/mtdna_network/review.py --size full
python examples/mtdna_network/study.py --figures docs/guide/img/mtdna_network
python examples/mtdna_network/review.py --figures docs/guide/img/mtdna_network
```

Runs write to `runs/<size>/` (or `--output`); `--figures` reads `results/`.

| File | Content |
|---|---|
| `conditions.csv` | parameters and initial counts of the 39 source conditions, from the authors' parameter tables |
| `reference-variance.csv` | the authors' archived variance series (10,000 paths, every 10 days to day 500) |
| `fetch.py` | regenerates both files from the pinned upstream revision (`--online-summaries` for the series) |
| `manifest.json` | primary conditions, seed, report times and decision rule, fixed before the runs |
| `results/slopes.csv` | one row per condition: Equation 13 slope, measured slope, ratio and 99% interval |
| `results/variance.csv` | the early variance series of every condition |
| `results/nominal-long.csv` | a 96-path H:60 run to day 50 from the full campaign |
| `results/review-*.csv` | the review runs of the results page (`review.py`) |

The upstream repository is licensed under GPL-3.0. No upstream code is
used here; `conditions.csv` and `reference-variance.csv` are derived from
its data files. They are not covered by this repository's MIT license and
stay under the terms of their source; see
[THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md#reference-data-of-the-guide-case-studies). `fetch.py` regenerates both from the pinned revision.
