# Spatial mtDNA spread

Use the [guide](../../docs/guide/19-mtdna-spread.md) for the model and its
connection to Figures 1C and 2E. All times are in days and spatial distances
are in micrometers.

```sh
python examples/mtdna_spread/study.py --size quick
python examples/mtdna_spread/study.py --size draft --threads 4
python examples/mtdna_spread/study.py --figures docs/guide/img/mtdna_spread
```

The draft run repeats each stochastic Figure 1C curve 2,000 times. Its
`terminal-checks.csv` compares the neutral cases with the source and
`two-regions.csv` includes the preferential-elimination cases. The script
builds each model from individual birth, death and migration channels.
At the fixed left boundary, an incoming hop disappears into the
reservoir and outgoing hops sample its fixed copy number.

For the exact source-size run, use
`python examples/mtdna_spread/study.py --size full --fibers 800`.
Each completed fiber is checkpointed. The 500-unit model has a long horizon.
A two-day pilot ran in about 1.2 seconds on one worker. The full cohort
is many CPU-hours. Use `--fibers 1` for a single source-scale run.
The supplement's point-estimate discussion gives 3,000 copies, while the
Figure 2 caption gives 3,500. This example uses the caption.

```sh
python examples/mtdna_spread/wave.py --fibers 48
python examples/mtdna_spread/wave.py --fibers 96 --step 0.125 --output wave-fine
```

The step-size study uses Numba. It computes Poisson births and competing
binomial departures while holding the local feedback fixed over a step.
A three-site mean profile smoother is applied before interpolation.
Three fresh approximate cohorts accompany the archived benchmark.
The 0.125-day cohort gives 0.14617 µm/day at half height, with a
99% interval of 0.13016–0.16360. That interval only partly overlaps
the ±0.03 µm/day equivalence range centered on the archive.
The fresh 800-fiber exact ensemble is still outstanding.

`fetch.py` can regenerate the compact reference CSVs from the pinned
public archive. Run outputs go to `runs/` by default. `--figures`
reads the committed result tables unless `--output` specifies a run.
[Paper and parameters](https://doi.org/10.1073/pnas.2122073119)
