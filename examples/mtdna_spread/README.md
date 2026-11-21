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
99% interval of 0.13016–0.16360. The interval excludes the archived
0.18973 µm/day: an unresolved discrepancy. The fresh 800-fiber exact
ensemble is still outstanding.

`review.py` holds the review experiments on that gap (deterministic part,
copy number, step size, exact-versus-approximate noise, and the
equal-density setting of Figure 1C). Its rows are in
`results/review-experiments.csv`; its working files go to `runs/review`.
The cohorts take 4 to 15 minutes each on one core.

`fetch.py` can regenerate the compact reference CSVs from the pinned
public archive. That archive has no license file. `reference-1c.csv` and
`reference-2e.csv` are numerical extracts of its data; they are not
covered by this repository's MIT license and stay under the authors'
terms; see
[THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md#reference-data-of-the-guide-case-studies). No upstream code is included. Run outputs go to `runs/` by default. `--figures`
reads the committed result tables unless `--output` specifies a run.
[Paper and parameters](https://doi.org/10.1073/pnas.2122073119)

## De novo deletions and mutant types

The paper's model has no de novo mutations. `mutations.py` adds them as an
extension of this guide ([page 19b](../../docs/guide/19b-mtdna-spread-de-novo.md)):
32 regions with the Figure 2 parameters, all wild type at the start, a new
deletion type with probability `m` per replication, and `m` fitted so that
10% of regions hold at least 60% mutants after 3 years of 360 days
(and 10, 40 and 80 years for the denser mutant; the 40- and 80-year tables
were run on a second machine by a reader of the guide with the same
command, and at 40 years the calibration failed). Its
exact per-copy sampler follows every copy's type label; with Numba it runs
about 120 ns per event, without Numba it is far slower and only suits the
quick size.

```sh
python examples/mtdna_spread/mutations.py run --size quick
python examples/mtdna_spread/mutations.py sampler
python examples/mtdna_spread/mutations.py run --size full --output examples/mtdna_spread/results
python examples/mtdna_spread/mutations.py figure --figures docs/guide/img/mtdna_spread
```

`run` checkpoints every fiber under `runs/de-novo/<size>/` (`--work`), so an
interrupted run resumes. `run --workers N` simulates up to `N` fibers at the
same time in separate processes; every fiber has its own seed, so the
tables do not depend on `N`, apart from the measured `seconds`. The quick size uses a small model (4 regions of 20
copies, 120 days) in pure Python and takes a few seconds. `sampler` compares
the sampler with the CRN-craft model of the same reactions (`crn_model`) and
with a moment law. The saved tables are `results/de-novo-*.csv`, with the
3-, 10-, 40- and 80-year rows in this order; `figure` draws
`de-novo-types.png` from them.

### Fiber design for long lifespans

At 40 and 80 years whole fibers are affected or not, so the first design
counts far too many independent units. `run --design fibers` makes the
fiber the sampling unit in the fit and in the validation
([page 19b](../../docs/guide/19b-mtdna-spread-de-novo.md#a-design-for-long-lifespans),
manifest entry `de_novo_fibers`). The default, `--design regions`, is the
first design and writes the same tables as before.

```sh
python examples/mtdna_spread/mutations.py run --size quick --design fibers
python examples/mtdna_spread/mutations.py coverage --output examples/mtdna_spread/results
python examples/mtdna_spread/mutations.py run --size full --years 40 --case main --design fibers --workers 16 --work runs/f40 --output runs/f40
python examples/mtdna_spread/mutations.py run --size full --years 80 --case main --design fibers --workers 16 --work runs/f80 --output runs/f80
python examples/mtdna_spread/mutations.py length --years 80 --regions 128 --fibers 32 --m-from runs/f80/de-novo-fibers-rates.csv --reference runs/f80/de-novo-fibers-summary.csv --workers 16 --work runs/l128 --output runs/l128
```

The run prints how many fibers it will simulate and an estimate of the
CPU-hours (0.13 CPU-seconds per region and year with Numba, from timings on a
2.3 GHz Xeon): 288 fibers without a refit and 512 with it, about 13 to 24
CPU-hours at 40 years and 27 to 47 at 80 years; the 128-region control
about 12. Every fiber is checkpointed in `--work`, so a stopped run
continues where it stopped, and the tables do not depend on `--workers`
apart from the measured `seconds`. It writes `de-novo-fibers-rates.csv`,
`-fit.csv`, `-summary.csv`, `-regions.csv`, `-zones.csv` and
`-segments.csv`, with the columns of the first design plus
`fibers_affected`, `fibers_full`, `m_boot_low95`/`m_boot_high95` (the
bootstrap interval of `m` over fibers), `decision` and `calibration`; `m`
keeps all its digits. In these tables `wilson_low`/`wilson_high` assume
independent regions, and `band` and `m_low95`/`m_high95` are the first
design's rule and regression, kept for comparison only. `length` writes the
region rows to `de-novo-length.csv` and one summary row, in the same
layout, to `de-novo-length-summary.csv`. `coverage` tests the fiber design
on synthetic fibers without simulation (`results/review-fiber-fit-coverage.csv`).
The 40- and 80-year runs of this design and the control are outstanding.
