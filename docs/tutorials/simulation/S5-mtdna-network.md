# S5. Mitochondrial networks and the rate of drift

## What this chapter reproduces

A cell holds hundreds to thousands of mtDNA copies. When some copies carry a
mutation, the mutant fraction (heteroplasmy) of a cell drifts at random as
copies are degraded and replicated, and cells that started alike become
different. Mitochondria also fuse into networks and split again. A 2019
paper asks how this network state changes the speed of the drift:

> Aryaman J., Bowles C., Jones N.S., Johnston I.G. (2019). Mitochondrial
> network state scales mtDNA genetic dynamics. *Genetics* 212(4),
> 1429–1443.
> [doi:10.1534/genetics.119.302423](https://doi.org/10.1534/genetics.119.302423)
> (Equations 1–10 and 13, Figure 2E–H, Figure S3D)

Its answer, Equation 13, is that the variance of the mutant fraction between
cells grows at a rate proportional to the fraction of copies that are not
fused. Fused copies are protected from degradation, so they do not take part
in the turnover that drives the drift.

The authors published their simulator, parameter tables and the variance
summaries of their simulations
([MitoNetworksGenetics](https://github.com/ImperialCollegeLondon/MitoNetworksGenetics/tree/8f158fb7a61eaaa622d6e715d654c8bdfb44db89),
revision `8f158fb`). This chapter takes the parameter tables and the
summaries from there, writes the model in CRN-craft, and simulates it
afresh. No code was copied or translated from the authors' repository; it
is licensed under GPL-3.0 and was read for facts only. The paper's authors
did not write or review this chapter.

The result in three lines:

- In the three conditions chosen before the runs (E:0, F:0, H:60), the
  measured early slope divided by the Equation 13 slope is 0.986, 0.986 and
  1.010; all three 99% simultaneous intervals lie inside the declared band
  0.80–1.20 (numerical agreement). In all 39 source conditions, 38 intervals
  lie inside the band; the 39th (H:65) was rerun with a fresh seed and gave
  1.042 (0.925–1.187), so the first value was sampling noise.
- Over the first two days, this agreement does not test the network
  assumption of Equation 13: at time zero the slope equals Equation 13
  exactly for any fusion and fission rate. The archived 0–500 day
  simulations are the test; there the slowest switching case deviates by
  7%, and a fresh 2,000-path run cannot resolve a deviation of that size.
- Fresh simulations agree with the authors' archived variance series within
  their 99% intervals (slow switching, 0–500 days; nominal, 0–100 days).

## Words used in this chapter

- **Singleton**: a copy in an unfused mitochondrion. **Fused**: a copy in a
  network. Every copy is wild type (`W`) or mutant (`M`), so a cell has four
  counts: `WS`, `MS`, `WF`, `MF`.
- **Heteroplasmy** `h`: the mutant fraction `(MS + MF) / n`, with `n` the
  total copy number.
- **Singleton fraction** `f_s`: `(WS + MS) / n`.
- **Condition**: one row of the authors' parameter tables, named by figure
  panel and row index, for example `H:60` (the nominal parameters).
- **Slope ratio**: the measured growth rate of `Var(h)` divided by the
  Equation 13 rate.

## The model

Time is in days. The reactions and their hazards (rates per day of the
whole cell) are those of the paper's Equations 1–10 and of the authors'
simulator, with fused copies never degraded (the paper's nominal case):

| Process | Reaction | Hazard |
|---|---|---|
| Fusion of two identical singletons | `WS + WS -> WF + WF` | `γ·WS·(WS−1)/2` |
| Fusion of a singleton with a fused copy | `WF + WS -> WF + WF`, `MF + WS -> MF + WF` | `γ·WF·WS`, `γ·MF·WS` |
| Fusion of two different singletons | `WS + MS -> WF + MF` | `γ·WS·MS` |
| Fission | `WF -> WS` | `β·WF` |
| Replication of a singleton | `WS -> WF + WF` | `λ·WS` |
| Replication of a fused copy | `WF -> WF + WF` | `λ·WF` |
| Degradation (mitophagy), singletons only | `WS -> ∅` | `μ·WS` |

The mutant reactions are the same with `W` and `M` exchanged, including
`MS + MS`, `MF + MS` and `WF + MS`; there are 15 reactions in total. The
replication rate per copy has linear feedback on the copy number
(Equation 10):

$$\lambda = \max\{0,\ \mu + b\,[\kappa - (WS + WF) - \delta\,(MS + MF)]\}.$$

`μ` is the mitophagy rate (nominal 0.023 per day, a half-life of 30 days),
`β` the fission rate (33.12 per day), `γ` the fusion rate per pair
(0.0379 per day), `b` the feedback strength (1.24e-5 per day per copy),
`κ` sets the target copy number, and `δ` weights mutants in the feedback.
The model is neutral: mutants and wild type behave alike.

Equation 13 gives the variance after a short time `t`, starting at the
deterministic steady state with copy number `n`, singleton fraction `f_s`
and heteroplasmy `h`:

$$\mathrm{Var}(h) \simeq \frac{2\mu f_s}{n}\,h(1-h)\,t.$$

The paper derives it assuming many copies and fusion and fission much
faster than turnover. For the nominal state (`n` = 1,001, `h` = 0.300,
`f_s` = 0.467) it gives 4.50e-6 per day.

```python
import crn_craft as crn

model = crn.Model("identical_fusion")
model.add_species({"WS": 3, "WF": 0})
model.add_reaction(0.4, ["WS", "WS"], ["WF", "WF"], name="fusion")
law = model.fsp(1e-7, [1e-7], {"WS": 3, "WF": 4}, tolerance=1e-14).metadata
p = dict(zip(map(tuple, law["states"]), law["snapshots"][0]["probabilities"]))
print(round(p[(1, 2)] / 1e-7, 4))
```

```text
1.2
```

Three singletons form three unordered pairs, so the fusion hazard is
`0.4 × 3 = 1.2` per day; CRN-craft's mass action counts a pair of identical
reactants once, as the authors' simulator does.

## What is taken from the authors and what is simulated

| Item | Source |
|---|---|
| Parameters and initial counts of 39 conditions (`conditions.csv`) | the authors' parameter tables, read by `fetch.py` |
| Archived variance series, 10,000 paths each, every 10 days to day 500 (`reference-variance.csv`) | the authors' `online_stats_<index>.csv` files |
| Equation 13 slope | computed here from the initial counts, as the authors' analysis does |
| All slope ratios, intervals and fresh variance series | simulated here with CRN-craft 0.8.0, Static Direct |

The 39 conditions are all rows of Figure 2E (mitophagy rate, 4 rows), 2F
(copy number, 5), 2G (mutant weight `δ`, 9), the fusion-to-fission ratio
`R` at magnitude `M` = 1 (Figure 2H, 11 rows) and the magnitude `M` at
`R` = 1 (Figure S3D, 11 rows); `H:60` belongs to both sweeps.

## The protocol

Declared in `manifest.json` before the runs:

- Primary conditions E:0, F:0 and H:60; 2,000 independent paths each, seed
  `200521 +` row number.
- Report times 0, 0.5, 1, 1.5 and 2 days. The slope is the variance of `h`
  over paths regressed on time through the origin.
- Interval: whole-path bootstrap of the slope ratio, two-sided 99% with a
  Bonferroni split over the three conditions (tails 0.01/6).
- Decision: agreement if all three intervals lie inside 0.80–1.20.
- The other 36 conditions are descriptive. For the four fastest switching
  magnitudes (`β` ≥ 100 per day, H:62–H:65) the window is shortened by
  `100/β` and the cohort is 1,024 paths, to bound the event count.

[The results page](S5a-mtdna-network-results.md) shows the outcome, the
review runs, the reconstruction choices and the limits.

## Running the example

You need CRN-craft 0.8.0, NumPy, pandas and Matplotlib. Run from the
repository root:

```sh
python examples/mtdna_network/study.py --size quick
python examples/mtdna_network/study.py --size full
python examples/mtdna_network/review.py --size full
python examples/mtdna_network/study.py --figures docs/guide/img/mtdna_network
python examples/mtdna_network/review.py --figures docs/guide/img/mtdna_network
```

`--size quick` runs 256 paths of H:60 over 0.25 days in about one second.
`--size draft` runs the three primary and four more conditions (about three
minutes on one core). `--size full` runs all 39 conditions with
checkpoints; it was not timed for this page, and its event count suggests
20–30 minutes on one core. `review.py --size full` repeats the review runs
of the results page (about 11 minutes on one core). The
`--figures` commands redraw the plots from the saved CSV files. The
[example README](../../../examples/mtdna_network/README.md) lists the files.

[Tutorial S6](S6-mtdna-bottleneck.md) follows variance across scheduled
cell divisions.
