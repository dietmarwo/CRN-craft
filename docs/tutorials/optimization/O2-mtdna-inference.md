# O2. Case study: inferring mtDNA takeover times

## What this chapter reproduces

A single-cell experiment sees each cell once. It cannot follow how fast a
deletion mutant takes over a cell; it only sees, at each age, how many
cells carry how much mutant mtDNA. A 2026 paper asks whether such
snapshots are enough to estimate the mutation probability and the
selection advantage of the mutants, and then turns the advantage into the
time a mutant needs to take over a cell:

> Kowald A., Kirkwood T.B.L. (2026). Inferring accumulation times of
> mitochondrial DNA deletion mutants from cross-sectional single-cell data:
> methodological framework and validation. *npj Aging* 12.
> [doi:10.1038/s41514-026-00431-4](https://doi.org/10.1038/s41514-026-00431-4)
> (Equation 7, Figures 5 and 6, Tables 2 and 3)

This chapter is a reconstruction from the published text and figures. The
authors' programs and data were not available to it (their repository
returned server errors). The printed values of Tables 1 to 3 and the
readings of Figures 5 and 6 were checked against the preprint of the paper
(bioRxiv [10.64898/2026.01.30.702736](https://doi.org/10.64898/2026.01.30.702736), version 1); the review did not have the journal PDF. The paper's authors did
not write or review this chapter.

The result in three lines:

- With 40 synthetic cohorts of 10,000 cells per regime, the fit recovers
  the mutation probability within 7%, the selection advantage within 20%
  and the fraction of advantageous mutants within 18% of the truth (95%
  percentile intervals); every true value lies inside its interval. The
  fits search a box around the true values, which a real application would
  not know ([what that means](O2b-mtdna-inference-results.md#how-to-read-the-recovery)).
- In the baseline regime the reconstructed intervals are narrower than the
  published ones and the published means are shifted (selection advantage
  +7%, advantageous fraction +9%); in the strong regime both agree more
  closely. The cause is not known. The optimizer budget, the search box
  and the surrogate of the forward law are ruled out.
- At selection advantage 0.10, 9.1% of single founder mutants ever reach
  90% of the copies, and those that do need a median of 1,038 days. Only
  58% of them get there within the three years the data cover.

## Words used in this chapter

- **Snapshot**: the mutant fraction of each cell at one age, without
  knowing the same cell's past or future.
- **Threshold**: the paper summarizes a snapshot by the fraction of cells
  with *more than* 0%, 5%, 20% and 50% mutant copies. A cell with `M`
  mutants among `N` copies is above 5% when `20*M > N`, above 20% when
  `5*M > N` and above 50% when `2*M > N`; a cell exactly at a threshold is
  not above it.
- **Cohort**: one synthetic data set, 10,000 cells observed at the 37 ages
  0, 30, …, 1,080 days. In this reconstruction the same 10,000 cells are
  followed through all ages.
- **Takeover**: a mutant lineage reaching 90% of the copies of its cell
  (900 of 1,000, at least).
- **Regime**: one set of true parameters. *Baseline* is `p = 6.8e-6`,
  `s = 0.10`, `f = 0.80` (Figure 6, Table 2); *strong* is `p = 3e-5`,
  `s = 0.50`, `f = 0.20` (Table 3).

## The model in plain words

A cell holds a fixed number of mtDNA copies, 1,000. Copies are replaced one
at a time: one copy is degraded and one copy is replicated, so the number
stays the same. The copy that replicates is chosen in proportion to its
fitness, `1 + s` for an advantageous mutant and 1 for the others; the copy
that is degraded is chosen at random, independently. When a wild-type copy
replicates, its daughter is a new deletion mutant with probability `p`; a
fraction `f` of new mutants is advantageous, the others are neutral.
Mutants do not mutate again. With a half-life of 10 days, one day is 69.3
such replacements. This is the Moran model of the paper's Equation 7.

The pages of the chapter:

| Page | Question |
| --- | --- |
| [The model and Figure 6](O2a-mtdna-inference-model.md) | What exactly is simulated, how does CRN-craft express it, and how close are the published curves? |
| [Fitting snapshots and takeover times](O2b-mtdna-inference-results.md) | How are the parameters fitted, how well are they recovered, and what does that mean for takeover times? |
| [Choices, checks and limits](O2c-mtdna-inference-limits.md) | What did the reconstruction have to decide, how was it checked, and what does it not show? |

## CRN-craft and the local solver

The Moran step is six reactions that move one copy between the wild type
`W`, the neutral mutant `U` and the advantageous mutant `A`. CRN-craft
runs them with Static Direct and solves their probability law with finite
state projection (FSP); [the model page](O2a-mtdna-inference-model.md)
runs a seven-copy cell both ways and the check compares them with exact
matrix calculations.

The study itself runs on a small C++17 library local to the example
(`examples/kk26/moran.cpp`). It changes nothing in CRN-craft. It exists
for two reasons:

1. The fit needs the full probability law of a 1,000-copy cell, 501,501
   joint count states, at 37 ages, for every parameter set the optimizer
   tries. A specialized lattice solver does this in about three seconds.
2. The paper counts time in Moran steps, a discrete clock, while CRN-craft
   runs in continuous time. The library has both: the exact discrete law
   as reference and synthetic data, and a continuous-time solver for the
   fit. The two agree to 1.2e-5 in probability.

## Running the example

You need Python 3.10 or newer, CRN-craft 0.8.0, fcmaes, SciPy, Matplotlib
and a C++17 compiler (set `CXX` to choose one). The library is built in a
temporary directory the first time a command needs it. Run the commands
from the repository root:

```sh
python -m pip install 'crn-craft==0.8.0' 'fcmaes==2.3.0' scipy matplotlib
python examples/kk26/study.py bootstrap --size quick --threads 2
python examples/kk26/study.py fit examples/kk26/results/reconstruction/snapshots-baseline.csv
python examples/kk26/study.py figures
```

`bootstrap --size quick` simulates and fits four baseline cohorts of 2,000
cells in about ten seconds. `fit` fits the first saved baseline cohort and
prints the result as JSON. `figures` redraws the two figures from the
saved results without simulating. The
[example README](../../../examples/kk26/README.md) lists all commands and
output files; [the results page](O2b-mtdna-inference-results.md#rerunning-the-study)
gives the commands of the full study.

See also: [Re-creating a published study](../simulation/S2-published-study.md) (the
same fixed copy number, with a pair rule),
[Comparing mtDNA mechanisms](../simulation/S3-mtdna-mechanisms.md) (the ATP model that
generates panel A of Figure 6),
[Fitting](O1-fitting.md).

[Tutorial S4](../simulation/S4-mtdna-spread.md) follows mtDNA expansion along a muscle fiber.
