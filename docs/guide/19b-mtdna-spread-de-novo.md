# De novo deletions and the number of mutant types

Part of [chapter 19](19-mtdna-spread.md).

**This page is an extension made in this guide. It is not a result of the
paper, and no author of any of the papers cited here wrote or reviewed it.**
The question was suggested by an author of the papers behind chapters 16
and 17 (personal communication, November 2026).

## The question

The paper of chapter 19 starts every simulation with mutants already
present and follows their spread; it has no de novo mutations
([chapter 19](19-mtdna-spread.md#what-this-chapter-reproduces)).
Chapters [16](16c-fig1cd.md) and [17](17-mtdna-mechanisms.md) rest on an
argument that needs de novo mutations. In a short-lived animal, a mechanism
of clonal expansion must take over an affected cell within the lifespan. If
it needs a high deletion rate to do so, the affected cells hold many
different mutants, whereas the experiments these papers cite find
essentially one deletion type per affected cell. Kowald & Kirkwood (2013, *Aging Cell* 12:728–731,
[doi:10.1111/acel.12098](https://doi.org/10.1111/acel.12098)) and
(2014, *PNAS* 111:2972–2977,
[doi:10.1073/pnas.1314970111](https://doi.org/10.1073/pnas.1314970111))
used this to compare mechanisms: chapter 17b finds about 35 mutant types
per affected cell for neutral drift at 3 years and about 1.26 for
transcription-linked replication. The paper of chapter 19 uses the same
argument against neutral models: their high de novo rates "produce an
unrealistically high mutational diversity, at odds with the observed
clonality" (section "Low Mutation Rates Can Yield Large Mutant Loads"). It
does not simulate the diversity of its own model.

This page asks the same question of stochastic survival of the densest:
with the paper's Figure 2 parameters, no mutants at the start and new
deletions at replication, how many different mutant types does an affected
region hold after 3 years, and after a human lifespan?

## The model

The chain of [chapter 19](19-mtdna-spread.md#the-model) with the Figure 2
parameters: death rate `mu` = 0.07/day, feedback `c` = 0.0002/day, `Nss` =
3,500 copies, `delta` = 0.2, hop rate `gamma` = 0.1/day per copy and
directed edge, regions of 30 µm. One reaction is added:

- **De novo deletion.** At every replication, of a wild-type or a mutant
  template, the daughter is a new deletion type with probability `m`;
  otherwise it copies its template's type. Every new type gets a label never
  used before (infinite alleles). A deletion of a mutant is a new type.
- **All mutant types are alike.** Each has the paper's density parameter
  `delta` = 0.2. The model does not know which part of the genome a
  deletion removes.

Every region starts with 3,500 wild-type copies and no mutant. Time runs
in days; a year has 360 days, as in chapter 16. The protocol named 3 years
(1,080 days), the case of a short-lived animal, as the essential lifespan,
and 10, 40 and 80 years as lifespans to add if the compute budget allowed.
All three cases are run at 3 years; the denser mutant is also run at 10,
40 and 80 years (3,600, 14,400 and 28,800 days).

**Geometry.** A fiber is a chain of 32 regions (960 µm) with closed ends: a
copy in an end region has one neighbor and hops only inward, as in the
authors' archived C code. A copy hops at 0.2 per day in total, so over
1,080 days its position spreads with a standard deviation of about
√(0.2 × 1,080) ≈ 15 regions: within 3 years, the copies of one region mix
with about half the chain. A longer fiber would add regions whose
neighborhoods are like those of the central regions here; it would change
little per region, add more regions per fiber, and raise the number of
types in a whole segment roughly in proportion to its length. A run with 64
regions tests this on the [results page](19c-mtdna-spread-de-novo-results.md#results).
At 40 and 80 years the spread is √(0.2 × 14,400) ≈ 54 and √(0.2 × 28,800)
≈ 76 regions, longer than the chain: the copies of every region mix with
the whole fiber, and the fiber, not the region, becomes the unit that
varies. The [results page](19c-mtdna-spread-de-novo-results.md#forty-and-eighty-years)
shows what this does to the calibration.

**Units.** A **region** is affected when at least 60% of its copies are
mutants (`5*M >= 3*(W+M)`, as in chapter 17); its types are the different
mutant labels among its copies. A **zone** is a maximal run of adjacent
affected regions in one fiber; its types are counted among all copies of
the zone. A **segment** is the whole 32-region fiber; its types are counted
among all its copies. The 10% target applies to regions. A region is not a
cell in the sense of chapters 16 and 17: it is the domain of one nucleus in
a long fiber, it holds 3,500 copies instead of 1,000, and it exchanges
copies with its neighbors.

**Three cases.** *Denser mutant* is the model above. *Equal density* sets
`delta` = 1, the neutral control: same chain, same hopping, no survival of
the densest. *Isolated regions* sets `gamma` = 0: every region is a closed
population of 3,500 copies, the closest analog of a cell in chapters 16
and 17.

## The simulator

Every copy in region `i`, wild type or mutant, of any type, is copied at
the same rate `lambda_i = max(0, mu + c * (Nss - W_i - delta * M_i))`, dies
at `mu` and hops at `gamma` per neighbor. The type only matters for the
counts that set `lambda_i`. So an exact path needs no rate per type: pick a
region with probability proportional to its total rate, pick one of its
copies uniformly, and let that copy be copied, die or hop. `mutations.py`
does exactly this (Gillespie's direct method over regions) and keeps one
label per copy. It stores the labels of a region with the mutants first,
so the run without labels draws the same random numbers and gives the same
counts; the check compares both runs exactly.

CRN-craft can write this model directly: one species `mtDNA` with the
attributes `clone` and `region`, `new_id()` for a new deletion, a hop rule
that changes `region`, and per-region observables with `reduce="distinct"`
for the types (`crn_model` in `mutations.py`). Only the `tree` method runs
species created during a run, and it recomputes every rule after every
event. With 32 regions of 3,500 copies it took 9.6 µs per event on the
two-core host of this study; the example-local sampler takes about 0.12 µs
with Numba. A 3-year fiber has about 56 million events: about 9 minutes on
CRN-craft, 7 seconds in the sampler. The CRN-craft model serves as the
reference for the sampler on a small chain.

### How the sampler was checked

Targets and pass rules were written down before these runs; the rows are
in `results/de-novo-sampler.csv`.

- **Against CRN-craft.** A small chain (3 regions, `Nss` = 20, `c` =
  0.01/day, `m` = 0.02, 60 days, the other rates as above): 20,000 paths
  from the sampler and 20,000 from the CRN-craft model on `tree`. The
  largest gap between the two empirical distribution functions must stay
  below the two-sample Dvoretzky–Kiefer–Wolfowitz bound at error
  probability 0.001, 0.0288, as in
  [chapter 17c](17c-mtdna-mechanisms-limits.md#how-the-reconstruction-was-checked).
  Observed: 0.0068 (wild type in region 0), 0.0042 and 0.0047 (mutants in
  regions 0 and 1), 0.0041 (types in region 0) and 0.0023 (types in the
  fiber). Types in the fiber are counted from CRN-craft's final instances:
  its `distinct` observable counts instances, that is pairs of type and
  region.
- **A moment law.** With `c` = 0 every copy is copied at exactly `mu`, so a
  copy is replaced by a new type at rate `m * mu` and the expected mutant
  count is `L * Nss * (1 - exp(-m * mu * t))`; the expected number of
  deletions is `m * mu * L * Nss * t`. 4,000 paths of 4 regions of 50
  copies (`m` = 0.05, 100 days): z = −0.41 and −0.01.
- **Exact laws in the check.** Labels on and off give identical counts and
  event totals; each region holds one label per mutant copy, every label
  is a type that was made; births minus deaths equal the change of the
  copy number; hops alone change nothing.

The sampler has no time step, so the stepping error of `wave.py` does not
arise here.

## The protocol

Written down before any calibration run (and recorded in the example's
`manifest.json`, entry `de_novo`), in the order of [chapter 17b](17b-mtdna-mechanisms-results.md#the-protocol),
with fewer units because a fiber costs seconds rather than a millisecond:

- **Target**: 10% of regions affected at the end of the lifespan (1,080
  days for 3 years). A prevalence interval
  is *within* the band if it lies inside 9.5–10.5%, *compatible* if it
  intersects it, and *outside* otherwise.
- **Fit**: a start value of `m` from one probe fiber, then 4 fibers (128
  regions) at 0.9, 0.97, 1.04 and 1.12 times it, extended until the points
  enclose 10%; then 8 fibers (256 regions) at 1/1.02, 1 and 1.02 times the
  logistic inverse at 10%. The fitted `m` inverts a logistic regression of
  the affected indicator on `log10(m)` through all fit regions. Only the
  affected indicator enters the fit, never a type count.
- **Validation**: 32 new fibers (1,024 regions) at the frozen `m`, with
  their own seeds (a hash of the manifest seed, case, size, lifespan, phase,
  `m` and fiber index). Prevalence gets a 95% Wilson interval, which treats regions as
  independent, and a 95% percentile bootstrap interval over whole fibers
  (1,999 resamples), which does not. Means per region, zone and segment get
  the fiber bootstrap.
- **Refit rule**: a validation outside the band is kept. One refit is
  allowed: stage 2 again with new seeds around the first `m`, a fit through
  all fit points, and a second validation cohort. Both are reported.
- **Statement rule**: a 95% interval of the mean types per affected region
  below 2 reads "about one type, like the transcription mechanism"; above
  10 reads "many types, like drift"; otherwise "intermediate". Fewer than
  20 affected regions: no statement.

Regions of one fiber share copies, so they are far from independent. The
Wilson interval is therefore too narrow, and with only 32 fibers the
prevalence of a cohort scatters more than the band allows; the bootstrap
over fibers shows by how much. The same holds for the fit, whose logistic
regression also treats regions as independent. The longer the lifespan, the
stronger the effect: at 40 and 80 years whole fibers are affected or not,
and the protocol, written for 3 years, no longer pins the prevalence.

## A design for long lifespans

At 40 and 80 years the first protocol has about as many independent units
as fibers, 32 per cohort, while its fit and its band count 1,024 regions
([results page](19c-mtdna-spread-de-novo-results.md#forty-and-eighty-years)).
A second design, `mutations.py run --design fibers`, makes the fiber the
sampling unit everywhere. It was written down before any run of it, in
`manifest.json` (entry `de_novo_fibers`) and here. The first design, its
code path (the default, `--design regions`) and its saved tables stay as
they are.

**What stays the same.** The target is still 10% of regions with at least
60% mutants at the end of the lifespan, so the rows of both designs can be
compared. The model, the parameters, the 32-region chain with closed ends,
the three cases, the year of 360 days, the units, the statement rule and
the seed rule are those above.

**Fit.** Stage 1 simulates 16 fibers at each of 0.92, 0.97, 1.03 and 1.09
times a start value of `m`. At a logistic slope of 20 to 30 per decade of
`m` this spans about 5% to 19% affected regions. The grid is extended by
factors of 1.25, as above, until it encloses 10%. Stage 2 simulates 32
fibers at 1/1.03, 1 and 1.03 times the stage-1 estimate. The start value
is the last fitted `m` of the first design for the same case and lifespan
(`results/de-novo-rates.csv`) where there is one, and otherwise the probe
fiber of the first design. The point estimate is the same logistic
regression of the region indicator on `log10(m)`. Its 95% interval comes
from a bootstrap over fibers: within each fit point, whole fibers are
resampled, the curve is refitted and inverted at 10%, 1,999 times. A
resample whose fit has no positive slope counts at both ends of the
interval, so it can only widen it.

**Validation.** 128 new fibers at the frozen `m`. The cohort is
*compatible* if the 95% bootstrap interval over fibers of the share of
affected regions contains 10% and at least 10 fibers hold an affected
region. Otherwise it is *outside* (the interval excludes 10%) or has *too
few affected fibers*. There is no band of ±0.5 points: with 128 fibers that
behave as single units, the prevalence of a cohort cannot be pinned that
closely. The Wilson interval is still reported, and it assumes independent
regions. The statement on the number of types needs at least 20 affected
regions in at least 10 fibers; otherwise none is made. This excludes the
one-fiber case of the first design at 40 years.

**Refit.** One refit is allowed if the first cohort is not compatible:
three new stage-2 points of 32 fibers around the `m` fitted through all
points, a fit through all points, and a new cohort of 128 fibers. If that
one is not compatible either, the calibration failed. All rows are kept.

**Seeds and tables.** The phases have their own names (`fibers-stage1`,
`fibers-validation` and so on), so no fiber of the first design is reused
and no checkpoint is shared. The tables are `de-novo-fibers-*.csv`, with
the columns of the first design and, in addition, the fibers with an
affected region (`fibers_affected`), the fibers affected along their whole
length (`fibers_full`), the bootstrap interval of `m` over fibers
(`m_boot_low95`, `m_boot_high95`) and the label of the cohort
(`decision`). The columns `m_low95`, `m_high95` (regression, regions
independent) and `band` (the first design's rule) are kept for comparison
only.

**Size and cost.** Without a refit the design simulates 4 × 16 + 3 × 32 +
128 = 288 fibers, with one 512. A 32-region fiber costs about 170
CPU-seconds at 40 years and 330 at 80 years with Numba on the machine that
made the first 40- and 80-year runs. That is about 13 to 24 CPU-hours at 40
years and 26 to 47 at 80 years. The quick size (4, 8 and 32 fibers of the
small model, at least 3 affected fibers) runs in seconds in pure Python; the
draft size (2, 4 and 16 fibers, at least 3) takes minutes at 3 years.

**Longer-chain control.** Whether 32 regions are long enough is tested
once: 32 fibers of 128 regions at the `m` of the fiber design, 80 years,
the denser mutant (`mutations.py length`, about 12 CPU-hours). If the 95%
intervals of types per affected region of the 32- and the 128-region
cohorts do not overlap, the result reads "chain length matters", otherwise
"no difference seen".

Next: [results and limits](19c-mtdna-spread-de-novo-results.md).
