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
interval, so it can only widen it. A fit through all points without a
positive slope, or whose inverse lies outside the range of the fit points,
ends the calibration as failed; at stage 1 such a fit puts stage 2 at the
geometric mean of the stage-1 values. (This rule was added after a run of
the small quick-size model extrapolated far beyond its fit points, before
any run of the draft or full size.)

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

**Fallback cohort.** A calibration that ends with a failed fit (at the
first or the second attempt) has no validation cohort and so no count of
types. Such a run still simulates one cohort of 128 fibers, at the fit
point of the last stage 2 whose observed share of affected regions is
closest to 10% (on a tie, the smaller `m`), with its own phase
(`fibers-fallback`) and seeds. Its summary row reads "calibration failed,
fallback cohort" and reports the same measures as a validation cohort; the
statement rule applies unchanged. The calibration stays failed, and the row
is never compatible: a fallback cohort measures the types at whatever
prevalence it happens to have, not at a calibrated 10%. (This was added
before any run of the draft or full size.)

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
years and 27 to 47 at 80 years. A fallback cohort takes the place of a
validation cohort and adds no fibers to these counts. The quick size (4, 8 and 32 fibers of the
small model, at least 3 affected fibers) runs in seconds in pure Python; the
draft size (2, 4 and 16 fibers, at least 3) takes minutes at 3 years.

**Longer-chain control.** Whether 32 regions are long enough is tested
once: 32 fibers of 128 regions at the `m` of the fiber design, 80 years,
the denser mutant (`mutations.py length`, about 12 CPU-hours). If the 95%
intervals of types per affected region of the 32- and the 128-region
cohorts do not overlap, the result reads "chain length matters", otherwise
"no difference seen". If either cohort has fewer than 10 affected fibers,
no statement is made.

### How well the design calibrates: synthetic fibers

Before any simulated run, the fit and the validation rule were tried on
fibers drawn from a known curve, with no simulation
(`mutations.py coverage`, about 3 minutes; targets and rules written
down first). The share of affected regions follows a logistic curve in
`log10(m)` with a slope of 25 per decade and 10% at a true `m`. In the
*whole fibers* scenario a fiber is affected along all 32 regions or not at
all, as at 40 and 80 years. In the *partly clustered* scenario the regions
of a fiber share a random probability (intra-fiber correlation 0.2, a
design effect of about 7, as at 3 and 10 years). Every repetition runs the
design's own code at the full size: 16 and 32 fibers per fit point, 128 per
cohort, 1,999 resamples, a start value up to 7% from the true `m`. 400
repetitions per scenario; the intervals are 95% Wilson intervals of the
shares (`results/review-fiber-fit-coverage.csv`):

<!-- mtdna-spread-fiber-coverage -->
| Experiment, quantity | Whole fibers | Partly clustered |
| --- | --- | --- |
| T4 fit interval: coverage of true m, fiber bootstrap | 0.810 (0.769–0.845), undercovers | 0.907 (0.875–0.932) |
| T4 fit interval: coverage of true m, regions independent | 0.235 (0.196–0.279) | 0.560 (0.511–0.608) |
| T4 fit interval: median width / true m, fiber bootstrap | 0.126 | 0.0498 |
| T4 fit interval: median width / true m, regions independent | 0.0158 | 0.0178 |
| T4 fit interval: failed point fits (no positive slope or inverse outside the fit points) | 43 | 10 |
| T4 fit interval: share of failed bootstrap refits | 0.0193 | 0.0016 |
| T5 validation rule at the true m: fiber interval of prevalence contains 0.10 | 0.943 (0.915–0.961) | 0.935 (0.906–0.955) |
| T5 validation rule at the true m: at least 10 affected fibers | 0.838 (0.798–0.870) | 1.000 (0.990–1.000) |
| T5 validation rule at the true m: compatible | 0.820 (0.779–0.855) | 0.935 (0.906–0.955) |
| T6 end to end, R0: done at attempt 1 | 0.660 (0.612–0.705) | 0.858 (0.820–0.888) |
| T6 end to end, R0: done at attempt 2 | 0.180 (0.145–0.221) | 0.117 (0.090–0.153) |
| T6 end to end, R0: failed calibrations | 0.160 (0.127–0.199) | 0.025 (0.014–0.045) |
| T6 end to end, R0: coverage of true m, fiber bootstrap of the last fit | 0.848 (0.809–0.879) | 0.948 (0.921–0.965) |
| T6 end to end, R0: median error of m in done calibrations (%) | 1.52 | 0.79 |
| T6 end to end, R1: done at attempt 1 | 0.660 (0.612–0.705) | 0.858 (0.820–0.888) |
| T6 end to end, R1: done at attempt 2 | 0.212 (0.175–0.255) | 0.120 (0.092–0.156) |
| T6 end to end, R1: failed calibrations | 0.128 (0.098–0.164) | 0.022 (0.012–0.042) |
| T6 end to end, R1: coverage of true m, fiber bootstrap of the last fit | 0.897 (0.864–0.924) | 0.950 (0.924–0.967) |
| T6 end to end, R1: median error of m in done calibrations (%) | 1.52 | 0.788 |
| T6 decision: R1 replaces R0 (last-fit coverage >= 0.90, no more failures) | 0, R0 kept | – |
| T7 fallback cohort, R0: fit-failed calibrations with a fallback cohort | 40 | 1 |
| T7 fallback cohort, R0: true prevalence at the fallback point, median | 0.0893 | 0.095 |
| T7 fallback cohort, R0: true prevalence at the fallback point, 2.5% quantile | 0.0436 | 0.095 |
| T7 fallback cohort, R0: true prevalence at the fallback point, 97.5% quantile | 0.218 | 0.095 |
| T7 fallback cohort, R0: share of fallback cohorts with at least 10 affected fibers | 0.750 (0.598–0.858) | 1.000 (0.207–1.000) |
<!-- /mtdna-spread-fiber-coverage -->

**The fit interval (T4).** When whole fibers are affected or not, the
interval of the first design's regression, which treats regions as
independent, contains the true `m` in only 23.5% of the repetitions; its
median width is 1.6% of `m`. The bootstrap over fibers is eight times
wider and contains the true `m` in 81%. That is below the 90% the rule
asked for, so by the rule declared before the run it *undercovers*. Most of
the gap comes from the 43 repetitions (11%) whose fit through all points
fails, counted as not covered: in 39 of them the inverse lies just outside
the range of the fit points. Among the other 357, the interval covers in
324 (91%); this split was not declared in advance. With partial clustering
the bootstrap covers in 91%, the regression in 56%.

**The validation rule at the true `m` (T5).** With whole fibers, a cohort
of 128 fibers at exactly the right `m` is compatible in 82% of the
repetitions: its prevalence interval contains 10% in 94%, but in 16% fewer
than 10 fibers hold an affected region (12.8 are expected). So about one
cohort in five triggers the refit even at the right `m`.

**End to end (T6).** The whole procedure, fit, validation and refit, ends
with a compatible cohort in 84% of the repetitions with whole fibers (66%
at the first cohort) and in 97.5% with partial clustering. A done
calibration misses the true `m` by 1.5% (median) with whole fibers and 0.8%
with partial clustering. After the result of T4, one change was tested,
written down before its run: revision R1 lets a failed fit at attempt 1
lead to the refit instead of ending the calibration. R1 lowers the failed
calibrations with whole fibers from 16% to 13% and raises the coverage of
the last fit from 85% to 90% (0.8975). The rule asked for at least 0.90, so
R1 was not adopted, and the full runs use the design as declared (R0). Its
interval of `m` is honest only to about 85% with whole fibers; the
validation cohort, not the fit interval, decides whether `m` is accepted.

**What a fallback cohort is worth (T7).** In the same end-to-end runs of
R0 (rule written down before the run), 40 of the 400 calibrations with
whole fibers end with a failed fit and get a fallback cohort; the other 24
failures have two cohorts that are not compatible. At the fallback points
the true share of affected regions has a median of 8.9% (2.5% to 97.5%
range: 4.4% to 21.8%), and 75% of the fallback cohorts have at least 10
affected fibers, so the statement rule can apply to them. With partial
clustering one calibration ends this way. A fallback cohort therefore
measures the types at a prevalence near, but not at, 10%.

Next: [results and limits](19c-mtdna-spread-de-novo-results.md).
