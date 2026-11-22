# De novo deletions: results and limits

Part of [tutorial S4](S4-mtdna-spread.md); the model and the protocol are on
[page S4b](S4b-mtdna-spread-de-novo.md).

**This page is an extension made in this guide. It is not a result of the
paper, and no author of any of the papers cited here wrote or reviewed it.**

## Results

The validation cohorts, 32 fibers of 32 regions each, at 1,080 days (3
years) and, for the denser mutant only, at 3,600, 14,400 and 28,800 days
(10, 40 and 80 years). The 40- and 80-year rows come from a second machine
([below](#forty-and-eighty-years)); at 40 years the calibration failed.
Tutorial S3's rows are its 20,000-cell validation cohorts at the same ages
of 360 days, for comparison of the layout and the order of magnitude only (the
table is generated from `results/de-novo-summary.csv` and tutorial S3's
`comparison.csv`, and checked by `validation/check_mtdna_spread.py`):

<!-- mtdna-spread-de-novo -->
| Age | Mechanism and unit | m per replication | Affected / units | Prevalence (95% Wilson) | Mutant types (95% bootstrap) | Target band |
| ---: | --- | --- | ---: | --- | --- | --- |
| 3 y | Survival of the densest, 32 regions (this guide), first fit | 0.0113 | 232 / 1024 regions | 22.66% (20.20–25.32%) | 1487 (1445–1516) | outside |
| 3 y | Survival of the densest, 32 regions (this guide), refit | 0.01117 | 82 / 1024 regions | 8.01% (6.50–9.83%) | 1456 (1418–1487) | compatible |
| 3 y | Equal density, 32 regions (this guide), first fit | 0.01133 | 146 / 1024 regions | 14.26% (12.25–16.53%) | 798 (778–813) | outside |
| 3 y | Equal density, 32 regions (this guide), refit | 0.01133 | 114 / 1024 regions | 11.13% (9.35–13.21%) | 816 (795–832) | compatible |
| 3 y | Survival of the densest, isolated regions (this guide), first fit | 0.009855 | 104 / 1024 regions | 10.16% (8.45–12.16%) | 271 (267–275) | compatible |
| 10 y | Survival of the densest, 32 regions (this guide), first fit | 0.002966 | 37 / 1024 regions | 3.61% (2.63–4.94%) | 849 (763–866) | outside |
| 10 y | Survival of the densest, 32 regions (this guide), refit | 0.003 | 83 / 1024 regions | 8.11% (6.59–9.94%) | 817 (767–856) | compatible |
| 40 y | Survival of the densest, 32 regions (this guide), first fit | 0.0004732 | 0 / 1024 regions | 0.00% (0.00–0.37%) | – | outside |
| 40 y | Survival of the densest, 32 regions (this guide), refit | 0.0004663 | 32 / 1024 regions | 3.12% (2.22–4.38%) | 266 (one fiber, no interval) | outside |
| 80 y | Survival of the densest, 32 regions (this guide), first fit | 0.0001375 | 152 / 1024 regions | 14.84% (12.80–17.15%) | 114 (106–121) | outside |
| 80 y | Survival of the densest, 32 regions (this guide), refit | 0.000134 | 111 / 1024 regions | 10.84% (9.08–12.89%) | 108 (95–119) | compatible |
| 3 y | Neutral drift, one cell (tutorial S3) | 0.00758 | 2122 / 20000 cells | 10.61% (10.19–11.04%) | 34.72 (34.51–34.93) | compatible |
| 3 y | Faster copying of shorter mtDNA, one cell (tutorial S3) | 0.00608 | 1994 / 20000 cells | 9.97% (9.56–10.39%) | 29.30 (29.11–29.51) | within |
| 3 y | Transcription-linked replication, one cell (tutorial S3) | 5.1e-06 | 1976 / 20000 cells | 9.88% (9.47–10.30%) | 1.26 (1.24–1.28) | compatible |
| 10 y | Neutral drift, one cell (tutorial S3) | 0.00142 | 1953 / 20000 cells | 9.77% (9.36–10.18%) | 9.46 (9.35–9.57) | compatible |
| 10 y | Faster copying of shorter mtDNA, one cell (tutorial S3) | 0.000599 | 1831 / 20000 cells | 9.16% (8.76–9.56%) | 4.92 (4.84–5.01) | compatible |
| 10 y | Transcription-linked replication, one cell (tutorial S3) | 1.42e-06 | 2049 / 20000 cells | 10.24% (9.83–10.67%) | 1.07 (1.06–1.08) | compatible |
| 40 y | Neutral drift, one cell (tutorial S3) | 0.000171 | 1945 / 20000 cells | 9.72% (9.32–10.14%) | 2.21 (2.16–2.25) | compatible |
| 40 y | Faster copying of shorter mtDNA, one cell (tutorial S3) | 3.25e-05 | 2125 / 20000 cells | 10.62% (10.21–11.06%) | 1.26 (1.24–1.28) | compatible |
| 40 y | Transcription-linked replication, one cell (tutorial S3) | 3.23e-07 | 1966 / 20000 cells | 9.83% (9.43–10.25%) | 1.01 (1.01–1.01) | compatible |
| 80 y | Neutral drift, one cell (tutorial S3) | 6.78e-05 | 1935 / 20000 cells | 9.68% (9.27–10.09%) | 1.50 (1.47–1.53) | compatible |
| 80 y | Faster copying of shorter mtDNA, one cell (tutorial S3) | 1.13e-05 | 2009 / 20000 cells | 10.04% (9.64–10.47%) | 1.10 (1.08–1.11) | within |
| 80 y | Transcription-linked replication, one cell (tutorial S3) | 1.57e-07 | 1985 / 20000 cells | 9.93% (9.52–10.35%) | 1.01 (1.00–1.01) | within |
<!-- /mtdna-spread-de-novo -->

![A: share of affected regions against m in the 3-year fit runs, with the fitted logistic curves and the 10% target. B: mean mutant types per affected region against age, both on log scales, in the last validation cohort of each case and lifespan with 95% bootstrap intervals over fibers: the denser mutant at 3, 10, 40 and 80 years, the 40-year point open because it rests on one fiber of a failed calibration; equal density and isolated regions at 3 years; tutorial S3's three mechanisms in one cell at the same four ages](../../guide/img/mtdna_spread/de-novo-types.png)

**The fit.** The share of affected regions rises from 0 to 1 over a narrow
range of `m`: for the denser mutant from 0 of 128 regions at `m` = 0.0102
to 128 of 128 at 0.0127. The regression's 95% interval for the refit `m`
is 0.01112–0.01121, but it treats regions as independent and is too narrow
for the same reason as the Wilson interval below. The rise is steep
because the mutant fraction of a region is set mostly by mutation
pressure, which is nearly the same everywhere. Without selection, a wild-type copy is replaced by a new
deletion at rate `m * mu`, so after time `t` the expected mutant fraction
is about `1 - exp(-m * mu * t)`. With equal density at `m` = 0.01133 this
gives 0.575, and the refit cohort's mean is 0.576. With the denser mutant
at `m` = 0.01117 it gives 0.570, and the cohort's mean is 0.581: survival
of the densest raises the mutant fraction by about 0.01 in 3 years. The
regions scatter around the mean with a standard deviation of about 0.011
within a fiber (0.014 with equal density), and the fiber means scatter by
about 0.008 (0.014). 10% of regions above 0.6 is a tail of this narrow
distribution.

**The 10% target.** At 3 years, the first validation cohorts of the two coupled cases
were outside the band (22.7% and 14.3% affected), and the predeclared
refit was used for both. The refit cohorts are compatible with the band
(8.0% and 11.1%), and so is the isolated case (10.2%). None lies within
it. The Wilson interval treats the 1,024 regions as independent; the
bootstrap over fibers, below, is two to three times wider, and it is the
honest one: with 32 fibers the prevalence of a cohort cannot be pinned to
±0.5 percentage points. The mean number of types hardly moves between
first fit and refit (1,487 and 1,456; 798 and 816), because it depends on
`m` smoothly, not through the 60% threshold.

**The number of types.** At 3 years, every case gives many types per affected region:
about 1,460 for the denser mutant, about 820 with equal density and about
270 in isolated regions, against 35 for neutral drift and 1.26 for
transcription-linked replication in one cell of tutorial S3. By the
predeclared rule all three read "many types, like drift". No affected
region held fewer than 230 types. The largest single type held 0.4% of
the copies of an affected region with the denser mutant, 0.6% with equal
density and 3.8% in isolated regions; in tutorial S2's drift cells it is
about 16%.

**What density contributes.** The denser mutant gives 1.78 times as many
types per affected region as equal density (95% bootstrap 1.73–1.84), but
the cause is the copy number, not clonal expansion. An affected region of
the denser mutant holds about 6,760 copies instead of 3,500, and the
number of mutant copies per type is about the same in both cases (2.85 and
2.64). Within 3 years, survival of the densest does not make any one type
large: it acts on the mutant as a whole, and all mutant types share its
density.

**What hopping contributes.** Isolated regions need a lower `m` (0.0099
instead of 0.0112), because without exchange between neighbors the regions
scatter more (standard deviation 0.056 instead of 0.011 within a fiber),
and they hold fewer, larger types (16 copies per type). Hopping spreads every lineage
over about 15 regions in 3 years, which mixes many small lineages into
each region.

**Zones and segments.**

<!-- mtdna-spread-de-novo-units -->
| Case, cohort | Fiber 95% interval of prevalence | Fibers with affected regions (all 32 affected) | Mutant copies per type | Largest type, share of copies | Zones | Types per zone | Types per 32-region segment |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| Denser mutant, 3 y, first fit | 14.7–31.4% | 26 (0) | 2.83 (2.78–2.90) | 0.39% (0.37–0.41%) | 73 | 2174 (1980–2391) | 9095 (9034–9161) |
| Denser mutant, 3 y, refit | 4.6–11.9% | 19 (0) | 2.85 (2.77–2.93) | 0.41% (0.38–0.44%) | 49 | 1741 (1601–1871) | 8937 (8882–8990) |
| Equal density, 3 y, first fit | 8.3–20.3% | 21 (0) | 2.71 (2.66–2.78) | 0.62% (0.59–0.67%) | 63 | 1058 (944–1190) | 5061 (5037–5084) |
| Equal density, 3 y, refit | 6.3–17.0% | 19 (0) | 2.64 (2.58–2.71) | 0.59% (0.55–0.63%) | 53 | 1056 (949–1172) | 5049 (5028–5071) |
| Isolated, 3 y, first fit | 8.7–11.7% | 32 (0) | 16.32 (16.11–16.50) | 3.77% (3.53–4.02%) | 95 | 297 (282–314) | 7401 (7357–7449) |
| Denser mutant, 10 y, first fit | 0.6–7.5% | 6 (0) | 4.87 (4.76–5.54) | 0.72% (0.68–0.77%) | 16 | 1007 (817–1095) | 3035 (2999–3070) |
| Denser mutant, 10 y, refit | 3.4–13.8% | 13 (0) | 5.13 (4.89–5.42) | 0.82% (0.76–0.88%) | 34 | 990 (888–1101) | 3107 (3077–3139) |
| Denser mutant, 40 y, first fit | 0.0–0.0% | 0 (0) | – | – | 0 | – | 611 (598–623) |
| Denser mutant, 40 y, refit | 0.0–9.4% | 1 (1) | 16.37 (one fiber, no interval) | 4.10% (one fiber, no interval) | 1 | 688 (one fiber, no interval) | 590 (575–606) |
| Denser mutant, 80 y, first fit | 4.9–26.5% | 8 (2) | 38.87 (36.73–40.91) | 10.41% (7.28–14.41%) | 17 | 148 (127–196) | 202 (192–212) |
| Denser mutant, 80 y, refit | 1.9–22.3% | 7 (2) | 44.01 (35.72–48.58) | 6.93% (5.34–9.20%) | 14 | 141 (110–189) | 202 (194–211) |
<!-- /mtdna-spread-de-novo-units -->

The fiber column counts the fibers with at least one affected region and,
in brackets, those affected along all 32 regions. An affected zone holds
more types than one of its regions (about 1,740
for the denser mutant at 3 years), and a whole 32-region segment holds
about 8,900 mutant types. No 3-year segment was affected as a whole (its
pooled mutant fraction stayed below 60%), so the segment rows describe all
segments, not affected ones. In isolated regions a zone is only a run of
neighbors that share no copies.

**Ten years.** For the denser mutant the same protocol at 3,600 days needed
a refit as well (first cohort 3.6%, refit 8.1% affected). The fitted `m`
is 0.0030, about a quarter of the 3-year value, and an affected region
holds about 820 mutant types (95% bootstrap 767–856; no affected region
fewer than 628), with 5.1 mutant copies per type and a largest type of
0.8% of the copies (at most 1.4%). Over 10 years survival of the densest
raises the mean mutant fraction more than over 3: the refit cohort's mean
is 0.568, while mutation pressure alone, `1 - exp(-m * mu * t)`, gives
0.530. That gain is spread over hundreds of types. Tutorial S3's drift cell
holds 9.5 types at 10 years. Each copy now travels about 27 regions, so a
32-region segment mixes almost completely; its 3,100 types are about 4
times the count of one region.

**Chain length.** Eight fibers of 64 regions at the refit `m` (0.01117,
`results/de-novo-length.csv`) gave 79 of 512 regions affected, with 1,558
types per affected region (95% bootstrap over fibers 1,509–1,599) and 2.65
mutant copies per type. The prevalence differs from the 32-region refit
cohort (8.0%), within the scatter between cohorts seen above; the number of
types stays in the same range. The end regions of the 32-region chains had
the same mean mutant fraction as the central ones (0.581 for regions 0 and
31, 0.581 for regions 8 to 23 in the refit cohort of the denser mutant), so
the closed ends do not shape the result.

### Forty and eighty years

The 40- and 80-year runs of the denser mutant were run on a second machine
by a reader of the guide with the command of
[Running the experiment](#running-the-experiment) (`run --size full
--years 40` and `--years 80`, `--case main`, with several workers), from this
repository at or after the commit that added `--workers`. The operating
system and package versions of that machine are not known. The six
tables of each run are appended unchanged to the saved tables. Two of their fibers
were simulated again on the review host (Linux x86-64, Python 3.13, NumPy
2.3.5, Numba 0.68.1) and agree bit for bit: fiber 28 of the second 80-year
validation cohort (15 affected regions in 6 zones) and fiber 21 of the
second 40-year validation cohort, the one with affected regions. Their 32
region rows (wild type, mutants, types, largest type), their zones and
their segment rows are identical. The 40-year `m` is stored with six
digits only; its full value was rebuilt from the run's probe fiber,
simulated again here, and the saved fit counts, and every fit value of `m`
it implies equals the saved one to the six digits stored. The other fibers were not rerun.

**80 years.** The first cohort was outside the band (152 of 1,024 regions,
14.8%). The refit is compatible with it: at `m` = 1.34 · 10⁻⁴, 111 of 1,024
regions are affected (10.8%, Wilson 9.1–12.9%). An affected region holds
108 mutant types (95% bootstrap over fibers 95–119; no affected region
fewer than 81), with 44 mutant copies per type and a largest type of 6.9%
of its copies (at most 11.6%). By the predeclared rule this reads "many
types, like drift". Tutorial S3's drift cell holds 1.50 types at 80 years,
the transcription cell 1.01. Four of the 32 segments have a pooled mutant
fraction of at least 60%, and a segment holds about 200 types.

**40 years: the calibration failed.** Both cohorts lie outside the band
after the one allowed refit: 0 of 1,024 regions are affected at `m` =
4.73 · 10⁻⁴, and 32 of 1,024 (3.1%, Wilson 2.2–4.4%) at 4.66 · 10⁻⁴. By the
predeclared rule the 40-year run does not reach the target, so it says
nothing about regions at a prevalence of 10%. The 32 affected regions of
the second cohort are one whole fiber. They hold 266 types each (214 to
297), 16 mutant copies per type, and a largest type of 4.1% of their
copies. With one fiber there is no interval: a bootstrap resample either
contains that fiber or has no affected region. This guide therefore makes
no statement on the number of types at 40 years. The saved summary row
reads "many types" only because the run applies the rule to that
degenerate interval.

**Why the calibration breaks down.** At these ages a fiber is affected
more or less as a whole. Within a fiber the mutant fractions of the regions
scatter with a standard deviation of 0.011 to 0.015 at every age, while the
fiber means scatter by 0.008 at 3 years, 0.019 at 10, 0.045 at 40 and
0.075 at 80 years (refit cohorts). The copies of a region spread over the
whole chain ([page S4b](S4b-mtdna-spread-de-novo.md#the-model)), and a
front of the denser mutant, once formed, travels 1.9 to 2.7 mm in 40 years
at the 0.13 to 0.19 µm/day of [tutorial S4](S4-mtdna-spread.md), more than
the 960 µm of a fiber. Of the 32 fibers of the 80-year refit, 2 are
affected along all 32 regions, 5 in part (1, 1, 2, 15 and 28 regions) and
25 not at all; no 3- or 10-year cohort has a fiber affected along its whole
length. The prevalence of a cohort is then close to a count of fibers. Its
variance over fibers is 26 to 33 times the binomial variance of 1,024
independent regions at 40 and 80 years, against 5 to 11 times at 3 and 10
years: the effective sample is about 31 to 39 units, the number of fibers,
not regions. So the fiber interval of the 80-year refit (1.9–22.3%) is
five times wider than its Wilson interval. The fit points are not even
ordered in `m`: at 40 years, 77 of 256 regions are affected at `m` =
4.78 · 10⁻⁴, 0 at 4.88 · 10⁻⁴ and 39 at 4.98 · 10⁻⁴; the 80-year refit
points give 44, 4 and 63 of 256 at 1.348, 1.375 and 1.402 · 10⁻⁴. Each
point is 8 fibers, and the logistic regression, which treats the 256
regions as independent, is far too confident: its 95% interval for the
refit `m` is about ±1% at 40 and ±2% at 80 years, no wider than the 2% step
between fit points whose prevalences differ tenfold.

**Expansion on top of mutation pressure.** Mutation pressure alone,
`1 - exp(-m * mu * t)`, gives a mean mutant fraction of 0.375 at 40 years
and 0.237 at 80 years at the refit `m`, against 0.57 at 3 years. The refit
cohorts' means over all regions are 0.504 and 0.513. At these ages a region
reaches 60% only through the expansion of the denser mutant on top of
mutation pressure. The types grow larger with age: 2.9 mutant copies per
type at 3 years, 5.1 at 10, 16 at 40 (one fiber) and 44 at 80. Their
number per affected region falls, from about 1,460 to 820, 270 (one fiber)
and 110, but stays far above one.

**The rate this needs.** [Page 19](S4-mtdna-spread.md#what-this-chapter-reproduces)
quotes the paper's estimate of the de novo rate for human fibers, 4.1 · 10⁻⁸
to 1.6 · 10⁻⁶ per replication. The 80-year `m` that brings 10% of regions to
60% mutants, 1.34 · 10⁻⁴, is about 84 to 3,300 times that (the failed
40-year values, about 4.7 · 10⁻⁴, about 290 to 11,000 times). The
10%-of-regions target forces a rate far above the paper's estimate, so the
many types at 80 years partly reflect the target. With the paper's rate,
affected regions would be much rarer, and this experiment does not show how
many types they would hold.

## What this shows and what it does not

What follows, for this model with these parameters: at 3 years, de novo
deletions bring 10% of regions to at least 60% mutants only at a deletion
probability of about 0.011 per replication, and then an affected region
holds well over a thousand mutant types, none of which exceeds 1% of its
copies. Survival of the densest does not change this: the neutral control
gives as many types per mutant copy. The affected regions of this
experiment look like the drift cells of tutorials S2 and S3, only more so,
and unlike the essentially clonal affected cells found in experiments. In
3 years the mechanism has no time to act on a single new deletion: a new
type starts as one copy among 3,500 to 6,800 and, in about 75 copy
generations, grows only by drift. At 10 years the picture is the same with
about 820 types per affected region. At 80 years an affected region holds
about 110 types of about 44 copies each, none above 12% of its copies:
survival of the densest makes the types larger and lifts the mutant
fraction well above mutation pressure, but no single type takes over. This
needs a deletion probability 84 to 3,300 times the paper's estimate.

What the experiment does not show:

- **Not a test of the paper.** The paper makes no claim about short-lived
  animals at 3 years. Its estimate of the de novo rate concerns human fibers
  over decades, and it is based on the takeover probability of one founder
  copy. The Figure 2 parameters were chosen by the authors to model the
  wave in rhesus monkey fibers, not for a short-lived animal. Other
  parameters, in particular a smaller `Nss` (the paper's Figure 4 shows
  that a founder's takeover probability depends on `Nss`), could change the
  result.
- **A region is not a cell.** Tutorials S2 and S3 classify cells of 1,000
  copies; here a region of 3,500 to 6,800 copies exchanges copies with its
  neighbors. A muscle fiber section classified as COX-negative in an
  experiment may cover one region or many. The zone and segment rows show
  how the count grows with the unit.
- **10, 40 and 80 years only for the denser mutant, and no result at 40
  years.** Tutorial S3's drift cells hold 2.2 and 1.5 mutant types at 40 and
  80 years. The long lifespans show that survival of the densest makes the
  types larger (2.9, 5.1, 16 and 44 mutant copies per type at 3, 10, 40 and
  80 years) while their number stays far above one; whether equal density
  or isolated regions would differ there is not run. At 40 years the
  calibration failed, and at 80 years the refit is compatible with the band
  only by the Wilson interval, which treats regions as independent.
- **The protocol does not suit long lifespans.** At 40 and 80 years whole
  fibers are affected or not, so 32 fibers are about 32 units. A design
  with fibers as the sampling unit, 128 fibers per cohort and a
  128-region control is declared on
  [page S4b](S4b-mtdna-spread-de-novo.md#a-design-for-long-lifespans). Its
  40- and 80-year runs and the control are outstanding; this page has no
  result from them.
- **The rate is far above the paper's estimate.** The 10% target sets `m`;
  at 80 years it is 84 to 3,300 times the paper's estimate for human
  fibers. The number of types at the paper's rate is not known from this
  experiment.
- **All deletion types are alike.** Real deletions differ in size and
  perhaps in density; one with a larger density advantage than the others
  would not be modeled here.
- **The 10% target is met only loosely.** The prevalence scatters between
  cohorts more than the ±0.5-point band, because regions in a fiber are
  correlated. At 3 years the type count barely depends on this (1,456 to
  1,487 types for `m` between 0.01117 and 0.01130); at 80 years the two
  cohorts give 114 and 108.
- **Fewer units than the papers.** 1,024 regions in 32 fibers per cohort,
  instead of 1,000 independent cells per run in the papers and 20,000 in
  tutorial S3. The intervals above account for the clustering by
  resampling fibers.

The parameters, the boundary conditions, the chain length, the unit
definitions, the year of 360 days and the treatment of mutant templates
are choices of this guide, given with their reasons on
[page S4b](S4b-mtdna-spread-de-novo.md#the-model).

## Running the experiment

| Size | Command | What it does | Time on one core |
| --- | --- | --- | --- |
| quick | `python examples/mtdna_spread/mutations.py run --size quick` | the whole protocol on a small model (4 regions of 20 copies, 120 days), pure Python | 3 s |
| sampler | `python examples/mtdna_spread/mutations.py sampler` | the sampler against CRN-craft and the moment law | 30 s with Numba |
| draft | `python examples/mtdna_spread/mutations.py run --size draft` | the Figure 2 parameters, 2 fibers per fit point and 8 validation fibers per case | not run; about 6 to 10 min with Numba (estimate from the full run) |
| full | `python examples/mtdna_spread/mutations.py run --size full --output examples/mtdna_spread/results` | the 3-year runs of this page, with the refits | 28 CPU-minutes with Numba (main 15, equal density 10, isolated 3) |
| full, 10 years | `python examples/mtdna_spread/mutations.py run --size full --years 10 --case main` | the 10-year runs (133 fibers with the refit) | 50 CPU-minutes with Numba |
| full, 40 years | `python examples/mtdna_spread/mutations.py run --size full --years 40 --case main --workers N` | the 40-year runs (133 fibers with the refit); the saved rows were run on a second machine by a reader of the guide | about 3 CPU-hours with Numba (estimate: 84 to 106 s per fiber on the review host) |
| full, 80 years | `python examples/mtdna_spread/mutations.py run --size full --years 80 --case main --workers N` | the 80-year runs (141 fibers with the refit); the saved rows were run on a second machine by a reader of the guide | about 6.5 CPU-hours with Numba (estimate: 164 s for one fiber on the review host) |
| length | `python examples/mtdna_spread/mutations.py length --m 0.01116696641687773 --regions 64` | eight 64-region fibers | 2 min |
| fiber design, quick | `python examples/mtdna_spread/mutations.py run --size quick --design fibers` | the [fiber design](S4b-mtdna-spread-de-novo.md#a-design-for-long-lifespans) on the small model, pure Python | 9 s |
| coverage | `python examples/mtdna_spread/mutations.py coverage --output examples/mtdna_spread/results` | the fiber design on synthetic fibers (page S4b) | 3 min |
| fiber design, 40 years | `python examples/mtdna_spread/mutations.py run --size full --years 40 --case main --design fibers --workers 16 --work DIR --output DIR` | 288 fibers, 512 with the refit; outstanding | 13 to 24 CPU-hours (estimate) |
| fiber design, 80 years | `python examples/mtdna_spread/mutations.py run --size full --years 80 --case main --design fibers --workers 16 --work DIR --output DIR` | 288 fibers, 512 with the refit; outstanding | 27 to 47 CPU-hours (estimate) |
| 128-region control | `python examples/mtdna_spread/mutations.py length --years 80 --regions 128 --fibers 32 --m-from DIR/de-novo-fibers-rates.csv --reference DIR/de-novo-fibers-summary.csv --workers 16 --output DIR2` | 32 fibers of 128 regions after the 80-year fiber run; outstanding | about 12 CPU-hours (estimate) |

The draft and full sizes need Numba in practice: without it the same
sampler runs in pure Python, roughly 100 times slower. Every fiber is
checkpointed under `runs/de-novo/<size>/` (`--work`), so an interrupted run
continues where it stopped. Add `--workers N` to `run` to simulate up to
`N` fibers at the same time in separate processes, for example
`run --size full --years 80 --case main --workers 16`; every fiber has its
own seed, so the tables do not depend on `N`, apart from the measured
`seconds`. The times in the table are
CPU time on one core. `mutations.py figure --figures FOLDER` redraws the
figure from the saved tables; add `--results RUN_FOLDER` to draw it from
the tables of a run of your own (the folder given to `run --output`).

The `run` command writes the tables of one lifespan; the saved tables hold
the 3-, 10-, 40- and 80-year rows in this order. A 10-year fiber has about
185 million events and took 23 s; a 40-year fiber has 730 to 890 million
and an 80-year fiber about 1.4 billion. Outstanding: equal density and
isolated regions at 10, 40 and 80 years (at 10 years about 35 and 15
CPU-minutes, at 40 and 80 years a few CPU-hours each), and the fiber-design
runs at 40 and 80 years with the 128-region control (no results yet).

The cost of the fiber design is an estimate of 0.13 CPU-seconds per region
and year with Numba, from the machine that made the first 40- and 80-year
runs (about 170 and 330 s per 32-region fiber). One 40-year fiber took 83
s on the two-core review host, about half of that. `run --design fibers`
prints the number of fibers and this estimate at its start. With
`--workers 16` the stages run in 18 rounds of 16 fibers without a refit and
32 with it, about 1.7 and 3 hours at 80 years.

## Exercises

1. Run `mutations.py run --size quick` and open `runs/de-novo/quick/de-novo-summary.csv`.
   Then make the regions larger: pass `nss=200` in `SIZES["quick"]["model"]`.
   How do the fitted `m` and the number of types per affected region change?
2. From `results/de-novo-regions.csv`, compute for each case the share of
   the largest type, `largest / (wild + mutant)`, in affected regions. Why
   is it about ten times larger in isolated regions than in the chain,
   although the deletion probability is almost the same?

See also: [tutorial S3b](S3b-mtdna-mechanisms-results.md), [tutorial S2c](S2c-fig1cd.md#result).
