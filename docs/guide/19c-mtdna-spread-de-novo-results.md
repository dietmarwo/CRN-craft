# De novo deletions: results and limits

Part of [chapter 19](19-mtdna-spread.md); the model and the protocol are on
[page 19b](19b-mtdna-spread-de-novo.md).

**This page is an extension made in this guide. It is not a result of the
paper, and no author of any of the papers cited here wrote or reviewed it.**

## Results

The validation cohorts, 32 fibers of 32 regions each, at 1,080 days (3
years) and, for the denser mutant only, at 3,600 days (10 years).
Chapter 17's rows are its 20,000-cell validation cohorts at 3 and 10 years
of 360 days, for comparison of the layout and the order of magnitude only (the
table is generated from `results/de-novo-summary.csv` and chapter 17's
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
| 3 y | Neutral drift, one cell (chapter 17) | 0.00758 | 2122 / 20000 cells | 10.61% (10.19–11.04%) | 34.72 (34.51–34.93) | compatible |
| 3 y | Faster copying of shorter mtDNA, one cell (chapter 17) | 0.00608 | 1994 / 20000 cells | 9.97% (9.56–10.39%) | 29.30 (29.11–29.51) | within |
| 3 y | Transcription-linked replication, one cell (chapter 17) | 5.1e-06 | 1976 / 20000 cells | 9.88% (9.47–10.30%) | 1.26 (1.24–1.28) | compatible |
| 10 y | Neutral drift, one cell (chapter 17) | 0.00142 | 1953 / 20000 cells | 9.77% (9.36–10.18%) | 9.46 (9.35–9.57) | compatible |
| 10 y | Faster copying of shorter mtDNA, one cell (chapter 17) | 0.000599 | 1831 / 20000 cells | 9.16% (8.76–9.56%) | 4.92 (4.84–5.01) | compatible |
| 10 y | Transcription-linked replication, one cell (chapter 17) | 1.42e-06 | 2049 / 20000 cells | 10.24% (9.83–10.67%) | 1.07 (1.06–1.08) | compatible |
<!-- /mtdna-spread-de-novo -->

![A: share of affected regions against m in the fit runs, with the fitted logistic curves and the 10% target. B: mutant types per affected region in the last validation cohort of each case, on a log scale, with chapter 17's three-year values for neutral drift and transcription-linked replication](img/mtdna_spread/de-novo-types.png)

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

**The 10% target.** The first validation cohorts of the two coupled cases
were outside the band (22.7% and 14.3% affected), and the predeclared
refit was used for both. The refit cohorts are compatible with the band
(8.0% and 11.1%), and so is the isolated case (10.2%). None lies within
it. The Wilson interval treats the 1,024 regions as independent; the
bootstrap over fibers, below, is two to three times wider, and it is the
honest one: with 32 fibers the prevalence of a cohort cannot be pinned to
±0.5 percentage points. The mean number of types hardly moves between
first fit and refit (1,487 and 1,456; 798 and 816), because it depends on
`m` smoothly, not through the 60% threshold.

**The number of types.** Every case gives many types per affected region:
about 1,460 for the denser mutant, about 820 with equal density and about
270 in isolated regions, against 35 for neutral drift and 1.26 for
transcription-linked replication in one cell of chapter 17. By the
predeclared rule all three read "many types, like drift". No affected
region held fewer than 230 types. The largest single type held 0.4% of
the copies of an affected region with the denser mutant, 0.6% with equal
density and 3.8% in isolated regions; in chapter 16's drift cells it is
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
| Case, cohort | Fiber 95% interval of prevalence | Mutant copies per type | Largest type, share of copies | Zones | Types per zone | Types per 32-region segment |
| --- | --- | ---: | ---: | ---: | --- | --- |
| Denser mutant, 3 y, first fit | 14.7–31.4% | 2.83 (2.78–2.90) | 0.39% (0.37–0.41%) | 73 | 2174 (1980–2391) | 9095 (9034–9161) |
| Denser mutant, 3 y, refit | 4.6–11.9% | 2.85 (2.77–2.93) | 0.41% (0.38–0.44%) | 49 | 1741 (1601–1871) | 8937 (8882–8990) |
| Equal density, 3 y, first fit | 8.3–20.3% | 2.71 (2.66–2.78) | 0.62% (0.59–0.67%) | 63 | 1058 (944–1190) | 5061 (5037–5084) |
| Equal density, 3 y, refit | 6.3–17.0% | 2.64 (2.58–2.71) | 0.59% (0.55–0.63%) | 53 | 1056 (949–1172) | 5049 (5028–5071) |
| Isolated, 3 y, first fit | 8.7–11.7% | 16.32 (16.11–16.50) | 3.77% (3.53–4.02%) | 95 | 297 (282–314) | 7401 (7357–7449) |
| Denser mutant, 10 y, first fit | 0.6–7.5% | 4.87 (4.76–5.54) | 0.72% (0.68–0.77%) | 16 | 1007 (817–1095) | 3035 (2999–3070) |
| Denser mutant, 10 y, refit | 3.4–13.8% | 5.13 (4.89–5.42) | 0.82% (0.76–0.88%) | 34 | 990 (888–1101) | 3107 (3077–3139) |
<!-- /mtdna-spread-de-novo-units -->

An affected zone holds more types than one of its regions (about 1,740
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
0.530. That gain is spread over hundreds of types. Chapter 17's drift cell
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

## What this shows and what it does not

What follows, for this model with these parameters: at 3 years, de novo
deletions bring 10% of regions to at least 60% mutants only at a deletion
probability of about 0.011 per replication, and then an affected region
holds well over a thousand mutant types, none of which exceeds 1% of its
copies. Survival of the densest does not change this: the neutral control
gives as many types per mutant copy. The affected regions of this
experiment look like the drift cells of chapters 16 and 17, only more so,
and unlike the essentially clonal affected cells found in experiments. In
3 years the mechanism has no time to act on a single new deletion: a new
type starts as one copy among 3,500 to 6,800 and, in about 75 copy
generations, grows only by drift. At 10 years the picture is the same with
about 820 types per affected region.

What the experiment does not show:

- **Not a test of the paper.** The paper makes no claim about short-lived
  animals at 3 years. Its estimate of the de novo rate concerns human fibers
  over decades, and it is based on the takeover probability of one founder
  copy. The Figure 2 parameters were chosen by the authors to model the
  wave in rhesus monkey fibers, not for a short-lived animal. Other
  parameters, in particular a smaller `Nss` (the paper's Figure 4 shows
  that a founder's takeover probability depends on `Nss`), could change the
  result.
- **A region is not a cell.** Chapters 16 and 17 classify cells of 1,000
  copies; here a region of 3,500 to 6,800 copies exchanges copies with its
  neighbors. A muscle fiber section classified as COX-negative in an
  experiment may cover one region or many. The zone and segment rows show
  how the count grows with the unit.
- **40 and 80 years are not run, and 10 years only for the denser
  mutant.** At 40 and 80 years chapter 17's drift cells hold 2.2 and 1.5
  mutant types, and survival of the densest has more time to act. Whether
  it then produces clonal regions is open; see
  [Running the experiment](#running-the-experiment) for the commands and
  costs.
- **All deletion types are alike.** Real deletions differ in size and
  perhaps in density; one with a larger density advantage than the others
  would not be modeled here.
- **The 10% target is met only loosely.** The prevalence scatters between
  cohorts more than the ±0.5-point band, because regions in a fiber are
  correlated. The type count barely depends on this (1,456 to 1,487 types
  for `m` between 0.01117 and 0.01130).
- **Fewer units than the papers.** 1,024 regions in 32 fibers per cohort,
  instead of 1,000 independent cells per run in the papers and 20,000 in
  chapter 17. The intervals above account for the clustering by
  resampling fibers.

The parameters, the boundary conditions, the chain length, the unit
definitions, the year of 360 days and the treatment of mutant templates
are choices of this guide, given with their reasons on
[page 19b](19b-mtdna-spread-de-novo.md#the-model).

## Running the experiment

| Size | Command | What it does | Time on one core |
| --- | --- | --- | --- |
| quick | `python examples/mtdna_spread/mutations.py run --size quick` | the whole protocol on a small model (4 regions of 20 copies, 120 days), pure Python | 3 s |
| sampler | `python examples/mtdna_spread/mutations.py sampler` | the sampler against CRN-craft and the moment law | 30 s with Numba |
| draft | `python examples/mtdna_spread/mutations.py run --size draft` | the Figure 2 parameters, 2 fibers per fit point and 8 validation fibers per case | not run; about 6 to 10 min with Numba (estimate from the full run) |
| full | `python examples/mtdna_spread/mutations.py run --size full --output examples/mtdna_spread/results` | the 3-year runs of this page, with the refits | 28 CPU-minutes with Numba (main 15, equal density 10, isolated 3) |
| full, 10 years | `python examples/mtdna_spread/mutations.py run --size full --years 10 --case main` | the 10-year runs (133 fibers with the refit) | 50 CPU-minutes with Numba |
| length | `python examples/mtdna_spread/mutations.py length --m 0.01116696641687773 --regions 64` | eight 64-region fibers | 2 min |

The draft and full sizes need Numba in practice: without it the same
sampler runs in pure Python, roughly 100 times slower. Every fiber is
checkpointed under `runs/de-novo/<size>/` (`--work`), so an interrupted run
continues where it stopped. `mutations.py figure --figures
docs/guide/img/mtdna_spread` redraws the figure from the saved tables.

The `run` command writes the tables of one lifespan; the saved tables hold
the 3-year rows followed by the 10-year rows. Outstanding: 10 years with
equal density and isolated regions (about 35 and 15 CPU-minutes), and 40
and 80 years (`--years 40`, `--years 80`). A 10-year fiber has about 185
million events and took 23 s; 40 and 80 years cost about 4 and 8 times as
much per fiber, so with a refit about 3.5 and 7 CPU-hours for the denser
mutant, and roughly 15 to 20 CPU-hours for all three cases at both ages.

## Exercises

1. Run `mutations.py run --size quick` and open `runs/de-novo/quick/de-novo-summary.csv`.
   Then make the regions larger: pass `nss=200` in `SIZES["quick"]["model"]`.
   How do the fitted `m` and the number of types per affected region change?
2. From `results/de-novo-regions.csv`, compute for each case the share of
   the largest type, `largest / (wild + mutant)`, in affected regions. Why
   is it about ten times larger in isolated regions than in the chain,
   although the deletion probability is almost the same?

See also: [chapter 17b](17b-mtdna-mechanisms-results.md), [chapter 16c](16c-fig1cd.md#result).
