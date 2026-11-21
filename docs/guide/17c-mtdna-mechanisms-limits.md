# Choices, checks and limits

Part of [chapter 17](17-mtdna-mechanisms.md).

## Reconstruction choices

The papers describe their models in words, equations and figures. Where
they leave something open, the reconstruction had to choose. Three points
that the papers leave open were confirmed by an author of the papers
(personal communication, November 2026): mutant templates can mutate, a
year has 360 days, and a cell whose ATP reaches zero remains in that state
for the rest of the simulation, with fixed ATP, wild-type and mutant
counts. All choices are listed here; the sensitivity runs keep the fitted
rate and change one choice (1,000 new cells each, in
`results/reconstruction/comparison.csv`); for a confirmed point the
sensitivity shows what the other convention would give.

| Choice | What the papers say | Reconstruction | Effect where measured |
| --- | --- | --- | --- |
| Maximal start probability `r` of the size model | Named (`replS`), no value | 1 per hour; `C` recomputed so the wild-type cell holds 1,000 copies | `r` = 0.1: wild type displaced after 831,793 instead of 831,783 h; 97 of 1,000 cells affected at 3 years, 30.41 lineages (29.30 with 20,000 cells at `r` = 1) |
| Copy count in the size model's start probability | All copies | Free plus busy copies after that hour's degradation, one probability per hour | not varied |
| Degradation probability per hour | Size paper: `ln(2)/halfL`. PNAS: not stated | Size model `ln(2)/240`; ATP model `1 - exp(-ln(2)/240)` | 0.14% apart, not varied |
| Order within the ATP model's hour | Lists degradation, replication, mutation and ATP production | In that order: replication only by copies that survived the hour; ATP relaxes for one hour with the new counts held fixed | 0.5-h step: 92 of 1,000 affected, 1.27 lineages |
| A cell whose ATP reaches zero | PNAS: "The system collapses". Confirmed by an author: the cell remains in that state for the rest of the simulation, with fixed ATP, wild-type and mutant counts | The sampler stops the cell at that moment and keeps its counts and lineages; at the end of the lifespan it is classified and its lineages counted like any other cell. CRN-craft's `lower_bound=0` (`stopped`, counts kept, ATP 0 at later reports) does the same | Nearly all affected transcription cells are such cells (1,904 of 1,976 at 3 years). Without them only 72 of 18,096 cells are affected, with 1.03 lineages |
| Which templates can produce a deletion | Not stated in the papers. Confirmed by an author: mutants can mutate | Daughters of wild-type and mutant templates, in all three mechanisms (drift, size, transcription), in the reference sampler and in the native models (`new_lineage`, `mutate_M`) | Only wild-type templates (ATP model): 92 of 1,000 affected and 1.04 lineages instead of 1.26 |
| Year length | Not stated in the papers. Confirmed by an author: 360 days; chapter 16 uses 360 days too. The size paper's printed 19.33 years fits 365 days (open question, see [limits](#limits)) | 360 days: lifespans of 1,080, 3,600, 14,400 and 28,800 days | 365 days, 3 years: 108, 109 and 111 affected; 34.44, 29.59 and 1.32 lineages. 80 years: 111, 97 and 80; 1.51, 1.04 and 1.01 (drift, size, transcription) |
| 60% threshold | Text of both papers: "at least 60%". Size paper figure captions and chapter 16: more than 60% | At least 60% (`5*M >= 3*(W+M)`); the strict count is saved as well | The strict rule removes 38 drift and 2 size cells at 3 years, no transcription cell |
| What a "type" is | PNAS: different mutants. Size paper text: different mtDNAs. 2013 table: matches with wild type included | Mutant lineages; the wild type is added for the 2013 and PNAS drift values | Adds one type to almost every affected 3-year cell |
| Drift time course | 2013 program steps hourly | Continuous time, the exact pair rule of chapter 16c | not varied |
| How the rate is adjusted | Adjusted until 10 ± 0.5% of 1,000 cells are affected | Logistic fit, then an independent 20,000-cell cohort ([protocol](17b-mtdna-mechanisms-results.md#the-protocol)) | — |

The ATP equations themselves are not a choice. They are Equations 1 to 3
of the PNAS paper, and the formula for `c` in its Table 1 evaluates to the
same 1.303.

Two of these points carry the transcription result, and an author of the
papers confirmed both (personal communication, November 2026): a cell
whose ATP runs out keeps its counts and is counted at the end of the
lifespan, and mutant templates produce new deletions. With both, the
reconstruction gives 1.26 lineages at three years, the value printed in
the PNAS paper. With mutations from wild-type templates only it would give
1.04 (1,000 cells), and the affected cells would still be almost all cells
that ran out of ATP.

## How the reconstruction was checked

`validation/check_kk14.py` runs in CI and in a few seconds:

- **Endpoints**: the inclusive and strict threshold in integers, empty
  cells, and both type counts on hand-made endpoints.
- **Deterministic targets**: the ATP collapse (247.41 days ± 0.1), the
  wild-type displacement in hours (94.95 ± 0.05 printed years, converted
  with 365 days), the busy fractions, mesh convergence, and the 2-day,
  10-h and 70% cases of the size paper's Figures 3 to 5.
- **Exact count law of drift**: 16,000 eight-copy cells from the count
  sampler, the lineage sampler and CRN-craft's `tree` pair rule, each
  against the exact distribution from the matrix exponential. The largest
  empirical CDF gap must stay below the Dvoretzky–Kiefer–Wolfowitz bound
  for 16,000 cells at error probability 0.001 (0.0154). Observed: 0.0038,
  0.0039 and 0.0069.
- **Hourly size queue**: 16,000 eight-copy cells from the sampler against
  16,000 from an independent per-copy Python implementation of the three
  hourly tasks, for wild-type count, mutant count and lineage count. The
  two-sample bound is twice the one-sample bound at 0.0005 (0.0322).
  Observed: 0.0081, 0.0077 and 0.0071.
- **Replay and marking**: the same cells with one and three threads, and
  with and without lineage marking, are identical; the native hybrid run
  with the lineage callback matches the run without it.
- **Saved evidence**: the sampler hash in the manifest matches
  `reference.cpp`, all 12 fits are refitted on the 360-day year and the
  six 365-day sensitivity rows are present, the 42 reference rows are
  present, the three-year decision rule replays, the table of
  [the results page](17b-mtdna-mechanisms-results.md#results) equals the
  table generated from the saved CSV files, and the deterministic size
  rows on [the models page](17a-mtdna-mechanisms-models.md#one-cell-deterministic)
  equal `single-cell.json`.

The check does not rerun the tissue study and compares no simulated cell
with a saved value; its sampling checks use bounds. The saved results were
produced on one Intel Xeon x86-64 host (2.30 GHz, two cores) with GCC
13.3.0 and Python 3.13. They have not been replayed on a second host. A
byte-for-byte replay needs the same C++ standard
library and C math library: `std::uniform_int_distribution` and
`std::binomial_distribution` are implementation-defined, so libc++ (macOS)
draws different, equally valid cells.

## Limits

- The original programs (Mathematica for the size paper, Java for the
  PNAS paper) were not available. Where the text is silent and no author
  confirmed the point, the table above is a choice, not a fact about the
  published runs.
- Open question for the authors: the size paper prints 19.33 years for a
  10-h wild-type copying time. The reconstruction gives 169,333 h, which
  is 19.33 years of 365 days but 19.60 years of 360 days, and its other
  three deterministic times also fit 365 days at least as well
  ([models](17a-mtdna-mechanisms-models.md#one-cell-deterministic)).
  Perhaps the size paper's deterministic program used 365 days while the
  PNAS tissue program used 360. The tissue comparison with the size
  paper's Figures 6 and 7 assumes 360 days; the 365-day sensitivity above
  shows how much that matters.
- The values of PNAS Figure 3 and the size paper's Figures 6 and 7 were
  read from the figures with an interval of ±2 pixels. In panel B of the
  PNAS figure the markers at 40 and 80 years overlap; they are read from
  their visible edges, and ±2 pixels is a quarter of a type there. The
  size paper has an erratum (doi:10.1016/j.jtbi.2014.02.008); its text was
  not available to the reconstruction.
- The size paper is not self-consistent in one place: its Results give
  almost 160 years for the 70%-size mutant, its Discussion almost 150 years
  for a mutant 30% shorter. The reconstruction gives 158 years.
- The published values have sampling errors that the papers do not report.
  The reconstruction's intervals describe its own cells only.
- The full tier ran with one seed schedule; the draft tier with another,
  and on a 365-day year. Both are kept.
- Matching the published numbers does not validate a mechanism
  biologically. It shows that the models as reconstructed behave as the
  papers report.

## Exercises

1. Switch from the inclusive to the strict 60% boundary. Compare the
   affected counts and lineage means in the saved
   `comparison.csv` (columns `strict_negative` and `strict_mutant_types`).
2. Fit the size model again at a **2-day** half-life. Predict the direction
   of the change in lineage diversity before viewing the fit. Figure 7 of
   the size paper reports **1.16 types per affected cell** at this turnover.
   The saved refit gives **1.12 mutant lineages** and 1,850 affected cells
   out of 20,000 (9.25%, an interval of 8.86–9.66%).
3. Compare results with a 365-day year at 3 and 80 years, at the same
   mutation probability. In the ATP model, check what happens when
   mutations are permitted only from wild-type templates or when the
   time step is halved. `compare --size quick --ages 3 --sensitivities`
   runs all of these with 128 cells each in seconds, `--size draft` with
   1,000.

See also: [Reproducing a published simulation](08-reproducibility.md#reproducing-a-published-simulation),
[Limits and lessons of chapter 16](16d-limits-lessons.md).
