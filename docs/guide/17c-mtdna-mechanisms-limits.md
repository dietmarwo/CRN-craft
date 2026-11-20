# Choices, checks and limits

Part of [chapter 17](17-mtdna-mechanisms.md).

## Reconstruction choices

The papers describe their models in words, equations and figures. Where
they leave something open, the reconstruction had to choose. All choices
are listed here; the sensitivity runs keep the fitted three-year rate and
change one choice (1,000 new cells each, in
`results/reconstruction/comparison.csv`).

| Choice | What the papers say | Reconstruction | Effect where measured |
| --- | --- | --- | --- |
| Maximal start probability `r` of the size model | Named (`replS`), no value | 1 per hour; `C` recomputed so the wild-type cell holds 1,000 copies | `r` = 0.1: wild type displaced after 94.954 instead of 94.952 years; 99 of 1,000 cells affected at 3 years, 28.86 lineages |
| Copy count in the size model's start probability | All copies | Free plus busy copies after that hour's degradation, one probability per hour | not varied |
| Degradation probability per hour | Size paper: `ln(2)/halfL`. PNAS: not stated | Size model `ln(2)/240`; ATP model `1 - exp(-ln(2)/240)` | 0.14% apart, not varied |
| Order within the ATP model's hour | Lists degradation, replication, mutation and ATP production | In that order: replication only by copies that survived the hour; ATP relaxes for one hour with the new counts held fixed | 0.5-h step: 100 of 1,000 affected, 1.26 lineages |
| A cell whose ATP reaches zero | "The system collapses" | The cell stops at that moment with its counts, and is classified like any other cell. CRN-craft's `lower_bound=0` does the same | Nearly all affected transcription cells are such cells (1,863 of 1,944 at 3 years). Without them only 81 of 18,137 cells are affected, with 1.06 lineages |
| Which templates can produce a deletion in the ATP model | Not stated | Daughters of wild-type and mutant templates | Only wild-type templates: 101 of 1,000 affected and 1.01 lineages instead of about 1.25 |
| Year length | Not stated. The size paper's 19.33 years for a 10-h copying time fits a 365-day year (19.60 years with 360 days); chapter 16 uses 360 days | 365 days | 360 days, 3 years: 100, 96 and 97 affected; 33.44, 29.41 and 1.26 lineages. 80 years: 115, 104 and 105; 1.48, 1.09 and 1.01 (drift, size, transcription) |
| 60% threshold | Text of both papers: "at least 60%". Size paper figure captions and chapter 16: more than 60% | At least 60% (`5*M >= 3*(W+M)`); the strict count is saved as well | The strict rule removes 36 drift and 5 size cells at 3 years, no transcription cell |
| What a "type" is | PNAS: different mutants. Size paper text: different mtDNAs. 2013 table: matches with wild type included | Mutant lineages; the wild type is added for the 2013 and PNAS drift values | Adds one type to almost every affected 3-year cell |
| Drift time course | 2013 program steps hourly | Continuous time, the exact pair rule of chapter 16c | not varied |
| ATP equations | Only as images in the PNAS article; parameters in Table 1 | The equations of [the models page](17a-mtdna-mechanisms-models.md#transcription-linked-replication), from the text and Table 1 | Deterministic collapse after 247 days (paper: about 250) |
| How the rate is adjusted | Adjusted until 10 ± 0.5% of 1,000 cells are affected | Logistic fit, then an independent 20,000-cell cohort ([protocol](17b-mtdna-mechanisms-results.md#the-protocol)) | — |

Two of these choices carry the transcription result. With collapsed cells
counted as affected and mutant templates able to produce new deletions,
the reconstruction gives 1.25 lineages at three years, close to the
published 1.26. Without mutations from mutant templates it gives 1.01, and
the affected cells are almost all cells that ran out of ATP. The papers do
not say which convention their programs used.

## How the reconstruction was checked

`validation/check_kk14.py` runs in CI and in a few seconds:

- **Endpoints**: the inclusive and strict threshold in integers, empty
  cells, and both type counts on hand-made endpoints.
- **Deterministic targets**: the ATP collapse (247.41 days ± 0.1), the
  wild-type displacement (94.95 years ± 0.05), the busy fractions, mesh
  convergence, and the 2-day, 10-h and 70% cases of the size paper's
  Figures 3 to 5.
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
  `reference.cpp`, all 12 fits and 42 reference rows are present, the
  three-year decision rule replays, and the table of
  [the results page](17b-mtdna-mechanisms-results.md#results) equals the
  table generated from the saved CSV files.

The check does not rerun the tissue study and compares no simulated cell
with a saved value; its sampling checks use bounds. The saved results were
produced on an AMD x86-64 host with GCC 14.2. A rerun of the full tier on
an Intel x86-64 host with GCC 13.3 and Python 3.13 reproduced every saved
table byte for byte. A byte-for-byte replay needs the same C++ standard
library and C math library: `std::uniform_int_distribution` and
`std::binomial_distribution` are implementation-defined, so libc++ (macOS)
draws different, equally valid cells.

## Limits

- The original programs (Mathematica for the size paper, Java for the
  PNAS paper) were not available. Where the text is silent, the table
  above is a choice, not a fact about the published runs.
- The PNAS equations are available only as images in the online article,
  and the values of PNAS Figure 3 and the size paper's Figures 6 and 7 were
  read from the figures with an interval of ±2 pixels. The size paper has
  an erratum (doi:10.1016/j.jtbi.2014.02.008); its text was not
  available to the reconstruction.
- The size paper is not self-consistent in one place: its Results give
  almost 160 years for the 70%-size mutant, its Discussion almost 150 years
  for a mutant 30% shorter. The reconstruction gives 158 years.
- The published values have sampling errors that the papers do not report.
  The reconstruction's intervals describe its own cells only.
- The full tier ran with one seed schedule; the draft tier with another.
  Both are kept.
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
   The saved refit gives **1.11 mutant lineages** and 2,014 affected cells
   out of 20,000.
3. Compare results with a 360-day year at 3 and 80 years, at the same
   mutation probability. In the ATP model, check what happens when
   mutations are permitted only from wild-type templates or when the
   time step is halved. `compare --size quick --ages 3 --sensitivities`
   runs all of these with 128 cells each in seconds, `--size draft` with
   1,000.

See also: [Reproducing a published simulation](08-reproducibility.md#reproducing-a-published-simulation),
[Limits and lessons of chapter 16](16d-limits-lessons.md).
