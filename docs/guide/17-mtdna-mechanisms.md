# 17. Case study: comparing mtDNA mechanisms

## What this chapter reproduces

Old post-mitotic cells such as muscle fibers and neurons accumulate
mtDNA deletion mutants. In an affected cell, the mutant copies usually
descend from a single deletion. Two papers from 2014 used this observation
to test three candidate mechanisms. Each mechanism was given the mutation
rate that makes 10% of simulated cells COX-negative at the end of a
virtual lifespan of 3, 10, 40 or 80 years of 360 days. Each mechanism was
then asked how many different mutants an affected cell contains:

> Kowald A., Kirkwood T.B.L. (2014). Transcription could be the key to the
> selection advantage of mitochondrial deletion mutants in aging.
> *PNAS* 111(8):2972–2977.
> [doi:10.1073/pnas.1314970111](https://doi.org/10.1073/pnas.1314970111)
> ([PMC3939916](https://pmc.ncbi.nlm.nih.gov/articles/PMC3939916/), Figures 2 and 3)
>
> Kowald A., Dawson M., Kirkwood T.B.L. (2014). Mitochondrial mutations
> and ageing: Can mitochondrial deletion mutants accumulate via a size
> based replication advantage? *Journal of Theoretical Biology*
> 340:111–118. [doi:10.1016/j.jtbi.2013.09.009](https://doi.org/10.1016/j.jtbi.2013.09.009)
> (Figures 2 to 7)

The drift baseline comes from the study re-created in
[chapter 16](16-published-study.md) (Kowald & Kirkwood 2013, *Aging Cell*
12:728–731). This chapter is a reconstruction from the published text and
figures. The original programs were not available to it. The papers'
authors did not write or review it.

The result in three lines:

- At three years, an affected cell holds 34.7 mutant lineages under neutral
  drift, 29.3 under the size advantage and 1.26 under transcription-linked
  replication (20,000 new cells per mechanism after the fit). The
  PNAS paper reports 1.26 for the transcription mechanism and about 30 for
  the others.
- The difference shrinks with lifespan. At 80 years, every mechanism gives
  1.0 to 1.5 lineages per affected cell.
- The transcription result depends on two points that the papers do not
  state and that an author of the papers confirmed (personal
  communication, November 2026): a cell whose ATP runs out stays as it is
  for the rest of the simulation and is counted, and mutant templates
  mutate too. Both, and the 360-day year, are explained on
  [the last page](17c-mtdna-mechanisms-limits.md).

## Words used in this chapter

- **COX-negative** (or *affected*): a cell whose respiratory chain fails
  because too much of its mtDNA is mutant. Both papers use 60% mutant
  copies as the threshold. A cell with `M` mutant and `W` wild-type copies
  is affected when `M >= 0.6 * (W + M)`, computed in integers as
  `5*M >= 3*(W+M)`. A cell with no copies at all is not classified.
- **Lineage**: all descendants of one simulated deletion event. Every
  deletion creates a new lineage, also when the template was already a
  mutant. A faithful copy belongs to its template's lineage.
- **Types**: the papers count "different mtDNAs" or "different mutants"
  per affected cell. This chapter reports the number of surviving mutant
  lineages. The drift table of the 2013 paper also counts the wild type as
  one type when any wild-type copy survives
  ([chapter 16c](16c-fig1cd.md#result) settled this). Where a published
  number uses that convention, the comparison uses it too, and the result
  files give both counts.

## The three mechanisms in plain words

| Mechanism | What makes a deletion mutant spread | Page |
| --- | --- | --- |
| Neutral drift | Nothing. Copies degrade at random and survivors replace them, so a mutant spreads only by chance. | [17a](17a-mtdna-mechanisms-models.md#neutral-drift) |
| Size advantage | A half-size mutant spends 1 h instead of 2 h being copied. Only free copies can start a new replication, so the mutant is available more often. | [17a](17a-mtdna-mechanisms-models.md#size-advantage) |
| Transcription-linked replication | Mutants start replication 1.5 times as often as wild type. ATP that only wild-type copies produce slows all replication. | [17a](17a-mtdna-mechanisms-models.md#transcription-linked-replication) |

The pages of the chapter:

| Page | Question |
| --- | --- |
| [The three models](17a-mtdna-mechanisms-models.md) | What exactly is simulated, and how does CRN-craft express it? |
| [Fitting and the comparison](17b-mtdna-mechanisms-results.md) | How is the 10% target fitted, and what do new cells show? How close are the published values? |
| [Choices, checks and limits](17c-mtdna-mechanisms-limits.md) | What did the reconstruction have to decide, how was it checked, and what does it not show? |

## CRN-craft and the reference sampler

The single-cell models run on CRN-craft: the drift pair rule on the
`tree` method, a staged version of the size model on `static_direct`, and
the ATP model on `hybrid` and `thinning` (the
[ATP chapter](12c-mtdna-atp.md) introduces both). The tissue study, about
420,000 simulated cells, runs on a small C++17 sampler local to the
example (`examples/kk14/reference.cpp`). It changes nothing in CRN-craft.

The reference sampler exists for two reasons:

1. Two of the published programs are hourly loops, not reaction networks.
   In the size model a copy is busy for exactly 2 h; a reaction network can
   only approximate a fixed delay with a chain of exponential stages. The
   transcription program applies degradation, replication, mutation and
   ATP production once per hour. The sampler follows those loops, so its
   differences from the papers are differences in the reconstructed
   model, not in the time stepping.
2. Speed. On a two-core x86-64 host, a three-year drift cell took about
   16 ms on CRN-craft's `tree` method (76,000 events) and 0.75 ms in the
   reference sampler. A three-year ATP cell took about 15 ms on
   `thinning` and 1 ms in the hourly sampler. The calibration needs mostly
   10- to 80-year cells. The full study took 136 minutes of sampling with
   one thread on a two-core Intel x86-64 host with the reference sampler,
   and would take roughly 20 times as long natively.

For neutral drift both run the same continuous-time process, and the
validation script checks both against its exact count law
([checks](17c-mtdna-mechanisms-limits.md#how-the-reconstruction-was-checked)).

## Running the example

You need Python 3.10 or newer, CRN-craft 0.8.0, SciPy, Matplotlib and a
C++17 compiler (set `CXX` to choose one). The compiler builds the reference
sampler in a temporary directory the first time a command needs it. Run
the commands from the repository root:

```sh
python -m pip install 'crn-craft==0.8.0' scipy matplotlib
python examples/kk14/study.py single-cell --size quick
python examples/kk14/study.py compare --size quick --ages 3 --threads 2
python examples/kk14/study.py figures --results examples/kk14/results/reconstruction
```

`single-cell` solves the two deterministic single-cell models in a few
seconds. `compare --size quick` simulates 128 new cells per mechanism at
the saved three-year mutation rates; its intervals are wide. `figures`
redraws the saved comparison without simulating. The
[example README](../../examples/kk14/README.md) lists all commands and
output files; [the results page](17b-mtdna-mechanisms-results.md#rerunning-the-study)
gives the commands of the full study.

See also: [Case study: mtDNA mutants](12-mtdna.md),
[Variable population: ATP](12c-mtdna-atp.md),
[Re-creating a published study](16-published-study.md).

Continue with [chapter 18](18-mtdna-inference.md) to infer mutation and
selection parameters and estimate takeover times from cell snapshots.
