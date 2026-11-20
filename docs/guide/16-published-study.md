# 16. Case study: re-creating a published study (Kowald & Kirkwood 2013)

Author: Axel Kowald

This tutorial re-creates all four panels of Figure 1 of a published
stochastic simulation study with the Python interface of CRN-craft (0.8):

> Kowald A., Kirkwood T.B.L. (2013). Mitochondrial mutations and aging:
> random drift is insufficient to explain the accumulation of mitochondrial
> deletion mutants in short-lived animals. *Aging Cell* 12(4):728–731.
> [doi:10.1111/acel.12098](https://doi.org/10.1111/acel.12098) (open access)

The tutorial was written for CRN-craft 0.7.1. The maintainers adapted the
text and the code to CRN-craft 0.8 (this version is a draft for the
author's review);
[the last page](16d-limits-lessons.md#what-changed-since-071) lists what
that changed.

The tutorial is a worked example of how to take a model from a paper's
Methods section and turn it into CRN-craft simulations, and how to check
the result against the published numbers. It covers several features:

| Panel | Question | CRN-craft features | Page |
| --- | --- | --- | --- |
| 1A | Fraction of COX-negative cells over 120 years | `crn.Model` builder, `static_direct`, `compile_loop`, batched ensembles | [Fig 1A](16a-fig1a.md) |
| 1B | Mutation rate needed for 5/10/15 % COX-negative cells | parameter sweep with independent seeds, logistic regression | [Fig 1B](16b-fig1b.md) |
| 1C/1D | Number of mtDNA types and size of the largest clone | species with attributes (`add_species(attributes=)`), `crn.copy`, `crn.new` with `new_id()`, `reduce="distinct"` and `reduce="max_count"` observables, tree method | [Fig 1C/1D](16c-fig1cd.md) |

The [last page](16d-limits-lessons.md) explains the limits (event,
sampling, instance, memory and output budgets) that a long run like this
one meets, which of them still need attention in 0.8, and the lessons of
the study.

## The biological model

The paper asks whether pure random drift of mtDNA deletion mutants, with no
selective advantage, can explain why COX-negative cells (cells without
functional respiratory chains) accumulate during aging. Its model of a
post-mitotic cell is:

- A cell holds N = 1000 mtDNA molecules.
- Every molecule is degraded with a half-life of 10 days. Degradation is
  first order with δ = ln 2 / 240 h⁻¹.
- After each degradation a randomly chosen survivor replicates, so N stays
  constant (a Moran process).
- Every replication, of a wild-type or of a mutant template, produces a new
  deletion mutant with probability P_mut. The copy of a mutant is always a
  mutant.
- A cell is COX-negative when more than 60 % of its mtDNA is mutant.

The paper advances time in hourly steps. We simulate the same process
exactly in continuous time (Gillespie), with time measured in hours and a
360-day year. Every page of the case study starts with these constants:

```python
import math
import os

import numpy as np
import crn_craft as crn

nMolecules   = 1000
halfLifeDays = 10.0
degRate      = math.log(2)/(24.0*halfLifeDays)                   # 1/h
coxThreshold = 0.6
hoursPerYear = 24.0 * 360.0
nThreads     = max(1, (os.cpu_count() or 1) - 2)
```

A result does not depend on the number of threads
([chapter 3](03-running.md#seeds-and-threads)), so `nThreads` only sets the
speed.

## The script

The complete script is
[`examples/kk13/kowald_kirkwood_2013.py`](https://github.com/dietmarwo/CRN-craft/blob/main/examples/kk13/kowald_kirkwood_2013.py).
Each figure has its own function: `KowaldKirkwood13Fig1A()`,
`KowaldKirkwood13Fig1B()` and `KowaldKirkwood13Fig1CD()`. The model
functions of the pages, `buildModel(pmut)` and `buildCloneModel(pmut)`, are
those of the script. Run it from the repository root:

```sh
python examples/kk13/kowald_kirkwood_2013.py --size quick                # all parts, small cell counts
python examples/kk13/kowald_kirkwood_2013.py --size draft                # the 0.8 tables of these pages
python examples/kk13/kowald_kirkwood_2013.py --size full                 # the author's cell counts
python examples/kk13/kowald_kirkwood_2013.py --size full --part 1A       # one part: 1A, 1B or 1CD
python examples/kk13/kowald_kirkwood_2013.py --size draft --figures      # plots from the stored results
```

| Size | Fig 1A cells | Fig 1B fine sweep | Fig 1C/1D cells per combination |
| --- | --- | --- | --- |
| `quick` | 200 per P_mut, 10 years | 200 per point, lifespans 3 and 10 y | 100, 3 years only |
| `draft` | 3,000 / 3,000 / 15,000 | 1,000 per point | 200 |
| `full` | 3,000 / 3,000 / 15,000 | 3,000 per point | 3,000 |

The result files of a run go to `examples/kk13/results/<size>/`, one file
per P_mut, lifespan or combination, together with a JSON file that records
the cell counts, seeds, events, run time and CRN-craft version. A unit that
already has its file is skipped, so a long run can be stopped and
continued. `--figures` writes the figures to `docs/guide/img/kk13/`, and
`examples/kk13/make_tables.py` fills the 0.8 tables of these pages from
the result files.

Each result page shows two runs: the author's run with CRN-craft 0.7.1 and
the full cell counts, and a draft run with CRN-craft 0.8 and the `draft`
sizes (about 25 minutes on two cores). The draft run has larger
statistical errors for Fig 1C/1D; its purpose is to show that 0.8 gives the
same results.

The code blocks of the pages are reduced versions of the script's code: a
few hundred cells and short horizons, so that each runs in seconds. Each
page says which script call gives the full result.

See also: [Fig 1A](16a-fig1a.md),
[Case study: mtDNA mutants](12-mtdna.md) (a different mtDNA model, with
mutant classes and a replication advantage),
[Reproducing a published simulation](08-reproducibility.md#reproducing-a-published-simulation).

[Chapter 17](17-mtdna-mechanisms.md) reuses the drift baseline to compare
fixed replication delays and transcription-linked replication.
[Chapter 18](18-mtdna-inference.md) fits a Moran model of the same fixed
copy number to cell snapshots.
