# Kowald & Kirkwood (2013), Figure 1

Author: Axel Kowald. The script of [tutorial S2](../../docs/tutorials/simulation/S2-published-study.md),
which re-creates Figure 1 of Kowald A., Kirkwood T.B.L. (2013), *Aging Cell*
12(4):728–731, [doi:10.1111/acel.12098](https://doi.org/10.1111/acel.12098).

| File | What it is |
| --- | --- |
| `kowald_kirkwood_2013.py` | the study: `--size quick`, `draft` or `full`, `--part 1A`, `1B` or `1CD`, `--figures`, `--threads` |
| `make_tables.py` | the 0.8 tables of the chapter from the result files (`--write`, `--check`) |
| `results/quick/` | the quick size, compared exactly by `validation/check_kk13.py` |
| `results/draft/` | the draft run the chapter shows |

Measured on two cores (`--threads 2`): `quick` 16 s; `draft` about 25
minutes (Fig 1A 40 s with `compile_loop`, Fig 1B 5 minutes, Fig 1C/1D 17
minutes). `full` uses the author's cell counts; its Fig 1C/1D took about
2.6 hours on 18 threads in the author's run.

Every unit (one P_mut of Fig 1A, one lifespan of Fig 1B, one combination of
Fig 1C/1D) writes its own CSV and a JSON with the cell counts, seeds,
events, run time and CRN-craft version, and is skipped when its CSV exists:
an interrupted run continues where it stopped. Delete a unit's CSV to run
it again.
