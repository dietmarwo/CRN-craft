# C7 saved results

All files come from runs in review on one thread, with the crn-diffsol
binary of the same revision, fixtures under `DATA/C7` (SHA-256 checked) and
seeds 9971 + start (searches), 4983 (archived designs checked), 4984
(interior parity designs) and 7 (bootstrap).

| File | Origin |
| --- | --- |
| `terminal.csv` | `run.py C7 --preset study` (ka kinetics): terminal hypervolume per arm and start |
| `trajectory.csv` | same run: quartiles of the hypervolume after each evaluation, per arm |
| `fronts.csv` | same run: nondominated designs of each start, physical units |
| `parity.csv` | same run: native and Summit-RK45 objective errors against a tight SciPy solve, 24 designs |
| `stats.json` | same run and `--kinetics pinned`: checks, medians, Mann–Whitney test, bootstrap interval |
| `terminal-pinned.csv` | `run.py C7 --preset study --kinetics pinned`; identical to the implementer's 8-thread run |
| `archive-groups.csv` | derived from Summit's `snar_results.csv`: terminal hypervolume of the archived runs by strategy and transform (median, quartiles, runs) |
| `archive-front.csv` | derived from Summit's `pareto_front_snar.csv`: its nondominated designs (STY, E-factor) |

The figures `img/C7-*.png` are drawn from these tables by
`publications.reactor.draw`.
