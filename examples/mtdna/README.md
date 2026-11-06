# mtDNA examples

Runnable examples for the [user manual](../../docs/mtdna-user-manual.md), the
[optimization tutorial](../../docs/mtdna-optimization-tutorial.md) and the
[performance analysis](../../docs/mtdna-performance-analysis.md). Run them
from the repository root with the `crn-craft` Python package installed. The
fitting scripts also need `fcmaes` and `scipy`.

Every Python block of the manual runs on its own, and
`python validation/check_manual.py --binary target/release/crn-craft` runs
them all in fresh interpreters and compares their output with the manual.

| Script | Manual section |
| --- | --- |
| `01_first_run.py` | 3, a first ensemble |
| `02_reproducibility.py` | 4, seeds and threads |
| `03_trackers.py` | 6, trackers |
| `04_event_log.py` | 7, event logs |
| `04_callback.py` | 7, live callbacks |
| `05_methods.py` | 8, choosing a method |
| `05_structure.py` | 8, the `mtdna_direct` fast path |
| `06_clones.py` | 9, the clone model |
| `07_atp.py` | 10, the ATP model |
| `08_sweep.py` | 11, parameter sweeps |
| `09_variants.py` | 12, model variants |
| `10_budgets.py` | 13, budgets |
| `11_cli_request.py [DIR [CRN_CRAFT]]` | 14, the command line |
| `12_builder.py` | the `rules.Model` builder |
| `13_continuous_atp.py` | continuous variables and `hybrid` (replaces `07_atp.py` in the ATP chapter) |
| `fit_takeover.py [EVALUATIONS]` | tutorial, section 2 |
| `fit_heteroplasmy.py [EVALUATIONS]` | tutorial, section 3 |
| `noise_and_parallel.py` | tutorial, sections 4 and 5 |
| `performance.py SECTION [CRN_CRAFT]` | performance analysis |
