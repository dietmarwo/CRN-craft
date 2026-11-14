# mtDNA examples

Runnable examples for the [case study of the user guide](../../docs/guide/12-mtdna.md), the
[fitting chapter of the guide](../../docs/guide/13-fitting.md) and the
[performance analysis](../../docs/evidence/mtdna-performance-analysis.md). Run them
from the repository root with the `crn-craft` Python package installed. The
fitting scripts also need `fcmaes` and `scipy`. They build the model in
Python and use the top-level module `crn_craft`, as the guide does; with
`CRN_CRAFT_SMOKE=1` they run with few cells and three evaluations.

The numbered scripts are the complete programs of the user manual of
release 0.7.0, which the guide replaced. They use the module
`crn_craft.rules` and the model files in `models/rules/`, so they show the
file form of the models, which the guide teaches on
[page 14a](../../docs/guide/14a-mtdna-model-files.md). The guide's own
blocks, which build the models in Python, are checked by
`validation/check_guide.py`.

| Script | Topic, and where the guide covers it |
| --- | --- |
| `01_first_run.py` | a first ensemble; [12a](../../docs/guide/12a-mtdna-takeover.md) |
| `02_reproducibility.py` | seeds and threads; [3](../../docs/guide/03-running.md#seeds-and-threads) |
| `03_trackers.py` | trackers; [12b](../../docs/guide/12b-mtdna-events.md) |
| `04_event_log.py` | event logs; [12b](../../docs/guide/12b-mtdna-events.md#mutation-events-in-a-log) |
| `04_callback.py` | live callbacks; [12b](../../docs/guide/12b-mtdna-events.md#stop-at-the-first-advantaged-mutant) |
| `05_methods.py` | choosing a method; [12d](../../docs/guide/12d-mtdna-methods.md) |
| `05_structure.py` | the `mtdna_direct` fast path; [14a](../../docs/guide/14a-mtdna-model-files.md#the-fast-path) |
| `06_clones.py` | the clone model; [14a](../../docs/guide/14a-mtdna-model-files.md#every-lineage-with-the-clone-model) |
| `13_continuous_atp.py` | continuous variables and `hybrid`; [12c](../../docs/guide/12c-mtdna-atp.md) |
| `14_independent_ode.py` | independent ODEs; [9](../../docs/guide/09-continuous-variables.md#independent-variables) |
| `07_atp.py` | the stored ATP block and `thinning`; [14a](../../docs/guide/14a-mtdna-model-files.md#the-stored-atp-block-and-thinning) |
| `08_sweep.py` | parameter sweeps; [12e](../../docs/guide/12e-mtdna-large-studies.md#parameter-sweeps) |
| `12_builder.py` | building models in Python; [12](../../docs/guide/12-mtdna.md#the-model) |
| `09_variants.py` | model variants; [12a](../../docs/guide/12a-mtdna-takeover.md#variants) |
| `10_budgets.py` | budgets; [3](../../docs/guide/03-running.md#budgets) |
| `11_cli_request.py [DIR [CRN_CRAFT]]` | the command line; [12e](../../docs/guide/12e-mtdna-large-studies.md#from-the-command-line) |
| `fit_takeover.py [EVALUATIONS]` | guide, chapter 13 |
| `fit_heteroplasmy.py [EVALUATIONS]` | guide, chapter 13 |
| `noise_and_parallel.py` | guide, chapter 13 |
| `performance.py SECTION [CRN_CRAFT]` | performance analysis |
