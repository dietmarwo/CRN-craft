# CRN-craft user manual: mtDNA mutant dynamics

The manual is now part of the [user guide](guide/01-install.md). What it
taught in general (running, results, methods, the command line) is in the
guide's chapters; the mtDNA models are the guide's
[case study](guide/12-mtdna.md), and their model files, the fast path and
the clone model are on the page
[The mtDNA model files](guide/14a-mtdna-model-files.md). The table says
where each section went.

Scripts written for release 0.6.0 against `crn_craft.rules` keep running.
The manual of that release is kept unchanged as
[validation/compat/manual-0.6.0.md](https://github.com/dietmarwo/CRN-craft/blob/main/validation/compat/manual-0.6.0.md),
and CI runs every block of it; [Migrating](guide/migrating.md) lists what
changed. The manual as it was in release 0.7.0 is in the repository at the
tag
[v0.7.0](https://github.com/dietmarwo/CRN-craft/blob/v0.7.0/docs/mtdna-user-manual.md).

| Manual section | Now |
| --- | --- |
| Introduction | [12, What this case study covers](guide/12-mtdna.md#what-this-case-study-covers) |
| 1. Installation | [1, Install](guide/01-install.md); [14a, The five model files](guide/14a-mtdna-model-files.md#the-five-model-files) |
| 2. The models | [12, The model](guide/12-mtdna.md#the-model); [14a, The five model files](guide/14a-mtdna-model-files.md#the-five-model-files) |
| The three-class model | [12, Parameters and rates](guide/12-mtdna.md#parameters-and-rates) |
| 3. A first ensemble | [12a, A first ensemble](guide/12a-mtdna-takeover.md#a-first-ensemble); [4, What a result holds](guide/04-results.md#what-a-result-holds) |
| Quantiles and the median | [4, Quantiles and the median](guide/04-results.md#quantiles-and-the-median); [12a, The mean hides the typical cell](guide/12a-mtdna-takeover.md#the-mean-hides-the-typical-cell) |
| Saving results | [4, Save a run in one file](guide/04-results.md#save-a-run-in-one-file) |
| 4. Seeds, threads and reproducibility | [3, Seeds and threads](guide/03-running.md#seeds-and-threads); [4, Single trajectories](guide/04-results.md#single-trajectories) |
| 5. Takeover, thresholds and stop conditions | [6, Stop conditions](guide/06-stops-trackers-logs.md#stop-conditions); [12a, Takeover fraction and time to threshold](guide/12a-mtdna-takeover.md#takeover-fraction-and-time-to-threshold) |
| 6. Event-level bookkeeping with trackers | [6, Trackers](guide/06-stops-trackers-logs.md#trackers); [12b, Class appearances and extinctions](guide/12b-mtdna-events.md#class-appearances-and-extinctions) |
| 7. Exact event records: logs and callbacks | [6, Event logs](guide/06-stops-trackers-logs.md#event-logs); [6, Callbacks](guide/06-stops-trackers-logs.md#callbacks); [12b, Mutation events in a log](guide/12b-mtdna-events.md#mutation-events-in-a-log); [12b, Stop at the first advantaged mutant](guide/12b-mtdna-events.md#stop-at-the-first-advantaged-mutant) |
| 8. Choosing a method | [7, Which method for which model](guide/07-methods.md#which-method-for-which-model); [12d, The same cells under every method](guide/12d-mtdna-methods.md#the-same-cells-under-every-method) |
| Selectors for Static Direct | [7, Static Direct and its selectors](guide/07-methods.md#static-direct-and-its-selectors) |
| The opt-in fast path `mtdna_direct` | [14a, The fast path](guide/14a-mtdna-model-files.md#the-fast-path); [reference: the fast path](reference/methods.md#the-fast-path-mtdna_direct); [14, A model document for three strains](guide/14-model-files.md#a-model-document-for-three-strains) |
| Optional: compiled event loop | [5, Compiled models](guide/05-speed.md#compiled-models); [12d, A compiled model](guide/12d-mtdna-methods.md#a-compiled-model); [reference: compiled event loops](reference/methods.md#compiled-event-loops) |
| 9. Individual mutant lineages: the clone model | [14a, Every lineage with the clone model](guide/14a-mtdna-model-files.md#every-lineage-with-the-clone-model) |
| 10. Variable population: continuous variables | [12c, ATP as a continuous variable](guide/12c-mtdna-atp.md#atp-as-a-continuous-variable); [9, A variable between events](guide/09-continuous-variables.md#a-variable-between-events) |
| Independent ODEs | [9, Independent variables](guide/09-continuous-variables.md#independent-variables) |
| Stored models: the ATP block and `thinning` | [14a, The stored ATP block and thinning](guide/14a-mtdna-model-files.md#the-stored-atp-block-and-thinning) |
| 11. Parameter sweeps | [5, Parameter sweeps](guide/05-speed.md#parameter-sweeps); [12e, Parameter sweeps](guide/12e-mtdna-large-studies.md#parameter-sweeps) |
| 12. Building models in Python | [2, A second model with two reactants](guide/02-first-model.md#a-second-model-with-two-reactants); [2, The builder calls](guide/02-first-model.md#the-builder-calls); [12, The model](guide/12-mtdna.md#the-model) |
| Larger models | [7, Larger models](guide/07-methods.md#larger-models) |
| 13. Model variants | [12a, Variants](guide/12a-mtdna-takeover.md#variants) |
| 14. Budgets, wall time and cancellation | [3, Budgets](guide/03-running.md#budgets); [12e, Budgets for long cells](guide/12e-mtdna-large-studies.md#budgets-for-long-cells) |
| 15. The command line and JSON requests | [12e, From the command line](guide/12e-mtdna-large-studies.md#from-the-command-line); [reference: contracts](reference/contracts.md) |
| 16. Rust | [developer: Rust API](developer/rust-api.md) |
| 17. Checking results | [12e, Checks behind these results](guide/12e-mtdna-large-studies.md#checks-behind-these-results) |
| 18. Reference | [reference: Python API](reference/python-api.md); [Migrating](guide/migrating.md); [reference: results](reference/results.md) |
| 19. Troubleshooting | [Troubleshooting](guide/troubleshooting.md) |
