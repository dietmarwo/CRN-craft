# Summary

[CRN-craft](index.md)
[How do I ...?](howto.md)
[How CRN-craft compares](comparison.md)

# User guide

- [1. Install](guide/01-install.md)
- [2. A first model](guide/02-first-model.md)
- [3. Running ensembles](guide/03-running.md)
- [4. Results](guide/04-results.md)
- [5. Speed: threads, parameter sweeps and compiled models](guide/05-speed.md)
- [6. Stop conditions, trackers and event logs](guide/06-stops-trackers-logs.md)
- [7. Choosing a method](guide/07-methods.md)
- [8. Reproducing runs and published results](guide/08-reproducibility.md)
- [9. Continuous variables](guide/09-continuous-variables.md)
- [10. The command line, request files and archives](guide/10-command-line.md)
- [11. Analysis, import and export](guide/11-analysis-import-export.md)
- [12. Case study: mtDNA mutants](guide/12-mtdna.md)
  - [Takeover in a fixed population](guide/12a-mtdna-takeover.md)
  - [Mutant classes: trackers, event logs and callbacks](guide/12b-mtdna-events.md)
  - [Variable population: ATP](guide/12c-mtdna-atp.md)
  - [Methods and speed](guide/12d-mtdna-methods.md)
  - [Large studies: sweeps, budgets and the command line](guide/12e-mtdna-large-studies.md)
- [13. Fitting parameters to data](guide/13-fitting.md)
- [14. Model files and species with attributes](guide/14-model-files.md)
  - [The mtDNA model files: clones, the fast path and stored models](guide/14a-mtdna-model-files.md)
- [15. Using CRN-craft from R](guide/15-from-r.md)
- [16. Case study: re-creating a published study](guide/16-published-study.md)
  - [Fig 1A: COX-negative cells over time](guide/16a-fig1a.md)
  - [Fig 1B: mutation rates for 5, 10 and 15 %](guide/16b-fig1b.md)
  - [Fig 1C/1D: clones in COX-negative cells](guide/16c-fig1cd.md)
  - [Limits and lessons](guide/16d-limits-lessons.md)
- [17. Case study: comparing mtDNA mechanisms](guide/17-mtdna-mechanisms.md)
  - [The three models](guide/17a-mtdna-mechanisms-models.md)
  - [Fitting and the comparison](guide/17b-mtdna-mechanisms-results.md)
  - [Choices, checks and limits](guide/17c-mtdna-mechanisms-limits.md)
- [18. Case study: inferring mtDNA takeover times](guide/18-mtdna-inference.md)
  - [The model and Figure 6](guide/18a-mtdna-inference-model.md)
  - [Fitting snapshots and takeover times](guide/18b-mtdna-inference-results.md)
  - [Choices, checks and limits](guide/18c-mtdna-inference-limits.md)
- [19. Spatial spread of mtDNA deletions](guide/19-mtdna-spread.md)
- [20. Mitochondrial networks and genetic drift](guide/20-mtdna-network.md)
- [21. The mtDNA bottleneck during development](guide/21-mtdna-bottleneck.md)
- [22. Protein distributions from four reactions](guide/22-two-stage-expression.md)
- [23. Fitting single-cell RNA distributions](guide/23-yeast-rna-fit.md)
- [Troubleshooting](guide/troubleshooting.md)
- [Migrating from earlier releases](guide/migrating.md)

# Reference

- [Model language](reference/model-language.md)
- [Methods and their options](reference/methods.md)
- [Python API](reference/python-api.md)
- [Contracts and versions](reference/contracts.md)
- [Results](reference/results.md)
- [Archives](reference/archives.md)
- [Diagnostics](reference/diagnostics.md)
- [Analysis report](reference/analysis.md)
- [GillesPy2 export](reference/gillespy2-export.md)

# Developer

- [Architecture](developer/architecture.md)
- [Rust API](developer/rust-api.md)
- [Extending and checks](developer/extending.md)
- [Specialized kernels and generated loops](developer/specialized-kernels.md)
- [Context for AI agents](developer/ai-context.md)
- [Releasing](developer/releasing.md)

# Evidence

- [Technical report](evidence/technical-report.md)
- [mtDNA performance analysis](evidence/mtdna-performance-analysis.md)
- [Comparison with GillesPy2](evidence/comparison-gillespy2.md)
- [Comparison with Catalyst and JumpProcesses](evidence/comparison-catalyst-jumpprocesses.md)
  - [The comparison before acceleration](evidence/comparison-catalyst-jumpprocesses-before-acceleration.md)
- [Ten further Julia-gap iterations](evidence/comparison-julia-gap10.md)
- [Julia at one and eight threads](evidence/comparison-julia-8-threads.md)
- [Parameter optimization against rebop](evidence/comparison-parameter-optimization.md)
- [Model and algorithm cards](evidence/model-cards/README.md)
  - [Direct SSA](evidence/model-cards/direct.md)
  - [Rate laws](evidence/model-cards/rate-laws.md)

# History

- [Plans, decisions and campaign logs](history/index.md)
