# Summary

[CRN-craft](index.md)
[Tutorials](tutorials.md)
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
- [12. Model files and species with attributes](guide/12-model-files.md)
  - [The mtDNA model files: clones, the fast path and stored models](guide/12a-mtdna-model-files.md)
- [13. Using CRN-craft from R](guide/13-from-r.md)
- [Troubleshooting](guide/troubleshooting.md)
- [Migrating from earlier releases](guide/migrating.md)

# Simulation tutorials

- [S1. Case study: mtDNA mutants](tutorials/simulation/S1-mtdna.md)
  - [Takeover in a fixed population](tutorials/simulation/S1a-mtdna-takeover.md)
  - [Mutant classes: trackers, event logs and callbacks](tutorials/simulation/S1b-mtdna-events.md)
  - [Variable population: ATP](tutorials/simulation/S1c-mtdna-atp.md)
  - [Methods and speed](tutorials/simulation/S1d-mtdna-methods.md)
  - [Large studies: sweeps, budgets and the command line](tutorials/simulation/S1e-mtdna-large-studies.md)
- [S2. Case study: re-creating a published study](tutorials/simulation/S2-published-study.md)
  - [Fig 1A: COX-negative cells over time](tutorials/simulation/S2a-fig1a.md)
  - [Fig 1B: mutation rates for 5, 10 and 15 %](tutorials/simulation/S2b-fig1b.md)
  - [Fig 1C/1D: clones in COX-negative cells](tutorials/simulation/S2c-fig1cd.md)
  - [Limits and lessons](tutorials/simulation/S2d-limits-lessons.md)
- [S3. Case study: comparing mtDNA mechanisms](tutorials/simulation/S3-mtdna-mechanisms.md)
  - [The three models](tutorials/simulation/S3a-mtdna-mechanisms-models.md)
  - [Fitting and the comparison](tutorials/simulation/S3b-mtdna-mechanisms-results.md)
  - [Choices, checks and limits](tutorials/simulation/S3c-mtdna-mechanisms-limits.md)
- [S4. Spatial spread of mtDNA deletions](tutorials/simulation/S4-mtdna-spread.md)
  - [Choices, checks and limits](tutorials/simulation/S4a-mtdna-spread-limits.md)
  - [Extension: de novo deletions](tutorials/simulation/S4b-mtdna-spread-de-novo.md)
  - [Extension: results and limits](tutorials/simulation/S4c-mtdna-spread-de-novo-results.md)
- [S5. Mitochondrial networks and genetic drift](tutorials/simulation/S5-mtdna-network.md)
  - [Results, choices and limits](tutorials/simulation/S5a-mtdna-network-results.md)
- [S6. The mtDNA bottleneck during development](tutorials/simulation/S6-mtdna-bottleneck.md)
  - [The model and the inference](tutorials/simulation/S6a-mtdna-bottleneck-model.md)
  - [Results](tutorials/simulation/S6b-mtdna-bottleneck-results.md)
  - [Choices, checks and limits](tutorials/simulation/S6c-mtdna-bottleneck-limits.md)
- [S7. Protein distributions from four reactions](tutorials/simulation/S7-two-stage-expression.md)

# Optimization tutorials

- [O1. Fitting parameters to data](tutorials/optimization/O1-fitting.md)
- [O2. Case study: inferring mtDNA takeover times](tutorials/optimization/O2-mtdna-inference.md)
  - [The model and Figure 6](tutorials/optimization/O2a-mtdna-inference-model.md)
  - [Fitting snapshots and takeover times](tutorials/optimization/O2b-mtdna-inference-results.md)
  - [Choices, checks and limits](tutorials/optimization/O2c-mtdna-inference-limits.md)
- [O3. Fitting single-cell RNA distributions](tutorials/optimization/O3-yeast-rna-fit.md)
  - [Data and model](tutorials/optimization/O3a-yeast-rna-fit-model.md)
  - [Fits and predictions](tutorials/optimization/O3b-yeast-rna-fit-results.md)
  - [Choices, checks and limits](tutorials/optimization/O3c-yeast-rna-fit-limits.md)
- [Ten optimization lessons with fcmaes](tutorials/optimization/lessons.md)

# Reference

- [Model language](reference/model-language.md)
- [Methods and their options](reference/methods.md)
- [Reproducible experiments](reference/experiments.md)
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
