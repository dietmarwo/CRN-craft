# Tutorials

The tutorials come in two groups. Pick the group by what you want to do.

- **[Simulation tutorials](#simulation-tutorials)** build a model from a
  paper or a question, simulate it and compare the outcome with a published
  figure, an exact law or archived data. No optimizer is needed to follow
  them.
- **[Optimization tutorials](#optimization-tutorials)** put an optimizer
  from [fcmaes](https://github.com/dietmarwo/fast-cma-es) around a
  simulation: they fit parameters to data, select models, design
  experiments or search for a trade-off, and then check what the optimizer
  returned.

The [user guide](guide/01-install.md), chapters 1 to 11, 14 and 15,
explains the interface that both groups use. Read chapters 1 to 4 first.

Every tutorial is runnable from this repository, saves the tables behind
its figures, and says which of its numbers agree with the source, which
are inconclusive and which disagree. A tutorial that misses a published
value reports the miss; none is tuned until it matches.

## Simulation tutorials

| Tutorial | Source | What it does | Outcome |
| --- | --- | --- | --- |
| [12. mtDNA mutants](guide/12-mtdna.md) | a teaching model | Takeover of a cell by mutant mtDNA: ensembles, trackers, event logs, a continuous ATP variable, methods and large studies | a worked example, not a reproduction |
| [16. Random drift of mtDNA mutants](guide/16-published-study.md) | Kowald & Kirkwood 2013, *Aging Cell* | Re-creates Figure 1A–D: COX-negative cells over 120 years, the mutation rate for 5, 10 and 15%, clones per cell | contributed by an author of the paper |
| [17. Three mtDNA mechanisms](guide/17-mtdna-mechanisms.md) | Kowald & Kirkwood 2014, *PNAS*; Kowald, Dawson & Kirkwood 2014, *J. Theor. Biol.* | Drift, size advantage and transcription-linked replication, each calibrated to 10% affected cells, then compared by mutant lineages per cell | 1.26 lineages at three years, as published; two unstated model points confirmed by an author |
| [19. Spatial spread of mtDNA deletions](guide/19-mtdna-spread.md) | Insalata et al. 2022, *PNAS* | Two coupled regions and a traveling wave along a muscle fiber; an extension with de novo deletions | Figure 1C agrees; the wave speed is an unresolved discrepancy (0.13–0.16 against 0.190 µm/day) |
| [20. Mitochondrial networks and drift](guide/20-mtdna-network.md) | Aryaman et al. 2019, *Genetics* | Variance of the mutant fraction under fusion and fission, against Equation 13 and the archived simulations | agreement in 38 of 39 conditions, the 39th after a rerun |
| [21. The mtDNA bottleneck](guide/21-mtdna-bottleneck.md) | Johnston et al. 2015, *eLife* | Six developmental phases with divisions; approximate Bayesian computation against copy numbers and variances | Figure 3 reproduced up to 50 days; the copy number at day 100 is 23 times below the published curve |
| [22. Protein distributions](guide/22-two-stage-expression.md) | Shahrezaei & Swain 2008, *PNAS* | Exact stationary laws of two-stage expression against the negative-binomial approximation | Figure 1D matches if its unit is bits, not nats; the paper does not state the unit |

## Optimization tutorials

Three chapters of the guide fit models to data with fcmaes:

| Tutorial | Source | What is optimized | Outcome |
| --- | --- | --- | --- |
| [13. Fitting parameters to data](guide/13-fitting.md) | the model of chapter 12 | Two rates of a stochastic model against single-cell data: Monte Carlo noise, common random numbers, staged refinement, held-out checks | a worked example |
| [18. Inferring mtDNA takeover times](guide/18-mtdna-inference.md) | Kowald & Kirkwood 2026, *npj Aging* | Mutation probability, selection advantage and advantageous fraction from cell snapshots | every true value inside its 95% interval; the baseline intervals are narrower than published, cause unknown |
| [23. Fitting single-cell RNA distributions](guide/23-yeast-rna-fit.md) | Munsky et al. 2018, *PNAS* | A four-state gene model against 59,557 yeast cells, by means and by whole histograms | the published scores reproduce; the searches start near the published vector, so they test the route, not a blind fit |

Ten lessons live with their code in
[`tutorials/optimization`](../tutorials/optimization/README.md) and run
through one command, `python tutorials/optimization/run.py C1 … C10`.
[Chapter 24](guide/24-publication-tutorials.md) introduces them. C1 to C4
teach the route on synthetic data; C5 to C10 each reconstruct a published
optimization study.

| Lesson | Source | What is optimized | Outcome |
| --- | --- | --- | --- |
| [C1. Calibration](../tutorials/optimization/C1-calibration.md) | synthetic data | Rates of a stochastic model, with reusable preparation | method lesson |
| [C2. Sensitivities](../tutorials/optimization/C2-sensitivities.md) | synthetic data | Rates of an ODE model with derivatives, then transferred to the stochastic model | method lesson |
| [C3. Oscillators](../tutorials/optimization/C3-oscillators.md) | synthetic circuits | Several objectives and a diversity archive over oscillating circuits | method lesson |
| [C4. Rules](../tutorials/optimization/C4-rules.md) | exact aggregation laws | Kernels of a rule-based model from finite samples | method lesson |
| [C5. RNA count calibration](../tutorials/optimization/C5-rna-counts.md) | Poovathingal & Gunawan 2010, *BMC Bioinformatics* | Three rates from count distributions | the three rates are recovered within 3% from 10,000 trajectories; from 100 trajectories two of four data sets miss one rate by a factor 1.7 along a flat valley of the score |
| [C6. Stimulus model selection](../tutorials/optimization/C6-stimulus-fit.md) | Gómez-Schiavon et al. 2017, *Genome Biology* | Seven hypotheses fitted to 376 measured cells, ranked by AIC and BIC | both criteria select the published model; the fitted rates lie within one posterior standard deviation of the paper's |
| [C7. Virtual reactor design](../tutorials/optimization/C7-reactor.md) | Felton et al. 2021, *Chemistry–Methods* | Yield and waste of a reactor in 50 evaluations | at 50 evaluations MODE is no better than random search and far below the paper's best strategy; with 500 it passes that strategy's median |
| [C8. Sampling design](../tutorials/optimization/C8-experiment-design.md) | Fox & Munsky 2019, *PLoS Comput. Biol.* | Sampling times and a stimulus pulse by Fisher information | the sampling optimum and the ranking of the two published designs in 97 of 100 parameter samples are reproduced, the archived matrices entry by entry |
| [C9. STAT5 calibration](../tutorials/optimization/C9-stat5.md) | Boehm et al. 2014, *J. Proteome Res.* | Nine parameters of a signaling model against a PEtab problem | the native model reproduces the archived likelihood within 2e-6 nats; a search with profiled noise scales reaches it, plain retries stop 9 nats above |
| [C10. Pathway calibration](../tutorials/optimization/C10-pathway.md) | Moles, Mendes & Banga 2003, *Genome Research* | Thirty-six parameters over sixteen experiments | the objective is verified; searches of up to 20,000 evaluations reach costs of 62 to 667 against a published 0.0013 after 2.8 million |

## How the tutorials were made

[Chapter 16](guide/16-published-study.md) was contributed by an author of
the paper it re-creates. The other tutorials were written by AI agents
directed by the maintainer: an agent read the paper, rebuilt the model,
ran it and compared the numbers; other agents then reviewed the work
without having seen it made. Unless a chapter says otherwise, the papers'
authors did not write or review it. The [README](../README.md#what-it-is-for)
describes this way of working and the two others the project is built for.
