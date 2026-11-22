# Ten optimization lessons with fcmaes

Tutorials [O1](O1-fitting.md), [O2](O2-mtdna-inference.md) and
[O3](O3-yeast-rna-fit.md) fit one model each. The ten lessons of
[`tutorials/optimization`](../../../tutorials/optimization/README.md) cover the
other things an optimizer does with a simulation: calibration with
reusable preparation, derivatives, several objectives, model selection,
experiment design and fits over many conditions. They run through one
command and keep their code, saved tables and figures next to the lesson
pages.

```sh
python tutorials/optimization/run.py C7 --output /tmp/crn-c7-smoke
```

Every lesson has a `smoke` preset that runs the whole route in seconds and
a `study` preset with the budgets of its page. A run writes the model, the
numerical checks, every proposal the optimizer made and a `summary.json`.
`passed` in that file means that the lesson's declared checks passed. It
does not mean that the optimizer found the true parameters.

## Method lessons

| Lesson | Question |
| --- | --- |
| [C1. Calibration](../../../tutorials/optimization/C1-calibration.md) | Does reusable preparation preserve a complete fit? |
| [C2. Sensitivities](../../../tutorials/optimization/C2-sensitivities.md) | Are derivatives correct, rates identifiable and ODE fits transferable to the stochastic model? |
| [C3. Oscillators](../../../tutorials/optimization/C3-oscillators.md) | Which trade-offs and behaviors survive independent paths? |
| [C4. Rules](../../../tutorials/optimization/C4-rules.md) | Can finite-path fits of exact aggregation laws tell kernels apart? |

## Publication studies

Each study checks its objective against the paper or the paper's archive
before any search, then runs fcmaes and verifies the finalists.

| Lesson | Paper | Objective | Outcome |
| --- | --- | --- | --- |
| [C5. RNA count calibration](../../../tutorials/optimization/C5-rna-counts.md) | Poovathingal & Gunawan 2010, [doi:10.1186/1471-2105-11-414](https://doi.org/10.1186/1471-2105-11-414) | weighted distance of cumulative count distributions | recovered within 3% from 10,000 trajectories; a flat valley of the score limits recovery from 100 |
| [C6. Stimulus model selection](../../../tutorials/optimization/C6-stimulus-fit.md) | Gómez-Schiavon et al. 2017, [doi:10.1186/s13059-017-1297-9](https://doi.org/10.1186/s13059-017-1297-9) | joint count likelihood of seven hypotheses, AIC and BIC | the published model is selected |
| [C7. Virtual reactor design](../../../tutorials/optimization/C7-reactor.md) | Felton et al. 2021, [doi:10.1002/cmtd.202000051](https://doi.org/10.1002/cmtd.202000051) | space–time yield against E-factor, by hypervolume | at 50 evaluations no better than random search; ahead of the paper's best strategy with 500 |
| [C8. Sampling design](../../../tutorials/optimization/C8-experiment-design.md) | Fox & Munsky 2019, [doi:10.1371/journal.pcbi.1006365](https://doi.org/10.1371/journal.pcbi.1006365) | smallest eigenvalue of the Fisher information | the sampling optimum and the 97-of-100 design ranking are reproduced |
| [C9. STAT5 calibration](../../../tutorials/optimization/C9-stat5.md) | Boehm et al. 2014, [doi:10.1021/pr5006923](https://doi.org/10.1021/pr5006923) | Gaussian likelihood with fitted noise scales | the archived likelihood is reproduced; a search with profiled noise scales reaches it in three of five seeds |
| [C10. Pathway calibration](../../../tutorials/optimization/C10-pathway.md) | Moles, Mendes & Banga 2003, [doi:10.1101/gr.1262503](https://doi.org/10.1101/gr.1262503) | weighted residuals over 16 experiments, 36 parameters | objective verified; the searches stay far from the published value at 0.7% of its budget |

Three of these outcomes are misses or partial results, and the lessons say
why. In C7 the paper's protocol of 50 evaluations leaves an evolutionary
optimizer four selection steps. In C10 the published value took 2.8
million evaluations; a local least-squares fit from a perturbed start
reaches the integration floor, so the objective is right and the gap is
budget. In C5 the generating rates are recovered when the search runs long
enough and the data are plentiful, and the remaining error at 100
trajectories follows a valley of the objective.

## Engine features the studies use

- **Probability sensitivities.**
  `crn_craft.sensitivities.forward_sensitivities` advances a projected law
  and its derivatives with respect to parameters together, including the
  sink, and `fisher_information` forms the Fisher matrix over observed
  bins ([reference](../../reference/experiments.md#probability-sensitivities-and-fisher-information)).
  C8 uses them.
- **Retained state sets.** `prepare_fsp(model, bounds, retained=states)`
  keeps a chosen set of count tuples; transitions to other states go to
  the sink
  ([reference](../../reference/experiments.md#retained-state-sets)).
- **Concentration trajectories.** The optional `crn-diffsol concentrations`
  command prepares a deterministic model with real-valued states, compartment
  volumes and time-dependent inputs once and then solves batches of
  candidates read from standard input
  ([interface](../../../integrations/diffsol/CONCENTRATIONS.md)). C7, C9 and
  C10 use it.
