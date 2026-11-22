# Reproducing published optimization studies

The six publication studies join a physical or stochastic model to a
bounded optimization objective. Each uses an external paper or archived
calculation to check numerical results before running the optimizer.

| Study | Optimization objective | Comparison |
|---|---|---|
| RNA counts | Weighted cumulative-distribution discrepancy | Simulation-based recovery and projected count probabilities |
| Stimulus model | Joint count likelihood and model selection | Four observed snapshots and seven hypotheses |
| Virtual reactor | Space–time yield and E-factor | Fifty proposals per run and two-dimensional hypervolume |
| Experimental design | Smallest Fisher-matrix eigenvalue | A sampling-interval sweep and a pulse design |
| STAT5 observations | Gaussian likelihood with fitted noise | The archived PEtab predictions and measurement tables |
| Three-step pathway | Weighted residual sum across conditions | Sixteen independent simulated trajectories |

The [publication tutorial index](https://github.com/dietmarwo/CRN-craft/tree/main/tutorials/optimization#publication-studies-c5c10)
links to the runnable lessons and fixture setup.

For timed designs, `crn_craft.sensitivities.forward_sensitivities` advances
the law and its forward derivatives together, including the FSP sink.
A prepared graph can take a `retained` set of count tuples in lowered species
order. Transitions outside that set contribute to the sink.

For concentration trajectories, the optional `crn-diffsol concentrations`
command prepares named formulas once and reads candidate batches over
standard input. The C7, C9 and C10 models provide full working examples.
