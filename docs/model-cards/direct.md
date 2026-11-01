# Direct SSA

**Status:** run contract only; planned M2, not implemented.

Reference: D. T. Gillespie, *Exact stochastic simulation of coupled chemical reactions* (1977), [DOI](https://doi.org/10.1021/j100540a008).

The future baseline sums finite non-negative propensities, draws a positive exponential wait and selects a reaction proportional to its propensity. Zero total propensity holds the state through the final save time. Count updates are checked. Inputs must define a time-homogeneous well-mixed count process; budgets prevent uncontrolled work but do not prove mathematical non-explosion.

Implementation must specify uniform endpoints, selection rounding, time-progress errors and stream mapping. Before any speed claim: sampler known-answer vectors, analytic marginal checks, save-time boundary tests, reactant/conservation checks, 1/2/8-thread numerical reproducibility and budget/cancellation failure tests. Record dynamic-range loss rather than silently assigning missing probability to the final channel.
