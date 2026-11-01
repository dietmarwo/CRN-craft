# Statistical verification policy, recorded before simulation

Date: 1 November 2026. M0 has no stochastic output, and no statistical solver tests have been run.

M2 fixtures will predeclare model bytes, initial laws, times, marginals, independent replicate count, fixed seeds, hypotheses and family-wise alpha before first execution. Initial analytic CDF target: n=100000 per fixture, suite alpha=0.001, DKW with Bonferroni allocation across declared tests. For m tests, epsilon = sqrt(log(2m/alpha)/(2n)). A union bound permits dependence between times/species. Analytical distribution errors and FSP numerical errors are separate terms.

Pure decay uses a binomial law. Immigration–death from zero uses a Poisson law. Nonzero deterministic initial populations add surviving initial particles; do not substitute a Poisson law indiscriminately. Zero-variance cases require exact/appropriate handling.

DSMTS settings and Z/Y values are retained in M3. No binomial model of the count of correlated pointwise failures is authorized. Freeze a justified simultaneous or calibrated family rule before the campaign. Fixed seeds make a regression repeatable; rotating seeds assess false-alarm behavior without changing the gate after a failure.

Each statistic names its sample count, estimator and uncertainty method. Only independent replicate-level samples count toward ensemble SE. Variance and quantiles do not inherit a mean-SE formula. Censored first-passage times stay censored. Partial completion is tracked and cannot silently change the sampling population.

For backend comparisons, predeclare the hypothesis family and combined uncertainty method. For speed, pair repeated comparable workloads and report variability, actual work and environment. For oscillation classification, freeze held-out topologies, parameters, seed sets, burn-in and horizon. For tau-leaping, assess bias relative to uncertainty, not noisy pointwise monotonicity.

Keep failed rows. Revisions to a test's mathematical validity require a documented decision and a newly declared campaign; never quietly loosen a threshold to pass.
