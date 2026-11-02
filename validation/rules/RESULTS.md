# Rules engine K1: exact-law gates

Protocol: `docs/decisions/rules-k0-k1-2026-11-02.md` (declared at `fe0ca4e`,
implementation details and power laws fixed at `1c5df66` before the run).
Script: `validation/check_rules.py`. Report: `k1-gates.json` (fixed-work
evidence; identical under `python -O`). Root seed 20261102, 100,000 replicates
per fixture, four worker threads.

**Result: 47 of 47 checks passed, no failure retained.**

| Fixture | Exact law | Marginals | Largest ECDF distance |
|---|---|---|---|
| Immigration–death (k_in 10, k_out 1) | Poisson, mean 10(1 − e^{−t}) | X at t = 0.5, 1, 3 | 0.00161 |
| Immigration–dimerization (k_in 20, k2 0.1), rules engine | master equation on 0–120, leakage 4.8e−106 | X at t = 0.5, 1, 3 | 0.00377 |
| Same model, network engine (`ssa`, Direct) | same | X at t = 0.5, 1, 3 | 0.00198 |
| mtDNA, fixed population (K 20, δ 1, sw 1, sm 1.5, μ 0.1, fracAdv 0.5) | master equation on (wt, sw mutants, sm mutants), 1,771 states | wild types, mutants, advantaged mutants at t = 1, 3, 10 | 0.00322 |
| Aggregation, constant kernel, 12 monomers | master equation on the 77 partitions of 12 | clusters, largest cluster at t = 0.02, 0.1, 0.5 | 0.00256 |
| Aggregation, multiplicative kernel, 12 monomers | same | same | 0.00306 |

DKW bound for one family of 30 comparisons at alpha 0.001:
ε = sqrt(ln(2·30/0.001)/(2·100,000)) = 0.00742. The largest distance over all
thirty marginals is 0.00377, about half the bound.

Power: every wrong law is rejected by the same samples (largest distance):
immigration–death k_in = 10.5 (0.059), dimerization k2 = 0.11 (0.084), mtDNA
sm = 1.35 (0.062), aggregation k = 1.1 (0.085 constant, 0.098 multiplicative).
The supplementary mtDNA check also rejects the law with the "natural"
competition term D = w + sw·a + sm·b (0.055), so the engine reproduces
Kowald's D = w + sm·(a + b) specifically, not merely an mtDNA-like model.

Other checks: four schema snapshots equal the CLI's generated schemas; each
model passes `rules-check` and validates against the model schema; every
request and result validates against its schema; a 2,000-replicate ensemble of
the multiplicative kernel is identical with one and with four threads.

Run time on the 2-vCPU development container: about two minutes for the whole
script (rules runs, network run, exact laws, JSON parsing); not a benchmark.

Limitations: the fixtures are small (at most 20 molecules or 12 monomers)
because their laws are computed exactly; large populations, many clones and
long horizons are covered by unit tests only for accounting, not for law. The
DKW bound is conservative for discrete laws. Validation does not cover
time-dependent rates (not implemented) or performance.
