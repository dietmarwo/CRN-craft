# Python builder and whole-rule propensities: validation

Script: `validation/check_builder.py` (installed wheel; rebop 0.10.0 optional).
Reports: `validation/builder/builder-checks.json` (first run) and
`validation/builder/builder-checks-static-direct-defaults.json` (second run,
rules-run/0.6 defaults), each identical under `python -O`.

## Protocol (declared before the first run)

Every arm runs 4,000 replicates with `static_direct` on the installed wheel,
or 4,000 rebop 0.10.0 runs. Each test is a z test on a mean at a fixed time:
against the exact mean with the exact variance where the law is Poisson,
against the exact mean with the sample variance otherwise, and between two
arms with Welch's standard error sqrt(v1/n1 + v2/n2). One family of 18 tests
at family alpha 0.001, Bonferroni: |z| <= z* = 4.03 (the 1 - 0.001/36
standard normal quantile). Rows that rebop cannot run (not installed) are
recorded as skipped and still count in the family size. No test is changed
after a run; failed rows are kept.

| # | Target | Arms (root seed) | Exact value or comparison |
|---|---|---|---|
| 1 | Birth–death, numeric rates k = 10, g = 1, A(0) = 0: mean A(10) | builder (201) | Poisson, 10 (1 − e^−10) |
| 2 | Same with string propensities `"k"`, `"g * A"` | builder (202) | same |
| 3 | Single-enzyme Michaelis–Menten with numeric rates: 5 enzymes, binding kon = 2 (pseudo-first order), koff = 0.5, kcat = 1: mean P(20) | builder (203) | 5 kcat π (t − (1 − e^−λt)/λ), λ = kon + koff + kcat, π = kon/λ |
| 4 | Michaelis–Menten flux as a string propensity `"vmax * S / (Km + S)"` on `S -> S + P`, S = 50, vmax = 4, Km = 25: mean P(10) | builder (204) | Poisson, 10 vmax S/(Km + S) |
| 5–10 | SIR, S(0) = 490, I(0) = 10, beta = 5e-4, gamma = 0.1: mean S and I at t = 10, 20, 40 | hand-written mass-action JSON (101) vs builder with string propensities `"beta * S * I"`, `"gamma * I"` (102) | Welch |
| 11–16 | Same SIR | hand-written JSON (101) vs rebop with numeric rates (rng 1001) | Welch |
| 17–18 | Dimerization A + A -> 0, A(0) = 100: mean A at t = 0.5, 1 | builder with numeric rate 2k = 0.02 (301) vs rebop with numeric rate k = 0.01 (rng 1002) | Welch; checks the documented factor 2 for identical reactants |

Also recorded as checks without a statistical test: the builder's JSON for
the SIR and birth–death models equals the hand-written documents, and a
builder model passes `rules.load`'s validation.

## Results

First run, 2026-11-06, on the 2-vCPU development container with rebop
0.10.0. **18 of 18 tests and 3 of 3 checks passed; no failed row to
retain.** The largest |z| is 2.62 against z* = 4.03.

| # | Test | Estimate | Expected or other arm | Standard error | z | Outcome |
|---|---|---|---|---|---|---|
| 1 | `birth_death_numeric_mean_A_10` | 9.9665 | 9.9995 | 0.0500 | -0.66 | pass |
| 2 | `birth_death_string_mean_A_10` | 10.1305 | 9.9995 | 0.0500 | +2.62 | pass |
| 3 | `michaelis_menten_mass_action_mean_P_20` | 56.3627 | 56.3265 | 0.0995 | +0.36 | pass |
| 4 | `michaelis_menten_propensity_mean_P_10` | 26.5963 | 26.6667 | 0.0816 | -0.86 | pass |
| 5 | `sir_json_vs_builder_S_10` | 440.2482 | 440.2058 | 0.4282 | +0.10 | pass |
| 6 | `sir_json_vs_builder_I_10` | 38.0868 | 38.2302 | 0.3293 | -0.44 | pass |
| 7 | `sir_json_vs_builder_S_20` | 323.0242 | 321.9948 | 1.1098 | +0.93 | pass |
| 8 | `sir_json_vs_builder_I_20` | 91.1787 | 91.8852 | 0.6105 | -1.16 | pass |
| 9 | `sir_json_vs_builder_S_40` | 115.6617 | 115.0447 | 0.8035 | +0.77 | pass |
| 10 | `sir_json_vs_builder_I_40` | 87.9950 | 87.5242 | 0.3557 | +1.32 | pass |
| 11 | `sir_json_vs_rebop_S_10` | 440.2482 | 439.9885 | 0.4314 | +0.60 | pass |
| 12 | `sir_json_vs_rebop_I_10` | 38.0868 | 38.4008 | 0.3330 | -0.94 | pass |
| 13 | `sir_json_vs_rebop_S_20` | 323.0242 | 321.9558 | 1.1154 | +0.96 | pass |
| 14 | `sir_json_vs_rebop_I_20` | 91.1787 | 91.8782 | 0.6127 | -1.14 | pass |
| 15 | `sir_json_vs_rebop_S_40` | 115.6617 | 115.1573 | 0.8065 | +0.63 | pass |
| 16 | `sir_json_vs_rebop_I_40` | 87.9950 | 87.5355 | 0.3531 | +1.30 | pass |
| 17 | `dimer_builder_2k_vs_rebop_k_A_0.5` | 50.1440 | 50.1925 | 0.1197 | -0.41 | pass |
| 18 | `dimer_builder_2k_vs_rebop_k_A_1` | 33.4950 | 33.5565 | 0.1034 | -0.59 | pass |

Rows 17–18 show the documented convention: a builder rate of 2k on `A + A`
reproduces rebop's k. With rate k instead, A(1) would be near
100/(1 + 100·0.01) = 50 rather than 33.5 (deterministic limit), a difference
of more than a hundred standard errors; this alternative was not run.

## Second run: rules-run/0.6 Static Direct defaults

Since rules-run/0.6, `method="static_direct"` in Python runs Ziggurat waits
and sparse trackers (no trackers here), so the same protocol and seeds draw
other samples. Rerun 2026-11-06 on the 2-vCPU development container with
rebop 0.10.0, unchanged in tests, sizes, seeds and alpha. **18 of 18 tests
and 3 of 3 checks passed.** The largest |z| is 1.54 against z* = 4.03.
The first run stays above as recorded; the explicit options
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}`
reproduce its samples.

| # | Test | Estimate | Expected or other arm | Standard error | z | Outcome |
|---|---|---|---|---|---|---|
| 1 | `birth_death_numeric_mean_A_10` | 10.0585 | 9.9995 | 0.0500 | +1.18 | pass |
| 2 | `birth_death_string_mean_A_10` | 10.0698 | 9.9995 | 0.0500 | +1.40 | pass |
| 3 | `michaelis_menten_mass_action_mean_P_20` | 56.2602 | 56.3265 | 0.0979 | -0.68 | pass |
| 4 | `michaelis_menten_propensity_mean_P_10` | 26.5587 | 26.6667 | 0.0816 | -1.32 | pass |
| 5 | `sir_json_vs_builder_S_10` | 440.4183 | 440.0230 | 0.4278 | +0.92 | pass |
| 6 | `sir_json_vs_builder_I_10` | 38.1003 | 38.4633 | 0.3325 | -1.09 | pass |
| 7 | `sir_json_vs_builder_S_20` | 323.1363 | 321.6965 | 1.1079 | +1.30 | pass |
| 8 | `sir_json_vs_builder_I_20` | 91.2095 | 92.0060 | 0.6073 | -1.31 | pass |
| 9 | `sir_json_vs_builder_S_40` | 115.2625 | 115.0720 | 0.7953 | +0.24 | pass |
| 10 | `sir_json_vs_builder_I_40` | 87.8738 | 87.3298 | 0.3539 | +1.54 | pass |
| 11 | `sir_json_vs_rebop_S_10` | 440.4183 | 439.9885 | 0.4335 | +0.99 | pass |
| 12 | `sir_json_vs_rebop_I_10` | 38.1003 | 38.4008 | 0.3347 | -0.90 | pass |
| 13 | `sir_json_vs_rebop_S_20` | 323.1363 | 321.9558 | 1.1148 | +1.06 | pass |
| 14 | `sir_json_vs_rebop_I_20` | 91.2095 | 91.8782 | 0.6103 | -1.10 | pass |
| 15 | `sir_json_vs_rebop_S_40` | 115.2625 | 115.1573 | 0.7915 | +0.13 | pass |
| 16 | `sir_json_vs_rebop_I_40` | 87.8738 | 87.5355 | 0.3538 | +0.96 | pass |
| 17 | `dimer_builder_2k_vs_rebop_k_A_0.5` | 50.3335 | 50.1925 | 0.1226 | +1.15 | pass |
| 18 | `dimer_builder_2k_vs_rebop_k_A_1` | 33.5490 | 33.5565 | 0.1048 | -0.07 | pass |

Limitations: the tests compare means at three or fewer times, not whole
distributions; the exact-law gates of the rules engine
(`validation/rules/RESULTS.md`) cover distributions for rate rules, and the
Rust tests check that `propensity: "k * A"` reproduces `rate: "k"` path by
path under every method.
