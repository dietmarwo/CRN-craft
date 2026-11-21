# The model and the inference

## Six phases, 36 divisions

The model follows one cell lineage from the fertilized egg (time 0) to a
mature egg cell. Wild-type copies `W` and mutant copies `M` follow the
same reactions, because the paper's fitted model has no selection:

| Reaction | Rate per copy (per hour) |
|---|---|
| `W -> 2 W`, `M -> 2 M` (replication) | `λ_i` |
| `W -> ∅`, `M -> ∅` (degradation) | `ν_i` |

Here `i` is the phase. Times are in hours inside the code and in dpc in
the tables and plots.

| Phase | Time | Divisions | Free parameters |
|---|---|---|---|
| 1 | 0 to 203 h (8.46 dpc) | 29, every 7 h | `λ_1`, `ν_1` |
| 2 | 203 to 315 h (13.125 dpc) | 7, every 16 h | `λ_2`, `ν_2` |
| 3, 4, 5 | lengths `τ_3`, `τ_4`, `τ_5` in days | none | `λ_i`, `ν_i`, `τ_i` |
| 6 | until the end | none | none (`λ_6 = 0`; `ν_6 = 0` in the stored run) |

At each division the tracked daughter receives each copy with
probability one half: `W' ~ Binomial(W, 1/2)` and `M' ~ Binomial(M, 1/2)`.
An observation made exactly at a division time sees the daughter. The
initial copy number `m_0` is free, and the initial heteroplasmy is fixed
at `h_0 = 0.2`, as in the paper.

The example runs one lineage by joining CRN-craft segments
(`native_lineage` in `model.py`): it simulates a phase with
`static_direct`, draws the binomial partition in Python, and starts the
next segment from the sampled counts.

The published example also exposes `protocol_lineages` in the same module.
It uses one scheduled count run and reports the daughter immediately after
each partition. A comparison with the phase recurrence checks the mean
and variance through the ordered division times. [Scheduled experiments](../reference/experiments.md#scheduled-stochastic-simulation)


## Moments without simulation

The inference never simulates. Within a phase the expected copy number
grows as `exp((λ - ν) t)`, and a division halves it. So the log mean
copy number is exact and piecewise linear in time.

For the spread, let `z` be the squared coefficient of variation of the
total copy number. Within a phase `z` grows by `(λ + ν) ∫ dt / E(m)`.
A binomial division adds `1 / (2 E(m'))`, where `m'` is the daughter's
count. Both formulas are exact for this linear model. If wild-type and
mutant copies are independent, expanding `h = M / (W + M)` to first
order around its mean gives `V'(h) ≈ z`. The appendix of the paper also
truncates its heteroplasmy variance at first order. The review measured
the error of this approximation against exact simulation (see
[Results](21b-mtdna-bottleneck-results.md#how-good-is-the-first-order-variance)).

## The distance and the priors

The paper's appendix defines the distance between a parameter set `θ`
and the data `D` as the sum of two parts:

- `Σ_i [log E_θ(m, t_i) - log E_D(m, t_i)]²` over copy-number rows;
- `1000 Σ_j [V'_θ(h, t_j) - V'_D(h, t_j)]²` over heteroplasmy rows.

The first sum runs over the 34 copy-number rows of Table 1, the second
over its 40 heteroplasmy rows. The example uses natural logarithms; the
paper does not name the base. A parameter set is accepted when
`ρ ≤ 40`, the Figure 3 threshold. Parameter sets whose mean copy number
reaches 500,000 at any time are rejected, as in the paper. The priors
are uniform, as in the appendix:

| Parameter | Prior |
|---|---|
| `λ_i`, `ν_i`, phases 1–5 | each in [0, 1] per hour |
| `m_0` | (0, 1,000,000] copies |
| `τ_3`, `τ_4`, `τ_5` | each in [0, 50] days |

The code stores the net growth `λ_i - ν_i` and `ν_i`. This linear change
of variables has Jacobian one, so the prior stays uniform.

## Observation ages

Table 1 gives the earliest ages from the number of cells, using the 7 h
and 16 h division times (its footnote 1). It prints 7 h, 14 h and 21 h as
0.29, 0.58 and 0.88 dpc. Read literally, 0.29 dpc (6.96 h) lies just
before the first division. These six rows (three each from two studies)
count cells that have already divided, so the example places them at
7/24, 14/24 and 21/24 dpc, where the model sees the daughter. The
mature egg cells of the next generation have no age in Table 1. The paper
models them at 100 dpc, and so does the example.

## The sampler

The paper runs a Metropolis chain with small normal steps and accepts a
step when its distance is at most the threshold. The example targets the
same posterior, the prior restricted to `ρ ≤ 40`, with a different
proposal:

1. A differential evolution search (SciPy, seed 61 + ensemble) finds a
   starting point with `ρ ≤ 40`.
2. A single-coordinate random walk of 160,000 steps spreads 16 starting
   points around it.
3. Sixteen walkers then move by differential evolution Markov chain
   steps: walker `i` proposes `x_i + γ (x_a - x_b) + ε`, with `a` and `b`
   two other walkers, `γ` from a fixed cycle of eight values, and `ε`
   small normal noise. The proposal is symmetric, so the
   Metropolis–Hastings ratio is the prior ratio times the indicator. With
   a uniform prior the step is accepted exactly when it stays inside the
   prior box and below the threshold.

The first quarter of each ensemble is discarded. The plan's target was at
least 1,000 effective draws across four independent ensembles. For each
ensemble the effective size is computed for eight observables: log copy
number and `V'` at 10, 35, 70 and 100 dpc. Each uses the walker-averaged
series, with autocorrelation summed to the first negative pair, scaled by
the ratio of the variance over all walkers to the variance of that
average. The smallest of the eight counts per ensemble, and the four
counts are summed. The stored files do not show whether the run length
of 360,000 sweeps was fixed before the run.

## The independent mouse groups

Figure 4 of the paper reports 339 new single-cell heteroplasmy
measurements in 15 groups (9 litters, 24 to 61 dpc) from a mouse line
that carries two mtDNA types. The example recomputes each group's mean,
variance and `V'` from the cell values and checks them against the
spreadsheet's summary rows. It does not fit them. For each group it
predicts two ranges at the group's age:

- the 95% posterior range of the population `V'`;
- the 95% range of the sample `V'` of a group of that size. For this it
  draws cells from a zero-or-geometric founder law matched to the exact
  mean and the first-order variance, starting at the group's observed
  mean heteroplasmy.

The paper did something different. It fitted half of the groups (every
other group ordered by age) with a distance weight of 3,000. It then
compared transformed single-cell values of the other half with a 95%
range. Neither the transformation nor that split is reconstructed here.
All 15 groups are therefore a new validation of the Table 1 fit. The
manifest flags litters 1424 and 1406, but the fit uses no group, so the
flag does not separate training from test data.
