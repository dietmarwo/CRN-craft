# Exact methods: Static Direct, selectors and next reaction

This page describes the exact simulation methods of the one engine, the event
selectors of Static Direct and the random stream they share. The model
language is in the [model reference](rules-language.md), the other methods
(`tree`, `mtdna_direct`, `tau`, `binomial`, `thinning`, `hybrid`) are in
[mtDNA ensembles](mtdna-rules.md), and the cost of each is in the
[performance analysis](mtdna-performance-analysis.md). An exact method samples
the model's process with finite-precision arithmetic and Monte Carlo
uncertainty. It is not a biological validation claim.

## The methods

| `method` | Exact | Scope | Sampler identity (prefix `chacha8-u53/` or `chacha8-u53-ziggurat/`) |
| --- | --- | --- | --- |
| `static_direct` | yes | closed models: a fixed set of attribute vectors, any number of reactant molecules in `model/1.0` | `rules-static-direct`, with a suffix for the selector and the trackers |
| `next_reaction` | yes | the same models | `rules-next-reaction/1` |
| `tree` (the default) | yes | every model without an ATP block or continuous variables; instances may change identity | `rules-direct/2` |

Static Direct lowers the model to channels once, at preparation: species,
channels, the changes each channel makes and the channels whose rates each
change affects. Each event then recomputes only the affected rates, draws a
waiting time and selects a channel. A zero total rate is absorbing: the state
is held to the horizon without consuming random numbers. Report values are
right-continuous, so an event exactly at a report time is included.

Three rules hold for every selector:

- A positive channel must be at least 2⁻⁵³ of the total rate. Below that, a
  53-bit uniform cannot select it exactly, and the replicate fails with
  `CRN-R031` rather than dropping the channel.
- A non-finite or negative rate fails the replicate (`CRN-R030`). There is no
  fallback to the last channel.
- Under `run/1.0` a positive wait too small to move the clock is accepted;
  1,000,000 such steps in a row fail with `CRN-R030`. Every `rules-run`
  version fails at the first one.

## Selectors

Under `crn-craft.run/1.0`, `method.selector` of `static_direct` chooses how
the next channel is selected. Every selector samples the same law; each uses
the random stream differently, so the same seed gives different paths, and the
sampler identity names the selector.

| `selector` | How it selects | Sampler suffix |
| --- | --- | --- |
| `linear` | a scan of the channels in order, as in every `rules-run` version | none |
| `composition_rejection` | channels grouped by the binary exponent of their rate; a group by its mass, a member by rejection | `-cr` |
| `sorting` | a scan in an order that moves each fired channel one place forward | `-sorting` |
| `rssa_cr` | candidates by composition–rejection from upper bounds of the rates over population intervals, accepted against lower bounds or the exact rate | `-rssa-cr` |
| `auto` (the default in `run/1.0`) | one of the others, chosen at preparation | the chosen one |

`auto` chooses RSSA-CR from 8 channels when no rate reads two or more
species, and from 128 channels when the species that start with molecules
average fewer than 8 molecules. Otherwise it chooses composition–rejection
from 128 channels, and otherwise the linear scan. It never chooses `sorting`.
It decides from the initial counts: a model whose populations later grow or
shrink far from them can end up with a slower selector. The thresholds are in
`describe()["selectors"]["auto"]`.

The `rules-run` versions run the linear scan and refuse the field
(`CRN-R020` at `/schema`). A request that the Python helpers write as
`run/1.0` for a `crn_craft.rules` call names `"selector": "linear"`, so
adding quantiles never changes the paths. The top-level `crn_craft.ensemble`
leaves the choice to `auto`. Generated event loops and rates kernels run the
linear scan only, and refuse any other selector, also when `auto` would pick
one (`CRN-R020` at `/method/selector`).

Selectors are ports of published algorithms: composition–rejection (Slepoy,
Thompson and Plimpton, J. Chem. Phys. 128, 205101 (2008)), sorting direct
(McCollum et al., Comput. Biol. Chem. 30 (2006), 39–49) and RSSA-CR
(Thanh, Priami and Zunino, J. Chem. Phys. 141, 134116 (2014), with
composition–rejection for the candidate draw).

## Waiting times

`waiting_time` is `ziggurat` or `inverse`. Ziggurat draws exponential waits
from tables adapted from the Julia ecosystem (see
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md)); the inverse transform is
`-ln(u) / total`. From `rules-run/0.6` and in `run/1.0`, `static_direct`
defaults to Ziggurat waits and to sparse integration of eligible count
trackers. `{"kind": "static_direct", "waiting_time": "inverse",
"sparse_trackers": false}` reproduces the Static Direct of release 0.3.0.

## The random stream

Every replicate draws from its own ChaCha8 stream. The key is the SHA-256 of
the domain `crn-craft.rules/0.1` and the root seed; the replicate
index is the stream id. A 53-bit uniform is the high 53 bits of a 64-bit word
times 2⁻⁵³. Results therefore depend on the seed and the method, never on the
thread count, and a replicate's path is the same in a larger ensemble.

Statistical equivalence of the selectors is checked against exact
distributions: master-equation solutions and the finite state projection, with
predeclared tests and a recorded family error rate (`validation/merge/`).
