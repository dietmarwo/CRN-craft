# Alternative exact methods

This page described the network engine's alternative exact methods
(DirectCR, SortingDirect, RSSA-CR and Gibson–Bruck next reaction, with
inverse or Ziggurat waits) and their `ssa-run` and `ssa-result` contracts.
Those contracts and their random stream (xoshiro256++) were removed in
release 0.7.0. The methods themselves live on in the one engine.

| 0.6.0 | 0.7.0 |
| --- | --- |
| `aggregator: "direct"`, `"cached_direct"` | `method: "static_direct"` with `selector: "linear"` |
| `aggregator: "direct_cr"` | `selector: "composition_rejection"` |
| `aggregator: "sorting_direct"` | `selector: "sorting"` |
| `aggregator: "rssa_cr"` | `selector: "rssa_cr"` |
| `aggregator: "next_reaction"` | `method: "next_reaction"` |
| `waiting_time: "inverse"` or `"ziggurat"` | the same field of `static_direct` |
| `fast: true` or `false` | no such switch |
| `crn.ssa_spec`, `simulate_ssa`, `prepare_ssa` | `crn_craft.ensemble(..., method="static_direct", selector=...)` |
| `crn-craft ssa MODEL REQUEST` | `crn-craft run MODEL REQUEST` with a `run/1.0` request |

Results are not byte-identical to those of 0.6.0, because the random stream
changed; they sample the same laws. See [exact methods](direct-ssa.md) for
the selectors, their sampler identities and the `auto` rule, and the
[changelog](../CHANGELOG.md) for the full migration tables.
