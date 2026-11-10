# crn-craft-rules

Experimental stochastic simulation of dynamic, attributed species for
[CRN-craft](https://github.com/dietmarwo/CRN-craft). Species types carry
numeric attributes; instances appear and disappear during a run. Rule
schemata over zero to two reactant patterns, observables and stop conditions
describe models whose species cannot be listed in advance, such as
mitochondrial DNA clones, polymers or clusters. A rule gives a rate per
reactant tuple or, since rules-model/0.3, its whole `propensity`; continuous
variables (ODEs between events that rates read) are simulated by the
approximate `hybrid` method, which integrates rules-model/0.4 independent
variables once per run. Exact-law gates cover a fixed-population mtDNA
model, immigration–death, immigration–dimerization, aggregation kernels and,
for `hybrid`, continuous variables with closed-form laws.

```toml
[dependencies]
crn-craft-rules = "0.6"
```

```rust
use crn_craft_rules::{Control, execute, prepare};

// Validate the exact model and request bytes once; errors are coded diagnostics.
let prepared = prepare(&model_bytes, &request_bytes).expect("valid model and request");
// Run the ensemble; replicates that hit a limit are reported in the result.
let result = execute(&prepared, &Control::default()).expect("result fits the output budget");
println!("{:?}: {} events", result.termination, result.events);
```

The [API documentation](https://docs.rs/crn-craft-rules) starts with a
complete, runnable example.

- **Contracts:** `crn-craft.rules-model/0.1`, `0.2` (the legacy ATP block),
  `0.3` (whole-rule `propensity` and `continuous` variables) and `0.4`
  (`independent` continuous variables),
  `crn-craft.rules-check/0.1`, `crn-craft.rules-run/0.1` to `0.6` (0.5 adds
  `hybrid` with `rtol` and `atol`; 0.6 makes Ziggurat waits and sparse
  count trackers the `static_direct` defaults) and
  `crn-craft.rules-result/0.1` to `0.3` (0.3 records the method and the
  continuous columns). Stored 0.1–0.3 models and 0.1–0.6 requests replay
  byte-identically. JSON Schemas come from
  `crn_craft_rules::schema`.
- **Methods:** exact `tree` (dynamic instances), `static_direct`,
  `mtdna_direct` (an opt-in fast path for one model structure, with a
  Static Direct fallback), `next_reaction` and `thinning` (legacy ATP block);
  approximate `tau`, `binomial` and `hybrid` (continuous variables,
  Dormand–Prince 5(4) at the request's tolerances).
- **Observation:** trackers and bounded event logs inside the result, and an
  `EventObserver` interface for live callbacks on the worker threads.
- **Generated code:** `generate_rules_kernel` and
  `generate_rules_loop` emit Rust source for a Static Direct model (the
  loop for every rate expression, string rate and observable Static Direct
  runs);
  `execute_with_rules_kernel` and `execute_with_rules_loop` run it with
  results byte-identical to the interpreter. From Python,
  `crn_craft.rules.compile_loop` builds a generated loop against this crate
  and runs it as a separate program.

Completed results depend only on the input bytes, the seed and the sampler
identity, not on the thread count.

The native CLI ([`crn-craft-cli`](https://crates.io/crates/crn-craft-cli))
runs rules models with `crn-craft rules-check` and `crn-craft rules`, and the
Python package [`crn-craft`](https://pypi.org/project/crn-craft/) exposes
them as `crn_craft.rules`, with a `rules.Model` builder and pandas results.
Further reading:

- [Rules language](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md)
- [mtDNA ensembles](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-rules.md)
  and the [mtDNA user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md)
- [Generated rules kernels and event loops](https://github.com/dietmarwo/CRN-craft/blob/main/docs/specialized-kernels.md)

MIT licensed.
