# crn-craft-rules

Experimental exact stochastic simulation of dynamic, attributed species for
[CRN-craft](https://github.com/dietmarwo/CRN-craft). Species types carry
numeric attributes; instances appear and disappear during a run. Rule
schemata over zero to two reactant patterns, observables and stop conditions
describe models whose species cannot be listed in advance, such as
mitochondrial DNA clones, polymers or clusters. Exact-law gates cover a
fixed-population mtDNA model, immigration–death, immigration–dimerization
and aggregation kernels.

```toml
[dependencies]
crn-craft-rules = "0.2"
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

- **Contracts:** `crn-craft.rules-model/0.1` and `0.2` (ATP coupling),
  `crn-craft.rules-check/0.1`, `crn-craft.rules-run/0.1` to `0.4` and
  `crn-craft.rules-result/0.1` and `0.2`. JSON Schemas come from
  `crn_craft_rules::schema`.
- **Methods:** exact `tree` (dynamic instances), `static_direct`,
  `mtdna_direct`, `next_reaction` and `thinning` (ATP); approximate `tau` and
  `binomial`.
- **Observation:** trackers and bounded event logs inside the result, and an
  `EventObserver` interface for live callbacks on the worker threads.
- **Generated code (Rust only):** `generate_rules_kernel` and
  `generate_rules_loop` emit Rust source for a Static Direct model;
  `execute_with_rules_kernel` and `execute_with_rules_loop` run it with
  results byte-identical to the interpreter.

Completed results depend only on the input bytes, the seed and the sampler
identity, not on the thread count.

The native CLI ([`crn-craft-cli`](https://crates.io/crates/crn-craft-cli))
runs rules models with `crn-craft rules-check` and `crn-craft rules`, and the
Python package [`crn-craft`](https://pypi.org/project/crn-craft/) exposes
them as `crn_craft.rules`. Further reading:

- [Rules language](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md)
- [mtDNA ensembles](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-rules.md)
  and the [mtDNA user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md)
- [Generated rules kernels and event loops](https://github.com/dietmarwo/CRN-craft/blob/main/docs/specialized-kernels.md)

MIT licensed.
