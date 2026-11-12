# crn-craft

The engine of [CRN-craft](https://github.com/dietmarwo/CRN-craft):
reproducible stochastic simulation and analysis of well-mixed chemical
reaction networks, from fixed species lists to species with attributes that
appear and disappear during a run. Models, requests and results are
versioned JSON documents with strict input checks, explicit budgets and
replayable seeds.

- One model format, `crn-craft.model/1.0`: species types with optional
  numeric attributes, rules with rate or whole-rule propensity expressions,
  observables, continuous variables and higher-order patterns. Network
  models (`crn-craft.model/0.1`) convert with
  `crn_craft::rules::convert::from_network`, and the `crn-craft.rules-model`
  0.1–0.4 and `crn-craft.rules-run` 0.1–0.6 documents of earlier releases
  still run and replay byte for byte.
- Exact simulation in `crn_craft::rules`: Static Direct with the event
  selectors `linear`, `composition_rejection`, `sorting` and `rssa_cr`, or
  `auto` (the default of `crn-craft.run/1.0`), Gibson–Bruck
  `next_reaction`, the `tree` method for species that appear and
  disappear, and a certified fast path for one mtDNA model structure.
  Explicitly approximate `tau`, `binomial` and `hybrid`. One ChaCha8 stream
  per replicate, so results do not depend on the thread count; means,
  variances and quantiles per report time.
- On the flat view of a closed model: sparse finite state projection,
  conservation laws and structure analysis, ODE right-hand sides and
  Jacobians, and GillesPy2 export. Sealed `crn-craft.archive/1.0` result
  directories with a verifier.
- Opt-in generated Static Direct event loops (Rust source), byte-identical
  to the interpreter.

```toml
[dependencies]
crn-craft = "0.6"
```

The command-line interface is the
[`crn-craft-cli`](https://crates.io/crates/crn-craft-cli) crate
(`cargo install crn-craft-cli --locked`), and the Python package is
[`crn-craft` on PyPI](https://pypi.org/project/crn-craft/). The crate
[`crn-craft-rules`](https://crates.io/crates/crn-craft-rules) is the
deprecated 0.6.0 name of `crn_craft::rules`, kept as a re-export.
Documentation, contracts, example models and the fixture-based test suite are
in the [repository](https://github.com/dietmarwo/CRN-craft). A successful
numerical run is not biological validation. MIT licensed; see
THIRD_PARTY_NOTICES.md for the Ziggurat tables adapted from Julia.
