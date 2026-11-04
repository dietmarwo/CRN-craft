# crn-craft-rules

Experimental exact stochastic simulation of dynamic, attributed species for
[CRN-craft](https://github.com/dietmarwo/CRN-craft): species types with numeric
attributes, instances created and removed during a run, rule schemata over
zero to two reactant patterns, observables and stop conditions (contracts
`crn-craft.rules-model/0.1` and `crn-craft.rules-run/0.1`). Exact-law gates
cover Kowald's fixed-population mtDNA model and aggregation kernels.

```toml
[dependencies]
crn-craft-rules = "0.2"
```

The native CLI ([`crn-craft-cli`](https://crates.io/crates/crn-craft-cli))
runs rules models with `crn-craft rules-check` and `crn-craft rules`. The
language reference is
[docs/rules-language.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md).
MIT licensed.
