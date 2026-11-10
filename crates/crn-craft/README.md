# crn-craft

The core library of [CRN-craft](https://github.com/dietmarwo/CRN-craft):
reproducible simulation and analysis of well-mixed chemical reaction networks.
Models, requests and results are versioned JSON documents with strict input
checks, explicit budgets and replayable seeds.

- Exact stochastic simulation: checked Direct and dependency-cached Direct
  archives, and `direct`, `cached_direct`, `direct_cr`, `sorting_direct` and
  `rssa_cr` with inverse or Ziggurat waits (`crn-craft.ssa-run/0.3`, where the
  JumpProcesses.jl-style `fast` policy is the default and `fast: false` the
  checked one), plus Gibson–Bruck `next_reaction` (`crn-craft.ssa-run/0.4`,
  checked policy).
- Sparse finite state projection, explicitly approximate Poisson leaping,
  network analysis and conservation laws.
- Reusable parameter templates for optimization loops.

```toml
[dependencies]
crn-craft = "0.6"
```

The command-line interface is the
[`crn-craft-cli`](https://crates.io/crates/crn-craft-cli) crate
(`cargo install crn-craft-cli --locked`), dynamic attributed species are
[`crn-craft-rules`](https://crates.io/crates/crn-craft-rules), and the Python
package is [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/).
Documentation, contracts, example models and the fixture-based test suite are
in the [repository](https://github.com/dietmarwo/CRN-craft). A successful
numerical run is not biological validation. MIT licensed; see
THIRD_PARTY_NOTICES.md for the Ziggurat tables adapted from Julia.
