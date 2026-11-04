# crn-craft-cli

The native command-line interface of
[CRN-craft](https://github.com/dietmarwo/CRN-craft), installed as `crn-craft`:

```sh
cargo install crn-craft-cli --locked
crn-craft describe
```

It checks, simulates, analyzes, exports and verifies reaction-network JSON
documents, runs the exact SSA methods, finite state projection and approximate
tau leaping, and runs the experimental rules engine for dynamic, attributed
species (`rules-check`, `rules`). Every command reads and writes versioned
JSON; `crn-craft describe` lists commands, methods and limits.

The libraries are [`crn-craft`](https://crates.io/crates/crn-craft) and
[`crn-craft-rules`](https://crates.io/crates/crn-craft-rules); the Python
package is [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/). See the
[repository](https://github.com/dietmarwo/CRN-craft) for the tutorial and
contracts. A successful numerical run is not biological validation. MIT
licensed.
