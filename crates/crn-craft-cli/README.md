# crn-craft-cli

The native command-line interface of
[CRN-craft](https://github.com/dietmarwo/CRN-craft), installed as `crn-craft`:

```sh
cargo install crn-craft-cli --locked
crn-craft describe
```

It checks, runs, analyzes, converts, exports and verifies CRN-craft JSON
documents. `crn-craft run MODEL REQUEST` runs a `crn-craft.run/1.0` request
on a `crn-craft.model/1.0` model with any exact or approximate method, or
computes a finite state projection; `--archive NEW_DIRECTORY` seals the run
and `crn-craft verify DIRECTORY` checks it. Network models of release 0.6.0
(`crn-craft.model/0.1`) convert with `crn-craft convert MODEL`. Every command
reads and writes versioned JSON; `crn-craft describe` lists the commands,
contracts, methods, selectors and limits, and `crn-craft describe KIND`
prints a JSON Schema. `rules` and `rules-check` remain aliases of `run` and
`check`.

The library is [`crn-craft`](https://crates.io/crates/crn-craft); the Python
package is [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/). See the
[repository](https://github.com/dietmarwo/CRN-craft) for the tutorial and
contracts. A successful numerical run is not biological validation. MIT
licensed.
