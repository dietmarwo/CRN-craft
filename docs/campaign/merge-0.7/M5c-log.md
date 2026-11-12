# M5c log: tau leaping on converted models

Task M5c of [the plan](../../plan-0.7.0-engine-merge.md): run the gates of
`validation/check_tau.py` against the rules `tau` method on converted
models, move the Poisson sampler that the rules engine calls to
`sampling.rs`, and measure the cost of `tau`'s channels² leap condition.

## Declared before running

In `validation/merge/RESULTS.md` ("M5c (tau leaping on converted models):
declaration", commit `f2c5110`), before any run: the gate (sizes, seeds and
tolerances of `check_tau.py`), the rules analog of the approximate
assessment, the budgets, the absence of a binomial gate and the
leap-condition measurement.

## What changed

- **`sampling.rs`**: `poisson_with` (the shared Poisson kernel: the
  multiplication method up to mean 16, Hörmann's PTRS above), its private
  `ptrs` and the constant `PTRS_MEAN` moved from `tau.rs` unchanged (same
  text, same `#[doc(hidden)] pub`). Three tests on the kernel itself, ported
  from `tau.rs`: the Poisson law on both sides of the switch (200,000 draws
  at five means, DKW band at family α = 0.001, means within five standard
  errors), the multiplication method's draws and words for small means, a
  few uniforms per variate for large means and errors of the uniform source
  passed through. They stay when M8 deletes `tau.rs`.
- **`tau.rs`**: keeps its wrapper `poisson` (mean check and poll) and its
  tests, and re-exports `crate::sampling::poisson_with` (hidden) until it
  is removed. **`rules/leap.rs`**: calls `crate::sampling::poisson_with`
  (same function; the leap code is otherwise untouched).
- **`validation/check_tau_rules.py`** (new; `check_tau.py` keeps testing the
  network command until M8): the `TauControl` network model converted with
  `crn-craft convert`, a `crn-craft.run/1.0` request with
  `{"kind": "tau", "epsilon": ε}` run by `crn-craft rules`, the envelopes of
  `check_tau.py` unchanged, and instead of 0.6.0's
  `assessment: approximate_ensemble`: termination `completed`, the sampler
  `chacha8-u53/rules-approximate-cgp-poisson/2` and the echoed method. The
  rows record events and random words (the rules result has no leap
  counts). CI runs it next to `check_tau.py` in the `reference-export` job,
  normally and under `-O` with a diff.
- **`validation/merge/measure_tau_leap_condition.py`**: the declared cost
  measurement.

## Tests written first

The three `sampling.rs` tests were written before the move and failed to
compile (`poisson_with` not found in `sampling`;
`/home/claude/m07/m5cdg/m5c-first-fail.log`); after the move they pass,
as do the `tau.rs` tests and the `tau` integration tests. The gate script
was written and its declaration committed before its first run; it passed
on the first run.

## Measured

| Gate (declared) | Rows | Failures |
| --- | ---: | ---: |
| `check_tau_rules.py`, ε = 0.1, 0.05, 0.025, n = 10,000 each | 3 | **0** |

CDF errors 0.0043, 0.0066, 0.0063 against envelopes 0.071, 0.046, 0.033;
mean errors 0.007, 0.064, 0.058 against 3.89, 2.07, 1.16 (RESULTS.md has
the rows). The network engine's `check_tau.py` in the same build: CDF
errors 0.0227, 0.0148, 0.0081. The rules `tau` takes far fewer leaps on
this model (1.8–2.1 random words per event; an exact step takes 2): its
species bound with g = 2 refuses a leap when ε · n / 2 < 1, which holds for
ε = 0.05 and 0.025 at these counts, so most steps are exact fallbacks.

**binomial**: `check_tau.py` has no binomial variant, so there is no
binomial gate; the rules `binomial` method keeps the laws checked by
`check_rules.py` and the rules tests.

**Leap-condition cost** (declared; W5 and W6 converted, 16 replicates,
1 thread, medians of 3):

| Model | Channels | `tau` ε = 0.05 | `static_direct` | Ratio | Words per event (`tau`) |
| --- | ---: | ---: | ---: | ---: | ---: |
| W5 | 510 | 3.23 s | 0.229 s | 14.1 | 2.0 |
| W6 | 992 | 4.96 s | 0.386 s | 12.9 | 2.0 |

At the declared ε neither model leaps (W5's counts of 4 make every channel
critical; W6's counts of 32 give ε · n / 2 = 0.8 < 1), so this measures
`tau`'s exact fallback, which recomputes every rate and every channel's
criticality at each step, not the channels² condition. Supplementary, not
declared, one run: W6 with ε = 0.1 (leaps attempted), 1 replicate: 9.97 s
against about 0.024 s per replicate for Static Direct (about 400×); each
attempted leap evaluates the 992 channels after one firing of each enabled
channel, twice, about 2 · 992² ≈ 2.0 million channel rates. Not optimized,
as the task says: the gate does not need it. Large networks should use the
exact methods, as the M4a log already advises.

## Checks

On the branch's final tree (M5d, M5g and M5c), release builds, under the
CPU lock, Python 3.12 in `venv-m5cdg` with the wheel built from the root.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 130 (3 new), `archives_v1` 17 (new), `archives` 11, `tau` 5, CLI 16 (1 new), I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26, `generated_kernel` 3), all others, doctests 37 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass after two links to `write()` were made unambiguous (only the known output-filename collision warning) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| Schemas against `crn-craft describe KIND` | 26 byte-identical (2 new: `archive-1.0`, `archive-check-1.0`); `analysis` and `analysis-v1` are no longer generated since M5b, as before |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | 22 blocks pass; `-O` identical |
| `check_contracts.py` | 103 of 103 (2 new schema rows); `-O` identical |
| `check_rules_python.py` | pass (new archive check); `-O` identical |
| `check_tau_rules.py`, `check_tau.py` | pass, 3 rows each; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 17 checks pass; `-O` identical |
| `validate_parameters.py --evidence … --skip msrv-tests` | 38 checks pass, `msrv-tests` skipped (the 1.85.0 toolchain download is refused by the host's network policy), so its summary is not `complete`: fmt (workspace, binding, worker), ruff, release build, worker build, clippy (workspace, binding, worker with all features), release tests, wheel, `check_contracts`, `check_alternative_contracts`, `check_python`, `check_ssa` (17 rows), `check_cme` (14), default and fast exact laws (210 rows each), each `-O` identical, `simulate-cached` reports equal, dependencies, installed, both parameter scripts |
| `check_jobs.py`, `check_parameter_binding.py`, `check_parameter_interface.py` | see the [M5g log](M5g-log.md#checks) |
| `integrations/jobs`, parameter worker builds with `--locked` | pass |

Not run: the full `mutation_check.py` (it rebuilds a copy of the tree; its
archive and tau anchors are unchanged and still unique, checked with the
script's own anchor test), `check_compiled_loop.py` (the generators are
unchanged; I3 holds by the fixtures), `check_mtdna.py`, `check_hybrid.py`,
`check_builder.py`, `check_tutorial.py` (paths these tasks do not touch).

## Deviations from the plan

- The migrated script is the new `check_tau_rules.py`, and `check_tau.py`
  stays on the old command until M8 (both in CI).
- The rules analog of the approximate assessment is the sampler identity
  and the echoed method (the rules result has no `assessment`).
- `tau.rs` re-exports the moved function (hidden) so that nothing that
  compiled against `crn_craft::tau::poisson_with` breaks before M8.

## Limitations

- On this gate model the rules `tau` mostly takes exact steps, so the gate
  shows that its laws stay within the envelopes, not how much it gains
  from leaping.
- `tau`'s exact fallback is O(channels) per event and its leap condition
  O(channels²) per attempted leap; neither was optimized.

## Migration notes

- `crn_craft::sampling::poisson_with` replaces
  `crn_craft::tau::poisson_with` (still re-exported, hidden, until M8).
- The rules `tau` method of `crn-craft.run/1.0` is the tau leaping of the
  unified engine; `crn-craft tau` and `tau-request/0.1` stay until M8.
