# Direct SSA in M2

`0.1.0-dev.4` implements a time-homogeneous, well-mixed stochastic count
process. Mass action, Hill and Michaelis–Menten use the existing validated
propensity evaluator and reactant-availability gates. FSP executes separately
via `fsp`. This is an exact **method**, with finite-precision arithmetic and
Monte Carlo uncertainty; it is not a biological validation claim.

```bash
cargo build --release --locked
./target/release/crn-craft simulate \
  models/examples/immigration-death.model.json \
  models/examples/immigration-death.run.json /tmp/crn-ssa-example
./target/release/crn-craft verify /tmp/crn-ssa-example
```

The output directory must be new and its parent must exist. `simulate` exits
0 only for a completed ensemble. An event/population/deadline/cancellation/
output stop produces exit 1 and, when finalization succeeds, an inspectable
archive. `verify` exits 0 for an internally consistent archive even if its
simulation is incomplete. Always inspect both `archive_complete` and
`simulation_complete`.

## Event and floating-point policy

The Direct method recomputes all rates in reaction declaration order. The
scheduler shares a known-valid count state across channel evaluations, so
state validation is not repeated for every channel. Checked reaction updates
remain transactional. Null reactions still consume waits, draws and events.

Sum the propensities in declaration order, draw an exponential wait and select
the first positive channel whose cumulative propensity exceeds the selection
threshold. There is no last-channel fallback. A zero total is absorbing: hold
the state to the horizon without consuming random numbers.

The sampler is `crn-craft.chacha8-direct/0.1`:

- `rand_chacha = 0.3.1`, `rand_core = 0.6.4`, default features disabled;
  full transitive pins are in Cargo.lock. No operating-system randomness.
- Key = SHA-256 of ASCII `crn-craft.chacha8/0.1`, one NUL byte, then the
  root seed as eight little-endian bytes. Set the ChaCha 64-bit stream ID to
  the replicate index and start at word position zero.
- Draw `next_u64`, take its high 53 bits and multiply by 2^-53. Selection uses
  [0,1); exponential sampling rejects zero, giving (0,1).
- Wait = `-libm::log(U) / total`, using pinned pure-Rust libm 0.2.15.
  Selection consumes a separate word only for an event within the horizon.

CRN-E042 still rejects a positive propensity that underflows. CRN-E061 rejects
nonfinite sums, positive additions lost in rounding and positive channel
probabilities below the 53-bit selection resolution. A threshold rounded to
the total is also an error. CRN-E062 rejects nonfinite/nonpositive waits or a
clock that cannot advance. These restrictions deliberately fail closed for
extreme dynamic ranges rather than claim arbitrary-range exact sampling.

Saved states include events at or before their time. One pending event is
retained while crossing save times; changing output density does not resample
it. An event exactly at the horizon is included. An event cap is a conservative
stop immediately after the allowed number of events unless the horizon was
already reached; the solver does not do additional work to discover later
absorption. `stop_time` is the latest certified event/sample time, or the
horizon for a completed path. A safety-limit violation is detected before
mutating the state or incrementing the event counter.

## Ensembles and statistics

Workers receive round-robin replicate IDs and return results through one
bounded channel per worker. The parent consumes IDs in ascending order. At
most a small fixed number of paths per worker are in flight; there is no
unbounded out-of-order buffer. Thread count is an allowance: use at most 64
workers and no more workers than replicates. Changing workers, saved output
or model metadata does not change the seed/stream mapping.

All requested replicates have records, including those never started after a
stop. Incomplete ensembles have **no ensemble statistics**. Reporting the
mean of only completed paths could introduce state-dependent selection bias.
The retained paths are the first requested IDs, including their available
saved-state prefixes; they are not a chosen successful subset.

Fully completed ensembles use ordered Welford updates. Requested `variance`
is the unbiased sample variance (denominator n-1). Requested `mean` includes
SE sqrt(sample variance/n) for n>1; a single replicate has no estimated SE.
No variance-SE or confidence-interval estimator is implied. Counts remain
exact i64; summary arithmetic is f64, including conversion of counts above
2^53. RNG and libm pins do not themselves establish cross-platform bit identity.

## Execution limits and cancellation

The unchanged run 0.1 budget supplies event cap, optional wall time, worker
allowance and maximum output bytes. The engine additionally caps requested
replicates at 1,000,000, conservatively reserved workspace at 256 MiB and
encoded archive contents at 256 MiB. Feasible combinations can be much smaller
than each individual maximum. Workspace admission accounts for bounded worker
paths/stacks, retained paths, reduction arrays and per-replicate evidence; it
is an allocation reservation, not a hard OS RSS limit. Model compilation and
JSON verification have their own size limits and allocation overhead.

The engine reserves room for input bytes and accounting before starting. It
counts serialized result bytes without allocating the serialization. If
retained numerical output does not fit, it removes trajectories/statistics,
records CRN-E064 and `output_limit`, and keeps every replicate's work and
termination record. Archive writes enforce the total file-content byte cap
again. Filesystem allocation units are not counted as content bytes.

`simulation::Control::cancel()` supports library callers. The CLI accepts
`--cancel-file PATH`; creating that file requests cooperative cancellation.
An already present file cancels before work starts. Polls occur at every
event, at most every 64 rate evaluations, and for every saved point. Each
poll reads the cancellation flag. Every 64th poll of a worker, starting with
its first, also reads the wall clock and, at most once per millisecond, checks
the cancel file; a found file sets the cancellation flag for all workers.
A worker keeps a stop once it has seen one, so its remaining replicates stay
unstarted. Without a file or wall-time limit, polls read only the flag.
Checking the file at every poll made small models about ten times slower (a
stat per poll against roughly 100 ns per event), and reading the clock at every
poll about doubled the cost of a wall-time limit. The wall clock starts before
worker setup and is polled during trajectory generation. Pending ordered collection, statistic reduction and finalization
may continue afterward; input parsing and model compilation precede it.
The process must still be allowed to finalize bounded evidence. This is not
a hard real-time deadline, signal handler or resumable checkpoint protocol.
Forced termination leaves an unmarked directory, which verification rejects.

## Rust API

`simulation::prepare(model_bytes, run_bytes)` validates exact input bytes and
execution admission. `execute(&prepared, &control)` returns numerical results
and separate provenance. `simulate_bytes` combines those two operations. The
byte API checks the model hash; only `archive::simulate_files` can validate
the run's relative filesystem reference. `archive::verify_archive` checks a
finished archive without running a solver. See [archive contracts](archive-contract.md).

## References and verification

Gillespie, 1977, [Exact stochastic simulation of coupled chemical reactions](https://doi.org/10.1021/j100540a008).
Bernstein, 2008, [ChaCha, a variant of Salsa20](https://cr.yp.to/chacha/chacha-20080128.pdf).
The installed rand_chacha 0.3.1 source documents its 64-bit counter/stream
layout; the independent scalar Python implementation tests 360 words across
nine pinned streams and multiple blocks. The acceptance plan and statistical
family are in [iteration-04-validation.md](decisions/iteration-04-validation.md).
