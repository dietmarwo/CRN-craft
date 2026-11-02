# Alternative exact SSA methods

This extension implements DirectCR, SortingDirect and RSSACR, plus separately
versioned inverse-transform and exponential Ziggurat waits. The original
`simulate` / `simulate-cached` commands and their model/run/result/archive
formats and random streams remain unchanged.

```python
import crn_craft as crn
model = crn.load_example("immigration-death")
request = crn.ssa_spec(
    model, aggregator="rssa_cr", waiting_time="ziggurat",
    end=100.0, save_times=[0.0, 100.0], replicates=64, seed=4242,
    max_events=1_000_000, max_sampling_steps=10_000_000,
)
prepared = crn.prepare_ssa(model, request)
result = prepared.simulate()
if result["termination"] != "completed":
    raise RuntimeError(result["diagnostics"])
print(result["statistics"])
```

`crn.simulate_ssa(model, request)` is the corresponding one-shot call. A prepared
handle is immutable: repeated calls restart the fixed seeds, and concurrent calls
own independent execution state. Prepare a new request to change inputs. The
same reuse interface exists for legacy execution: `crn.prepare(model, run,
cached=True).simulate()`.

The native and installed Python CLIs both provide:

```sh
crn-craft describe ssa-run
crn-craft describe ssa-result
crn-craft ssa model.json request.json [--cancel-file PATH]
```

The request schema is `crn-craft.ssa-run/0.1` or `0.2`, and the result is
`crn-craft.ssa-result/0.1` or `0.2` to match. These are standalone JSON results, not legacy sealed
archives. The command exits 0 only for a completed ensemble; incomplete runs
exit 1 and withhold ensemble statistics. Unknown fields, duplicate keys and
mismatched model references are rejected. `max_output_bytes` covers compact
result JSON including its trailing newline. If even minimal diagnostic evidence
cannot fit, execution returns a coded error instead of oversized evidence.

### The `fast` switch and ssa-run/0.2

`crn-craft.ssa-run/0.2` adds an optional `method.fast` boolean (default
`false`). With `fast: true`, every method keeps its selection algorithm but
emulates JumpProcesses.jl's numerical policy:

- random words come from Xoshiro256++, one state per replicate, seeded from a
  keyed ChaCha8 stream with its own domain (`crn-craft.xoshiro256pp-fast/0.1`),
  so results stay independent of the thread count and share no words with the
  default methods at the same root seed;
- DirectCR, SortingDirect and RSSACR keep a running total updated by
  `+= new - old` instead of the pairwise tree, recomputed from the leaves every
  65,536 updates or when it falls below 2^-30 of its recent peak; Direct and
  cached Direct keep their prefix sums by additions only;
- no lost-increment or 53-bit resolution checks; a positive rate below half an
  ulp of its running sum is then never selected, as in Julia;
- cooperative stops are polled every 256 event-loop iterations and before
  every rejection retry, instead of before every random word.

Selection remains exact for the totals in use (composition–rejection still
accepts with probability rate/bound, and SortingDirect redraws a threshold
beyond its scanned sum). Sampling budgets, population and event limits,
checked state updates, cancellation and output budgets still apply. The
sampler identity is `xoshiro256pp-u53/<method>/<wait>/fast-1`. The fast policy
has its own exact-law gate: `check_alternative_ssa.py --fast` runs the same 210
comparisons (100,000 paths per fixture, family alpha 0.001); it passes with a
largest CDF error of 0.00526 against the 0.00805 limit and rejects all 20 wrong
laws (`validation/fast/fast-laws.json`).

`fast` is off by default because the default's guarantees cannot be checked
after the fact: a result does not reveal that a tiny channel lost its
probability or that a total drifted. Use it for well-scaled models where speed
matters, and for benchmarks against tools that make the same choices.

**Migrating from ssa-run/0.1.** Nothing changes for existing requests: 0.1 is
still accepted, rejects `fast`, and returns `crn-craft.ssa-result/0.1` with the
same bytes as before. To use `fast`, set `schema` to `crn-craft.ssa-run/0.2`
and add `"fast": true` to `method`; the result is then
`crn-craft.ssa-result/0.2`, and its `method` echoes `fast: true` (omitted when
false). `crn.ssa_spec(..., fast=True)` writes such a request; without `fast`
it still writes 0.1. Discovery (`describe/0.7`) lists both request versions
under `alternative_exact_ssa`.

## Methods and numerical policy

| `aggregator` | Selection and updates |
| --- | --- |
| `direct` | Full ordered propensity sum and binary prefix selection |
| `cached_direct` | Dependency updates and resumed ordered prefix sums |
| `direct_cr` | Dependency updates, balanced totals, power-of-two composition/rejection |
| `sorting_direct` | Dependency updates, balanced totals, linear search; fired channel moves one place forward |
| `rssa_cr` | Bounded propensity envelopes, composition/rejection, and thinning; refresh only after leaving a population bracket |

`waiting_time` is `inverse` or `ziggurat`. Each pair has its own sampler identifier.
The new Direct/cached inverse variants retain legacy paths; alternative selection
or wait sampling changes random consumption. The Ziggurat table/algorithm is
adapted from Julia 1.13.3 under MIT terms in [THIRD_PARTY_NOTICES](../THIRD_PARTY_NOTICES.md).
No Julia runtime is required. Replays within each pair are independent of worker
count. SortingDirect resets its learned order per replicate.

DirectCR is an original variant in the composition–rejection family: a bucket
is chosen in proportion to its envelope mass, a member is chosen using unbiased
integer sampling, and rejection yields actual-rate proportions. For buckets of
up to 1024 members, one random word supplies both the member (Lemire's
multiply-and-reject) and the acceptance uniform (the top 53 bits of the low
product word, which given the member is uniform on a lattice of spacing at
most 2^-54); larger buckets draw them separately. This is sampler version 2 of
`direct_cr` and `rssa_cr`. The event wait
uses the actual total. It is not a port of Julia's internal data structures and
is not strictly constant-time: dependency updates cost O(log R), and selection
uses a binary search over cached bucket envelope sums, rebuilt only when a
bucket gains or loses a member. Balanced trees rebuild sums from children,
avoiding incremental add/subtract drift. Non-finite totals, lost positive additions and channels below
53-bit selection resolution fail explicitly. Tiny or overflowing waits also
fail. Integer state updates remain checked and transactional.

Why selection stays exact to 53 bits: bucket bounds are powers of two, so the
acceptance test `U·2^k < a` has probability exactly `a/2^k` at 53-bit
resolution, and the member index uses unbiased integer sampling (with the joint word, the
acceptance probability errs by at most 2^-53 + 2^-54). Bucket choice
and totals round at about 2^-53 relative, as in legacy Direct, and the wait
always uses the actual total. The guarantee that every positive channel is
resolvable comes from the tree's check that the smallest positive leaf is at
least 2^-53 of the total. The tree evaluates that check lazily: it keeps a floor
that never exceeds a current positive leaf, returns the total directly when the
floor is at least 2^-52 of a normal total (then no pairwise addition can lose a
positive operand, because that would need at most half an ulp, i.e. 2^-53 of a
sum no larger than the total), and otherwise evaluates the fully checked tree.
Decisions and totals are identical to checking every node. Two edge cases are
known: bucket bounds below 2^-1021 are subnormal, so the acceptance product
rounds (absolute error at most about 2^-53); and the top bucket (rates at or
above 2^1023) uses `f64::MAX`, which is not a power of two. Neither occurs in a
model whose total propensity is finite and whose rates are normal numbers.

SortingDirect scales its threshold by the pairwise tree total but scans in
learned order; the two sums can differ in the last bit. A threshold at or above
the sequential sum is redrawn (counted against the sampling budget), so accepted
thresholds are uniform below the sum the scan partitions.

One-variable Hill laws
memoize their checked value per integer regulator count below 8192 (first 128
Hill laws of a network, lazily allocated, at most 64 KiB each); failures are
never stored. Michaelis–Menten and mass action are evaluated directly.

RSSACR uses thinning of a process with the summed upper-bound rate. Every
candidate advances time, including a rejected candidate; saved states retain
the same right-continuous convention as legacy Direct. It accepts below a
certified lower bound without evaluating the propensity. Otherwise it evaluates
the current gated propensity and rejects or fires accordingly.

For a population x>0, a bracket has half-width `clamp(x/10,4,64)` (integer division),
clamped to nonnegative i64 bounds. Zero has bracket [0,0]. Bounds cover the actual
floating-point evaluator on the integer box:

- Normal elementary mass-action products with a safe coefficient range use
  monotone endpoints, including the minimum positive enabled-state rate.
- Hill and Michaelis–Menten rates enumerate the feasible integer regulator range
  (at most 129 values). Availability gates are separate. A consumed repressor can
  be disabled at zero and maximal at an interior state; endpoint-only bounds
  would be unsound in this case.
- Higher-order/subnormal/overflow-risk mass action pins every read coordinate to
  its current value. A numerical error at a neighboring bound state also triggers
  an exact point fallback when the actual current propensity is valid. Pinned
  coordinates refresh on every change. Current-state errors are never hidden.

When no lower bound guarantees a positive reaction, RSSACR checks current rates
after each accepted state change (stopping at the first positive rate). This
detects true absorption even if neighboring states keep the envelopes positive.

Population, event, output, workspace, cancellation and wall-time controls still
apply. `max_sampling_steps_per_replicate` counts **every u64 random word**, including
all zero-uniform retries, bucket/index/acceptance retries and Ziggurat tails.
It is not a number of accepted events. Exhaustion returns `sampling_limit` with
CRN-E067, so a high rejection rate cannot bypass the work budget. Bounds also poll
cooperative stops. Work fields report propensity evaluations, selection proposals,
bound updates, point fallbacks, rejected RSSA candidates and random words.

Bucket buffers shrink as memberships fall, so memory is bounded by live reaction
membership rather than every exponent bucket visited over a long trajectory.
Request validation
reserves selector/bound/reverse-reader workspace and rejects counter bounds that
could overflow u64. Admission still uses the legacy archive reserve (model and
request size plus a fixed envelope) when it checks `max_output_bytes` up front,
so a budget smaller than that reserve is rejected even if the compact result
would fit; the result itself is then checked once against the full budget.

New methods draw from the same per-path ChaCha8 streams as legacy Direct for
the same `root_seed`. Different methods consume the words differently, but
their early states are strongly correlated with legacy cached Direct at equal
seeds (about 0.6–0.8 at short times on a birth–death model). Use different
seeds when comparing methods statistically.

Known limitation: RSSACR passes the smallest positive rate anywhere in a box to
the resolution check. A very steep Hill law (for example repression with n=130
and the regulator at K) can then fail with CRN-E061 at t=0 even though the
current rates are resolvable, while `direct`, `cached_direct` and `direct_cr`
complete. Use one of those methods for such models.

## Validation and scope

The [predeclared protocol](decisions/ssa-acceleration-2026-11-02.md) fixes the
100,000-path, 210-marginal statistical family and the million-sample sampler gates.
Independent integer-box enumeration covers mixed availability/regulator directions,
consumed repression and exact point fallbacks. The implementation log retains
failed development checks and the final results. These methods support the same
four built-in rate laws as the engine, not arbitrary time-dependent expressions,
delays, jump–ODE coupling or Julia's full modeling language.

Algorithm references: Slepoy, Thompson & Plimpton (2008),
[doi:10.1063/1.2919546](https://doi.org/10.1063/1.2919546); McCollum et al. (2006),
[doi:10.1016/j.compbiolchem.2005.10.007](https://doi.org/10.1016/j.compbiolchem.2005.10.007);
Thanh, Priami & Zunino (2014), [doi:10.1063/1.4896985](https://doi.org/10.1063/1.4896985),
and Thanh, Zunino & Priami (2017), [doi:10.1109/TCBB.2016.2530066](https://doi.org/10.1109/TCBB.2016.2530066).
The Ziggurat algorithm is Marsaglia & Tsang (2000),
[doi:10.18637/jss.v005.i08](https://doi.org/10.18637/jss.v005.i08).
