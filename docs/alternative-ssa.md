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

The request schema is `crn-craft.ssa-run/0.1`, and the result is
`crn-craft.ssa-result/0.1`. These are standalone JSON results, not legacy sealed
archives. The command exits 0 only for a completed ensemble; incomplete runs
exit 1 and withhold ensemble statistics. Unknown fields, duplicate keys and
mismatched model references are rejected. `max_output_bytes` covers compact
result JSON including its trailing newline. If even minimal diagnostic evidence
cannot fit, execution returns a coded error instead of oversized evidence.

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
integer sampling, and rejection yields actual-rate proportions. The event wait
uses the actual total. It is not a port of Julia's internal data structures and
is not strictly constant-time: dependency updates cost O(log R), and selection
scans active exponent buckets. Balanced trees rebuild sums from children, avoiding
incremental add/subtract drift. Non-finite totals, lost positive additions and
channels below 53-bit selection resolution fail explicitly. Tiny or overflowing
waits also fail. Integer state updates remain checked and transactional.

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
Request validation reserves selector/bound/reverse-reader workspace and rejects
counter bounds that could overflow u64.

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
