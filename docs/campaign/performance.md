# Ten performance experiments

The functional campaign is complete before P01. P00 is its release binary;
GitHub 5a07334 remains the original comparison baseline. Each iteration records
source changes, a passing correctness gate, seven alternating paired timings,
exact numerical digests, acceptance/rejection and a local commit. Tests and
benchmarks are serial; heavy optional builds are finished before timing.

Frozen workloads and acceptance are in acceptance.md and benchmark_campaign.py.
Raw reports live in validation/campaign/P*.json. Rejected candidates retain their
diff and evidence; no rejected change is silently counted as an improvement.
Final tables and limitations are generated after P10. Shared-host ranges are
reported; these are workload measurements, not a universal solver ranking.

P01 accepted: normal-range zero/first-order mass-action fast paths, scaled fallback
for exceptional values. All exact digests match; paired median speedups 1.152x two
channels, 1.043x sparse 66, 1.066x circuit, 1.010x dense. Correctness/fmt/clippy pass.

P02 accepted: normal-range homo/heterodimer factors, with exact small-count cases
and scaled exceptional fallback. The secondary dimer control gives 1.176x speedup
with identical payloads. The four original workloads do not exercise these factors;
the serial remeasurement is retained separately from an exploratory timing that
overlapped the secondary control. No general speedup is claimed for those rows.

P03 accepted: trusted scheduler firings skip a redundant full-state scan while
keeping propensity enablement and transactional limit/overflow checks. Exact
payloads match; speedups range 1.008–1.046x across the four frozen controls.

P04 accepted: selected positive channels use one checked transactional update,
avoiding duplicate propensity evaluation. Public fire_reaction still validates
untrusted callers. Exact payloads match; paired timings are retained by workload.

P05 rejected as inconclusive: preallocating each saved-path vector passed every
correctness check, but paired speedups 0.984–1.015x showed no consistent benefit.
The original allocation behavior is restored; the candidate diff and all timings
remain available. This experiment is not counted as a speedup.

P06 accepted for latency: one-worker execution now reduces on the caller without
spawning a thread or exchanging each path through a channel. The ordered fold is
shared with the multi-worker route. Long-work timings are roughly neutral; the
21-pair short-path control measures 11.30x core speedup with identical results.

P07 accepted for the eight-worker short-path case: bounded queues hold four paths
instead of one, reducing blocking (1.166x on 21 paired short-path repeats). Workspace
reservation increases from three to six path buffers per worker when workers>1;
the 256 MiB admission cap still applies. One-worker frozen workloads are unchanged
in design and their small timing variations are retained, not claimed as gains.

P08 accepted: strict sum and minimum-positive-rate tracking share one pass. A
single minimum/total resolution comparison replaces per-channel divisions, with
the original ordered sum and lost-increment failure. Identical results; speedups
1.036x two channels, 1.311x sparse, 1.005x circuit and 1.122x dense.

P09 rejected: precomputed Hill logarithms passed exact checks, but the inline
cache enlarged every rate enum and regressed unrelated rates. A boxed cache
variant avoided that layout cost but gave only 1.023x on the circuit while
other medians were 0.972–0.995x. Both timing variants and the final candidate diff
are retained; the simpler P08 rate representation is restored.

P10 accepted: FSP normalizes immutable sparse rows once and reuses all product/
series buffers. The 401-state, horizon-40 control measures 1.187x process speedup;
the full FSP payload (including work/error fields) is identical. All 16 independent
FSP checks are rerun. SSA timing variations in this iteration are not attributed
to an FSP change.

Ten experiments are complete: eight accepted, P05 and P09 rejected with evidence.

Final compatibility correction: queue capacity falls back to one when the larger
reservation would exceed the cap but the legacy reservation fits. An actual
250000-save-point cancelled archive verifies under this path. The optional job
bridge also reserves one supplied CPU allowance for ordered reduction when using
multiple numerical workers. Stable clippy's Option unwrap suggestion is resolved
with an equivalent tuple pattern. All final numerical gates pass on this source.
The final comparison with GitHub main is final-compatible-performance.json;
final-performance-summary.json is derived from that report. Earlier final-vs-main
timings are retained as superseded measurements.
