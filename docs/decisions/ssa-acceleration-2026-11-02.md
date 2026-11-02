# Five-step exact SSA acceleration campaign: predeclared protocol

Date: 2 November 2026. GitHub base: baed0e4590f93e4c1dcc77cd7f039a2604a1594c.
The working branch also contains the preceding Julia comparison at dbb6809.
No solver is published or pushed as part of this campaign.

## Iterations

1. Measure existing whole-run, preparation and isolated sampling/rate costs.
   Expose the already-existing Rust prepare/execute split through Python.
2. Implement composition–rejection Direct using existing reaction dependencies.
3. Implement Sorting Direct; evaluate a versioned exponential Ziggurat sampler
   separately from inverse-transform sampling.
4. Implement rejection SSA with composition–rejection selection and independently
   checked bounds, including availability gates, rejected-work limits and stops.
5. Repeat the same-host Julia comparison and add larger connected networks.

Each iteration gets a release build, targeted correctness checks, measurements
where applicable, a log entry and a commit. Failed development runs are retained.
The old model/run/result/archive schemas, RNG vectors and Direct/cached results
must remain unchanged. New exact methods use separate request/result schemas
and explicit algorithm/sampler identities. Deterministic paths must agree across
worker counts within each fixed new method, but not across different samplers.

## Numerical design requirements

The new methods must retain validated inputs, full custom propensities, checked
transactional count updates, zero-rate absorption, save-time/event ordering,
population/event/output/workspace limits and cancellation/wall-time handling.
Every rejection loop, including exponential sampling, must have a bounded work
counter and cooperative stop checks. Incomplete ensembles withhold statistics.
The new exact command returns its own JSON result; legacy sealed Direct archives
and their verifier are preserved. Do not label new output as a legacy archive.

Use published algorithms, cite their papers, and retain any adapted MIT notices.
Do not change RNG or floating-point behavior underneath the existing methods.
Use deterministic balanced sums or an explicitly checked numerical policy for
new selection structures. No silently clipped negative totals or discarded
positive channels. Selectors must handle activation/deactivation, equal rates,
extreme dynamic ranges, null reactions, and rare rejection/tail branches.

RSSACR bounds must cover the actual gated rate evaluator on each integer state
box. For built-in one-variable custom laws, enumerating a bounded-width integer
interval is acceptable and avoids assuming perfect floating-point monotonicity.
For elementary mass action use monotone fast paths only where their arithmetic
is proved applicable; otherwise pin affected coordinates to their current
values and re-evaluate on changes. Unsupported or unrepresentable bounds must
produce an explicit failure or documented exact fallback, never an invalid
bound. Include gates involving a distinct species and repression involving its
own consumed regulator in independent bound tests.

## Correctness gates (fixed before experiments)

- Preserve the baseline Rust regressions, schema snapshots, RNG vectors and
  legacy Direct/cached path identity. Run fmt, clippy -D warnings, release tests,
  schema/CLI checks normally and with Python -O, and the dependency gate.
- Validate all new algorithm/sampler combinations against independent linear
  laws and the four existing finite-state/nonlinear CME oracles. Reuse 100,000
  paths per fixture, the established seeds and save times. Include Direct,
  cached Direct, DirectCR, SortingDirect and RSSACR with inverse and Ziggurat
  waits: 10 variants * (9 linear + 12 nonlinear marginals) = 210 checks.
  Use family alpha 0.001 and DKW sqrt(log(2*210/0.001)/(2*100000)). Retain
  power checks rejecting the wrong homodimer convention. Never loosen this gate.
- Additional categorical and exponential tests: 1,000,000 samples, fixed seeds;
  DKW bound sqrt(log(2*100/0.001)/(2*1000000)) for at most 100 CDF comparisons.
  Exponential mean/second-moment diagnostics are descriptive, not retuned gates.
- Independently enumerate small state boxes to check each reported bound against
  every enclosed gated propensity. Exercise activation and deactivation at zero,
  mixed gate/regulator directions and numerical edge cases. Test proposal/work
  exhaustion and cancellation inside rejected work, and thread replay.
- Python prepared and one-shot interfaces must return identical numerical
  results for the same request. CLI schemas must validate actual new results,
  reject unknown/duplicate fields and identify incomplete outputs honestly.

## Profiling and performance

Record isolated component timings separately from whole-run timings. Component
microbenchmarks include their public-API validation overhead and are diagnostic,
not a decomposition whose percentages may be added. Try hardware profiling only
if available; do not infer exact cycle attribution from end-to-end timings.

Retain W1–W4 unchanged: 64 paths, horizon 100, saves [0,100], seed 4242 and
one numerical worker. Five rotated/interleaved timing rounds, with all builds
and heavy checks finished first. Include baseline Direct/cached, each new
algorithm/wait-sampler pair, prepared API measurements, and the same five Julia
algorithms. Julia 1.13.3 / Catalyst 17.0.0 / JumpProcesses 10.0.1 stay pinned.
Use fresh copies of pristine Julia problems outside the Julia warm timer.
Record medians, minima/maxima, all observations and own-engine replay digests.
Keep one-shot and prepared timing boundaries explicit. Event counts are measured
for Rust; Julia event-normalized values remain labelled CRN-count proxies.

Add three synthetic but connected workloads: a 256-species nearest-neighbor
diffusion chain (510 channels), a heterogeneous 32-state fully connected
conversion network (992 channels), and a 32-compartment reversible binding and
substrate-diffusion network (96 species,126 channels). Use 16 paths, horizon 20,
saves [0,20], five rounds. Freeze model JSON/hashes before timing. Compare 512
validation paths on all seven workloads before timing: linear W1/W2 marginals
against exact Poisson laws, other terminal marginals against legacy cached SSA,
with the same union-bound DKW formula family alpha 0.01 and the actual declared
number of species/engine comparisons computed from the frozen fixture set.

Success means correct runnable implementations and an honest measured report.
Matching the fastest Julia result is a target, not an acceptance claim or reason
to relax numerical checks. Retain slow variants and explain their costs; report
any unsupported cases and limits. Ship a patch ZIP, implementation log, updated
comparison, raw evidence, pinned environments and reproduction instructions.
