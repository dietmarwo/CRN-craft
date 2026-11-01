# M2 acceptance targets declared before simulation

Date: 1 November 2026. Base: `bb2b9aed6129734bbe5856a709e0f2d1e32ca198`.

Scope: Direct SSA with pinned ChaCha8 streams and owned uniform/exponential/
selection mappings; ordered ensemble reduction, mean and unbiased sample
variance; event/population/time/output/workspace limits and cooperative
cancellation; versioned results and checksummed, completion-marked archives.
FSP and comparative performance claims remain later milestones.

## Deterministic gates

- Preserve the reviewed 100 Rust regressions, 58 contract, 218 network,
  309 structure and six GillesPy2 checks; preserve all 17 source mutants.
- Pin RNG known-answer vectors, seed/stream mapping and sampler endpoints.
  Compare the ChaCha8 stream with an independent scalar Python implementation.
  Exercise zero rejection, exponential inversion and categorical boundaries.
- Test zero rates, absorption, reaction availability, integer overflow,
  population safety stops, null reactions, finite total propensity, lost small
  positive increments and time stagnation. Preserve CRN-E042 underflow errors.
- Events at a save time are included; pending events are reused across save
  points. Changing save density must not change path endpoints or event counts.
- Completed fixed-work payloads and ordered mean/variance reductions must be
  identical at 1, 2 and 8 threads. Timing/thread provenance is separate.
- Account for every replicate, including never-started work. Incomplete runs
  carry no unconditional exact-ensemble statistics; partial trajectories and
  termination evidence remain inspectable.
- Enforce retained JSON bytes and bounded workspace; test below/at capacity.
  Cancellation/deadline polling occurs during event/rate work and output loops.
  Measure cancellation latency on an active long run; target under one second
  locally, report the measurement and do not claim a real-time guarantee.
- Archive tests: original input bytes, hashes, unique fixed paths, malformed/
  duplicate metadata, missing completion marker, corruption, symlinks, extra
  files, overwrite refusal and byte limits. A complete archive can record an
  incomplete simulation; completion of the write is a separate property.
- Add mutations targeting wait inversion, selection, saves, stream separation,
  event caps, variance and archive verification. A passing baseline is mandatory.
- Run fmt, clippy with warnings denied, release/locked tests on Rust 1.85,
  1.88 and available stable; Python normal/-O equality and dependency checks.
  Apply the delivered patch to a fresh baseline and rebuild/test it.

## Analytic campaign (frozen before first stochastic execution)

Use 100,000 independent replicates per fixture, root seeds 2026110401,
2026110402 and 2026110403 respectively. Save at 0, 0.25, 1 and 2 seconds.
Each fixture has one species X; all rates are stochastic count coefficients.

1. Pure decay: X(0)=20; X -> empty, coefficient 0.7.
   X(t) follows Binomial(20, exp(-0.7*t)).
2. Immigration/death: X(0)=0; empty -> X, coefficient 3;
   X -> empty, coefficient 0.8. X(t) follows
   Poisson(3/0.8 * (1-exp(-0.8*t))).
3. Nonzero immigration/death: same rates, X(0)=12.
   X(t) is the convolution of Binomial(12, exp(-0.8*t)) and the
   independent immigration Poisson law above.

Test all nine positive-time marginal CDFs with DKW and Bonferroni:
family alpha=0.001, m=9, epsilon=sqrt(log(2*m/alpha)/(2*100000)).
Compute exact-form CDFs by binomial/Poisson recurrences and convolution;
no normal approximation. Check the complete empirical integer support and
both tails. t=0 is exact. Use at most 10,000 events per replicate, no deadline
or population stop, and sufficient explicit output budget. Any incomplete
replicate fails the campaign. Preserve failed rows, seeds and fixture hashes;
do not relax epsilon after seeing output.

Mean and sample variance are checked independently from retained paths;
mean SE is sqrt(sample_variance/n), absent for n=1. No variance-SE claim.
Add small independent first-reaction/reference trajectory checks and archive
checker corruption probes; statistical agreement alone does not validate
the PRNG, event loop or evidence format.
