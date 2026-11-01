# Iteration 4 / M2 validation

Date: 1 November 2026. Base: `bb2b9aed6129734bbe5856a709e0f2d1e32ca198`.
Untouched baseline: 100 Rust tests and 58 contract checks passed locally.

| Gate | Executed result | Evidence |
|---|---|---|
| Rust 1.85.0 / 1.88.0 / stable 1.99.1 | 126 release tests each; fmt and clippy -D warnings | `rust-*-tests.log`, `rust-*-fmt.log`, `rust-*-clippy.log` |
| Contract/schema/CLI | 61 passed; normal/-O identical | `contracts*.json` |
| M1a network/rate oracle | 218 passed; unchanged fixture stream; normal/-O identical | `network*.json` |
| Reviewed M1b structure oracle | 309 passed; normal/-O identical | `structure*.json` |
| Actual GillesPy2 1.8.3 export | Six checks, 15,717 comparisons; normal/-O identical | `export*.json` |
| M2 analytic/trajectory/CLI suite | 17 checks; normal/-O identical | `ssa.json`, `ssa-optimized.json` |
| Independent scalar ChaCha8 | 360 words / nine streams | `rng*.json`, `../rng-vectors.json` |
| Source mutants | 26/26 caught after a passing baseline | `mutations.json` |
| Dependency / Python quality | Passed; dependency normal/-O identical | `dependencies*.json`, `python-*.log` |

Test breakdown: 100 reviewed tests, three new unit tests, 14 SSA integration
tests and nine archive tests. An oversized-file regression was added after the
mutation run; the implementation was unchanged and all final release suites
include it. Source tests exercise exact boundaries and failure paths, not just
agreement with statistical output.

## Frozen statistical family

The acceptance plan precedes simulations in
`docs/decisions/iteration-04-validation.md`. Three fixtures, 100,000 independent
paths each; seeds 2026110401/02/03; times 0, 0.25, 1 and 2. Pure decay is
binomial; zero-initial immigration/death is Poisson; nonzero-initial
immigration/death is a binomial/Poisson convolution. The nine positive-time
CDF comparisons use family alpha 0.001 and DKW/Bonferroni epsilon
0.006999331052635781. Largest observed error: 0.003334014071384228. All paths
completed and all nine checks passed. Total events: 4,263,447.

Each marginal's histogram and analytic error is retained in `ssa.json`, along
with exact model/run/result hashes. Full trajectories are generated and checked
in temporary archives; the report preserves reproducible sufficient evidence
for the marginal CDF/statistic checks without committing large path dumps.
Mean/variance/SE agree with separate integer reductions. A zeroed solver fails
the CDF checker for all three fixtures. No tolerance was relaxed.

The separate Python event-history reference matched 64 paths and event counts
at 1, 2 and 8 threads. The enzyme conservation checks also passed. Incomplete
runs have full replicate accounting and no unconditional ensemble statistics.
Normal/-O result reports are identical; timed stop measurements are separate.

## Cancellation, integrity and environment

Final active CLI cancellation latency: approximately 5.3 ms normally and
5.9 ms under `-O`, below the predeclared one-second local target. See
`cancellation-latency*.json`; differing event counts are expected and excluded
from fixed-work comparisons. The live wall-time stop passed as well.

Archive tests cover original input bytes, hashes/lengths, completion markers,
strict parsing, duplicate paths, symlinks, extra/missing files, false accounting
and statistics, overwrite refusal, output limits and oversized payloads. They
do not establish authenticity or secure concurrent access to a hostile directory.

`environment.json` records Rust/Python/package versions, CPU quota/affinity,
RNG source digests and the inspected sibling repository. RNG dependencies and
licenses are in `rng-dependencies.json`. Local execution was Linux x86_64;
GitHub workflows are configured, not claimed as a hosted run of this patch.

## Retained issues

The initial Python lint failure is retained in `development-python-lint-01.log`.
It concerned strict zip and loop-variable style, not numerical logic. Workspace
build-output failures are retained in `development-environment-*`: a zero-byte
non-executable test binary and a zero-length object mapping error. Clean
temporary target directories passed. No numerical or structural fixture failed.
`ssa-first.json` preserves the first successful campaign before the additional
enzyme acceptance row. `mutations-first.json` preserves the earlier successful
26-mutant run; the final gate includes the scheduler's validated-state fast path.

The mail patch is applied and rebuilt separately after committing source; its
post-commit tree hashes and fresh validation are in the downloadable ZIP.
The project now executes Direct SSA, but no controlled simulator performance
comparison or physical-core scaling result is claimed. M3 precedes M4a.

## Repository review, 1 November 2026

The delivered mail patch `0001-m2-direct-ssa-and-archives.patch` (SHA-256
`3fc39c3d…`) applied cleanly to `bb2b9ae` and produced tree `b977d47`. The
delivery's gates reproduced on Linux x86_64 (two CPUs) with Rust 1.97.0,
Python 3.13.18 and the pinned GillesPy2 1.8.3 stack:

- 126 release tests;
- 61 contract, 218 network, 309 structure, 6 reference-export and 17 SSA
  checks and the 360-word RNG check, each report byte-identical to the
  delivered one and under `-O`;
- the dependency gate;
- 26 of 26 mutants caught after a passing baseline.

Independent checks run in review:

| Check | Result |
|---|---|
| Exact master equation (`check_cme.py`, rule declared in `docs/decisions/iteration-04-review-cme.md` before the first run) | Dimerization with a null reaction, heterodimer binding, a Michaelis–Menten/Hill loop (n = 2.5) and immigration with homodimer annihilation, 100,000 paths each. All 12 DKW tests pass: largest error 0.00245 against ε = 0.00710. Laws with the falling-factorial convention (no 1/2 for 2A) are rejected for both homodimer fixtures, with errors up to 0.445 and 0.372. `review-cme.json` holds the report |
| Oracle accuracy | Uniformization against `scipy.linalg.expm`: largest difference 1.47e-12, above the declared 1e-12 agreement; this row failed and is retained. Against a 40-digit mpmath matrix exponential, uniformization errs by at most 6.0e-13 and expm by 3.4e-14, both nine orders of magnitude below ε |
| Thread invariance | The four fixtures with 2,000 replicates at 1, 3 and 8 threads give identical results apart from the run hash |
| Polling cost | Delivered code read the clock (when a wall-time limit was set) and stat'ed the cancel file at every poll, three or more per event. Median simulation time per event of seven runs, two-reaction birth/death model, delivered (fixed): no limit 95–110 ns in both, wall-time limit 234 ns (110), cancel file 1224 ns (106), both 1238 ns (113); with 66 reactions, cancel file 2202 ns (749). Fixed in review; event counts unchanged |
| Mutation probe of 24 rules in the delivered sampler, scheduler, archive and CLI code | 14 caught by the delivered tests. Survivors: lost positive increment accepted, zero wait accepted, one extra retained trajectory, archive verification without the byte budget, the unstarted-record rules, the initial-state check or the saved-state-after-stop check, the completed-record stop-time check, the pre-write self-check, and `simulate` exiting 0 for incomplete runs. New tests catch eight. The completed-record stop-time check is implied by the saved-state rule, and the pre-write self-check cannot fail for a correct engine |
| Mutation probe of the review's poller | Skipping the first clock/file read, dropping the kept stop and inverting the file interval are caught. A file found by one worker not cancelling the others survives: each worker still finds the file within a millisecond |

After the review:

- 133 release tests: 4 unit, 10 archive, 35 contract, 33 network, 18 SSA,
  22 structure, 11 CLI;
- all earlier Python reports unchanged byte for byte; the CME check passes 14
  rows, identical under `-O`;
- active CLI cancellation took 6.3 ms normally and 5.2 ms under `-O`, after
  about 723,000 and 690,000 events (the delivered binary managed about 66,000
  events in the same window on this host);
- the mutation gate catches 35 of 35 mutants after a passing baseline.
