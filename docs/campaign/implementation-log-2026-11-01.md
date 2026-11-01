# CRN-craft implementation and performance campaign

Date: 1 November 2026 (UTC). Repository: https://github.com/dietmarwo/CRN-craft.

The remaining milestone campaign is implemented and tested, followed by ten
performance experiments. Eight performance candidates were accepted; two were
rejected and reverted with their diffs and measurements retained. The deliverable
is a local patch against reviewed GitHub main, with evidence and an installable
Linux wheel. No branch, tag, release or package index was published.

## Starting point and execution record

GitHub main was fetched and inspected at
`5a0733444aa46c8cea2a5b61bd8601daf5890d79`, the merged iteration-4/M2 review.
The supplied iteration-4 review was read. M0 through M2 were already merged;
the reviewed plan identified M3 as next. A separate local campaign branch was used.
The final remote-main observation is retained as `github-main-final.log`.

`acceptance.md` was committed before numerical development. Each functional
iteration and each P01–P10 experiment has a separate local commit. Corrections
discovered by testing were committed afterwards. The ZIP's `APPLY.md` and
`PACKAGE-MANIFEST.json` give the exact base, final source commit and expected Git
tree. Use either the consolidated patch or the ordered series, once, on that base.

All evidence paths below are relative to `validation/campaign/` in the patched
repository. The ZIP also includes a copy under `evidence/campaign/`. Earlier failed
and superseded rows remain available; the final report names below identify the
accepted results. The last runtime-source commit is `70a5d64`; the numerical core
has no changes after the final compatibility correction, `a4cb4c4`.

## Functional implementation loop

| Stage | Implemented behavior | Executed acceptance |
|---|---|---|
| M3 | Reachable sparse FSP, explicit escaping sink, split uniformization, separate leakage/series/roundoff accounting, work/state/edge/memory/output/deadline/cancellation limits | 16 independent matrix-exponential, analytic and selected DSMTS checks; maximum per-state/sink discrepancy 4.07e-14 against the frozen 2e-12 threshold |
| M4a | Frozen Direct workloads and pinned same-host reference harness; separate process/core timing and recorded event/output work | Five paired dynamic-reference repetitions with rebop 0.10.0 and GillesPy2 1.8.3 NumPySSASolver; unavailable Julia and physical-core scaling are recorded |
| M5 | Isolated PyO3 binding and Python package; native-backed CLI, strict byte/hash handling, trajectory autocorrelation/spectrum and saved-grid censored first passage | CLI/binding equality and 1/2/8-worker identity; seven grouped checks; synthetic development/held-out controls and 16 three-gene paths; normal/-O agreement |
| M4b | Exact dependency-cached Direct, retaining ordered sums, channel selection and RNG mapping; external cooperative cancellation hook | Identical frozen trajectory digests and full analytic SSA/nonlinear CME reports for original and cached implementations |
| M6 | Native fcmaes worker registration for fit-moments/circuit-score, log coordinates, fixed candidate seed sets, real SDK lifecycle, bounded workers and failed/incomplete-objective rejection | 15 optimization/evaluate jobs, 16 manifests including a command job, active deadline, explicit failure accounting, disjoint-seed validation and 1/8-thread equality |
| M7 | Optional diffsol BDF, count-space kinetics with analytic state/parameter Jacobians, forward sensitivities, BiteOpt-to-SQP fitting and conservative topology screening | 16 analytic/DOP853/sensitivity comparisons; global/local synthetic fit; both screening topologies receive held-out SSA allocation |
| M8 | Provenance-preserving explicit SBML L3V2 and expanded BioNetGen .net subsets, exact rational propensity checks, explicit rate convention, coded rejection of unsupported semantics | 29 importer checks normally and under -O, including undefined semantics, unsafe XML, homodimer factors, coefficient overflow and positive underflow |
| M9 | Separate explicitly approximate adaptive Poisson-leap request/result, critical low-count SSA fallback, reject/halve/fallback handling, no negative-population clipping | 30,000 paths over epsilon 0.1/0.05/0.025; all mean/CDF envelopes pass; request, budget, cancellation, reproducibility and critical-reaction checks |
| M10 | ABI3 wheel, source distribution, clean installed CLI/API acceptance, resources/license inclusion, optional platform workflow and owner-controlled release preparation | Linux x86_64 wheel installed on Python 3.12; extracted sdist rebuilt and installed; seven grouped installed checks pass, including CLI export/legacy analysis/tau/schema selection |

M4b was selected from the measured full-propensity-evaluation bottleneck.
SortingDirect, NRM and RSSA variants were not added without a further workload
need. This is the demand-driven branch of the reviewed plan, not evidence that
every advanced aggregator is unnecessary.

FSP enumerates reachable states up to 100,000 states and 2,000,000 edges, charging
up to 200,000,000 work units with a conservative 256 MiB workspace cap. It splits
intervals so the uniformization mean is at most 16. Escaping probability goes to
the sink; it is not reflected or renormalized. A geometric series-tail majorant
is reported separately from ordinary floating-point mass residuals. Roundoff is
not certified, and the retained unnormalized moments do not bound unbounded
physical moments. DSMTS inputs were read externally at pinned revision
`e79fd58e94757fb39926c6373297c1f2a1845894`; their hashes are in the report. LGPL
reference files were not copied into the MIT engine.

The original eleven schema snapshots are byte-for-byte unchanged. Discovery is
now 0.5; FSP result and tau request/result have separate 0.1 schemas. Exact SSA
result/0.2 and archive/0.1 remain intact. Incomplete ensembles withhold statistics.
`original-schema-compatibility.json` records all eleven source hashes.

## Ten subsequent performance experiments

All ten experiments followed functional implementation. P00 is the functional
campaign baseline; GitHub main is the overall baseline. Frozen exact workloads
use seed 4242, 64 paths, horizon 100, saves at 0 and 100, one worker and seven
alternating paired repetitions. Every accepted numerical payload must match.
Short-path controls use 21 pairs; the FSP control has its own paired report.

| Experiment | Decision and evidence |
|---|---|
| P01 | Accepted zero/first-order mass-action fast paths with scaled exceptional fallback. Per-step paired gains 1.010–1.152x across the four controls. |
| P02 | Accepted checked elementary second-order fast paths. The dimer control improves 1.176x; the four original workloads do not establish this change's benefit because they lack these factors. |
| P03 | Accepted skipping redundant full-state validation on trusted scheduler firings, preserving checked transactional updates. Measured gains 1.008–1.046x. |
| P04 | Accepted avoiding duplicate selected-propensity evaluation. Public reaction firing retains validation. Exact payloads and gates pass. |
| P05 | Rejected trajectory preallocation: paired gains 0.984–1.015x are inconclusive. Original allocation restored; candidate diff and timings retained. |
| P06 | Accepted caller-thread reduction for one-worker execution. Short-path core latency improves 11.30x; this is not a process-startup or long-path speedup claim. |
| P07 | Accepted capacity-four bounded queues for the eight-worker short-path case: 1.166x. Larger reservation falls back to legacy capacity one when needed to preserve archive admission. |
| P08 | Accepted fused ordered propensity sum/minimum-resolution tracking. Measured gains 1.005–1.311x, with the same ordered sum and failure behavior. |
| P09 | Rejected two Hill-log caches: inline enum growth regressed other rates; boxing yielded only 1.023x on the circuit and 0.972–0.995x elsewhere. Both variants remain in evidence; original representation restored. |
| P10 | Accepted pre-normalized sparse FSP rows and reused product/series buffers: 1.187x process speed on the 401-state/horizon-40 control, with the complete FSP payload identical. |

The P07 compatibility correction was tested with an actual cancelled archive
containing 250,000 save points. Multi-worker optional jobs reserve one supplied
CPU allowance for ordered reduction. The final numerical gates were rerun after
these corrections.

Final measurements compare original GitHub-main Direct with final cached Direct:

| Frozen workload | Median paired core speedup | Paired bootstrap 95% interval | Before → after ns/event |
|---|---:|---:|---:|
| Two channels | 1.551x | 1.496–1.621x | 103.28 → 66.48 |
| Sparse 66 channels | 4.338x | 4.274–4.703x | 1169.89 → 266.22 |
| Three-gene circuit, 12 channels | 2.929x | 2.893–3.058x | 568.04 → 190.49 |
| Dense 66-channel control | 1.601x | 1.539–1.634x | 1181.11 → 746.42 |

All event counts and numerical digests match. These overall gains include M4b;
they are not the product of ten isolated claimed improvements. Raw paired
core/process times and ranges are in `final-compatible-performance.json`;
`final-performance-summary.json` uses 20,000 bootstrap resamples with seed
20261101. The intervals describe seven pairs on this shared host, not general
hardware uncertainty. Earlier `final-vs-main.json` timings are superseded.

The host reports AMD EPYC 9V74, nine available virtual CPUs and an eight-CPU-equivalent
quota. It does not establish eight dedicated physical cores. Benchmarks were
serialized. Rust 1.88 release builds were used for these timings.

Pinned external timings are intentionally qualified: rebop/CRN core ratios in
the final reference remeasurement are 0.573 on two channels and 6.036 on sparse
66 channels. Thus CRN remains about 1.74x slower on the tiny reference workload
and is about 6.04x faster on the sparse one. These rows include Python/xarray
overhead for rebop and ordered statistics for CRN, with different sample streams
and uninstrumented reference event counts. They precede the memory-admission
correction, which does not alter the measured one-worker numerical path. The
GillesPy2 comparison uses NumPySSASolver, not its compiled C++ solver. There is
no universal simulator-ranking claim. See `m4a-references.json` and
`final-rebop-reference.json`.

## Final verification

| Gate | Result and final evidence |
|---|---|
| Native Rust | 142 tests pass on Rust 1.85, 1.88 and stable; fmt and clippy with warnings denied pass. Use `final-rust-185-*`, `final-compatible-*` for 1.88, and `final-rust-stable-*`. An earlier partial `final-rust-188-tests.log` is superseded. |
| Fresh patch application | A clean checkout applied the source series, matched the expected tree and passed all 142 tests, clippy and release build: `source-patch-application.json`, `fresh-*.log`. Final packaging additionally checks both complete patch routes. |
| Mutation resistance | Baseline passes and all 43 mutants are caught: `final-mutations.json`. |
| Independent native acceptance | Contracts 67, network 219, structure 309, analytic SSA 17, nonlinear CME 14, FSP 16 and tau three statistical rows pass: `verified-*.json`. All seven checker families also pass under -O with identical reports. |
| Exact cached SSA | Original/cached SSA and CME reports agree byte-for-byte, normally and under -O: `cached-quality-status.json` and `cached-*.json`. Each implementation's declared campaign has 700,000 analytic/CME paths; repeated gates are not extra independent evidence. |
| RNG and export | Nine streams/360 RNG words pass; GillesPy2 export passes six checks and 15,717 independent propensity comparisons, normally and under -O. |
| Optional integrations | `final-jobs.json` passes; four Rust-generated typed schemas and ten example inputs validate. `final-ode.json` passes 16 comparisons. Optional builds/clippy pass. |
| Python and packaging | Seven Python groups and 29 importer checks pass normally/-O. Clean wheel and rebuilt-sdist installations pass seven groups: `parity-installed.json`, `parity-sdist-installed.json`; source identity is recorded in `parity-source-identity.json`. Binding fmt/clippy pass. |
| Dependency/schema policy | Default-workspace dependency audit and its five mutation probes pass; optional dependencies remain isolated. All original schema bytes are unchanged. |

FSP's maximum discrepancy is 4.07e-14. Tau CDF discrepancies are 0.0238184,
0.0089739 and 0.0071184 at epsilon 0.1, 0.05 and 0.025, within the predeclared
approximation-plus-sampling envelopes. No negative populations occur. These
small controls do not establish stiff or arbitrary-network tau accuracy.

The global/local fit uses synthetic birth/death truth (3.0, 0.8), 64 BiteOpt
calls and active-set SQP. Loss falls from 0.5901 to 4.97e-20, with 94 native
propagations including derivative checks; independent analytic loss is
6.65e-19. `m7-global-local.json` records work and solver status. The ring and
Vilar controls each receive eight held-out SSA paths; many classifications are
ambiguous at the measured horizon. `m7-screen.json` reports these outcomes.

Failures that led to corrections remain in the record: invalid initial FSP
fixtures, a package-source glob omitting nested files, a version-specific archive
tamper fixture, unfixed random probes in SciPy's norm estimate, SDK discovery and
budget assumptions, SQP's `active-set` spelling, stable-clippy issues and importer
literal underflow. No failed performance candidate is counted as an improvement.
Final workflow inspection also corrected required FSP/tau binary arguments.
Hosted GitHub Actions were configured but not executed during this local campaign.

## Dependencies, reproduction and scope limits

The default core and CLI keep their reviewed third-party dependency set and
Rust 1.85 minimum. Optional pins are PyO3 0.29.0; diffsol 0.18.0 with Rust 1.89;
fcmaes source `583f33a9a815ede0e1115c8ddd69029a874593f9`; worker SDK
`602d31afacf8bd2ecb0b0e941d18bd8e2821c1e3`; tested runner
`cc6b03d594f73d1f3b38f088e7e3bd8663814d75`. The Python fit uses fcmaes-rust
0.2.1 and SQP-craft 0.17.0. Exact lockfiles are included.

After applying the patch, representative local verification commands are:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
cargo build --release --workspace --locked
python -m pip install -r validation/requirements.txt -r validation/reference-requirements.txt
python validation/check_dependencies.py
python validation/check_contracts.py --output contracts.json
python validation/check_network.py --output network.json
python validation/check_structure.py --output structure.json
python validation/check_ssa.py --output ssa.json
python validation/check_cme.py --output cme.json
python validation/check_ssa.py --command simulate-cached --output cached-ssa.json
python validation/check_cme.py --command simulate-cached --output cached-cme.json
python validation/check_fsp.py --binary target/release/crn-craft --output fsp.json
python validation/check_tau.py --binary target/release/crn-craft --output tau.json
python validation/mutation_check.py --output mutations.json
```

Pass `--dsmts PATH` to include the external source-file comparisons after checking
out the recorded DSMTS revision. Repeat Python checkers with `python -O` to compare
reports. Mutations build disposable source copies. Independent ODE/worker checks
take explicit binary paths; see their `--help` and M6/M7 notes.

For package verification, install the supplied wheel in a fresh Python 3.12
environment, then run `python validation/check_installed.py --output installed.json`.
Rebuild with maturin 1.9.6 using `maturin build --release --locked --out dist`
and `maturin sdist --out dist`. The supplied ABI3 wheel targets Linux x86_64,
glibc 2.34+ and Python 3.10+; it was executed on Python 3.12. Pinned optional
fitting extras require newer Python according to their package metadata.

Remaining external validation is explicit: Julia is unavailable here; macOS and
Windows wheel jobs have not been run; dedicated physical-core and home-cluster
measurements need the corresponding hardware. Importers implement documented
numerical subsets rather than full SBML/BioNetGen conformance or rule expansion.
The optional ODE path rejects unsupported gates and singular domains. Tau is an
experimental approximation whose small-mean product Poisson sampler still costs
work proportional to expected firings. Biological detector accuracy, general
95%/99% sensitivity targets and parameter identifiability are not established by
these synthetic controls. Publishing remains an owner-controlled later action.
