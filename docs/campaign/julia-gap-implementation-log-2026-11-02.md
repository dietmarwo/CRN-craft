# Julia-gap campaign: implementation log (2 November 2026)

Protocol: [decisions/julia-gap-2026-11-02.md](../decisions/julia-gap-2026-11-02.md).
Evidence: `validation/julia-gap/`. Host: 2 vCPU Intel Xeon (AVX2/AVX-512),
Linux, Rust 1.97 stable, one numerical worker. Base: the imported acceleration
patch (`3184f65`, from `baed0e4`).

## Import and review

The acceleration patch (Astra) applied cleanly as one commit; build, 181 tests,
clippy and the contract checkers passed. Three independent reviews covered the
engine core, RSSACR bounds and Ziggurat, and the Python/validation/docs side.
Confirmed findings and their handling:

| Finding | Severity | Handling |
|---|---|---|
| SortingDirect fails (CRN-E061) when the ordered scan sum rounds one ulp below the tree total | bug, ~1e-16 per event | Fixed: redraw the threshold (counted); unit test with a one-ulp case |
| Alternative-SSA output budget reduced twice (legacy archive allowance, then own check) | medium | Fixed: payload limit is a scheduler parameter; exact-size test |
| Python `ssa` CLI exits 0 for rejected or incomplete runs; rejection is a message string | medium-high | Fixed: native-style diagnostics, exit 1; checker covers both CLIs |
| New methods share legacy streams at equal seeds; Julia-comparison validation counts identical-by-construction rows | medium | Documented; `compare_julia_accelerated.py` now uses distinct seeds per engine |
| RSSACR fails on very steep Hill laws (n≈130) at t=0 | robustness | Documented as a known limitation, with working alternatives |
| Exactness argument undocumented; top bucket bound `f64::MAX` | low | Documented in `alternative-ssa.md` |
| Changelog, discovery version, CLI milestone, notices pointer to a missing file, W6 reuse claim | low | Fixed in docs and strings |

No review found a law error: 210/210 exact-law comparisons, 3.1e9 bounded
integer states and the Ziggurat tables were rechecked independently.

## Why CRN-craft was slower

Reading JumpProcesses.jl 10.0.1, the Catalyst-generated code and callgrind
profiles of the imported engine gave the causes listed in the protocol: an
out-of-line stop poll before every random word, ChaCha8 on its SSE2 fallback,
Hill laws evaluated through `log1p`/`exp`/`log`, strict per-channel checks in
every ordered sum, three fields per rate-tree node and a full ancestor walk per
changed reaction, composition weights recomputed with divisions on every
attempt, and three random words per composition attempt. Julia maintains
incremental sums without drift or resolution checks and draws from Xoshiro.

## Iterations and decisions

Exact iterations had to keep numerical digests (results with work counters
removed) on all 35 rows; A08 changes paths and passed the statistical gate.

| ID | Change | First round | B re-test (targets by ≥5% fewer instructions) | Kept |
|---|---|---|---|---|
| A01 | Poll per loop iteration; inline poller | reject | 18 targets, median 1.105×; W4 `cached_direct` 0.963× over 10 rounds | reverted |
| A02 | ChaCha8 AVX2 dispatch (`rand_chacha/std`) | reject | 10 targets, median 1.082× | **yes** |
| A03 | Hill memo by regulator count | **accept** (1.99×) | 5 targets, median 1.987× | **yes** |
| A04 | Skip redundant RSSACR evaluation | reject | 1 target (W6 1.051×); W4 `cached_direct` 0.960× | reverted |
| A05 | Cached composition envelope prefix | reject | 10 targets, median 1.090× | **yes** |
| A06 | Sums-only rate tree, lazily checked | reject | 13 targets, median 1.308× | **yes** |
| A07 | Additions-only ordered sums | reject | W5/W6 cached rows ≤3/5 faster (sum is latency-bound) | reverted |
| A09 | Batched ancestor rebuild | reject | 5 targets, median 1.221×; W4 `rssa_cr` 0.946×, W7 `direct_cr` 0.949× | reverted |
| A10 | Keep small bucket buffers | reject | no target row | reverted |
| A08 | One word for member and acceptance | — | 8 targets, median 1.035×; 210/210 laws, max CDF error 0.00440 (limit 0.00805) | **yes** |

The first round rejected everything but A03 because its target sets were too
broad (they included rows a change cannot reach). The B re-test was declared
and committed before its timings; it selects targets from deterministic
instruction counts and keeps every threshold. Rejected iterations were
reverted. The two fix commits failed the no-regression rule only on W4
`cached_direct` (0.965×, identical instruction count, code they do not touch)
and were kept because they correct wrong results.

W4 `cached_direct` moved by 3–7% between unrelated builds with unchanged
instruction counts, which points to code-layout sensitivity of its hot loop.
It decided three rejections (A01, A04 and the fix commits). A future protocol
should randomize or average over code layout, or decide on instruction counts
with a separate wall-clock sanity check.

## Results

Cumulative paired speedup from the import to the final tree
(`final-cumulative.json`, five rotated rounds, all 5/5 faster except where
noted): W3 2.0–2.4× on every method (Hill memo); outside W3, DirectCR
1.25–1.91×, SortingDirect 1.07–1.73×, RSSACR 1.22–1.42× (W6 4/5), and legacy
and cached Direct 1.00–1.05× (W5 legacy 3/5, W4 and W5 cached 4/5).

Projected against the comparison report's Julia medians (the Julia package
hosts are blocked here, so this multiplies the report's CRN medians by the
measured speedups; it is not a Julia measurement):

| Workload | Fastest projected CRN | Projected s | Fastest Julia | Julia s | CRN / Julia |
|---|---|---|---|---|---|
| W1 | legacy cached | 0.0749 | Direct | 0.0511 | 1.47× (was 1.54×) |
| W2 | DirectCR | 0.1272 | DirectCR | 0.0968 | 1.31× (was 1.89×) |
| W3 | SortingDirect | 0.3278 | SortingDirect | 0.2944 | 1.11× (was 2.39×) |
| W4 | RSSACR | 0.1946 | RSSACR | 0.1061 | 1.83× (was 2.24×) |
| W5 | DirectCR | 0.0768 | DirectCR | 0.0477 | 1.61× (was 2.20×) |
| W6 | RSSACR | 0.0614 | RSSACR | 0.0396 | 1.55× (was 2.21×) |
| W7 | RSSACR | 0.0187 | RSSACR | 0.0133 | 1.41× (was 1.79×) |

Geometric mean 1.46× (from 2.02×). The declared target of 1.4× was missed.

## What still separates the engines

Final callgrind profiles (path functions only) of each workload's fastest
method:

- **Stop polling, 6–14%.** The out-of-line poll before every random word
  remains because A01 was reverted. Julia has no cancellation checks.
- **Per-leaf ancestor walks, 13–34%** (`RateTree::set`; W5 DirectCR 34%,
  W3 SortingDirect 21%). A09 removed most of this but was reverted. Julia's
  DirectCR updates group sums in O(1) without drift control.
- **Composition sampling and membership, 6–15%** on DirectCR/RSSACR.
- **Rate evaluation inside RSSACR bounds, about 18% on W4**: Hill and
  Michaelis–Menten bounds enumerate up to 129 integer states; Julia evaluates
  two endpoints of user-asserted monotone laws.
- **Ordered Direct sums** (W1 legacy cached, W5/W6 cached): a sequential
  floating-point addition chain per event, kept for byte-identical legacy
  paths; Julia's Direct recomputes too, but its SortingDirect and DirectCR
  avoid ordered sums entirely.
- **Random words**: ChaCha8 (now AVX2) against Xoshiro256++, plus a counted
  budget check per word.
- **SortingDirect order** resets per path for worker-count replay; Julia keeps
  the learned order across a batch.

Some of the gap is policy: checked transactional state updates, explicit
resolution and loss checks (now cheap in the tree), cooperative cancellation
and sampling budgets, and replay that is independent of the worker count.

## Validation of the final tree

`cargo fmt`, `clippy -D warnings`, release tests, the dependency gate,
`check_alternative_ssa.py` (210/210, normal and `-O` byte-identical),
`check_alternative_contracts.py` and `check_contracts.py` (normal and `-O`
identical) passed locally; the full CI matrix runs on the pull request.
Limitations: no Julia rerun on this host; one host and one CPU vendor; the
projection assumes relative speedups transfer to the AMD host and to the
prepared-Python boundary of the original report.
