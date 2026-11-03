# CRN-craft fast generic rules implementation log

> Contributed with the patch and kept as written apart from this note. The
> series was built on the pre-review head `6edc100` (tree `861d295`), not on
> `main`; it was applied there and merged with `main` (`f255ee9`). Since then
> Static Direct reports CRN-R031 for a sub-resolution channel, not CRN-R030 as
> stated below. The [review](fast-generic-rules-review-log-2026-11-03.md)
> changed the generated-kernel fingerprint, observer calls and the
> sparse-tracker identity, and added Rust tests for the generated path.

3 November 2026. Implemented S1–S6, including conditional S5. The patch
was tested on local mtDNA base `6edc100c24ce1996ad386d2de25dd6c07b18945e` with rustc 1.88.0 (6b00bc388 2025-06-23). A fresh
GitHub fetch timed out and the merge commit `f255ee9` named by the plan
was absent from the available clone. Patch applicability to that GitHub
merge remains unverified. The verified result tree is `abf845ce89f71830aab3676143dbbb3187612190` at
`915d3a7a16a14f07c9fdada1a7df79f6101948ca`.

## Incremental results

| Milestone | Implementation and gate |
| --- | --- |
| S1 | Cached rule factors and common vacancy bases. Re-evaluated only dependent channels and kept fixed-order total, resolution and draw behavior. Updated live instance accounting from each committed delta. Static Direct, Next Reaction, ATP thinning and the leap rate sweeps share rate caching. Default results passed raw replay. |
| S2 | Refreshed affected projections in declaration order. Deferred certified finite report-only projections. Compiled deterministic predicates to register bytecode with short-circuit branches and a fused guarded fraction comparison. Stop invalidation uses conservative observable dependencies. Potentially overflowing projections retain eager checks. Default results passed raw replay. |
| S3 | Added optional Ziggurat waits and constant-span count-tracker integration to rules-run/0.4. The default sampler and output schema remain stable. The new sampler passed CME, worker replay and tracker reconstruction gates. |
| S4 | Compared both engine pairs at one and eight workers, alternating order with three repetitions. The mtDNA gap exceeded the 1.5× threshold, triggering S5. |
| S5 | Added ahead-of-time Rust generation for closed models with constant or recognized vacancy rates and weighted sums. It emits shared vacancy expressions, ordered projection sums, static dependencies and a topology fingerprint. The generated mtDNA kernel used the normal Static Direct scheduler and replayed its output byte for byte. |
| S6 | Added `docs/specialized-kernels.md` and linked it from `AGENTS.md`. The report-at-event rule follows the existing right-continuous contract. The new generic path handles an overflowing wait as reaching the horizon. |

## Timing and decision

These runs used the same host and an eight-CPU quota. The mtDNA fixture
uses 512 cells, seed 17, a 25,920-hour horizon, 720-hour reports, the 90%
takeover stop and two scalar count trackers. The sparse fixture has 48
species, 64 cells and a 100-hour horizon. S1 and S2 rows are single
observations. The first S1 variant was slower on mtDNA and the refined
variant returned approximately to baseline. The 1.6× mtDNA target was
missed. The 3× sparse target was met in the refined measurement. S2
improved mtDNA but fell short of the additional 1.4× target.

| Step | mtDNA, one worker | Speedup | Sparse, one worker | Speedup |
| --- | ---: | ---: | ---: | ---: |
| S1 initial | 13.046 → 13.359 s | 0.98× | 1.161 → 0.362 s | 3.21× |
| S1 refined | 12.922 → 12.922 s | 1.00× | 1.067 → 0.350 s | 3.05× |
| S2 against S1 | 13.434 → 10.819 s | 1.24× | 0.323 → 0.304 s | 1.06× |

S4 below reports median child-process wall time. The mtDNA generic row
uses Ziggurat and sparse trackers, with retained paths on both sides.
The sparse row compares the cached inverse Static Direct path with the
updated Next Reaction path. Raw evidence also retains the earlier
comparison against baseline Next Reaction.

| Workload | Workers | Generic Static Direct | Reference | Reference time | Generic / reference |
| --- | ---: | ---: | --- | ---: | ---: |
| mtdna | 1 | 9.170 s | Specialized mtDNA | 2.753 s | 3.33× |
| sparse | 1 | 0.302 s | Next Reaction | 0.150 s | 2.01× |
| mtdna | 8 | 1.489 s | Specialized mtDNA | 0.440 s | 3.39× |
| sparse | 8 | 0.055 s | Next Reaction | 0.028 s | 1.95× |

The S5 comparison retained count trackers on both sides and used the
same sampler, request and output bytes. Each cell in this table is a
median of three alternating runs. The ahead-of-time probe build took
0.79 s using a warm Cargo cache. The first
replay attempt exposed different JSON field ordering between the
standalone probe and CLI. Serializing both through the CLI's ordered
JSON value representation made the raw-byte gate comparable.

| Workers | Interpreter | Generated | Generated speedup |
| --- | ---: | ---: | ---: |
| 1 | 9.149 s | 8.267 s | 1.11× |
| 8 | 1.536 s | 1.394 s | 1.10× |

The generated speedup was modest and did not close the gap to the
specialized kernel. Keep the generated route opt-in and
`mtdna_direct` available for fixed mtDNA throughput. Generation binds
the exact model bytes, including parameters, so binding a different
candidate requires regenerating the source. Other rate and projection
families continue through generic Static Direct.

## Validation

- 228 Rust release tests passed. Workspace formatting and Clippy
  with warnings denied passed. The expression comparison covered
  512 seeded expressions across four input states, signed zero, NaN,
  short-circuit branches and division by zero. The deterministic
  predicate language rejects random `choice`.
- 34 one/eight-worker default replay cases
  matched raw CLI output, stderr and exit status against the local
  baseline. These include overflow, unresolved positive propensities,
  instance budget and undefined stop diagnostics. The local baseline
  uses CRN-R030 for unresolvable totals. The patch preserves that code.
- The new sampler's independent K=20 master-equation gate used 32,768
  paths, family alpha 0.001 and the existing bounds. Its largest CDF
  distance was 0.00609 against 0.01306. Count integrals
  reconstructed from a full event log agreed within 1e-12 relative.
- Python wheel, callback identity and dependency checks passed. The CLI,
  rules, mtDNA and Python gates agreed under normal Python and Python
  `-O`. The generated runner also rejected a changed model fingerprint.

Raw results are in `validation/mtdna/fast-generic-rules.json` and
`validation/mtdna/gates-fast-generic.json`. The ZIP has a single patch,
an ordered commit series, instructions and checksums.
