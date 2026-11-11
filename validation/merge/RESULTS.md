# 0.7.0 engine merge: declarations and results

Task M0c of `docs/plan-0.7.0-engine-merge.md`. This file was committed with
`models.json` **before** any reference or baseline run. Only two sizing probes
ran before it (per-path cost of `direct` and `cached_direct` at the timing
horizons, and FSP feasibility of the bounds below); they are listed in
`docs/campaign/merge-0.7/M0c-log.md` and produced no stored result.

Later tasks append their rows to the results section at the end. Failed rows
are kept. No tolerance below is changed after a result is seen.

## Old-engine build (0.6.0)

Built from branch `merge/M0c` (0.6.0 code plus the plan) with the CI
commands, from the repository root:

```sh
cargo build --release --workspace --locked
maturin build --release --locked --out dist
```

- CLI: `/home/claude/m07/bin060/crn-craft` (copy of `target/release/crn-craft`).
- Wheel: `/home/claude/m07/bin060/crn_craft-0.6.0-cp310-abi3-manylinux_2_34_x86_64.whl`.
- Python: `/home/claude/m07/venv060` (`uv venv -p 3.12`; the wheel, numpy
  2.3.5, pandas 3.1.1, xarray 2026.9.0, jsonschema 4.26.0, rebop 0.10.0,
  ruff 0.11.13).
- Kinetics printer: `validation/merge/kinetics-reference` (a standalone Cargo
  package against the 0.6.0 `crn-craft` crate), binary copied to
  `/home/claude/m07/bin060/kinetics-reference`.

Both scripts take the CLI and the Python interpreter as arguments, so the
0.6.0 engine can be rerun on this host after the network engine is deleted.

## Models

`models.json` is the machine-readable form of this table. Every model is a
`crn-craft.model/0.1` file; its SHA-256 is in `models.json`. Report times
always start at 0 and end at `end` (the 0.6.0 engine requires both).

| Model | File | Species | Reactions | Law `end` | Report times | Replicates | Chunk | FSP |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| W1 | `validation/ssa-acceleration/models/W1.model.json` | 1 | 2 | 10 | 0, 0.5, 1, 2, 5, 10 | 200,000 | 10,000 | X ≤ 250 |
| W2 | `validation/ssa-acceleration/models/W2.model.json` | 33 | 66 | 20 | 0, 1, 2, 5, 10, 20 | 200,000 | 10,000 | no |
| W3 | `validation/ssa-acceleration/models/W3.model.json` | 6 | 12 | 30 | 0, 2, 5, 10, 20, 30 | 200,000 | 10,000 | no |
| W4 | `validation/ssa-acceleration/models/W4.model.json` | 8 | 66 | 20 | 0, 1, 2, 5, 10, 20 | 200,000 | 10,000 | every Xi ≤ 1000 |
| W5 | `validation/ssa-acceleration/models/W5.model.json` | 256 | 510 | 2 | 0, 0.25, 0.5, 1, 2 | 200,000 | 5,000 | no |
| W6 | `validation/ssa-acceleration/models/W6.model.json` | 32 | 992 | 2 | 0, 0.25, 0.5, 1, 2 | 200,000 | 10,000 | no |
| W7 | `validation/ssa-acceleration/models/W7.model.json` | 96 | 126 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | no |
| enzyme-mechanism | `models/examples/enzyme-mechanism.model.json` | 4 | 3 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | S ≤ 50, E ≤ 10, ES ≤ 10, P ≤ 50 |
| immigration-death | `models/examples/immigration-death.model.json` | 1 | 2 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | A ≤ 150 |
| michaelis-menten | `models/examples/michaelis-menten.model.json` | 2 | 1 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | S ≤ 50, P ≤ 50 |
| repressilator | `models/examples/repressilator.model.json` | 6 | 12 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | no |
| vilar-noise | `models/examples/vilar-noise.model.json` | 9 | 16 | 10 | 0, 1, 2, 5, 10 | 200,000 | 10,000 | no |

All twelve models get 200,000 replicates. The law horizons are shorter than
the timing horizons where 200,000 replicates at the timing horizon would not
fit the budget: at the timing horizons the reference would cost about
19,000 CPU seconds (W3 1,400, W4 2,700, W5 4,400, W6 8,600, W7 900). With the horizons
above the estimated cost is about 3,300 CPU seconds of simulation (W6 860,
W4 570, W7 440, vilar-noise 440, W5 440, W3 270, W2 155, repressilator 100,
the rest under 30), plus trajectory output and moment accumulation, which fits
the budget of about 90 CPU minutes on this two-CPU host. Every horizon still
contains every reaction firing many times per path: W5 and W6 at `end = 2`
have about 2,000 and 1,400 events per path; W3 covers the first oscillations
of the repressilator, the example `repressilator` its first 10 time units.
The example horizons are those of their run files (`end = 10`); vilar-noise
has no run file and uses 10 hours, which covers the first activator peak.

### Seeds

- 0.6.0 reference: chunk k (0-based, `replicates / chunk` chunks) uses root
  seed `1000000 + k`, the same for every model.
- Engine under test (law gate): chunk k uses root seed `2000000 + k`, with the
  same chunk size. A law-gate run of any method or selector uses these seeds;
  a failing run is not repeated with other seeds.
- Timing: root seed 4242 (as in `validation/compare_julia_accelerated.py`).

## Reference method and statistics

The reference law is that of 0.6.0 `direct` (checked policy: ChaCha8 words,
inverse waits, `fast: false`). It is computed with `cached_direct` under the
same policy, whose paths are bit-identical to `direct` (alternative-ssa.md:
"The new Direct/cached inverse variants retain legacy paths") and which costs
up to 4.5 times less on W5 and W6. `references.py` verifies the identity for
every model: the first `min(chunk, 1000)` replicates of chunk 0 are rerun
with `direct` and their saved states and per-replicate event counts must be
equal; otherwise the script stops.

Each chunk is an `ssa-run/0.3` request with `keep_trajectories = chunk`, so
every replicate's state at every report time is available. Per cell (species,
report time) the script accumulates the exact integer power sums
S1 … S4 of the n = 200,000 counts and derives, with exact rational arithmetic
rounded once to `f64`:

- mean x̄ = S1/n and its standard error se = sqrt(v/n);
- unbiased variance v = (S2 − S1²/n)/(n − 1);
- fourth central moment m4 = Σ(x − x̄)⁴/n;
- the standard error of the variance (fourth-moment estimate)
  sev = sqrt(max(0, (m4 − v²·(n − 3)/(n − 1))/n)).

The engine under test must compute the same four quantities with the same
formulas from its own n replicates.

## FSP references

For W1, W4, enzyme-mechanism, immigration-death and michaelis-menten:
`crn-craft fsp` (`run/0.1`, method `fsp`) with the bounds of the table,
`max_states` as in `models.json`, tolerance 1e-12, at the law report times.
The full FSP result (states, probabilities, error terms, raw moments) is
stored. The law moments (mean, variance, fourth central moment) are computed
from the retained probabilities normalized by their sum. A stored FSP
reference requires leakage + mass defect ≤ 1e-9 at every report time;
otherwise the bounds were wrong and the task stops. The other seven models
exceed the 100,000-state limit.

**Self-check of the references.** The 0.6.0 SSA reference is compared with
the FSP law by the FSP comparisons of the law gate below (n = 200,000):
89 cells, 178 comparisons, 0 failures expected. A failure stops the task.

## Kinetics and analysis references

- `crn-craft analyze MODEL` output (bytes and exit code) for every model.
- `kinetics::Kinetics` (0.6.0) right-hand side, state Jacobian and parameter
  Jacobian at the initial state and at a probe state (count of species i =
  initial count + i + 1), with the model's parameter values, for every model
  `Kinetics::compile` accepts; refusals are stored with their diagnostics.
  Expected refusals: W4 (Michaelis–Menten law consuming more than its
  substrate), W5 (256 species), W6 (992 reactions), W7 (96 species).

## Law gate (M2, M4a, M4b, M4c, after M8)

For every model of the run, the engine under test simulates the model's
`replicates` at its `end` and `report_times` with the gate seeds and computes
x̄, v, se, sev per cell as above. Per cell:

1. mean against the reference: |x̄ − x̄_ref| ≤ 5·sqrt(se² + se_ref²);
2. variance against the reference: |v − v_ref| ≤ 5·sqrt(sev² + sev_ref²);
3. where an FSP reference exists, mean against the FSP law (mean μ, variance
   σ², fourth central moment μ4): |x̄ − μ| ≤ 5·sqrt(σ²/n);
4. and variance against the FSP law:
   |v − σ²| ≤ 5·sqrt(max(0, μ4 − σ⁴·(n − 3)/(n − 1))/n).

A bound of zero demands equality (report time 0, where every count is the
initial count). The FSP comparisons use the exact law's standard errors, so a
rare event that one sample does not contain cannot fail them spuriously.

Size of a full run (all twelve models): 2,318 cells, 4,636 comparisons
against the reference and 178 against FSP, N = 4,814. Of these, 940 are at
time 0 (deterministic). Each of the 3,874 stochastic comparisons fails
under a correct engine with probability about 5.7e-7 (two-sided normal tail at
5 standard errors), so about 0.002 failures are expected per full run:
**0 failures are expected**. Any failure is investigated and reported; it is
not rerun with another seed, and the factor 5 is not changed. Runs over a
subset of models (M2: the examples, W1, W2, W3, W6; M4a: W5, W7; M4b: W4)
apply the same rules to their cells.

## Speed gate (M4d)

Workloads W1–W7 at the sizes of `validation/acceleration_workloads.py`
(`end` 100 with 64 replicates for W1–W4, `end` 20 with 16 replicates for
W5–W7), and the example models at `end = 10` with 10,000 replicates
(enzyme-mechanism, immigration-death, michaelis-menten), 1,000
(repressilator) and 256 (vilar-noise); save times 0 and `end`, one thread,
root seed 4242, no kept trajectories, statistics mean and variance.

- **0.6.0 methods** (21): `direct`, `cached_direct`, `direct_cr`,
  `sorting_direct`, `rssa_cr`, each with `fast` false and true and with
  inverse and Ziggurat waits (`ssa-run/0.3`), and `next_reaction`
  (`ssa-run/0.4`, inverse waits, `fast: false`, its only policy).
- **Boundary (gated):** wall time of one `crn-craft ssa MODEL REQUEST`
  process, from start to exit, stdout written to a file. One untimed warm-up,
  then three timed runs; the row's time is their median. All runs of one row
  must give byte-identical output.
- **Boundary (informative):** the Python prepared call
  (`crn_craft.prepare_ssa(...)`, then the native `simulate_json`), timed
  without preparation and without JSON decoding, also median of three after a
  warm-up.
- **0.6.0 best** per model: the smallest median over the 21 methods.
- **Gate:** for each of the twelve models, the unified engine's best exact
  method (smallest median over its methods and selectors, same sizes, same
  thread count, same boundary, its own run command) takes at most
  **1.25 ×** the 0.6.0 best. Timings on a shared virtual machine drift
  between sessions, so M4d reruns `baseline.py` with the kept 0.6.0 build in
  the same session, under the CPU lock, and compares against that rerun; the
  rerun's non-timing fields must equal `baseline-0.6.0.json`. All timings run
  under `flock /home/claude/craft5/cpu.lock`.

## Selector rule (M4c)

Composition–rejection is ported. From `baseline-0.6.0.json` (gated CLI
boundary), for each workload W1–W7 let t_cr, t_sort and t_rssa be the
smallest median of `direct_cr`, `sorting_direct` and `rssa_cr` over their four
variants (`fast` × waits). Sorting is ported if t_cr / t_sort ≥ 1.3 on at least
one workload; RSSA-CR is ported if t_cr / t_rssa ≥ 1.3 on at least one
workload. Otherwise the method is dropped and the M4c log says so. The example
models do not count. The committed file decides; M4c does not remeasure.

## M0c gate

- `baseline-0.6.0.json` and `reference/` exist and are committed; under 20 MB
  together.
- Determinism: `references.py` rerun on immigration-death,
  michaelis-menten, enzyme-mechanism and W1 (declared sizes) gives
  byte-identical files, normally and under `python -O`. `baseline.py` rerun
  on W1, W7 and immigration-death gives identical non-timing fields, normally
  and under `python -O`.
- Reference self-check: 0 failures among the 178 FSP comparisons.
- ruff format and ruff check pass on the scripts.

## Results

### M0c (2026-11-11, 0.6.0 engine, this host)

Host: Intel Xeon @ 2.10 GHz, 2 logical CPUs, Linux 6.18, 7.8 GiB. Every
run under `flock /home/claude/craft5/cpu.lock`, with the kept 0.6.0 build.

**References** (`references.py`, 2 worker threads): 39 min 57 s wall,
74.5 CPU minutes, 1.74e10 events, 52 files (513 KB).

| Check | Declared | Result |
| --- | --- | --- |
| `cached_direct` = `direct` paths (1,000 replicates per model) | identical | identical for 12 of 12 models |
| FSP leakage + mass defect | ≤ 1e-9 | largest 6.1e-13 |
| Self-check, SSA reference against FSP law | 178 comparisons, 0 failures | 178 comparisons, **0 failures**, largest \|z\| 2.08 |
| Kinetics refusals | W4, W5, W6, W7 | W4, W5, W6, W7 (`CRN-E080`) |
| `crn-craft analyze` | every model | 12 reports, exit 0 |

Self-check per model (comparisons, failures, largest |z|): W1 12, 0, 2.01;
W4 96, 0, 1.67; enzyme-mechanism 40, 0, 2.08; immigration-death 10, 0, 1.49;
michaelis-menten 20, 0, 1.08.

**Baseline** (`baseline.py`): 18 min 21 s wall, 252 rows (12 models × 21
methods). Best 0.6.0 method per model, gated CLI boundary, median of three
(the full table, including the best checked method, is in the M0c log):

| Model | Best method | Seconds | Events |
| --- | --- | ---: | ---: |
| W1 | `cached_direct/ziggurat/fast` | 0.0478 | 1,272,658 |
| W2 | `direct_cr/ziggurat/fast` | 0.0942 | 1,235,619 |
| W3 | `cached_direct/ziggurat/fast` | 0.2224 | 3,656,543 |
| W4 | `rssa_cr/ziggurat/fast` | 0.1062 | 1,227,674 |
| W5 | `direct_cr/ziggurat/fast` | 0.0677 | 325,267 |
| W6 | `rssa_cr/ziggurat/fast` | 0.0630 | 229,211 |
| W7 | `rssa_cr/ziggurat/fast` | 0.0209 | 108,058 |
| enzyme-mechanism | `cached_direct/ziggurat/fast` | 0.0765 | 460,645 |
| immigration-death | `cached_direct/ziggurat/fast` | 0.0724 | 626,862 |
| michaelis-menten | `direct/ziggurat/fast` | 0.0612 | 307,345 |
| repressilator | `sorting_direct/ziggurat/fast` | 0.2867 | 4,194,928 |
| vilar-noise | `cached_direct/ziggurat/fast` | 1.0809 | 14,224,888 |

**Selector rule (M4c).** t_cr / t_sort per workload W1–W7: 1.158, 0.746,
**1.315**, 1.090, 0.367, 0.598, 0.890. t_cr / t_rssa: 0.856, 0.822, 0.874,
**8.030**, 0.938, **3.351**, 1.071. Sorting reaches 1.3 on W3 and RSSA-CR on
W4 and W6: **both are ported**. The sorting margin is small (1.315 against
1.3); at the informative Python boundary it is 1.384 on W3.

**M0c gate.**

| Check | Result |
| --- | --- |
| Outputs committed, total size | 1.05 MB |
| `references.py` rerun (immigration-death, michaelis-menten, enzyme-mechanism, W1), normal | 24 of 24 files byte-identical |
| The same under `python -O` | 24 of 24 files byte-identical |
| `baseline.py` rerun (W1, W7, immigration-death), non-timing fields, normal | 63 of 63 rows identical (host and tools too) |
| The same under `python -O` | 63 of 63 rows identical |
| ruff 0.11.13 format and check on `validation/merge` | pass |
| Workspace checks (fmt, clippy, release tests incl. I1/I3, rustdoc, manual I2 normal and `-O`) | pass |

### M2 (converted models, `static_direct`): declaration

Written and committed before the gate run.

- Engine under test: branch `merge-0.7` with the M2 converter, through
  Python (the wheel built from the checkout) with 2 worker threads:
  `check_laws.py --python ... --method '{"kind": "static_direct"}'
  --expect-refused W4,W5,W7`. The method runs in `crn-craft.rules-run/0.6`,
  so with Ziggurat waits (sampler `chacha8-u53-ziggurat/rules-static-direct/1`).
- Models: all twelve rows of `models.json`. The converter must refuse W4
  (order 8), W5 (256 species) and W7 (96 species) with `CRN-R016`, and
  convert the other nine.
- Size: the nine converted models have 510 cells, 1,020 comparisons against
  the reference and 82 against FSP laws (W1, enzyme-mechanism,
  immigration-death, michaelis-menten): N = 1,102, of which 204 are at t = 0.
  At about 5.7e-7 per stochastic comparison, about 5e-4 failures are
  expected: **0 failures are expected**. A failure is investigated once for a
  conversion error, not rerun with other seeds.
- Output: `validation/merge/laws/static_direct.json` (every row, failed rows
  included), summarized below.
