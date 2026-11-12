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

### M2 (converted models, `static_direct`): results

2026-11-11, this host, under the CPU lock, engine at `d5fb778` (wheel
`_native` SHA-256 `77d90cd7…`), 2 threads: 18 min 36 s wall, 34.7 CPU
minutes, 1.56e10 events. Report: `validation/merge/laws/static_direct.json`.

| Model | Cells | Comparisons (reference) | Failures | max \|z\| | Comparisons (FSP) | Failures | max \|z\| |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 | 6 | 12 | 0 | 2.45 | 12 | 0 | 2.28 |
| W2 | 198 | 396 | 0 | 3.39 | – | – | – |
| W3 | 36 | 72 | 0 | 1.88 | – | – | – |
| W4 | refused (`CRN-R016`, order 8), as declared | | | | | | |
| W5 | refused (`CRN-R016`, 256 species), as declared | | | | | | |
| W6 | 160 | 320 | 0 | 3.33 | – | – | – |
| W7 | refused (`CRN-R016`, 96 species), as declared | | | | | | |
| enzyme-mechanism | 20 | 40 | 0 | 1.78 | 40 | 0 | 2.16 |
| immigration-death | 5 | 10 | 0 | 1.55 | 10 | 0 | 1.81 |
| michaelis-menten | 10 | 20 | 0 | 1.90 | 20 | 0 | 2.20 |
| repressilator | 30 | 60 | 0 | 1.88 | – | – | – |
| vilar-noise | 45 | 90 | 0 | 2.94 | – | – | – |
| **Total** | 510 | 1,102 (N), 204 at t = 0 | **0** | 3.39 | | | |

**Gate passed:** 0 failures among 1,102 comparisons, as declared; the
refusals are exactly W4, W5 and W7. W3 and repressilator convert to the same
bytes and run with the same seeds, so their cells at t = 0, 2, 5 and 10 are
the same comparisons (the 0.6.0 references share seeds in the same way): 48
of the 1,102 comparisons repeat others, 36 of them stochastic.

### M4a (W5 and W7, `static_direct`): declaration

Written and committed before the gate run.

- Engine under test: branch `merge-0.7` with the M4a limits (wheel built
  from the checkout), through Python with 2 worker threads:
  `check_laws.py --python ... --method '{"kind": "static_direct"}'
  --models W5,W7 --label m4a-static_direct`. The converted models are
  `crn-craft.model/1.0`, so the requests are `crn-craft.run/1.0` (Ziggurat
  waits, sampler `chacha8-u53-ziggurat/rules-static-direct/1`).
- Models: W5 (256 species, 510 reactions) and W7 (96 species, 126
  reactions), the rows of `models.json`; both must convert (no expected
  refusal).
- Size: W5 1,280 cells, W7 480 cells: 1,760 cells, 3,520 comparisons
  against the reference (no FSP law), of which 704 at t = 0. At about
  5.7e-7 per stochastic comparison, about 1.6e-3 failures are expected:
  **0 failures are expected**. A failure is investigated once, not rerun
  with other seeds.
- Output: `validation/merge/laws/m4a-static_direct.json` (every row, failed
  rows included). The wall time per model is printed on stderr and recorded
  in the M4a log; it is not gated here (the speed gate is M4d).

### M4a (W5 and W7, `static_direct`): results

2026-11-11, this host, under the CPU lock, engine at `66560f4` (wheel
`_native` SHA-256 `32d62d00…`), `crn-craft.run/1.0`, sampler
`chacha8-u53-ziggurat/rules-static-direct/1`, 2 threads: 5 min 43 s wall
for both models. Report: `validation/merge/laws/m4a-static_direct.json`.

| Model | Cells | Comparisons (reference) | Failures | max \|z\| | Events (gate) | Events (0.6.0 reference) | Wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W5 | 1,280 | 2,560 | 0 | 3.88 | 407,994,800 | 408,007,617 | 221.2 s |
| W7 | 480 | 960 | 0 | 2.91 | 676,470,648 | 676,431,647 | 121.9 s |
| **Total** | 1,760 | 3,520 (N), 704 at t = 0 | **0** | 3.88 | | | |

**Gate passed:** 0 failures among 3,520 comparisons, as declared; both
models convert. The wall time covers conversion, simulation with every path
kept and the decoding in Python (2 worker threads); it is informative only.

### M4b (higher-order reactions): declaration

Written and committed before the FSP references were computed and before
any run of the checks below.

**Models** (`validation/merge/higher-order/models.json`, network models
`crn-craft.model/0.1`, converted by the engine under test):

| Model | Reaction | Initial | c | Report times | FSP bounds, max states |
| --- | --- | --- | ---: | --- | --- |
| third-order-decay | 3A → ∅, mass action | A = 40 | 0.002 | 0, 0.5, 1, 2, 5, 10 | A ≤ 40, 41 |
| trimolecular | 2A + B → C, mass action | A = 30, B = 20, C = 0 | 0.002 | 0, 0.25, 0.5, 1, 2, 5 | A ≤ 30, B ≤ 20, C ≤ 20, 100 |

The parameters were chosen with an independent dense-matrix solution of the
two chains (not the engine under test) so that every report time after 0
has between 5 and 9 states of probability at least 1e-3.

**References.** `higher_order_references.py` runs the kept 0.6.0 CLI
(`/home/claude/m07/bin060/crn-craft fsp`, `run/0.1`, tolerance 1e-12) and
stores the full FSP results and their moments under
`validation/merge/reference/higher-order/`, with the M0c rule: leakage +
mass defect ≤ 1e-9 at every report time, otherwise the task stops.

**Exact-law check** (`check_higher_order.py`, Python with the wheel under
test, 2 threads): each model with `static_direct` and `next_reaction`
(`crn-craft.run/1.0`), n = 200,000 replicates per model and method in
chunks of 50,000, chunk k with root seed 3,000,000 + k. Comparisons per
model, method and report time:

1. mean and variance of each species against the FSP law (formulas 3 and
   4 of the law gate, factor 5);
2. each joint state with FSP probability p ≥ 1e-3 (normalized by the
   retained mass): |p̂ − p| ≤ 5·sqrt(p(1 − p)/n);
3. the remaining states together, with the same bound.

The independent solution gives 53 comparisons per method for
third-order-decay (12 moments, 35 states, 6 remainders) and 83 for
trimolecular (36, 41, 6): N = 272 for both methods, of which 24 are at
t = 0 and demand equality. The count follows from the stored references;
the M4b log reports it. About 1.5e-4 failures are expected under a correct
engine: **0 failures are expected**. Before M4b is implemented the
converter refuses both models (order 3), so the check fails first; that run
is recorded below as the failing-first run, not as a gate run. A failure
after the implementation is investigated once, not rerun with other seeds.

**Law gate for W4** (`check_laws.py --python ... --method '{"kind":
"static_direct"}' --models W4 --label m4b-static_direct`, gate seeds of
`models.json`): 48 cells, 96 comparisons against the reference and 96
against the FSP law, N = 192, of which 32 at t = 0; **0 failures are
expected**.

**Regression run after M4a and M4b** (`check_laws.py --python ... --method
'{"kind": "static_direct"}' --label m4-static_direct`, all twelve models,
no `--expect-refused`): every model must convert; N = 4,814 comparisons
(4,636 against the references, 178 against FSP laws), 940 at t = 0;
**0 failures are expected**. Wall time per model is printed on stderr and
recorded in the M4b log (not gated; the speed gate is M4d).

### M4b (higher-order reactions): results

**References** (2026-11-11, kept 0.6.0 CLI, `higher_order_references.py`):
third-order-decay 14 reachable states, trimolecular 16; leakage + mass
defect at most 5.0e-13 and 5.5e-13 (limit 1e-9). Rerun under `python -O`
into another directory: all five files byte-identical. The means agree with
the independent dense-matrix solution used to choose the parameters
(third-order-decay at t = 0.5: 25.248 both).

**Failing-first run** (engine at `39aecef`, M4a wheel, before the M4b
implementation): both models refused by the converter (`CRN-R016`, order
3), exit 1, 0 comparisons. Report:
`validation/merge/laws/higher-order-before-m4b.json`. Not a gate run.

**Exact-law check** (engine at `1f66966`, wheel `_native` SHA-256
`85673a83…`, `crn-craft.run/1.0`, 2 threads): 12 s wall. Report:
`validation/merge/laws/higher-order.json`; the run under `python -O` wrote
a byte-identical report. The converted models hold
`{"label": "x", "type": "A", "count": 3}` and
`{"label": "x", "type": "A", "count": 2}, {"label": "y", "type": "B"}`.

| Model | Method | Sampler | Comparisons | at t = 0 | Failures | max \|z\| |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| third-order-decay | `static_direct` | `chacha8-u53-ziggurat/rules-static-direct/1` | 53 | 4 | 0 | 2.83 |
| third-order-decay | `next_reaction` | `chacha8-u53/rules-next-reaction/1` | 53 | 4 | 0 | 2.56 |
| trimolecular | `static_direct` | `chacha8-u53-ziggurat/rules-static-direct/1` | 83 | 8 | 0 | 1.77 |
| trimolecular | `next_reaction` | `chacha8-u53/rules-next-reaction/1` | 83 | 8 | 0 | 2.48 |
| **Total** | | | 272 (N) | 24 | **0** | 2.83 |

**Passed:** N = 272 as declared (53 + 83 per method), 0 failures.

**Law gate for W4** (engine at `1f66966`, wheel `85673a83…`, `run/1.0`,
2 threads). A first run under the lock was stopped by a host restart at
chunk 17 of 20, before any comparison was written; it was rerun with the
declared seeds. Report: `validation/merge/laws/m4b-static_direct.json`.

| Model | Cells | Comparisons (reference) | Failures | max \|z\| | Comparisons (FSP) | Failures | max \|z\| | Events | Wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W4 | 48 | 96 | 0 | 1.94 | 96 | 0 | 1.72 | 742,510,499 (0.6.0 reference 742,514,249) | 558.1 s |

**Passed:** N = 192, 32 at t = 0, 0 failures, as declared. The converted
W4 has 66 rules: 33 sources and 33 higher-order rules of eight count-1
patterns with a whole (Michaelis–Menten) propensity, so this gate exercises
the whole-propensity path of `Multi`; the combination count is exercised by
the exact-law check above.

### M4a and M4b: regression run over all twelve models

Engine at `1f66966`, wheel `85673a83…`, `{"kind": "static_direct"}` in
`run/1.0` (Ziggurat waits), 2 threads, no expected refusal: 32 min 6 s wall.
Report: `validation/merge/laws/m4-static_direct.json`.

| Model | Cells | Comparisons (reference) | Failures | max \|z\| | Comparisons (FSP) | Failures | max \|z\| | Wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 | 6 | 12 | 0 | 2.45 | 12 | 0 | 2.28 | 13.2 s |
| W2 | 198 | 396 | 0 | 3.39 | – | – | – | 67.0 s |
| W3 | 36 | 72 | 0 | 1.88 | – | – | – | 158.4 s |
| W4 | 48 | 96 | 0 | 1.94 | 96 | 0 | 1.72 | 532.8 s |
| W5 | 1,280 | 2,560 | 0 | 3.88 | – | – | – | 212.5 s |
| W6 | 160 | 320 | 0 | 3.33 | – | – | – | 253.1 s |
| W7 | 480 | 960 | 0 | 2.91 | – | – | – | 112.0 s |
| enzyme-mechanism | 20 | 40 | 0 | 1.78 | 40 | 0 | 2.16 | 2.9 s |
| immigration-death | 5 | 10 | 0 | 1.55 | 10 | 0 | 1.81 | 2.2 s |
| michaelis-menten | 10 | 20 | 0 | 1.90 | 20 | 0 | 2.20 | 2.3 s |
| repressilator | 30 | 60 | 0 | 1.88 | – | – | – | 58.3 s |
| vilar-noise | 45 | 90 | 0 | 2.94 | – | – | – | 510.2 s |
| **Total** | 2,318 | 4,814 (N), 940 at t = 0 | **0** | 3.88 | | | | |

**Passed:** N = 4,814 as declared, 0 failures, every model converts. The
power sums and event counts of every model are identical to the earlier
gate runs (M2 for the nine M2 models, M4a for W5 and W7, the W4 gate
above): the 1.0 limits, the index-based lowering and the higher-order code
changed no number of these runs. Wall times include conversion, simulation
with every path kept and decoding in Python; they are informative (the
speed gate is M4d).

### M4e (quantile statistics): declaration

Written and committed before the gate run. Script:
`validation/merge/check_quantiles.py --binary target/release/crn-craft`.

- Model: `models/rules/immigration-death.rules.json` (k_in = 10, k_out = 1,
  X(0) = 0) with only its schema string rewritten to
  `crn-craft.model/1.0`. Exact law at time t: Poisson with mean
  λ(t) = 10 (1 − e^−t); at the late report time t = 10,
  λ = 9.999546.
- Request: `crn-craft.run/1.0`, `end` 10, report times 0 and 10,
  100,000 replicates, root seed 20261111, 2 threads, statistics mean and
  variance, `outputs.quantiles` [0.025, 0.5, 0.975], no kept paths.
- Methods: `static_direct` (run/1.0 defaults), `tree`, `next_reaction`.
- Poisson quantile: Q(p) = the smallest integer k with F(k) ≥ p, F the
  Poisson distribution function of mean λ(10). Q(0.025) = 4,
  Q(0.5) = 10, Q(0.975) = 17. Margins: F(3) = 0.01034, F(4) = 0.02926;
  F(9) = 0.45799, F(10) = 0.58310; F(16) = 0.97297, F(17) = 0.98573. The
  smallest margin, 0.975 − F(16) = 0.00203, is 4.1 standard errors of the
  empirical distribution function at n = 100,000 (0.00049).
- **Gate:** for each method and p, the reported quantile at t = 10 lies
  within one count of Q(p): |q̂ − Q(p)| ≤ 1. 9 comparisons. A failure
  needs the empirical distribution function to miss F by at least 0.0178
  at Q(p) − 2 or Q(p) + 1 (closest: p = 0.975, F(18) − 0.975 = 0.01782;
  then p = 0.025, 0.025 − F(2) = 0.02223), at least 36 standard errors:
  **0 failures are expected**. A failure is investigated, not rerun with
  another seed. At t = 0 every quantile must be exactly 0.
  (Corrected before the gate run: the first committed text gave this
  margin as 0.024 and 48 standard errors, from F(15) at p = 0.975 only.)
- Informative, not gated: whether q̂ equals Q(p) exactly; the mean
  and variance at t = 10 against λ (|x̄ − λ| in standard errors).
- Contract checks in the same script (exact, no tolerance): the request
  validates against `schemas/run-1.0.schema.json` and each result against
  `schemas/result-1.0.schema.json`; the result lists the quantiles in
  request order with values shaped like `statistics.mean`; the same request
  as `crn-craft.rules-run/0.6` is refused with `CRN-R020` at `/schema`.

Measurement (no target): the added time and memory of quantiles 0.025,
0.5 and 0.975 for 4,000 cells × 37 report times × 6 observables
(`models/rules/mtdna-three-class-cme.rules.json`, `static_direct`, the
37-point grid of `check_compiled_loop.py`, 2 threads), by
`validation/merge/measure_quantiles.py`: wall time of `simulate()`, median
of five, and peak resident memory of a fresh process, with and without
quantiles, under the CPU lock.

### M4e (quantile statistics): results

2026-11-11, this host, under the CPU lock, engine as committed in `bea803a`
(CLI `target/release/crn-craft` of the task worktree), 2 threads. Report:
`validation/merge/laws/quantiles-immigration-death.json`.

| Method | p | q̂ at t = 10 | Q(p) | \|q̂ − Q\| | Passed |
| --- | ---: | ---: | ---: | ---: | --- |
| `static_direct` | 0.025 | 4 | 4 | 0 | yes |
| `static_direct` | 0.5 | 10 | 10 | 0 | yes |
| `static_direct` | 0.975 | 17 | 17 | 0 | yes |
| `tree` | 0.025 | 4 | 4 | 0 | yes |
| `tree` | 0.5 | 10 | 10 | 0 | yes |
| `tree` | 0.975 | 17 | 17 | 0 | yes |
| `next_reaction` | 0.025 | 4 | 4 | 0 | yes |
| `next_reaction` | 0.5 | 10 | 10 | 0 | yes |
| `next_reaction` | 0.975 | 17 | 17 | 0 | yes |

**Gate passed:** 9 comparisons, 0 failures, every quantile equal to the
Poisson quantile; every quantile at t = 0 is 0. Informative: the means at
t = 10 lie 0.73, −1.46 and −1.23 standard errors from λ (static_direct,
tree, next_reaction). Contract checks: requests and results valid against the
committed 1.0 schemas, quantiles in request order shaped like the mean, and
the `rules-run/0.6` refusal (`CRN-R020` at `/schema`).

Run history, failed run kept: the first started run was cut off by a server
restart before it wrote anything. The first complete run gave the same 9
rows (byte-identical) and failed one contract check, the `rules-run/0.6`
refusal: the engine refused correctly (exit 1, `CRN-R020` at `/schema`), but
the script read the CLI's output as a bare list of diagnostics, while the CLI
prints `{"success": false, "diagnostics": [...]}`. The parser was fixed (no
tolerance or row changed) and the script rerun: all checks pass, and the
report is byte-identical under `python -O`.

**Measurement** (`validation/merge/measure_quantiles.py`, report
`validation/merge/quantile-cost.json`): 4,000 cells of the three-class CME
model, 37 report times, 6 observables, `static_direct`, 2 threads, about
1.30e9 events per ensemble; three fresh processes per configuration, five
timed `simulate()` calls each.

| Configuration | Median `simulate()` | Range | Peak RSS (three processes) |
| --- | ---: | ---: | --- |
| no quantiles | 66.21 s | 63.62–68.60 s | 88,808 / 88,896 / 88,812 KiB |
| quantiles 0.025, 0.5, 0.975 | 64.31 s | 61.77–66.69 s | 95,756 / 95,992 / 96,048 KiB |

Added memory: 7,180 KiB (median of the peaks), against a buffer of
4,000 × 37 × 6 × 8 = 7,104,000 bytes (6,938 KiB) plus a 32,000-byte sort
column. Added time: −1.90 s, inside the run-to-run spread of about ±3 s; the
quantile work (222 columns of 4,000 values sorted once, three evaluations
each) is not measurable against 65 s of simulation.

### M5a (finite state projection on the flat view): declaration

Written and committed before the first run of the checks below. The engine
under test solves `crn-craft.run/1.0` requests with method `fsp` through
`crn-craft rules MODEL REQUEST` (`crn_craft::fsp::v1`), on models converted
by the engine under test (`crn-craft convert`).

**Reference gate** (`validation/merge/check_fsp_references.py --binary
target/release/crn-craft`, standard library only). For each stored 0.6.0
FSP reference (W1, W4, enzyme-mechanism, immigration-death and
michaelis-menten from M0c; third-order-decay and trimolecular from M4b):
convert the network model, request the same bounds (by species name),
`max_states` and tolerance (1e-12) as the reference, with `end` and the
report times of the reference's save grid, and compare by species name and
state:

1. the projected state set equals the reference's (7 checks);
2. every probability: |p − p_ref| ≤ 1e-10 absolute (11,232 comparisons:
   W1 251 × 6, W4 1,001 × 6, enzyme-mechanism 506 × 5, immigration-death
   151 × 5, michaelis-menten 51 × 5, third-order-decay 14 × 6,
   trimolecular 16 × 6 states × report times);
3. the sink probability at every report time: |leak − leak_ref| ≤ 1e-10
   (38 comparisons).

N = 11,277 checks; **0 failures are expected**. Rates are evaluated by
expressions, so bits may differ (Hill laws, mass action over pairs); all
seven models are mass action or Michaelis–Menten. The report lists the
largest difference per model, and as information (not gated) the
uniformization rate, the edge count and the largest moment difference. The
report must be identical normally and under `python -O`.

**Migrated `check_fsp.py`** (`validation/check_fsp_rules.py --binary
target/release/crn-craft --dsmts DSMTS_CHECKOUT`): every case of
`validation/check_fsp.py` (the four CME fixtures F1–F4 of `check_cme.py`
and DSMTS 001-01 and 002-01) on the converted model through the new
command, with the old script's tolerances: distribution error against an
independent matrix exponential ≤ 2e-12 at every report time, series tail
bound ≤ 1.01e-14, mass balance |Σp + discarded − 1| < 2e-12, DSMTS
analytic mean error < 1e-9 and variance error < 1e-7, and the upstream
rounded DSMTS means < 5.1e-6 from the analytic law. `check_fsp.py` stays
on the network `crn-craft fsp` command until the network engine is
removed (M8). The report must be identical normally and under `python -O`.

### M5a (finite state projection on the flat view): results

Engine of the M5a commit, release build, this host, under the CPU lock.

**Reference gate** (`validation/merge/check_fsp_references.py`, report
`validation/merge/fsp/references.json`, identical under `python -O`):
**0 failures** in 11,278 checks. The declaration counted 38 sink
comparisons and N = 11,277; the seven references have 6 + 6 + 5 + 5 + 5 +
6 + 6 = 39 report times, so the run made 39 sink comparisons and N =
11,278. The checks themselves are as declared.

| Model | States × report times | State set | Largest probability difference | Largest sink difference | λ, edges (new = reference) |
| --- | --- | --- | ---: | ---: | --- |
| W1 | 251 × 6 | equal | 0 | 0 | 350, 501 |
| W4 | 1,001 × 6 | equal | 0 | 0 | 196.0588…, 2,001 |
| enzyme-mechanism | 506 × 5 | equal | 1.4e-17 | 0 | 6, 1,365 |
| immigration-death | 151 × 5 | equal | 0 | 0 | 34, 301 |
| michaelis-menten | 51 × 5 | equal | 0 | 0 | 3.5714…, 50 |
| third-order-decay | 14 × 6 | equal | 0 | 0 | 19.76, 13 |
| trimolecular | 16 × 6 | equal | 0 | 0 | 17.4, 15 |

Six of seven models reproduce the 0.6.0 probabilities bit for bit; in
enzyme-mechanism the mass action over pairs (E + S) differs in the last
bits (the M2 tolerance of 4 ε relative), largest raw-moment difference
2.3e-13.

**Migrated `check_fsp.py`** (`validation/check_fsp_rules.py` with the
pinned DSMTS checkout `e79fd58`): 16 of 16 checks pass, normally and under
`python -O` with identical reports. Largest distribution error against the
independent matrix exponential 4.49e-14 (F4), the same as `check_fsp.py`
on the network command; every row equals the old script's row except F3
(Hill and Michaelis–Menten), 5.83e-16 against 6.11e-16. `check_fsp.py` on
the network `crn-craft fsp` command of this build gives a report
byte-identical to the 0.6.0 binary's (it includes the SHA-256 of every
result), so the old command's output is unchanged.

### M5b (network analysis on the flat view): declaration

Written and committed before the first run of the checks below. The engine
under test analyzes a model with `crn-craft analyze MODEL`, which reads a
network `crn-craft.model/0.1` through the converter and a rules model or
`crn-craft.model/1.0` directly, and reports `crn-craft.analysis/1.0`.

**Reference gate** (`validation/merge/check_analysis_references.py --binary
target/release/crn-craft`, standard library only). For each of the twelve
models of `models.json`, two runs of `crn-craft analyze`:

1. on the network model file (the input of the 0.6.0 reference);
2. on the model converted by the engine under test (`crn-craft convert`).

Each parsed report is compared with the stored 0.6.0 report
(`reference/<name>.analysis.json.gz`) leaf by leaf: every object key, every
list length and every value. Ignored: `schema` (0.6.0
`crn-craft.analysis/0.2`, now `crn-craft.analysis/1.0`, checked to be that
string), and in run 2 `model_sha256`, which is the hash of the converted
bytes there (checked to be that hash). The exit code must be 0 as in the
reference. Values must be equal, with one declared exception: an element of
`propensities` may differ when the 0.6.0 channel is a Hill law (relative
difference ≤ 1e-14 · max(1, n)) or mass action over two molecules (relative
difference ≤ 4 ε), the rate tolerances M2 declared and measured for the
converted laws. Every such difference is listed in the report; a difference
outside it, or in any other leaf, is a failure. N is the number of leaves
compared, counted by the script. **0 failures are expected.** The report,
`validation/merge/analysis/references.json`, holds no timings and must be
identical normally and under `python -O`.

**Migrated `check_network.py` and `check_structure.py`**: the same random
network streams (seeds 20261101 and 20261103, unchanged fixture hashes),
example models and special cases, analyzed by `crn-craft analyze` of the
engine under test on the network files, with every structural and numerical
check and tolerance of the 0.6.0 scripts. Cases whose expected outcome
changes with the converter are counted and listed in the reports: a network
without reactions is refused (`CRN-R016`, a rules model needs a rule); a
propensity that overflows f64 fails with `CRN-R030` (0.6.0: `CRN-E042`); a
positive propensity that underflows to 0 is reported as 0 (0.6.0: refused
with `CRN-E042`); its tolerance check against the exact value applies as
before. Both reports must be identical normally and under `python -O`.

### M5f (GillesPy2 export on the flat view): declaration

Written and committed before the first run of the checks below. The engine
under test exports with `crn-craft export MODEL gillespy2`, which reads a
network `crn-craft.model/0.1` through the converter and a rules model or
`crn-craft.model/1.0` directly. A rule whose rate is a constant and whose
channel takes at most two molecules exports as GillesPy2 mass action, as in
0.6.0; any other channel exports as a custom propensity string.

**Export equivalence** (`validation/merge/check_export_references.py
--binary target/release/crn-craft --reference-binary BINARY_0_6_0`,
standard library only), over the twelve models of `models.json`, the five
examples and the fixture of `validation/check_export.py`:

1. every model the 0.6.0 binary exports successfully exports again, and its
   `artifact` (source, aliases, notes) is byte-identical to the 0.6.0
   artifact;
2. every model the 0.6.0 binary refused with `CRN-E050` (Hill and
   Michaelis–Menten laws, reactant order above two) now exports; each
   custom propensity string, evaluated with Python floats (the operators
   GillesPy2 1.8.3 keeps: `+ - * /` and `pow`), equals the engine's
   propensity (`crn-craft analyze MODEL STATE`) within the tolerance of
   `check_export.py` (relative 2e-12 plus absolute 8 × 2^-1074) at four
   states per model: the initial counts, all counts 0, all counts 7, and
   count i = initial count + i + 1 of species i;
3. no other export outcome changes (refusals keep their code).

N is the number of artifacts and propensities compared, counted by the
script. **0 failures are expected.** The report,
`validation/merge/export/references.json`, holds no timings and must be
identical normally and under `python -O`.

**Migrated `check_export.py` and `check_importers.py`**: every case and
tolerance of the 0.6.0 scripts. `check_export.py` keeps its fixture, its
GillesPy2 1.8.3 run and its rate comparisons; the two models it expected to
be refused (Michaelis–Menten, repressilator) now export custom propensities
and are compared like the others. The importers return the converted
`crn-craft.model/1.0` document; `check_importers.py` reads its rate
constants from the converted rules. Both reports must be identical
normally and under `python -O`.

### M5b (network analysis on the flat view): results

Engine of the M5b commit, release build, this host, under the CPU lock.

**Reference gate** (`validation/merge/check_analysis_references.py`, report
`validation/merge/analysis/references.json`, identical under `python -O`):
**0 failures** in 194,388 compared leaves (twelve models, each analyzed
from the network file and from its converted model). No propensity
differs: the declared exception for Hill laws and mass action over pairs
was not needed at the initial states.

| Model | Leaves (network file / converted) | Failures | Listed propensity differences |
| --- | ---: | ---: | ---: |
| W1 | 45 / 44 | 0 | 0 |
| W2 | 1,005 / 1,004 | 0 | 0 |
| W3 | 239 / 238 | 0 | 0 |
| W4 | 4,963 / 4,962 | 0 | 0 |
| W5 | 11,514 / 11,513 | 0 | 0 |
| W6 | 74,784 / 74,783 | 0 | 0 |
| W7 | 3,698 / 3,697 | 0 | 0 |
| enzyme-mechanism | 134 / 133 | 0 | 0 |
| immigration-death | 45 / 44 | 0 | 0 |
| michaelis-menten | 53 / 52 | 0 | 0 |
| repressilator | 239 / 238 | 0 | 0 |
| vilar-noise | 481 / 480 | 0 | 0 |

(The converted run has one leaf fewer: its `model_sha256` is checked
separately.)

**Migrated `check_network.py`** (report schema `network-validation/0.2`,
identical under `python -O`): 219 of 219 checks, the same fixture stream
(`e6a64d11…`) and the same 219 cases as the 0.6.0 script on the 0.6.0
binary. 1,030 rates evaluated (0.6.0: 1,028; the two underflow cases are now
evaluated), largest relative error 4.36e-15 (0.6.0: 2.37e-14). Changed
expectations, as declared: 21 random networks without reactions refused
with `CRN-R016`; `mass_2000_1000_1` fails with `CRN-R030`;
`hill_repression_1_5e-324` and `mm_underflow` report 0.

**Not declared before the run**: one case, `hill_activation_1_1e+300`
(x = 1, k = vmax = 1e300, n = 2, exact value 1e-300), fails the 0.6.0
tolerance: the converted law `basal + vmax / (1 + (k / R)^n)` gives 0 where
`(k / R)^n` overflows. This is the limitation M2 declared and measured for
the converted Hill laws (absolute floor vmax · 2^-1020 where the induced
term underflows; there the converted propensity is basal). The migrated
script accepts a Hill propensity that misses the 0.6.0 tolerance only when
it equals the basal rate exactly and lies within M2's floor, and lists each
such case (`hill_floor_cases`: this one case). No other tolerance changed.

**Migrated `check_structure.py`** (report schema
`structure-validation/0.2`, identical under `python -O`): 309 of 309
checks, the same fixture stream (`1373b12d…`) and cases as 0.6.0. 37 of the
300 random networks have no reactions and are refused with `CRN-R016` as
declared, so fewer rays and bounds are checked than in 0.6.0: 211 extreme
rays (0.6.0: 373), 349 basis-bounded and 363 cone-bounded species (0.6.0:
511 and 525); the cone adds 14 bounds over the basis in both.

### M5f (GillesPy2 export on the flat view): results

Engine of the M5f commits, release build, this host, under the CPU lock;
reference binary: the 0.6.0 build of M0c (`/home/claude/m07/bin060`).

**Export equivalence** (`validation/merge/check_export_references.py`,
report `validation/merge/export/references.json`, identical under `python
-O`): **0 failures**. The nine models 0.6.0 exported (W1, W2, W5, W6, W7,
enzyme-mechanism, immigration-death, vilar-noise and the `check_export.py`
fixture) export byte-identical artifacts. The four models 0.6.0 refused with
`CRN-E050` (W3, W4, michaelis-menten, repressilator) now export 40 custom
propensities; at the four declared states each, 160 evaluations, every one
equal to the engine's propensity bit for bit (largest relative difference
0).

**Migrated `check_export.py`** (GillesPy2 1.8.3, report schema
`reference-export-validation/0.2`, identical under `python -O`): 9 of 9
checks; 15,756 rate comparisons, largest relative error 1.35e-16, the grid
fixture exact as in 0.6.0. The Michaelis–Menten and repressilator examples
now export and pass the rate checks; a GillesPy2 SSA smoke run of the
Michaelis–Menten custom propensity conserves S + P; the two refusal rows
use a rate with a function and a stop condition (`CRN-E050`).

**Migrated `check_importers.py`**: 49 of 49 checks (the 48 of 0.6.0 and the
new rejection of an SBML model without reactions), identical under
`python -O`.

### M5e (ODE right-hand side and Jacobians on the flat view): declaration

Written and committed before the first run of the checks below. The engine
under test is `crn_craft::kinetics::Kinetics` of the M5e commits, which
reads every model through its flat view (`Kinetics::from_bytes`; a network
`crn-craft.model/0.1` is converted first) and evaluates the rate and
propensity expressions with dual numbers.

**Reference gate** (`crates/crn-craft/tests/kinetics_references.rs`, run
with `cargo test --release -p crn-craft --test kinetics_references`; with
`CRN_KINETICS_REPORT=PATH` it also writes its report, committed as
`validation/merge/kinetics/references.json`, without timings). For each of
the twelve models of `models.json` and its M0c reference
`reference/<name>.kinetics.json`:

1. the model file's SHA-256 equals the reference's `model_sha256`;
2. acceptance: the engine under test accepts the network model file if and
   only if 0.6.0 did; a refusal is `CRN-E080`, as the reference's (W4, W5,
   W6, W7 are expected to be refused and have no values);
3. for each accepted model, run 1 on the network model file and run 2 on
   the model converted by the engine under test (`crn-craft.model/1.0`
   bytes from `convert::from_network`): the species labels, parameter names
   and parameter values equal the reference's `species`, `parameters` and
   `parameter_values` exactly, and at each stored point (the initial and
   the probe state, the reference's `y` and parameter values) every entry
   of the right-hand side (n), the state Jacobian (n × n) and the parameter
   Jacobian (n × m) is compared with the reference value r:
   - r ≠ 0: |x − r| ≤ 1e-12 · |r| (relative);
   - **r = 0 (either sign): |x| ≤ 1e-12 · M**, where M is the largest
     absolute value of the same array (the right-hand side, the state
     Jacobian or the parameter Jacobian) at the same point of the same
     reference; when M = 0, x must be 0. Zero references are thus compared
     on the scale of their array, since a zero that comes from cancelling
     terms (W2's probe: 3 − 0.2 · 15) can be off by rounding. The report
     counts the zero references reproduced exactly (x = ±0) and lists every
     zero reference that is not;
   - an evaluation that fails where the reference has values is a failure.

N = 3,126 values per run (W1 4, W2 2,244, W3 168, enzyme-mechanism 64,
immigration-death 8, michaelis-menten 20, repressilator 168, vilar-noise
450), 6,252 in both runs, plus the twelve acceptance checks. **0 failures
are expected.** The report lists, per model and run, the largest relative
difference over the nonzero references, the number of values that are bit
for bit equal, and the zero-reference counts. As information (not gated):
mass action of order 0 and 1 is expected to be bit for bit equal; pairs,
Michaelis–Menten (derivatives) and Hill laws (evaluated as the converted
fraction instead of the 0.6.0 logistic form) may differ in the last bits.

**Downstream.** `integrations/diffsol` (its own workspace and lockfile)
builds with `--locked` and its tests pass (`cargo test --release --locked`
in that directory, with a temporary target directory).
`validation/check_ode.py --binary <crn-diffsol>` passes with its own
thresholds, unchanged, normally and under `python -O` with identical
reports.

### M5e (ODE right-hand side and Jacobians on the flat view): results

Engine of the M5e commits, release build, this host, under the CPU lock.

**Reference gate** (`crates/crn-craft/tests/kinetics_references.rs`,
report `validation/merge/kinetics/references.json`): **0 failures** in
6,252 values (3,126 per run) and 12 acceptance checks, as declared. The
network file (run 1) and the converted `model/1.0` (run 2) gave the same
numbers:

| Model | Values per run | Bit for bit equal | Largest relative difference (r ≠ 0) | Zero references (reproduced exactly) |
| --- | ---: | ---: | ---: | ---: |
| W1 | 4 | 4 | 0 | 0 (0) |
| W2 | 2,244 | 2,244 | 0 | 2,113 (2,113) |
| W3 | 168 | 155 | 2.6e-16 | 106 (106) |
| enzyme-mechanism | 64 | 64 | 0 | 25 (25) |
| immigration-death | 8 | 8 | 0 | 1 (1) |
| michaelis-menten | 20 | 16 | 1.7e-16 | 4 (4) |
| repressilator | 168 | 155 | 2.6e-16 | 106 (106) |
| vilar-noise | 450 | 450 | 0 | 356 (356) |

W4, W5, W6 and W7 are refused with `CRN-E080`, as in 0.6.0: W4's 33
Michaelis–Menten propensities take molecules of species they do not vanish
with; W5, W6 and W7 exceed the size cap. Every zero reference was
reproduced exactly, so the declared scale rule for zero references never
applied. The 30 values per run that differ in the last bits are the Hill
laws (evaluated as the converted fraction instead of the logistic form) and
the Michaelis–Menten derivatives.

**Downstream.** `integrations/diffsol` builds with `--locked`; its tests
pass, 10 of 10 (2 new: a converted model solves bit for bit like its network
model, and an expression rate of a `model/1.0` matches its closed form).
`validation/check_ode.py`: 50 of 50 checks, identical reports normally and
under `python -O`. Against the last stored 50-row report
(`validation/campaign/review-ode-jobs-ode.json`), 31 rows are identical and
19 rows (repressilator, the Hill fixtures, the Michaelis–Menten long
horizon) differ in the trailing digits of their errors (largest: the
repressilator DOP853 error 7.72221e-10 against 7.72216e-10); every
threshold holds with the same margin. CI runs neither `check_ode.py` nor the
diffsol build or tests.
