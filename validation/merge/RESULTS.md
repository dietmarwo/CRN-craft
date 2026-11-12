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

### M4c (selectors): declaration

Written and committed before any run of the checks below. Mechanical smoke
runs of the scripts (300 replicates of the exact-law check, 2,000 events of
the benchmark, 2,000 replicates of the law gate on W1 and W4) checked only
that the scripts run; they are listed in the M4c log and decide nothing.
Engine under test: branch `merge-0.7` at the M4c implementation, CLI
`target/release/crn-craft` and the wheel built from the checkout, under the
CPU lock.

**1. Selector-threshold benchmark** (`selector_threshold.py --binary ...`,
output `validation/merge/selector-threshold.json`). Five families of
`crn-craft.model/1.0` models:

| Family | Network | Molecules per species | Rates |
| --- | --- | ---: | --- |
| `ring` | n species, X_i → X_(i+1 mod n) | 100 | 1 |
| `ring_small` | the same | 4 | 1 |
| `dense` | m species, every X_i → X_j, i ≠ j (n = m(m − 1)) | 32 | 1 |
| `michaelis` | the ring | 100 | whole propensity vmax X_i / (km + X_i), vmax 1, km 10 |
| `two_reads` | the ring | 100 | whole propensity X_i X_(i+1) / (1 + X_(i+1)) |

Grid: n = 8, 16, 32, 64, 128, 256, 512, 1,024, 2,048, 4,096 channels for the
ring families; m = 3, 4, 6, 8, 12, 16, 23, 32, 46, 64 (n = 6 … 4,032) for
`dense`. Every point runs `linear`, `composition_rejection`, `sorting` and
`rssa_cr` (`crn-craft.run/1.0`, Ziggurat waits, one replicate, one thread,
root seed 4242, end = 1,000,000 / a0 with a0 the initial total propensity).
A point's time is the median wall time of three `crn-craft rules` processes
after a warm-up, divided by the run's events. The rule, fixed from the file
without change:

- N_cr: the smallest grid n such that composition–rejection is faster than
  the linear scan at every point with at least n channels of every family;
- N_rssa: the smallest grid n such that RSSA-CR is faster than both the
  linear scan and composition–rejection at every such point of every family
  except `two_reads`;
- `auto` = `rssa_cr` when every rule's rate reads at most one species (no
  RSSA-CR bound is a point by construction) and n ≥ N_rssa; else
  `composition_rejection` when n ≥ N_cr; else `linear`. Sorting is never
  picked. If no grid n qualifies, the branch never applies.

**2. Exact-law check per selector** (`check_selectors.py --binary ...`, 2
threads, report `validation/merge/laws/selectors-exact.json`), the structure
of `validation/check_alternative_ssa.py`: the fixtures pure-decay,
immigration-zero and immigration-nonzero (closed-form laws, the seeds of their
run files) and the CME fixtures F1 dimerization, F2 heterodimer, F3
Michaelis–Menten with Hill recycling and F4 open annihilation of
`validation/check_cme.py` (laws by uniformization, their seeds), converted by
the engine under test. Each runs with the four selectors and both waiting
times (8 variants), n = 100,000 replicates, every path kept. Per variant,
fixture and positive report time: the largest CDF distance of the marginal
to the exact law is at most ε = sqrt(ln(2 · 168 / 0.001) / (2n)) = 0.00798
(168 tests, family α = 0.001); the reported mean and variance equal the
paths' moments to 2e-12 relative; for F1 and F4 every variant rejects the
falling-factorial law at some time (16 power checks). Expected failures:
about 0.001 for the family; **0 failures are expected**. A failure is
investigated once, not rerun with other seeds.

**3. Law gate per selector** (`check_laws.py --python ... --method
'{"kind": "static_direct", "selector": S}' --label m4c-S`) for S = `linear`,
`composition_rejection`, `sorting`, `rssa_cr`: all twelve models, the gate
seeds of `models.json`, 2 threads. Per selector N = 4,814 comparisons (4,636
against the references, 178 against FSP laws), 940 at t = 0; about 0.002
failures expected per selector: **0 failures are expected** for each. For
`linear`, the power sums and event counts of every model must also equal
those of `laws/m4-static_direct.json` (the M4b regression run, the same
linear scan under `crn-craft.run/1.0` before the selector existed). Failed
rows are kept; a failure is investigated once, then reported.

**4. Informative timings** (`compare_speed.py --old /home/claude/m07/bin060/crn-craft
--new target/release/crn-craft`): W1–W7 at the sizes of the speed gate, the
0.6.0 best method rerun in the same session, and the unified engine with
each selector and with `auto`; medians of three after a warm-up. No target;
the speed gate is M4d.

### M4c (selectors): results

**1. Selector-threshold benchmark** (2026-11-11, this host, under the CPU
lock; report `validation/merge/selector-threshold.json`). Nanoseconds per
event (median of three processes, warm-up excluded, process start and
lowering included), the fastest selector of each point last:

| Family | Channels | linear | composition_rejection | sorting | rssa_cr | Fastest |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `ring` | 8 | 77 | 109 | 98 | 70 | `rssa_cr` |
| `ring` | 16 | 79 | 119 | 107 | 72 | `rssa_cr` |
| `ring` | 32 | 103 | 117 | 124 | 75 | `rssa_cr` |
| `ring` | 64 | 133 | 132 | 154 | 83 | `rssa_cr` |
| `ring` | 128 | 212 | 140 | 202 | 80 | `rssa_cr` |
| `ring` | 256 | 344 | 143 | 297 | 85 | `rssa_cr` |
| `ring` | 512 | 647 | 153 | 499 | 100 | `rssa_cr` |
| `ring` | 1,024 | 1,414 | 207 | 856 | 107 | `rssa_cr` |
| `ring` | 2,048 | 2,550 | 203 | 1,529 | 132 | `rssa_cr` |
| `ring` | 4,096 | 4,977 | 290 | 3,138 | 204 | `rssa_cr` |
| `ring_small` | 8 | 74 | 160 | 102 | 196 | `linear` |
| `ring_small` | 16 | 84 | 167 | 112 | 190 | `linear` |
| `ring_small` | 32 | 112 | 172 | 127 | 192 | `linear` |
| `ring_small` | 64 | 141 | 165 | 156 | 207 | `linear` |
| `ring_small` | 128 | 213 | 174 | 210 | 203 | `composition_rejection` |
| `ring_small` | 256 | 349 | 193 | 326 | 210 | `composition_rejection` |
| `ring_small` | 512 | 652 | 197 | 470 | 241 | `composition_rejection` |
| `ring_small` | 1,024 | 1,268 | 210 | 875 | 236 | `composition_rejection` |
| `ring_small` | 2,048 | 2,477 | 263 | 1,565 | 269 | `composition_rejection` |
| `ring_small` | 4,096 | 4,834 | 338 | 3,113 | 352 | `composition_rejection` |
| `dense` | 6 | 85 | 177 | 134 | 106 | `linear` |
| `dense` | 12 | 113 | 237 | 185 | 110 | `rssa_cr` |
| `dense` | 30 | 163 | 330 | 273 | 124 | `rssa_cr` |
| `dense` | 56 | 229 | 425 | 369 | 128 | `rssa_cr` |
| `dense` | 132 | 393 | 631 | 557 | 145 | `rssa_cr` |
| `dense` | 240 | 570 | 762 | 793 | 161 | `rssa_cr` |
| `dense` | 506 | 998 | 1,140 | 1,244 | 196 | `rssa_cr` |
| `dense` | 992 | 1,554 | 1,462 | 1,789 | 238 | `rssa_cr` |
| `dense` | 2,070 | 3,005 | 2,279 | 3,235 | 309 | `rssa_cr` |
| `dense` | 4,032 | 5,686 | 3,347 | 5,256 | 414 | `rssa_cr` |
| `michaelis` | 8 | 131 | 167 | 155 | 78 | `rssa_cr` |
| `michaelis` | 16 | 141 | 165 | 163 | 78 | `rssa_cr` |
| `michaelis` | 32 | 183 | 174 | 184 | 80 | `rssa_cr` |
| `michaelis` | 64 | 197 | 185 | 218 | 82 | `rssa_cr` |
| `michaelis` | 128 | 270 | 205 | 264 | 108 | `rssa_cr` |
| `michaelis` | 256 | 422 | 200 | 387 | 94 | `rssa_cr` |
| `michaelis` | 512 | 711 | 228 | 542 | 102 | `rssa_cr` |
| `michaelis` | 1,024 | 1,316 | 259 | 955 | 160 | `rssa_cr` |
| `michaelis` | 2,048 | 2,547 | 332 | 1,683 | 160 | `rssa_cr` |
| `michaelis` | 4,096 | 5,126 | 472 | 3,163 | 239 | `rssa_cr` |
| `two_reads` | 8 | 182 | 212 | 207 | 399 | `linear` |
| `two_reads` | 16 | 232 | 228 | 221 | 360 | `sorting` |
| `two_reads` | 32 | 197 | 251 | 250 | 352 | `linear` |
| `two_reads` | 64 | 239 | 241 | 256 | 379 | `linear` |
| `two_reads` | 128 | 316 | 257 | 320 | 435 | `composition_rejection` |
| `two_reads` | 256 | 456 | 277 | 416 | 393 | `composition_rejection` |
| `two_reads` | 512 | 805 | 280 | 641 | 427 | `composition_rejection` |
| `two_reads` | 1,024 | 1,448 | 320 | 937 | 449 | `composition_rejection` |
| `two_reads` | 2,048 | 2,733 | 422 | 1,736 | 548 | `composition_rejection` |
| `two_reads` | 4,096 | 5,075 | 623 | 3,214 | 831 | `composition_rejection` |

Rule, applied as declared: **N_cr = 512** (composition–rejection lost to the
linear scan at `dense` n = 506, 1,140 against 998 ns, and won at every point
from 512 on); **N_rssa: none** (RSSA-CR was the fastest selector at every
point of `ring` and `michaelis` and at every `dense` point from n = 12 on,
but slower than composition–rejection at every `ring_small` point, by 2–4%
at the two largest: 269 against 263 ns at n = 2,048 and 352 against 338 ns
at n = 4,096). So `auto` picks composition–rejection from
512 channels on and the linear scan below, and never sorting or RSSA-CR:
`AUTO_COMPOSITION_REJECTION = 512`, `AUTO_RSSA_CR = None` in
`crates/crn-craft/src/rules/select.rs`.

Run history: a server restart stopped the first run after 27 of 50 points
(`ring`, `ring_small`, `dense` up to n = 506), each written when it finished.
The script gained `--resume` (finished points are kept, the others run; no
setting changed), and the remaining 23 points ran after the restart, one
family per command. The first 27 points ran on a CLI built from the M4c
sources before `cargo fmt` and an MSRV-compatible `next_up` helper, which
changes only the widening of higher-order channels (none in these models);
the other 23 on the CLI rebuilt from `eb336a6`.

**2. Exact-law check per selector** (2026-11-11, CLI built at `f1b7dac`, 2
threads, under the CPU lock, 119 s; report
`validation/merge/laws/selectors-exact.json`, byte-identical under
`python -O`): **168 of 168 CDF tests pass**, ε = 0.00798; every reported
mean and variance equals the paths' moments; all 16 power checks reject the
falling-factorial law (smallest largest distance 0.369).

| Selector | Waits | Tests | Failures | Largest CDF distance |
| --- | --- | ---: | ---: | ---: |
| `linear` | inverse | 21 | 0 | 0.00298 |
| `linear` | ziggurat | 21 | 0 | 0.00366 |
| `composition_rejection` | inverse | 21 | 0 | 0.00374 |
| `composition_rejection` | ziggurat | 21 | 0 | 0.00283 |
| `sorting` | inverse | 21 | 0 | 0.00358 |
| `sorting` | ziggurat | 21 | 0 | 0.00408 |
| `rssa_cr` | inverse | 21 | 0 | 0.00499 |
| `rssa_cr` | ziggurat | 21 | 0 | 0.00403 |

**3. Law gate per selector** (2026-11-11/12, wheel `_native` SHA-256
`1cad0420…` built at `f1b7dac`, Python backend, 2 threads, gate seeds,
under the CPU lock; reports `validation/merge/laws/m4c-<selector>.json`).
The runs were split into pieces of under ten minutes with `check_laws.py
--state` (stored chunks), because of server restarts; a report is written
when every chunk of every model is stored, so per-model wall times are not
recorded. Failures (comparisons) and the largest |z| per model:

| Model | `linear` | `composition_rejection` | `sorting` | `rssa_cr` |
| --- | --- | --- | --- | --- |
| W1 | 0 (24), 2.45 | 0 (24), 1.86 | 0 (24), 1.21 | 0 (24), 1.80 |
| W2 | 0 (396), 3.39 | 0 (396), 3.21 | 0 (396), 3.30 | 0 (396), 2.98 |
| W3 | 0 (72), 1.88 | 0 (72), 1.90 | 0 (72), 2.90 | 0 (72), 1.78 |
| W4 | 0 (192), 1.94 | 0 (192), 1.48 | 0 (192), 1.15 | 0 (192), 2.10 |
| W5 | 0 (2560), 3.88 | 0 (2560), 3.83 | 0 (2560), 3.51 | 0 (2560), 3.94 |
| W6 | 0 (320), 3.33 | 0 (320), 2.89 | 0 (320), 3.35 | 0 (320), 2.71 |
| W7 | 0 (960), 2.91 | 0 (960), 3.49 | 0 (960), 2.78 | 0 (960), 3.26 |
| enzyme-mechanism | 0 (80), 2.16 | 0 (80), 2.39 | 0 (80), 2.25 | 0 (80), 2.99 |
| immigration-death | 0 (20), 1.81 | 0 (20), 2.61 | 0 (20), 1.30 | 0 (20), 2.39 |
| michaelis-menten | 0 (40), 2.20 | 0 (40), 1.99 | 0 (40), 2.20 | 0 (40), 1.59 |
| repressilator | 0 (60), 1.88 | 0 (60), 1.90 | 0 (60), 2.90 | 0 (60), 1.78 |
| vilar-noise | 0 (90), 2.94 | 0 (90), 2.26 | 0 (90), 1.98 | failed in chunk 1, see below |
| **Total** | **0** (4814), 3.88 | **0** (4814), 3.83 | **0** (4814), 3.51 | **0** (4724), 3.94 |

- `linear`, `composition_rejection`, `sorting`: **gate passed**, 0 failures
  among N = 4,814 comparisons each (940 at t = 0).
- `linear`: the power sums, chunk event counts and samplers of all twelve
  models equal those of `laws/m4-static_direct.json` (the M4b regression
  run), as declared.
- `rssa_cr`: **gate failed for vilar-noise**; the other eleven models pass
  (0 failures among 4,724 comparisons). Chunk 1 (root seed 2,000,001) of
  vilar-noise did not complete: replicate 6298 stopped at t = 7.1624 with
  `CRN-R030` "clock failed to advance" (`/time`), so `check_laws.py`
  stopped (statistics withheld). Investigated once
  (`laws/m4c-rssa_cr-vilar-noise-failure.json`): the RSSA-CR upper total
  was 9,729.05 (exact total 3,353.06), and the Ziggurat wait unit / 9,729
  was below half an ulp of t (4.4e-16), so t + wait rounded to t. The
  linear scan has the same strict check, but RSSA-CR draws a wait for
  every trial, accepted or rejected, against a larger total, so a
  sub-ulp wait is more frequent (vilar-noise is the example where 0.6.0's
  RSSA-CR needed about 5.7 trials per event). The 0.6.0 network engine
  behaves the same way: its RSSA-CR trial loop calls
  `sampling::event_time(time, wait)?` before the composition draw and the
  acceptance test, and fails the replicate (`CRN-E062`) when the clock does
  not advance, also for a trial that would be rejected. Decision: keep the
  0.6.0 semantics (no engine change); the failed row is kept; the
  remaining 18 chunks were not run, and the factor and seeds are
  unchanged.

**4. Informative timings** (`compare_speed.py`, 2026-11-12, CLI built at
`4fe800d`, the kept 0.6.0 CLI rerun in the same session, one thread, the
sizes of the speed gate, medians of three; report
`validation/merge/speed-m4c.json`). Seconds, and in parentheses the ratio to
the same-session 0.6.0 best:

| Model | 0.6.0 best | `linear` | `composition_rejection` | `sorting` | `rssa_cr` | `auto` (picks) |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| W1 | `cached_direct/ziggurat/fast` 0.0499 | 0.0657 (1.32) | 0.0914 (1.83) | 0.0814 (1.63) | 0.0805 (1.61) | 0.0668 (1.34), `linear` |
| W2 | `direct_cr/ziggurat/fast` 0.0947 | 0.1574 (1.66) | 0.1328 (1.40) | 0.1513 (1.60) | 0.1143 (1.21) | 0.1577 (1.67), `linear` |
| W3 | `cached_direct/ziggurat/fast` 0.2156 | 0.4307 (2.00) | 0.6968 (3.23) | 0.5384 (2.50) | 0.5156 (2.39) | 0.4146 (1.92), `linear` |
| W4 | `rssa_cr/ziggurat/fast` 0.1052 | 1.7719 (16.84) | 2.2711 (21.59) | 2.2011 (20.92) | 0.2486 (2.36) | 1.8672 (17.75), `linear` |
| W5 | `direct_cr/ziggurat/fast` 0.0654 | 0.2262 (3.46) | 0.0942 (1.44) | 0.1792 (2.74) | 0.0831 (1.27) | 0.2252 (3.45), `linear` |
| W6 | `rssa_cr/ziggurat/fast` 0.0750 | 0.4131 (5.50) | 0.4878 (6.50) | 0.4448 (5.93) | 0.0823 (1.10) | 0.3511 (4.68), `composition_rejection` |
| W7 | `rssa_cr/ziggurat/fast` 0.0197 | 0.0312 (1.58) | 0.0307 (1.55) | 0.0418 (2.12) | 0.0209 (1.06) | 0.0302 (1.53), `linear` |

No target (the speed gate is M4d). The best selector is `rssa_cr` on W2, W4,
W5, W6 and W7 and `linear` on W1 and W3; against the 1.25 gate of M4d the
best ratios are W1 1.32, W2 1.21, W3 1.92 (`auto`, which runs the linear
scan), W4 2.36, W5 1.27, W6 1.10, W7 1.06.
Timings on this host varied between sessions by up to about 30% per row (an
earlier run of the same script, before `4fe800d`, gave for example W6
`composition_rejection` 0.379 s and W1 `linear` 1.16).

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

### M4f (per-event cost, RSSA-CR clock and threshold, `auto`): declaration

Written and committed before any run of the checks below. Engine under test:
branch `merge-0.7` at the M4f commits (clock fix `5d99986`, RSSA-CR
threshold reuse `dccc6da` and the result-preserving changes between and
after them), CLI `target/release/crn-craft` and the wheel built from the
checkout, under the CPU lock.

**1. RSSA-CR law gate, rerun as a new row.** Two M4f changes alter RSSA-CR
results under `crn-craft.run/1.0`: a positive wait below half an ulp of the
time is accepted (at most 1,000,000 consecutive steps; the M4c failure of
vilar-noise was such a wait), and the acceptance test reuses the
composition's threshold instead of drawing another uniform. So the whole
RSSA-CR gate runs again: `check_laws.py --python ... --method '{"kind":
"static_direct", "selector": "rssa_cr"}' --label m4f-rssa_cr`, all twelve
models, the gate seeds of `models.json` (chunk k: root seed 2,000,000 + k),
2 threads, report `validation/merge/laws/m4f-rssa_cr.json`. N = 4,814
comparisons (4,636 against the references, 178 against FSP laws), 940 at
t = 0; **0 failures are expected**. The M4c rows (eleven models passed,
vilar-noise failed) stay as recorded. A failure is investigated once, not
rerun with other seeds.

**2. The other selectors are unchanged.** Every other M4f change keeps the
arithmetic and the random words (byte-identical outputs on 156 requests:
twelve models × every selector, both waiting times, `next_reaction`, kept
paths and instances). The clock change differs from M4c only where a wait
leaves the rounded time unchanged, which failed a replicate in M4c; every
chunk of the M4c `linear`, `composition_rejection` and `sorting` gates
completed, so no such wait occurred there. Check: chunks 0 and 1 of every
model are rerun for each of the three selectors with `check_laws.py
--state` into new state directories, and their stored power sums S1–S4,
event counts and samplers must equal the M4c chunks
(`/home/claude/m07/m4c/state-<selector>`) exactly (72 chunks).

**3. Exact-law check per selector** (`check_selectors.py`, as declared in
M4c: 168 DKW tests, ε = 0.00798, 16 power checks, family α = 0.001), rerun
on the M4f engine because RSSA-CR changed; report
`validation/merge/laws/selectors-exact-m4f.json`. **0 failures are
expected**; the M4c report stays.

**4. `auto` benchmark** (`selector_rule.py`, output
`validation/merge/selector-rule.json`). Features of a lowered network: the
channel count n; `pinned`, whether some rule's rate reads two or more
lowered species (RSSA-CR bounds of such rules are points); the population
P̄, the mean initial count of the lowered species with a positive initial
count (∞ when every species starts at zero: such counts carry no scale).
Synthetic families of `crn-craft.model/1.0` models, each at the M4c grid of
channel counts (8 … 4,096; `dense` m = 3 … 64):

| Family | Network | Molecules per species | Rates | pinned |
| --- | --- | --- | --- | --- |
| `ring4`, `ring32`, `ring256` | X_i → X_(i+1 mod n) | 4, 32, 256 | 1 | no |
| `dense4`, `dense32`, `dense256` | every X_i → X_j, i ≠ j | 4, 32, 256 | 1 | no |
| `michaelis` | the ring | 100 | vmax X_i / (km + X_i), vmax 1, km 10 | no |
| `two_reads` | the ring | 100 | X_i X_(i+1) / (1 + X_(i+1)) | yes |

Each point runs `linear`, `composition_rejection` and `rssa_cr` as in M4c
(one replicate, one thread, root seed 4242, Ziggurat waits, about 10⁶
events, median of three processes after a warm-up, nanoseconds per event).
The target models are W1–W7 and the five examples, timed by `compare_speed.py
--selectors linear,composition_rejection,sorting,rssa_cr` at the speed-gate
sizes (report `validation/merge/speed-m4f-selectors.json`, after the
benchmark, same engine).

Rule family, with parameters N_small, N_large, N_cr ∈ {1, 2, 4, 8, 16, 32,
64, 128, 256, 512, 1,024, 2,048, 4,096, ∞} and P_cut ∈ {8, 16, 32}:
`auto` = `rssa_cr` if not pinned and n ≥ (N_small if P̄ < P_cut else
N_large); else `composition_rejection` if n ≥ N_cr; else `linear` (never
`sorting`). The parameters are chosen by this order of criteria, each
breaking the ties of the one before: (a) fewest target models whose pick is
more than 1.10 × slower than their fastest of the four selectors; (b) the
smallest largest ratio pick / fastest over the target models; (c) the
smallest largest ratio pick / fastest of the three measured selectors over
the synthetic points; (d) the largest N_cr, then N_large, then N_small;
(e) P_cut = 16, then 8, then 32. The chosen rule and every target model's
pick and ratio are recorded; no target is expected to miss (a) but a miss is
reported, not hidden. The M4c result (`selector-threshold.json`) stays as
the earlier row.

**5. Timing evidence for M4f** is the M4d speed gate (declared in M0c and
above under "Speed gate (M4d)"); its unified methods are the four selectors
(and `auto`) of Static Direct and `next_reaction`, all interpreted. The
compiled event loop of the unified engine (`compile_loop`, linear scan) is
timed at the Python boundary for information only.

### M4d (law gate of the run/1.0 default, worker speed rerun): declaration

Written and committed before either run below. Engine: `merge-0.7` at
`a21c358` (engine sources of `57cd574`, the commit of the passing speed-gate
rerun), wheel built from it with `maturin build --release --locked` and
installed in the engine venv, under the CPU lock.

**1. Law gate of the run/1.0 default.** `check_laws.py --python
<engine venv> --method '{"kind": "static_direct"}' --label m4d-auto --state
<dir>`, all twelve models, the gate seeds of `models.json`, 2 threads,
report `validation/merge/laws/m4d-auto.json`. The request names no selector,
so `auto` (the M4f rule) decides per model: `rssa_cr` for W2 and W4–W7,
`linear` for the others. N = 4,814 comparisons; **0 failures are expected**.
Also expected, because `auto` runs the selector it picks with the same random
words and every change since those gates keeps results: each model's
sampler, chunk event counts and power sums S1–S4 equal those of
`laws/m4f-rssa_cr.json` (models where `auto` picks `rssa_cr`) or of
`laws/m4c-linear.json` (models where it picks `linear`). A difference counts
as a failure of this check and is investigated once. The Python backend is
used because it stores finished chunks (`--state`), so a server restart does
not lose a run and each piece stays under ten minutes; the CLI backend has no
such state. The wheel and `target/release/crn-craft` are release builds of
the same sources, and I1/I2 tie the CLI to the same results.

**2. Worker speed rerun** (the M5g comparison after M4f). Script
`validation/merge/compare_parameter_workers.py`, the migrated benchmark
worker and the fixtures from branch `merge/M5cdg` at `97a60c3` (not yet
merged into `merge-0.7`), extracted into `/home/claude/m07/m4d/workers/`; the
migrated worker built in release mode with its own lockfile against the
`crn-craft` crate of this checkout (`a21c358`). Old workers: the M5g
binaries `/home/claude/m07/m5cdg/old-worker/dynamic` and `generated` (0.6.0
network engine, cached Direct, Ziggurat waits, `fast: true`). Protocol as in
the M5g declaration on that branch: the three frozen cases (`birth_death`,
`sparse33`, `binding12`), one worker process per route and repetition,
`init`, 128 training evaluations at the same log-uniform points for every
route and one held-out evaluation, 3 repetitions with rotated route order,
medians of the wall time of the 129 calls, 1 thread, under the CPU lock;
within each engine the routes must return identical moments. The migrated
worker's `crn-craft.run/1.0` request names no selector, so `auto` decides
(`linear` for `birth_death` and `binding12`, `rssa_cr` for `sparse33`, 66
channels starting empty). **Target:** for each case, new `fresh` / old
`fresh` ≤ 1.25 and new `reuse` / old `reuse` ≤ 1.25 (six ratios; the gate
passes only if all six do). New `reuse` / old `generated` is recorded for
information. The M5g row on `merge/M5cdg` stays as the earlier row. A miss
is recorded with its numbers and hot spots, not optimized further (owner
decision for M4d).

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

### M5g (parameter binding and workers): declaration

Written and committed before the first run of the measurement below. The
checks of `validation/check_jobs.py` (objectives equal an independent
recomputation within relative 1e-12, and its other rows) and of the new
`validation/check_parameter_binding.py` (exact equalities) keep their own
rules; they are not statistical gates.

**Worker speed comparison** (`validation/merge/compare_parameter_workers.py`;
no target, recorded). Workers: the 0.6.0-engine benchmark worker built from
this tree before the M5g change (the `dynamic` binary for the routes `fresh`
and `reuse`, the `generated` binary for `generated`; network model, cached
Direct, Ziggurat waits, `fast: true`, the worker defaults) and the migrated
worker (routes `fresh` and `reuse` on the rules engine: Static Direct,
Ziggurat waits, the model converted to `crn-craft.model/1.0`). For each of
the three frozen cases of `validation/parameter-optimization/fixtures.json`
(`birth_death`, `sparse33`, `binding12`; 32 training paths with seed 77201,
4,096 held-out paths with seed 880431, the frozen report times): one worker
process per route and repetition, `init`, then 128 training evaluations at
points drawn log-uniformly in [lower, upper] from
`numpy.random.default_rng(20261112 + case index)` (the same points for
every route), then one held-out evaluation at the true values. Recorded per
route: the sums of the replies' `bind_seconds` and `simulation_seconds` and
the wall time of the 129 calls seen from Python; 3 repetitions, route order
rotated per repetition, medians reported. 1 worker thread, under the CPU
lock, this host. Checked as part of the run (not timing): within each
engine the routes return identical moments at every point (new `fresh` =
new `reuse`; old `fresh` = old `reuse` = old `generated`).

### M5c (tau leaping on converted models): declaration

Written and committed before the first run. Engine under test: the rules
method `tau` of `crn-craft.run/1.0` (`{"kind": "tau", "epsilon": ε}`), run
by `crn-craft rules` on the network model of `validation/check_tau.py`
converted with `crn-craft convert`.

**Gate** (`validation/check_tau_rules.py`; the gates of `check_tau.py` with
the same tolerances): the model `TauControl` (X starts at 30, birth 20,
death 0.5 per molecule), end 2, report times 0 and 2, for ε = 0.1, 0.05 and
0.025 with root seeds 9101, 9102 and 9103, n = 10,000 replicates, every
path kept, 1 thread. Exact law at t = 2: Binomial(30, e⁻¹) convolved with
Poisson(40 (1 − e⁻¹)). For each ε: the Kolmogorov distance between the
empirical CDF of the final X and the exact CDF is at most
ε / 2 + sqrt(ln(6 / 0.001) / (2n)); the absolute difference between the
sample mean and the exact mean is at most ε · mean + 4.5 · sqrt(variance / n);
the result's termination is `completed`; its sampler is the approximate
`chacha8-u53/rules-approximate-cgp-poisson/2` and its `method` echoes the
request (the rules analog of 0.6.0's `approximate_ensemble` assessment).
**0 failures expected (3 rows).** The report is identical normally and under
`python -O`. Budgets (not gates): 10,000 events and 10⁸ sampling steps per
replicate. `check_tau.py` has no binomial variant, so `binomial` gets no
gate here. A failing row is investigated once; if it persists, the task
stops and reports.

**Leap-condition cost** (recorded, no target; `tau`'s propensity condition
evaluates every channel after one firing of every enabled channel,
O(channels²) per attempted leap): W5 and W6, converted, end 20, report times
0 and 20, 16 replicates, root seed 4242, 1 thread; `tau` with ε = 0.05
against `static_direct` (Ziggurat waits); wall time of `crn-craft rules`
(process), medians of 3, under the CPU lock; events per replicate recorded.

### M5g (parameter binding and workers): results

Release builds of this branch, this host, under the CPU lock.

**Worker speed comparison** (`validation/merge/parameter-workers.json`):
wall time of the 129 calls (128 training evaluations and the held-out one),
median of 3, seconds.

| Case | old fresh | old reuse | old generated | new fresh | new reuse | new fresh / old fresh | new reuse / old reuse | new reuse / old generated |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| birth_death | 0.181 | 0.178 | 0.209 | 0.271 | 0.296 | 1.49 | 1.66 | 1.42 |
| sparse33 | 3.126 | 3.023 | 3.060 | 5.227 | 5.326 | 1.67 | 1.76 | 1.74 |
| binding12 | 0.067 | 0.045 | 0.041 | 0.075 | 0.074 | 1.11 | 1.65 | 1.81 |

Bind seconds (sum over the 129 calls): old `reuse` 0.002–0.024, old `fresh`
0.013–0.154, new `fresh` 0.016–0.138, new `reuse` 0.013–0.129: the rules
binding re-serializes and prepares the whole model, so it costs about what a
fresh preparation costs. Simulation dominates (sparse33: 2.9 s old, 5.0–5.1 s
new): the new routes run Static Direct with the linear selector of this
branch against 0.6.0's cached Direct with `fast: true`. Within each engine
the routes returned identical moments at every point, as declared. The old
generated route was no faster than the old reusable route on these models.

**Checks**: `check_jobs.py` passes 15 of 15 rows (14 objectives equal the
independent recomputation with relative difference 0.0), normally and under
`python -O` with identical reports; `check_parameter_binding.py` passes (18
bindings, 36 worker evaluations, 12 refused inputs, 12 concurrent bindings);
`check_parameter_interface.py` (network path) passes (180 bindings), both
with identical reports under `python -O`.

### M5c (tau leaping on converted models): results

Release build of `0007ea9`, this host, under the CPU lock.

**Gate** (`validation/check_tau_rules.py`, report
`validation/merge/tau/check_tau_rules.json`): **0 failures in 3 rows**;
identical reports normally and under `python -O`.

| ε | CDF error | CDF envelope | mean error | mean envelope | events | random words | passed |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.1 | 0.00426 | 0.07086 | 0.00729 | 3.888 | 736,991 | 1,293,328 | yes |
| 0.05 | 0.00660 | 0.04586 | 0.06371 | 2.072 | 735,205 | 1,572,298 | yes |
| 0.025 | 0.00629 | 0.03336 | 0.05809 | 1.164 | 736,709 | 1,493,418 | yes |

For comparison (not a gate), `check_tau.py` on the network engine in the
same build: CDF errors 0.0227, 0.0148, 0.0081 with 73,776, 278,259 and
730,139 leaps. The rules `tau` leaps far less on this model (about two
random words per event): its propensity condition and its species bound
(g = 2, which refuses a leap when ε · n / 2 < 1) send most steps to the
exact fallback, so its laws are close to the exact ones.

**Leap-condition cost** (`validation/merge/tau/leap-condition-cost.json`):

| Model | Channels | Method | Wall time (s, median of 3) | Events per replicate | Random words per replicate | tau / static_direct |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| W5 | 510 | tau, ε = 0.05 | 3.228 | 20,437 | 40,876 | 14.1 |
| W5 | 510 | static_direct | 0.229 | 20,420 | 41,536 | |
| W6 | 992 | tau, ε = 0.05 | 4.962 | 14,365 | 28,731 | 12.9 |
| W6 | 992 | static_direct | 0.386 | 14,357 | 29,208 | |

At the declared ε = 0.05 neither model leaps (two words per event, as an
exact step): W5's counts of 4 make every channel critical, and W6's counts
of 32 give ε · n / 2 = 0.8 < 1, so the species bound refuses a leap before
the propensity condition runs. The 13–14× cost is therefore the exact
fallback inside `tau`, which recomputes every channel's rate and the
criticality of every channel at each step, not the channels² condition.
Supplementary row, not declared and run once: W6 with ε = 0.1 (where
ε · n / 2 = 1.6, so leaps are attempted), 1 replicate: 9.97 s, 14,043
events, 561,290 random words, against about 0.024 s per replicate for
Static Direct (about 400×). Each attempted leap evaluates the 992 channels
at the state after one firing of each enabled channel, twice (criticality
and change pass): about 2 · 992² ≈ 2.0 million channel rates per attempt.

### M4f (per-event cost, RSSA-CR clock and threshold, `auto`): results

**1. RSSA-CR law gate rerun** (2026-11-12, wheel `_native` SHA-256 `877542df…` of the M4f sources from `dccc6da` on (as recorded in the report), Python backend, 2 threads, gate seeds, under the CPU lock, seven pieces of under ten minutes with `--state`, 05:39–06:44; report `validation/merge/laws/m4f-rssa_cr.json`): **0 failures among N = 4,814 comparisons** (940 at t = 0), largest |z| 3.88. Failures (comparisons) and the largest |z| per model:

| Model | `rssa_cr`, M4f | `rssa_cr`, M4c (earlier row) |
| --- | --- | --- |
| W1 | 0 (24), 1.86 | 0 (24), 1.80 |
| W2 | 0 (396), 3.08 | 0 (396), 2.98 |
| W3 | 0 (72), 2.18 | 0 (72), 1.78 |
| W4 | 0 (192), 1.77 | 0 (192), 2.10 |
| W5 | 0 (2560), 3.88 | 0 (2560), 3.94 |
| W6 | 0 (320), 2.87 | 0 (320), 2.71 |
| W7 | 0 (960), 3.22 | 0 (960), 3.26 |
| enzyme-mechanism | 0 (80), 2.42 | 0 (80), 2.99 |
| immigration-death | 0 (20), 1.57 | 0 (20), 2.39 |
| michaelis-menten | 0 (40), 2.03 | 0 (40), 1.59 |
| repressilator | 0 (60), 2.18 | 0 (60), 1.78 |
| vilar-noise | 0 (90), 2.36 | failed in chunk 1 (`CRN-R030`) |
| **Total** | **0** (4814), 3.88 | **0** (4724), vilar-noise not compared |

vilar-noise completes every chunk, including chunk 1 (root seed 2,000,001), where M4c stopped. The M4c rows stay as recorded above.

**2. The other selectors are unchanged** (2026-11-12, same wheel, `check_laws.py --state` into new state directories, chunks 0 and 1 of every model; W5's chunks hold 5,000 replicates, so `--replicates 20000` ran its chunks 0–3): **78 of 78 chunks identical** to the M4c chunks (root seed, sampler, event count and power sums S1–S4), 26 per selector for `linear`, `composition_rejection` and `sorting` (`/home/claude/m07/m4f/unchanged/compare-<selector>.json`). Declared were 72 chunks (two per model); the four W5 chunks cover the same 20,000 replicates as two chunks of the other models.

**3. Exact-law check per selector** (2026-11-12, `target/release/crn-craft` of the M4f sources before the `auto` change, 2 threads, under the CPU lock; report `validation/merge/laws/selectors-exact-m4f.json`, byte-identical under `python -O`): **168 of 168 CDF tests pass**, ε = 0.00798; every reported mean and variance equals the paths' moments; all 16 power checks reject the falling-factorial law (smallest largest distance 0.370).

| Selector | Waits | Tests | Failures | Largest CDF distance (M4f) | (M4c) |
| --- | --- | ---: | ---: | ---: | ---: |
| `linear` | inverse | 21 | 0 | 0.00298 | 0.00298 |
| `linear` | ziggurat | 21 | 0 | 0.00366 | 0.00366 |
| `composition_rejection` | inverse | 21 | 0 | 0.00374 | 0.00374 |
| `composition_rejection` | ziggurat | 21 | 0 | 0.00283 | 0.00283 |
| `sorting` | inverse | 21 | 0 | 0.00358 | 0.00358 |
| `sorting` | ziggurat | 21 | 0 | 0.00408 | 0.00408 |
| `rssa_cr` | inverse | 21 | 0 | 0.00464 | 0.00499 |
| `rssa_cr` | ziggurat | 21 | 0 | 0.00407 | 0.00403 |

The three unchanged selectors give the M4c distances (same fixtures, seeds and words); `rssa_cr` changed with the threshold reuse.

**4. `auto` benchmark** (2026-11-12, CLI built at `bbbef17`, under the CPU lock; report `validation/merge/selector-rule.json`). Nanoseconds per event (median of three processes after a warm-up, process start and lowering included), the fastest of the three measured selectors, and the pick of the fitted rule with its ratio to the fastest:

| Family | Channels | linear | composition_rejection | rssa_cr | Fastest | `auto` (ratio) |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `ring4` | 8 | 59 | 152 | 169 | `linear` | `linear` (1.00) |
| `ring4` | 16 | 77 | 141 | 170 | `linear` | `linear` (1.00) |
| `ring4` | 32 | 114 | 195 | 181 | `linear` | `linear` (1.00) |
| `ring4` | 64 | 154 | 150 | 164 | `composition_rejection` | `linear` (1.03) |
| `ring4` | 128 | 260 | 171 | 188 | `composition_rejection` | `rssa_cr` (1.10) |
| `ring4` | 256 | 481 | 176 | 198 | `composition_rejection` | `rssa_cr` (1.13) |
| `ring4` | 512 | 872 | 187 | 198 | `composition_rejection` | `rssa_cr` (1.06) |
| `ring4` | 1,024 | 1,706 | 199 | 204 | `composition_rejection` | `rssa_cr` (1.03) |
| `ring4` | 2,048 | 3,394 | 243 | 247 | `composition_rejection` | `rssa_cr` (1.01) |
| `ring4` | 4,096 | 6,985 | 335 | 342 | `composition_rejection` | `rssa_cr` (1.02) |
| `ring32` | 8 | 60 | 130 | 96 | `linear` | `rssa_cr` (1.58) |
| `ring32` | 16 | 96 | 126 | 95 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 32 | 108 | 142 | 95 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 64 | 159 | 138 | 91 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 128 | 268 | 144 | 95 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 256 | 477 | 157 | 100 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 512 | 961 | 178 | 112 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 1,024 | 1,794 | 195 | 126 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 2,048 | 3,501 | 253 | 186 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring32` | 4,096 | 7,023 | 316 | 285 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 8 | 73 | 126 | 82 | `linear` | `rssa_cr` (1.12) |
| `ring256` | 16 | 78 | 155 | 98 | `linear` | `rssa_cr` (1.26) |
| `ring256` | 32 | 144 | 179 | 106 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 64 | 163 | 126 | 87 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 128 | 267 | 218 | 141 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 256 | 505 | 153 | 111 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 512 | 943 | 162 | 106 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 1,024 | 1,718 | 173 | 113 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 2,048 | 3,617 | 280 | 193 | `rssa_cr` | `rssa_cr` (1.00) |
| `ring256` | 4,096 | 6,975 | 309 | 234 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 6 | 71 | 165 | 177 | `linear` | `linear` (1.00) |
| `dense4` | 12 | 85 | 257 | 193 | `linear` | `linear` (1.00) |
| `dense4` | 30 | 180 | 382 | 191 | `linear` | `linear` (1.00) |
| `dense4` | 56 | 213 | 456 | 179 | `rssa_cr` | `linear` (1.19) |
| `dense4` | 132 | 382 | 663 | 193 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 240 | 589 | 804 | 204 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 506 | 1,138 | 1,149 | 224 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 992 | 2,255 | 1,748 | 235 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 2,070 | 4,275 | 1,989 | 310 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense4` | 4,032 | 8,099 | 2,708 | 334 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 6 | 66 | 149 | 93 | `linear` | `linear` (1.00) |
| `dense32` | 12 | 93 | 214 | 145 | `linear` | `rssa_cr` (1.56) |
| `dense32` | 30 | 190 | 283 | 138 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 56 | 197 | 422 | 145 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 132 | 400 | 609 | 128 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 240 | 650 | 724 | 224 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 506 | 1,200 | 954 | 187 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 992 | 2,148 | 1,412 | 259 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 2,070 | 4,226 | 1,857 | 373 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense32` | 4,032 | 8,229 | 2,547 | 435 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 6 | 61 | 139 | 89 | `linear` | `linear` (1.00) |
| `dense256` | 12 | 92 | 195 | 87 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 30 | 131 | 264 | 94 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 56 | 197 | 365 | 89 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 132 | 435 | 607 | 89 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 240 | 627 | 787 | 98 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 506 | 1,313 | 1,089 | 125 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 992 | 2,247 | 1,576 | 127 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 2,070 | 4,376 | 1,927 | 143 | `rssa_cr` | `rssa_cr` (1.00) |
| `dense256` | 4,032 | 8,315 | 2,663 | 192 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 8 | 126 | 143 | 72 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 16 | 128 | 132 | 69 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 32 | 149 | 131 | 73 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 64 | 259 | 161 | 99 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 128 | 304 | 167 | 80 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 256 | 543 | 180 | 99 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 512 | 959 | 184 | 101 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 1,024 | 1,884 | 209 | 109 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 2,048 | 3,679 | 265 | 183 | `rssa_cr` | `rssa_cr` (1.00) |
| `michaelis` | 4,096 | 7,326 | 441 | 228 | `rssa_cr` | `rssa_cr` (1.00) |
| `two_reads` | 8 | 138 | 198 | 312 | `linear` | `linear` (1.00) |
| `two_reads` | 16 | 163 | 236 | 310 | `linear` | `linear` (1.00) |
| `two_reads` | 32 | 190 | 256 | 303 | `linear` | `linear` (1.00) |
| `two_reads` | 64 | 250 | 223 | 314 | `composition_rejection` | `linear` (1.12) |
| `two_reads` | 128 | 354 | 221 | 362 | `composition_rejection` | `composition_rejection` (1.00) |
| `two_reads` | 256 | 593 | 246 | 417 | `composition_rejection` | `composition_rejection` (1.00) |
| `two_reads` | 512 | 1,013 | 274 | 375 | `composition_rejection` | `composition_rejection` (1.00) |
| `two_reads` | 1,024 | 1,968 | 343 | 546 | `composition_rejection` | `composition_rejection` (1.00) |
| `two_reads` | 2,048 | 3,990 | 372 | 590 | `composition_rejection` | `composition_rejection` (1.00) |
| `two_reads` | 4,096 | 7,576 | 603 | 852 | `composition_rejection` | `composition_rejection` (1.00) |

Target models (`compare_speed.py --selectors linear,composition_rejection,sorting,rssa_cr`, same CLI, speed-gate sizes, seconds, median of three; report `validation/merge/speed-m4f-selectors.json`):

| Model | Channels | P̄ | pinned | `linear` | `composition_rejection` | `sorting` | `rssa_cr` | Fastest | `auto` picks (ratio) |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | --- |
| W1 | 2 | ∞ | no | 0.0515 | 0.0875 | 0.0683 | 0.0734 | `linear` | `linear` (1.00) |
| W2 | 66 | ∞ | no | 0.1997 | 0.1339 | 0.1501 | 0.1115 | `rssa_cr` | `rssa_cr` (1.00) |
| W3 | 12 | 5.0 | no | 0.2798 | 0.5510 | 0.3729 | 0.4124 | `linear` | `linear` (1.00) |
| W4 | 66 | 20.0 | no | 1.1032 | 1.3205 | 1.1238 | 0.1130 | `rssa_cr` | `rssa_cr` (1.00) |
| W5 | 510 | 4.0 | no | 0.3178 | 0.0832 | 0.1794 | 0.0747 | `rssa_cr` | `rssa_cr` (1.00) |
| W6 | 992 | 32.0 | no | 0.5061 | 0.3239 | 0.4059 | 0.0708 | `rssa_cr` | `rssa_cr` (1.00) |
| W7 | 126 | 13.3 | no | 0.0351 | 0.0301 | 0.0318 | 0.0224 | `rssa_cr` | `rssa_cr` (1.00) |
| enzyme-mechanism | 3 | 30.0 | no | 0.0478 | 0.1039 | 0.0699 | 0.1249 | `linear` | `linear` (1.00) |
| immigration-death | 2 | ∞ | no | 0.0559 | 0.0883 | 0.0569 | 0.0858 | `linear` | `linear` (1.00) |
| michaelis-menten | 1 | 50.0 | no | 0.0435 | 0.0559 | 0.0404 | 0.0595 | `sorting` | `linear` (1.08) |
| repressilator | 12 | 5.0 | no | 0.3417 | 0.6841 | 0.4683 | 0.4954 | `linear` | `linear` (1.00) |
| vilar-noise | 16 | 1.0 | no | 1.0273 | 2.5791 | 1.8835 | 6.3203 | `linear` | `linear` (1.00) |

Rule, chosen by the declared criteria: **N_small = 128, N_large = 8, N_cr = 128, P_cut = 8**. (a) no target model is more than 1.10 × slower than its fastest selector; (b) the largest target ratio is 1.08 (michaelis-menten, one channel, `linear` against `sorting`); (c) the largest grid ratio is 1.58 (`ring32` at 8 channels, where the linear scan is fastest; N_large = 16 would give `michaelis` at 8 channels to the linear scan, 1.75 × slower than RSSA-CR). So `auto` = `rssa_cr` when no rule's rate reads two or more species and n ≥ 8 (n ≥ 128 when P̄ < 8); else `composition_rejection` from 128 channels; else `linear`: `AUTO_RSSA_CR_SMALL = Some(128)`, `AUTO_RSSA_CR_LARGE = Some(8)`, `AUTO_COMPOSITION_REJECTION = Some(128)`, `AUTO_POPULATION_CUT = 8.0` in `crates/crn-craft/src/rules/select.rs` (commit `9c4b5f9`). The M4c rule (composition–rejection from 512 channels, else linear) stays as the earlier row.

### M4d (speed gate): results

Declared in M0c ("Speed gate (M4d)" above). Both runs: `baseline.py` with the kept 0.6.0 CLI (`/home/claude/m07/bin060/crn-craft`) and venv in four parts of three models (each under ten minutes), each part checked with `--check-against baseline-0.6.0.json`: **63 of 63 rows per part with equal non-timing fields, no difference**, in both runs; then `compare_speed.py --old <0.6.0 CLI> --new target/release/crn-craft --baseline <the four parts>` (the unified engine's four selectors, `auto` and `next_reaction`; gated CLI boundary, median of three after a warm-up), all under the CPU lock, one session per run. Seconds; ratio = unified best / 0.6.0 best of the same session.

**Run 1** (2026-11-12, 07:31–07:57, CLI of `9c4b5f9`; reports `validation/merge/baseline-0.6.0-m4d-run1.json`, `validation/merge/speed-m4d-run1.json`): **gate failed on W5 (1.28)**; the other eleven models were within 1.25.

| Model | 0.6.0 best (same session) | s | `linear` | `composition_rejection` | `sorting` | `rssa_cr` | `auto` (picks) | `next_reaction` | Unified best | s | Ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: |
| W1 | `direct/ziggurat/fast` | 0.0654 | 0.0504 | 0.1231 | 0.0648 | 0.0970 | 0.0557 (`linear`) | 0.1150 | `linear` | 0.0504 | 0.77 |
| W2 | `direct_cr/ziggurat/fast` | 0.0894 | 0.1927 | 0.1250 | 0.1474 | 0.1378 | 0.1032 (`rssa_cr`) | 0.1391 | `auto` | 0.1032 | 1.15 |
| W3 | `cached_direct/ziggurat/fast` | 0.2254 | 0.2651 | 0.7067 | 0.4446 | 0.4378 | 0.2799 (`linear`) | 0.4490 | `linear` | 0.2651 | 1.18 |
| W4 | `rssa_cr/ziggurat/fast` | 0.1046 | 0.8528 | 1.2250 | 1.2388 | 0.1124 | 0.1306 (`rssa_cr`) | 1.2278 | `rssa_cr` | 0.1124 | 1.08 |
| W5 | `direct_cr/ziggurat/fast` | 0.0602 | 0.2992 | 0.1138 | 0.1807 | 0.0771 | 0.0770 (`rssa_cr`) | 0.0832 | `auto` | 0.0770 | **1.28, failed** |
| W6 | `rssa_cr/inverse/fast` | 0.0641 | 0.4997 | 0.3064 | 0.3824 | 0.0653 | 0.0610 (`rssa_cr`) | 0.3144 | `auto` | 0.0610 | 0.95 |
| W7 | `direct_cr/ziggurat/fast` | 0.0202 | 0.0368 | 0.0294 | 0.0286 | 0.0195 | 0.0196 (`rssa_cr`) | 0.0224 | `rssa_cr` | 0.0195 | 0.97 |
| enzyme-mechanism | `direct/ziggurat/fast` | 0.0761 | 0.0551 | 0.1170 | 0.0680 | 0.1080 | 0.0505 (`linear`) | 0.0714 | `auto` | 0.0505 | 0.66 |
| immigration-death | `cached_direct/ziggurat/fast` | 0.0695 | 0.0488 | 0.0906 | 0.0573 | 0.0857 | 0.0511 (`linear`) | 0.0679 | `linear` | 0.0488 | 0.70 |
| michaelis-menten | `cached_direct/ziggurat/fast` | 0.0575 | 0.0381 | 0.0634 | 0.0596 | 0.0640 | 0.0481 (`linear`) | 0.0416 | `linear` | 0.0381 | 0.66 |
| repressilator | `cached_direct/ziggurat/fast` | 0.2676 | 0.3152 | 0.5714 | 0.5168 | 0.4846 | 0.3302 (`linear`) | 0.4796 | `linear` | 0.3152 | 1.18 |
| vilar-noise | `cached_direct/ziggurat/fast` | 1.0723 | 0.9841 | 2.5166 | 1.8086 | 6.3688 | 1.0370 (`linear`) | 1.9928 | `linear` | 0.9841 | 0.92 |

**Fix round** (one, as the plan allows): W5 under callgrind (RSSA-CR, 510 channels, 4 molecules per species): 2.1 trials per event, half of the composition draws rejected because the group upper bounds (8) are exact powers of two, each candidate evaluated exactly, about 10 random words (ChaCha8) per event. Changes that keep the arithmetic and the words were limited to the composition sampler: the ceilings of the active groups are now computed when the cumulative masses are rebuilt, not per draw (`57cd574`, byte-identical on the 156 harness requests and in the I1 and generated-loop replays; about 1 % fewer instructions on W5). Fewer words or trials per event would change RSSA-CR's results or the random-number consumption and were not attempted.

**Run 2, the rerun** (2026-11-12, 08:11–08:33, CLI of `57cd574`; reports `validation/merge/baseline-0.6.0-m4d.json`, `validation/merge/speed-m4d.json`): **gate passed**, every model within 1.25.

| Model | 0.6.0 best (same session) | s | `linear` | `composition_rejection` | `sorting` | `rssa_cr` | `auto` (picks) | `next_reaction` | Unified best | s | Ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: |
| W1 | `direct/ziggurat/fast` | 0.0482 | 0.0635 | 0.0873 | 0.0649 | 0.1048 | 0.0511 (`linear`) | 0.1188 | `auto` | 0.0511 | 1.06 |
| W2 | `direct_cr/ziggurat/fast` | 0.0906 | 0.2091 | 0.1281 | 0.1382 | 0.1031 | 0.1144 (`rssa_cr`) | 0.1571 | `rssa_cr` | 0.1031 | 1.14 |
| W3 | `cached_direct/ziggurat/fast` | 0.2248 | 0.2713 | 0.6758 | 0.4233 | 0.3962 | 0.2707 (`linear`) | 0.4548 | `auto` | 0.2707 | 1.20 |
| W4 | `rssa_cr/ziggurat/fast` | 0.1119 | 0.7341 | 1.2303 | 1.3746 | 0.1620 | 0.1084 (`rssa_cr`) | 1.2262 | `auto` | 0.1084 | 0.97 |
| W5 | `direct_cr/ziggurat/fast` | 0.0614 | 0.3211 | 0.0945 | 0.1736 | 0.0770 | 0.0685 (`rssa_cr`) | 0.0776 | `auto` | 0.0685 | 1.12 |
| W6 | `rssa_cr/inverse/fast` | 0.0639 | 0.5309 | 0.3188 | 0.3820 | 0.0930 | 0.0599 (`rssa_cr`) | 0.3106 | `auto` | 0.0599 | 0.94 |
| W7 | `cached_direct/ziggurat/fast` | 0.0185 | 0.0397 | 0.0456 | 0.0315 | 0.0216 | 0.0231 (`rssa_cr`) | 0.0274 | `rssa_cr` | 0.0216 | 1.17 |
| enzyme-mechanism | `cached_direct/ziggurat/fast` | 0.0731 | 0.0683 | 0.1034 | 0.0622 | 0.1220 | 0.0670 (`linear`) | 0.0788 | `sorting` | 0.0622 | 0.85 |
| immigration-death | `cached_direct/ziggurat/fast` | 0.0687 | 0.0555 | 0.0881 | 0.0744 | 0.0917 | 0.0474 (`linear`) | 0.0762 | `auto` | 0.0474 | 0.69 |
| michaelis-menten | `direct/ziggurat/fast` | 0.0561 | 0.0397 | 0.0579 | 0.0392 | 0.0672 | 0.0368 (`linear`) | 0.0426 | `auto` | 0.0368 | 0.66 |
| repressilator | `cached_direct/ziggurat/fast` | 0.2809 | 0.3178 | 0.6362 | 0.4214 | 0.4785 | 0.3281 (`linear`) | 0.4793 | `linear` | 0.3178 | 1.13 |
| vilar-noise | `cached_direct/ziggurat/fast` | 1.0514 | 1.0986 | 2.8116 | 1.6707 | 6.1344 | 0.9451 (`linear`) | 2.0337 | `auto` | 0.9451 | 0.90 |

**Timing noise.** The gate takes the smallest median among the unified methods, and `auto` runs the same computation as the selector it picks (outputs equal apart from the method name), so each model has two timings of one computation. In run 2 they differ by up to 1.55 × (W6 0.0599 against 0.0930 s, W4 0.1084 against 0.1620 s), and the 0.6.0 best of W1 moved from 0.0654 s (run 1) to 0.0482 s (run 2). With the slower of the two timings, run 2's ratios would be W1 1.32, W2 1.26, W3 1.21, W4 1.45, W5 1.26, W6 1.46, W7 1.24, immigration-death 0.81, michaelis-menten 0.70, repressilator 1.17, vilar-noise 1.04 (enzyme-mechanism 0.85, best `sorting`). Pooled over both runs (each side's smaller median over the two runs; run 1 used the CLI before the 1 % fix): W1 1.05, W2 1.15, W3 1.18, W4 1.04, W5 1.14, W6 0.94, W7 1.05, enzyme-mechanism 0.69, immigration-death 0.69, michaelis-menten 0.66, repressilator 1.18, vilar-noise 0.90. The 0.6.0 side is a minimum over 21 methods with close variants, so it profits from the same effect. The gate passes by its declared rule; for W2, W3, W5, W7 and the repressilator (1.12–1.20) the margin is within this host's run-to-run noise. W5's change from 1.28 to 1.12 is mostly such noise (the fix saves about 1 % of its instructions).

**Compiled event loops** (informative, 2026-11-12 08:45, wheel of `57cd574`, `/home/claude/m07/m4d/compiled.py` under the CPU lock with `CRN_CRAFT_SOURCE` at the checkout; output `/home/claude/m07/m4d/compiled.json`). Python boundary: `rules.ensemble(...)` at the speed-gate sizes with `selector: linear` (the scan the generated loops run), `compile_loop`, `loop.verify(run)`, then the median of three `run.simulate()` calls after a warm-up, interpreted and compiled; 0.6.0 best = the smallest Python-boundary median of the same-session baseline rerun. Seconds:

| Model | Interpreted (`linear`) | Compiled loop | Compiled / interpreted | 0.6.0 best (Python) | s | Compiled / 0.6.0 best |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| W1 | 0.0477 | 0.0353 | 0.74 | `cached_direct/ziggurat/fast` | 0.0440 | 0.80 |
| W2 | 0.1915 | 0.1727 | 0.90 | `direct_cr/ziggurat/fast` | 0.0882 | 1.96 |
| W3 | 0.2627 | 0.1808 | 0.69 | `sorting_direct/ziggurat/fast` | 0.2502 | 0.72 |
| W4 | 0.7850 | refused: source over 1 MiB (`CRN-R020`) | | `rssa_cr/ziggurat/fast` | 0.1014 | |
| W5 | 0.2916 | refused: source over 1 MiB (`CRN-R020`) | | `direct_cr/inverse/fast` | 0.0527 | |
| W6 | 0.5022 | refused: source over 1 MiB (`CRN-R020`) | | `rssa_cr/ziggurat/fast` | 0.0468 | |
| W7 | 0.0306 | 0.0350 | 1.15 | `cached_direct/ziggurat/fast` | 0.0154 | 2.27 |
| enzyme-mechanism | 0.0454 | 0.0370 | 0.81 | `cached_direct/ziggurat/fast` | 0.0411 | 0.90 |
| immigration-death | 0.0453 | 0.0410 | 0.90 | `sorting_direct/ziggurat/fast` | 0.0403 | 1.02 |
| michaelis-menten | 0.0389 | 0.0309 | 0.79 | `sorting_direct/ziggurat/fast` | 0.0266 | 1.16 |
| repressilator | 0.3319 | 0.2080 | 0.63 | `cached_direct/ziggurat/fast` | 0.2604 | 0.80 |
| vilar-noise | 1.0012 | 0.6126 | 0.61 | `sorting_direct/ziggurat/fast` | 1.3030 | 0.47 |

The compiled loops run the linear scan, so on W2 and W7, where RSSA-CR is the best interpreted method, they are slower than the 0.6.0 best; W4–W6 exceed the generated-source limit. Build seconds (cache hits after the first): W1 49.3 (first build), the others 0.1–6.4.

**Remaining hot spots** (callgrind, CLI of `57cd574`; for a follow-up release):

- W1 (2 channels, `linear`, 1.06): the two ChaCha8 words per event (wait and selection; about 3.7 ns per word against 1.5 ns for 0.6.0's xoshiro256++) and the Ziggurat wait.
- W3 and the repressilator (12 channels, Hill rates, `linear`, 1.20 and 1.13): one mispredicted scan exit per event (scans of more than four rates keep the early exit), the running-sum resume after each event, the Hill programs of the affected channels, two ChaCha8 words per event.
- W2, W5 and W7 (`rssa_cr`, 1.14, 1.12, 1.17): trials per event (W5 2.1), composition rejection (W5: about half of the draws, because the group upper bounds are exact powers of two), exact evaluation of every candidate, about 10 words per event on W5.
- W4 (`rssa_cr`, 0.97): eight-species count updates and interval tests per event; range refreshes on 7 % of the events.
- Levers that change results or the build, left to a follow-up: RSSA-CR interval widths or group ceilings that are not exact powers of two, prefix-sum selection as a new selector, a cheaper generator under a new sampler identity (D4 fixes ChaCha8), and LTO or one codegen unit (which also changes the 0.6.0 comparison build).

### M4d (law gate of the run/1.0 default, worker speed rerun): results

**1. Law gate of the run/1.0 default** (declared above; Python backend,
wheel of `57cd574`'s engine sources, `--state` pieces under the CPU lock;
report `validation/merge/laws/m4d-auto.json`): **passed**: 4,814 comparisons, 0 failures, largest |z| 3.88 (W5). `auto` picked `rssa_cr` for W2 and W4–W7 and `linear` for the other seven models, as declared. Sampler, chunk layout, event counts and power sums S1–S4 equal those of `laws/m4f-rssa_cr.json` (the `rssa_cr` picks) and `laws/m4c-linear.json` (the `linear` picks) for all twelve models. The `linear` selector keeps the sampler identity `chacha8-u53-ziggurat/rules-static-direct/1` of rules-run Static Direct.

**2. Worker speed rerun** (declared above; run 2026-11-12 09:26–09:28 under
the CPU lock; report `validation/merge/parameter-workers-m4d.json`). Wall
time of the 129 calls, median of 3 repetitions, seconds:

| Case | old fresh | old reuse | old generated | new fresh | new reuse | new fresh / old fresh | new reuse / old reuse | new reuse / old generated (information) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| birth_death | 0.1687 | 0.1684 | 0.1494 | 0.1944 | 0.1848 | 1.15 | 1.10 | 1.24 |
| sparse33 | 2.9089 | 2.8939 | 2.9794 | 3.1464 | 3.0428 | 1.08 | 1.05 | 1.02 |
| binding12 | 0.0520 | 0.0343 | 0.0334 | 0.0538 | 0.0531 | 1.04 | **1.55, failed** | 1.59 |

**Target failed** on one of the six gated ratios: `binding12` new `reuse`
/ old `reuse` = 1.55 (the other five: 1.04–1.15). Of the 0.019 s gap,
binding is about 0.011 s (new 0.0137 s, about 106 µs per call, against
0.0024 s: `Prepared::bind_parameters` re-serializes and prepares the whole
model per call) and simulation about 0.008 s (32 short paths per call, where
per-call run preparation and per-path setup outweigh the events). Within
each engine the routes returned identical moments at every point. Recorded,
not optimized (owner decision for M4d); follow-up: a bind path that patches
the parameter values of an already prepared run.

### M10c (GillesPy2 comparison on the unified engine): declaration

Committed before the first run. The measurement repeats
`docs/comparison-gillespy2.md` on the unified engine. It is informative:
there is no speed target.

- **Script:** `validation/compare_gillespy2.py` as ported in M8a, unchanged:
  workloads W1–W4 of `benchmark_campaign.models()`, 64 paths over t in
  [0, 100], reports at 0 and 100, one thread, root seed 4242; the agreement
  checks use 512 paths per simulator with the script's fixed seeds (777,
  778, 779).
- **Build:** the CLI and the wheel of the branch `merge/M10c` at this
  commit, `cargo build --release --workspace --locked` and
  `maturin build --release --locked` from the root; the wheel installed in a
  virtual environment with the package set of the merge suite (Python
  3.12.3, GillesPy2 1.8.3, NumPy 2.3.5, SciPy 1.18.1, matplotlib 3.11.2,
  SCons 4.11.1); g++ 13.3.
- **Runs:** four parts, `--only W1` to `--only W4`, each under the CPU lock,
  with the script's defaults: five interleaved repeats of every CRN-craft
  row and of GillesPy2's `SSACSolver`, three repeats of its
  `NumPySSASolver`; a row's time is the median. A part that is interrupted is
  rerun whole. The four reports are merged by workload into
  `validation/merge/gillespy2-comparison-0.7.json`; their other fields must
  be equal. The figure is plotted from the merged report with
  `--plot-only`.
- **Rows:** the CLI process (`crn-craft run`, wall time of the process), the
  Python calls `crn.ensemble(...).simulate()` (preparation included) for the
  default Static Direct request (`auto` selector, Ziggurat waits), the
  checked configuration (`linear`, inverse waits) and each selector with
  Ziggurat waits, GillesPy2's C++ solver (built once; compile time reported
  separately) and its NumPy solver.
- **Agreement rule:** the script reports 20 two-sample Kolmogorov–Smirnov
  p-values (up to three species per workload, for the checked and the
  default configuration) and 8 pooled chi-square p-values against the exact
  Poisson laws of W1 and W2 (both simulators, both configurations). With a
  family-wise error rate of 0.01 over these 28 tests (Bonferroni), a p-value
  below 0.01 / 28 = 3.6e-4 is a failed row. On W4 every species must stay
  equal within each path, in both simulators. Failed rows are kept, and
  nothing is rerun with other seeds.
