# M5g log: parameters and workers

Task M5g of [the plan](../../history/plan-0.7.0-engine-merge.md): the network
`ParameterTemplate` is replaced, for the downstream workers, by the rules
engine's parameter binding (`rules::Prepared::bind_parameters`, Python
`PreparedEnsemble.bind_parameters`); `integrations/jobs` and
`validation/parameter-optimization/worker` move to it with converted
models.

## What changed

- **`integrations/jobs`** (own workspace and lockfile, unchanged lock):
  - The model input is a `crn-craft.model/1.0` document (or any rules
    model). A `crn-craft.model/0.1` network model is refused with a message
    that names `crn-craft convert`, because the fit-data columns change
    meaning (below).
  - Per thread allowance, one prepared Static Direct `run/1.0` request per
    seed of `seed_set` (`Templates`, built lazily and cached); each
    candidate calls `bind_parameters` with exp(x) of `parameter_names` and
    runs `rules::execute`. The job's cancellation and deadline stop a run
    through a watcher thread that cancels the rules `Control` (the request
    of a prepared template cannot carry the remaining wall time).
    Budgets: `max_events` as before, no sampling or instance cap
    (`u64::MAX`), 256 MiB output.
  - `crn-craft.fit-data/0.1` columns, and `amplitude_species`, refer to the
    result's observables, sorted by name (for a converted model: the
    species names, sorted). The observables come from a cancelled run of
    the first template when the job is prepared, which also validates the
    model and the requests then, as the old `CompiledNetwork::compile` did.
  - Examples: `examples/fit.model.json` (the former inline immigration-death
    model, converted) and `examples/circuit.model.json` (converted); the
    columns of `examples/circuit.data.json` reordered from
    m1, p1, m2, p2, m3, p3 to m1, m2, m3, p1, p2, p3. Problem version `"3"`;
    `about` describes the inputs and the engine; `describe` lists the
    `crn-craft.model/1.0` input schema (`schemas/crn-craft.model-1.0.schema.json`
    replaces the 0.1 copy).
  - Tests: the amplitude penalty on `PathSeries` (ported), the examples are
    converted models whose columns match the observables, and a bound
    candidate equals a fresh preparation with the same values.
- **`validation/check_jobs.py`** (migrated in place; the worker has no old
  path left): the independent recomputation edits the model's parameters
  and prepares a fresh rules run (`rules.prepare(...).simulate()`), reads
  the statistics frames row-major and the amplitude observable from
  `Result.paths`; the broken ring sets the converted `tx3` propensity
  `alpha0 + alpha / (1 + (p2 / K)^n)` to its basal rate `alpha0` (0.6.0 set
  `vmax` to 0); report schema `crn-craft.jobs-validation/0.3`.
- **`validation/parameter-optimization/worker`** (own lockfile, unchanged):
  routes `fresh` (the network model converted with
  `rules::convert::from_network`, parameters set, prepared afresh) and
  `reuse` (`bind_parameters` on the prepared phase template), Static Direct
  with `wait` (default Ziggurat), moments reported in the model's type
  order (the network species order of a converted model). The init fields
  `algorithm` and `fast` and the command `validate_kernel` are gone (unknown
  field or command). The `generated` feature and its four modules
  `crn_*.rs` are removed (the deletions landed in `f2c5110` through the
  index; `bc7b5e1` completes the worker change). `rebop-macro` is
  unchanged.
- **Generated kernels: mapped to the interpreted path.** The old route ran
  a compiled `ssa::PropensityKernel` that read the bound coefficients at
  run time. The rules generators (`generate_rules_kernel`,
  `generate_rules_loop`) compile parameter values into the source and
  their fingerprint covers the model bytes, so `execute_with_rules_kernel`
  and `execute_with_rules_loop` refuse every rebound candidate
  (`CRN-R020` at `/model`); a parameter search would need one compilation
  per candidate. Both CRN routes are therefore interpreted. On these
  models the old generated route was no faster than the old reusable route
  (table below).
- **Benchmark scripts**: `generate_parameter_workers.py` writes only the
  rebop modules (identical output); `parameter_optimization.py` has the
  routes `crn_fresh` and `crn_reuse`, builds the `dynamic` and
  `rebop_macro` workers, hashes the rules Python and binding sources, reads
  the gate `evidence/parameter-binding.json` and compares the traces of two
  CRN routes; `report_parameter_optimization.py` keeps the 0.6.0 route
  names of the committed evidence it renders (output byte-identical with
  the old and the new script). The worker README says that the committed
  evidence is the 0.6.0 measurement.
- **Checks**: new `validation/check_parameter_binding.py` (the rules path:
  worker strictness, for each fixture and waiting time three candidates
  whose binding equals a fresh preparation in Python and whose worker
  `fresh` and `reuse` replies are identical and equal the Python binding,
  moments in type order, template not mutated, invalid values refused in
  Python and in the worker, concurrent bindings, cancellation, budget
  stops, thread independence, and the generated loop source changing with
  the parameter values). `validation/check_parameter_interface.py` keeps
  testing the network `ParameterTemplate` Python path until M6a/M8 and
  lost its worker parts (they moved). `validation/validate_parameters.py`
  builds the worker, runs both scripts normally and under `-O`, and gained
  `--evidence DIR` (default: the committed evidence directory), `--resume`
  and `--skip NAME` (recorded; the summary is then not `complete`).
- **`validation/merge/compare_parameter_workers.py`**: the declared speed
  comparison.

## Tests written first

- `integrations/jobs`: the three tests were written with the code; run
  against the network-model example files, the two new ones failed
  (`/home/claude/m07/m5cdg/m5g-jobs-first-fail.log`), then passed with the
  converted examples.
- `check_parameter_binding.py` was written before its first run against the
  migrated worker; its first run failed on an infinite value that JSON
  cannot carry to the worker (test bug; the case became a wrong value
  count).

## Measured

Declared in `validation/merge/RESULTS.md` (commit `f2c5110`) before the
run; report `validation/merge/parameter-workers.json`. Wall time of 129
worker calls (128 training evaluations, one held-out), median of 3, 1
thread, this host:

| Case | old fresh | old reuse | old generated | new fresh | new reuse | new reuse / old reuse | new reuse / old generated |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| birth_death | 0.181 s | 0.178 s | 0.209 s | 0.271 s | 0.296 s | 1.66 | 1.42 |
| sparse33 | 3.126 s | 3.023 s | 3.060 s | 5.227 s | 5.326 s | 1.76 | 1.74 |
| binding12 | 0.067 s | 0.045 s | 0.041 s | 0.075 s | 0.074 s | 1.65 | 1.81 |

The difference is the simulation (Static Direct with this branch's linear
selector against 0.6.0's cached Direct with `fast: true`); the rules binding
costs about a fresh preparation (bind seconds 0.013–0.129 over 129 calls,
against 0.002–0.024 for the old template). M4c/M4f, not on this branch,
change the Static Direct selector and per-event cost; rerun the comparison
after the merge.

## Checks

| Check | Result |
| --- | --- |
| `cargo test --release --locked` in `integrations/jobs` | 3 of 3 pass |
| `cargo build --release --locked` in `integrations/jobs` and in the worker | pass, lockfiles unchanged |
| worker `cargo clippy --all-features --all-targets --locked -- -D warnings` | pass (rebop 0.10.0 fetched) |
| `check_jobs.py` (runner `optimization-run` built from the pinned optimization-jobs `12f9d57`) | 15 of 15 rows, 17 jobs, 18 manifests; 14 objectives equal the recomputation (relative difference 0.0); `-O` identical; rerun on the final tree: the same report |
| `check_parameter_binding.py --worker` | pass: 18 bindings, 36 worker evaluations, 12 refusals, 12 concurrent bindings; `-O` identical |
| `check_parameter_interface.py` | pass: 180 bindings, 10 refusals, 12 concurrent bindings; `-O` identical |
| `validate_parameters.py --evidence /home/claude/m07/m5cdg/params-evidence --skip msrv-tests` | 38 checks pass, `msrv-tests` skipped (see Limitations); details in the [M5c log](M5c-log.md#checks) |

The other checks of the branch are in the [M5c log](M5c-log.md#checks).

## Limitations

- The rules binding re-serializes and prepares the whole model for each
  candidate; it is a convenience and a validation of names, not a fast
  path. Generated rules code cannot serve a parameter search (above).
- The jobs worker watches the job context every 2 ms to cancel a run; a
  run stops at the engine's next poll after that.
- `validate_parameters.py`'s `msrv-tests` (`cargo +1.85.0`) could not run
  here: the toolchain download (`static.rust-lang.org`) is refused by the
  host's network policy. CI tests 1.85.0.
- The committed optimization evidence, `docs/comparison-parameter-optimization.md`
  and `docs/reusable-parameters.md` still describe the 0.6.0 routes; the
  full optimization run (`parameter_optimization.py --run`) was not
  repeated (it needs the final gates in the committed evidence directory).

## Deviations from the plan

- The plan names `validate_parameters.py` as a gate script; it is the
  workspace's whole final-gate orchestrator (old-engine scripts included),
  run here with `--evidence` outside the tree and `msrv-tests` skipped.
- `check_parameter_interface.py` stays on the old path (as M5a did with
  `check_fsp.py`); the migrated checks are the new
  `check_parameter_binding.py`.
- The jobs worker refuses network models instead of converting them,
  because the fit-data column order would silently change.
- The generated route is removed rather than mapped to rules generators
  (reason above).

## Migration notes

- Jobs (`crn-job-worker`, problems `app.crn.fit-moments` and
  `app.crn.circuit-score`, version 3): pass `crn-craft.model/1.0` models
  (`crn-craft convert model.json`); order fit-data columns like the
  result's observables, sorted by name; `amplitude_species` names an
  observable. Objective values differ from version 2 (another engine and
  random stream).
- Parameter benchmark worker: `init` takes `wait` but no longer
  `algorithm` or `fast`; engines `fresh`, `reuse`, `rebop_macro`; no
  `generated` engine, feature or `validate_kernel` command. Benchmark
  routes `crn_fresh` and `crn_reuse` replace `crn_*_fast`.
- Rust users of `parameters::ParameterTemplate` move to
  `rules::prepare(model, request)?.bind_parameters(&values)` with a
  converted model; the network API stays until M8.
