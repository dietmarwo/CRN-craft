# Catalyst / JumpProcesses reference benchmark

This optional environment supports `validation/compare_julia.py`. It does not
change the CRN-craft engine or its default dependencies. The comparison uses the
frozen W1–W4 fixtures from `validation/benchmark_campaign.py`.

## Reproduce

From the repository root, with Rust and Julia 1.13.3 available:

```sh
cargo build --release --workspace --locked
python3 -m venv .venv-julia-comparison
. .venv-julia-comparison/bin/activate
python -m pip install -r validation/julia/requirements.txt
python -m pip install .
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=validation/julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
python validation/compare_julia.py \
  --julia "$(command -v julia)" \
  --binary target/release/crn-craft \
  --output validation/campaign/julia-comparison-rerun.json
```

Use `--depot /absolute/path` for a separate Julia depot. Set `JULIA_DEPOT_PATH`
to the same path for the preceding `Pkg.instantiate()` call. For a shorter
integration check, append `--smoke --rounds 2` and choose another output file.
Smoke runs use 16 validation paths and two timed paths; their times are not the
published benchmark.

The full run performs 512-path statistical checks before five interleaved
64-path timing rounds. Builds and other heavy workloads should finish before
timing. Each numerical engine has one worker. The Python process communicates
with a persistent Julia process; Julia reports its own elapsed batch time, so
JSON IPC and process startup are excluded from warm Julia times.

`Project.toml` and `Manifest.toml` pin the registered Julia packages.
`runtime-provenance.json` records the official Julia binary and verified digest.
`python-environment.txt` is an observed environment record, including the local
wheel's digest; install from source as above instead of reusing its temporary
wheel path. The native build used rustc/cargo 1.88.0 and maturin 1.9.6.

## Model and RNG details

The Python harness emits `generated_models.jl` from the frozen fixtures.
Catalyst constructs the actual `JumpProblem`s, specialized mass-action jumps
and dependency information. Custom laws use the full-propensity `=>` arrow.
The W4 translation omits availability gates only on its invariant equal-count
class; it is not a general-purpose CRN-craft-to-Catalyst exporter.

Within a batch, a single seeded Xoshiro stream advances through serial Julia
solves. A pristine problem is copied before, and outside, each batch timer to
reset SortingDirect's learned order. The CRN-craft RNG is different; seeds do
not imply identical paths between engines. Repeated timed batches must replay
their own numerical digests. Native Direct and cached Direct must agree exactly.

The predeclared acceptance protocol is
`docs/decisions/julia-comparison-2026-11-02.md`. All final samples, validation
distances, medians, ranges and hashes are retained in the JSON report. Julia
stderr is written beside each output. No statistical cutoff is fitted to the
observed samples.

## Development evidence

The retained `smoke-01` failure was an unsupported `maxiters` keyword passed to
JumpProcesses 10's `SSAStepper`. It was removed. The retained `smoke-02` failure
was the harness inferring an overly narrow Julia dictionary value type, then
rejecting a scalar elapsed time; the result dictionary now explicitly allows
heterogeneous values. Neither was a numerical solver failure. Both occurred
before the full benchmark and did not change the acceptance thresholds.
The subsequent `smoke-03` and full runs passed. The full run recorded all 239
distribution checks and 140 timing observations; its Julia stderr was empty.

`crn-tests.log`, `crn-clippy.log`, `crn-fmt.log`, the two `contracts*.json` files,
and `dependencies.log` record checks of the unmodified native engine. Julia
package installation and compilation are recorded in `installation.log`.
