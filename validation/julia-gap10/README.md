# Ten further Julia-gap iterations — evidence and reproduction

Protocol: `docs/decisions/julia-gap10-2026-11-02.md`.
Base: `add6cca0fb146e8a40418432e963219514086024`, merged PR #7.
The preceding A01–A10 work is already in this base. This is C01–C10.

`C01.json` through `C10.json` retain every paired observation, numerical digest,
explicit target flag, acceptance decision, parent/candidate binary hash and build
cost. Matching `.patch`, `.build*.log`, `.tests*.log` and `.timing.log` files keep
candidate code and failed development runs. Rejected candidate patches are for
review; **do not apply them as part of the final patch**. Only C05, C06 and C09
remain in the delivered implementation.

- `final/`: release/MSRV/clippy/schema/interface/legacy/exact-law gate outputs.
- `comparison-validation.json`: final-engine independent-seed comparison checks.
- `W*.samples.json.gz`: all retained validation terminal samples, lossless, with
  hashes in the validation report.
- `comparison.json`: baseline/final/Julia timing observations and replay hashes.
- `all-method-timings.csv`: all 41 engine/API rows per workload, including slower
  choices, with all five observations and their median/minimum/maximum.
- `summary.json`: summary computed from those observations.
- `source-investigation.json`: hashes of inspected pinned Julia sources and host.
- `evidence-audit.json`: final source/report/sample hash verification, matching
  normal/-O outputs, and final native/wheel hashes.
- `regression-diagnostics.json`: independent follow-up pairs for every exposed
  prepared-row regression above 5%, plus native localization samples. Primary
  measurements and candidate decisions remain unchanged.

The full statistical gate uses 100,000 paths on each of seven fixtures for every
algorithm/wait combination: 7 million paths per pass, 210 marginals. The existing
family-alpha .001 DKW threshold and wrong-convention power checks are unchanged.
Normal and Python -O reports must match. Julia comparison validation uses distinct
seeds by engine index, 512 paths, 6,946 checks and family alpha .01. Boundaries,
counts and limitations are detailed in the new comparison document.

## Build and validate

Use Rust 1.88 (minimum 1.85), Python 3.12 and a clean environment. Install the
captured third-party Python dependencies, then build/install this source wheel:

```sh
python -m pip install -r validation/julia-gap10/python-requirements.txt
python -m pip install maturin==1.9.6
cargo build --release --workspace --examples --locked
maturin build --release --locked --out dist
python -m pip install --force-reinstall --no-deps dist/*.whl
cargo test --release --workspace --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
python validation/check_alternative_ssa.py --output final-laws.json
python -O validation/check_alternative_ssa.py --output final-laws-optimized.json
cmp final-laws.json final-laws-optimized.json
python validation/check_contracts.py --output contracts.json
python validation/check_alternative_contracts.py --binary target/release/crn-craft --output alternative-contracts.json
```

`validate_gap10.py` records the executed full gate list; its target/toolchain paths
are task-local and should be adjusted for another workspace. The independent
checkers themselves accept explicit binary/output paths. Default engine
requirements have not changed.

## Baseline wheel and direct Julia rerun

Build the unmodified base in a separate clean worktree and extract its wheel:

```sh
git worktree add --detach /tmp/crn-gap10-baseline add6cca0fb146e8a40418432e963219514086024
cd /tmp/crn-gap10-baseline
maturin build --release --locked --out /tmp/crn-gap10-baseline-wheels
python -m zipfile -e /tmp/crn-gap10-baseline-wheels/CRN_BASELINE_WHEEL.whl /tmp/crn-gap10-baseline-package
```

Replace the wheel placeholder with the actual generated name. Return to the final
checkout. Instantiate the retained Julia project with Julia 1.13.3. `JULIA_BIN`
and `JULIA_DEPOT` below are task-specific local paths:

```sh
"$JULIA_BIN" --project=validation/ssa-acceleration/julia -e 'using Pkg; Pkg.instantiate()'
python validation/compare_julia_gap10.py --julia "$JULIA_BIN" --depot "$JULIA_DEPOT" --binary target/release/crn-craft --baseline-package /tmp/crn-gap10-baseline-package --validation-only --output validation/julia-gap10/comparison-validation.json
# Finish every other heavy check before the next command.
python validation/compare_julia_gap10.py --julia "$JULIA_BIN" --depot "$JULIA_DEPOT" --binary target/release/crn-craft --baseline-package /tmp/crn-gap10-baseline-package --validated validation/julia-gap10/comparison-validation.json --output validation/julia-gap10/comparison.json
python validation/report_gap10.py
```

The baseline is loaded into an isolated process through PYTHONPATH; its IPC is
outside its measured prepared-call timer. This avoids mixing native extension
versions in one interpreter. The final harness checks all twelve baseline/final
numerical results for identical paths, while allowing explained work counters
to decrease. Julia uses pristine batch copies outside its timer. See the report
for the application-boundary differences. This is a direct same-host measurement,
not a projection from earlier hosts.

## Diagnostic regression follow-up

`validation/diagnose_gap10_regressions.py` records the six exposed prepared rows
in five fresh alternating pairs, two repeats per side. It then measures W4
DirectCR/inverse and DirectCR/Ziggurat natively at the base and retained C05,
C06 and C09 builds. Its `/tmp/crn-gap10-*` paths refer to the campaign's isolated
baseline package and preserved `bench_ssa` binaries; adjust them to rebuilt
artifacts when reproducing. Re-running this diagnostic overwrites its own JSON,
never the primary comparison report. The native Ziggurat pass was executed as
an additional standalone pass after the initial diagnostic; the script includes
both passes for future reproduction. Native sequential samples are diagnostic,
not a replacement for the paired primary experiment.
