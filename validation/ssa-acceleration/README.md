# Five-step SSA acceleration evidence

Protocol: `docs/decisions/ssa-acceleration-2026-11-02.md`.
Workload models and hashes were frozen before measurement in `models/`.
The implementation log in `docs/campaign/` explains failures and later fixes.

The two `final-laws*.json` files are the complete 210-marginal acceptance gate
on the delivered engine. Their bytes must match under normal and optimized
Python. Earlier `step4-laws*` reports precede the explicit absorption check and
remain development evidence. No numerical threshold was changed.

`julia-validation.json` contains all 6,946 distribution checks for 17 engines on
seven workloads, plus propensity/stoichiometry probes. The `W*.samples.json.gz`
files preserve all 512 final samples per engine per model, with SHA256 identities
in the report. `superseded-validation/` is an intentionally interrupted earlier
validation, retained separately and never used as a completed comparison.

`julia-comparison.json` contains the five warm timing observations, minima,
medians, maxima, event/work counters and own-engine replay hashes. One-shot and
prepared CRN API boundaries are separate. Builds and heavy checks finish before
these timings. Historical pre-acceleration evidence remains in
`validation/campaign/julia-comparison.json`.

## Reproduce

Use a release build and an installed current Python wheel. The optional Python
environment needs NumPy, SciPy, jsonschema and ruff. The executed environment is
pinned by the captured versions and requirements file; numerical validation does
not require Julia. The default engine has no new dependency.

```sh
python -m pip install -r validation/ssa-acceleration/python-requirements.txt
python -m pip install maturin==1.9.6
maturin build --release --locked --out dist
python -m pip install --force-reinstall --no-deps dist/*.whl
cargo build --release --workspace --locked
cargo test --release --workspace --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
python validation/check_alternative_ssa.py --output final-laws.json
python -O validation/check_alternative_ssa.py --output final-laws-optimized.json
cmp final-laws.json final-laws-optimized.json
python validation/check_alternative_contracts.py --binary target/release/crn-craft --output alternative-contracts.json
python validation/check_contracts.py --output contracts.json
python -O validation/check_contracts.py --output contracts-optimized.json
cmp contracts.json contracts-optimized.json
python validation/check_dependencies.py
```

Instantiate the pinned Julia project using Julia 1.13.3. `JULIA_BIN` and
`JULIA_DEPOT` below are task-specific local paths supplied by the user:

```sh
"$JULIA_BIN" --project=validation/ssa-acceleration/julia -e 'using Pkg; Pkg.instantiate()'
python validation/compare_julia_accelerated.py --julia "$JULIA_BIN" --binary target/release/crn-craft --depot "$JULIA_DEPOT" --validation-only --output validation/ssa-acceleration/julia-validation.json
python validation/compare_julia_accelerated.py --julia "$JULIA_BIN" --binary target/release/crn-craft --depot "$JULIA_DEPOT" --validated validation/ssa-acceleration/julia-validation.json --output validation/ssa-acceleration/julia-comparison.json
```

The validation-reuse option requires matching source hashes, including Rust and
Python engine source, and matching frozen workload hashes. The second invocation
still performs fresh per-engine warm-ups and validates every timing replay.
Julia's model construction/pristine copying and JSON IPC are outside its warm
timer. CRN's prepared timer includes result serialization/Python decoding;
one-shot additionally includes parsing and preparation. No archive writing is
timed. Julia event-normalized figures, if shown, use CRN event-count proxies.
