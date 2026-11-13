# M5d log: archives of rules runs

Task M5d of [the plan](../../history/plan-0.7.0-engine-merge.md): seal a rules run
in a new directory as `crn-craft.archive/1.0` and verify it, reusing the
manifest and path-safety code of `archive.rs`; `verify` keeps reading the
network engine's 0.1 archives until M8.

## What changed

- **`crates/crn-craft/src/archive/v1.rs`** (new, `crn_craft::archive::v1`,
  a child module of `archive.rs`, so the old file keeps its content and the
  new code reuses its private helpers: `write_new` (exclusive create and
  sync), `read_regular` (bounded regular files, no symlinks), `encoded`
  (byte-limited compact JSON), the `CRN-E065` constructor).
  - Files: `model.json` and `request.json` (exact bytes; the request keeps
    its own `model.path`), `result.json` (the result exactly as
    `crn-craft rules MODEL REQUEST` prints it: compact JSON with sorted keys
    and a newline), `manifest.json` (`crn-craft.archive/1.0`: path, length
    and SHA-256 of the three payloads, in that order) and `COMPLETE` (the
    manifest's SHA-256 and a newline, written last; the directory is synced
    after it on Unix).
  - `run_files(model, request, output, control)`: the request file's
    `model.path` must name the model file below the request's directory
    (`CRN-R020` at `/model/path`, as `crn-craft rules`), then `write`.
    `write(model_bytes, request_bytes, output, control)`: a simulation
    request is prepared first (input errors create no directory), then the
    directory is created (`CRN-E065` if it exists), the inputs are written,
    the run executes, the result, manifest and marker follow. An FSP
    request (`run/1.0` with method `fsp`) is solved first and written as a
    whole, so a failed solve leaves nothing.
  - `verify(directory)`: exactly the five regular files within their size
    limits (16 MiB inputs, 1 GiB result, 8 KiB manifest, 65-byte marker),
    the marker, a strict manifest naming each payload once, lengths and
    digests; the inputs must still prepare; the result must fit the
    request's `max_output_bytes`, parse strictly and be byte for byte the
    engine's serialization of itself; then the relations between request
    and result: identity (version, both hashes, sampler, seed, types,
    rules, observables, report times, replicate count, method, continuous
    columns, taken from a **cancelled execution of the archived inputs**,
    which simulates nothing and uses the engine's own construction), every
    record (ids, times, budgets, rule events summing to events, unstarted
    and complete records, diagnostics of halted ones), the aggregates, the
    ensemble termination and the exact list of ensemble diagnostics
    (`CRN-R050`, the inputs' `CRN-R051` warnings, `CRN-R045`), kept paths
    (which replicates, lengths, and the initial row when the first report
    time is 0, taken from one probe replicate of at most one event), final
    instances and observations, and statistics (present exactly for a
    complete ensemble that asked for them, estimators, quantile
    probabilities, dimensions, finite values, standard errors equal to
    sqrt(variance / samples), quantiles nondecreasing in p, and with every
    path kept the mean, variance and quantiles recomputed from the paths bit
    for bit). An FSP result is checked for its hashes, the flat view's
    labels and initial state, the bounds, the report times and every
    snapshot's dimensions.
  - Report `crn-craft.archive-check/1.0` (`ArchiveReport`): `valid`,
    `archive_complete`, `simulation_complete` (termination `completed`, or
    an FSP result), `result` (the result's schema string), `termination`
    (rules termination; null for FSP), both hashes, `diagnostics`.
  - `is_archive(directory)`: the manifest's schema string, read leniently,
    is `crn-craft.archive/1.0`; without a readable manifest, a
    `request.json` without `run.json`.
- **`archive.rs`**: `pub mod v1;`, `AnyReport` (untagged 0.1 or 1.0 report)
  and `verify_any`, which dispatches with `v1::is_archive`. The 0.1 writer
  and verifier are unchanged (the mutation anchors still match).
- **CLI**: `crn-craft rules MODEL REQUEST [--archive NEW_DIRECTORY]
  [--cancel-file PATH]` (options in any order, each once; anything else is
  `CRN-E031`). With `--archive` it prints the archive-check/1.0 report and
  exits 0 only for a valid archive of a completed run, like the old
  `simulate`. `crn-craft verify DIRECTORY` calls `verify_any`. The help text
  and `describe` list the kinds `archive-1.0` and `archive-check-1.0`.
- **Schemas**: `schemas/archive-1.0.schema.json` and
  `schemas/archive-check-1.0.schema.json` (from `crn_craft::schema`);
  `check_contracts.py` snapshots them. Discovery: `rules::describe()`
  gains `contracts_1_0.archive` (command, verify command, manifest and
  report versions, schema kinds, Python names). `rules/mod.rs` otherwise
  changes only `mod quantile` to `pub(crate) mod quantile` (the verifier
  recomputes quantiles with the engine's function).
- **Python**: native `rules_archive_json(model, request, output,
  cancel_file=None)` and `verify_json` through `verify_any`;
  `crn_craft.rules.archive(model, request, directory, *, cancel_file=None)`,
  `PreparedEnsemble.archive(directory, *, cancel_file=None)` and
  `crn_craft.rules.verify_archive(directory)`, each returning the report as
  a dictionary. `python -m crn_craft verify` reads both versions.

## Python: what 0.6.0 offered and the choice made

0.6.0 had no archive function in the Python API: archives were written and
checked by `python -m crn_craft simulate|simulate-cached|verify`, which call
the native `archive_json(model_path, run_path, output, cached, cancel_file)`
and `verify_json(directory)` and return the report. Mirrored: a function
that runs the request and writes the archive, returning the report
(`rules.archive`, and `PreparedEnsemble.archive` for the exact bytes a
prepared run holds), and a verify function returning the report.
`Result.archive(directory)` was not added: a `Result` holds decoded frames,
not the engine's result bytes, so it would have to run the request again,
which does not reproduce a cancelled run, a wall-time stop or a replicate
stopped by a callback; a `bind_parameters` candidate has only reconstructed
model bytes. `PreparedEnsemble.archive` refuses such a candidate with
`ValueError`. A `python -m crn_craft rules ... --archive` command is left to
M7, which defines the unified Python command line.

## Tests written first

- `crates/crn-craft/tests/archives_v1.rs` (17 tests) against stub entry
  points returning a failed report: 17 of 17 failed
  (`/home/claude/m07/m5cdg/m5d-first-fail.log`). Ported from
  `tests/archives.rs`: inputs preserved and overwrite refused; interrupted
  or tampered archives (no marker, extra file, a directory, `{}` as result,
  `../model.json` in a resealed manifest); digests checked independently
  of JSON (zero marker, an equal-length string change); refreshed hashes
  that hide false accounting or statistics (25 mutations, each with its
  expected message); duplicate manifest entries and duplicate keys
  (`CRN-E002`), a 0.1 manifest version (`CRN-E003`); cancellation and
  event limits (valid archives of incomplete runs); symlinked payload and
  root; bad model path, output reservation, hash and syntax leave no
  directory; an oversized payload is refused before reading; partial
  evidence (started flag, missing diagnostic, a path longer than its
  record, a wrong initial row, padding beyond the output budget); a large
  report grid (100,000 report times, 8 threads, cancelled). New: the byte
  API writes the same files; non-canonical result bytes and unknown fields;
  every method (`tree`, `static_direct` with inverse waits,
  `next_reaction`, `tau`, `binomial`, `hybrid`) under `run/1.0` and
  `rules-run/0.6`; FSP archives and six FSP tamper cases; `verify_any` on
  both versions; schema snapshots and discovery.
- `crates/crn-craft-cli/tests/cli.rs`:
  `rules_archives_runs_and_verify_reads_both_archive_versions` failed
  against the old CLI (`m5d-cli-first-fail.log`), then passed.
- `validation/check_rules_python.py`: `check_archives` (prepared-run
  archive equals the simulated result, the input bytes, `rules.archive`
  equals `PreparedEnsemble.archive`, a re-serialized model refused at
  `/model/sha256` without a directory, an existing directory refused,
  tampering detected, an FSP archive, a cancelled run, a 0.1 archive read
  through `rules.verify_archive`, a candidate refused).

## Checks

See the [M5c log](M5c-log.md#checks) for the checks of the whole branch
(fmt, clippy, all tests including I1 and I3, rustdoc, ruff, I2 and the
Python scripts); all pass.

## Deviations from the plan

- The archive option sits on the existing `crn-craft rules` command (the
  unified `run` is M7); the verify command is the existing `verify`, which
  dispatches on the archive version.
- No `provenance.json`: the plan's archive holds model, request, result and
  manifest, and the rules result has no timing, host or executable fields
  that would need one. The completion marker of 0.1 is kept.
- The request file is named `request.json` (0.1: `run.json`), which also
  tells the versions apart without a manifest.
- The output budget applies to the result (the rules contract's
  `max_output_bytes`), not to the whole directory as in 0.1.
- 0.1's "saved state lies beyond observation time" relation is not ported:
  the rules engine fills report times before the drawn next event time and
  records the last event's time as `stop_time` for a halted replicate, so a
  kept row can lie after `stop_time` legitimately.
- FSP requests can be archived too (not in the plan).
- The verifier derives the expected identity from a cancelled execution and
  the initial row from a one-event probe instead of duplicating the
  engine's header logic, so it follows the selector strings of M4c without
  a change.

## Limitations

- Verification prepares the inputs three times (validation, cancelled
  header, probe) and holds the parsed result and its re-serialization in
  memory: about three times the result size for a 1 GiB result.
- Observations (trackers and event logs) are checked for order, replicate
  and tracker count only; their contents are not recomputed. Final
  instances are checked for count, type and attributes, not against the
  kept paths. Without every path kept, the statistics are checked for
  consistency only.
- Integrity checks are not signatures; a forger who recomputes a
  consistent result for the same inputs is not detected.
- `docs/archive-contract.md` still describes 0.1 only; M10b documents 1.0.

## Migration notes

- New: `crn-craft rules MODEL REQUEST --archive NEW_DIRECTORY`,
  `crn_craft::archive::v1::{run_files, write, verify, is_archive}`,
  `crn_craft::archive::{verify_any, AnyReport}`, contracts
  `crn-craft.archive/1.0` and `crn-craft.archive-check/1.0`
  (`describe archive-1.0`, `describe archive-check-1.0`), Python
  `crn_craft.rules.archive`, `PreparedEnsemble.archive`,
  `crn_craft.rules.verify_archive`.
- `crn-craft verify` and `python -m crn_craft verify` read both versions;
  for a 1.0 archive they print `crn-craft.archive-check/1.0`.
- Unchanged until M8: `simulate`, `simulate-cached`, `archive::simulate_files`,
  `archive::verify_archive` and the 0.1 contracts.
