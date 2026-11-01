# Result 0.2 and archive 0.1

The previous result 0.1 was a reserved envelope with no simulator emitter.
`describe result` now emits the actual result **0.2** schema. The unchanged
reserved schema remains available as `describe result-v1`. Model, run, check,
state, analysis 0.2, legacy analysis 0.1 and export 0.1 shapes are unchanged.
Discovery advances to **0.4** and advertises implemented Direct SSA.
Historical preflight/analysis report flags remain unchanged; use `describe`
for the current executable capability catalog.

## Files and completion

| File | Meaning |
|---|---|
| `model.json` | Exact consumed model bytes |
| `run.json` | Exact consumed run bytes, including its original relative model path |
| `result.json` | Result 0.2: numerical payload, accounting, termination and assessment |
| `provenance.json` | Compiler/package/lockfile/executable identities, OS/architecture, workers, elapsed time and workspace reservation |
| `manifest.json` | Archive 0.1; four unique fixed payload names, byte lengths and SHA-256 digests |
| `COMPLETE` | SHA-256 of the exact manifest bytes, followed by LF; written last |

Output is a new directory, never an overwrite or merge with existing data.
Each file is created exclusively and synchronized. On Unix, the directory is
synchronized after the marker. There are no symlinks or nested manifest paths.
The verifier rejects additional files, missing files, duplicates, unknown
fields/versions, mismatched lengths/digests, invalid completion markers and
oversized content. Strict duplicate-key and depth handling is retained;
input documents keep 16 MiB limits while result bytes have the archive cap.

This adapts the durable-write and checksum conventions inspected in
`optimization-jobs` at `cc6b03d594f73d1f3b38f088e7e3bd8663814d75`,
`crates/optimization-job-api/src/storage.rs`. No job protocol or worker
dependency is introduced. M6 can integrate the current shared runner later.
No third-party implementation code was copied.

The archive marker certifies that the **write** finished. A cancelled, failed
or budget-limited simulation may still have a complete archive. If a process
is killed before finalization, the directory has no valid completion marker;
it may contain only input snapshots or partial files. There is no resume or
per-event recovery journal in M2. Integrity checks do not provide signatures,
authenticity, biological validity or protection against hostile concurrent
filesystem replacement.

## Result relationships

Result 0.2 records model/run hashes, sampler version, seed, species order,
save times, requested/completed replicates, aggregate executed events,
per-replicate records, retained trajectories, optional statistics and
diagnostics. Every replicate ID from zero to requested-1 occurs once.
Each record reports whether it started, its stop reason, executed events,
observation time, available save-point count and optional error diagnostic.

`completed` requires all paths to reach the horizon. `exact_ensemble` is
reserved for this completed Direct result; it describes the stochastic method,
not zero sampling/rounding error. Numerical failures produce `failed_checks`;
other stops produce `incomplete`. An output-limited result is incomplete even
when all simulation work finished. It has accounting but no retained numerical
payload or unconditional statistics. Aggregate stop priority is failure,
cancellation, wall time, population limit, event limit, completion; output
removal takes precedence as `output_limit`. Individual reasons remain visible.

Statistics are present only for a fully completed ensemble when requested.
Verification checks dimensional/finite/domain constraints and consistency
between records, totals, stop reasons, sample counts, assessment and requested
outputs. When **all** trajectories are retained, it recomputes the same ordered
mean/variance/SE and rejects inconsistent summaries. When only some paths are
retained, full numerical recomputation is unavailable; checksums and structural
consistency are not a substitute for independent statistical validation.

Numerical reproducibility comparisons exclude run hashes when run bytes differ
only by worker allowance, and exclude `provenance.json`. The consumed run hash
must still identify those exact bytes. Retained `run.json` is an original
snapshot, so its original relative path need not equal the archive filename
`model.json`; replay it through the byte API or reconstruct a new checked run
reference without mislabelling the changed run as the original.

## CLI results

`simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]` returns an
`archive-check/0.1` report and exit 0 only for a completed simulation.
`verify DIRECTORY` returns the same report and exit 0 for a valid archive,
including valid incomplete evidence. Read `valid`, `archive_complete`,
`simulation_complete`, `termination` and `assessment` together.

Input/preflight failures create no archive. A write or finalization failure
returns a coded diagnostic and may leave an unmarked directory. The caller
must not reinterpret that directory as a completed result.
