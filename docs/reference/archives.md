# Archives

An archive is the sealed record of one run: the exact documents consumed and
the result, with hashes. `crn-craft run MODEL REQUEST --archive NEW_DIRECTORY`
and `PreparedEnsemble.archive` (Python) write one. `crn-craft verify DIRECTORY`
and `crn_craft.verify_archive` check one without simulating. The contracts
are `crn-craft.archive/1.0` (the manifest) and `crn-craft.archive-check/1.0`
(the report of writing or verifying). The schemas are `crn-craft describe
archive-1.0` and `archive-check-1.0`.

Archives of the removed network engine (`archive/0.1`, from `simulate`) are
not read by this release. Verify them with the release that wrote them.

## Files and completion

| File | Content |
| --- | --- |
| `model.json` | the exact model bytes |
| `request.json` | the exact request bytes, with their original `model.path` |
| `result.json` | the result as `crn-craft run` prints it: compact JSON with sorted keys and a trailing newline. It is `crn-craft.result/1.0` (or `rules-result`) for a simulation, and `crn-craft.fsp-result/1.0` for an `fsp` request |
| `manifest.json` | `crn-craft.archive/1.0`: length and SHA-256 of the three payloads |
| `COMPLETE` | the SHA-256 of the manifest bytes and a line feed, written last |

The output is a new directory, never an overwrite or a merge. Each file is
created exclusively. The verifier rejects extra or missing files, duplicates,
unknown fields or versions, mismatched lengths or digests, an invalid
completion marker and oversized content. Input documents keep the 16 MiB
limit of every contract.

The marker certifies that the write finished. A canceled or budget-limited
simulation can have a complete archive. If a process is killed before the
marker, the directory is not an archive. There is no resume. Integrity
checks are not signatures; they do not prove authenticity or protect
against a hostile process that replaces files concurrently.

## What verification checks

`verify` checks that the directory has exactly these regular files, that the
marker matches the manifest, the lengths and digests, that the model and
request still prepare, that the result parses strictly and is the engine's
own serialization, the hashes inside the result, and the relations between
request and result: identity, replicate accounting, terminations, kept
paths and statistics. It does not run the simulation again. To check that a
run reproduces, run the archived request and compare the result bytes.

## The report

`archive-check/1.0` has these fields, to be read together.

| Field | Meaning |
| --- | --- |
| `valid` | the archive was written, or verified, without error |
| `archive_complete` | every file, the manifest and the marker are present and consistent |
| `simulation_complete` | the archived result is complete: a simulation that ended `completed`, or a finite state projection |
| `termination` | the ensemble termination; null for a projection and for an invalid archive |
| `result` | the result's schema, for example `crn-craft.result/1.0` |
| `model_sha256`, `run_sha256` | the hashes of the archived bytes |
| `diagnostics` | why the archive is invalid (`CRN-E065`) |

`crn-craft run --archive` exits 0 only for a completed simulation.
`verify` exits 0 for any valid archive, including one that records an
incomplete simulation. Input or preflight failures create no archive. A
failure while writing returns a coded diagnostic and may leave a directory
without a marker, which must not be read as a result.

## Replaying an archive

`request.json` keeps its original `model.path`, which need not equal the
archive file name `model.json`. To replay, copy `model.json` next to a request
whose `model.path` names it and whose `model.sha256` is unchanged, or build a
new request and say that it is new. Numerical comparisons use `result.json`;
nothing in it depends on timing, thread count or the host.

## Provenance

The archive format adapts the durable-write and checksum conventions of
`optimization-jobs` at `cc6b03d594f73d1f3b38f088e7e3bd8663814d75`
(`crates/optimization-job-api/src/storage.rs`). No job protocol or worker
dependency is involved, and no implementation code was copied.
