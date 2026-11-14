# M6b log: `Result.save(path)` (FR1b)

Task M6b of [the plan](../../history/plan-0.7.0-engine-merge.md), the second part of
the feature request FR1: `result.save(file_name)` writes the model details
and the simulation results to one file, markdown or JSON by extension. Done
early, on the branch `merge/M6b` from `merge-0.7` after M1, in parallel with
the engine tasks. Python only: no Rust, schema, binding or CLI change.

## What changed

- `python/crn_craft/report.py` (new): the markdown and JSON writers,
  `save(result, path, *, paths)` and `load_report(path)`.
  - Formats by extension, case-insensitive: `.md`, `.markdown`, `.json`;
    any other extension raises `ValueError` naming the three.
  - Atomic write: a temporary file `.<name>.<random>.tmp` in the target
    directory (created with `O_EXCL` and mode `0o666`, so the umask applies
    as for `open`), written, flushed, `fsync`ed and renamed with
    `os.replace`. On any error the temporary file is removed and the error
    raised, so a missing directory or a target that is a directory raises
    `OSError` and leaves nothing behind. An existing file is replaced.
  - Deterministic bytes: no timestamp, host name or timing; numbers through
    `repr`; JSON strings ASCII-escaped; fixed key order; LF endings; one
    trailing newline.
  - `load_report` parses strictly (duplicate keys and `NaN`/`Infinity`
    literals raise `ValueError`) and checks `schema ==
    "crn-craft.report/1.0"`.
- `python/crn_craft/rules.py`:
  - `Result.model` and `Result.request`: the documents that ran, parsed anew
    on each access (so a caller cannot alter what `save` writes). `None` for
    a `Result` built directly from a metadata document and columns.
  - `Result.save(path, *, paths=True)`, returning `pathlib.Path(path)`; the
    docstring states the relation to archives.
  - `rules.load_report`, also exported as `crn_craft.load_report`
    (`python/crn_craft/__init__.py`).
  - `PreparedEnsemble` keeps `_ran = (model bytes, request bytes, bound)`
    and gives it to every result it creates (`_result`), on the interpreted
    path, with `simulate(compiled=...)` and with an attached loop.
    `CompiledLoop.verify(run)` (`_rules_loop.py`) now returns
    `run._result(...)`, so its result knows the documents too.
- `validation/check_rules_python.py`, `validation/check_installed.py`: the
  tests below. `validation/fixtures/sir-report.md` (new): the expected SIR
  report.
- `.github/workflows/python-release.yml`: `validation/fixtures/**` added to
  the pull-request path filter, since `check_installed.py` now reads the
  fixture.
- `docs/campaign/merge-0.7/M6b-manual-draft.md`: the manual section "Saving
  a report" with the SIR example and its real output, the JSON layout
  table, the reference rows and the R line, for M10b (I2 freezes the manual
  until M10).

### The documents of each path

Read before writing the code:

| Path | Documents today | What results report |
| --- | --- | --- |
| `prepare(model, request)` and `ensemble(...)` | `PreparedEnsemble._documents`, the exact `(model, request)` bytes | those bytes, parsed |
| `simulate(compiled=loop)`, attached loop, `loop.verify(run)` | the loop runs `_documents[1]` against its own model bytes, which the request's hash pins to `_documents[0]` | the same exact bytes |
| `bind_parameters(overrides)` | none in Python: the engine (`Prepared::bind_parameters`) clones the parsed model, sets the parameters, serializes it with serde, and binds a clone of the request to the hash of those bytes | reconstructed, see below |

For a candidate the engine's bytes are not reachable from Python without a
binding change. The candidate's documents are therefore rebuilt in Python
from the template's: the model with the overrides applied and all
parameters as floats sorted by name (the engine's `BTreeMap<String, f64>`),
and the template's request with `model.sha256` replaced, on access, by the
result's `model_sha256` (the hash of the bytes the engine ran). A candidate
of a candidate applies its overrides on top. Nothing that runs changes: the
native call and its arguments are the same as before (`dict(overrides)` is
built once and passed on), and the candidate's `_documents` stays `None`,
so compiled loops still refuse candidates exactly as before. The rebuilt
documents equal what the engine ran apart from fields that serde writes
with default values (for example `"attributes": {}` on products). The test
below checks that they run to the candidate's numbers: it prepares them
again and compares sampler, records, rule events, statistics and paths
exactly.

## Markdown layout decisions

The content rules of the plan are followed: sections, their order, their
lines, the observable order and the number rules. The line layout differs
in one point, decided by the coordinator after the first version (see
"Fenced sections" below). Where the plan left a choice open, this is what
the code does:

- **Withheld statistics.** `## Summary statistics` stays as a heading,
  followed by the single line `statistics withheld (<reason>)`; nothing
  else of the section. The reason is `ensemble termination: <termination>`
  (from `metadata["termination"]`, for example `event_limit`), or `no
  statistics requested` for a completed run whose request asked for none.
- **Method line.** `method = <kind>` followed by `, <option> = <value>` for
  each option in request order, for example `method = hybrid, rtol = 1e-06,
  atol = 1e-09`. The request's method is used, else the result's; a request
  without a method (rules-run/0.1, the tree method) prints `method = tree`.
- **Species types** (models with attributes): `mtDNA: clone, s, wt`; a type
  without attributes is its name alone.
- **Start values** of attributed entries: `mtDNA {clone=0, s=sw, wt=1} =
  w0`, attributes in the type's declaration order, values and counts as
  written (expressions stay expressions). Several entries of an
  attribute-free type are summed when all counts are numbers, else joined
  with ` + `.
- **Reactions.** A rule is written as a reaction line only when it has no
  fields besides `name`, `reactants`, `rate`/`propensity` and `products`,
  its patterns only `label`, `type` (and a future `count`, repeated), its
  products only `type` (and an empty `attributes`), and every type involved
  is attribute-free. Otherwise `<name>: <compact JSON>`. So a rule over an
  attributed type is JSON even without `where`, `factor` or `copy`.
- **Continuous variables.** One line per variable: `x: initial = 0,
  derivative = 2, lower_bound = 0, on_bound = clamp`, then `independent =
  true` and any other field, in that order.
- **Two sections the plan does not list**, written only for models that
  need them, so builder models and the SIR example are unaffected:
  `## Observables` (each observable as `name: <compact JSON>`) when some
  observable is not a plain count of the type it is named after, because
  otherwise the trajectories of a model such as the three-class model have
  no definition in the report; and `## ATP coupling` (the `atp` block as
  compact JSON) for models with one.
- Trajectories and statistics use the plan's observable order; ids are
  sorted ascending even if a result lists kept paths in another order.

### Fenced sections (decision after the first version)

The plan's layout writes the lines of a section directly under its heading,
without blank lines. A markdown renderer (GitHub, for example) joins such
lines into one paragraph, so the first version rendered each section as a
single run-on line. Decision: a markdown file must render as intended, and
a fenced code block keeps the user's line layout intact. So the data lines
of every section are in a fenced block: a line with three backticks and
`text` before them, a line with three backticks after them. Headings stay
outside the fences; there is a blank line before each heading (except the
first), after it, and around each fence. `Name: <model name>` stays a plain
paragraph under `# Model`, and `# Results` has no lines of its own.

Fenced: Species types, Reactions, Start values, Parameter values,
Observables, Stop condition, Continuous variables, ATP coupling,
Simulation, Terminations, rule_events, Single trajectories, the lines of
Summary statistics itself (`samples = N`, or the `statistics withheld
(...)` line) and each statistics subsection. The coordinator's list did
not name Species types and the Summary statistics lines; they are fenced
too under its rule "the data lines of every section". Content, order,
number rules and deterministic bytes are unchanged; the file still ends
with one newline (after the closing fence of the last section).

No data line can contain a backtick (names are identifiers, and the engine
has validated every expression before the run), so three backticks always
close the fence.

Tests first: the markdown parser of `check_rules_python.py` now checks the
layout (blank lines around headings and fences, exactly one ```` ```text ````
fence per section with data, no text outside a fence except the `Name:`
paragraph), the SIR model-section check uses the fenced text, and
`validation/fixtures/sir-report.md` was converted from the committed
version by a separate throwaway script (wrap each section's lines in a
fence, blank lines around). Both scripts failed against the unfenced
writer (`blank lines around the heading '## Reactions'`; `the SIR report
differs from .../sir-report.md`), and passed after the writer change; the
writer reproduced the converted fixture byte for byte. Five layout
mutations (an unfenced Species types section, no blank line before a
fence, none after a fence, two fences in one section, an unfenced
`samples` line) were each caught.

## JSON layout decisions

- Keys and nesting as in the plan; objects one member per line with two
  spaces of indentation, lists of numbers inline.
- Numbers in `result` follow the markdown number rule (integral floats
  below 2^53 without `.0`); `model` and `request` use the standard JSON
  forms, so `load_report(p)["model"] == result.model` and
  `["request"] == result.request` hold exactly, types included.
- Non-finite values anywhere are the strings `"nan"`, `"inf"`, `"-inf"`.
- `statistics.quantiles` is `[]` when the result has none; a statistic the
  result does not hold (`mean_standard_error` without `variance`) is
  `null`.

## Quantiles (M4e) and the fixture

Quantile statistics do not exist yet. The writers read
`result.statistics["quantiles"]` (a dict from probability to a DataFrame
shaped like `statistics["mean"]`) when present and ignore
`statistics["median"]` (it is `quantiles[0.5]`). Markdown titles: `###
median` for 0.5, else `### <100 p>% percentile` with the shortest exact
decimal (`Decimal(repr(p)) * 100`, normalized: 2.5, 25, 97.5), after
`### sample_variance`, ascending. Tested with results whose statistics the
test extends by hand with `np.quantile(result.paths.values, p, axis=0)`.

`validation/fixtures/sir-report.md` was created **without quantiles**.
**M4e must regenerate it** with quantiles 0.025, 0.5 and 0.975 (add
`quantiles=(0.025, 0.5, 0.975)` to the `rules.ensemble` call in
`check_installed.check_report`, or switch it to the top-level API whose
default these are, after M6a) and review the diff: it must add exactly the
three subsections `### 2.5% percentile`, `### median`, `### 97.5%
percentile`, each with its fenced lines, after `### sample_variance` (and
new hashes if the request version changes).

The fixture was written test-first and independently of the
implementation: a throwaway script built the expected text from the plan's
rules and the result's raw metadata (its own number formatting, its own
section assembly). The implementation, written afterwards, reproduced it
byte for byte on the first run. After the fence decision the fixture was
converted by another separate script, and the changed writer again
reproduced it byte for byte.

`check_installed.py` compares the report with the fixture byte for byte
except the one line `crn-craft = <version>`, which must name the installed
version. Without that substitution, every version bump (0.7.0 in M11) would
need a fixture edit in the release-prep commit.

## Tests written first, and how they failed first

Both scripts were extended before `report.py` existed and run against the
unchanged package (venv `/home/claude/m07/venv-m6b`, wheel of the M1 state):

- `check_rules_python.py` failed at the first new check:
  `AttributeError: 'Result' object has no attribute 'model'`.
- `check_installed.py` failed with `FileNotFoundError` for the fixture, and,
  once the independently derived fixture was in place,
  `AttributeError: 'Result' object has no attribute 'save'`.

New checks:

- `check_installed.py` (`rules_result_report`): the SIR model of the plan
  (10 replicates, seed 1, every path kept, `static_direct`, stop `I == 0`)
  saved as markdown equals `validation/fixtures/sir-report.md`; saving
  twice gives identical bytes; `save` returns the `pathlib.Path`; the JSON
  report loads with `crn_craft.load_report` and holds the model and request
  that ran.
- `check_rules_python.py` ("result model and request of every path",
  "markdown and JSON reports"):
  - `result.model`/`request` from `ensemble` equal the builder's document
    and the request bytes; a fresh dictionary per access; `None` on a bare
    `Result`, whose `save` raises `ValueError` and writes nothing.
  - `bind_parameters`: the bound values (and those of a candidate of a
    candidate) appear, the request carries the candidate's `model_sha256`,
    and preparing the reported documents again reproduces the candidate's
    sampler, records, rule events, statistics and paths exactly.
  - The compiled path, with a stand-in `CompiledLoop` whose `execute` runs
    the interpreter: `simulate(compiled=loop)` and `loop.verify(run)` give
    results that know their documents.
  - JSON round trip: header, `model == result.model`, `request ==
    result.request`, provenance, report times, observables in result order,
    terminations against `records`, mean rule events against the exact
    integer sums of `rule_events`, every statistic against its DataFrame,
    quantiles in ascending order, paths against `result.paths` (NaN where
    NaN), all compared exactly; the file is strict JSON (no `NaN` literal,
    no duplicate key) with one trailing newline; integral result numbers
    without `.0`, documents in their standard form.
  - Markdown: section order, simulation lines, terminations, mean rule
    events, every trajectory line and every statistics line parsed back and
    compared exactly, observable order (`S, I, R` for the SIR builder model;
    result order for the three-class model; `A, x` for a hybrid model), the
    SIR model section text.
  - `paths=False` gives `"paths": null` and no trajectory section;
    `paths=1` raises `TypeError`.
  - Withheld statistics (an event-limited run): `"statistics": null` and
    the markdown line naming the termination; `nan` in partial paths.
  - Quantiles 0.975, 0.025, 0.5, 0.25 added by hand: titles `2.5%`, `25%`,
    `median`, `97.5%` in that order; JSON `p` ascending.
  - The three-class mtDNA model: `## Species types`, rules as name and
    compact JSON, attributed start values, parameters, trajectories of two
    kept cells. A hybrid model: `## Continuous variables`, `0 -> A, x`, no
    parameter section, the method line with tolerances.
  - Extensions: `.MD`, `.MARKDOWN` (replacing an existing file) and `.Json`
    work; `.txt`, none, `.md.bak`, `.jsonl` raise `ValueError` naming the
    supported ones and write nothing.
  - A failed write (directory missing; target is a directory) raises
    `OSError` and leaves no temporary file.
  - Saving changes neither `metadata` nor any frame of the result.
  - `load_report` refuses another schema, a list and a document without
    `schema`; `crn_craft.load_report is rules.load_report`.

To see that the new checks discriminate, each of these mutations of the
implementation was run against them (in a scratch process, not committed),
and each was caught: result order instead of declaration order; a write
without temporary-file cleanup; a wrong percentile title; `NaN` literals in
JSON; an empty withheld reason; integral floats with `.0`; a dropped
trajectory; bound candidates reporting the template's documents;
`loop.verify` and `simulate(compiled=...)` without documents. (The first
version of the withheld check accepted an empty reason; it was tightened to
require the termination name.)

## Checks

Run with the task venv (`/home/claude/m07/venv-m6b`, the M1 wheel plus the
changed Python files copied into `site-packages/crn_craft/`):

| Check | Result |
| --- | --- |
| `python validation/check_rules_python.py` and `python -O ...` (from `validation/`) | passed; reports identical |
| `python -I .../validation/check_installed.py` and `-I -O`, run from a temporary directory outside the source tree | passed, 17 checks (new: `rules_result_report`); reports identical |
| `python validation/check_manual.py --binary /home/claude/fb2/engine/target/release/crn-craft --require-cargo` with `CRN_CRAFT_SOURCE=/home/claude/m07/wt-m6b`, normally and under `-O` (I2) | passed, 22 of 22 blocks match the manual, including both compiled-loop blocks; reports identical; 147 s and 103 s |
| `python validation/check_compiled_loop.py --require-cargo` with `CRN_CRAFT_SOURCE` set, normally and under `-O` (exercises the changed `CompiledLoop.verify`) | passed, 50 rows; reports identical; the loop program came from the build cache (the crate sources equal M1's) |
| `python -m ruff format --check --exclude '*.md' .` and `python -m ruff check .` (ruff 0.11.13, as pinned in `validation/requirements.txt`) | 112 files formatted; all checks passed |
| Python 3.10 grammar (`ast.parse(..., feature_version=(3, 10))` under Python 3.11) for every changed Python file | parsed |
| `cargo fmt --all -- --check` | passed |

Not run: `cargo clippy`, `cargo test --release --workspace --locked`
(I1, I3) and `cargo doc`. The task changes no file under `crates/`,
`bindings/`, `schemas/`, `tests/` and no Cargo file
(`git diff --stat 23c7456 HEAD -- crates bindings tests schemas Cargo.toml
Cargo.lock` is empty), so they cannot change; running them would have
built a second full `target/` on the shared machine. The coordinator's
runs on `merge-0.7` cover them.

After the fence change (markdown section helper in `report.py`, the
`save` docstring, tests, fixture, docs): `check_rules_python.py` and
`check_installed.py` (outside the tree, `-I`) were rerun normally and under
`-O` and passed with identical reports; ruff format and check are clean.
`check_manual.py` and `check_compiled_loop.py` were not rerun: they do not
call `save`, and nothing else changed.

Measured (no target declared, before the fence change, which adds four
short lines per section): `Result.save` of a three-class CME ensemble with
4,000 cells x 37 report times x 6 observables, every path kept, median of
three under the CPU lock:

| Call | Time | File size |
| --- | --- | --- |
| `save("r.md")` | 0.633 s | 3,463,221 bytes |
| `save("r.json")` | 1.010 s | 3,213,890 bytes |
| `save("r.md", paths=False)` | 0.003 s | 16,142 bytes |
| `save("r.json", paths=False)` | 0.003 s | 19,826 bytes |

## Limitations

- Bound candidates' documents are rebuilt in Python (see above), not read
  from the engine. If exact bytes are wanted, a binding accessor on the
  native prepared object (re-serializing its source model and request)
  would provide them; M6a, which changes the bindings anyway, could add
  it. Exact equality with the engine's bytes is not
  claimed, only that the documents run to the same numbers.
- The manual draft shows the report in a four-backtick fence because the
  report contains fences; `check_manual.py` must learn longer fences before
  M10b adds that block (noted in the draft).
- The markdown report shows the model's `description` nowhere (the JSON
  report holds the full model).
- The fixture depends on the request version and hashes that
  `crn_craft.rules` writes; M4e (quantiles) and any change of the request
  `rules.ensemble` builds regenerate it.

## Deviations from the plan

- The fixture has no quantiles yet (coordinator instruction); M4e
  regenerates it.
- `check_installed.py` compares the fixture byte for byte except the
  `crn-craft = <version>` line, which is compared with the installed
  version.
- Bound candidates' documents are reconstructed (above).
- Two conditional markdown sections the plan does not list (`## Observables`,
  `## ATP coupling`), and the interpretation of withheld statistics as the
  heading plus the one line.
- Fenced sections (coordinator decision): the data lines of every section
  are in a ```` ```text ```` fence, with blank lines before each heading and
  around each fence, instead of the plan's bare lines. Reason: a markdown
  file must render as intended, and fences keep the user's line layout
  intact; the plan's layout rendered each section as one paragraph.
- The manual text is in `M6b-manual-draft.md` for M10b, as the plan
  schedules the docs in M10b and I2 freezes the manual.

## Migration notes

Additive only; no existing behavior or result byte changes.

- New: `Result.model`, `Result.request`, `Result.save(path, *,
  paths=True)`, `crn_craft.rules.load_report(path)` and
  `crn_craft.load_report(path)`; the JSON report format
  `crn-craft.report/1.0` (Python-level, no schema file).
- `CompiledLoop.verify(run)` returns the same result as before, which now
  also knows its documents.
- Internal: `PreparedEnsemble.__init__` gained a keyword-only `ran`.
