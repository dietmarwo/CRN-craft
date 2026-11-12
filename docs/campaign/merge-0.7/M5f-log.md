# M5f log: import and export

Task M5f of [the plan](../../plan-0.7.0-engine-merge.md): the importers
return converted `crn-craft.model/1.0` documents, and the GillesPy2 export
(`export.rs`) works on the `FlatView`: a constant rate exports as mass
action, anything else as a custom propensity string.

## What changed

- **`export.rs`** rewritten on the view (`FlatView::of_any_model`, so it
  reads `model/0.1` through the converter, the rules models and
  `model/1.0`; the view of `FlatView::with_every_type`, so a converted
  network keeps every species and reaction). One GillesPy2 species `S{i}`
  per species label, one reaction `R{c}` per channel:
  - a channel whose rule has a constant `rate` (reads no observable) and
    that takes at most two molecules is GillesPy2 mass action with constant
    rate × pattern `factor` (halved for a homodimer; `CRN-E051` when the
    halving is inexact), written exactly as 0.6.0 wrote it;
  - any other channel (an observable-dependent rate, a whole-rule
    `propensity`, more than two molecules) is a custom propensity string:
    the rule's expression with parameter values inlined, each observable as
    the weighted sum of its species, `^` as `**`, times the falling
    factorial n (n − 1) … / k! per reactant (`rate` rules) or times one
    availability factor `x/(x + 1e-300)` per molecule taken (whole-rule
    propensities: exactly 1 while the count suffices, else 0);
  - refused with `CRN-E050` (no partial source): comparisons, logic and
    functions in the expression (GillesPy2's Python solvers have no `exp`
    and its C++ solver no `min`/`max`/`abs`; its expression printer drops
    comparisons), observables reduced other than by `sum`, more than 64
    molecules of one species in a custom propensity, a negative constant
    rate, and a `stop` condition (GillesPy2 cannot stop a trajectory, as it
    could not keep the old population limits).
  - The report stays `crn-craft.export/0.1` (same shape). Its diagnostics
    are the errors of reading and lowering the model plus `CRN-E050` and
    `CRN-E051`; model paths are the rules model's (`/rules/i/rate`, for a
    converted network the reaction's index). Reaction aliases are the rule
    names, `rule[j]` for the channels of a rule with several. A fifth note
    describes custom propensities when the source has one.
- **`rules/expr.rs`**: `Expr::to_python` writes an expression as Python
  source (every operation parenthesized, numbers in the shortest round-trip
  form) or names its first unsupported construct; `python_literal`.
  **`rules/flat_view.rs`**: `FlatView::python_rate(rule)` maps the slots
  (parameters inlined, observables as sums over the view's species).
- **Discovery**: the core `reference_exporters` entry now says
  `"subset": "mass_action_and_custom_propensities"`,
  `"maximum_mass_action_order": 2`, `"maximum_custom_order": 64`.
  `schemas/export.schema.json` regenerated (descriptions only).
- **`python/crn_craft/importers.py`**: `import_sbml` and `import_bng_net`
  build and validate the network model exactly as before, then convert it
  with `crn_craft.convert` and return a `crn-craft.import/1.0` report whose
  `model` is the converted `crn-craft.model/1.0` document as a dictionary
  (the converter's canonical bytes parsed in their key order). The
  provenance fields (`source_format`, `source_sha256`, conventions,
  `source_parameters`, `species_patterns`, `groups`, `parameter_policy`)
  are unchanged. A model the converter refuses (a network without
  reactions) raises `ImportFailure` (`CRN-E090`) with its diagnostics. The
  `crn-craft-import` script writes the new report.
- **Validation**: `validation/merge/check_export_references.py` (export
  equivalence against the 0.6.0 binary, declared first) and its report;
  `check_export.py`, `check_importers.py`, `check_contracts.py` and
  `check_installed.py` migrated in place (the `export` command and the
  importers themselves changed). CI keeps running `check_export.py` in the
  `reference-export` job with `validation/reference-requirements.txt`; no
  workflow change was needed.

## Return value for M6a

`read_sbml(path, convention=...)` and `read_bngl(...)` (M6a) can wrap
`import_sbml(raw, convention=...)["model"]` in a `Model`: the `model`
entry is a complete `crn-craft.model/1.0` dictionary, and the rest of the
report is the provenance M6a may keep as model metadata. The importers keep
their keyword arguments and their `ImportFailure` (a `ValueError`).

## Tests written first

Run against the 0.6.0 export (logs in
`/home/claude/m07/m5bf/first-fail-export.log`):

- `crates/crn-craft/tests/export.rs` (8 tests): 6 failed. They check that a
  network model keeps its 0.6.0 artifact byte for byte (the expected source
  was printed by the 0.6.0 binary) and that its converted model exports the
  same source; the custom propensities of Michaelis–Menten, Hill repression
  and third-order mass action; an attributed `model/1.0` with labels,
  `factor`, weighted observables, a rule with two channels and a
  whole-rule propensity of count 2; every `CRN-E050` refusal (and that 64
  molecules are accepted); converter refusals and `CRN-R017`; the `CRN-E051`
  path. The two that passed are ports that hold for both exporters
  (homodimer halving and aliases; contract and count boundary).
- `tests/contracts_v1.rs`: export of `model/1.0` and rules models succeeds
  (failed first).
- `crates/crn-craft-cli/tests/cli.rs`: Michaelis–Menten exports, the ATP
  model is refused with `CRN-R017` (failed first); `export` left the list of
  commands that refuse `model/1.0`.

## Measured

Declared in `validation/merge/RESULTS.md` before the runs; rows there.

| Check | Result |
| --- | --- |
| Export equivalence (twelve models and the `check_export.py` fixture) | 9 artifacts byte-identical to 0.6.0; 40 custom propensities × 4 states = 160 evaluations equal to the engine bit for bit; **0 failures**; `-O` identical |
| `check_export.py` (GillesPy2 1.8.3) | 9 of 9; 15,756 rate comparisons, largest relative error 1.35e-16; `-O` identical |
| `check_importers.py` | 49 of 49 (one new rejection); `-O` identical |
| `check_contracts.py` | 101 of 101, every example exports; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 17 checks pass; `-O` identical |

## Deviations from the plan

- The importers return the report with the converted model, not a bare
  model dictionary or `Model` objects; see "Return value for M6a".
- A constant rate of a channel that takes more than two molecules is a
  custom propensity: GillesPy2's mass action stops at order two.
- Refusals listed above (`stop`, functions, comparisons, non-sum
  observables, order above 64, negative constants) are not in the plan;
  they are what GillesPy2 cannot evaluate in both solvers.
- The export report no longer repeats the structural warnings
  (`CRN-W001`–`W003`) of the network validation; `analyze` reports them.

## Limitations

- Custom propensities are evaluated by GillesPy2 in Python or C++ floating
  point, not by the engine; a division by a zero count (for example the
  converted Hill activation `vmax / (1 + (k / R)^n)` at R = 0) raises
  `ZeroDivisionError` when Python evaluates it with integer or float
  operands. The tested models do not divide by a count.
- GillesPy2 1.8.3 parses custom propensities with a deprecated `ast`
  visitor; under `python -W error::DeprecationWarning` building such a
  model fails inside GillesPy2.
- A network without reactions can no longer be imported (the converter
  refuses it).

## Migration notes

- `crn-craft export MODEL gillespy2` reads every model version; Hill and
  Michaelis–Menten laws, reactions of order above two and rules models
  export (custom propensities). `CRN-E051` is reported at `/rules/i/rate`.
- `crn_craft.importers.import_sbml` and `import_bng_net` (and
  `crn-craft-import`) return `crn-craft.import/1.0` with a
  `crn-craft.model/1.0` model; read rate constants from
  `report["model"]["rules"][i]["rate"]` (a string such as `"0.1"`).
- Discovery's `reference_exporters` entry changed (above).
