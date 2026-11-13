# M6b manual draft: "Saving a report"

User-facing text for the user manual (`docs/mtdna-user-manual.md`), written
in task M6b of [the plan](../../history/plan-0.7.0-engine-merge.md) for task M10b to
merge. Invariant I2 freezes the manual until M10, so the text is kept here.
Both blocks and their outputs were produced with the M6b implementation on
the 0.6.0 wheel through `crn_craft.rules`; the first output equals
`validation/fixtures/sir-report.md` byte for byte.

## Notes for M10b

- Rewrite the blocks on the top-level API (`import crn_craft as crn`,
  `crn.Model`, `crn.ensemble`) and rerun them. Through the top-level API the
  model is `crn-craft.model/1.0` and the request `crn-craft.run/1.0` (M6a),
  and quantiles 0.025, 0.5 and 0.975 are on by default (M4e). So the hashes,
  the request schema and the statistics change, and the report gains the
  subsections `### 2.5% percentile`, `### median` and `### 97.5% percentile`
  after `### sample_variance`. Paste the new output from the run; it must
  equal the regenerated `validation/fixtures/sir-report.md` (M4e) when the
  arguments are the same.
- The report prints the line `crn-craft = <version>`. `check_manual.py`
  compares it literally, so the release-prep commit that changes the
  version (M11) must update this output too, or the block should print the
  report without that line. `check_installed.py` already substitutes the
  installed version for that one line of the fixture.
- The printed report contains ```` ```text ```` fences, so the manual shows
  it in a four-backtick fence (````` ````text `````), as below. Today
  `check_manual.py` recognizes only three-backtick fences
  (`FENCE = re.compile(r"^```(\w*)\s*$")`) and ends a block at any line
  that starts with three backticks, so it would read the report's own
  fences as blocks. M10b must teach it fences of any length (an opening
  run of three or more backticks, closed by a run of the same length)
  before adding this block, or print the report without its fences.
- Suggested place: a new subsection at the end of section 3 ("A first
  ensemble") or a short section of its own after section 17 ("Checking
  results"), plus the reference rows below in section 18.
- `docs/using-from-r.md`: one line, for example after the statistics
  example: `result$save("sir.md")  # one markdown report: model, run, results`.

## Text for the manual

### Saving a report

`result.save(path)` writes the model, the request and the results to one
file. The extension chooses the format: `.md` or `.markdown` for a
markdown report to read, `.json` for the same content as JSON
(`crn-craft.report/1.0`). It returns the path as a `pathlib.Path`.

```python
"""An SIR ensemble saved as one markdown report."""

from pathlib import Path
import tempfile

from crn_craft import rules

sir = rules.Model("sir")
sir.add_species({"S": 999, "I": 1, "R": 0})
sir.add_parameter({"beta": 2e-4, "gamma": 0.01})
sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
sir.set_stop("I == 0")
result = rules.ensemble(
    sir,
    end=250.0,
    report_times=[0, 50, 100, 150, 200, 250],
    replicates=10,
    seed=1,
    keep_paths="all",
    method="static_direct",
).simulate()
with tempfile.TemporaryDirectory() as folder:
    path = result.save(Path(folder) / "sir.md")  # or result.save("sir.md")
    print(path.read_text(), end="")
```

````text
# Model

Name: sir

## Reactions

```text
S + I -> I + I, beta * S * I
I -> R, gamma * I
```

## Start values

```text
S = 999
I = 1
R = 0
```

## Parameter values

```text
beta = 0.0002
gamma = 0.01
```

## Stop condition

```text
I == 0
```

# Simulation

```text
seed = 1
end = 250
replicates = 10
report_times = 0, 50, 100, 150, 200, 250
method = static_direct
sampler = chacha8-u53-ziggurat/rules-static-direct/1
crn-craft = 0.6.0
model_sha256 = 3a56995648ee284605c18b4e13aaed5afcc8abf5136fdd9c475c2fa403195575
run_sha256 = 837c92783f25c7a1202fb90c0d1ef8014dae7d5657a4117417125f48913d968e
```

# Results

## Terminations

```text
completed = 10
```

## rule_events (mean per replicate)

```text
infection = 999
recovery = 877.7
```

## Single trajectories for S, I, R

```text
S000: 999, 221, 0, 0, 0, 0
I000: 1, 705, 583, 357, 218, 135
R000: 0, 74, 417, 643, 782, 865
S001: 999, 26, 0, 0, 0, 0
I001: 1, 825, 511, 298, 176, 114
R001: 0, 149, 489, 702, 824, 886
S002: 999, 32, 0, 0, 0, 0
I002: 1, 802, 500, 298, 185, 112
R002: 0, 166, 500, 702, 815, 888
S003: 999, 65, 0, 0, 0, 0
I003: 1, 830, 538, 321, 195, 106
R003: 0, 105, 462, 679, 805, 894
S004: 999, 104, 0, 0, 0, 0
I004: 1, 788, 566, 338, 210, 129
R004: 0, 108, 434, 662, 790, 871
S005: 999, 45, 0, 0, 0, 0
I005: 1, 797, 513, 313, 197, 121
R005: 0, 158, 487, 687, 803, 879
S006: 999, 448, 0, 0, 0, 0
I006: 1, 513, 585, 347, 204, 119
R006: 0, 39, 415, 653, 796, 881
S007: 999, 141, 0, 0, 0, 0
I007: 1, 750, 558, 343, 209, 135
R007: 0, 109, 442, 657, 791, 865
S008: 999, 102, 0, 0, 0, 0
I008: 1, 776, 549, 326, 195, 124
R008: 0, 122, 451, 674, 805, 876
S009: 999, 63, 0, 0, 0, 0
I009: 1, 795, 520, 321, 204, 128
R009: 0, 142, 480, 679, 796, 872
```

## Summary statistics

```text
samples = 10
```

### mean

```text
S: 999, 124.69999999999999, 0, 0, 0, 0
I: 1, 758.1, 542.3000000000001, 326.2, 199.29999999999998, 122.30000000000001
R: 0, 117.20000000000002, 457.7, 673.8, 800.7, 877.7
```

### mean_standard_error

```text
S: 0, 40.423329347736264, 0, 0, 0, 0
I: 0, 29.556142433605157, 9.70458081984425, 6.343500610861484, 3.932909581693657, 3.0841890704401083
R: 0, 12.462209541917785, 9.70458081984425, 6.343500610861485, 3.932909581693662, 3.08418907044011
```

### sample_variance

```text
S: 0, 16340.455555555556, 0, 0, 0, 0
I: 0, 8735.655555555553, 941.7888888888888, 402.40000000000015, 154.67777777777772, 95.12222222222219
R: 0, 1553.0666666666668, 941.7888888888888, 402.4000000000002, 154.67777777777815, 95.12222222222229
```
````

The report has three parts. Each section's lines are in a fenced `text`
block, so a markdown viewer such as GitHub shows them one per line in a
monospace font, and an editor shows the same lines:

- **Model**: one line per reaction (`reactants -> products, rate`; a rate
  marked `(mass action)` is a constant per reactant combination, anything
  else is the whole propensity), the start values of every species, the
  parameter values sorted by name, and the stop condition and continuous
  variables when the model has them. A model with attributes also lists its
  species types, and a rule that uses attributes, `where`, `factor` or
  `copy` is printed as its name and its JSON. Observables that are not
  plain species counts are listed with their definitions.
- **Simulation**: seed, horizon, replicates, report times, the method with
  its options, the sampler, the package version and the SHA-256 hashes of
  the model and request that ran.
- **Results**: how many cells ended in each termination, the mean number
  of events per rule, the kept trajectories (one line per cell and
  observable, `S003` is S in cell 3; values after a budget stop are `nan`)
  and the summary statistics with one line per observable. Observables
  appear in the order the species were declared, then continuous variables
  (models whose observables are not all species counts use result order,
  sorted by name). When the ensemble withholds its statistics (section 14),
  the report says `statistics withheld (ensemble termination: ...)`.

Numbers are written in their shortest round-trip form, so reading them back
gives the same floats; integral values have no `.0`. The same result always
gives the same bytes: a report has no timestamp, host name or timing.
`result.save(path, paths=False)` leaves the trajectories out, which keeps
the file small for large ensembles. An existing file is replaced. The report
is written to a temporary file next to it and renamed, so a failed write
raises `OSError` and leaves no partial file.

The JSON report holds the same information for programs.
`crn_craft.load_report(path)` (also `rules.load_report`) reads it back and
checks its schema:

```python
"""The same report as JSON, read back with load_report."""

from pathlib import Path
import tempfile

import crn_craft
from crn_craft import rules

sir = rules.Model("sir")
sir.add_species({"S": 999, "I": 1, "R": 0})
sir.add_parameter({"beta": 2e-4, "gamma": 0.01})
sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
sir.set_stop("I == 0")
result = rules.ensemble(
    sir,
    end=250.0,
    report_times=[0, 50, 100, 150, 200, 250],
    replicates=10,
    seed=1,
    keep_paths="all",
    method="static_direct",
).simulate()
with tempfile.TemporaryDirectory() as folder:
    report = crn_craft.load_report(result.save(Path(folder) / "sir.json"))
print(list(report))
print(report["model"] == sir.to_dict(), report["request"] == result.request)
body = report["result"]
print(body["observables"], body["terminations"], body["rule_events_mean"])
print("mean I:", body["statistics"]["mean"]["I"])
print("replicate 3, I:", body["paths"]["values"]["I"][3])
```

```text
['schema', 'crn_craft_version', 'model', 'request', 'result']
True True
['I', 'R', 'S'] {'completed': 10} {'infection': 999, 'recovery': 877.7}
mean I: [1, 758.1, 542.3000000000001, 326.2, 199.29999999999998, 122.30000000000001]
replicate 3, I: [1, 830, 538, 321, 195, 106]
```

| Key | Content |
| --- | --- |
| `schema`, `crn_craft_version` | `"crn-craft.report/1.0"` and the package version that wrote it |
| `model`, `request` | the model and request documents that ran, equal to `result.model` and `result.request` |
| `result.model_sha256`, `run_sha256`, `sampler` | provenance, as in `result.metadata` |
| `result.report_times`, `observables` | the report grid and the observables in result order (sorted names, then continuous variables) |
| `result.terminations` | `{termination: count}` for each termination that occurred |
| `result.rule_events_mean` | `{rule: mean events per cell}`, in rule order |
| `result.statistics` | `null` when withheld; else `samples`, then `mean`, `mean_standard_error`, `sample_variance` as `{observable: [value per report time]}` (`null` if not computed), and `quantiles`, a list of `{"p": p, "values": {observable: [...]}}` in ascending `p` (empty without quantiles) |
| `result.paths` | `null` without kept paths or with `paths=False`; else `replicates` (the kept cell ids) and `values`, `{observable: [[value per report time] per cell]}` |

Numbers in `result` follow the markdown rule (integral values without
`.0`); `model` and `request` are written as the documents were. Non-finite
values are the strings `"nan"`, `"inf"` and `"-inf"`, so the file is strict
JSON. No schema file is published for this Python-level format; this table
documents it.

A report is a readable record of one run. It is not an archive: an archive
is the sealed, verifiable record of a run for replay. To rerun a report's
model and request, set `request["model"]["sha256"]` to the SHA-256 of the
model bytes you pass (`hashlib.sha256(crn_craft.model_bytes(model))
.hexdigest()` for a dictionary) and prepare them again with `rules.prepare`.

`result.model` and `result.request` are the documents themselves, as new
dictionaries. For a candidate from `run.bind_parameters(...)` they hold the
bound parameter values, and `request["model"]["sha256"]` is the candidate's
`model_sha256`.

### Reference rows (section 18)

For the table "The result, `run.simulate()`":

| Attribute | Content |
| --- | --- |
| `model`, `request` | the model and request documents that ran, as new dictionaries; for a `bind_parameters` candidate with the bound values |
| `save(path, *, paths=True)` | writes the model, the request and the results to one markdown (`.md`, `.markdown`) or JSON (`.json`, `crn-craft.report/1.0`) file; returns the `pathlib.Path` |

And under the module functions: `crn_craft.load_report(path)` (also
`rules.load_report`) reads a JSON report and checks its schema.
