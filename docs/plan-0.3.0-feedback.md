# CRN-craft 0.3.0: plan for the first user feedback

Six findings came from a researcher who tried the 0.2.0 Python interface on mtDNA mutant models. This plan fixes every one in small increments.

- Each increment is a commit with its own tests. A fix starts with a test that fails before it.
- All increments land in one pull request.
- The result is release 0.3.0.

Owner decisions (2026-11-06):

- `mtdna_direct` stays as a documented opt-in.
- The ATP coupling is generalized in 0.3.0.
- pandas becomes a required dependency of the Python package.

## Findings and decisions

| # | Finding | Cause | Decision |
| --- | --- | --- | --- |
| 1 | `mtdna_direct` looks hard-coded for one model | It is a fast path certified at preparation for the three-class structure. It falls back to generic Static Direct with `CRN-R051`, and is about 3.9× faster for that model. | Keep it as an opt-in. The docs lead with the generic methods and state plainly what the fast path requires. |
| 2 | Need values for every replicate; DataFrames instead of ndarrays | Per-replicate series exist (`keep_paths`), but only as raw metadata with no accessor or documentation. Results are dicts of ndarrays. | Results become pandas objects: statistics, records, rule events, paths, trackers and event logs. Add `keep_paths="all"` and `to_xarray()`. |
| 3 | `event_log={"max_entries": 50, …}` fails with `CRN-R020: observation reservation exceeds budget` | The reservation is a worst-case bound. It uses the global product limit instead of the model's, and the default output budget is 64 MiB. The message does not say what is needed. | Base the bound on the model's actual product and attribute counts. Make the error state the bytes needed, the budget and the fix. Check the default budget. |
| 4 | Manual chapter 7: `names` is undefined in `watch` | The manual shows an excerpt of `04_event_log.py` without its first lines. | Make every manual block run on its own, and run all of them in CI. Expose rule names on prepared runs and results, so callbacks need no lookup table. |
| 5 | Where is dA/dt = f·wt − v1·A − v2·N in the ATP model file? | It is a fixed `atp` block with that ODE form hard-coded. | Add generic continuous variables: ODEs written as expressions that rule rates may read, simulated by integrating the hazards. The ATP model becomes an ordinary model file; the legacy block stays readable. |
| 6 | Building models in JSON is cumbersome; wanted is rebop-style `add_reaction(rate, reactants, products)` | There is no builder API. | Add a Python `Model` builder. A number is a mass-action constant; a string is the full propensity. Observables are derived from species. The engine gains a whole-rule propensity form. |

## Increments

### I1. Self-contained manual examples and rule names (finding 4)

- **Failing-first test:** execute every `python` block of `docs/mtdna-user-manual.md` in a fresh interpreter, in order, from the repository root, with timeouts. Today the `watch` block fails with `NameError`.
- **Fix:** every block defines what it uses.
- **Rule names:**
  - add `PreparedEnsemble.rule_names` and `Result.rule_names`;
  - callback events keep the rule index; the manual shows `run.rule_names[event["rule"]]`.
- **CI:** add the manual check to the Python job.

### I2. Event-log budget (finding 3)

- **Failing-first test:** reproduce the reported request with the Python default budget: the three-class model, `event_log={"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]}`, at the replicate counts the manual uses. It fails today with `CRN-R020`.
- **Tighter bound:** compute the per-entry bound from the model's largest product count and attribute count instead of `MAX_PRODUCTS`, so it stays an upper bound on what the result can hold. Test that the bound still covers measured output for a log filled to `max_entries`.
- **Actionable error:** the `CRN-R020` message names the reservation (replicates × entries × bytes per entry, plus trackers), the budget, and the smallest `max_output_bytes` that would pass.
- **Default budget:** check whether the Python default is too small for the manual's documented uses. Change it only with a stated reason.

### I3. pandas results (finding 2)

- **Dependency:** `pandas` is required in `pyproject.toml`; `xarray` is optional.
- **`Result.statistics`:** a dict of DataFrames, one per statistic: index = report time, columns = observables.
- **`Result.records`:** a DataFrame, one row per replicate, indexed by replicate id. `termination` is categorical with names.
- **`Result.rule_events`:** a DataFrame with rule-name columns.
- **`Result.paths`:** a DataFrame with a (replicate, time) MultiIndex and one column per observable.
  - `keep_paths="all"` or `True` keeps every replicate.
  - `Result.to_xarray()` returns a Dataset with dims (replicate, time).
- **`tracker_columns(i)` and `event_columns(id)`:** return DataFrames.
- **NumPy access:** stays available through `.to_numpy()`. Breaking changes get migration notes in the CHANGELOG.
- **Tests:**
  - shapes, index values and dtypes;
  - identity with the raw metadata;
  - `keep_paths="all"`;
  - `-O` behavior;
  - `check_installed.py`.

### I4. Whole-rule propensity and the `Model` builder (finding 6)

- **Engine** (`crn-craft.rules-model/0.3`): a rule may give `propensity` (the total rate expression) instead of `rate` (per reactant tuple).
  - Allowed only when the reactants carry no `where` filter and the expression reads no reactant label.
  - The propensity is the expression's value when every reactant type has enough molecules, else 0. It must be finite and nonnegative, or the replicate fails with a coded diagnostic.
  - Firing picks distinct reactant molecules uniformly.
- **Methods and mtDNA path:** all methods support it, or refuse with a coded diagnostic; `mtdna_direct` falls back with `CRN-R051`.
- **Schema:** regenerated, with a new contract version only for the model.
- **Python `crn_craft.rules.Model`:**
  - `Model(name)`, `add_parameter(name, value)`, `add_species(name, initial=0)`, `add_reaction(rate, reactants, products, name=None)`, `set_initial(...)`, `to_dict()` / `to_json()`;
  - `ensemble()` accepts a `Model`;
  - a numeric rate is mass action, with the engine's counting of identical reactants documented and compared with rebop's;
  - a string rate is the full propensity, as in rebop;
  - observables are all species;
  - species are attribute-free types (attributes stay available through JSON).
- **Tests:**
  - against analytic results: birth–death stationary mean and Michaelis–Menten mean flux;
  - builder JSON equals hand-written JSON;
  - a rebop-style SIR model gives the same mean within statistical tolerance as the equivalent mass-action JSON model.

### I5. Generic continuous variables (finding 5)

- **Model contract** (rules-model/0.3): an optional `continuous` list of variables. Each has:
  - `name`, `initial` (expression);
  - `derivative`: an expression over parameters, observables and continuous variables;
  - optional `lower_bound` with `on_bound`: `"stop"` or `"clamp"`.
- **Rates:** rule rates and propensities may read continuous variables.
- **Method `hybrid`:**
  - between events, integrate the ODE together with the cumulative total hazard using a fixed-tolerance adaptive embedded Runge–Kutta method (a published method, cited);
  - an event fires when the integrated hazard reaches an exponential threshold, located by dense output and root finding;
  - tolerances are declared in the request, recorded in the result, and the method is labeled as approximate at that tolerance;
  - results carry the continuous variables as extra observable columns.
- **ATP equivalence:** the ATP model is rewritten as an ordinary rules-model/0.3 file with continuous variable `ATP`. Validate it against the existing exact `thinning` method on the legacy file: takeover-time distribution and final ATP distribution, with predeclared two-sample tests and fixed seeds and sizes. Failed rows are retained.
- **Analytic checks:** a decaying variable driving a death rule (survival probability in closed form); a constant derivative (linear rate, time-inhomogeneous Poisson count); a bound stop.
- **Legacy files:** the legacy `atp` block and `thinning` keep working unchanged, so stored 0.2 models replay byte-identically.

### I6. Method guidance (finding 1)

- **Manual and guides:** the user manual, `docs/mtdna-rules.md`, `describe` and the Python docstrings lead with the generic methods (`static_direct`, `next_reaction`, `tree`, `tau`, `binomial`, `hybrid`).
- **Opt-in fast paths:** `mtdna_direct` and the Rust generated loops are described as opt-in. The text says what each requires, that it is certified per model at preparation, that it falls back otherwise, and what it gains.
- **Performance analysis:** it says which numbers depend on the fast path.

### I7. Documentation, examples and release preparation

- **User manual:**
  - a new chapter on the `Model` builder;
  - chapter 10 rewritten for generic continuous variables;
  - DataFrame results throughout;
  - all examples rerun and outputs refreshed.
- **Other docs:** README, `python/README.md`, `docs/rules-language.md` and `docs/mtdna-rules.md`.
- **Version 0.3.0:** bump the workspace, the bindings, `pyproject.toml` and the five lockfiles. Add a CHANGELOG entry with highlights and migration notes. Update `check_installed.py` for the builder, DataFrames, the event-log budget and `hybrid`.
- **Checks:**
  - fmt, clippy `-D warnings`, release tests;
  - `check_contracts.py` (normal and `-O`);
  - `check_rules.py`;
  - the manual examples;
  - ruff;
  - packaging dry runs.
- **Release:** the owner tags `v0.3.0` (RELEASING.md).

## Order and ownership

Two streams run in parallel and are then merged into one integration branch for one PR:

- **Engine:** I2 → I4 → I5, with `docs/rules-language.md`.
- **Python and docs:** I1 → I3 → I6.

I7 runs on the integration branch after both streams. No increment changes the results of an existing request: stored 0.1–0.4 rules runs replay byte-identically.

## Status (2026-11-06)

All seven increments are done on the branch `release/0.3.0`, each with its
own tests; the [CHANGELOG](../CHANGELOG.md) entry for 0.3.0 lists the
changes and the migration notes.

| Increment | Status |
| --- | --- |
| I1. Self-contained manual examples and rule names | done |
| I2. Event-log budget | done |
| I3. pandas results | done |
| I4. Whole-rule propensity and the `Model` builder | done |
| I5. Generic continuous variables | done |
| I6. Method guidance | done |
| I7. Documentation, examples and release preparation | done |

Known limits:

- `hybrid` is about 11 times slower than `thinning` on the ATP model
  (400 replicates on two threads: about 62 s against 6 s), since every event
  needs a Dormand–Prince step.
- There is no tolerance study yet: the error of `hybrid` has not been
  measured as a function of `rtol` and `atol`; the gates use the defaults.
- A bound crossing that dips below the bound and returns within one
  integration step can be missed.
- Stop conditions and trackers cannot read continuous variables; a stop on
  a variable is written as its `lower_bound`.
- The large ATP equivalence row (4,000 replicates per arm) has p = 0.036
  for ATP exhaustion (416 against 476 stopped) against the predeclared level
  0.0025. It passes, but it is the first place for a bigger independent run
  ([results](../validation/hybrid/RESULTS.md)).
