# Rust API

The crate `crn-craft` is the engine behind the Python package and the command
line. Its module `crn_craft::rules` holds the functions below. The minimum
supported Rust version is 1.85. The crate `crn-craft-rules` is a deprecated
facade that re-exports the module, so `crn_craft_rules::X` still compiles
(see [Migrating from 0.6.0](#migrating-from-060)). The documentation of every
public item is in the crate's rustdoc, which CI builds with warnings denied;
this page is the map.

## Entry points

| Task | API |
| --- | --- |
| Check a model | `rules::check_model`, `rules::check_model_file` |
| Check any document, with a request | `check_any`, `check_any_files` |
| Prepare and run | `rules::prepare`, `rules::execute`, `rules::solve`, `rules::solve_files` |
| Observe events | `rules::execute_with_observer`, `rules::observe::EventObserver` |
| Parameter sweeps | `Prepared::bind_parameters` |
| Generated event loops | `rules::generate_rules_loop`, `rules::execute_with_rules_loop` |
| Convert a network model | `rules::convert::from_network` |
| Analysis | `analyze_bytes`, `analyze_files`, `rules::FlatView` |
| Sealed archives | `archive::v1::run_files`, `archive::v1::write`, `archive::v1::verify` |
| Finite state projection | `fsp::v1::solve`, `fsp::v1::solve_files` |
| GillesPy2 export | `export::export_gillespy2` |
| Deterministic kinetics | `kinetics::Kinetics` |
| Discovery | `describe`, `schema`, `rules::schema` |

Report-style functions never fail: they return a report with a `valid` or
`success` flag. In a stochastic simulation, event, sampling, instance and
wall-time limits and cancellation are typed outcomes in each replicate's
record, not errors. Cancel a running call through `rules::Control`.

## Prepare and execute

`prepare` validates exact model and request bytes, and `execute` runs them.
The request binds the model by the SHA-256 of its exact bytes. A prepared run
can be executed again with independent execution state, and
`Prepared::bind_parameters` derives a candidate with other parameter values
for a sweep or an optimizer. A candidate is prepared again from the
re-serialized bound model and request, so its hashes differ from those of the
files.

The program below is complete. It builds the model and the request as bytes,
runs the ensemble and runs two candidates of a parameter sweep.

```rust
use crn_craft::rules::{Control, Termination, execute, prepare};
use std::collections::BTreeMap;

fn main() -> Result<(), String> {
    // Immigration at rate k_in, and death at rate k_out per molecule.
    let model = br#"{
      "schema": "crn-craft.model/1.0", "name": "immigration_death",
      "parameters": {"k_in": 10.0, "k_out": 1.0},
      "types": [{"name": "X"}], "initial": [],
      "observables": {"X": {"type": "X"}},
      "rules": [
        {"name": "immigration", "reactants": [], "rate": "k_in",
         "products": [{"type": "X"}]},
        {"name": "death", "reactants": [{"label": "x", "type": "X"}],
         "rate": "k_out", "products": []}]}"#;
    let request = format!(
        r#"{{
      "schema": "crn-craft.run/1.0",
      "model": {{"path": "model.json", "sha256": "{}"}},
      "method": {{"kind": "static_direct"}},
      "time": {{"end": 5.0, "report_times": [0.0, 5.0]}},
      "ensemble": {{"replicates": 200, "root_seed": 7}},
      "outputs": {{"statistics": ["mean", "variance"], "keep_paths": 0, "keep_instances": 0}},
      "budget": {{"max_events_per_replicate": 100000,
                 "max_sampling_steps_per_replicate": 1000000, "max_instances": 1,
                 "threads": 2, "max_output_bytes": 1048576}}}}"#,
        crn_craft::sha256(model)
    );

    let prepared = prepare(model, request.as_bytes()).map_err(|d| format!("{d:?}"))?;
    let result = execute(&prepared, &Control::default()).map_err(|d| format!("{d:?}"))?;
    assert_eq!(result.termination, Termination::Completed);
    let mean = result.statistics.as_ref().and_then(|s| s.mean.as_ref()).unwrap();
    // The exact law at t = 5 is Poisson with mean 10 (1 - e^-5), about 9.93.
    assert!((mean[1][0] - 9.93).abs() < 1.5);
    println!("{} events, mean X at t = 5: {:.2}", result.events, mean[1][0]);

    // A sweep: the same request and seed with another k_in.
    let values = BTreeMap::from([("k_in".to_string(), 20.0)]);
    let candidate = prepared.bind_parameters(&values).map_err(|d| format!("{d:?}"))?;
    let swept = execute(&candidate, &Control::default()).map_err(|d| format!("{d:?}"))?;
    let swept_mean = swept.statistics.as_ref().and_then(|s| s.mean.as_ref()).unwrap();
    assert!(swept_mean[1][0] > mean[1][0]);
    Ok(())
}
```

## Observing events

`execute_with_observer(&prepared, &control, &observer)` calls an
`EventObserver` for every committed event on the native workers, concurrently
for different replicates, so the observer must be `Sync`. It returns
`ObserverAction::Stop` to end that replicate: the replicate is recorded as
`stopped` at the event's time and holds its state to the end of the report
grid. An observer that sets `ENABLED` to false compiles away, and results stay
byte-identical to those of `execute`. An observer must not panic.

```rust
use crn_craft::rules::observe::{Event, EventObserver, ObserverAction};
use crn_craft::rules::{Control, execute_with_observer, prepare};
use std::sync::atomic::{AtomicU64, Ordering};

/// Counts events and stops a replicate once X reaches 5 copies.
struct StopAtFive {
    events: AtomicU64,
}

impl EventObserver for StopAtFive {
    fn on_event(&self, event: &Event<'_>) -> ObserverAction {
        self.events.fetch_add(1, Ordering::Relaxed);
        if event.observables[0] >= 5.0 {
            ObserverAction::Stop
        } else {
            ObserverAction::Continue
        }
    }
}

fn main() -> Result<(), String> {
    let model = br#"{
      "schema": "crn-craft.model/1.0", "name": "immigration_death",
      "parameters": {"k_in": 10.0, "k_out": 1.0},
      "types": [{"name": "X"}], "initial": [],
      "observables": {"X": {"type": "X"}},
      "rules": [
        {"name": "immigration", "reactants": [], "rate": "k_in",
         "products": [{"type": "X"}]},
        {"name": "death", "reactants": [{"label": "x", "type": "X"}],
         "rate": "k_out", "products": []}]}"#;
    let request = format!(
        r#"{{
      "schema": "crn-craft.run/1.0",
      "model": {{"path": "model.json", "sha256": "{}"}},
      "method": {{"kind": "static_direct"}},
      "time": {{"end": 5.0, "report_times": [0.0, 5.0]}},
      "ensemble": {{"replicates": 16, "root_seed": 1}},
      "outputs": {{"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0}},
      "budget": {{"max_events_per_replicate": 100000,
                 "max_sampling_steps_per_replicate": 1000000, "max_instances": 1,
                 "threads": 2, "max_output_bytes": 1048576}}}}"#,
        crn_craft::sha256(model)
    );
    let prepared = prepare(model, request.as_bytes()).map_err(|d| format!("{d:?}"))?;
    let observer = StopAtFive { events: AtomicU64::new(0) };
    let result = execute_with_observer(&prepared, &Control::default(), &observer)
        .map_err(|d| format!("{d:?}"))?;
    assert_eq!(observer.events.load(Ordering::Relaxed), result.events);
    assert_eq!(result.stopped_replicates, result.requested_replicates);
    Ok(())
}
```

## Generated event loops

*Opt-in.* `generate_rules_loop(&prepared)` emits Rust source for the whole
Static Direct event of one prepared model, for every model that Static Direct
runs (expression and string rates included). It requires a `static_direct`
request and at most 1 MiB of source, and refuses anything else with
`CRN-R020`. Parameters are compiled in, so another parameter set needs another
loop. Compile the source into your program and run it with
`execute_with_rules_loop`, which runs only a loop generated for the same model
bytes. Results are byte-identical to those of `execute`. Observers are not
supported on this path.

On the three-class example it is about 2.1 times as fast as interpreted Static
Direct (timed from Python with `compile_loop`, 12 November 2026). On networks
with many constant-rate channels it can be slower, so measure first
([details](specialized-kernels.md#generated-rules-event-loops)). From Python,
`crn.compile_loop` builds such a loop into a separate program and runs it
([methods reference](../reference/methods.md#compiled-event-loops)).

## Migrating from 0.6.0

The two engines of 0.6.0 became one in 0.7.0. For Rust:

- `crn_craft_rules::X` is now `crn_craft::rules::X`. The old paths still
  compile through the deprecated facade crate.
- Struct literals gain fields: `Method::StaticDirect { .., selector: None }`,
  `Outputs { .., quantiles: None }`, `Statistics { .., quantiles: None }` and
  `RulesCheck { .., run_sha256: None }`.
- `crn_craft::describe()` is `describe/1.0`.
- `Kinetics::compile(&contracts::Model)` is now `Kinetics::from_bytes`, which
  reads every model version. The public fields `Kinetics::model` and
  `Kinetics::network` are gone; use `species()`, `initial()` and
  `parameter_values()`.
- The analysis functions take a `&crn_craft::rules::FlatView`.
- `parameters::ParameterTemplate` is now
  `rules::prepare(model, request)?.bind_parameters(&values)`.
- `crn_craft::tau::poisson_with` is now `crn_craft::sampling::poisson_with`.

The contracts that changed are in [contracts](../reference/contracts.md), and
the Python side of the migration is in the
[migration page](../guide/migrating.md).

## How the examples on this page are checked

The `rust` blocks of this page are complete programs. They are checked as
doctest-style examples: `validation/check_rust_docs.py` writes each block as
the `main.rs` of a scratch crate that depends on the workspace crate
`crn-craft` by path, builds it with `cargo` and runs it, so a block that no
longer compiles or whose assertions fail is reported. The scratch crate is
kept outside the repository. Run it under the CPU lock on the shared host.
