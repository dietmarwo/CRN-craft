# Fixed-topology parameter optimization

`prepare_parameters` validates a model and SSA request once. Each binding
supplies the complete named-parameter vector, while sharing immutable species,
stoichiometry, dependencies and reverse readers. Rate arguments and Hill memos
are resolved afresh; custom-bound groups are rebuilt when parameters split or
join identical laws. All path workspaces are private to an execution.

```python
import crn_craft as crn

model = crn.load_example("immigration-death")
request = crn.ssa_spec(model, fast=True, aggregator="cached_direct",
                       waiting_time="ziggurat", end=10.0, replicates=32,
                       save_times=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0], seed=77201)
template = crn.prepare_parameters(model, request)
print(template.parameter_names)  # ('birth', 'death')
bound = template.bind({"birth": 3.0, "death": 0.8})
result = bound.simulate()
assert result["termination"] == "completed"
# Full consumed snapshots permit ordinary, independently prepared replay:
replay = crn.simulate_ssa(bound.model_bytes, bound.request_bytes)
assert replay == result
```

For an optimizer, call `template.simulate(values)` with a numeric vector in
`parameter_names` order. `bind(values, seed=...)` overrides the root seed for an
independent held-out evaluation. Names/default values/topology identity are
read-only. Mappings must contain exactly the declared parameter names; unused
parameters must still be supplied. Boolean/string/nonfinite/domain-invalid
values fail. Binding errors leave the template and prior bindings intact.
Initial populations, topology, literal coefficients, rate-law kinds, limits and
run settings are fixed. Prepare a new template to change those settings.

The Rust API is `parameters::ParameterTemplate::{prepare, bind, bind_named}`;
`BoundParameters` owns the prepared SSA inputs and exact model/request byte
snapshots. Ordinary `ssa::execute` runs a binding. Template binding still clones
and validates the candidate model and serializes auditable snapshots; it avoids
rebuilding compiled topology/graphs, rather than eliminating all preparation.
Request/result schemas and archive contracts are unchanged. Default checked
semantics stay in force unless the supplied request explicitly selects `fast`.
Completed streams/statistics remain independent of the worker count. Per-path
RNG, selector, range-cache and population state are reset for every execution.

## Optional generated propensities

```sh
cargo run --release --locked --example generate_kernel -- MODEL.json kernel.rs
```

Or obtain the same source with `crn.mass_action_kernel_source(model)`. Compile
it as a Rust module in an application depending on `crn-craft`. The emitted
`GeneratedKernel` implements `ssa::PropensityKernel`. Call
`ssa::execute_with_kernel(&bound.prepared, &control, &kernel)` from Rust.
Named parameter values remain runtime inputs, so the compiled code is reused
across optimizer candidates. A topology fingerprint guards accidental use with
a different model. Trusted handwritten implementations of the trait carry the
same numerical-equivalence obligation as other custom Rust simulator code.

Generation supports mass-action models. Zero/first/dimer/pair evaluation is
specialized; exceptional arithmetic and higher orders use the existing checked
evaluator. Custom Hill/MM laws use the reusable dynamic interface and are
rejected by this generator. Generated kernels retain the core SSA scheduler,
bounds, checked population updates, limits, cancellation and reduction. This is
propensity specialization, not a generated whole solver. The core performs no
runtime compilation or dynamic loading and adds no dependency or unsafe code.

For a complete example and measured compilation/optimization/held-out costs,
see [the reproduction guide](../validation/parameter-optimization/README.md) and
[the full comparison](comparison-parameter-optimization.md). Existing
legacy M6 job applications keep their current execution path; applications can
adopt the new Rust/Python interface without changing the job protocol.
