# Parameter sweeps and optimization loops

The old parameter templates (`prepare_parameters`, `ParameterTemplate.bind`,
generated propensity kernels) are replaced by one call on a prepared run:
`bind_parameters`. It is the interface for sweeps and for optimizers that
call the simulator many times with new parameter values.

| 0.6.0 | 0.7.0 |
| --- | --- |
| `crn.prepare_parameters(model, request)` | `crn.ensemble(model, ...)`, or `crn.prepare(model, request)`: the prepared run is the template |
| `template.bind({name: value})`, `template.simulate(values)` | `run.bind_parameters({name: value}).simulate()` |
| `template.parameter_names` | the model's `parameters` |
| `crn.mass_action_kernel_source(model)` | `crn.loop_source(model)` and `crn.compile_loop(model)`: generated event loops, see [the manual](mtdna-user-manual.md#optional-compiled-event-loop) |
| Rust `parameters::ParameterTemplate::{prepare, bind, bind_named}` | `rules::prepare(model, request)?.bind_parameters(&values)` |

## What binding does

`bind_parameters` takes a mapping from parameter names to numbers. It returns
a new prepared run with those values in the model and everything else equal:
the same request, the same seed, the same method. The candidate is validated
and lowered again, which takes about 0.1 to 0.3 ms on a model of a dozen
species. The template and earlier candidates stay intact.

- Only the named parameters change. A name that the model does not declare, or
  a value that is not a finite real number, is refused with `CRN-R020`.
- Initial populations, topology, rules, the request and its method are fixed.
  To change them, prepare a new run.
- The same seed gives every candidate the same random streams (common random
  numbers), so differences between candidates are not noise from different
  seeds alone. Pass another `seed` to `ensemble` for an independent held-out
  evaluation.
- The candidate's `model_sha256` and `run_sha256` hash the re-serialized bound
  model and request. A result names its parameters through them.
- A candidate runs interpreted. A compiled loop serves only the parameter set
  it was compiled for, so `bind_parameters` returns a run without it.
- Results do not depend on the number of threads.

The [user manual](mtdna-user-manual.md#11-parameter-sweeps) runs a sweep over
two parameters. The [optimization tutorials](../tutorials/optimization/README.md)
use the call inside fitting loops, with the comparison in
[the parameter optimization report](evidence/comparison-parameter-optimization.md).

## Limits

Binding prepares the whole model again for each candidate. The 0.6.0
templates reused the compiled topology, so on small, short runs a sweep is
slower than it was (1.5 times on one benchmark; see the known limitations of
the changelog). Reusing the prepared topology is left to a later release.
