# Concentration trajectories

`crn-diffsol concentrations MODEL` reads a `crn-craft.concentration-model/0.1`
document, prepares it once and then answers requests: one
`crn-craft.concentration-request/0.1` document per line on stdin, one
`crn-craft.concentration-result/0.1` document per line on stdout, flushed
after each line, in request order. The process ends with exit code 0 at the
end of its input. `crn-diffsol concentrations describe` prints the three
schemas, generated from the Rust types; the files in `schemas/` are that
output, pretty-printed, and a test compares them. The count-based
`crn-craft.ode-request/0.1` path (`crn-diffsol MODEL REQUEST`) is separate
and unchanged.

## Model

A model names its states and parameters (ASCII letters, digits and `_`, not
starting with a digit, at most 64 bytes, all distinct, `time` reserved) and
gives, in the expression language of the rules engine:

- `initial`: one expression per state, of the parameters only;
- `assignments`: optional named intermediates, in order; each may read
  `time`, the parameters, the states and earlier assignments;
- `derivatives`: one expression per state, the time derivative of the state;
- `observations`: optional readouts, with the scope of the derivatives.

Limits: 1 to 64 states, at most 128 parameters, 128 assignments and 64
readouts, 4,096 bytes per expression, 1 MiB per model file. A refused model
ends the process with a `CRN-E080` diagnostic on stderr.

**Units and meaning.** The solver attaches no units. A state is whatever the
model's derivative is the rate of; in the tutorials it is a concentration
(mol/L, or the paper's concentration unit), and time is in the unit of the
rate constants. A flux is an assignment: amount per volume per time when it
is added to a concentration derivative directly. Compartments are explicit:
for an amount flux J between compartments of volumes V₁ and V₂, write
`dC1/dt = −J/V1` and `dC2/dt = J/V2` with the volumes as parameters; then
`V1·C1 + V2·C2` is conserved. The solver does not check that a volume is
positive. A readout is the value of its expression at a requested time.

**Negative states.** The expressions see `max(x, 0)` for every state `x`,
and the Jacobian treats a negative state as constant, so a trial step below
zero cannot drive a rate expression outside its domain. A reported state
below `−(atol + rtol · max_t x)` fails the case.

## Request and result

A request holds 1 to 32 cases. A case gives the parameter vector (finite,
the model's length), 2 to 201 strictly increasing finite times starting at
0, and the tolerances `rtol` in [1e-13, 1e-2] and `atol` in [1e-15, 1e-2].
The integrator is diffsol's variable-order BDF with these per-step local
error tolerances (`atol` for every state); the global error is larger,
in the checks about five times `rtol` relative. States are interpolated at
the requested times; at most 20,000 steps are taken per case.

The result has either `cases`, one entry per case in order, or `error` for
a line refused as a whole. A case entry is `{"values": {"states": …,
"observations": …}}`, matrices with one row per time, or
`{"error": {"code": …, "message": …}}`. Codes:

| Code | Meaning |
| --- | --- |
| `CRN-E001` to `CRN-E003` | the line is not valid JSON, has a duplicate key, or has a wrong type, an unknown field or another schema |
| `CRN-E004` | a line over 1 MiB, or not 1 to 32 cases |
| `CRN-E015` | times, tolerances or parameters out of their domain, or a negative or non-finite initial state |
| `CRN-E080` | the integration failed (step size, nonlinear solver, step limit), a reported state is negative beyond the tolerance, or a readout is not finite |

A failed case does not stop the batch or the process, so a fitting client
can record the candidate as failed. The same request gives the same bytes.
There is no time limit inside the process; a client that needs one stops
the process (the Python client in `tutorials/optimization/publications/`
does).
