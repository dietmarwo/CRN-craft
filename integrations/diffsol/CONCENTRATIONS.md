# Concentration trajectories

`crn-diffsol concentrations MODEL` prepares the model once. It reads one
request per line and writes one result per line. Run `crn-diffsol concentrations
describe` for the three generated schemas.

State values and times are continuous. A model gives ordered parameter and
state names, an initial expression per state, optional named assignments,
and one derivative expression per state. Assignments may read `time`,
parameters, states and earlier assignments. Observation expressions use the
same scope. Supply compartment-volume factors in the derivative formulas.
The solver evaluates continuous Jacobian products through the expression
derivatives. Supported trajectories have nonnegative concentrations and
finite readouts.

Each request supplies up to 32 cases with a parameter vector, 2 to 201
strictly increasing times starting at zero, and relative and absolute
tolerances. A case returns state and observation matrices in time order
or its numerical diagnostic. Results use `crn-craft.concentration-result/0.1`.
The `ode-request/0.1` count kinetics command remains available through
`crn-diffsol MODEL REQUEST`.
