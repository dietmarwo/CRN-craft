# Optional deterministic path: numerical gates (declared before the campaign)

Pinned diffsol 0.18.0 BDF. Default Rust core remains MSRV 1.85; this optional
workspace declares 1.89, to be checked against resolved dependency metadata.
Use rtol 1e-10 / atol 1e-12 on linear and small nonlinear fixtures.

- Analytic immigration/death states and raw-rate sensitivities: absolute error
  <= 1e-7 at t = 0, 0.5, 1, 2.
- Nonlinear enzyme and three-gene deterministic solutions against SciPy DOP853:
  max scaled error <= 1e-5 at t = 0, 0.5, 1, 2, 5.
- Sensitivities against independent central differences (h = 1e-4 max(1,p)):
  max scaled error <= 1e-4. Also compare gradients in the positive interior.
- Explicitly reject unsupported gates, singular Hill boundaries and negative
  states. No deterministic trajectory is labelled an exact stochastic mean.
- Fit synthetic linear time-course data with fcmaes-rust followed by SQP-craft;
  re-evaluate the final loss and forward-sensitivity gradient independently.
- A deterministic screen may only prioritize topology candidates. All candidates
  retain a nonzero stochastic scoring allocation. Show this on the noise-driven
  control and ring/broken-ring fixtures, with withheld seeds accounted separately.
