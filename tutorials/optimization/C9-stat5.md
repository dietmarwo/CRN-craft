# C9 · Fit a deterministic STAT5 model from PEtab tables

The pinned PEtab v1 problem supplies 48 measurements, three observations,
six estimated rates and three estimated Gaussian noise scales. Its SBML
model has fractional initial concentrations, two compartments, and a
time-varying input. The concentration model translates rates into each
compartment's concentration units. The reference calculation solves the
original reaction amounts via libSBML and SciPy.
[Paper](https://doi.org/10.1021/pr5006923) · [Archived problem](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/49900951a5594a757ad5e651684c7f0a59272881/problems/Boehm_JProteomeRes2014/v1)

```bash
python tutorials/optimization/publications/sources.py C9
python tutorials/optimization/run.py C9 --preset study --output results/c9
```

Verify the archived negative log-likelihood 138.2219997062 before running
DE and CMA-ES. Retain all Gaussian normalization terms because the noise
standard deviations vary. The score reader supports one condition, one
model, fixed and estimated parameters, linear-scale measurements and
Gaussian noise from one named parameter per row. The SBML-to-CRN mapping
is explicit for this case. `--profile-evaluations` adds a comparison that
solves each noise scale analytically for a given kinetic vector.
`--joint-refine-evaluations` then reopens all nine variables in CMA-ES.
