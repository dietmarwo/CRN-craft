# C5 · Calibrate stochastic RNA counts

The three-rate case study generates RNA distributions under stochastic
promoter switching. Use the nominal rates in Table 1 and reactions in
Supplementary Table S1. Regenerated observations preserve the half-minute
grid and rounded normal error with variance 0.25 in the study preset.
[Poovathingal and Gunawan (2010)](https://doi.org/10.1186/1471-2105-11-414)

```bash
python tutorials/optimization/run.py C5 --preset study --output results/c5
python tutorials/optimization/run.py C5 --preset smoke --source-budget --source-evaluations 90 --output results/c5-ssa
```

The faster arm projects the count law and applies the rounded observation
error analytically. The default sparse-data analysis uses one empirical
midrank per tied integer threshold and evaluates the corresponding model
mid-distribution. `--published-rank` selects the individually ranked
Equation 9 reconstruction. `--source-budget` additionally uses 10,000 SSA
prediction paths per candidate with the original population size and a
120,030-call allowance, including the initial population. The second command
above is a short source-method pilot.

Compare generated-path seeds, rates, objective values and fresh-path scores
in `summary.json`. `cdf-contour.png` uses ten trajectories and fixes the first
rate. The source dataset seeds remain unavailable. Match the nominal vector
and compare recovery across regenerations with Table 1.
