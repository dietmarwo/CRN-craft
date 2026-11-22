# C10 · Fit 36 kinetic parameters across 16 conditions

Generate eight time series at each published substrate/product combination.
The reconstruction uses the starting concentrations and six-second sampling
from the later Table 1 case study. It uses the pathway's third transcript
in the last translation equation. `--printed-equation-12` selects the
second transcript in the 2003 printed numerator. Both variants preserve
the explicit source choices in their reconstruction files.
[2003 paper and Equation 15](https://doi.org/10.1101/gr.1262503) · [Later reconstruction, Table 1](https://citeseerx.ist.psu.edu/document?doi=f05547c0006562be54acfe8d0ceaa917d693a5ad&repid=rep1&type=pdf)

```bash
python tutorials/optimization/run.py C10 --preset study --evaluations 2000 --output results/c10-2k
python tutorials/optimization/run.py C10 --preset study --evaluations 20000 --resume-from results/c10-2k --output results/c10-20k
```

Each candidate scores all 16 conditions. Normalize each state by its largest
measured value in that condition. Compare actual evaluations and state fit
quality between DE and CMA-ES. Both the publication-scale bounds and a
narrower, declared decade box around the published parameter scales are
reported. The study has checkpoints at 2,000, 20,000 and 200,000 calls.
Replayed proposals consume their original portion of the budget without
requiring the simulations to be run again. Unsuccessful numerical evaluations
and independently rejected finalists remain visible in the trace and summary.

The printed best score of 0.0013 is a historical reference. The reconstructed
time grid and initial conditions come from the later study, so the report
compares its own objectives on common inputs and also checks a denser grid.
