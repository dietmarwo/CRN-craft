# C8 · Optimize a sampling schedule with Fisher information

Use the observed probability law to measure information about fitted rate
parameters. Sum state probabilities and derivatives into measured bins
before constructing Fisher information. The burst experiment uses an
RNA-only readout. The toggle design uses both counts.
[Paper, Figures 5 and 7](https://doi.org/10.1371/journal.pcbi.1006365) · [Author archive](https://github.com/Munsky/FSPFIM_2018/tree/2e2f76fdeeb0d9290b48aa981210422d254b07c2)

```bash
python tutorials/optimization/publications/sources.py C8
python tutorials/optimization/run.py C8 --preset smoke --output results/c8-smoke
python tutorials/optimization/run.py C8 --preset study --native-samples 3 --native-evaluations 10 --native-retries 2 --native-pair-samples 100 --threads 7 --cache results/fim-cache --output results/c8-paper
```

Check the burst maximum near 6.1 min. The 50-point grid in the author script
and the approximate 60-point grid described in the article are both tested.
For the toggle, reconstruct the 672 archive scores and compare the two
published designs across 100 parameter samples. The robust objective
averages the smallest eigenvalues of individual Fisher matrices.

Native runs carry probability and tangent vectors through fixed-parameter
segments. The published pulse law and the archived code have different
constants. `--pulse-law archive-code` selects the constants in the code,
while the default uses Equation 28. In the archive coordinate order the
stimulus acts on its first species. The native timing uses four evenly
spaced observations, following the experiment schematic. Check the sink,
marginal cutoff and source convention before interpreting comparisons.
`--native-grid --cache PATH` recomputes all designs and reuses completed
sample evaluations on subsequent runs. A complete native search with all 100
samples can be expensive. The report records exactly which samples and
search budgets were used.
