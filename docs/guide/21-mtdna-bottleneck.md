# 21. Inheriting an mtDNA bottleneck

## What this chapter reproduces

A fertilized mouse egg holds about 200,000 copies of mitochondrial DNA
(mtDNA). Over the next week the copy number per cell falls to about
1,000, and it recovers only in the developing egg cells of the next
generation. Cells that start with the same mutant fraction end up with
different mutant fractions. This chapter asks what copy numbers and the
spread of mutant fractions, measured together, reveal about how that
spread arises:

> Johnston I.G., Burgstaller J.P., Havlicek V., Kolbe T., Rülicke T.,
> Brem G., Poulton J., Jones N.S. (2015). Stochastic modelling, Bayesian
> inference, and new in vivo measurements elucidate the debated mtDNA
> bottleneck mechanism. *eLife* 4:e07464.
> [doi:10.7554/eLife.07464](https://doi.org/10.7554/eLife.07464)
> (Table 1, Figure 3A–B, Figure 4 source data)

The paper is open access under the Creative Commons Attribution license.
The chapter takes the measurements of Table 1 and the individual cell
measurements of Figure 4 from the paper. It reads the published Figure 3
curves from the figure by eye. Everything else is simulated or
calculated here. The authors' programs were not used. The paper's authors
did not write or review this chapter.

The result in three lines:

- Approximate Bayesian computation at the paper's distance threshold 40
  reproduces Figure 3 up to 50 days post-conception (dpc). There, every
  reconstructed median lies inside the published outer envelope.
- After 50 dpc the data constrain copy number weakly. The reconstructed
  copy-number band at 100 dpc spans more than three orders of magnitude,
  and its median is 23 times below the published curve.
  [Results](21b-mtdna-bottleneck-results.md)
- The fitted model, which never sees the 15 new mouse groups of
  Figure 4, puts 12 of their 15 measured variances inside the range that
  a sample of that group's size would show.

## Words used in this chapter

- **Heteroplasmy** `h`: the mutant fraction of the mtDNA in one cell.
- **Normalized variance** `V'(h) = Var(h) / (E(h) (1 - E(h)))`: the
  variance of `h` across cells, divided by its largest possible value for
  that mean. It does not depend on the mean to first order.
- **dpc**: days post-conception.
- **Phase**: one of six stretches of development with its own
  replication and degradation rates.

## One division, by hand

At a division each copy enters the tracked daughter with probability one
half, independently of the others. A mother with two mutant copies
therefore passes on zero, one or two of them with these probabilities:

```python
from math import comb

probabilities = [comb(2, copies) * 0.5**2 for copies in range(3)]
print(probabilities)
```

```text
[0.25, 0.5, 0.25]
```

Between divisions each copy is copied at rate `λ` and removed at rate
`ν`. Thirty-six such divisions, with turnover in between, produce the
bottleneck. [The model and the inference](21a-mtdna-bottleneck-model.md)

## Run it

| Size | What it runs | Time, one core with Numba |
|---|---|---|
| `quick` | one ensemble, 16 walkers, 200 sweeps | about 9 s |
| `draft` | four ensembles, 15,000 sweeps each | about 15 s |
| `full` | four ensembles, 360,000 sweeps each (the stored results) | about 1 min |

```sh
python examples/mtdna_bottleneck/study.py --size quick
python examples/mtdna_bottleneck/study.py --size full
python examples/mtdna_bottleneck/study.py --size full --free-final-degradation
python examples/mtdna_bottleneck/review.py
python examples/mtdna_bottleneck/study.py --figures docs/guide/img/mtdna_bottleneck
python validation/check_mtdna_bottleneck.py
```

Numba is optional. Without it the quick size runs in about the same
time, but the draft and full sizes are too slow to be practical. The
check needs neither Numba nor a stored run. It verifies the exact laws,
the reference files, and that every table on the results page equals the
stored result files.

The chapter has three parts:

1. [The model and the inference](21a-mtdna-bottleneck-model.md)
2. [Results](21b-mtdna-bottleneck-results.md)
3. [Choices, checks and limits](21c-mtdna-bottleneck-limits.md)

[Example and result files](../../examples/mtdna_bottleneck/README.md) ·
[Chapter 22](22-two-stage-expression.md) turns to the shape of a different
stochastic count distribution.
