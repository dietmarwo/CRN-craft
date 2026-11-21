# Two-stage stochastic gene expression

The [guide](../../docs/guide/22-two-stage-expression.md) gives the four
reactions and the limit yielding a negative-binomial protein count.

```sh
python examples/two_stage_expression/study.py --size quick
python examples/two_stage_expression/study.py --size full --threads 4
python examples/two_stage_expression/review.py
python examples/two_stage_expression/study.py --figures docs/guide/img/two_stage_expression
```

`review.py` solves the truncated joint master equation for its
stationary vector, independently of the generating function, and writes
`results/review-stationary-check.csv`. It also recomputes the Figure 1D
divergence under each open convention and writes
`results/review-figure1d-conventions.csv`. The guide page's tables are
rendered from the result files by `tables.py`, and the chapter check
compares them.

The stationary sweeps use a generating function by integrating over
the times of mRNA arrivals and their exponential lifetimes. It gives
the finite-time protein marginal from zero counts. To approximate
stationarity, the run doubles the horizon, doubles inversion support
and refines quadrature. It records residual tail mass and checks the
analytical mean and Fano factor. The joint-state native FSP comparison
for Figure 1B uses an early report at three protein lifetimes and a
bounded lattice. The stationary table uses the marginal generating
function because the larger joint projection exceeds the native work
limit. All four Figure 1B–C histograms have independent native samples.

Figure 1D's method description leaves the direction, binning, tail
handling and unit of the KL divergence open. This example reports the
reaction law relative to the negative-binomial law, grouping counts above
150. `divergence.csv` stores it in nats together with total variation;
the figure shows it in bits. In bits the published values agree within
the figure's reading error; in nats they are about 1.47 times larger. The
figure image was not available to the review, so an axis-calibration
error in `digitize.py` is not ruled out.
[Paper](https://doi.org/10.1073/pnas.0803850105)
