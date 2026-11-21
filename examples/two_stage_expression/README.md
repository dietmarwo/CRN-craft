# Two-stage stochastic gene expression

The [guide](../../docs/guide/22-two-stage-expression.md) gives the four
reactions and the limit yielding a negative-binomial protein count.

```sh
python examples/two_stage_expression/study.py --size quick
python examples/two_stage_expression/study.py --size full --threads 4
python examples/two_stage_expression/study.py --figures docs/guide/img/two_stage_expression
```

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

Figure 1D's method description leaves the direction and tail handling
of the KL divergence open. This example reports reaction law relative
to the negative-binomial law, grouping counts above 150. Results include
both this KL divergence and total variation. Open markers in the plot
come from the digitized paper. The decay-rate sweep reproduces the
decrease in approximation error. The ordinates differ most in the
broad-burst case at the smallest decay-rate ratios.
[Paper](https://doi.org/10.1073/pnas.0803850105)
