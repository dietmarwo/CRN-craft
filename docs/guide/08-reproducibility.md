# 8. Reproducing runs and published results

Two questions come up in every study. Can I get my own result again, and can
I get the result of a paper? The first has a strict answer: the same model,
request and seed give the same bytes. The second is statistical, and this
chapter gives a recipe for it. The setup block builds the immigration–death
model of [chapter 2](02-first-model.md). The blocks of this chapter continue
from it:

```python
import json

import numpy as np
import crn_craft as crn

model = crn.Model("immigration_death")
model.add_parameter({"birth": 4.0, "death": 0.2})
model.add_species("A", 0)
model.add_reaction("birth", [], ["A"], name="birth")
model.add_reaction("death * A", ["A"], [], name="death")
```

## What makes a run reproducible

Replicate *i* of a run draws from its own ChaCha8 random stream, keyed by
the root seed, with *i* as the stream number. Replicates never share random
numbers, and their results are combined in replicate order. So the numbers
do not depend on the number of threads or on the machine:

```python
digests = set()
for threads in [1, 2, 8]:
    result = crn.ensemble(model, end=50.0, report_times=[0.0, 50.0], replicates=64,
                          threads=threads, keep_paths=4, seed=42,
                          method="static_direct").simulate()
    numbers = dict(result.metadata)
    numbers.pop("run_sha256")  # the run's own hash covers `threads`
    digests.add(json.dumps(numbers, sort_keys=True))
assert len(digests) == 1
print("identical across 1, 2 and 8 threads")
```

```text
identical across 1, 2 and 8 threads
```

Sums of rates are computed in a fixed order, so there is no rounding that
depends on scheduling. A result also names its model and its request by the
SHA-256 of their exact bytes. If either changes by one character, the hash
changes, and a result cannot be attributed to a model that was edited later.

## What to record for a paper

Five items identify a run: the package version, the model hash, the request
hash, the seed and the sampler. The result holds them all:

```python
result = crn.ensemble(model, end=30.0, report_times=[0.0, 30.0], replicates=400,
                      seed=42, method="static_direct").simulate()
for key in ("sampler", "root_seed", "model_sha256", "run_sha256"):
    print(f"{key:13s} {result.metadata[key]}")
saved = result.save("run.json", paths=False)
print(sorted(crn.load_report(saved)))
```

```text
sampler       chacha8-u53-ziggurat/rules-static-direct/1
root_seed     42
model_sha256  9977266823169a5efba1fd1937ade26f112c0f46ae6eb39f751e1bb2c8fa4168
run_sha256    26fa8850983792c48980481e5b20b929918d916e2bb4e7faca005e22e0de2668
['crn_craft_version', 'model', 'request', 'result', 'schema']
```

`result.save("run.json")` writes one file with the model, the request, the
results and the package version, which is what a supplement needs
([chapter 4](04-results.md#save-a-run-in-one-file)). Anyone with that file
and the same release can run the model again and compare. An archive goes
one step further and seals the files with a manifest that a verifier checks
([chapter 10](10-command-line.md#run-archive-and-verify)). In the text of a
paper, name the release, the method and the seed.

## Reproducing a published simulation

A paper reports a model and summary numbers: means, standard deviations,
perhaps a figure. The worked example is model 001-01 of the SBML discrete
stochastic models test suite (Evans, Gillespie and Wilkinson, Bioinformatics
24, 285–286, 2008, [doi:10.1093/bioinformatics/btm566](https://doi.org/10.1093/bioinformatics/btm566)).
It is a birth–death process: `X → 2 X` at rate 0.1 per molecule and `X → 0`
at rate 0.11, starting from 100 molecules. Its mean is `100 exp(-0.01 t)`
and its variance `2100 exp(-0.01 t) (1 - exp(-0.01 t))`.

```python
published = crn.Model("birth_death")
published.add_species("X", 100)
published.add_reaction(0.1, ["X"], ["X", "X"], name="birth")
published.add_reaction(0.11, ["X"], [], name="death")

n, times = 10_000, [10.0, 20.0, 30.0, 40.0, 50.0]
result = crn.ensemble(published, end=50.0, report_times=times, replicates=n,
                      seed=1, threads=2, method="static_direct").simulate()
t = np.array(times)
mu = 100 * np.exp(-0.01 * t)  # the published mean
var = 2100 * np.exp(-0.01 * t) * (1 - np.exp(-0.01 * t))  # the published variance
mean = result.statistics["mean"]["X"].to_numpy()
s2 = result.statistics["sample_variance"]["X"].to_numpy()
z = np.sqrt(n) * (mean - mu) / np.sqrt(var)
y = np.sqrt(n / 2) * (s2 / var - 1)
for row in zip(t, mu, mean, z, var, s2, y):
    print("t {:2.0f}  mean {:5.2f}  simulated {:5.2f}  Z {:5.2f}   "
          "variance {:5.1f}  simulated {:5.1f}  Y {:5.2f}".format(*row))
```

```text
t 10  mean 90.48  simulated 90.48  Z -0.01   variance 180.8  simulated 180.3  Y -0.22
t 20  mean 81.87  simulated 81.64  Z -1.34   variance 311.7  simulated 310.6  Y -0.25
t 30  mean 74.08  simulated 74.00  Z -0.40   variance 403.2  simulated 399.6  Y -0.64
t 40  mean 67.03  simulated 67.09  Z  0.26   variance 464.1  simulated 462.7  Y -0.21
t 50  mean 60.65  simulated 60.82  Z  0.75   variance 501.2  simulated 497.7  Y -0.50
```

The steps, for any paper:

1. **Translate the model.** Check what a rate constant means in the source.
   CRN-craft counts distinct combinations of molecules, so `A + A` with
   constant `k` fires at `k A (A - 1) / 2`. A source whose propensity is
   `k A (A - 1)`, or whose deterministic rate is `k A²`, corresponds to
   `2 k` here. Convert concentrations to molecule counts with the volume.
   If in doubt, write the rate as a string, which is the whole propensity
   ([chapter 2](02-first-model.md#rates-as-numbers-and-as-strings)).
2. **Do not expect the paths.** A seed does not carry over between tools,
   and not between methods or selectors of one tool. Reproduce the
   distribution: the same means and variances within sampling error.
3. **Compare with standard errors.** With the published mean μ and standard
   deviation σ, `Z = sqrt(n) (mean - μ) / σ` is close to a standard normal
   number if the simulation samples the published law. `Y = sqrt(n/2)
   (s² / σ² - 1)` does the same for the variance, less sharply, because
   counts are not normally distributed. The test suite's user guide defines
   these two statistics. It treats `Z` beyond ±3 as evidence of an error,
   and, because `Y` rests on a normal approximation, only `Y` beyond ±5.
   Above, every `Z` and `Y` lies well within ±3. A value beyond the range
   at one of many report times happens by chance now and then; a row of
   them is a discrepancy.
4. **Choose the number of replicates.** The standard error of the mean is
   `σ / sqrt(n)`. To see a difference `d` at three standard errors you need
   `n = (3 σ / d)²`: about 12,000 replicates to resolve 1% of the mean at
   t = 50 here, where σ is 22.4.
5. **Cross-check.** Run a second exact method (`method="next_reaction"`).
   For a small system, compare with the exact distribution from
   `Model.fsp` ([chapter 7](07-methods.md#exact-distributions-with-finite-state-projection)),
   which has no sampling error. Or export the model to another simulator
   ([chapter 11](11-analysis-import-export.md#export-to-gillespy2)).
6. **Plot beside the published figure.** Export the statistics as in
   [chapter 4](04-results.md#export-and-plot).

If the paper used an approximate method with a fixed time step, an exact run
will differ from it by the error of that method. Compare with the exact law
first; [chapter 7](07-methods.md#approximate-leaping) and
[the case study](12d-mtdna-methods.md#matching-an-hourly-scheme) show how to
match such a scheme.

## What does not reproduce bit for bit

- **Another method or selector.** Each has its own sampler name and uses the
  random numbers differently. The law is the same, the paths are not.
- **Another release with a new sampler name.** A release that changes how
  random numbers are used changes the sampler name, so old seeds are never
  silently reused for other paths. Results of the network engine of releases
  before 0.7.0 do not replay bit for bit
  ([Migrating](migrating.md#what-replays-and-what-does-not)).
- **Another request.** Report times and kept paths do not change the
  trajectories, but they change the request hash.

## What a run does not show

A check against an exact law, a published table or a second simulator
confirms that the simulator solves the model. It does not confirm that the
model describes the biology.
[Results and claims](../reference/results.md#claims) says what a result does
and does not claim.

See also: [Results](04-results.md),
[The command line, request files and archives](10-command-line.md),
[Technical report](../evidence/technical-report.md).
