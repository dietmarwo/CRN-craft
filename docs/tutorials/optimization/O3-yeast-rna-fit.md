# O3. Case study: fitting yeast single-cell RNA distributions with fcmaes

## What this chapter reproduces

Cells sampled at the same time after a stress do not all hold the same
number of RNA molecules. Some hold none, others hold a hundred. A 2018
paper asks what a model learns when it is fitted to these whole count
distributions instead of only to their averages:

> Munsky B., Li G., Fox Z.R., Shepherd D.P., Neuert G. (2018). Distribution
> shapes govern the discovery of predictive models for gene regulation.
> *PNAS* 115(29):7533–7538.
> [doi:10.1073/pnas.1804060115](https://doi.org/10.1073/pnas.1804060115)
> (Figure 2A–D and Figure 4D–E)

The measurements are described in a data paper:

> Li G., Neuert G. (2019). Multiplex RNA single molecule FISH of inducible
> mRNAs in single yeast cells. *Scientific Data* 6:94.
> [doi:10.1038/s41597-019-0106-6](https://doi.org/10.1038/s41597-019-0106-6)

The study uses the processed count arrays and the fitted parameter vectors
that the authors published in their repository
([MunskyGroup/Munsky_et_al_2018 at revision 1aa0b5c](https://github.com/MunskyGroup/Munsky_et_al_2018/tree/1aa0b5ce25d21484c7c4e8a894477b1719209ce7)).
Their MATLAB code was read for facts only: which replicates and time
points the fit used, how the parameter vector is laid out. Nothing was
copied or translated from it. The equations were checked against the
bioRxiv version of the paper and the first part of its supplement
([10.1101/154401](https://doi.org/10.1101/154401)); the journal
supplement was not reachable. The paper's authors did not write or review
this chapter.

The result in four lines:

- The same four-state model was fitted to the measured means and to the
  measured count histograms. Scored on whole histograms, the fit to means
  costs 2.39690 nats per cell and the fit to histograms 2.03845 (59,557
  cells). Both scores reproduce from the saved parameters.
- Both fits searched a box centered on the authors' published histogram
  fit. The histogram fit stays close to that vector; the mean fit ends at
  the box edge. So the comparison shows how the two objectives rank models
  near the published one. It does not show what each objective finds from
  no prior knowledge, and the held-out scores are not clean held-out tests
  ([why](O3b-yeast-rna-fit-results.md#what-the-comparison-shows)).
- The spatial model was not refitted. The best point of its 3 × 3 check
  grid is the authors' published vector.
- The transcription-site prediction recomputes the authors' prediction from
  their published vector. It agrees with their model curve to 0.0024 in
  active fraction. Against the measured traces it misses by up to 0.42 at
  0.4 M.

## Words used in this chapter

- **Count distribution, histogram**: for one dose, replicate and sampling
  time, the number of cells with 0, 1, 2, … RNA molecules.
- **Nats per cell**: the negative log likelihood of a histogram divided by
  its number of cells. Smaller is better. A difference of 0.01 nats per
  cell over 59,557 cells is 596 nats.
- **Archived vector**: a parameter vector saved in the authors' repository
  (`Fit_Results/.../Best_Pars_AC.mat`). There are four: means and
  histogram fits of the nonspatial and of the spatial model.
- **Replicate**: one biological time course, named by its date. 0.2 M has
  two (20140425, 20140429), 0.4 M has three (201301004, 20131022,
  20160509).

## Score a distribution

A count histogram keeps the cells without RNA as well as the cells with
many. If a bin holds `n_k` cells and the model gives it probability `p_k`,
the bin adds `-n_k log p_k` to the negative log likelihood.

```python
from math import log

counts = [7, 2, 1]
probabilities = [0.6, 0.3, 0.1]
score = -sum(n * log(p) for n, p in zip(counts, probabilities, strict=True))
print(f"{score / sum(counts):.4f} nats per cell")
```

```text
0.8286 nats per cell
```

A fit to means sees only the average of each histogram. Two models with
the same means can give very different histograms, and the histogram score
tells them apart.

The pages of the chapter:

| Page | Question |
| --- | --- |
| [Data and model](O3a-yeast-rna-fit-model.md) | Which cells are used, what is the model, and how is it solved? |
| [Fits and predictions](O3b-yeast-rna-fit-results.md) | What do the two fits give, what does the optimizer comparison show, and how good is the transcription-site prediction? |
| [Choices, checks and limits](O3c-yeast-rna-fit-limits.md) | What did the reconstruction have to decide, how was it checked, and what is not verified? |

## Running the example

You need Python 3.11 or newer, CRN-craft 0.8.0, NumPy, SciPy, pandas,
fcmaes 2.2.0, Matplotlib and a C++17 compiler for the fits (set `CXX` to
choose one). The compiled solver is cached in the system's temporary
directory. Run the commands from the repository root:

```sh
python examples/yeast_rna_fit/study.py --size quick
python examples/yeast_rna_fit/study.py --figures docs/guide/img/yeast_rna_fit
python validation/check_yeast_rna_fit.py
```

`--size quick` fits 0.2 M at five sampling times with 80 evaluations per objective in
about three seconds; it shows the workflow, not a result. `--figures`
redraws the three figures from the saved result files. The check needs
about a second and no network. The larger runs and their cost are listed
on [the results page](O3b-yeast-rna-fit-results.md#rerunning-the-study);
the [example README](../../../examples/yeast_rna_fit/README.md) lists the
files.

See also: [Protein distributions from four reactions](../simulation/S7-two-stage-expression.md)
(the count distribution of a simpler model),
[Fitting](O1-fitting.md) (fcmaes with CRN-craft).
