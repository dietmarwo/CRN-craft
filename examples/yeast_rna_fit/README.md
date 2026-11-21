# Fit single-cell RNA count distributions

Follow [chapter 23](../../docs/guide/23-yeast-rna-fit.md) for the data and
fit comparison. The examples use NumPy, SciPy, pandas, fcmaes and a C++17
compiler. The C++ example solves the total-RNA master equation with a
banded BDF2 scheme. The joint nuclear/cytoplasmic solver uses SciPy's sparse
BDF integrator. The native CRN-craft example and validation cover the
constant-signal limit.

```sh
python examples/yeast_rna_fit/study.py --size quick
python examples/yeast_rna_fit/study.py --size draft --workers 4 --cross-validate --budget-check
python examples/yeast_rna_fit/study.py --size full --workers 4 --budget-check
python examples/yeast_rna_fit/study.py --size full --spatial-only --spatial-fit --workers 4
python examples/yeast_rna_fit/study.py --figures docs/guide/img/yeast_rna_fit
```

For the first two fits, all 13 nonspatial coefficients are eligible to
vary within the shared multiplicative box around the authors' count
fit. The same first-moment and likelihood expressions score every
restart. The fit uses fixed empirical sampling covariance and an
integer-count variance floor at low variance. The paper's simulation
also includes intermediate moment-based fits. The spatial exercise
uses 15 coefficients, with initiation in state 3 and cytoplasmic decay
eligible to vary on the checked sparse-solve grid. The spatial result
table gives both the interpolated and directly evaluated objectives.
The fit keeps the best directly evaluated candidate when the interpolant
overstates an improvement.
The spatial run reads the delay from the saved parameter vector and
uses the same forcing in the mature and independent site propagations.

The independent prediction uses the numerical STL1 transcription-site
traces saved in the authors’ figure files and the paper's elongation rate of 63 nt/s. These intensities
are read only after the maturation model parameters are fixed. The
study works with archived single-cell count tuples grouped by dose,
replicate and time. It retains both the reported total and the reported
compartment counts. Where they differ, the respective objectives use
the corresponding source variable. The audit identifies 235 such cells. The published high-dose count
is 39,943, which the transcription-site archive matches. The mature
count file used by the authors’ analysis scripts contains 40,132
high-dose rows. The fit uses the latter file as archived.

Use `fetch.py` to rebuild the source tables from the pinned authors'
archive. It records the SHA-256 of each source MATLAB file. Use the
processed arrays as data and the paper and supplement for equations.
The tests fetch nothing.

[Paper](https://doi.org/10.1073/pnas.1804060115) · [Author archive](https://github.com/MunskyGroup/Munsky_et_al_2018/tree/1aa0b5ce25d21484c7c4e8a894477b1719209ce7)
