# mtDNA mechanism comparison

[Guide chapter 17](../../docs/guide/17-mtdna-mechanisms.md) reproduces the
matching-prevalence experiment of the 2014 papers. Run the commands from
the CRN-craft root with Python 3.10 or newer and a C++17 compiler on Linux
or macOS. Set `CXX` to choose the compiler.

```sh
python -m pip install 'crn-craft==0.8.0' scipy matplotlib
python examples/kk14/study.py single-cell --size quick
python examples/kk14/study.py compare --size quick --ages 3 --threads 2
python examples/kk14/study.py compare --size draft --ages 3 --threads 4 --sensitivities
python examples/kk14/study.py native --threads 2
python examples/kk14/study.py figures --results examples/kk14/results/reconstruction
```

`quick` uses the supplied rates and 128 new cells. It is a smoke example.
`draft` refits the selected conditions and simulates a held-out cohort of
1,000 cells at each fit. `full` uses a narrower calibration grid and an
independent 20,000-cell validation cohort. It includes the ancillary
sensitivities. The paper-size repeat uses 1,000 cells with the rate frozen.
Use distinct `--results` directories for simultaneous invocations. Use
`--resume` after interruption. Completed batches are replaced atomically
and are reused only for the matching sampler code, scientific configuration,
phase and replicate indices. Thread count does not change the simulated
result. Record the compiler and C++ library for a byte-for-byte replay.

```sh
python examples/kk14/study.py compare --size full --ages 3 10 40 80 \
    --threads 8 --seed 20261121 --resume
python validation/check_kk14.py --output kk14-validation.json
```

| Output | What it records |
| --- | --- |
| `fit.csv` | Count data at each calibration probability, including the broad bracket and full-tier local batches |
| `rates.csv` | Frozen probabilities and local logistic intervals for the 12 mechanism–age fits and the sensitivity cases |
| `comparison.csv` | Prevalence, Wilson intervals, both threshold rules, lineage means, bootstrap intervals, largest-mutant fraction and biological or numerical stops |
| `diversity-3y.csv` | One row per observed clone count and its frequency among affected cells |
| `reference-differences.csv` | Numeric differences from the printed and digitized sources, with source reading intervals and both conventions for drift |
| `manifest.json` | Scientific configuration, sampler, seed bases, work sizes and the three-year decision rule |
| `batches/`, `timings.json` | Local per-cell arrays for restart and wall times by cohort, produced by `compare` |
| `single-cell.json`, `native-*.csv`, `native-marks.json` | Deterministic benchmarks and small native solver exercises |

The checked-in [reconstruction results](results/reconstruction/comparison.csv)
use the full tier. Run `figures --results examples/kk14/results/reconstruction
--images docs/guide/img/kk14` to regenerate the guide's tissue figures.
The app reads and writes the example data only. The reference processes
are in `reference.cpp`. The fixed-delay rule makes the hourly scheduler
part of the model. The ATP calculation uses an explicitly declared
reconstruction of the article's algorithm. `native.py` exercises the
CRN-craft pair, static stage and ATP models alongside the study-specific
samplers.

## Reference data

- The PNAS Figure 3 markers were digitized from the 901 × 432 px figure
  embedded in the journal PDF. After reading the tick positions, the
  colored segments were extrapolated to each marker center. Rates are
  read in log space. The source intervals allow ±2 pixels in the ordinate.
  The printed three-year ATP lineage mean of 1.26 takes precedence.
  [Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3939916/)
  · [PDF](https://europepmc.org/articles/PMC3939916?pdf=render)
- The size paper includes the original bars in Figure 6. Figure 7 gives
  the 2-day half-life sensitivity. Bar tops or marker centers use the same
  ±2 px rule. The printed 1.16 is retained at its reported precision.
  The journal published an erratum. Its bibliographic record was available
  during the reconstruction.
  [Article](https://doi.org/10.1016/j.jtbi.2013.09.009)
  · [Erratum](https://doi.org/10.1016/j.jtbi.2014.02.008)
- The 2013 study's 10% column is printed in its Table 1. The matching
  type count uses the convention that includes any surviving wild type.
  The [chapter 16 table](../../docs/guide/16c-fig1cd.md)
  gives the earlier data. The endpoints in chapter 17 also retain the
  mutant-only count. [Paper](https://doi.org/10.1111/acel.12098)

A digitization interval gives the reading precision of the image.
Sampling error in the published mean is unavailable. The bootstrap and
Wilson intervals in the reconstruction describe its own simulated cells.
