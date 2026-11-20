# mtDNA mechanism comparison

[Guide chapter 17](../../docs/guide/17-mtdna-mechanisms.md) reproduces the
matching-prevalence experiment of the 2014 papers. Its
[last page](../../docs/guide/17c-mtdna-mechanisms-limits.md#reconstruction-choices)
lists every reconstruction choice in one table. Run the commands from
the CRN-craft root with Python 3.10 or newer and a C++17 compiler on Linux
or macOS. Set `CXX` to choose the compiler.

```sh
python -m pip install 'crn-craft==0.8.0' scipy matplotlib
python examples/kk14/study.py single-cell --size quick
python examples/kk14/study.py compare --size quick --ages 3 --threads 2
python examples/kk14/study.py compare --size draft --ages 3 --threads 4 --sensitivities
python examples/kk14/study.py native --threads 2
python examples/kk14/study.py native --threads 2 --cells 4000 --stage-seed 900000 --atp-seed 77000
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
result. A byte-for-byte replay needs the same C++ standard library and C
math library: `std::uniform_int_distribution` and
`std::binomial_distribution` are implementation-defined. The saved full
tier (AMD x86-64, GCC 14.2) was reproduced byte for byte on Intel x86-64
with GCC 13.3. The sampler is compiled with `-ffp-contract=off`, so FMA
targets do not change the rounding.

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
| `single-cell.json`, `native-*.csv`, `native-marks.json` | Deterministic benchmarks and small native solver exercises; `native-size-4000.csv` and `native-atp-4000.csv` come from `native --cells 4000` with the seeds above |
| `draft-*.csv`, `draft-manifest.json` | The earlier draft tier (seed 20261120), kept with its four rows outside the target band |

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
  embedded in the journal PDF. The axes are calibrated on the tick
  positions (panel A: 10⁻² at row 52.0 and 10⁻⁷ at row 370.0; panel B: 40
  at row 32.0 and 0 at row 345.5). A marker center is the center of its
  disk or ring, or the centroid of its triangle; a marker that is partly
  hidden is placed from its visible edge. Rates are read in log space. The
  source intervals allow ±2 pixels in the ordinate. The printed three-year
  ATP lineage mean of 1.26 takes precedence.
  Four panel B rows (drift and size advantage at 3 and 10 years) were
  first read 2 to 3.5 pixels too low, by extrapolating the steep colored
  segments to the marker. They were read again from the marker outlines
  after the full-tier run. The corrected readings agree with the sources
  the PNAS figure cites: 35.7 printed in the 2013 table and 29.8 in the
  size paper's Figure 6B. No simulated value depends on this file;
  `manifest.json` records its hash after a replay from the saved batches,
  and `draft-manifest.json` keeps the hash of the file at the draft run.
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
