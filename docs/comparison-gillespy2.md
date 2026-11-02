# CRN-craft compared with GillesPy2

Measured on 1 November 2026 with CRN-craft 0.2.0-dev.1 (main 88677c0) and
GillesPy2 1.8.3, Python 3.13 and g++ 13.3, on one shared 2-vCPU cloud host.
Script: [`validation/compare_gillespy2.py`](../validation/compare_gillespy2.py);
raw report: [`validation/campaign/gillespy2-comparison.json`](../validation/campaign/gillespy2-comparison.json).

## Summary

On the four frozen benchmark networks, CRN-craft's `simulate-cached` beats
GillesPy2's compiled C++ solver (`SSACSolver`) on every network:

- 1.8 times faster on a two-channel network;
- 3.3 times faster on a sparse 66-reaction network;
- 1.1 times faster on the 12-reaction repressilator;
- 1.1 times faster on a dense 66-reaction Michaelis–Menten network.

Against GillesPy2's pure Python solver (`NumPySSASolver`), CRN-craft is 43 to
126 times faster. GillesPy2's C++ solver also compiles each model before its
first run, which takes 1.5 to 2.1 seconds.

![Time per SSA event for CRN-craft simulate-cached, CRN-craft simulate and GillesPy2's C++ solver on W1–W4](images/gillespy2-time-per-event.png)

## What was measured

| ID | Network | Reactions | Rate laws |
|---|---|---|---|
| W1 | immigration–death | 2 | mass action |
| W2 | 33 independent immigration–death species | 66 | mass action |
| W3 | three-gene repressilator | 12 | Hill repression and mass action |
| W4 | eight species born and removed together | 66 | mass action and gated Michaelis–Menten |

Each run simulates 64 paths over t from 0 to 100, with output at 0 and 100, on
one thread. These are the frozen workloads of the performance campaign. Each
row below is the median of five interleaved runs (three for the NumPy solver).

How each column was timed:

- **CRN-craft:** the wall time of the Python call `crn.simulate(model, run)`,
  which includes JSON encoding and decoding. The command line gives the same
  times within about 5%.
- **GillesPy2:** the wall time of `model.run(solver=..., number_of_trajectories=64)`
  with a solver object built beforehand. For the C++ solver this excludes the
  one-time compilation, which is reported separately.

Time per event uses CRN-craft's event count for both tools. Both simulate the
same stochastic process, so the expected number of events is the same.

## Results

Median seconds for 64 paths:

| | W1 | W2 | W3 | W4 |
|---|---|---|---|---|
| CRN-craft `simulate-cached` | **0.061** | **0.190** | **0.553** | **0.696** |
| CRN-craft `simulate` | 0.063 | 0.613 | 1.223 | 0.800 |
| GillesPy2 `SSACSolver` (C++) | 0.110 | 0.634 | 0.621 | 0.756 |
| GillesPy2 `NumPySSASolver` | 4.16 | 23.8 | 26.6 | 30.0 |
| GillesPy2 C++ compile (once per model) | 2.1 | 1.7 | 1.5 | 1.6 |
| Events (CRN-craft) | 1,268,165 | 1,234,576 | 3,663,350 | 1,225,923 |

Time per event, in nanoseconds:

| | W1 | W2 | W3 | W4 |
|---|---|---|---|---|
| CRN-craft `simulate-cached` | 48 | 154 | 151 | 568 |
| CRN-craft `simulate` | 49 | 496 | 334 | 652 |
| GillesPy2 C++ | 87 | 513 | 170 | 617 |
| GillesPy2 NumPy | 3,278 | 19,292 | 7,259 | 24,449 |

How many times longer GillesPy2 takes than CRN-craft:

| | W1 | W2 | W3 | W4 |
|---|---|---|---|---|
| GillesPy2 C++ / CRN-craft `simulate-cached` | 1.81 | 3.34 | 1.12 | 1.09 |
| GillesPy2 C++ / CRN-craft `simulate` | 1.76 | 1.03 | 0.51 | 0.95 |
| GillesPy2 NumPy / CRN-craft `simulate-cached` | 68 | 126 | 48 | 43 |

Run-to-run spread was small. Every timed run of every column lies within 8%
of its median, and most within 4%. The largest deviations were 7.6% for
CRN-craft `simulate` on W1, 8.0% for one CLI run on W4 and 6.4% for GillesPy2's
C++ solver on W3.

## Reading the results

- **The two compiled solvers use different per-event strategies.** In
  GillesPy2 1.8.3 (`SSASolver.cpp`), after each event the C++ solver:
  - recomputes only the propensities the fired reaction affects;
  - sums all propensities again;
  - picks the reaction by a linear scan;
  - loops over all species to apply the update.

  CRN-craft's `simulate-cached` also recomputes only the affected
  propensities. In addition, it resumes its running sums from the first
  changed reaction and finds the reaction by binary search over those sums.
  CRN-craft's plain `simulate` recomputes every propensity after every event.
  That is why plain `simulate` is about twice as slow as GillesPy2 on the
  repressilator, whose Hill laws are expensive, while `simulate-cached` is
  faster on all four networks.
- **The largest gap is on the sparse network (W2).** Each reaction there
  changes one species and therefore one propensity. GillesPy2 then spends most
  of an event summing and scanning 66 propensities and 33 species. CRN-craft
  re-sums only from the changed reaction onward (half the channels on average)
  and finds the reaction by binary search in about six steps.
- **CRN-craft does more work per event and is still faster.** On every event it
  checks that no positive propensity is lost to rounding or too small to be
  selected. It also guarantees byte-identical results across threads and
  between its two implementations, and seals results in checksummed archives.
  GillesPy2 makes none of these checks.
- **Compilation matters for short jobs.** GillesPy2's C++ solver compiles each
  model before its first run, in 1.5 to 2.1 seconds. That is about 20 times
  the W1 run itself. CRN-craft needs no per-model compilation.
- **The NumPy solver** is 43 to 126 times slower than CRN-craft. That is as
  expected for a solver that evaluates propensities in Python.

## The two tools simulate the same process

Before timing, both simulators ran 512 paths of each network with fixed,
different seeds:

- **Final counts agree.** Two-sample Kolmogorov–Smirnov tests and mean z-scores
  compared the final counts of up to three species per network.
  - W1–W3: five of the seven KS p-values are above 0.05.
  - The two exceptions: W2 species `X2` (mean z 3.14, KS p 0.018) and W3
    species `m2` (KS p 0.038).
  - W4's three species are identical within each path, so all three give
    KS p 0.34.
- **The exact law settles the borderline cases.** W1 and W2 have Poisson laws
  at t = 100, with means 100 and 15. Against those exact laws, both simulators
  pass:
  - CRN-craft: pooled chi-square p 0.17 on W1 and 0.71 on W2, and the root mean
    square of the 33 W2 mean z-scores is 1.03.
  - GillesPy2: pooled chi-square p 0.20 on W1 and 0.81 on W2, with a root mean
    square z of 1.09.
  - The W2 `X2` difference is therefore chance. A separate run of CRN-craft
    with six more seeds and 2,000 paths each gave a root mean square z of 0.99
    over 198 species means.
- **W4:** all eight species stay equal in every path of both simulators, as the
  network's structure requires.

## How the models were translated

- **Mass action:** as GillesPy2 mass-action reactions. W1–W4 have only zeroth-
  and first-order reactants, so no combinatorial convention differs.
- **Hill repression:** as the propensity expression
  `basal + vmax / (1 + pow(x/K, n))`.
- **Michaelis–Menten:** as `vmax * S / (Km + S)`.
- **Reactant gates.** CRN-craft fires a custom-rate reaction only when its
  reactants are present. GillesPy2 expressions have no comparisons, so a gate
  is written as `x / (x + 1e-300)`, which is exactly 0 for x = 0 and 1
  otherwise.
- **W4 is timed without gates.** Its eight species always stay equal, so the
  gates on its Michaelis–Menten deaths never change a propensity. The
  validation confirms this. With the gates, GillesPy2's C++ solver takes 7.08
  seconds instead of 0.76, about 10 times longer. We consider that cost an
  artifact of the translation rather than of GillesPy2, so the comparison
  uses the faster form. The gated time is kept in the raw report.

## Limitations

- **Host:** one shared 2-vCPU cloud host. Absolute times differ on other
  machines, and the ratios may too.
- **Random streams** are independent, so this compares cost for the same
  process, not identical paths.
- **Event counts:** GillesPy2's are not instrumented. Its time per event uses
  CRN-craft's count, and the 64-path counts of the two tools differ by sampling
  noise.
- **What the timings include:** GillesPy2's include building its result
  objects and, for the C++ solver, launching its compiled program and parsing
  its output. CRN-craft's include JSON encoding and decoding in Python.
- **Scope:** only exact SSA on one thread. CRN-craft can also divide
  replicates among threads without changing results. GillesPy2 also offers
  tau-leaping, hybrid and ODE solvers, which are not compared here.
- **Workload size:** the networks are small to medium (up to 66 reactions).
  Larger networks favor dependency-based updates, as W2 shows.

## Reproducing

```sh
python -m pip install -r validation/reference-requirements.txt
maturin build --release --locked --out dist && python -m pip install dist/*.whl
cargo build --release --workspace --locked
cd validation
python compare_gillespy2.py --binary ../target/release/crn-craft \
    --output gillespy2.json --figure gillespy2.png
```

The script needs a C++ compiler for GillesPy2's `SSACSolver` and takes about
ten minutes, mostly in the NumPy solver.
