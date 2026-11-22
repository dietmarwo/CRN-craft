# C9 saved results

All files except `profiled-seeds.csv` are the `tables/` folder of one study
run, rerun in review on one thread:

```bash
python tutorials/optimization/run.py C9 --preset study \
  --profile-evaluations 4000 --joint-refine-evaluations 2400 \
  --data-dir DATA --ode-binary BINARY --output OUT
```

Seeds: plain retries 47425 to 47432 (DE) and 48425 to 48432 (CMA-ES), profiled
DE 7458, profiled CMA-ES 7459, joint CMA-ES 7460 (seed offset 0). Each
candidate is scored in the same order whatever the thread count; this rerun
reproduces the trace files and every verified NLL of the implementer's study
record (three worker threads) bit for bit. It took about 100 s.

| File | Content |
|---|---|
| `nll.csv` | NLL per signal and total: archived simulations, nominal vector (native and libSBML routes), best plain retry, profiled fit, joint refinement; `calls` counts objective calls including earlier phases |
| `fits.csv` | Verified NLL and the nine estimated coordinates (log10) of the nominal vector, the eight plain retries, the profiled fit and the joint refinement |
| `measurements.csv` | The 48 measurements with the archived simulation and our values at the same times |
| `curves.csv` | Our signals on 161 equally spaced times from 0 to 240 min |
| `verification.csv` | Native model against the libSBML/SciPy route at the nominal vector and two off-nominal vectors; NLL at the tight and at the search tolerance |
| `profiled-seeds.csv` | The profiled route (with joint refinement) at seed offsets 0 to 4 and once at half budget, rerun in review with `--retries 1 --evaluations 100 --profile-evaluations P --joint-refine-evaluations J --seed-offset S` |

Measurements, archived simulations and the PEtab nominal vector come from the
pinned Benchmark-Models-PEtab problem `Boehm_JProteomeRes2014` (revision
49900951a5594a757ad5e651684c7f0a59272881). The figures `img/C9-*.png` are
drawn from these tables by `publications.stat5.draw`.
