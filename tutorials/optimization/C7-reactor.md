# C7 · Compare MODE on a virtual reactor

Represent the published simulator as a deterministic concentration model.
The objectives are space–time yield and E-factor. Use residence time,
reagent equivalents, inlet concentration and temperature as the independent
variables of the emulator. The inputs and outputs refer to the numerical
benchmark. [Summit paper, Figure 3](https://doi.org/10.1002/cmtd.202000051) · [Pinned emulator and analysis](https://github.com/sustainable-processes/summit/tree/1de682d05e97adcfb96cd8376e876cef2d6160d3/experiments/snar_benchmark)

```bash
python tutorials/optimization/publications/sources.py C7
python tutorials/optimization/run.py C7 --preset study --output results/c7
```

The archived results record zero measurement noise. Inside fcmaes, the
objectives are negative yield and positive E-factor, both minimized.
Compute hypervolume against the archived reference point `(-2957, 10.7)`.
The study gives MODE and seeded random proposals 50 evaluations in each of
20 starts. All initialization evaluations count toward that total.

`pareto-and-hypervolume.png` plots the two objectives in their original units
and summarises the observed hypervolume trajectories. The report contains
matching archived strategy and transformation summaries for comparison.
