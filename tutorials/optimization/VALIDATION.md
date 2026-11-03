# Validation record

C1–C4 are simulator-owned Python and CLI lessons; they change no core SSA, rule,
ODE or optimizer algorithm. C3 uses a bounded reference example inside
fcmaes-rust's oscillator workspace whose backend is the local ReBop
compatibility implementation.

[results/verified.json](results/verified.json) holds the numbers below for one
complete smoke and one complete study run of the committed lessons (C1 with
`--include-checked`, `--threads 2`). [record.py](record.py) writes it from the
runs' summaries and refuses incomplete or failed runs:

```bash
python tutorials/optimization/record.py --smoke SMOKE_DIR --study STUDY_DIR \
  --output tutorials/optimization/results/verified.json
```

SSA runs use the default `fast` policy (`crn-craft.ssa-run/0.3`, sampler
`xoshiro256pp-u53/...`); C1's checked arm uses `fast=False` (`chacha8-u53/...`).
The two policies draw different paths, so their rows differ and are never
pooled.

| Check | Result (smoke; study) |
|---|---|
| C1 fresh and reusable candidate streams, ordinary replay, held-out moments | Equal streams and replay for every arm and both policies; largest held-out mean error 2.0 standard errors (smoke), 2.2 (study), against a gate of 6 |
| C1 training loss of the fitted versus the true rates | Fitted below truth in most arms, e.g. study binding, checked: 0.012–0.014 against 0.108; the fits follow the shared training sample |
| C2 integration and log-gradient audits | Gradient error 9.1e-10, state error 6.5e-11, tightened 3.2e-12 |
| C2 local identifiability | Full-observation condition number 2.6 at the truth; a terminal-product twin (0.74, 0.84) matches P(T) to 2e-13 with full-trajectory RMSE 0.029 |
| C2 fits | Multistart and global-local SQP recover (1.0, 0.8); BiteOpt alone 12% rate error (smoke), 0.4% (study) |
| C3 generator, cross-backend control, MODE/archive and held-out paths | Generator error below 1e-14; control within 1.4 (smoke) and 1.7 (study) standard errors; archives occupy 3 of 16 cells in both presets |
| C3 oscillation behind the selections | Quality-optimal Pareto design: 75% of training and 31% of held-out paths oscillatory (smoke), 94% and 78% (study); archive elites at most 6% |
| C4 mass and four-marginal law gate per kernel | Largest held-out CDF distance 0.030 (smoke) and 0.011 (study) against bounds 0.066 and 0.023 |
| C4 fitted rate versus exact-law reference | Constant kernel 0.87 (smoke) and 1.17 (study) against 1.2; multiplicative 0.60, the start (smoke), and 0.64 (study) against 0.67 |
| Fifteen Python unit and negative checks | Pass normally and under `python -O` |

C1's start and box are multiples of the true rates, and one training sample
serves all optimizer seeds; `truth_training_loss` and `best_is_start` make the
consequences visible. C2's terminal-product rank of one holds by dimension; the
twin pair shows what it means. C3's archive is sparse and its period cells
beyond about 33.5 time units cannot hold oscillatory paths (see the resolution
note in C3); completion does not establish that a niche contains a robust
oscillator. The Vilar control is only a short conservation test of a separate
model.

Timing is descriptive and depends on host load. Independent training and
validation seeds and native work counts establish the numerical protocol, not
statistical power for a scientific discovery. Clone inference, learned
topology proposals, stock MO/QD workers and generated-kernel performance
studies are follow-up work.
