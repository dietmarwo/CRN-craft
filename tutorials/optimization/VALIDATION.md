# Validation record

C1–C4 are simulator-owned Python and CLI lessons; they change no core SSA, rule,
ODE or optimizer algorithm. C3 uses a bounded reference example inside
fcmaes-rust's oscillator workspace whose backend is the local ReBop
compatibility implementation.

[results/verified.json](results/verified.json) holds the numbers below for one
complete smoke and one complete study run of the committed lessons (C1 with
`--include-next-reaction`, `--threads 2`), made on the unified 0.7.0 API. [record.py](record.py) writes it from the
runs' summaries and refuses incomplete or failed runs:

```bash
python tutorials/optimization/record.py --smoke SMOKE_DIR --study STUDY_DIR \
  --output tutorials/optimization/results/verified.json
```

SSA runs are Static Direct ensembles of `crn-craft.run/1.0` (sampler
`chacha8-u53-ziggurat/rules-static-direct/1`); C1's second arm uses the
`next_reaction` method (`chacha8-u53/rules-next-reaction/1`). The two methods
draw different paths, so their rows differ and are never pooled. The
0.6.0 record (xoshiro256++ `fast` and ChaCha8 `checked` policies) is replaced:
the random stream changed with the merge, so the sampled numbers below differ
from it; the deterministic ones (C2 audits and ODE fits, C4, which ran on the
rules engine before) are unchanged, and every threshold is unchanged.

| Check | Result (smoke; study) |
|---|---|
| C1 fresh and reusable candidate streams, ordinary replay, held-out moments | Equal streams and replay for every arm and both methods; largest held-out mean error 2.3 standard errors (smoke), 2.8 (study), against a gate of 6 (0.6.0: 2.0 and 2.2) |
| C1 training loss of the fitted versus the true rates | Fitted below truth in most arms, e.g. study binding, next-reaction: 0.005–0.013 against 0.025; the fits follow the shared training sample |
| C2 integration and log-gradient audits | Gradient error 9.1e-10, state error 6.5e-11, tightened 3.2e-12 |
| C2 local identifiability | Full-observation condition number 2.6 at the truth; a terminal-product twin (0.74, 0.84) matches P(T) to 2e-13 with full-trajectory RMSE 0.029 |
| C2 fits | Multistart and global-local SQP recover (1.0, 0.8); BiteOpt alone 12% rate error (smoke), 0.4% (study) |
| C3 generator, cross-backend control, MODE/archive and held-out paths | Generator error below 1e-14; control within 1.1 (smoke) and 1.8 (study) standard errors; archives occupy 2 of 16 cells in both presets (0.6.0: 3) |
| C3 oscillation behind the selections | Quality-optimal Pareto design: 75% of training and 62% of held-out paths oscillatory (smoke), 81% and 66% (study); archive elites at most 15% |
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

## C5–C10 publication studies

C5–C10 reconstruct published optimization studies on the same engine and the
same `run.py` driver; they add the finite state projection with retained sets,
forward sensitivities on a prepared projection and the persistent
concentration solver of the diffsol workspace, and change no SSA, rule or
optimizer algorithm. [verified.json](results/verified.json) records C1–C4 only: the
studies of C5–C10 take hours on eight threads, and the record keeps one smoke
and one study run of the chapters that finish in minutes.

Each chapter keeps its compact tables in `results/C#/`, and `results/C#/README.md`
names for every file where it came from: the study run of the implementer on
eight threads (kept as delivered, not rerun in review), or a run repeated in
review with its command, seed and thread count. The numbers of a chapter are
on its own page, next to the published value; this page repeats none of them.
Figures are drawn only from these tables.

The smoke preset of each chapter checks the workflow with a small budget; it
does not reproduce a published number. It runs the gates that are fixed in the
chapter's module before any search:

| Chapter | Smoke gates (check names in `summary.json`) |
|---|---|
| C5 | FSP against SSA with the count noise of the data, the two-rank score, the search ran |
| C6 | the 376 archived cells (pinned files, SHA-256), the joint likelihood at a fixed vector |
| C7 | the archived reference noise, outlet parity of the virtual reactor, a hand-computed hypervolume, exactly 50 evaluations |
| C8 | probability sensitivities against finite differences and a Fisher fixture, the bursting-interval peak near 6.1 minutes, the archived grid and its 97-of-100 ranking |
| C9 | the archived Gaussian likelihood and the nominal signal at the published parameters |
| C10 | the sixteen conditions against an independent trajectory, a hand-computed weighted loss |

`record.py --publications DIR --output results/publications-smoke.json` writes
the check names and verdicts of complete smoke runs of C5–C10 (one folder per
chapter, `C5/` to `C10/`) and the package versions; it refuses a missing,
failed, non-smoke or mixed-version run, and it records no number, path or
timing.

Every search arm writes `NAME.csv` (the proposals in evaluation order with
their scores and error text) and `NAME.run.json` (seed, optimizer, population,
bounds, context, backend and its version, failure count, replayed rows and
the elapsed time). The CSV files and the run.json fields other than
`elapsed_seconds`, `workers`, `replayed` and the requested budget are
deterministic for a fixed seed; elapsed times depend on the host and are kept
out of the traces. `--resume-from` replays the completed proposals of an
earlier, completed output and continues the same optimizer sequence; it
refuses a trace that is cut, edited, incomplete or from another seed,
optimizer, population, bounds or context. A run that stopped before it wrote
its `run.json` cannot be resumed. Whether a context string names the objective
completely is the chapter's duty: the replay trusts the saved scores.

CI runs the unit checks (normally and under `python -O`) and the smoke presets of
C1, C4, C5 and C6; C2, C3 and C7–C10 need the diffsol binary, the reference binary
of another repository or a large archive and are run by hand.

Timing is descriptive and depends on host load. Independent training and
validation seeds and native work counts establish the numerical protocol, not
statistical power for a scientific discovery. Clone inference, learned
topology proposals, stock MO/QD workers and generated-kernel performance
studies are follow-up work.
