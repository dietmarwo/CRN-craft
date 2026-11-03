# Validation and implementation record

Implemented C1–C4 as simulator-owned Python/CLI tutorials. No core SSA, rule,
ODE or optimizer algorithm was changed. C3 adds a bounded reference example
inside fcmaes-rust's existing oscillator workspace; its backend is explicitly
the local ReBop compatibility implementation.

The [verification record](results/verified.json) records dependency versions,
executed gates and source hashes. Full raw command artifacts, study summaries,
logs and source-base identities are included in the combined tutorial bundle.
Native CRN/fcmaes/SQP bindings were built from the recorded sources; LP-IP-craft
0.29.0 used its published native wheel. Builds and reports are host observations,
not release certification.

| Check | Result |
|---|---|
| C1 smoke and study, fresh/bound traces and ordinary replay | Passed; study has 12 checked-policy arms |
| C2 independent integration, log-gradient and rank audits | Passed; gradient error about 9.1e-10, tightened state error about 3.2e-12 |
| C3 generator, common backend scoring, MODE/archive and held-out paths | Passed; smoke/study archive occupied 3/2 of 16 cells |
| C4 mass and four-marginal law gate per kernel | Passed; largest study CDF distance about 0.0114 versus bound 0.0234 |
| Eight independent/negative Python unit checks | Passed normally and under `python -O` |
| Core fmt, clippy with warnings denied, release workspace tests | Passed with locked dependencies |
| Schema/CLI contract checker | 78 passed normally and under `python -O` |
| Core dependency gate | Passed |
| Command-job execution and optimized-Python smoke route | Passed for all four chapters |

The smoke and larger study presets were both executed. C1's optional fast arm
has its own recorded check; it is never combined with the checked-policy timing
arm. C3's full larger preset was also run under `python -O` after selecting the
final declared grid. Its archive is intentionally sparse: completion does not
establish that every niche contains a robust oscillator. The Vilar control is
only a short, separately identified conservation test.

Timing is descriptive. The larger teaching runs overlapped repository checks;
their wall times must not be used as a performance benchmark. Independent
training/validation seeds and native work counts establish the numerical
protocol, not statistical power for a scientific discovery. Optional clone
inference, learned topology proposals, new stock MO/QD workers and generated
kernel performance studies are follow-up work.
