# How CRN-craft compares

This page compares CRN-craft 0.7.0 with GillesPy2 1.8.3 and rebop 0.10.0, the
Python packages a newcomer most likely knows, and names a few other tools.
Statements about other tools were checked against their documentation,
package pages or source on 13 November 2026 ([Sources](#sources)). Speed
numbers come only from measurements committed to this repository.

## In short

CRN-craft fits a study that needs many replicates per run, reproducible
ensembles, long exact runs, a model worth compiling, species with attributes,
or a continuous variable that rates read. Another tool is the better choice
for full SBML with events and rules, ODE or hybrid ODE–stochastic solvers,
spatial models or a graphical interface (GillesPy2 with StochSS), rule-based
models with very many states (BioNetGen and NFsim), or a model compiled into
your own Rust program (rebop's macro). GillesPy2 and BioNetGen have long
histories; CRN-craft's first release was in November 2026.

## Features

| | CRN-craft 0.7.0 | GillesPy2 1.8.3 | rebop 0.10.0 |
| --- | --- | --- | --- |
| Model input in Python | a builder: `add_species`, `add_parameter`, `add_reaction`, `set_stop`, `add_ode`; a rate is a number (mass action) or an expression string ([chapter 2](guide/02-first-model.md)) | `Model` with `Species`, `Parameter` and `Reaction`; a rate constant (mass action) or a propensity string | `Gillespie().add_reaction(rate, reactants, products)`; a number (mass action) or an expression string, with parameters passed to `run` |
| Model input from files | its own JSON model files; a strict subset of SBML Level 3 Version 2 and of BioNetGen `.net` that refuses what it cannot represent exactly ([chapters 11](guide/11-analysis-import-export.md#import-sbml-and-bionetgen), [14](guide/14-model-files.md)) | SBML import and export (needs libSBML); JSON | none; SBML is on its list of planned features |
| Exact methods | Static Direct with five selectors (`linear`, `composition_rejection`, `sorting`, `rssa_cr`, `auto`), next reaction, a tree method for species that appear and disappear ([chapter 7](guide/07-methods.md)) | the direct method, in C++ (`SSACSolver`) and NumPy (`NumPySSASolver`) | one SSA; other SSA algorithms are planned |
| Approximate methods | tau leaping and binomial splitting; the sampler name of the result says "approximate" | tau leaping (C++ and NumPy); chemical Langevin (`CLESolver`) | none; tau leaping is planned |
| Deterministic and hybrid | `hybrid` for models with continuous variables; no ODE solver for the network in the wheel ([chapter 9](guide/09-continuous-variables.md)) | ODE solvers (C++ and NumPy); `TauHybridSolver`, with species that switch between discrete and continuous | none; hybrid models are planned |
| Exact distributions | finite state projection with separate error terms ([chapter 7](guide/07-methods.md#exact-distributions-with-finite-state-projection)) | not documented | not documented |
| Ensembles and parallelism | one call for any number of replicates on native threads; statistics and quantiles as pandas data frames ([chapter 3](guide/03-running.md#many-trajectories-in-one-call)) | `number_of_trajectories` in one `run`; the C++ SSA runs them one after another in one process; ensemble mean and standard deviation | one trajectory per `run` call; parallelization is planned |
| Reproducibility | one ChaCha8 stream per replicate, results independent of the thread count, model and request named by SHA-256, sealed archives with a verifier ([chapters 8](guide/08-reproducibility.md), [10](guide/10-command-line.md)) | a `seed` argument | a seed or NumPy generator per call |
| Compilation from Python | `crn.compile_loop(model)`, opt-in; needs a Rust toolchain; cached; byte-identical results ([chapter 5](guide/05-speed.md#compiled-models)) | the C++ solvers compile each model with g++ and SCons; NumPy solvers otherwise | the Python API compiles nothing; the Rust macro `define_system!` compiles a mass-action model into a Rust program |
| Species with attributes | in model files, with new identities at run time; the builder cannot write them yet ([chapter 14](guide/14-model-files.md)) | not documented | not documented; rule modeling is listed as "probably not to come" |
| Continuous variables | ODEs between events that rates read; independent variables integrated once per run | rate rules and assignment rules (hybrid solver) | not documented |
| Event-level output | stop conditions, trackers (first passage, extremes, integrals), bounded event logs, Python callbacks ([chapter 6](guide/06-stops-trackers-logs.md)) | SBML-style events with triggers, delays and assignments; live output during a run | every reaction with `nb_steps=0`; events are listed as "probably not to come" |
| Results | pandas and xarray; mean, variance, standard error, quantiles; markdown and JSON reports ([chapter 4](guide/04-results.md)) | trajectories as NumPy arrays; matplotlib and Plotly plots; CSV | an xarray dataset per run |
| Languages | Python, a command line, Rust; R through reticulate ([R](using-from-r.md)) | Python | Python and Rust |
| Numeric rate on `A + A` | k·n(n−1)/2 ([model reference](reference/model-language.md#counting-identical-reactants-and-rebop)); exports GillesPy2 source with the constant halved | k·n(n−1) at volume 1 | k·n(n−1) |
| License | MIT | GPL-3.0 | MIT |
| What it has that CRN-craft lacks | | SBML events and rules, ODE, hybrid and Langevin solvers, the StochSS web interface, spatial models through the sibling project SpatialPy | a Rust macro that compiles a model into your program |

## Two things to try first

**Many trajectories as one ensemble.** `crn.ensemble(model, ...).simulate()`
runs any number of replicates on native threads and returns their means,
variances and quantiles; you write no loop ([chapter 3](guide/03-running.md#many-trajectories-in-one-call)).
Results do not depend on the thread count; two threads ran 2,000 mtDNA cells
2.1 times as fast as one (release 0.2.0, see Speed).

**A compiled model while you stay in Python.** `loop = crn.compile_loop(model)`
builds a native event loop once and caches it; `run.simulate(compiled=loop)`
gives byte-identical results ([chapter 5](guide/05-speed.md#compiled-models)).
On the mtDNA model it ran 2.14 times as fast as the interpreter (12 November
2026). It needs a Rust toolchain; the gain depends on the model.

## Speed

All runs: a shared two-CPU cloud host, under a CPU lock. Timings vary (most
GillesPy2 runs within 10% of their median, single runs up to 34%; rebop
repeats up to 1.6 times), so this page quotes ratios: other tool over CRN-craft.

**GillesPy2, current engine** (12 November 2026; unified engine of 0.7.0 at
`cec91ed`, whose build reports 0.6.0; GillesPy2 1.8.3, Python 3.12.3, g++ 13.3;
64 paths to t = 100, one thread, medians of five runs). `SSACSolver` against
CRN-craft's default on the networks W1 to W4: 2.9, 10.6, 3.6, 9.3; against the
linear scan with inverse waits: 2.8, 5.6, 3.7, 1.3. `NumPySSASolver`: 134 to
358. With the selectors `composition_rejection` and `sorting`, which `auto`
does not pick there, CRN-craft was slower on W4 (0.88, 0.93). The C++ compile
(2.1 to 2.9 s per model) is excluded; 28 of 28 agreement tests passed. Not measured:
tau, ODE and hybrid solvers, threads ([evidence](evidence/comparison-gillespy2.md)).

**rebop, Python API, current engine** (13 November 2026; crn-craft 0.7.0 and
rebop 0.10.0 from PyPI, Python 3.12.3, Rust 1.97.0; 64 paths to t = 100, one
thread, medians of three repeats; 8 of 8 agreement tests passed). W1, two
reactions: 1.21 against the default, 1.50 against the compiled loop; these
two are not resolved against each other. W2, 66 reactions: 25.0 against the
default, 13.4 against the linear scan; the compiled loop, linear scan only, ran
0.95 times as fast as the interpreted one. Not measured: rebop's macro,
W3 and W4 (rates that are not mass action), threads
(`validation/merge/rebop-comparison-0.7.json`; [`RESULTS.md`](../validation/merge/RESULTS.md), "rebop on 0.7.0").

**rebop's macro, removed engine** (2 November 2026; CRN-craft 0.2.0.dev1, the
network engine that 0.7.0 removed; rebop 0.10.0; 54 optimization runs).
**rebop was faster:** its macro beat CRN-craft by 3.44 times on birth–death and
1.72 on binding12, and tied on sparse33 (0.94 to 1.00); with one build the
totals differed by 2.5%. rebop's Python expression route took 15 times as long
as CRN-craft on birth–death. Not repeated on 0.7.0 ([evidence](evidence/comparison-parameter-optimization.md)).

**Compiled models** (12 November 2026). mtDNA three-class model file, 200 cells, two threads:
compiled 2.14 times as fast as interpreted Static Direct, `mtdna_direct` 2.84 times
(`validation/merge/manual-ratios/timings.json`).
Twelve networks, linear scan: compiled over interpreted time 0.61 to 0.90 on
eight, 1.15 on W7; W4 to W6 refused (source over 1 MiB). First build 49.3 s,
later 0.1 to 6.4 s ([`RESULTS.md`](../validation/merge/RESULTS.md), "Compiled event loops").

**Threads** (release 0.2.0, `mtdna_direct`, 2,000 cells): two threads 2.1 times as fast as
one; 20 threads on eight CPUs 5.8 to 6.1 times ([analysis](evidence/mtdna-performance-analysis.md#3-threads)).
**Julia** (Catalyst.jl, JumpProcesses.jl): compared only with the 0.6.0 engine; Julia was not on the
0.7.0 host ([comparison](evidence/comparison-catalyst-jumpprocesses.md), [ten iterations](evidence/comparison-julia-gap10.md), [eight threads](evidence/comparison-julia-8-threads.md)).

## Other tools

**StochPy** 2.5 (BSD-3-Clause) has the direct, first reaction and next reaction
methods and tau leaping, reads its own `.psc` files and SBML (optional
libSBML), and ships DSMTS models with a test call. **BioNetGen and NFsim** are
the reference for rule-based models; NFsim simulates models with very many or
unbounded states without generating the network. They are the nearest
relatives of species with attributes, and CRN-craft imports a subset of
expanded `.net` files ([chapter 11](guide/11-analysis-import-export.md#import-sbml-and-bionetgen)).
**JumpProcesses.jl** (Julia, SciML) provides SSAs and couples jumps to ODEs and
SDEs. There is no measurement against StochPy or NFsim.

## Where CRN-craft is weaker

- No deterministic ODE solver for a network in the wheel; no spatial models;
  no graphical interface.
- SBML and BioNetGen import cover a strict subset: no events, rules, fast
  reactions or non-integer stoichiometry.
- A compiled model needs a Rust toolchain, runs the linear scan only, refuses
  sources over 1 MiB, and can be slower than interpreted `rssa_cr`.
- Tau leaping is far slower than Static Direct on networks with about a
  thousand reactions; `bind_parameters` prepares the model again per candidate.
- The builder cannot write species with attributes yet. In the one measured
  optimization, rebop's macro was faster. First release: November 2026.

## Sources

Checked on 13 November 2026. CRN-craft cells: run or read on the 0.7.0 wheel.

- **GillesPy2 1.8.3** (19 September 2024, current): [PyPI](https://pypi.org/project/gillespy2/),
  [documentation](https://gillespy2.readthedocs.io/), [README](https://github.com/StochSS/GillesPy2);
  source of the 1.8.3 package: `core/model.py` (SBML, `get_best_solver`), `core/reaction.py` (`A + A`),
  `core/results.py`, `solvers/cpp`, `solvers/numpy`, `ssa_cpp_solver/SSASolver.cpp` (trajectory loop).
  StochSS and SpatialPy: [stochss.org](https://stochss.org).
- **rebop 0.10.0** (24 August 2026): [PyPI](https://pypi.org/project/rebop/) (README with planned
  and "probably not to come" features), [crates.io](https://crates.io/crates/rebop); source of the
  0.10.0 package: `gillespie.py` (`add_reaction`, `run`, the `A + A` convention).
- **StochPy 2.5** (18 June 2025): [PyPI](https://pypi.org/project/stochpy/). **BioNetGen, NFsim:**
  READMEs of [bionetgen](https://github.com/RuleWorld/bionetgen) and [nfsim](https://github.com/RuleWorld/nfsim).
  **JumpProcesses.jl:** [README](https://github.com/SciML/JumpProcesses.jl).
