# Julia comparison: protocol declared before measurements

Date: 2 November 2026. CRN-craft base: baed0e4590f93e4c1dcc77cd7f039a2604a1594c.

Reuse W1–W4 from `validation/benchmark_campaign.py`, as in the GillesPy2
comparison: 64 paths, horizon 100, saves 0 and 100, one numerical worker.
Build the current release CLI and Python binding from this exact base.
Install official Julia 1.13.3, verify its published SHA-256, resolve registered
Catalyst/JumpProcesses into an isolated environment and retain its Manifest.

Compare CRN-craft Direct and cached Direct with Catalyst-generated JumpProblems
using JumpProcesses Direct, SortingDirect, DirectCR, NRM and RSSACR, all through
SSAStepper. Do not use a hand-written replacement SSA. Catalyst must generate the
mass-action specialization and dependency graphs. Custom rates are full
propensities. For W4 only, omit redundant reactant gates on the proven invariant
class X0=...=X7, as in the GillesPy2 comparison. No approximate solver is timed.

Before performance interpretation, check translated stoichiometry and propensity
values against CRN-craft, including zero/low counts, Hill rates, equal-state W4
and a separate homodimer convention control. Relative tolerance is 2e-12.

Validate 512 independent paths per workload and Julia algorithm plus CRN-craft.
Check all 34 W1/W2 species against their exact final Poisson laws; compare all
six W3 species and the single independent W4 coordinate with CRN-craft using
two-sample empirical CDF distance. A conservative family alpha=0.01 covers
34*6 + 7*5 = 239 comparisons. Use union-bound DKW envelopes, without assuming
independent species: sqrt(log(2*239/0.01)/(2*512)) for one-sample checks and
twice this bound for two-sample checks (union bound over both empirical CDFs
uses log(4*239/0.01) instead of log(2*239/0.01)). Means/z-scores are descriptive.
All runs must finish at t=100 with nonnegative integer counts; W4 equality must
hold. Retain failures and do not loosen thresholds after seeing results.

Report Julia package loading, model/problem construction and first execution
separately from warm execution. Run five timing rounds, rotate/interleave the
seven implementations, and keep all raw times, medians, ranges and numerical
digests. Builds and heavy checks finish before timing. Reset each Julia timing
batch from an unsolved problem copy outside the timer, because SortingDirect
retains its learned ordering. Use independent Julia and CRN RNG streams, seed
4242 for timing and 777 for validation; do not expect cross-engine trajectories
to match. Warm-up must not alter the measured initial state.

Primary comparison: warm Python native-API call versus Julia batch solve plus
ordered mean/variance reduction, both including per-trajectory solve setup.
Describe serialization and object-creation differences. Julia process startup,
imports, model generation, deepcopy of pristine batch problems and IPC JSON are
outside its warm timer. Include cold costs separately. If actual Julia event
counts cannot be established without perturbing the solver, label CRN-count
normalization as a proxy and keep ensemble seconds primary.

Five timing repetitions are a small shared-host sample, not a universal ranking.
Show every Julia algorithm and the fastest observed option per workload; do not
compare only an unfavorable Julia baseline. These four small/medium models do
not establish behavior at thousands of reactions or on multiple physical cores.
If a package or algorithm fails, retain the failure and qualify coverage.
