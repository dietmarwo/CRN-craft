# Rules engine K0/K1: implementation log

Date: 2 November 2026. Branch `rules-k0-k1`, based on `reusable-parameters`
(Astra's parameter-binding import, under review separately). Protocol:
`docs/decisions/rules-k0-k1-2026-11-02.md`. Motivation: Axel Kowald's
mitochondrial mutant-accumulation models, where every mutant is a new clone,
and his request for an engine that is not tied to that one model.

## What was built

- **Crate `crn-craft-rules`** (MIT, no unsafe code, depends only on the core
  crate's strict parser and diagnostics plus pinned workspace dependencies).
  - `model.rs`: the four contracts and the compiler. Every name and
    expression is checked before execution; patterns are deduplicated into
    classes by comparing compiled expressions after label resolution, so
    `a.size` and `b.size` count as identical factors.
  - `expr.rs`: a small Pratt parser and evaluator (arithmetic, comparisons,
    short-circuit logic, `min max exp log sqrt abs floor ceil if`; random
    draws only in product attributes). Scopes are enforced at compile time.
  - `tree.rs`: growable pairwise sum trees with an optional minimum-positive
    track for the 2⁻⁵³ resolution check, pair trees (pattern weight, pair
    weight and smallest factor per node) for identical-pattern pairs, and max
    trees for `max`/`min` observables.
  - `engine.rs`: one replicate of the direct method. Per type, an arena of
    instances keyed by attribute bits (equal products merge, empty instances
    are freed); per pattern class, a sum tree of h·n or, for identical pairs,
    a pair tree; per observable, a sum or max tree. All rule propensities are
    recomputed after each event; the global rate is evaluated only when the
    reactant weight is positive. Products equal to a reactant cancel its
    consumption.
  - `lib.rs`: run request validation and admission limits, parallel
    ensembles over keyed ChaCha8 streams (domain `crn-craft.rules/0.1`)
    reduced in replicate order (Welford statistics, only kept series stored),
    records, output budget, file entry points, schemas and the discovery
    entry.
- **CLI**: `rules-check MODEL`, `rules MODEL REQUEST [--cancel-file PATH]`,
  `describe rules-model|rules-check|rules-run|rules-result`. The native CLI's
  `describe` reports `describe/0.8` with a `dynamic_species` entry; the core
  library and the Python package, which has no rules commands, stay at 0.7.
- **Models** in `models/rules/`: immigration–death, immigration–dimerization,
  Kowald's fixed-population mtDNA model, aggregation with constant and with
  multiplicative kernel.
- **Validation**: 21 Rust tests in the crate (parsing and scope errors, strict input,
  merging and `new_id`, mass conservation, `distinct`/`max`/`min`, stop
  conditions, every budget, cancellation, numerical failures including the
  resolution check, pair semantics, cross-type rules, thread-count replay,
  generated schemas, exact pair selection, the review regressions);
  `validation/check_rules.py` with the declared exact-law
  gates; contract checks extended to the four schemas and the five models; a
  new CI job runs the gates normally and under `python -O`.
- **Docs**: `docs/rules-language.md` (language reference), README section,
  architecture table, diagnostics, changelog, results.

## Results

All 47 checks passed in both gate runs (`validation/rules/RESULTS.md`): the 30
declared marginals are within the DKW bound 0.00742 (largest distance 0.00377
with sampler version 1, 0.00342 with version 2); each fixture's wrong law and
the supplementary "natural denominator" mtDNA law are rejected; the network
engine agrees with the same dimerization law; one and four threads give
identical ensembles. Both reports are identical under `python -O`.

Rough single-thread speed on the 2-vCPU container (not a benchmark), sampler
version 2: Kowald's model at K = 1000 runs about 4.3 million events per
second; aggregation of 10,000 monomers about 3.9 million (constant kernel) and
1.8–2.4 million (multiplicative kernel, including past gelation).

## Independent review and second gate run

A separate agent reviewed the engine without seeing these conclusions. It
confirmed the sampled law (pair, cross-type, cancellation, slot reuse, report
times, stop, budgets, thread independence; about 55 strict-input cases) and
found, each reproduced through the CLI:

1. Memory grew with replicates × report times × observables (all series were
   kept for the statistics). Now: ordered streaming reduction with a bounded
   window of replicates in flight, Welford statistics, only kept series stored,
   and an admission check (2 GiB reservation, output lower bound) before any
   run.
2. A model with a large attribute and no `new_id()` failed on its first event
   ("new_id() exceeds 2^53"). Now the limit applies only to models that call
   `new_id()`, at compile time and when an identity is handed out.
3. At 2⁵³, `new_id()` silently returned the same identity. Now an error.
4. An overflowing pair propensity became 0 (NaN through `max`), and
5. an overflowing observable total became `null` in a "completed" result.
   Both now fail with `CRN-R030`.
6. `new_id()` can collide with identities written by product constants; the
   protocol's wording promised more than the code does. Documented.
7. A duplicate-label diagnostic had an escaped path. Fixed.

It also showed that ½[(Σ h n)² − Σ h² n] loses accuracy when one instance
dominates and that the rejection loop then needs about h_max/2 tries per
event. Identical pairs now use pair trees with exact selection and no
rejection (sampler version 2). Smaller items: NaN conditions are errors,
`min`/`max` keep NaN, `choice` evaluates only its chosen branch, model names
are identifiers, tokens are separated by ASCII whitespace only, and the docs
explain `rng_words` past the budget and `peak_instances` at the limit.

The declared gates were rerun unchanged on version 2: 47/47, identical under
`python -O`. The first report is kept as `k1-gates-sampler1.json`.

## Decisions during implementation

Recorded in the protocol before the gate run: the global factor is the field
`rate`; observables gained `distinct`, `max` and `min`; a `rules-check`
contract; products equal to a reactant cancel its consumption; identical
pairs were first drawn by rejecting the same molecule (replaced after the
review, see above); the wrong laws were chosen from exact-law differences
only.

## Problems found and fixed during development

- `created_instances` counted a kept reactant twice (removed, then re-created
  by its `copy` product). Fixed by cancelling equal reactant/product pairs,
  which also saves the re-evaluation of factors.
- Test expectations: unknown fields are `CRN-E003`, not `CRN-E001`; a stop
  test looked at a report time before the stop; the output-limit test now
  derives its budget from the measured result sizes.
- The first discovery change bumped the core `describe` to 0.8 and would have
  made the Python CLI advertise commands it lacks; the native CLI now adds the
  rules entry itself.
- The exact-law script first required Poisson weights to sum to 1 within
  1e-13, which rounding over 2,400 terms exceeds; it now checks 1e-9 (the
  truncated tail is below 1e-25).

## Checks run

`cargo fmt --all -- --check`; `cargo clippy --workspace --all-targets --locked
-- -D warnings`; `cargo test --release --workspace --locked`;
`validation/check_contracts.py` normally and under `-O` (identical, no
failures);
`validation/check_dependencies.py` normally and under `-O` (identical; the
crate is added to the allowed workspace packages); `validation/check_rules.py`
normally and under `-O` (47/47, identical); `ruff format --check` and `ruff
check` on the new script.

## Limitations

- MSRV 1.85 was not compiled locally (rustup could not reach
  static.rust-lang.org from this container); clippy's MSRV lint passed and CI
  covers 1.85.0 and 1.88.0.
- No time-dependent rates or ODE coupling, so Kowald's variable-population
  model with ATP is not expressible yet (K2). No compartments, no Python
  builder, no counts inside expressions (largest clone *size* needs them).
- Two different patterns over one type are rejected. The pair resolution
  check is conservative (smallest factor squared against the pair total) and
  can reject extreme factor ranges that would still be selectable.
- All rule propensities are recomputed per event: fine for a few rules,
  wasteful for hundreds. No benchmark against Kowald's Java code yet (running
  it needs his permission).
- Validated laws are small exact systems; larger systems are covered for
  accounting by unit tests, not for law.
