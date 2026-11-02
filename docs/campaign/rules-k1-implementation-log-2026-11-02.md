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
    track for the 2⁻⁵³ resolution check, and max trees for `max`/`min`
    observables.
  - `engine.rs`: one replicate of the direct method. Per type, an arena of
    instances keyed by attribute bits (equal products merge, empty instances
    are freed); per pattern class, trees of h·n (and h²·n for identical
    pairs); per observable, a sum or max tree. All rule propensities are
    recomputed after each event; the global rate is evaluated only when the
    reactant weight is positive. Products equal to a reactant cancel its
    consumption.
  - `lib.rs`: run request validation, ordered parallel ensembles over keyed
    ChaCha8 streams (domain `crn-craft.rules/0.1`), records, statistics,
    output budget, file entry points, schemas and the discovery entry.
- **CLI**: `rules-check MODEL`, `rules MODEL REQUEST [--cancel-file PATH]`,
  `describe rules-model|rules-check|rules-run|rules-result`. The native CLI's
  `describe` reports `describe/0.8` with a `dynamic_species` entry; the core
  library and the Python package, which has no rules commands, stay at 0.7.
- **Models** in `models/rules/`: immigration–death, immigration–dimerization,
  Kowald's fixed-population mtDNA model, aggregation with constant and with
  multiplicative kernel.
- **Validation**: 17 Rust tests (parsing and scope errors, strict input,
  merging and `new_id`, mass conservation, `distinct`/`max`/`min`, stop
  conditions, every budget, cancellation, numerical failures including the
  resolution check, pair semantics, cross-type rules, thread-count replay,
  generated schemas); `validation/check_rules.py` with the declared exact-law
  gates; contract checks extended to the four schemas and the five models; a
  new CI job runs the gates normally and under `python -O`.
- **Docs**: `docs/rules-language.md` (language reference), README section,
  architecture table, diagnostics, changelog, results.

## Results

All 47 checks of the gate run passed (`validation/rules/RESULTS.md`): the 30
declared marginals are within the DKW bound 0.00742 (largest distance
0.00377); each fixture's wrong law and the supplementary "natural denominator"
mtDNA law are rejected; the network engine agrees with the same dimerization
law; one and four threads give identical ensembles. The report is identical
under `python -O`.

Rough single-thread speed on the 2-vCPU container (not a benchmark): Kowald's
model at K = 1000 runs about 3.3 million events per second; aggregation of
10,000 monomers about 2.4–2.9 million.

## Decisions during implementation

Recorded in the protocol before the gate run: the global factor is the field
`rate`; observables gained `distinct`, `max` and `min`; a `rules-check`
contract; products equal to a reactant cancel its consumption; identical
pairs are drawn by rejection of the same molecule; the wrong laws were chosen
from exact-law differences only.

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
`validation/check_contracts.py` normally and under `-O` (78/78, identical);
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
- Two different patterns over one type are rejected. Identical pairs use
  rejection, which slows down when one instance dominates the weight.
- All rule propensities are recomputed per event: fine for a few rules,
  wasteful for hundreds. No benchmark against Kowald's Java code yet (running
  it needs his permission).
- Validated laws are small exact systems; larger systems are covered for
  accounting by unit tests, not for law.
