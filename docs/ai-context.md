# Agent context: CRN-craft M0

CRN-craft models well-mixed chemical reaction networks. This version supports **validation and capability discovery only**. No command runs SSA, FSP, ODEs, optimization or numerical network analysis.

Start by running `crn-craft describe`; use `describe model` and `describe run` for canonical schemas. Select an example from `models/examples`. After editing a model, compute SHA-256 on its exact bytes and update its run reference. Keep the model in or below the run directory. Run `check MODEL RUN` and inspect every diagnostic; zero exit status means input preflight passed, not that a simulation ran.

Identifiers use ASCII letters/digits/underscore, beginning with a letter or underscore, up to 64 bytes. Species/parameter names must be disjoint. Reactions have positive integer stoichiometry entries; omit zero entries. Parameter references are exact strings. Mass-action coefficients are stochastic count coefficients; Ω is not an automatic rescaling instruction. Use counts and seconds. Hill/MM rates are complete propensities. Read docs/model-cards/rate-laws.md before translating another simulator's input.

Read docs/diagnostics.md for codes. Paths are RFC 6901 pointers; `""` is the document root. A syntax/type error normally stops that document's parsing. Semantic checks collect multiple errors in stable order. Unavailable methods must not be simulated by invented results.

Useful current prompt:

> Create a two-species conversion model by adapting an example. Describe the count-rate convention, bind its run file to the model bytes, run the preflight checker, and explain the diagnostics. State explicitly that no numerical simulation has been performed.

Next implementation prompt:

> Implement M1a from docs/implementation-plan-2026-11-01.md. Compile the model into indexed stoichiometry and dependencies, implement checked propensity evaluation and reactant-availability gating, and add analyze. Test homodimers and heterodimers against independently calculated values and compare 200 seeded small networks with an independent oracle. Deliver tests, evidence and a log; do not add SSA before this gate passes.
