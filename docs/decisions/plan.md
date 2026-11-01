# M0 design decisions

1 November 2026. This delivery follows the selected CRN-craft name, Rust core, later minimal Python interface, and M7 → M8 → demand-driven M9 order.

- First iteration is M0; no SSA or numerical analysis is represented as implemented.
- Local source/ZIP delivery is sufficient; publishing or creating a remote is not a prerequisite.
- Parameter and species names are disjoint to keep later expressions/exporters unambiguous.
- Use tagged rate/method variants with unknown fields denied; custom laws are full propensities.
- Separate safety population limits from FSP truncation bounds. Finite truncation can approximate an unbounded physical state space.
- Explicit save times avoid ambiguous cumulative step rounding.
- Model identity is exact raw-byte SHA-256; normalized representations never replace its meaning.
- `check_files` also verifies path identity; `check_bytes` cannot verify a filesystem reference.
- Defer expression AST support until M1b needs it. JSON trees cannot contain cycles; test unknown reference syntax and depth limits instead.
- Reserve a result envelope without emitting or claiming verified results. Version its full implementation in M2.
- Advance a useful Direct-based circuit workflow after M3/M4a without requiring all large-network aggregators first.
- The default installed toolchain is Rust 1.88.0; the explicit stable alias resolves to Rust 1.99.1. Test and record them separately, alongside MSRV 1.85.0. Hosted CI is configured but is not claimed as executed here.
