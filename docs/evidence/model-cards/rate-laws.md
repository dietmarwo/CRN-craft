# Rate-law contracts

> **History (note of 12 November 2026).** This card describes the rate laws
> of the network model `crn-craft.model/0.1` as the network engine of
> releases 0.1.0 to 0.6.0 evaluated them. Since release 0.7.0 such a model is
> read through the converter: mass action becomes a rule `rate`, and the Hill
> and Michaelis–Menten laws become whole-rule `propensity` expressions in a
> stable form ([architecture](../../developer/architecture.md#network-models-and-the-converter)).
> A negative or non-finite rate is `CRN-R030`; `CRN-E042` is retired.

**Status:** argument validation and numerical evaluation implemented in M1a (`0.1.0-dev.2`). No stochastic time integration is implemented.

| Kind | Complete propensity a(x) | Argument domains |
|---|---|---|
| mass_action | c × product binomial(xᵢ, νᵢ) | c ≥ 0; finite |
| hill_activation | basal + vmax × xⁿ/(kⁿ+xⁿ) | basal,vmax ≥ 0; k,n > 0; finite |
| hill_repression | basal + vmax × kⁿ/(kⁿ+xⁿ) | basal,vmax ≥ 0; k,n > 0; finite |
| michaelis_menten | vmax × x/(km+x) | vmax ≥ 0; km > 0; finite |

All custom laws denote full channel propensities. Reactant availability is a separate gate; no reaction may consume unavailable particles. Counts and regulatory thresholds are in molecule counts, time in seconds. Hill/MM are effective models, not elementary reaction mechanisms. The explicit enzyme example uses three elementary mass-action channels and is distinct from the reduced Michaelis–Menten example.

At fixed Ω, stored coefficients do not rescale automatically. The mass-action ODE/Ω convention is specified in the revised plan. A `2A` channel uses `c*x*(x-1)/2`, including its factorial divisor. Empty reactant sides have an empty product of one. If any reactant is unavailable, every law returns zero **before** evaluating its formula. Catalysts with equal reactant/product coefficients still gate a reaction.

## Numerical implementation and limits

The core pins [libm 0.2.15](https://docs.rs/libm/0.2.15/libm/) for pure-Rust scaling and transcendental operations. Mass action combines the coefficient and binomial factors in a scaled mantissa/exponent product, so a large binomial need not overflow before multiplication by a tiny coefficient. Hill fractions use logarithms and a stable logistic form, with `log1p` near x=k. Very small fractions are combined with log(vmax) before exponentiation. Michaelis–Menten evaluates `vmax * (x/(km+x))`, avoiding an overflowing `vmax*x` intermediate. Basal terms are added after the induced Hill term.

Evaluated propensities are f64. Overflow, non-finite/negative values and an entirely underflowed strictly positive rate fail with CRN-E042. Exact zero coefficients and unavailable-reactant gates return zero. Subnormal positive rates are allowed. Normal rounding can discard a tiny induced contribution when a positive basal term remains. Integer counts and updates stay exact i64; conversion to f64 for rate formulas can round counts above 2^53. This is not a correctly rounded or cross-platform identity guarantee.

Validation targets were fixed in `docs/decisions/iteration-02-validation.md`. Tests cover exact homodimer/heterodimer values for counts 0–50, high-count/tiny-coefficient combinations, extreme Hill arguments, steep thresholds, MM overflow avoidance, subnormals, all gates and transactional failures. An independent Python oracle uses integer combinations, Fraction linear algebra and 100-digit Decimal formulas, including 200 seeded networks. See `validation/iteration-02/RESULTS.md` for results and the retained development failure.

For interoperability, [Catalyst's API documentation](https://docs.sciml.ai/Catalyst/stable/api/core_api/) describes combinatorial-rate conventions. An exporter must explicitly translate its chosen target convention; matching a law's name alone does not establish coefficient equivalence. No external simulator code is copied or linked here.
