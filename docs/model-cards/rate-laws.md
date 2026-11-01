# Rate-law contracts

**Status:** argument validation implemented in M0; numerical evaluation planned for M1.

| Kind | Complete propensity a(x) | Argument domains |
|---|---|---|
| mass_action | c × product binomial(xᵢ, νᵢ) | c ≥ 0; finite |
| hill_activation | basal + vmax × xⁿ/(kⁿ+xⁿ) | basal,vmax ≥ 0; k,n > 0; finite |
| hill_repression | basal + vmax × kⁿ/(kⁿ+xⁿ) | basal,vmax ≥ 0; k,n > 0; finite |
| michaelis_menten | vmax × x/(km+x) | vmax ≥ 0; km > 0; finite |

All custom laws denote full channel propensities. Reactant availability is a separate gate; no reaction may consume unavailable particles. Counts and regulatory thresholds are in molecule counts, time in seconds. Hill/MM are effective models, not elementary reaction mechanisms. The explicit enzyme example uses three elementary mass-action channels and is distinct from the reduced Michaelis–Menten example.

At fixed Ω, stored SSA coefficients do not rescale automatically. The mass-action ODE/Ω convention is specified in the revised plan. M1 acceptance includes exact `2A` and `A+B` small-count checks, independently evaluated values, parameter-domain checks, unavailable reactants and finite/overflow behavior. No such numerical evaluation is claimed in M0.
