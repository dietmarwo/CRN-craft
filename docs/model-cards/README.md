# Model/algorithm card template

Every numerical extension should record:

- Name, implementation status and supported method/version.
- Mathematical definition and primary references.
- Assumptions, units, domains, reactant gating and unsupported cases.
- Numerical strategy, overflow/dynamic-range handling and error reporting.
- Seed/reproducibility scope where relevant.
- Independent fixtures, predeclared tolerances and mutation tests.
- Performance workload and hardware when measured.
- Known limitations and untested conditions.

The two initial cards document **contracts and plans**. They do not imply implemented rate evaluation or SSA.
