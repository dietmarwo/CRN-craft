# Direct SSA

**Status:** implemented in M2 (`0.1.0-dev.4`).

Gillespie's Direct method for a time-homogeneous, well-mixed count process.
Existing mass-action/Hill/MM propensities and availability gates are reused;
count changes are checked and atomic. The solver recomputes ordered rates,
uses a positive exponential wait and a separate categorical draw, and treats
zero total propensity as absorption. It includes events exactly on save times.

Random streams: pinned rand_chacha 0.3.1 / rand_core 0.6.4, root-seed key and
replicate-index stream; owned 53-bit mappings and libm 0.2.15 inversion.
No clipped rates, lost positive channels or last-channel selection fallback.
CRN-E042 underflow and CRN-E061/E062 resolution/time errors remain explicit.

Completed fixed-work outputs use ordered reduction and are checked at 1/2/8
threads. Incomplete ensembles have accounting and partial retained paths,
but no unconditional mean/variance. Limits and cancellation do not prove
non-explosion. Exact-method sampling still has Monte Carlo and floating-point
limitations. Cross-platform numerical identity and comparative throughput
are not established by this iteration.

See [Direct SSA details](../direct-ssa.md), [archive contracts](../archive-contract.md)
and [M2 acceptance targets](../decisions/iteration-04-validation.md).
Reference: [Gillespie (1977)](https://doi.org/10.1021/j100540a008).
