# C4 — Infer aggregation rates with the rules engine

Dynamic species require different semantics from a fixed reaction network.
This lesson uses the native `crn-craft run` command with `crn-craft.run/1.0`
requests, a small independent finite-state reference and a bounded scalar
search. It does not pass dynamic species through the Python builder or the
parameter-binding interface. See [rules.py](rules.py).

## A system small enough to check independently

Start with six monomers. An unordered pair of clusters merges into one cluster
of their combined mass. Under the constant kernel, each pair merges at rate
`k`; under the multiplicative kernel, the rate is `k * size_a * size_b`.
The total mass remains six and at most five aggregation events can occur.

Enumerate all integer partitions of six. Build the transition generator by
considering each unordered cluster pair, including multiplicities when clusters
have equal sizes. A SciPy matrix exponential gives the exact distribution of
cluster count and largest-cluster size at any requested time. Probability
conservation and nonnegative probabilities are independently checked. This
reference does not reuse CRN's sampling or existing uniformization checker.

Synthetic target distributions use the constant kernel with `k=1.2`, observed
at times 0.15 and 0.4. Fit a single log rate separately for the constant and
multiplicative models, with bounds 0.05 to 3. The loss is mean-square CDF error,
not just an error in mean cluster count.

## Execute the bounded workflow

Build the core CLI as described in [setup](README.md), then run:

```bash
python tutorials/optimization/run.py C4 --output /tmp/c4-smoke
python tutorials/optimization/run.py C4 --preset study --output /tmp/c4-study
```

Use `--binary /absolute/crn-craft` when necessary. Each native call has a
20-second outer timeout, at most ten events per replicate, sixteen live
instances, 10,000 sampling steps and a 16 MiB output limit. The full request,
model, model hash and result are retained for every call.

The smoke fit requests 16 calls per kernel with 64 training paths; the study
requests 64 with 256 paths. Final validation uses 1,024 or 8,192 independent
paths at held-out times 0.25 and 0.7. Read the actual call counts and total
native path/event work rather than multiplying only the requested budget.

## Separate implementation error from model mismatch

For every path, the mass observable must remain six. At the fitted rate, compare
each held-out empirical CDF with the independent law for **that same kernel**.
A predeclared Dvoretzky–Kiefer–Wolfowitz bound uses family error 0.001 over four
marginal checks (two observables at two times). A union bound does not require
the observables or times to be independent.

This check validates sampling of the selected model. A different number in the
summary, `exact_heldout_prediction_mse`, compares its exact predictions against
the generating constant kernel. Each kernel also gets an independent reference:
`exact_law_optimal_rate` and `exact_law_optimal_loss` minimize the same CDF loss
against the exact law of that kernel instead of finite paths (about 1.2 with
zero loss for the constant kernel, about 0.67 with loss 7.1e-3 for the
multiplicative one). A fitted rate far from its reference, or `best_is_start`,
means the search budget, not the kernel, decided the ranking. A wrong kernel can pass the simulator law
check while making poor predictions. The tutorial reports the observed
ranking; it does not make successful identification a condition for passing.

## Prompt and follow-up

> Infer an aggregation rate from synthetic distributions, then fit the wrong
> kernel. Verify mass and the finite-state law at held-out times under fixed
> resource limits. Explain separately whether the simulator is correct and
> whether the selected model predicts the observations.

Try fewer observation times or noisier empirical targets. A clone-dynamics
extension should start from the existing fixed-population fixture with its own
law and biological assumptions; it is not included in this aggregation lesson.
Time-dependent rates, ODE coupling and counts inside rule expressions remain
outside the current K1 feature set. Read [rules semantics](../../docs/rules-language.md)
and the existing [exact-law evidence](../../validation/rules/RESULTS.md).
