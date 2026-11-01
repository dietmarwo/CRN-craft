# M4a: same-host Direct baseline

The frozen workloads and paired timing harness are in `validation/benchmark_campaign.py`.
The baseline binary was built from GitHub main 5a07334, not from an older delivered
patch. Model generation, 64 streams, time horizon 100 and two output points are fixed.
Raw event counts, numerical digests, core times and CLI process times are retained.

`compare_references.py` compares complete same-size elementary ensembles against
rebop 0.10.0 dynamic Python API and GillesPy2 1.8.3 NumPy SSA, five paired repeats.
Reference streams differ; event counts are not instrumented there. Reference times
include Python result construction, whereas CRN core time includes ordered statistics
but excludes process/archive startup. These limitations are in each report row.

Baseline rebop/CRN core ratios: 0.387 for two channels (CRN is about 2.58x slower),
1.274 for 66 sparse channels (CRN is about 1.27x faster). These are workload-specific.
No Julia executable is installed, so Julia quality/cost is an unavailable row.
The container advertises nine virtual cores and a quota of eight CPU equivalents;
there is no evidence that eight dedicated physical cores are available.

Next measured bottleneck: Direct still evaluates every propensity per event.
M4b will preserve exact sum/selection order and cache only unaffected propensities.
