# M4b: dependency-cached Direct

`simulate-cached MODEL RUN NEW_DIRECTORY` and `simulation::execute_cached` cache
only propensities unaffected by the last reaction. M1's reviewed dependency graph
includes custom-rate regulators and reactant gates. Rates are still summed and
selected in the original declaration order, so pinned random words and every
completed numerical payload agree with ordinary Direct. No run/result/archive
contract changes are required for this exact implementation variant.

The numerical/library/archive/CLI tests pass. Additional tests compare all four
model families and 1/2/8 workers, and exercise a cooperative integration checkpoint.
The seven paired repetitions in m4b-cached-benchmark.json require identical
numerical digests before comparing time. Median core speedups over main 5a07334:
1.11x two channels, 3.10x sparse 66, 2.35x circuit 12, 1.38x dense 66. These are
shared-host observations with raw ranges, not general complexity guarantees.

The external checkpoint hook invokes the existing worker's cancellation/deadline
check at the amortized polling stride. Ordinary Direct remains the default and
its published sampler/archive semantics remain unchanged.
