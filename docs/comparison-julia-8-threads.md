# CRN-craft fast mode versus Julia: eight threads

> **Measured before the 0.2.0 release.** These measurements come from
> development builds of 2 November 2026 and have not been repeated on
> release 0.2.0 (4 November 2026). The network SSA engine changed in
> between: `fast` became the default of `ssa` requests
> (`crn-craft.ssa-run/0.3`), the Ziggurat waiting-time sampler was
> refactored so the rules engine can share it, and
> `crn-craft.ssa-run/0.4` added the `next_reaction` aggregator, which
> this comparison does not include. Requests of the measured versions
> return the same results as before; their timings on 0.2.0 may differ.

> **Import note (2 November 2026).** This report was integrated as delivered.
> Its harness (`validation/compare_eight_threads.py`, `validation/eight-thread/`)
> and raw evidence (`validation/eight-thread/evidence/comparison.json`) were not
> part of the delivery; the hashes under *Provenance* identify them. On import,
> every summary median, range and 1→8 ratio (224 rows) was recomputed from the
> individual observations, as were the selected methods, the fast/Julia ratios,
> the 6/7 and 3/7 counts, the geometric means 0.69× and 1.04× and the client
> latency ratios; all matched. The one-worker comparison at the end is derived
> from the same observations on import.

Measured 2 November 2026 (UTC), on the same host. CRN-craft GitHub main `c16fcc55654d0867a58af48d48a8e860efce8896`; numerical engine unchanged.

At eight numerical workers, across seven equally weighted workloads, the geometric mean of **fast-mode time / Julia time** is **0.69× for the original batches** and **1.04× for 256-path batches**. Fast mode has the lower best median on 6/7 and 3/7 workloads respectively. A ratio above 1 means Julia is faster; below 1 means CRN-craft is faster. This compares CRN's prepared Python API with Julia's prepared ensemble solve; the Julia boundary omits some result-contract work. These are fresh measurements, not extrapolations from serial results.

Each workload selects the lowest median among ten CRN fast configurations (five aggregators × inverse/Ziggurat waiting times) and five Julia algorithms. This describes the best observed choices on these fixtures, not an automatic method selector or a universal speedup. Lower-median counts include near ties; small differences are within the observed timing variation. All measurements are retained below.

## Original batches

Times are batch milliseconds, median of five warm observations. The 1→8 column compares the **same selected algorithm**, batch and seed at one and eight workers. Higher scaling is better.

| Model | Paths | Fast choice | Fast 8 ms | Fast 1→8 | Julia choice | Julia 8 ms | Julia 1→8 | Fast / Julia |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | 64 | cached_direct/ziggurat | 11.799 | 5.14× | Direct | 9.481 | 5.56× | 1.24× |
| W2 | 64 | direct_cr/ziggurat | 19.584 | 5.42× | DirectCR | 22.115 | 4.49× | 0.89× |
| W3 | 64 | cached_direct/ziggurat | 50.042 | 5.12× | SortingDirect | 51.929 | 5.78× | 0.96× |
| W4 | 64 | rssa_cr/ziggurat | 23.506 | 5.86× | RSSACR | 23.538 | 4.41× | 1.00× |
| W5 | 16 | direct_cr/ziggurat | 11.027 | 5.07× | RSSACR | 26.200 | 2.19× | 0.42× |
| W6 | 16 | rssa_cr/ziggurat | 9.428 | 6.13× | RSSACR | 26.748 | 1.51× | 0.35× |
| W7 | 16 | direct_cr/ziggurat | 3.902 | 5.48× | DirectCR | 8.511 | 2.06× | 0.46× |

## 256-path throughput batches

Times are batch milliseconds, median of five warm observations. The 1→8 column compares the **same selected algorithm**, batch and seed at one and eight workers. Higher scaling is better.

| Model | Paths | Fast choice | Fast 8 ms | Fast 1→8 | Julia choice | Julia 8 ms | Julia 1→8 | Fast / Julia |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | 256 | cached_direct/ziggurat | 44.973 | 4.98× | Direct | 34.044 | 6.32× | 1.32× |
| W2 | 256 | direct_cr/ziggurat | 79.099 | 5.79× | DirectCR | 79.466 | 4.69× | 1.00× |
| W3 | 256 | cached_direct/ziggurat | 187.343 | 5.70× | SortingDirect | 202.774 | 5.86× | 0.92× |
| W4 | 256 | rssa_cr/ziggurat | 92.985 | 5.28× | RSSACR | 86.834 | 4.88× | 1.07× |
| W5 | 256 | direct_cr/ziggurat | 161.852 | 5.67× | DirectCR | 174.464 | 4.80× | 0.93× |
| W6 | 256 | rssa_cr/ziggurat | 143.023 | 5.37× | RSSACR | 141.267 | 4.74× | 1.01× |
| W7 | 256 | rssa_cr/ziggurat | 46.370 | 6.36× | RSSACR | 41.810 | 5.27× | 1.11× |

## What eight threads means here

| Setting | Value |
| --- | --- |
| CPU | AMD EPYC 9V74 80-Core Processor |
| Platform | Linux-6.18.57-x86_64-with-glibc2.39 |
| CPU affinity | 0, 1, 2, 3, 4, 5, 6, 7, 8 |
| Cgroup CPU quota | 800000/100000 µs = 8 CPU equivalents |
| CRN-craft | 0.2.0.dev1; Rust 1.88.0 release, maturin wheel |
| Python | 3.12.14 (main, Aug 25 2026, 14:00:49) [Clang 22.1.3 ] |
| NumPy | 2.3.5 |
| Julia / Catalyst / JumpProcesses / SciMLBase | 1.13.3 / 17.0.0 / 10.0.1 / 3.58.0 |
| BLAS / OMP threads | 1 / 1 |
| Julia numerical threads, serial / parallel | 1 / 8 |
| Julia GC threads, serial / parallel | 1 / 8 |

Both implementations parallelize **independent trajectories**, not the reaction loop of one trajectory. CRN-craft creates eight scoped native workers per prepared batch and combines results in trajectory order. Julia uses `EnsembleThreads()` with eight numerical threads; the serial baseline uses `EnsembleSerial()`. Mutable JumpProblem state belongs to a task, with `safetycopy=false` and the pinned SciMLBase task-copy implementation. Each Julia trajectory receives an indexed Xoshiro RNG. BLAS is not the source of SSA parallelism. A matching runtime probe also finds one Julia interactive thread; it coordinates requests and is not a ninth numerical worker. The quota is eight CPU equivalents; dedicated physical cores and an otherwise idle host are not guaranteed.

The concurrency setup follows the [JumpProcesses ensemble guidance](https://docs.sciml.ai/JumpProcesses/stable/tutorials/ensembles_and_problem_reuse/) and [SciMLBase ensemble interface](https://docs.sciml.ai/SciMLBase/stable/interfaces/Ensembles/), with behavior checked against the pinned installed source.

Julia output callbacks observed thread IDs [2, 3, 4, 5, 6, 7, 8, 9]. Across the eight-worker observations, individual batches used 8–8 distinct output threads. These are measured callback thread IDs, not CPU indices; no reaction-loop instrumentation was added.

## Workloads and timing boundaries

| Model | Network | Species | Reactions | Horizon | Original paths |
| --- | --- | --- | --- | --- | --- |
| W1 | Birth/death | 1 | 2 | 100.0 | 64 |
| W2 | 33 independent birth/death pairs | 33 | 66 | 100.0 | 64 |
| W3 | Repressilator | 6 | 12 | 100.0 | 64 |
| W4 | Dense shared-species Michaelis–Menten | 8 | 66 | 100.0 | 64 |
| W5 | 256-species diffusion chain | 256 | 510 | 20.0 | 16 |
| W6 | Dense 32-species conversion | 32 | 992 | 20.0 | 16 |
| W7 | 32-compartment binding/diffusion | 96 | 126 | 20.0 | 16 |

Models, initial populations, rate conventions, horizons and save times `[0, end]` are identical across implementations. W4 uses the existing Catalyst translation valid on its equal-population invariant; explicit propensity and stoichiometry probes check that mapping. Timed output keeps means/variances and no retained trajectory history. Different RNGs and algorithmic random-number consumption mean event histories differ between implementations; this is equal simulation specification, not equal event traces.

CRN timing wraps the **prepared Python call**, including native worker creation, simulation, checks, statistics, result serialization and Python decoding. Julia timing wraps the **prepared ensemble solve**, including task problem copies, scheduling, trajectory completion/save/nonnegativity checks and ordered terminal statistics. It excludes IPC, final JSON serialization, replay digest calculation and the initial pristine batch deepcopy. Julia therefore has a narrower timing boundary; the table below also reports its measured client latency and pristine-copy cost for the selected methods. Both exclude model construction, installation, compilation and JIT warmup.

CRN's public result additionally carries per-trajectory event/work records, standard errors, and statistics at both saved times. The Julia driver reduces terminal counts to means/variances and checks the two saved times; it does not construct the full CRN result contract. Read the primary comparison as prepared API latency versus prepared Julia solve latency, rather than an isolated kernel speed measurement.

| Batch | Julia choice | Solve ms | Client ms | Pristine copy ms | Fast / Julia client |
| --- | --- | --- | --- | --- | --- |
| W1/original | Direct | 9.481 | 10.743 | 0.643 | 1.10× |
| W1/throughput | Direct | 34.044 | 35.761 | 0.645 | 1.26× |
| W2/original | DirectCR | 22.115 | 25.500 | 2.301 | 0.77× |
| W2/throughput | DirectCR | 79.466 | 81.994 | 2.175 | 0.96× |
| W3/original | SortingDirect | 51.929 | 53.668 | 1.259 | 0.93× |
| W3/throughput | SortingDirect | 202.774 | 204.381 | 1.038 | 0.92× |
| W4/original | RSSACR | 23.538 | 25.379 | 1.552 | 0.93× |
| W4/throughput | RSSACR | 86.834 | 89.246 | 1.852 | 1.04× |
| W5/original | RSSACR | 26.200 | 37.908 | 9.967 | 0.29× |
| W5/throughput | DirectCR | 174.464 | 190.148 | 10.500 | 0.85× |
| W6/original | RSSACR | 26.748 | 35.431 | 8.053 | 0.27× |
| W6/throughput | RSSACR | 141.267 | 149.295 | 9.044 | 0.96× |
| W7/original | DirectCR | 8.511 | 12.927 | 3.581 | 0.30× |
| W7/throughput | RSSACR | 41.810 | 47.382 | 3.576 | 0.98× |

Every configuration was warmed before timing. The two Julia processes received one full garbage collection after all warmups; no forced collections occurred between observations. Five rounds rotated method order and reversed it on alternate rounds. Numerical work ran sequentially between configurations; no build, package installation or validation job overlapped measured rounds. There are 2,240 timed observations: seven models × two batch sizes × two worker counts × sixteen configurations × five rounds. Timing seed is 4242.

| Julia workers | Process to all problems ready (s) | Package load (s) | Model/problem setup (s) |
| --- | --- | --- | --- |
| 1 | 67.776 | 8.263 | 50.181 |
| 8 | 63.524 | 7.026 | 48.703 |

These setup observations come from an already installed and precompiled Julia environment, once per worker configuration. They are excluded from warm times and do not constitute a cold-start or complete optimization-run comparison.

## Numerical checks and replay

All **6,514 predeclared terminal-marginal checks passed**, using 512 paths per model/engine for ten fast variants, five Julia algorithms and a legacy cached reference. Propensity/stoichiometry probes passed (8,870 rate values). Population checks and model-specific conservation checks passed. The baseline release suite passed 192 tests; formatting, warnings-denied Clippy, CLI/schema and alternative contract checks, and the dependency gate passed. Python contract/dependency checks also passed under `-O` with matching reports.

Family alpha is 0.01. W1/W2 compare against exact Poisson CDFs with bound 0.117261; other models compare against the independent reference using the conservative two-sample triangle bound 0.240224. These are finite-sample marginal checks with broad simultaneous bounds, not proof of joint trajectory-law equivalence, rare-event accuracy or optimization quality. No statistical thresholds were changed after observing samples.

| Model | CRN 1/8 replay | Julia 1/8 differing methods | Julia warm replay differences | Fresh SortingDirect |
| --- | --- | --- | --- | --- |
| W1 | 10/10 equal | none | none | equal / stable |
| W2 | 10/10 equal | SortingDirect | none | equal / stable |
| W3 | 10/10 equal | SortingDirect | none | equal / stable |
| W4 | 10/10 equal | SortingDirect | none | equal / stable |
| W5 | 10/10 equal | SortingDirect | none | equal / stable |
| W6 | 10/10 equal | SortingDirect | none | equal / stable |
| W7 | 10/10 equal | SortingDirect | none | equal / stable |

CRN one/eight validation compared full numerical results after removing only the run-spec hash; all timed CRN configurations replayed their warmup results. Julia one/eight replay checks compare the ordered terminal-count digest for 32 indexed paths per method/model. Standard SortingDirect may retain its learned reaction order within task-owned problems, so task partitioning or thread count can change seeded paths while preserving the intended distribution. `SortingDirect-fresh` sets `safetycopy=true`, creates fresh per-trajectory state and passed both cross-thread and warm replay. Its measurements are a diagnostic, excluded from best-method selection.

Fast mode was explicitly enabled with `method.fast=true` / `fast=True`. It uses per-replicate Xoshiro256++, running totals, and batched stop polling, while omitting the default lost-increment and 53-bit resolution checks. It retains checked state updates and budgets. These results do not extend validation to arbitrary population/rate scales. See the [fast-mode documentation](https://github.com/dietmarwo/CRN-craft/blob/c16fcc55654d0867a58af48d48a8e860efce8896/docs/alternative-ssa.md).

## Interpretation limits

Use the same-method scaling columns to judge whether a workload benefits from eight workers. Short batches can be dominated by worker creation, problem copies, reduction and serialization. Larger batches amortize some of that cost. The measurements alone do not isolate compiler code generation, sampling, cache behavior or scheduling as the cause of an implementation gap. Five observations on one shared host do not establish statistical significance for small differences. This report measures prepared SSA ensembles; it does not measure model recompilation for new parameters, full objective evaluation, optimizer convergence, final solution quality, or rebop.

## All measured configurations

Milliseconds below are **median [minimum–maximum]** across all five observations. A speedup below 1 means eight workers were slower than one. No observations were trimmed. `fresh` marks the excluded replay diagnostic.

### W1/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 80.954 [74.912–94.353] | 15.935 [15.481–23.950] | 5.08× |
| direct/ziggurat | 63.484 [59.435–77.177] | 12.568 [12.368–23.017] | 5.05× |
| cached_direct/inverse | 81.267 [73.612–93.709] | 15.116 [14.660–18.048] | 5.38× |
| cached_direct/ziggurat | 60.695 [55.351–68.692] | 11.799 [11.282–15.845] | 5.14× |
| direct_cr/inverse | 101.320 [94.607–106.293] | 19.759 [19.142–27.631] | 5.13× |
| direct_cr/ziggurat | 87.349 [74.736–98.297] | 16.278 [15.705–18.323] | 5.37× |
| sorting_direct/inverse | 75.618 [75.003–101.359] | 16.715 [16.021–19.850] | 4.52× |
| sorting_direct/ziggurat | 72.643 [64.856–75.517] | 13.088 [13.001–13.525] | 5.55× |
| rssa_cr/inverse | 101.751 [93.114–111.089] | 20.387 [20.201–20.651] | 4.99× |
| rssa_cr/ziggurat | 94.638 [79.089–108.786] | 17.179 [16.850–18.581] | 5.51× |
| Julia/Direct | 52.742 [50.477–58.265] | 9.481 [9.201–10.079] | 5.56× |
| Julia/SortingDirect | 58.351 [56.086–69.194] | 11.037 [10.794–11.173] | 5.29× |
| Julia/DirectCR | 86.247 [82.329–103.597] | 16.851 [16.462–17.462] | 5.12× |
| Julia/NRM | 72.144 [69.735–108.252] | 14.692 [14.568–14.908] | 4.91× |
| Julia/RSSACR | 74.628 [73.433–90.386] | 15.660 [14.624–17.838] | 4.77× |
| Julia/SortingDirect-fresh | 80.201 [78.782–86.060] | 20.229 [13.974–22.554] | 3.96× |

### W1/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 298.025 [286.303–330.424] | 67.952 [61.509–71.396] | 4.39× |
| direct/ziggurat | 238.679 [221.187–264.554] | 50.429 [49.545–77.067] | 4.73× |
| cached_direct/inverse | 292.005 [261.971–344.110] | 57.513 [57.007–70.551] | 5.08× |
| cached_direct/ziggurat | 223.838 [214.892–267.852] | 44.973 [44.827–51.277] | 4.98× |
| direct_cr/inverse | 397.518 [386.730–429.807] | 76.287 [75.353–80.478] | 5.21× |
| direct_cr/ziggurat | 332.934 [318.577–362.558] | 67.897 [62.568–79.558] | 4.90× |
| sorting_direct/inverse | 336.422 [303.272–395.037] | 65.222 [63.168–78.238] | 5.16× |
| sorting_direct/ziggurat | 268.671 [257.594–289.461] | 51.100 [50.348–62.867] | 5.26× |
| rssa_cr/inverse | 420.210 [390.921–439.156] | 79.575 [77.579–82.567] | 5.28× |
| rssa_cr/ziggurat | 400.485 [331.324–437.591] | 65.767 [64.306–66.666] | 6.09× |
| Julia/Direct | 215.211 [208.528–246.471] | 34.044 [34.004–37.760] | 6.32× |
| Julia/SortingDirect | 262.414 [234.518–272.975] | 41.163 [40.758–42.583] | 6.37× |
| Julia/DirectCR | 368.628 [338.852–421.025] | 64.062 [62.037–69.221] | 5.75× |
| Julia/NRM | 304.348 [280.861–342.737] | 55.158 [54.088–57.357] | 5.52× |
| Julia/RSSACR | 311.750 [298.391–326.129] | 55.828 [55.135–87.053] | 5.58× |
| Julia/SortingDirect-fresh | 344.062 [305.034–433.895] | 57.470 [54.547–61.608] | 5.99× |

### W2/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 572.165 [523.126–612.134] | 111.587 [109.712–115.599] | 5.13× |
| direct/ziggurat | 590.295 [516.710–604.992] | 111.165 [110.052–116.201] | 5.31× |
| cached_direct/inverse | 130.101 [127.839–142.111] | 23.048 [22.156–54.547] | 5.64× |
| cached_direct/ziggurat | 120.474 [115.122–123.578] | 20.266 [19.655–22.289] | 5.94× |
| direct_cr/inverse | 120.838 [110.662–131.004] | 22.577 [22.176–25.566] | 5.35× |
| direct_cr/ziggurat | 106.051 [97.625–137.387] | 19.584 [18.937–20.622] | 5.42× |
| sorting_direct/inverse | 165.169 [140.257–193.711] | 28.614 [28.315–40.088] | 5.77× |
| sorting_direct/ziggurat | 138.583 [132.705–150.330] | 25.601 [25.289–28.649] | 5.41× |
| rssa_cr/inverse | 139.690 [123.040–140.697] | 25.323 [24.511–33.459] | 5.52× |
| rssa_cr/ziggurat | 128.365 [122.270–133.848] | 22.760 [21.275–40.978] | 5.64× |
| Julia/Direct | 205.645 [190.967–256.647] | 50.342 [41.165–66.287] | 4.08× |
| Julia/SortingDirect | 129.163 [122.701–141.296] | 23.966 [20.192–45.983] | 5.39× |
| Julia/DirectCR | 99.262 [92.566–122.535] | 22.115 [20.050–43.076] | 4.49× |
| Julia/NRM | 114.280 [108.142–136.043] | 24.671 [21.987–36.127] | 4.63× |
| Julia/RSSACR | 107.166 [103.257–109.688] | 22.231 [20.769–25.328] | 4.82× |
| Julia/SortingDirect-fresh | 229.476 [203.791–247.889] | 46.731 [40.418–49.338] | 4.91× |

### W2/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 2320.830 [2180.435–2380.459] | 434.318 [424.954–437.646] | 5.34× |
| direct/ziggurat | 2436.165 [2158.145–2499.327] | 438.808 [403.194–454.907] | 5.55× |
| cached_direct/inverse | 519.580 [476.185–553.623] | 88.120 [86.794–142.378] | 5.90× |
| cached_direct/ziggurat | 500.682 [465.967–518.147] | 82.042 [75.982–128.663] | 6.10× |
| direct_cr/inverse | 458.670 [434.271–493.370] | 94.898 [86.208–154.906] | 4.83× |
| direct_cr/ziggurat | 458.207 [426.062–527.904] | 79.099 [75.463–96.195] | 5.79× |
| sorting_direct/inverse | 669.654 [605.892–826.487] | 115.100 [110.997–278.435] | 5.82× |
| sorting_direct/ziggurat | 587.029 [558.809–624.147] | 105.664 [98.310–116.815] | 5.56× |
| rssa_cr/inverse | 551.550 [521.507–571.079] | 141.512 [95.834–154.844] | 3.90× |
| rssa_cr/ziggurat | 501.758 [484.040–547.626] | 91.237 [83.294–104.877] | 5.50× |
| Julia/Direct | 838.674 [783.512–900.645] | 154.630 [152.778–180.921] | 5.42× |
| Julia/SortingDirect | 512.266 [490.753–585.900] | 83.752 [71.995–92.328] | 6.12× |
| Julia/DirectCR | 372.578 [369.521–439.800] | 79.466 [72.230–94.287] | 4.69× |
| Julia/NRM | 453.420 [437.986–520.169] | 87.215 [81.695–125.280] | 5.20× |
| Julia/RSSACR | 438.853 [411.715–486.784] | 80.515 [74.374–90.538] | 5.45× |
| Julia/SortingDirect-fresh | 927.311 [803.925–1047.520] | 173.269 [155.051–338.024] | 5.35× |

### W3/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 472.270 [457.364–511.499] | 101.999 [96.921–106.373] | 4.63× |
| direct/ziggurat | 419.119 [403.721–477.100] | 85.954 [84.409–90.764] | 4.88× |
| cached_direct/inverse | 290.885 [279.293–305.084] | 60.146 [59.207–63.383] | 4.84× |
| cached_direct/ziggurat | 256.420 [245.701–285.889] | 50.042 [43.735–52.354] | 5.12× |
| direct_cr/inverse | 427.373 [401.752–432.322] | 79.163 [76.899–91.561] | 5.40× |
| direct_cr/ziggurat | 412.854 [389.483–428.905] | 71.757 [70.659–74.255] | 5.75× |
| sorting_direct/inverse | 336.652 [321.559–357.900] | 64.919 [60.942–66.481] | 5.19× |
| sorting_direct/ziggurat | 292.186 [259.686–321.060] | 53.175 [50.194–55.657] | 5.49× |
| rssa_cr/inverse | 474.321 [445.618–523.088] | 85.057 [83.522–86.852] | 5.58× |
| rssa_cr/ziggurat | 479.083 [441.011–507.199] | 76.166 [73.912–79.718] | 6.29× |
| Julia/Direct | 443.516 [407.789–469.969] | 69.508 [68.415–73.123] | 6.38× |
| Julia/SortingDirect | 300.133 [270.360–314.997] | 51.929 [50.441–54.402] | 5.78× |
| Julia/DirectCR | 395.661 [375.732–428.115] | 71.970 [71.323–74.649] | 5.50× |
| Julia/NRM | 381.132 [358.005–441.471] | 76.915 [72.023–79.952] | 4.96× |
| Julia/RSSACR | 379.205 [351.623–411.250] | 61.329 [59.835–62.973] | 6.18× |
| Julia/SortingDirect-fresh | 338.162 [315.548–356.086] | 64.637 [56.038–96.438] | 5.23× |

### W3/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 1842.459 [1759.199–1879.961] | 392.099 [388.085–528.601] | 4.70× |
| direct/ziggurat | 1741.494 [1686.435–1791.829] | 337.265 [329.128–531.352] | 5.16× |
| cached_direct/inverse | 1249.766 [1173.299–1330.890] | 236.564 [233.404–335.546] | 5.28× |
| cached_direct/ziggurat | 1067.678 [974.582–1124.629] | 187.343 [185.956–277.186] | 5.70× |
| direct_cr/inverse | 1727.498 [1681.933–1784.405] | 326.516 [315.127–469.107] | 5.29× |
| direct_cr/ziggurat | 1652.911 [1621.978–1698.224] | 276.381 [270.891–384.750] | 5.98× |
| sorting_direct/inverse | 1276.874 [1216.977–1317.687] | 245.149 [244.572–346.535] | 5.21× |
| sorting_direct/ziggurat | 1172.548 [1083.649–1237.483] | 197.438 [194.547–265.565] | 5.94× |
| rssa_cr/inverse | 1851.137 [1775.213–1883.202] | 339.725 [337.309–455.825] | 5.45× |
| rssa_cr/ziggurat | 1742.077 [1642.443–1829.561] | 301.170 [287.773–375.351] | 5.78× |
| Julia/Direct | 1611.122 [1593.749–1838.227] | 259.779 [253.190–402.829] | 6.20× |
| Julia/SortingDirect | 1187.266 [1120.678–1358.255] | 202.774 [200.047–309.377] | 5.86× |
| Julia/DirectCR | 1705.417 [1591.585–1932.849] | 295.789 [277.339–400.562] | 5.77× |
| Julia/NRM | 1571.560 [1491.267–1799.565] | 279.431 [252.743–421.117] | 5.62× |
| Julia/RSSACR | 1498.597 [1446.717–1875.385] | 237.277 [225.971–348.195] | 6.32× |
| Julia/SortingDirect-fresh | 1380.642 [1286.240–1572.524] | 240.013 [235.046–373.029] | 5.75× |

### W4/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 812.878 [797.700–882.319] | 173.169 [146.493–238.427] | 4.69× |
| direct/ziggurat | 858.210 [751.395–889.787] | 170.916 [104.245–182.933] | 5.02× |
| cached_direct/inverse | 795.386 [717.028–821.848] | 146.720 [123.478–203.348] | 5.42× |
| cached_direct/ziggurat | 810.445 [716.843–884.282] | 141.278 [139.612–165.646] | 5.74× |
| direct_cr/inverse | 1074.184 [953.655–1194.969] | 202.392 [201.930–208.014] | 5.31× |
| direct_cr/ziggurat | 970.552 [948.843–1183.686] | 207.066 [203.813–222.372] | 4.69× |
| sorting_direct/inverse | 865.321 [823.182–999.286] | 169.527 [167.192–174.155] | 5.10× |
| sorting_direct/ziggurat | 863.109 [789.953–920.647] | 169.655 [144.761–171.021] | 5.09× |
| rssa_cr/inverse | 140.960 [128.243–175.097] | 27.541 [25.998–27.656] | 5.12× |
| rssa_cr/ziggurat | 137.658 [119.995–158.589] | 23.506 [23.143–23.979] | 5.86× |
| Julia/Direct | 247.168 [241.236–288.927] | 45.396 [42.672–47.401] | 5.44× |
| Julia/SortingDirect | 566.217 [552.531–724.941] | 121.305 [113.553–125.910] | 4.67× |
| Julia/DirectCR | 801.354 [723.940–951.448] | 168.183 [149.605–168.999] | 4.76× |
| Julia/NRM | 1037.474 [949.384–1083.842] | 196.516 [176.464–199.115] | 5.28× |
| Julia/RSSACR | 103.908 [100.921–124.085] | 23.538 [22.111–29.721] | 4.41× |
| Julia/SortingDirect-fresh | 662.319 [634.707–683.184] | 134.234 [127.129–136.447] | 4.93× |

### W4/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 3126.339 [3015.438–3347.120] | 678.070 [587.788–794.662] | 4.61× |
| direct/ziggurat | 3191.541 [3052.029–3348.195] | 690.657 [665.669–756.440] | 4.62× |
| cached_direct/inverse | 2977.135 [2955.794–3083.202] | 594.954 [557.450–623.817] | 5.00× |
| cached_direct/ziggurat | 3121.211 [2924.573–3281.113] | 559.079 [548.105–576.783] | 5.58× |
| direct_cr/inverse | 4066.863 [3937.939–4351.458] | 802.278 [740.044–859.246] | 5.07× |
| direct_cr/ziggurat | 4195.910 [4044.492–4289.080] | 798.160 [772.981–1127.049] | 5.26× |
| sorting_direct/inverse | 3393.164 [3284.807–3669.903] | 678.617 [624.548–763.259] | 5.00× |
| sorting_direct/ziggurat | 3396.561 [3285.756–3544.152] | 648.994 [632.675–921.376] | 5.23× |
| rssa_cr/inverse | 553.429 [530.155–592.078] | 106.506 [103.135–138.682] | 5.20× |
| rssa_cr/ziggurat | 490.796 [480.887–549.170] | 92.985 [90.088–145.587] | 5.28× |
| Julia/Direct | 1046.119 [974.380–1108.285] | 166.007 [163.185–234.393] | 6.30× |
| Julia/SortingDirect | 2335.999 [2143.770–2466.221] | 449.218 [446.947–457.666] | 5.20× |
| Julia/DirectCR | 3062.097 [3035.156–3475.544] | 663.218 [644.750–681.355] | 4.62× |
| Julia/NRM | 4015.300 [3854.586–4352.540] | 778.242 [774.178–827.711] | 5.16× |
| Julia/RSSACR | 423.411 [394.148–474.138] | 86.834 [83.932–114.270] | 4.88× |
| Julia/SortingDirect-fresh | 2703.768 [2645.240–3053.327] | 580.121 [491.468–703.275] | 4.66× |

### W5/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 1051.411 [1000.010–1097.007] | 219.799 [217.625–300.741] | 4.78× |
| direct/ziggurat | 1111.479 [1033.144–1212.513] | 228.820 [227.645–382.563] | 4.86× |
| cached_direct/inverse | 115.454 [107.813–128.873] | 16.978 [16.184–38.376] | 6.80× |
| cached_direct/ziggurat | 112.192 [105.466–124.260] | 17.691 [15.539–31.599] | 6.34× |
| direct_cr/inverse | 60.058 [57.004–72.017] | 11.830 [11.473–26.998] | 5.08× |
| direct_cr/ziggurat | 55.878 [54.709–71.778] | 11.027 [10.574–20.615] | 5.07× |
| sorting_direct/inverse | 170.129 [152.750–173.730] | 32.625 [31.863–53.617] | 5.21× |
| sorting_direct/ziggurat | 158.652 [152.148–164.855] | 32.137 [31.541–54.486] | 4.94× |
| rssa_cr/inverse | 71.492 [68.322–79.436] | 14.950 [13.773–22.311] | 4.78× |
| rssa_cr/ziggurat | 70.110 [64.212–72.232] | 14.326 [12.096–15.317] | 4.89× |
| Julia/Direct | 469.923 [440.595–508.824] | 140.542 [118.923–258.331] | 3.34× |
| Julia/SortingDirect | 127.806 [123.473–139.393] | 36.485 [33.110–45.197] | 3.50× |
| Julia/DirectCR | 49.409 [46.742–55.143] | 33.030 [26.306–53.661] | 1.50× |
| Julia/NRM | 69.215 [61.336–78.956] | 34.719 [28.417–38.105] | 1.99× |
| Julia/RSSACR | 57.335 [55.042–65.860] | 26.200 [25.065–50.266] | 2.19× |
| Julia/SortingDirect-fresh | 282.239 [265.175–354.841] | 61.662 [44.610–76.972] | 4.58× |

### W5/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 17148.952 [16484.494–18924.126] | 3335.624 [3306.018–3670.297] | 5.14× |
| direct/ziggurat | 18224.624 [17647.814–20280.132] | 3325.778 [3076.151–3615.691] | 5.48× |
| cached_direct/inverse | 1797.787 [1767.735–2028.849] | 254.835 [251.551–289.002] | 7.05× |
| cached_direct/ziggurat | 1811.226 [1720.067–2104.612] | 240.795 [234.199–330.670] | 7.52× |
| direct_cr/inverse | 965.080 [939.026–1128.173] | 177.015 [172.073–225.676] | 5.45× |
| direct_cr/ziggurat | 916.916 [869.939–1075.984] | 161.852 [160.176–204.097] | 5.67× |
| sorting_direct/inverse | 2713.174 [2490.487–2946.069] | 504.841 [500.666–531.835] | 5.37× |
| sorting_direct/ziggurat | 2586.358 [2449.080–3186.939] | 493.193 [460.602–542.000] | 5.24× |
| rssa_cr/inverse | 1211.469 [1169.686–1327.753] | 204.495 [203.182–206.794] | 5.92× |
| rssa_cr/ziggurat | 1104.397 [1014.532–1198.051] | 179.080 [178.604–183.785] | 6.17× |
| Julia/Direct | 7923.723 [7736.687–8878.019] | 1637.255 [1485.677–1668.184] | 4.84× |
| Julia/SortingDirect | 2182.426 [2087.607–2363.903] | 316.065 [282.214–354.406] | 6.90× |
| Julia/DirectCR | 837.709 [769.941–944.504] | 174.464 [158.855–195.796] | 4.80× |
| Julia/NRM | 1094.790 [1048.696–1206.486] | 202.204 [191.999–208.721] | 5.41× |
| Julia/RSSACR | 934.645 [917.701–1128.010] | 182.887 [178.984–208.638] | 5.11× |
| Julia/SortingDirect-fresh | 5166.795 [4830.355–5419.928] | 1030.794 [974.345–1073.831] | 5.01× |

### W6/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 1389.729 [1361.097–1529.029] | 299.093 [295.821–307.794] | 4.65× |
| direct/ziggurat | 1508.031 [1491.958–1532.563] | 315.587 [309.059–329.014] | 4.78× |
| cached_direct/inverse | 253.063 [230.001–273.744] | 40.939 [38.053–60.199] | 6.18× |
| cached_direct/ziggurat | 269.512 [234.324–299.174] | 39.260 [37.772–48.324] | 6.86× |
| direct_cr/inverse | 237.090 [223.370–313.615] | 48.103 [45.283–59.406] | 4.93× |
| direct_cr/ziggurat | 225.828 [214.504–244.112] | 45.029 [44.322–46.010] | 5.02× |
| sorting_direct/inverse | 306.652 [292.046–350.999] | 63.112 [61.463–67.005] | 4.86× |
| sorting_direct/ziggurat | 302.396 [294.576–325.262] | 62.245 [61.537–64.469] | 4.86× |
| rssa_cr/inverse | 51.035 [46.310–57.059] | 10.522 [9.723–16.662] | 4.85× |
| rssa_cr/ziggurat | 57.787 [46.288–70.613] | 9.428 [9.251–12.026] | 6.13× |
| Julia/Direct | 653.915 [644.351–790.811] | 155.618 [147.644–163.657] | 4.20× |
| Julia/SortingDirect | 227.661 [215.778–252.367] | 57.732 [42.335–177.664] | 3.94× |
| Julia/DirectCR | 217.342 [205.346–253.546] | 65.939 [53.219–70.338] | 3.30× |
| Julia/NRM | 314.898 [304.588–349.043] | 78.281 [68.848–160.226] | 4.02× |
| Julia/RSSACR | 40.381 [36.747–48.212] | 26.748 [18.763–39.494] | 1.51× |
| Julia/SortingDirect-fresh | 388.452 [381.743–599.450] | 78.991 [59.731–213.077] | 4.92× |

### W6/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 23066.909 [22747.589–24430.257] | 4450.573 [4129.645–4934.719] | 5.18× |
| direct/ziggurat | 24733.992 [23867.859–24976.558] | 4628.570 [4464.662–4711.015] | 5.34× |
| cached_direct/inverse | 4101.410 [3962.707–4285.939] | 603.424 [599.921–606.470] | 6.80× |
| cached_direct/ziggurat | 4170.929 [3896.213–4408.966] | 598.074 [587.916–617.313] | 6.97× |
| direct_cr/inverse | 3588.430 [3459.692–4052.781] | 722.429 [676.790–725.295] | 4.97× |
| direct_cr/ziggurat | 3457.772 [3443.194–4271.307] | 701.049 [666.472–707.462] | 4.93× |
| sorting_direct/inverse | 5064.538 [4924.830–5772.261] | 971.630 [919.457–1002.133] | 5.21× |
| sorting_direct/ziggurat | 5036.879 [4664.121–5458.743] | 978.483 [965.366–992.857] | 5.15× |
| rssa_cr/inverse | 833.940 [742.565–880.976] | 158.430 [147.932–164.519] | 5.26× |
| rssa_cr/ziggurat | 768.128 [747.057–922.326] | 143.023 [138.828–238.696] | 5.37× |
| Julia/Direct | 10502.531 [10256.248–11356.974] | 2184.773 [2000.755–2232.516] | 4.81× |
| Julia/SortingDirect | 3140.812 [2793.272–3487.370] | 503.036 [486.744–527.647] | 6.24× |
| Julia/DirectCR | 3414.713 [3219.504–4096.317] | 689.957 [632.817–702.828] | 4.95× |
| Julia/NRM | 5009.541 [4852.382–5708.608] | 929.535 [890.018–947.584] | 5.39× |
| Julia/RSSACR | 669.604 [622.152–887.396] | 141.267 [134.222–150.260] | 4.74× |
| Julia/SortingDirect-fresh | 6109.619 [5590.521–7738.590] | 1379.807 [1294.531–1490.332] | 4.43× |

### W7/original

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 95.632 [87.647–145.431] | 21.192 [20.288–21.725] | 4.51× |
| direct/ziggurat | 109.406 [107.508–120.743] | 21.948 [21.491–24.596] | 4.98× |
| cached_direct/inverse | 19.316 [19.138–28.348] | 4.285 [3.774–4.902] | 4.51× |
| cached_direct/ziggurat | 21.267 [20.454–23.346] | 4.115 [3.865–6.278] | 5.17× |
| direct_cr/inverse | 20.711 [17.719–30.119] | 4.308 [4.015–5.804] | 4.81× |
| direct_cr/ziggurat | 21.401 [19.828–25.307] | 3.902 [3.766–4.017] | 5.48× |
| sorting_direct/inverse | 23.096 [20.344–26.156] | 5.221 [4.586–5.523] | 4.42× |
| sorting_direct/ziggurat | 20.179 [18.544–25.568] | 5.568 [4.294–7.490] | 3.62× |
| rssa_cr/inverse | 19.753 [18.818–23.499] | 4.853 [3.807–43.302] | 4.07× |
| rssa_cr/ziggurat | 19.572 [18.260–22.470] | 5.185 [3.476–8.084] | 3.77× |
| Julia/Direct | 46.777 [40.204–48.099] | 15.582 [15.064–16.222] | 3.00× |
| Julia/SortingDirect | 18.697 [16.502–19.920] | 11.252 [8.572–22.865] | 1.66× |
| Julia/DirectCR | 17.517 [14.574–22.118] | 8.511 [8.074–14.465] | 2.06× |
| Julia/NRM | 19.446 [17.318–22.509] | 9.016 [8.248–11.077] | 2.16× |
| Julia/RSSACR | 14.018 [13.567–15.959] | 9.007 [7.513–17.389] | 1.56× |
| Julia/SortingDirect-fresh | 84.236 [62.668–130.362] | 19.974 [12.127–24.424] | 4.22× |

### W7/throughput

| Engine / waiting time | 1 worker ms | 8 workers ms | 1→8 |
| --- | --- | --- | --- |
| direct/inverse | 1502.941 [1399.340–1652.365] | 310.496 [305.444–325.811] | 4.84× |
| direct/ziggurat | 1643.409 [1538.161–1681.217] | 317.432 [313.481–318.978] | 5.18× |
| cached_direct/inverse | 322.512 [312.802–393.356] | 55.061 [53.173–55.961] | 5.86× |
| cached_direct/ziggurat | 322.623 [308.883–347.994] | 51.299 [49.936–60.843] | 6.29× |
| direct_cr/inverse | 320.264 [294.758–348.507] | 58.583 [57.143–59.107] | 5.47× |
| direct_cr/ziggurat | 310.143 [269.870–392.506] | 52.971 [52.705–55.625] | 5.85× |
| sorting_direct/inverse | 360.961 [320.789–392.935] | 67.013 [66.398–68.462] | 5.39× |
| sorting_direct/ziggurat | 332.116 [301.374–342.315] | 64.514 [62.931–92.456] | 5.15× |
| rssa_cr/inverse | 293.791 [270.086–310.602] | 52.416 [51.878–55.258] | 5.60× |
| rssa_cr/ziggurat | 294.754 [266.006–352.016] | 46.370 [46.086–52.984] | 6.36× |
| Julia/Direct | 697.242 [666.690–861.538] | 165.367 [159.013–169.905] | 4.22× |
| Julia/SortingDirect | 266.245 [257.731–286.625] | 44.608 [43.801–49.602] | 5.97× |
| Julia/DirectCR | 248.936 [230.685–255.264] | 54.135 [50.518–62.756] | 4.60× |
| Julia/NRM | 302.200 [292.814–310.942] | 59.799 [58.917–63.999] | 5.05× |
| Julia/RSSACR | 220.198 [209.275–227.299] | 41.810 [41.517–43.022] | 5.27× |
| Julia/SortingDirect-fresh | 1041.442 [978.436–1195.435] | 245.708 [236.191–381.871] | 4.24× |

## Individual observations

Ordered milliseconds for rounds 1–5. These values make every summary median and range independently checkable without an accompanying file.

### W1/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 80.954, 83.251, 74.912, 77.021, 94.353 | 15.935, 15.481, 16.690, 15.768, 23.950 |
| direct/ziggurat | 59.435, 63.484, 60.710, 72.876, 77.177 | 12.568, 13.117, 12.559, 12.368, 23.017 |
| cached_direct/inverse | 81.267, 93.709, 78.922, 73.612, 81.822 | 15.116, 15.203, 14.776, 14.660, 18.048 |
| cached_direct/ziggurat | 68.692, 62.492, 55.351, 59.273, 60.695 | 12.940, 11.282, 11.433, 11.799, 15.845 |
| direct_cr/inverse | 103.208, 101.320, 106.293, 99.472, 94.607 | 26.762, 19.142, 19.759, 19.689, 27.631 |
| direct_cr/ziggurat | 92.373, 74.736, 84.453, 98.297, 87.349 | 16.278, 16.720, 16.246, 15.705, 18.323 |
| sorting_direct/inverse | 75.239, 88.112, 75.003, 75.618, 101.359 | 16.859, 16.021, 16.715, 16.413, 19.850 |
| sorting_direct/ziggurat | 73.233, 75.517, 66.356, 64.856, 72.643 | 13.525, 13.088, 13.442, 13.065, 13.001 |
| rssa_cr/inverse | 93.114, 111.089, 109.656, 101.751, 93.383 | 20.455, 20.201, 20.254, 20.651, 20.387 |
| rssa_cr/ziggurat | 95.100, 108.786, 89.606, 79.089, 94.638 | 17.179, 18.581, 17.049, 17.333, 16.850 |
| Julia/Direct | 53.351, 50.725, 52.742, 50.477, 58.265 | 9.481, 9.271, 10.079, 9.201, 9.664 |
| Julia/SortingDirect | 58.351, 56.086, 62.398, 58.200, 69.194 | 10.794, 11.173, 11.037, 10.828, 11.147 |
| Julia/DirectCR | 86.247, 83.723, 89.460, 82.329, 103.597 | 16.951, 17.462, 16.770, 16.851, 16.462 |
| Julia/NRM | 71.502, 72.144, 75.918, 69.735, 108.252 | 14.568, 14.649, 14.716, 14.908, 14.692 |
| Julia/RSSACR | 85.722, 73.433, 74.628, 74.339, 90.386 | 14.624, 17.838, 14.661, 16.608, 15.660 |
| Julia/SortingDirect-fresh | 78.782, 80.201, 78.937, 81.691, 86.060 | 13.974, 22.409, 20.229, 22.554, 18.634 |

### W1/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 329.415, 298.025, 292.548, 286.303, 330.424 | 61.509, 71.396, 67.952, 64.371, 68.215 |
| direct/ziggurat | 241.255, 222.129, 238.679, 221.187, 264.554 | 49.545, 50.219, 59.972, 77.067, 50.429 |
| cached_direct/inverse | 344.110, 319.964, 281.189, 261.971, 292.005 | 57.007, 57.513, 57.453, 59.267, 70.551 |
| cached_direct/ziggurat | 227.723, 216.131, 267.852, 223.838, 214.892 | 51.277, 44.930, 44.827, 47.268, 44.973 |
| direct_cr/inverse | 429.807, 397.518, 386.730, 387.375, 412.028 | 80.478, 75.474, 75.353, 76.287, 79.435 |
| direct_cr/ziggurat | 327.214, 318.577, 332.934, 362.558, 350.123 | 67.897, 62.568, 64.772, 68.006, 79.558 |
| sorting_direct/inverse | 303.272, 336.422, 332.617, 340.938, 395.037 | 65.454, 63.266, 65.222, 63.168, 78.238 |
| sorting_direct/ziggurat | 285.377, 268.671, 257.594, 263.523, 289.461 | 51.100, 50.348, 59.247, 50.732, 62.867 |
| rssa_cr/inverse | 390.921, 420.210, 394.539, 424.345, 439.156 | 81.308, 77.579, 78.788, 79.575, 82.567 |
| rssa_cr/ziggurat | 412.920, 400.485, 331.324, 437.591, 395.817 | 66.666, 64.306, 66.043, 65.145, 65.767 |
| Julia/Direct | 208.528, 215.211, 210.410, 237.619, 246.471 | 34.031, 34.004, 34.044, 34.643, 37.760 |
| Julia/SortingDirect | 251.924, 234.518, 272.975, 272.018, 262.414 | 41.163, 40.775, 40.758, 41.786, 42.583 |
| Julia/DirectCR | 368.628, 342.998, 338.852, 421.025, 371.525 | 64.062, 62.037, 63.780, 68.000, 69.221 |
| Julia/NRM | 296.874, 280.861, 304.348, 342.737, 330.549 | 55.139, 55.158, 54.088, 57.357, 56.959 |
| Julia/RSSACR | 310.484, 298.391, 326.129, 318.882, 311.750 | 55.828, 55.479, 55.135, 87.053, 66.573 |
| Julia/SortingDirect-fresh | 344.062, 335.182, 350.031, 305.034, 433.895 | 54.547, 58.684, 54.998, 57.470, 61.608 |

### W2/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 612.134, 604.697, 564.000, 523.126, 572.165 | 109.712, 111.983, 111.587, 115.599, 110.321 |
| direct/ziggurat | 590.295, 592.542, 516.710, 604.992, 570.581 | 110.052, 111.165, 113.654, 116.201, 111.069 |
| cached_direct/inverse | 142.111, 127.839, 130.101, 130.733, 129.680 | 22.156, 22.878, 23.048, 54.547, 25.776 |
| cached_direct/ziggurat | 115.122, 121.539, 119.658, 123.578, 120.474 | 19.655, 19.918, 20.724, 22.289, 20.266 |
| direct_cr/inverse | 110.662, 112.573, 131.004, 120.838, 123.286 | 22.176, 22.229, 25.566, 23.266, 22.577 |
| direct_cr/ziggurat | 105.079, 97.625, 118.234, 106.051, 137.387 | 18.937, 19.584, 19.781, 19.171, 20.622 |
| sorting_direct/inverse | 165.169, 140.257, 168.220, 193.711, 162.069 | 28.315, 28.614, 35.466, 28.496, 40.088 |
| sorting_direct/ziggurat | 150.330, 138.583, 143.683, 132.705, 134.477 | 25.289, 28.649, 25.479, 25.601, 28.550 |
| rssa_cr/inverse | 139.926, 123.040, 140.697, 139.690, 136.228 | 24.511, 25.031, 25.323, 25.596, 33.459 |
| rssa_cr/ziggurat | 122.270, 125.648, 133.848, 128.365, 131.116 | 21.275, 40.978, 22.521, 22.760, 33.607 |
| Julia/Direct | 190.967, 197.457, 256.647, 205.645, 215.759 | 41.165, 66.287, 64.019, 50.342, 49.794 |
| Julia/SortingDirect | 122.701, 129.163, 123.328, 131.529, 141.296 | 23.966, 20.192, 42.615, 45.983, 21.371 |
| Julia/DirectCR | 99.262, 101.697, 92.566, 97.200, 122.535 | 20.581, 22.115, 23.590, 43.076, 20.050 |
| Julia/NRM | 108.142, 114.280, 111.531, 136.043, 119.196 | 23.498, 28.900, 36.127, 24.671, 21.987 |
| Julia/RSSACR | 107.166, 109.688, 104.993, 109.428, 103.257 | 20.769, 25.328, 24.530, 21.598, 22.231 |
| Julia/SortingDirect-fresh | 203.791, 247.889, 213.750, 229.476, 243.312 | 40.418, 43.539, 48.312, 46.731, 49.338 |

### W2/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 2180.435, 2270.905, 2320.830, 2336.408, 2380.459 | 435.766, 434.318, 437.646, 427.604, 424.954 |
| direct/ziggurat | 2436.165, 2442.797, 2499.327, 2253.095, 2158.145 | 438.808, 435.803, 443.525, 403.194, 454.907 |
| cached_direct/inverse | 553.623, 519.580, 520.843, 482.469, 476.185 | 88.120, 86.794, 142.378, 106.269, 87.630 |
| cached_direct/ziggurat | 506.034, 518.147, 500.682, 482.879, 465.967 | 82.042, 76.646, 128.663, 92.913, 75.982 |
| direct_cr/inverse | 458.670, 443.129, 469.970, 493.370, 434.271 | 98.331, 86.208, 154.906, 94.898, 90.444 |
| direct_cr/ziggurat | 426.062, 501.137, 458.207, 527.904, 444.166 | 77.806, 75.463, 96.195, 83.816, 79.099 |
| sorting_direct/inverse | 669.654, 608.074, 688.135, 826.487, 605.892 | 118.346, 110.997, 114.748, 115.100, 278.435 |
| sorting_direct/ziggurat | 624.147, 584.452, 587.029, 591.405, 558.809 | 116.815, 98.310, 98.367, 105.664, 106.597 |
| rssa_cr/inverse | 560.797, 521.507, 571.079, 551.550, 541.776 | 124.308, 95.834, 154.844, 141.512, 148.925 |
| rssa_cr/ziggurat | 495.834, 484.040, 541.181, 547.626, 501.758 | 104.877, 83.294, 90.270, 91.237, 92.650 |
| Julia/Direct | 795.845, 783.512, 900.645, 847.562, 838.674 | 180.921, 157.453, 152.778, 154.630, 153.543 |
| Julia/SortingDirect | 498.145, 490.753, 585.900, 512.266, 517.574 | 92.328, 74.392, 71.995, 86.804, 83.752 |
| Julia/DirectCR | 385.948, 372.578, 439.800, 369.521, 372.017 | 94.287, 72.230, 79.466, 73.049, 87.950 |
| Julia/NRM | 453.420, 437.986, 520.169, 480.355, 439.803 | 125.280, 87.215, 87.646, 81.695, 82.396 |
| Julia/RSSACR | 441.014, 411.715, 486.784, 438.853, 416.395 | 90.538, 74.374, 76.378, 85.391, 80.515 |
| Julia/SortingDirect-fresh | 933.174, 803.925, 1047.520, 869.297, 927.311 | 183.518, 164.702, 173.269, 155.051, 338.024 |

### W3/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 482.953, 468.070, 511.499, 457.364, 472.270 | 98.526, 106.373, 102.201, 96.921, 101.999 |
| direct/ziggurat | 438.240, 419.119, 419.037, 403.721, 477.100 | 84.409, 90.764, 85.954, 85.797, 86.263 |
| cached_direct/inverse | 305.084, 290.885, 298.394, 279.293, 290.419 | 59.207, 60.146, 59.926, 61.308, 63.383 |
| cached_direct/ziggurat | 253.141, 256.420, 285.889, 270.277, 245.701 | 43.735, 52.354, 49.214, 50.042, 50.307 |
| direct_cr/inverse | 401.752, 431.596, 407.720, 427.373, 432.322 | 78.393, 79.163, 80.911, 76.899, 91.561 |
| direct_cr/ziggurat | 389.483, 412.854, 416.412, 428.905, 404.653 | 74.255, 74.001, 70.659, 70.885, 71.757 |
| sorting_direct/inverse | 321.559, 325.509, 357.900, 336.652, 348.350 | 64.146, 64.919, 66.481, 65.623, 60.942 |
| sorting_direct/ziggurat | 295.577, 259.686, 321.060, 292.186, 266.240 | 53.193, 55.657, 51.778, 50.194, 53.175 |
| rssa_cr/inverse | 474.321, 445.618, 475.917, 469.166, 523.088 | 85.246, 86.852, 83.522, 85.057, 84.883 |
| rssa_cr/ziggurat | 507.199, 468.797, 494.536, 479.083, 441.011 | 79.718, 77.224, 73.912, 75.730, 76.166 |
| Julia/Direct | 427.019, 443.516, 407.789, 463.938, 469.969 | 73.123, 69.069, 70.834, 69.508, 68.415 |
| Julia/SortingDirect | 298.933, 300.133, 314.997, 270.360, 301.426 | 50.441, 51.929, 51.608, 52.668, 54.402 |
| Julia/DirectCR | 395.661, 428.115, 410.991, 375.732, 389.894 | 71.970, 71.323, 74.649, 72.769, 71.965 |
| Julia/NRM | 383.664, 441.471, 381.132, 367.221, 358.005 | 76.915, 79.952, 73.964, 72.023, 77.022 |
| Julia/RSSACR | 379.205, 411.250, 370.922, 351.623, 381.136 | 62.907, 61.088, 61.329, 62.973, 59.835 |
| Julia/SortingDirect-fresh | 335.214, 356.086, 338.162, 315.548, 343.550 | 60.235, 64.637, 96.438, 66.823, 56.038 |

### W3/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 1879.961, 1842.459, 1801.748, 1759.199, 1854.965 | 528.601, 399.766, 391.455, 388.085, 392.099 |
| direct/ziggurat | 1791.829, 1763.248, 1735.978, 1741.494, 1686.435 | 531.352, 337.265, 329.962, 329.128, 339.952 |
| cached_direct/inverse | 1318.244, 1249.766, 1198.997, 1173.299, 1330.890 | 335.546, 287.523, 236.564, 233.404, 234.359 |
| cached_direct/ziggurat | 1124.629, 1118.854, 974.582, 1064.670, 1067.678 | 277.186, 252.385, 186.227, 187.343, 185.956 |
| direct_cr/inverse | 1690.561, 1764.283, 1784.405, 1727.498, 1681.933 | 469.107, 385.928, 322.184, 326.516, 315.127 |
| direct_cr/ziggurat | 1675.009, 1649.329, 1652.911, 1621.978, 1698.224 | 384.750, 296.902, 276.381, 273.671, 270.891 |
| sorting_direct/inverse | 1317.687, 1258.763, 1276.874, 1216.977, 1289.176 | 328.130, 346.535, 244.661, 244.572, 245.149 |
| sorting_direct/ziggurat | 1229.747, 1083.649, 1172.548, 1095.165, 1237.483 | 264.058, 265.565, 197.438, 195.871, 194.547 |
| rssa_cr/inverse | 1822.453, 1883.202, 1879.859, 1851.137, 1775.213 | 432.753, 455.825, 337.309, 339.725, 337.592 |
| rssa_cr/ziggurat | 1829.561, 1742.077, 1642.443, 1764.126, 1662.603 | 375.351, 368.833, 287.773, 301.170, 287.977 |
| Julia/Direct | 1838.227, 1593.749, 1754.430, 1602.574, 1611.122 | 402.829, 253.190, 259.779, 313.700, 258.433 |
| Julia/SortingDirect | 1358.255, 1279.171, 1187.266, 1120.678, 1136.286 | 309.377, 226.125, 202.774, 200.047, 202.754 |
| Julia/DirectCR | 1932.849, 1705.417, 1720.363, 1629.758, 1591.585 | 327.000, 400.562, 295.789, 277.339, 287.013 |
| Julia/NRM | 1799.565, 1571.560, 1491.267, 1556.519, 1572.658 | 295.591, 421.117, 279.431, 252.743, 273.454 |
| Julia/RSSACR | 1875.385, 1535.817, 1498.597, 1476.284, 1446.717 | 247.746, 348.195, 237.277, 225.971, 231.902 |
| Julia/SortingDirect-fresh | 1572.524, 1371.952, 1384.064, 1286.240, 1380.642 | 247.648, 373.029, 235.046, 240.013, 235.788 |

### W4/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 860.449, 882.319, 797.700, 812.878, 811.956 | 177.734, 238.427, 173.169, 146.493, 170.926 |
| direct/ziggurat | 862.799, 889.787, 751.395, 858.210, 797.016 | 177.082, 182.933, 170.916, 168.414, 104.245 |
| cached_direct/inverse | 795.386, 821.848, 786.809, 798.062, 717.028 | 203.348, 149.935, 146.720, 145.044, 123.478 |
| cached_direct/ziggurat | 884.282, 871.490, 810.445, 800.314, 716.843 | 139.612, 165.646, 142.094, 139.710, 141.278 |
| direct_cr/inverse | 1108.697, 1194.969, 953.655, 1074.184, 960.960 | 201.930, 204.695, 208.014, 202.392, 202.167 |
| direct_cr/ziggurat | 1039.634, 1183.686, 970.552, 948.843, 954.403 | 207.066, 222.372, 206.036, 210.773, 203.813 |
| sorting_direct/inverse | 999.286, 938.718, 826.900, 865.321, 823.182 | 167.934, 169.527, 167.192, 170.421, 174.155 |
| sorting_direct/ziggurat | 920.647, 909.582, 863.109, 821.326, 789.953 | 169.333, 171.021, 169.655, 170.943, 144.761 |
| rssa_cr/inverse | 175.097, 140.960, 142.263, 128.243, 131.985 | 27.656, 27.566, 26.884, 27.541, 25.998 |
| rssa_cr/ziggurat | 158.589, 138.052, 137.658, 119.995, 122.727 | 23.979, 23.506, 23.143, 23.539, 23.416 |
| Julia/Direct | 288.927, 253.783, 247.168, 244.939, 241.236 | 45.700, 44.608, 47.401, 42.672, 45.396 |
| Julia/SortingDirect | 724.941, 561.291, 588.844, 566.217, 552.531 | 125.910, 123.901, 121.305, 113.960, 113.553 |
| Julia/DirectCR | 951.448, 809.815, 767.910, 801.354, 723.940 | 168.909, 168.183, 168.999, 166.498, 149.605 |
| Julia/NRM | 1083.842, 1050.093, 972.925, 1037.474, 949.384 | 199.115, 177.762, 196.516, 176.464, 196.805 |
| Julia/RSSACR | 117.885, 124.085, 100.921, 103.908, 101.880 | 29.721, 29.289, 22.111, 23.538, 22.637 |
| Julia/SortingDirect-fresh | 662.319, 683.184, 679.118, 634.707, 648.295 | 127.129, 127.481, 136.447, 134.234, 136.414 |

### W4/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 3347.120, 3037.279, 3126.339, 3015.438, 3128.129 | 768.544, 794.662, 587.788, 678.070, 658.123 |
| direct/ziggurat | 3348.195, 3107.707, 3248.738, 3191.541, 3052.029 | 756.440, 665.669, 682.770, 690.657, 727.306 |
| cached_direct/inverse | 2955.794, 3083.202, 2987.839, 2977.135, 2957.929 | 623.817, 557.450, 596.086, 587.942, 594.954 |
| cached_direct/ziggurat | 3281.113, 3234.772, 3121.211, 2924.573, 3002.283 | 576.783, 559.472, 559.079, 557.959, 548.105 |
| direct_cr/inverse | 4255.426, 4351.458, 4066.863, 3937.939, 4066.511 | 859.246, 743.342, 802.278, 808.588, 740.044 |
| direct_cr/ziggurat | 4195.910, 4244.088, 4289.080, 4044.492, 4072.561 | 1127.049, 772.981, 798.160, 832.171, 780.556 |
| sorting_direct/inverse | 3669.903, 3565.583, 3284.807, 3393.164, 3392.770 | 763.259, 624.548, 678.617, 681.275, 636.183 |
| sorting_direct/ziggurat | 3544.152, 3396.561, 3326.028, 3285.756, 3459.254 | 921.376, 632.675, 640.196, 657.044, 648.994 |
| rssa_cr/inverse | 592.078, 555.164, 530.155, 549.611, 553.429 | 138.682, 106.840, 103.135, 103.148, 106.506 |
| rssa_cr/ziggurat | 549.170, 526.877, 482.383, 480.887, 490.796 | 145.587, 91.276, 90.088, 94.333, 92.985 |
| Julia/Direct | 1108.285, 977.775, 1046.119, 974.380, 1104.679 | 234.393, 198.564, 163.185, 164.293, 166.007 |
| Julia/SortingDirect | 2466.221, 2266.400, 2335.999, 2143.770, 2337.716 | 457.666, 451.160, 446.947, 448.434, 449.218 |
| Julia/DirectCR | 3475.544, 3037.023, 3062.097, 3323.763, 3035.156 | 671.748, 653.935, 644.750, 681.355, 663.218 |
| Julia/NRM | 4352.540, 4289.458, 3951.394, 3854.586, 4015.300 | 793.681, 774.467, 827.711, 778.242, 774.178 |
| Julia/RSSACR | 474.138, 458.755, 403.683, 423.411, 394.148 | 85.523, 86.926, 114.270, 83.932, 86.834 |
| Julia/SortingDirect-fresh | 3053.327, 2857.223, 2663.403, 2645.240, 2703.768 | 580.121, 703.275, 583.761, 515.908, 491.468 |

### W5/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 1092.152, 1097.007, 1000.010, 1051.411, 1002.745 | 300.741, 220.211, 219.799, 218.954, 217.625 |
| direct/ziggurat | 1212.513, 1113.177, 1111.479, 1096.414, 1033.144 | 382.563, 228.820, 227.645, 229.973, 228.562 |
| cached_direct/inverse | 128.873, 115.454, 107.813, 116.952, 110.559 | 38.376, 16.978, 16.184, 16.609, 18.372 |
| cached_direct/ziggurat | 112.192, 124.260, 111.845, 120.872, 105.466 | 31.599, 18.461, 15.539, 16.349, 17.691 |
| direct_cr/inverse | 66.820, 57.004, 72.017, 60.058, 58.663 | 26.998, 11.862, 11.473, 11.676, 11.830 |
| direct_cr/ziggurat | 70.141, 54.912, 54.709, 71.778, 55.878 | 20.615, 14.538, 10.574, 10.913, 11.027 |
| sorting_direct/inverse | 170.129, 171.742, 157.649, 173.730, 152.750 | 53.617, 32.115, 31.863, 33.887, 32.625 |
| sorting_direct/ziggurat | 161.653, 164.855, 152.148, 158.652, 152.563 | 54.486, 31.590, 32.137, 32.479, 31.541 |
| rssa_cr/inverse | 72.162, 68.322, 79.436, 71.492, 68.680 | 22.311, 20.939, 14.360, 14.950, 13.773 |
| rssa_cr/ziggurat | 72.232, 72.015, 65.282, 70.110, 64.212 | 15.317, 15.181, 12.096, 14.326, 12.170 |
| Julia/Direct | 508.824, 440.595, 475.868, 469.923, 457.628 | 206.123, 140.542, 121.751, 258.331, 118.923 |
| Julia/SortingDirect | 139.393, 123.473, 127.806, 130.635, 127.090 | 40.050, 36.485, 33.110, 33.125, 45.197 |
| Julia/DirectCR | 55.143, 46.742, 49.085, 49.755, 49.409 | 53.661, 33.030, 37.916, 30.910, 26.306 |
| Julia/NRM | 71.284, 62.230, 61.336, 69.215, 78.956 | 36.499, 33.685, 28.417, 34.719, 38.105 |
| Julia/RSSACR | 63.371, 55.042, 56.608, 65.860, 57.335 | 50.266, 25.065, 37.990, 26.200, 26.032 |
| Julia/SortingDirect-fresh | 354.841, 282.239, 267.607, 321.156, 265.175 | 76.972, 54.467, 61.662, 44.610, 62.969 |

### W5/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 18924.126, 16484.494, 17148.952, 17299.481, 16863.756 | 3670.297, 3335.624, 3474.705, 3306.018, 3324.641 |
| direct/ziggurat | 20280.132, 17780.667, 17647.814, 18541.195, 18224.624 | 3615.691, 3510.829, 3325.778, 3311.795, 3076.151 |
| cached_direct/inverse | 2028.849, 1767.735, 1797.787, 1782.159, 1879.937 | 289.002, 257.029, 251.700, 254.835, 251.551 |
| cached_direct/ziggurat | 2104.612, 1744.337, 1811.226, 1817.207, 1720.067 | 330.670, 263.367, 240.795, 238.167, 234.199 |
| direct_cr/inverse | 1128.173, 967.675, 965.080, 952.121, 939.026 | 225.676, 178.462, 172.073, 173.040, 177.015 |
| direct_cr/ziggurat | 1075.984, 885.467, 963.331, 916.916, 869.939 | 204.097, 164.229, 161.852, 161.638, 160.176 |
| sorting_direct/inverse | 2946.069, 2490.487, 2731.070, 2713.174, 2516.185 | 531.835, 513.533, 504.841, 501.184, 500.666 |
| sorting_direct/ziggurat | 3186.939, 2449.080, 2601.046, 2586.358, 2577.897 | 542.000, 506.158, 493.193, 460.602, 487.954 |
| rssa_cr/inverse | 1327.753, 1228.588, 1194.885, 1169.686, 1211.469 | 205.335, 204.199, 204.495, 203.182, 206.794 |
| rssa_cr/ziggurat | 1198.051, 1169.128, 1014.532, 1104.397, 1060.407 | 180.460, 178.604, 179.080, 178.673, 183.785 |
| Julia/Direct | 8878.019, 8114.480, 7923.723, 7736.687, 7854.520 | 1546.674, 1637.255, 1485.677, 1668.184, 1653.600 |
| Julia/SortingDirect | 2363.903, 2115.525, 2182.426, 2305.897, 2087.607 | 354.406, 282.214, 321.569, 294.654, 316.065 |
| Julia/DirectCR | 944.504, 837.709, 796.891, 928.585, 769.941 | 195.796, 158.855, 190.169, 163.009, 174.464 |
| Julia/NRM | 1206.486, 1048.696, 1094.790, 1180.285, 1058.193 | 195.178, 202.204, 191.999, 208.721, 203.117 |
| Julia/RSSACR | 1064.052, 928.396, 934.645, 1128.010, 917.701 | 208.638, 178.984, 188.823, 179.031, 182.887 |
| Julia/SortingDirect-fresh | 5419.928, 5166.795, 4830.355, 5267.615, 4859.569 | 1009.524, 1067.981, 1030.794, 974.345, 1073.831 |

### W6/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 1369.111, 1529.029, 1361.097, 1389.729, 1430.517 | 295.821, 303.626, 307.794, 296.333, 299.093 |
| direct/ziggurat | 1508.031, 1491.958, 1532.563, 1510.545, 1507.257 | 309.499, 329.014, 317.753, 309.059, 315.587 |
| cached_direct/inverse | 273.744, 230.001, 253.063, 253.332, 249.816 | 40.939, 60.199, 39.694, 38.053, 42.013 |
| cached_direct/ziggurat | 269.512, 234.324, 299.174, 284.027, 263.122 | 39.260, 42.941, 48.324, 37.772, 37.999 |
| direct_cr/inverse | 226.424, 248.809, 313.615, 237.090, 223.370 | 48.103, 45.283, 59.406, 49.698, 45.427 |
| direct_cr/ziggurat | 225.828, 214.504, 240.710, 219.149, 244.112 | 45.214, 44.322, 46.010, 44.356, 45.029 |
| sorting_direct/inverse | 350.999, 331.642, 306.652, 292.046, 304.598 | 67.005, 61.463, 63.112, 64.359, 62.316 |
| sorting_direct/ziggurat | 325.262, 302.396, 294.576, 296.351, 318.431 | 62.245, 61.537, 64.469, 63.365, 61.907 |
| rssa_cr/inverse | 51.035, 47.708, 57.059, 46.310, 54.390 | 9.981, 9.723, 10.536, 10.522, 16.662 |
| rssa_cr/ziggurat | 57.828, 46.288, 57.787, 70.613, 54.102 | 9.345, 9.251, 9.428, 10.377, 12.026 |
| Julia/Direct | 650.715, 644.351, 790.811, 653.915, 692.799 | 155.618, 156.948, 148.153, 163.657, 147.644 |
| Julia/SortingDirect | 215.778, 219.712, 252.367, 229.460, 227.661 | 57.732, 48.967, 177.664, 42.335, 170.535 |
| Julia/DirectCR | 207.726, 217.342, 253.546, 205.346, 219.487 | 66.488, 53.219, 65.510, 65.939, 70.338 |
| Julia/NRM | 304.588, 344.930, 349.043, 313.741, 314.898 | 91.717, 160.226, 72.232, 78.281, 68.848 |
| Julia/RSSACR | 40.381, 43.108, 48.212, 40.120, 36.747 | 18.763, 39.494, 32.478, 21.495, 26.748 |
| Julia/SortingDirect-fresh | 388.452, 395.423, 599.450, 381.743, 386.877 | 168.831, 78.991, 68.959, 213.077, 59.731 |

### W6/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 24430.257, 23066.909, 22747.589, 24018.344, 22939.256 | 4934.719, 4450.573, 4129.645, 4473.062, 4259.901 |
| direct/ziggurat | 24509.848, 24830.829, 24976.558, 23867.859, 24733.992 | 4711.015, 4595.102, 4628.570, 4665.154, 4464.662 |
| cached_direct/inverse | 4285.939, 4101.410, 4039.155, 3962.707, 4219.454 | 606.470, 603.424, 604.130, 601.021, 599.921 |
| cached_direct/ziggurat | 4408.966, 4170.929, 3896.213, 3967.140, 4306.827 | 587.916, 598.074, 617.313, 596.220, 606.071 |
| direct_cr/inverse | 4052.781, 3692.933, 3526.067, 3459.692, 3588.430 | 676.790, 722.892, 715.855, 725.295, 722.429 |
| direct_cr/ziggurat | 4271.307, 3446.370, 3443.194, 3457.772, 3516.931 | 666.472, 703.189, 696.955, 707.462, 701.049 |
| sorting_direct/inverse | 5772.261, 4964.496, 5451.818, 4924.830, 5064.538 | 919.457, 941.506, 980.308, 1002.133, 971.630 |
| sorting_direct/ziggurat | 5335.212, 4664.121, 5036.879, 4784.073, 5458.743 | 978.483, 977.330, 992.857, 978.827, 965.366 |
| rssa_cr/inverse | 880.976, 759.408, 833.940, 742.565, 865.508 | 148.199, 162.033, 164.519, 158.430, 147.932 |
| rssa_cr/ziggurat | 879.845, 768.128, 922.326, 747.057, 767.220 | 141.634, 238.696, 147.402, 143.023, 138.828 |
| Julia/Direct | 11356.974, 10645.830, 10502.531, 10256.248, 10295.594 | 2223.903, 2184.773, 2232.516, 2089.236, 2000.755 |
| Julia/SortingDirect | 3140.812, 3487.370, 3222.931, 2793.272, 2990.835 | 503.036, 486.744, 526.778, 527.647, 496.463 |
| Julia/DirectCR | 4096.317, 3751.743, 3345.109, 3219.504, 3414.713 | 688.795, 632.817, 693.694, 702.828, 689.957 |
| Julia/NRM | 5708.608, 5009.541, 4852.382, 5055.557, 4968.774 | 947.584, 924.323, 938.048, 929.535, 890.018 |
| Julia/RSSACR | 887.396, 692.599, 664.428, 622.152, 669.604 | 141.267, 135.583, 144.496, 134.222, 150.260 |
| Julia/SortingDirect-fresh | 7738.590, 5919.567, 6300.938, 5590.521, 6109.619 | 1413.107, 1344.340, 1379.807, 1294.531, 1490.332 |

### W7/original

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 100.339, 145.431, 87.845, 95.632, 87.647 | 21.343, 21.725, 20.766, 20.288, 21.192 |
| direct/ziggurat | 109.406, 120.743, 109.872, 107.508, 107.959 | 24.596, 22.039, 21.948, 21.829, 21.491 |
| cached_direct/inverse | 28.348, 19.233, 19.138, 23.942, 19.316 | 4.317, 3.791, 4.902, 4.285, 3.774 |
| cached_direct/ziggurat | 21.267, 23.346, 20.913, 22.100, 20.454 | 6.278, 3.865, 4.115, 4.628, 3.895 |
| direct_cr/inverse | 30.119, 26.058, 17.719, 20.711, 17.768 | 4.559, 4.066, 4.308, 5.804, 4.015 |
| direct_cr/ziggurat | 23.891, 20.950, 25.307, 21.401, 19.828 | 3.902, 3.766, 4.017, 3.845, 4.004 |
| sorting_direct/inverse | 26.139, 23.096, 26.156, 20.344, 21.987 | 5.523, 4.789, 5.221, 5.468, 4.586 |
| sorting_direct/ziggurat | 25.568, 19.484, 18.544, 20.179, 21.319 | 6.429, 4.294, 7.490, 5.568, 4.623 |
| rssa_cr/inverse | 23.499, 19.753, 22.331, 18.818, 19.126 | 10.056, 3.828, 4.853, 43.302, 3.807 |
| rssa_cr/ziggurat | 19.845, 22.470, 19.572, 18.260, 19.526 | 8.084, 3.641, 5.185, 6.225, 3.476 |
| Julia/Direct | 48.099, 47.116, 46.777, 42.157, 40.204 | 15.597, 15.064, 16.222, 15.582, 15.372 |
| Julia/SortingDirect | 19.920, 19.892, 18.697, 16.634, 16.502 | 10.347, 8.572, 22.865, 17.244, 11.252 |
| Julia/DirectCR | 19.729, 22.118, 17.517, 16.963, 14.574 | 8.074, 8.511, 8.387, 10.123, 14.465 |
| Julia/NRM | 22.509, 20.188, 17.318, 18.151, 19.446 | 11.077, 8.539, 9.016, 10.411, 8.248 |
| Julia/RSSACR | 15.561, 15.959, 13.567, 14.018, 13.597 | 9.939, 17.389, 8.889, 9.007, 7.513 |
| Julia/SortingDirect-fresh | 72.825, 90.554, 62.668, 84.236, 130.362 | 19.974, 13.863, 23.076, 24.424, 12.127 |

### W7/throughput

| Engine / waiting time | 1 worker, five observations (ms) | 8 workers, five observations (ms) |
| --- | --- | --- |
| direct/inverse | 1590.762, 1470.634, 1502.941, 1652.365, 1399.340 | 310.096, 325.811, 311.579, 310.496, 305.444 |
| direct/ziggurat | 1643.409, 1661.368, 1574.602, 1681.217, 1538.161 | 318.978, 315.927, 313.481, 317.432, 318.402 |
| cached_direct/inverse | 353.182, 322.512, 312.802, 393.356, 320.747 | 55.961, 53.173, 55.535, 54.515, 55.061 |
| cached_direct/ziggurat | 308.883, 330.629, 312.776, 347.994, 322.623 | 60.843, 49.936, 59.098, 50.765, 51.299 |
| direct_cr/inverse | 348.507, 320.264, 294.758, 329.843, 306.123 | 59.107, 57.817, 58.583, 57.143, 58.632 |
| direct_cr/ziggurat | 310.143, 326.710, 269.870, 392.506, 290.371 | 53.595, 52.705, 52.971, 52.827, 55.625 |
| sorting_direct/inverse | 392.935, 348.159, 320.789, 383.987, 360.961 | 67.013, 66.398, 66.633, 67.173, 68.462 |
| sorting_direct/ziggurat | 301.374, 342.315, 319.813, 332.116, 340.081 | 92.456, 62.931, 64.514, 63.587, 66.728 |
| rssa_cr/inverse | 270.086, 293.791, 283.164, 310.602, 309.850 | 52.223, 53.065, 55.258, 51.878, 52.416 |
| rssa_cr/ziggurat | 294.754, 304.037, 266.006, 352.016, 279.748 | 46.144, 46.370, 47.671, 46.086, 52.984 |
| Julia/Direct | 697.242, 666.690, 687.416, 779.787, 861.538 | 169.905, 159.013, 165.367, 164.232, 169.113 |
| Julia/SortingDirect | 260.325, 257.731, 266.245, 277.927, 286.625 | 44.608, 49.602, 46.549, 43.942, 43.801 |
| Julia/DirectCR | 254.485, 238.716, 248.936, 230.685, 255.264 | 54.135, 51.757, 62.756, 62.310, 50.518 |
| Julia/NRM | 302.200, 310.942, 292.814, 303.731, 297.659 | 62.430, 59.799, 59.716, 58.917, 63.999 |
| Julia/RSSACR | 224.408, 214.461, 227.299, 209.275, 220.198 | 41.571, 41.810, 42.331, 43.022, 41.517 |
| Julia/SortingDirect-fresh | 1041.442, 978.436, 1195.435, 1008.967, 1063.807 | 245.605, 253.987, 245.708, 381.871, 236.191 |

## Provenance and implementation log

1. Fetched GitHub main and created an isolated worktree. Built the current fast-mode release wheel and CLI without changing the numerical engine.
2. Installed official Julia 1.13.3 with SHA-256 verification, instantiated the repository's pinned Catalyst/JumpProcesses environment, and completed precompilation.
3. Declared models, methods, thread counts, sample sizes, statistical bounds, timing boundaries and replay requirements before experiments.
4. Added a persistent Julia ensemble driver and Python comparison harness; ran the baseline gates, rate probes, distribution/conservation checks and replay checks.
5. Warmed all configurations, collected five rotated rounds and generated this Markdown directly from complete evidence. No remote push or publication was performed.

Evidence SHA-256 (`validation/eight-thread/evidence/comparison.json`): `611b0a0d78ac2b567bf148094139a0ab5a2147952783d1e8b355faf8f89c646f`.

| Model | Frozen JSON SHA-256 |
| --- | --- |
| W1 | `d2fa8a1f5bdc9976c0bce4f36e8e472690f2ecd1ca3c0f3a4cf642059a37b61e` |
| W2 | `d0f6a86b31cd697b821b5ed2bf384323c8eb1328868fe2274777bac0f5e73355` |
| W3 | `69f9ff6bfec8a473b70111caf6ba2d57e9ce6d363f537eac185f0f88097c8637` |
| W4 | `fb9204e42395eb3eeef76de4e2d54bf33611bdd43435eb7a4eb8fda652df8715` |
| W5 | `c8c2671550a2c67562262f1c2ea22da468bf5acbce6a21309094db3e9aba2681` |
| W6 | `4333e9086f36c5d713c16e7f8f12b43ff9163d2709459820bd92161a9c28f7e0` |
| W7 | `5b6c5f5e329e9fe7146b68cb0ac7d54ee503b64b0e91c6d796c438575edf7ca1` |

| Source | SHA-256 |
| --- | --- |
| Cargo.lock | `6e58d9e3d1c2e2bd36cbe37f30cfe1b98ea287a368fc50f13cbe920db993cc4b` |
| bindings/python/Cargo.lock | `9999d565a026401d196860ef5a23b2066e085bc0bb946ae0186033609538cab8` |
| validation/acceleration_workloads.py | `82a6479d06e61dd7b8c674cc876eea17b453b058e2709cdc29339bcb65540a9e` |
| validation/compare_eight_threads.py | `7437610e5a9bd84cf6e0652123758da1d59d81a6b9491d290115eb561e787e47` |
| validation/compare_julia.py | `24897a366971d2bc8746d507f348d32f7cca242303bf3ef0cfb9dd48e19a5524` |
| validation/eight-thread/ensemble_server.jl | `0f45d8ba29cc9b6127173cc7e744a795ed2823cb41406b6f239e6fa0015ce78a` |
| validation/eight-thread/generated_models.jl | `5caeb8a25625cccd4d0c56b185e6bf5002002e40f676452a37c1caf6dd6de8ee` |
| validation/ssa-acceleration/julia/Manifest.toml | `2b7c9970e8bce9f1ca0c8f322ca7e4c20869aa4ec1f67c3a89d9b9df16b54655` |
| validation/ssa-acceleration/julia/Project.toml | `4721301b2e017a177b376b00765982bf67a1e98beadf1b1ce6862fb38e904885` |

Reference API documentation: [JumpProcesses ensemble/problem reuse](https://docs.sciml.ai/JumpProcesses/stable/tutorials/ensembles_and_problem_reuse/), [SciMLBase ensembles](https://docs.sciml.ai/SciMLBase/stable/interfaces/Ensembles/). Thread-copy and RNG behavior was also checked against the installed pinned SciMLBase 3.58.0 source, not inferred from elapsed time.

## One worker, derived on import

The one-worker observations above also give a serial comparison, measured on
the same host and with the same boundaries. For each workload the lowest
one-worker median is selected independently (it can differ from the
eight-worker choice). This replaces the earlier serial *projection* of 1.02×
in [the Catalyst/JumpProcesses comparison](comparison-catalyst-jumpprocesses.md),
which scaled checked-default timings by fast/default speedups from another host.

| Batch | Fast choice | Fast 1 ms | Julia choice | Julia 1 ms | Fast / Julia |
| --- | --- | --- | --- | --- | --- |
| W1/original | cached_direct/ziggurat | 60.695 | Direct | 52.742 | 1.15× |
| W2/original | direct_cr/ziggurat | 106.051 | DirectCR | 99.262 | 1.07× |
| W3/original | cached_direct/ziggurat | 256.420 | SortingDirect | 300.133 | 0.85× |
| W4/original | rssa_cr/ziggurat | 137.658 | RSSACR | 103.908 | 1.32× |
| W5/original | direct_cr/ziggurat | 55.878 | DirectCR | 49.409 | 1.13× |
| W6/original | rssa_cr/inverse | 51.035 | RSSACR | 40.381 | 1.26× |
| W7/original | cached_direct/inverse | 19.316 | RSSACR | 14.018 | 1.38× |
| W1/throughput | cached_direct/ziggurat | 223.838 | Direct | 215.211 | 1.04× |
| W2/throughput | direct_cr/ziggurat | 458.207 | DirectCR | 372.578 | 1.23× |
| W3/throughput | cached_direct/ziggurat | 1067.678 | SortingDirect | 1187.266 | 0.90× |
| W4/throughput | rssa_cr/ziggurat | 490.796 | RSSACR | 423.411 | 1.16× |
| W5/throughput | direct_cr/ziggurat | 916.916 | DirectCR | 837.709 | 1.09× |
| W6/throughput | rssa_cr/ziggurat | 768.128 | RSSACR | 669.604 | 1.15× |
| W7/throughput | rssa_cr/inverse | 293.791 | RSSACR | 220.198 | 1.33× |

Geometric mean of fast / Julia at one worker: **1.15× for the original
batches** (fast lower on 1/7; range 0.85–1.38×) and
**1.12× for 256-path batches** (fast lower on 1/7; range
0.90–1.33×). With one worker Julia's best method is faster on six of seven
workloads. On the short original batches the selected fast methods scale
5.07–6.13× from one to eight workers and Julia's 1.51–5.78×, which reverses
the order; the 256-path batches are close to parity. The measurements do not
isolate the cause. The same limits apply as above: five observations per
configuration on one shared host, different timing boundaries, and best
observed choices rather than an automatic selector.
