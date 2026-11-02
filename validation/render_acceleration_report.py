"""Render the acceleration comparison from retained evidence without retiming."""

import json
from pathlib import Path
import re
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation/ssa-acceleration"
LABELS = {
    "legacy_direct": "Legacy Direct",
    "legacy_cached": "Legacy cached",
    "direct": "Direct",
    "cached_direct": "Cached",
    "direct_cr": "DirectCR",
    "sorting_direct": "SortingDirect",
    "rssa_cr": "RSSACR",
}
NAMES = {
    "W1": "Immigration–death",
    "W2": "Sparse independent channels",
    "W3": "Repressilator",
    "W4": "Dense synchronized custom rates",
    "W5": "Diffusion chain",
    "W6": "Dense conversion",
    "W7": "Compartment binding/diffusion",
}


def label(key):
    if key.startswith("Julia/"):
        return "Julia " + key.split("/")[1]
    bits = key.split("/")
    return (
        "CRN "
        + LABELS[bits[0]]
        + (
            " · " + bits[1]
            if len(bits) > 1 and bits[1] in ["inverse", "ziggurat"]
            else ""
        )
    )


def table(headers, rows):
    return (
        "\n".join(
            [
                "| " + " | ".join(headers) + " |",
                "|" + "|".join(["---"] * len(headers)) + "|",
            ]
            + ["| " + " | ".join(map(str, row)) + " |" for row in rows]
        )
        + "\n"
    )


def main():
    r = json.loads((DATA / "julia-comparison.json").read_text())
    laws = json.loads((DATA / "final-laws.json").read_text())
    if not r.get("complete") or not r["passed"] or laws["failed"]:
        raise RuntimeError("incomplete evidence")
    timings = r["timings"]
    names = list(timings)

    def median(name, key):
        return timings[name][key]["median_seconds"]

    def fastest(name, predicate):
        return min(
            (k for k in timings[name] if predicate(k)), key=lambda k: median(name, k)
        )

    best_rust = {n: fastest(n, lambda k: k.endswith("/prepared")) for n in names}
    best_julia = {n: fastest(n, lambda k: k.startswith("Julia/")) for n in names}
    best_new = {
        n: fastest(n, lambda k: k.endswith("/prepared") and not k.startswith("legacy_"))
        for n in names
    }
    rows = []
    for n in names:
        a, b = best_rust[n], best_julia[n]
        rows.append(
            [
                n,
                label(a),
                f"{median(n, a):.4f}",
                label(b),
                f"{median(n, b):.4f}",
                f"{median(n, a) / median(n, b):.2f}×",
                f"{median(n, 'legacy_cached/prepared') / median(n, a):.2f}×",
            ]
        )
    summary = table(
        [
            "Workload",
            "Fastest CRN prepared",
            "Seconds",
            "Fastest Julia",
            "Seconds",
            "CRN / Julia",
            "Cached / best CRN",
        ],
        rows,
    )
    medians = []
    keys = [
        k for k in timings["W1"] if k.endswith("/prepared") or k.startswith("Julia/")
    ]
    for group in [names[:4], names[4:]]:
        medians.append(
            table(
                ["Implementation"] + group,
                [
                    [label(k)]
                    + [
                        f"{median(n, k):.4f} [{timings[n][k]['min_seconds']:.4f}–{timings[n][k]['max_seconds']:.4f}]"
                        for n in group
                    ]
                    for k in keys
                ],
            )
        )
    boundary = table(
        [
            "Workload",
            "Best new method",
            "One-shot seconds",
            "Prepared seconds",
            "One-shot / prepared",
        ],
        [
            [
                n,
                label(best_new[n]),
                f"{median(n, best_new[n].replace('/prepared', '/one_shot')):.4f}",
                f"{median(n, best_new[n]):.4f}",
                f"{median(n, best_new[n].replace('/prepared', '/one_shot')) / median(n, best_new[n]):.2f}×",
            ]
            for n in names
        ],
    )
    workloads = table(
        ["ID", "Network", "Species", "Reactions", "Paths", "Horizon"],
        [
            [
                v["name"],
                NAMES[v["name"]],
                v["species"],
                v["reactions"],
                v["timing_paths"],
                v["end"],
            ]
            for v in r["workloads"]["rows"]
        ],
    )
    validation = table(
        ["Workload", "Checks", "Largest CDF distance", "Bound", "Result"],
        [
            [
                n,
                len(v["rows"]),
                f"{max(x['cdf_error'] for x in v['rows']):.5f}",
                f"{v['rows'][0]['cdf_bound']:.5f}",
                "Pass" if v["passed"] else "FAIL",
            ]
            for n, v in r["validation"].items()
        ],
    )
    workrows = []
    for n in names:
        old = timings[n]["cached_direct/inverse/prepared"]
        new = timings[n]["rssa_cr/ziggurat/prepared"]
        workrows.append(
            [
                n,
                f"{old['work']['propensity_evaluations']:,}",
                f"{new['work']['propensity_evaluations']:,}",
                f"{new['work']['bound_updates']:,}",
                f"{new['work']['rejected_candidates']:,}",
                f"{new['events']:,}",
            ]
        )
    worktable = table(
        [
            "Workload",
            "Cached rate evaluations",
            "RSSACR rate evaluations",
            "RSSACR bound updates",
            "Rejected RSSA candidates",
            "Accepted RSSACR events",
        ],
        workrows,
    )
    ns_event = table(
        ["Workload", "Best CRN ns / its event", "Julia ns / CRN-event proxy"],
        [
            [
                n,
                f"{1e9 * median(n, best_rust[n]) / timings[n][best_rust[n]]['events']:.0f}",
                f"{1e9 * median(n, best_julia[n]) / timings[n]['legacy_cached/prepared']['events']:.0f}",
            ]
            for n in names
        ],
    )
    setup = r["julia_setup"]
    tests = sum(
        map(
            int,
            re.findall(
                r"test result: ok\. (\d+) passed",
                (DATA / "final-tests.log").read_text(),
            ),
        )
    )
    gains = {
        n: median(n, "legacy_cached/prepared") / median(n, best_rust[n]) for n in names
    }
    largest = max(gains, key=gains.get)
    ratios = [median(n, best_rust[n]) / median(n, best_julia[n]) for n in names]
    wins = sum(median(n, best_rust[n]) < median(n, best_julia[n]) for n in names)
    # Scientific plot uses the exact medians and observed ranges, log scale.
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 5),
        gridspec_kw={"width_ratios": [4, 3]},
        layout="constrained",
    )
    series = [
        ("legacy_cached/prepared", "Legacy cached", "#64748b"),
        ("best_crn", "Best CRN prepared", "#0284c7"),
        ("best_julia", "Best Julia", "#d97706"),
    ]
    for ax, group in zip(axes, [names[:4], names[4:]], strict=True):
        for offset, (key, title, color) in enumerate(series):
            actual = [
                best_rust[n]
                if key == "best_crn"
                else best_julia[n]
                if key == "best_julia"
                else key
                for n in group
            ]
            vals = [median(n, k) for n, k in zip(group, actual, strict=True)]
            lo = [
                median(n, k) - timings[n][k]["min_seconds"]
                for n, k in zip(group, actual, strict=True)
            ]
            hi = [
                timings[n][k]["max_seconds"] - median(n, k)
                for n, k in zip(group, actual, strict=True)
            ]
            x = [i + (offset - 1) * 0.23 for i in range(len(group))]
            ax.errorbar(
                x, vals, yerr=[lo, hi], fmt="o", capsize=4, color=color, label=title
            )
        ax.set_xticks(range(len(group)), group)
        ax.set_yscale("log")
        ax.set_ylabel("Batch seconds · log scale")
        ax.grid(axis="y", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("64 paths · t=100" if group[0] == "W1" else "16 paths · t=20")
    axes[0].legend(loc="best", fontsize=9)
    fig.suptitle(
        "CRN-craft and Catalyst / JumpProcesses\nFive warm runs; dots are medians, whiskers are observed ranges",
        fontsize=13,
    )
    image = ROOT / "docs/images/julia-comparison-accelerated.png"
    fig.savefig(image, dpi=180)
    plt.close(fig)
    text = f"""# CRN-craft compared with Catalyst.jl and JumpProcesses.jl

Updated **2 November 2026**, after the five-step exact-SSA acceleration campaign.
CRN-craft **0.2.0-dev.1** now has DirectCR, SortingDirect, RSSACR, reusable Python
preparation and a separately versioned exponential Ziggurat sampler. Julia remains
**{setup["versions"]["Julia"]}**, Catalyst **{setup["versions"]["Catalyst"]}**, and
JumpProcesses **{setup["versions"]["JumpProcesses"]}**. Exact source hashes are
retained in the measurement report; the implementation is an unpushed patch over
GitHub base `baed0e4590f93e4c1dcc77cd7f039a2604a1594c`.

**The largest measured improvement over legacy cached Direct is {gains[largest]:.2f}×
on {largest}; W4 improves {gains["W4"]:.2f}×.** Julia's fastest measured method remains
{min(ratios):.2f}–{max(ratios):.2f}× faster than the fastest CRN option across these
seven workloads. The fastest observed CRN prepared median beats Julia on
**{wins} of 7** workloads. All slower variants
and timing observations remain in the tables and raw evidence. The rightmost
column compares the current legacy cached solver with the best CRN option at the
same prepared boundary. Ratios use unrounded medians; CRN / Julia above 1 means
Julia was faster.

{summary}
**Update after the Julia-gap campaign (projected, not remeasured).** A
follow-up campaign ([protocol](decisions/julia-gap-2026-11-02.md),
[log](campaign/julia-gap-implementation-log-2026-11-02.md)) kept five of ten
iterations under predeclared rules: AVX2 ChaCha8 dispatch, an exact Hill-rate
memo, a cached composition envelope, a sums-only lazily checked rate tree and a
joint member/acceptance word. Julia's package hosts were unreachable from that
session, so the table multiplies the CRN medians above by the paired speedups
measured on another host; it is not a Julia measurement, and the tables below
remain the measured comparison.

| Workload | Fastest projected CRN | Projected seconds | Fastest Julia | CRN / Julia |
|---|---|---|---|---|
| W1 | Legacy cached | 0.0749 | Direct | 1.47× |
| W2 | DirectCR | 0.1272 | DirectCR | 1.31× |
| W3 | SortingDirect | 0.3278 | SortingDirect | 1.11× |
| W4 | RSSACR | 0.1946 | RSSACR | 1.83× |
| W5 | DirectCR | 0.0768 | DirectCR | 1.61× |
| W6 | RSSACR | 0.0614 | RSSACR | 1.55× |
| W7 | RSSACR | 0.0187 | RSSACR | 1.41× |

The projected geometric-mean gap is 1.46× (measured 2.02× before). A direct
rerun of `validation/compare_julia_accelerated.py` remains the authoritative
check.

![Warm ensemble times after SSA acceleration](images/julia-comparison-accelerated.png)

The [earlier comparison](comparison-catalyst-jumpprocesses-before-acceleration.md)
is preserved. Its absolute times are historical, not observations combined with
this rerun. The [GillesPy2 report](comparison-gillespy2.md) uses the same W1–W4
models but was not rerun here, so this is not a three-way timing comparison.

## What was implemented

| Step | Delivered change | What the evidence supports |
|---|---|---|
| 1 | Component profiling and immutable `prepare` / `prepare_ssa` handles | Preparation is a small fraction of the original W1–W4 batches; reusable calls preserve results and support concurrent reuse |
| 2 | DirectCR with dependency updates and power-of-two composition/rejection | Removes ordered-prefix maintenance; balanced-tree and rejection costs can dominate small/dense workloads |
| 3 | SortingDirect and a versioned Ziggurat wait sampler | Both execute and replay; changing only the wait sampler does not guarantee an end-to-end gain |
| 4 | RSSACR with certified integer-box bounds, exact point fallbacks and explicit rejected-work limits | Avoids many actual-rate evaluations while preserving gates, absorption and checked state updates |
| 5 | Fresh Julia rerun and three larger connected networks | Shows where algorithm choice pays and where a measured gap remains |

This reimplements relevant algorithmic techniques in Rust; it does not recreate
Julia's compiler, symbolic model generator or complete solver ecosystem.
[Implementation details and API examples](alternative-ssa.md) explain the
published algorithm families, numerical policy and new contracts. The old
Direct/cached interfaces, schemas and byte identities are preserved. New methods
return `ssa-result/0.1` JSON, not a legacy sealed archive.

## Models and timing boundaries

{workloads}
W1–W4 are unchanged. W5 has 256 species, four initial molecules each, and
nearest-neighbor conversion at rate 0.5 in each direction. W6 has 32 species,
32 initial molecules each, with heterogeneous pairwise conversion rates
`0.005 * (1 + (7*i + 11*j) % 8)`. W7 has 32 compartments, each initially
A=24, B=12, C=4; A+B binds at 0.015, C dissociates at 0.5, and A diffuses
between adjacent compartments at 0.1. Frozen JSON and hashes accompany the
report. All channels remain separate.

All batches use one numerical worker, seed 4242 and saves at 0 and the horizon.
Julia uses one advancing Xoshiro stream per batch; CRN uses its pinned per-path
ChaCha8 streams. Matching seed numbers do not imply matching paths. Five rounds
use a deterministic shuffled order, then rotate/reverse it. No build or heavy
validation runs during the final warm rounds. No outlier is removed.

| Boundary | Included | Outside the warm timer |
|---|---|---|
| CRN one-shot | Request encoding, native validation/model preparation, simulation, standard metadata/statistics, result encoding and Python decoding | Model bytes/request dictionary creation; archive I/O |
| CRN prepared | Simulation, standard metadata/statistics, result encoding and Python decoding | Validation and preparation; archive I/O |
| Julia | Serial solve calls with SSAStepper, per-path initialization, saved-time/count checks, final statistics reduction | Catalyst construction/JIT setup, copy of pristine batch problem, JSON IPC |

These are application-level boundaries with different result/allocation policies,
not isolated equal-instruction kernels. Julia reuses its problem within each
batch, including SortingDirect's learned order, and starts each batch from a
pristine copy [1]. CRN resets ordering per replicate to retain worker-count
replay. Prepared calls produce identical results to one-shot calls.

Host: {r["environment"]["cpu_model"]}, Linux x86-64, affinity
{r["environment"]["cpu_affinity"]} logical CPUs and cgroup quota
`{r["environment"]["cpu_quota"]}`. Release Rust 1.88.0, Python 3.12.14,
NumPy {r["environment"]["numpy"]}; Julia and BLAS each use one thread.
This shared-host study has five observations per row, not confidence intervals.

## Warm results

Seconds per batch: **median [minimum–maximum]**. CRN rows below use prepared calls.

{medians[0]}
{medians[1]}
The best-new-method boundary comparison is below. All one-shot/prepared
observations for every method are in the raw report, including slower choices.

{boundary}
Event normalization is secondary. Julia event counts were not instrumented:
its following values use the legacy CRN count as an explicitly labelled proxy.

{ns_event}

## Why the results change

The first-step component measurements found roughly 3 ns per random word and
16–18 ns per inverse exponential draw, but preparation savings were small on
the original W1–W4 batches. On W6, RSSACR/Ziggurat took {median("W6", "rssa_cr/ziggurat/one_shot"):.4f} s one-shot
versus {median("W6", "rssa_cr/ziggurat/prepared"):.4f} s prepared, but RSSACR/inverse moved only from {median("W6", "rssa_cr/inverse/one_shot"):.4f} s to
{median("W6", "rssa_cr/inverse/prepared"):.4f} s, and one-shot calls were faster than prepared ones on W1, W2 and W4;
with five shared-host observations per row, the reuse benefit is not resolved. These diagnostic component timings include their own API overhead;
they are not additive percentages of total simulation time.

The new selectors use checked balanced sums and bounded sampling. They pay for
tree maintenance, accounting and rejection draws, so replacing a compact prefix
scan can make a small or dense model slower. RSSACR instead changes how often
actual propensities are evaluated. The following measured work counts compare
cached inverse sampling with RSSACR/Ziggurat at the same batch size; their random
paths and accepted event counts differ.

{worktable}
W4 is a particularly favorable RSSA case: every event changes all eight species
and all 33 custom death propensities, while bounded rejection can reuse their
envelopes. It is deliberately synthetic and has an equal-count invariant; its
behavior does not predict arbitrary dense-network performance. W5–W7 add
connected sparse, densely connected and compartmental structures but remain
synthetic. The full median tables, not an algorithm name alone, determine the
measured winner for each model.

No compiler-level attribution follows from these timings. Remaining differences
can include selection/summation design, bounds, checked count updates, RNG and
wait generation, generated rate code, metadata and API costs. This campaign
does not isolate the individual contribution of each difference.

## Correctness before speed

The protocol was committed before experiments. The final engine passed **{tests}
release tests**, strict clippy/fmt, legacy path/CME gates, 69 schema/CLI checks,
normal/-O native/Python/prepared/concurrent checks, and the unchanged default
dependency gate. Rust 1.85 MSRV results are retained separately.

Every one of the ten method/wait combinations ran 100,000 paths on each of the
three linear and four nonlinear fixtures: **7,000,000 paths per full pass**.
All **210 marginal comparisons** passed. With family alpha 0.001, the frozen DKW
limit is {laws["dkw_epsilon"]:.8f}; the largest observed error is
{max(v["max_cdf_error"] for v in laws["checks"]):.8f}. Wrong homodimer conventions
were rejected. Normal and Python -O final reports are byte-identical. Additional
million-sample categorical/exponential gates passed, as did forced tail/rejection
limits and independent enumeration of gated integer boxes. No tolerance was
relaxed after a failure.

The Julia comparison ran 512 validation paths per engine/model: **60,928 paths**
across 17 engines and seven workloads. W1/W2 use exact Poisson terminal laws;
other species are checked against legacy cached SSA. With K=6,946 and family
alpha 0.01, the one-sample bound is `sqrt(log(2K/alpha)/(2N))`, and the conservative
two-sample bound is `2*sqrt(log(4K/alpha)/(2N))`. This latter bound applies the
one-sample inequality to both empirical CDFs and uses the triangle inequality.
No independence between species is assumed for the union bound.

{validation}
The "Checks" column counts every engine/species comparison, including rows that
cannot fail. On W3–W7, three CRN engines (`legacy_direct`, `direct/inverse` and
`cached_direct/inverse`) reproduce the legacy cached reference sample exactly at
the shared seed, and `direct/ziggurat` reproduces `cached_direct/ziggurat`; their
distance is 0 by construction (W3 18 of 96 rows, W4 24 of 128, W5 768 of 4096,
W6 96 of 512, W7 288 of 1536). All engines also share the legacy random streams
at that seed, so even different methods are correlated with the reference, which
weakens a two-sample test that assumes independent samples. With N=512 the
two-sample gate detects only CDF differences above 0.24. For W3–W7 these checks
are therefore translation checks of Julia against legacy cached SSA; the
100,000-path exact-law gate above is the evidence for the new methods' laws.
Future runs should use distinct seeds per engine and report distinct-sample
comparisons separately.

Propensity and stoichiometry probes passed on all seven models. W5/W6 preserve
total molecule count; W7 preserves global A+C and each compartment's B+C; W4
retains equal counts. Every timed batch reproduces its own numerical digest.
Compressed sample files retain all individual terminal counts. These conservative
terminal checks are useful translation checks, not proofs about rare tails or
arbitrary model equivalence; the larger 100,000-path laws provide a separate gate.

Catalyst uses ordinary mass-action arrows with the combinatorial convention and
full-rate arrows for custom laws. As in the original comparison, W4's redundant
Julia availability gates are omitted only on its equal-count reachable class.
The Hill and Michaelis–Menten reference bounds are valid for these translated
single-variable monotone laws [2]. CRN's new bounds also handle the distinct-gate
and consumed-repressor cases tested separately. Arbitrary nonmonotone Julia rates
need sound user-supplied bounds; this benchmark does not authorize endpoint-only
bounds for them.

## Startup, scope and retained evidence

The fresh timing process took **{r["julia_process_to_ready_seconds"]:.2f} s** to its
ready response, including **{setup["load_seconds"]:.2f} s** of package imports and
**{setup["model_setup_seconds"]:.2f} s** of setup for all seven models and five
algorithms. These nested figures are not additive and are not per-model startup
costs. Installation/precompilation and Rust/wheel builds are excluded. A long-lived
Julia process amortizes setup. Individual setup/first-use observations are retained.

Catalyst/JumpProcesses still offers a richer symbolic/programmatic modeling
interface and additional exact algorithms such as NRM and Coevolve [2]. This
patch focuses on CRN-craft's four existing built-in rate laws and strict JSON
execution contracts. It does not add arbitrary time-dependent rates, hybrid
jump–ODE models, spatial solvers, distributed execution, GPU support or new FSP/
ODE algorithms. Existing CRN FSP, approximate leaps and optional fitting/ODE
integrations are outside these timings.

One shared host, one worker, seven models and five rounds cannot establish a
universal ranking, parallel scaling, memory scaling or production workload result.
Algorithm choice and bracket tuning can change outcomes. Ziggurat changes paths
and must keep its separate sampler identity. New exact-method JSON has no sealed
archive wrapper yet. Other-platform wheels were not tested.

Reproduce with [validation/ssa-acceleration/README.md](../validation/ssa-acceleration/README.md).
The [raw timing report](../validation/ssa-acceleration/julia-comparison.json)
contains **1,015 observations**; the [validation report](../validation/ssa-acceleration/julia-validation.json)
contains all 6,946 checks and hashes of the sample files, sources and pinned Julia
environment. The [implementation log](campaign/ssa-acceleration-implementation-log-2026-11-02.md)
retains development failures and the absorption follow-up. The intentionally
interrupted validation is kept separately and is not used for these conclusions.
The figure and this document can be regenerated with
`python validation/render_acceleration_report.py`.

Sources checked 2 November 2026; timing and validation values are our measurements:

1. [JumpProcesses problem reuse and SortingDirect replay](https://docs.sciml.ai/JumpProcesses/stable/tutorials/ensembles_and_problem_reuse/).
2. [Jump types, exact algorithms, dependency graphs, rate bounds and original papers](https://docs.sciml.ai/JumpProcesses/stable/jump_types/).
3. [Catalyst model semantics](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/math_models_intro/) and [full-rate DSL arrows](https://docs.sciml.ai/Catalyst/stable/model_creation/dsl_advanced/).
4. [Julia 1.13.3 exponential sampler source](https://github.com/JuliaLang/julia/blob/v1.13.3/stdlib/Random/src/normal.jl); adapted MIT notices accompany the patch.
"""
    (ROOT / "docs/comparison-catalyst-jumpprocesses.md").write_text(text)
    print(summary)
    print(
        "maximum ratio",
        max(median(n, best_rust[n]) / median(n, best_julia[n]) for n in names),
    )


if __name__ == "__main__":
    main()
