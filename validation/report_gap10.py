"""Render the ten-round decisions and direct Julia comparison, without retiming."""

import csv
import json
import math
from pathlib import Path
import re
import statistics
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "validation/julia-gap10"


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


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


def label(key):
    return (
        key.replace("Baseline/", "")
        .replace("/prepared", "")
        .replace("legacy_cached", "legacy cached")
        .replace("legacy_direct", "legacy Direct")
        .replace("rssa_cr", "RSSACR")
        .replace("direct_cr", "DirectCR")
        .replace("sorting_direct", "SortingDirect")
        .replace("cached_direct", "cached Direct")
        .replace("Julia/", "")
    )


def main():
    r = json.loads((FOLDER / "comparison.json").read_text())
    if not r.get("complete") or not r["passed"]:
        raise RuntimeError("comparison incomplete")
    timings = r["timings"]
    names = list(timings)
    selection = {}
    for name, rows in timings.items():
        groups = dict(
            baseline=[k for k in rows if k.startswith("Baseline/")],
            final=[
                k
                for k in rows
                if k.endswith("/prepared") and not k.startswith("Baseline/")
            ],
            julia=[k for k in rows if k.startswith("Julia/")],
        )
        selection[name] = {
            kind: min(keys, key=lambda k: rows[k]["median_seconds"])
            for kind, keys in groups.items()
        }
    data = []
    gaps = []
    basegaps = []
    gains = []
    for name in names:
        rows = timings[name]
        chosen = selection[name]
        base, final, julia = [
            rows[chosen[k]]["median_seconds"] for k in ["baseline", "final", "julia"]
        ]
        gains.append(base / final)
        gaps.append(final / julia)
        basegaps.append(base / julia)
        data.append(
            [
                name,
                label(chosen["baseline"]),
                f"{base:.5f}",
                label(chosen["final"]),
                f"{final:.5f}",
                f"{base / final:.3f}×",
                label(chosen["julia"]),
                f"{julia:.5f}",
                f"{final / julia:.3f}×",
            ]
        )
    summary = table(
        [
            "Model",
            "Baseline best",
            "s",
            "Final best",
            "s",
            "Gain",
            "Julia best",
            "s",
            "Final / Julia",
        ],
        data,
    )
    decisions = [
        json.loads((FOLDER / f"C{i:02}.json").read_text()) for i in range(1, 11)
    ]
    decision_table = table(
        ["Iteration", "Decision", "Target gain", "All-row gain"],
        [
            [
                v["iteration"],
                "Keep" if v["accepted"] else "Revert",
                f"{v['target_speedup']:.4f}×",
                f"{v['guard_speedup']:.4f}×",
            ]
            for v in decisions
        ],
    )
    rows = []
    for name, v in r["validation"].items():
        rows.append(
            [
                name,
                len(v["rows"]),
                f"{max(x['cdf_error'] for x in v['rows']):.5f}",
                f"{v['rows'][0]['cdf_bound']:.5f}",
                "Pass" if v["passed"] else "FAIL",
            ]
        )
    validation = table(
        ["Model", "Checks", "Largest CDF distance", "Fixed bound", "Result"], rows
    )
    tests = sum(
        map(
            int,
            re.findall(
                r"test result: ok\. (\d+) passed",
                (FOLDER / "final/release-tests.log").read_text(),
            ),
        )
    )
    laws = json.loads((FOLDER / "final/check_alternative_ssa.json").read_text())
    n_observations = sum(
        len(v["seconds"]) for rows in timings.values() for v in rows.values()
    )
    with (FOLDER / "all-method-timings.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "workload",
                "engine_boundary",
                "median_seconds",
                "minimum_seconds",
                "maximum_seconds",
                "all_seconds",
            ]
        )
        for name, rows in timings.items():
            for key, row in rows.items():
                writer.writerow(
                    [
                        name,
                        key,
                        row["median_seconds"],
                        row["min_seconds"],
                        row["max_seconds"],
                        json.dumps(row["seconds"]),
                    ]
                )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 5),
        gridspec_kw={"width_ratios": [4, 3]},
        layout="constrained",
    )
    for ax, models in zip(axes, [names[:4], names[4:]], strict=True):
        for offset, (kind, title, color) in enumerate(
            [
                ("baseline", "GitHub main", "#64748b"),
                ("final", "Final CRN", "#0284c7"),
                ("julia", "Fastest Julia", "#d97706"),
            ]
        ):
            entries = [timings[n][selection[n][kind]] for n in models]
            y = [v["median_seconds"] for v in entries]
            ax.errorbar(
                [i + (offset - 1) * 0.22 for i in range(len(models))],
                y,
                yerr=[
                    [v["median_seconds"] - v["min_seconds"] for v in entries],
                    [v["max_seconds"] - v["median_seconds"] for v in entries],
                ],
                fmt="o",
                capsize=4,
                color=color,
                label=title,
            )
        ax.set_xticks(range(len(models)), models)
        ax.set_yscale("log")
        ax.set_ylabel("Prepared batch seconds · log scale")
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("64 paths · t=100" if models[0] == "W1" else "16 paths · t=20")
    axes[0].legend(fontsize=9)
    fig.suptitle(
        "Ten further Julia-gap iterations\nSame host, five warm rounds; medians and observed ranges",
        fontsize=13,
    )
    plot = ROOT / "docs/images/julia-gap10.png"
    fig.savefig(plot, dpi=180)
    plt.close(fig)
    kept = [v["iteration"] for v in decisions if v["accepted"]]
    diagnostics = json.loads((FOLDER / "regression-diagnostics.json").read_text())
    regression_table = table(
        ["Prepared row", "Original baseline / final", "Diagnostic paired gain"],
        [
            [
                label(key),
                f"{v['original_speedup']:.3f}×",
                f"{v['diagnostic_speedup']:.3f}×",
            ]
            for key, v in diagnostics["python"].items()
        ],
    )
    native_table = table(
        ["Native build", "DirectCR / inverse seconds", "DirectCR / Ziggurat seconds"],
        [
            [
                build,
                f"{diagnostics['native'][build]['median_seconds']:.5f}",
                f"{diagnostics['native_ziggurat'][build]['median_seconds']:.5f}",
            ]
            for build in ["base", "C05", "C06", "C09"]
        ],
    )
    text = f"""# Ten further CRN-craft / Julia-gap iterations

Measured 2 November 2026. GitHub base: `add6cca0fb146e8a40418432e963219514086024`
(merged PR #7). Julia {r["julia_setup"]["versions"]["Julia"]}, Catalyst
{r["julia_setup"]["versions"]["Catalyst"]}, JumpProcesses
{r["julia_setup"]["versions"]["JumpProcesses"]}. This is a direct same-host rerun,
including a separately built GitHub-main Python wheel; it replaces no historical
measurement and uses no cross-host projection.

Ten candidates were implemented and tested in sequence. **{len(kept)} were kept:
{", ".join(kept)}.** The fastest final CRN median improves by up to
**{max(gains):.3f}×** over the fastest baseline choice. The geometric-mean
CRN/fastest-Julia gap is **{geometric(gaps):.3f}×**, versus
**{geometric(basegaps):.3f}×** for GitHub main on this same host. The aspirational
1.25× target was **{"met" if geometric(gaps) <= 1.25 else "missed"}**.

{summary}
![Same-host prepared comparison](images/julia-gap10.png)

Ratios use unrounded medians. Final / Julia >1 means Julia was faster.
Best methods are selected from twelve prepared CRN options and five Julia
algorithms. These are five observations on one shared host; small differences
are not a precise ranking. All methods, including slower variants, are retained
in [the complete timing CSV](../validation/julia-gap10/all-method-timings.csv)
and [raw report](../validation/julia-gap10/comparison.json).
The fastest W5 and W7 choices were respectively 1.7% and 1.4% slower than
baseline in this run. Julia remained faster on every workload.

## Exposed regressions and follow-up

The all-method run exposed six prepared rows with original baseline/final
ratios below 0.95. A separate diagnostic used five fresh alternating pairs,
two repeats per side, after warmup. It did not replace any primary observation
or change an iteration decision. Original ratios below use the ratio of the
original medians; diagnostic ratios use the median paired speedup.

{regression_table}
Five of six drops did not reproduce at the 5% level. W4 DirectCR/Ziggurat did:
0.878× throughput, approximately 13.9% more elapsed time. This is a retained
tradeoff outside the declared eleven-row gate, even though W4's best RSSACR
choice improved by 1.350×. Choose RSSACR for the tested W4 workload; these
changes do not establish a win for every algorithm or wait sampler.

To investigate attribution, the same W4 DirectCR paths were run through the
native prepared benchmark at the base and each retained commit (three samples
per binary, median below; inverse and Ziggurat tested separately):

{native_table}
These sequential native samples do not reproduce a final/base slowdown, and
the exact path digest and event count are identical across builds. The main
and diagnostic Python measurements include serialization and decoding, unlike
this native boundary. Generated-code layout, allocation and host variability
are possible explanations, not measured causes. Without a profiler, the
persistent Python-boundary regression remains unattributed. The patch retains
it explicitly; the original complete CSV includes every slower row.

[Diagnostic code](../validation/diagnose_gap10_regressions.py) and
[all diagnostic observations](../validation/julia-gap10/regression-diagnostics.json)
are included. They are follow-up evidence, not extra attempts to satisfy the
fixed candidate acceptance rule.

## Ten measured iterations

The [protocol](decisions/julia-gap10-2026-11-02.md) was committed before candidate
measurements. Each candidate passed the Rust suite and 35 numerical-identity
checks, then five paired rounds on eleven representative workload/method rows.
Targets and thresholds were not changed after seeing timings. Rejected changes
were reverted, with their diffs and observations kept.

{decision_table}
C05 avoids a third mass-action bound evaluation when the enabled lower corner
already gives the minimum positive rate. C06 shares successful custom-law ranges
only for bit-identical law parameters, regulator identity and feasible integer
interval; distinct availability gates are still checked. C09 compiles elementary
shapes into enum payloads, retaining the original rounded arithmetic and checked
general fallback. All existing stream/sampler identities remain in place.
C01/C02 polling changes, C03/C04 tree changes, C07 bucket lookup, C08 contiguous
SortingDirect and C10 Thin LTO did not satisfy the fixed rule on this host.

## Why Julia is faster

[The source investigation](campaign/julia-gap10-source-investigation.md) cites
the locally inspected pinned sources and their recorded hashes. Julia uses
compact mutable group sums, typed inlined mass-action loops, endpoint bounds for
appropriate laws, Xoshiro/Ziggurat, and learned SortingDirect order retained
between solves. CRN retains checked population updates, canonical sums and
resolution checks, cooperative stops, counted sampling budgets and worker-count
replay. Some mechanisms can be adopted without changing those guarantees;
source inspection alone does not establish that an implementation is faster.
The seven reverted experiments demonstrate that distinction.

No instruction profiler is installed here. Old callgrind percentages in the
attached log describe another host; they are not claimed as current attribution.
This campaign measures end-to-end effects and records evaluation/work counters.

## Validation before timing

The final engine passed **{tests} release tests** on Rust 1.88 and the Rust1.85
minimum, core/binding clippy with warnings denied, formatting, strict contract and
installed-package checks, legacy SSA/CME identities, and the dependency gate.
All **210 independent exact-law marginal comparisons** passed at 100,000 paths
per fixture/variant (7 million paths per full pass). Normal/-O reports match
byte for byte. The fixed family alpha is .001, epsilon
{laws["dkw_epsilon"]:.8f}; largest observed error
{max(v["max_cdf_error"] for v in laws["checks"]):.8f}. Power checks reject the
wrong homodimer convention. No acceptance tolerance was relaxed.

The Julia translation ran 512 paths per final engine/model with distinct seeds,
17 engines and seven workloads: 60,928 validation paths, 6,946 fixed comparisons.
W1/W2 use exact Poisson laws; other marginals use the legacy cached reference.
The same family-alpha .01 DKW rules are retained, including the conservative
two-sample triangle bound. Every timed baseline/final method also returns
identical numerical results; saved work counters are allowed to fall when work
is genuinely avoided.

{validation}
All model propensity/stoichiometric probes and conservation checks passed.
Compressed terminal samples and their hashes are retained. These finite checks
support the tested models; they do not establish arbitrary-model or rare-tail
accuracy. Known steep-Hill RSSACR precision limitations remain documented.

## Timing boundaries and limitations

The run includes **{n_observations:,} timing observations**: 24 final Rust API rows
(one-shot and prepared), twelve prepared baseline rows, and five Julia rows on
each of seven workloads over five rotated/reversed rounds. All builds and heavy
checks finished before warm timing; no observation was removed. Each engine's
replay digest remained stable. Baseline/final prepared calls include native
simulation, result serialization and Python decoding. The baseline lives in an
isolated process; its IPC is outside its own timer. Julia times serial solves,
saved-state checks and final reduction; problem construction, pristine batch
copy and IPC are outside its warm timer. No archive writing is timed.

One numerical worker; seed4242; original frozen W1–W4 (64 paths,t100) and larger
W5–W7 (16 paths,t20), saving initial/final states. These are application boundaries,
not identical allocation/check pipelines. Julia starts each batch from a pristine
problem and keeps its learned order inside the batch; CRN resets per replicate.
The host is {r["environment"]["cpu_model"]}, affinity
{r["environment"]["cpu_affinity"]}, quota `{r["environment"]["cpu_quota"]}`.
Julia's fresh timing-process setup took {r["julia_process_to_ready_seconds"]:.2f}s
for all seven models/five aggregators, not a per-model charge. Compiler/package
installation costs are outside warm timing. No multithread, GPU, other-host or
memory-scaling claim is made.

[Reproduction and raw evidence](../validation/julia-gap10/README.md),
[implementation log](campaign/julia-gap10-implementation-log-2026-11-02.md), and
[all final gate outputs](../validation/julia-gap10/final/summary.json) accompany
the patch. No remote was pushed and no package was published.
"""
    (ROOT / "docs/comparison-julia-gap10.md").write_text(text)
    summary_json = dict(
        kept=kept,
        geometric_gap=geometric(gaps),
        baseline_geometric_gap=geometric(basegaps),
        largest_best_method_speedup=max(gains),
        target_met=geometric(gaps) <= 1.25,
        rows=data,
        observations=n_observations,
        tests=tests,
        regression_diagnostics=diagnostics["python"],
    )
    (FOLDER / "summary.json").write_text(json.dumps(summary_json, indent=2) + "\n")
    print(json.dumps(summary_json, indent=2))


if __name__ == "__main__":
    main()
