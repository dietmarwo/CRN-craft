"""Render complete optimization costs/quality without retiming or dropping seeds."""

import json
import math
import re
import statistics
from parameter_optimization_models import ROOT, FOLDER
from parameter_optimization import ENGINES

LABELS = dict(
    crn_fresh_fast="CRN fresh, fast",
    crn_reuse_fast="CRN reusable, fast",
    crn_generated_fast="CRN generated, fast",
    rebop_expression="rebop expression",
    rebop_rebuild_lma="rebop rebuilt LMA",
    rebop_macro="rebop macro",
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


def median(rows, key):
    return statistics.median(r[key] for r in rows)


def main():
    report = json.loads((FOLDER / "evidence/optimization.json").read_text())
    if not report["complete"] or not report["passed"]:
        raise RuntimeError(
            "Complete successful observations required; failed evidence remains unchanged"
        )
    manifest = json.loads((FOLDER / "fixtures.json").read_text())
    summary, per_case, quality = [], [], []
    for fixture in manifest["rows"]:
        case = fixture["case"]
        rows = [r for r in report["rows"] if r["case"] == case]
        groups = {e: [r for r in rows if r["engine"] == e] for e in ENGINES}
        data = []
        for engine, observations in groups.items():
            data.append(
                [
                    LABELS[engine],
                    f"{median(observations, 'optimization_seconds'):.4f}",
                    f"{median(observations, 'heldout_seconds'):.4f}",
                    f"{median(observations, 'total_with_compilation_seconds'):.4f}",
                    f"{median(observations, 'exact_moment_loss'):.6f}",
                    f"{median(observations, 'heldout_loss'):.6f}",
                ]
            )
        per_case.append(
            "## "
            + case
            + "\n\n"
            + table(
                [
                    "Route",
                    "Optimizer s",
                    "Held-out s",
                    "First-run total s",
                    "Exact loss",
                    "Held-out loss",
                ],
                data,
            )
        )
        reuse, generated = groups["crn_reuse_fast"], groups["crn_generated_fast"]
        reuse_time, gen_time = (
            median(reuse, "optimization_seconds"),
            median(generated, "optimization_seconds"),
        )
        extra_build = (
            median(generated, "compilation_seconds")
            + median(generated, "generation_seconds")
            - median(reuse, "compilation_seconds")
        )
        saved = reuse_time - gen_time
        summary.append(
            dict(
                case=case,
                generated_optimizer_speedup=reuse_time / gen_time,
                reuse_over_fresh=median(
                    groups["crn_fresh_fast"], "optimization_seconds"
                )
                / reuse_time,
                macro_over_reuse=median(groups["rebop_macro"], "optimization_seconds")
                / reuse_time,
                conservative_bundle_break_even_runs=math.ceil(extra_build / saved)
                if saved > 0
                else None,
            )
        )
        for e, obs in groups.items():
            quality.append(
                [
                    case,
                    LABELS[e],
                    f"{min(r['exact_moment_loss'] for r in obs):.6f}–{max(r['exact_moment_loss'] for r in obs):.6f}",
                    f"{max(max(r['relative_parameter_errors']) for r in obs):.1%}",
                    sum(r["actual_evaluations"] for r in obs),
                ]
            )
    preparation = []
    for fixture in manifest["rows"]:
        rows = [r for r in report["rows"] if r["case"] == fixture["case"]]
        grouped = {
            e: [r for r in rows if r["engine"] == e]
            for e in ["crn_fresh_fast", "crn_reuse_fast"]
        }
        costs = {
            e: statistics.median(
                sum(v["bind_seconds"] for v in r["trace"]) / r["actual_evaluations"]
                for r in group
            )
            for e, group in grouped.items()
        }
        preparation.append(
            [
                fixture["case"],
                f"{costs['crn_fresh_fast'] * 1e6:.1f}",
                f"{costs['crn_reuse_fast'] * 1e6:.1f}",
                f"{costs['crn_fresh_fast'] / costs['crn_reuse_fast']:.2f}×",
            ]
        )
    campaign = []
    for engine in ENGINES:
        rows = [r for r in report["rows"] if r["engine"] == engine]
        build = rows[0]["compilation_seconds"] + rows[0]["generation_seconds"]
        campaign.append(
            [
                LABELS[engine],
                f"{build:.4f}",
                f"{sum(r['optimization_seconds'] for r in rows):.4f}",
                f"{sum(r['measured_run_seconds'] for r in rows) + build:.4f}",
            ]
        )
    tests = sum(
        map(
            int,
            re.findall(
                r"test result: ok\. (\d+) passed",
                (FOLDER / "evidence/final/release-tests.log").read_text(),
            ),
        )
    )
    calls = sum(r["actual_evaluations"] for r in report["rows"])
    builds = report["builds"]["workers"]
    text = (
        f"""# Reusable parameters, generated kernels and rebop: full optimization

Measured 2 November 2026, America/La_Paz. Base GitHub main `c16fcc5`, including
PR #8's optional `method.fast`. Actual engine versions: CRN-craft 0.2.0.dev1,
rebop {report["versions"]["rebop"]}, fcmaes-rust {report["versions"]["fcmaes-rust"]}.
This is a same-host experiment; it does not project the previous Julia timings.

All **{len(report["rows"])} complete optimization runs** finished: three fixed
models, six execution routes and three optimizer seeds. There were **{calls:,}
actual objective evaluations**. All nine three-route CRN optimizer trace
comparisons were identical, point for point and loss for loss. Compiled code
reused each topology across all candidates; coefficients remained runtime data.

The reusable interface reduced total optimizer time by about 7% on sparse33
and 1.73× on binding12 (per seed 1.69–1.77×); birth/death was essentially flat.
Three seeds with one repetition each cannot resolve smaller effects: on
birth/death, where fresh and reusable CRN simulate identical paths, their ratio
still ranged 0.85–1.16× across seeds. Generated propensities gave a median
1.18× on birth/death (per seed 0.99–1.27×), 0.92–0.97× on sparse33 (median
about 8.5% slower) and 0.73–1.01× on binding12; none of these is resolved.
Generation did not repay its bundle compilation in these 128-evaluation runs.
Rebop macros were faster than reusable CRN by 3.44× on birth/death (per seed
2.99–4.35×) and 1.72× on binding12 (1.64–1.77×), and approximately tied on
sparse33 (0.94–1.00×), using optimizer time. Across all nine runs per route,
complete cost including one build was 13.21s for reusable CRN, 17.53s for
generated CRN and 12.88s for rebop macros; the 2.5% difference between reusable
CRN and rebop macros is within the run-to-run variation. These aggregated
results also include different cross-engine search paths and final points.

The new parameter interface and optional generator are runnable. Numerical
speed and optimization quality are reported separately below; the generated
backend specializes propensities and retains the core solver.

"""
        + "\n".join(per_case)
        + "\n## Binding cost\n\nThe following per-candidate native preparation costs are medians of each run's mean; the optimizer tables also include simulation and application overhead.\n\n"
        + table(
            ["Case", "Fresh preparation µs", "Reusable binding µs", "Preparation gain"],
            preparation,
        )
        + """
Times/losses are medians over the three declared optimizer seeds. Lower is
better. “Optimizer s” includes all objective calls, transport, reduction and
optimizer overhead. “First-run total” includes process setup/teardown, source
generation, the entire relevant bundle's dependency-cache-warm compilation,
optimizer time, exact-oracle assessment and 4096-path held-out validation.
It excludes dependency installation/compilation. Every seed and timing
observation remains in [the raw report](../validation/parameter-optimization/evidence/optimization.json).

## Repeated-use cost

The following totals cover all nine runs per route, charging the bundle's build
once. They illustrate reuse across models/seeds; they are not per-model compiler
cost estimates. Existing applications can also amortize the generic Rust worker
build over arbitrary later models.

"""
        + table(
            [
                "Route",
                "Generation + build s",
                "Nine optimizations s",
                "Nine complete runs + one build s",
            ],
            campaign,
        )
        + """
## Quality and evaluation accounting

Targets use independent Poisson moments for the two linear systems and a
13-state generator/matrix exponential for reversible binding. Each objective
uses 32 trajectories with a fixed training seed; each final point receives 4096
new trajectories on an independent held-out seed. Exact loss evaluates the
fitted point's mathematical moments, without simulation noise. The objective
normalizes squared mean/variance residuals and excludes t=0. Parameter errors
below are the largest relative error over parameters and optimizer seeds.
A small moment loss does not guarantee parameter identifiability.

"""
        + table(
            [
                "Case",
                "Route",
                "Exact loss range",
                "Largest parameter error",
                "Actual objective calls",
            ],
            quality,
        )
        + """
All route/model median exact losses met the predeclared <=.05 aspiration.
Parameter recovery was less precise: worst errors reached 51.9% for CRN binding
and 91.2% for a rebop birth/death run. These are finite-budget fits to noisy
32-path moments; a low moment loss is insufficient evidence of accurate rate
identification. Held-out loss and parameter error remain separate results.

The exact-loss <=.05 aim was aspirational, not a reason to repeat or drop runs.
Cross-engine streams differ, so different final points need not imply a solver
bias or a superior optimizer. Three search seeds and these three synthetic
models support only the reported scope.

## Compilation and boundaries

"""
        + table(
            ["Worker bundle", "Measured rebuild s"],
            [[k, f"{v['compilation_seconds']:.4f}"] for k, v in builds.items()],
        )
        + f"""
Model-independent dependencies were primed consistently. Each worker package
was cleaned and rebuilt separately, with the same standard Cargo release
profile. The generated bundle contains the three model kernels and an arithmetic
guard. A single-model production build could have a different cost. No
clean-machine installation, LTO tuning or other-platform claim is made.

CRN fresh/reusable/generated all use cached Direct, Ziggurat and **fast=True**,
one numerical worker. Fast's state checks, work budgets, stop polling and
accounting remain. Rebop macros use generated model-specific state/reaction
code; rebop dynamic uses runtime expressions or rebuilt mass-action definitions.
Rebop macros specialize the full reaction code and state representation; our
AOT route specializes propensity evaluation within the existing checked solver.
Native startup also includes the common worker framework/reference seed-list
preparation; it is not the minimum overhead of an in-process Rust application.
Dynamic rebop explicitly uses sparse storage on the 33-species model. Both
Python routes construct per-trajectory xarray outputs; native routes communicate
with a persistent worker and return reduced moments. These are complete
application costs with different allocation/check pipelines. Native bind and
simulation times are retained for CRN; event counts are not, so the engines'
different search paths cannot be normalized by work. The CRN worker also
serializes and hashes each result. Ratios of optimization times are not
isolated event-kernel throughput measurements.

The host is {report["environment"]["cpu"]}, affinity allowance
{report["environment"]["affinity"]}, CPU quota `{report["environment"]["quota"] or "unavailable"}`.
All builds and heavy gates finished before optimizer timing. No multithread,
GPU, biological-network generalization or Julia performance claim is made.

## Validation and reproduction

The final engine passed **{tests} release tests** on Rust 1.88 and 1.85,
core/binding/worker clippy with warnings denied, strict schema/CLI, installed
Python, legacy SSA/CME and dependency gates. Both the checked-default and fast
families passed **210 exact-law comparisons each**, 100,000 paths per
fixture/variant, with the existing .001 family alpha and .008046120180018303
DKW epsilon unchanged. Normal/-O outputs match. Independent generated checks
cover 180 bindings/executions, 3150 arithmetic evaluations, cache refresh,
invalid inputs, cancellation, concurrency, worker replay and wrong topology.
The final gate logs record `check_parameter_interface.py` with SHA-256
`9f08cca6…`; the committed file (`8ee7cc7f…`) was edited afterwards and only
`parameter-validation.json` was refreshed with it. The evidence describes the
import commit `85a085b`; later commits change hashed sources.

[API guide](reusable-parameters.md),
[fixed protocol](decisions/reusable-parameters-2026-11-02.md),
[reproduction instructions](../validation/parameter-optimization/README.md) and
[implementation log](campaign/reusable-parameters-implementation-log-2026-11-02.md)
accompany the patch. Failed development outputs remain separate. The checked
mode remains the default. No remote was pushed or package published.
"""
    )
    (ROOT / "docs/comparison-parameter-optimization.md").write_text(text)
    save_path = FOLDER / "evidence/summary.json"
    save_path.write_text(
        json.dumps(
            dict(
                runs=len(report["rows"]),
                objective_calls=calls,
                tests=tests,
                cases=summary,
            ),
            indent=2,
        )
        + "\n"
    )
    print(save_path.read_text())


if __name__ == "__main__":
    main()
