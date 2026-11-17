"""Re-creating Figure 1 of Kowald & Kirkwood (2013) with CRN-craft.

Author: Axel Kowald. Adapted to CRN-craft 0.8 by the CRN-craft maintainers.
MIT license, like the rest of the repository.

Paper: Kowald A., Kirkwood T.B.L. (2013). Mitochondrial mutations and aging:
random drift is insufficient to explain the accumulation of mitochondrial
deletion mutants in short-lived animals. Aging Cell 12(4):728-731.
doi:10.1111/acel.12098

Guide: docs/guide/16-published-study.md. Run from the repository root:

    python examples/kk13/kowald_kirkwood_2013.py --size quick           # seconds
    python examples/kk13/kowald_kirkwood_2013.py --size draft           # the guide's tables
    python examples/kk13/kowald_kirkwood_2013.py --size full            # the author's cell counts
    python examples/kk13/kowald_kirkwood_2013.py --size draft --part 1B
    python examples/kk13/kowald_kirkwood_2013.py --figures              # plots only

Every unit of work (one P_mut of Fig 1A, one lifespan of Fig 1B, one
combination of Fig 1C/1D) writes its own result file under
examples/kk13/results/<size>/ and is skipped when that file exists, so a long
run can be stopped and continued.
"""

import argparse
import json
import math
import os
import platform
import time
from pathlib import Path

import numpy as np

import crn_craft as crn

nMolecules = 1000
halfLifeDays = 10.0
degRate = math.log(2) / (24.0 * halfLifeDays)  # 1/h
coxThreshold = 0.6
hoursPerYear = 24.0 * 360.0
nThreads = max(1, (os.cpu_count() or 1) - 2)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIFESPANS = [3, 10, 40, 80, 120]
TARGETS = [0.05, 0.10, 0.15]
FIG1A_PMUT = [1e-4, 1e-5, 1e-6]

# Cell counts and horizons per size. "full" is the author's run (about
# 2.6 h on 18 threads for Fig 1C/1D); "draft" gives the guide's tables in
# well under an hour on two cores; "quick" is the CI check.
SIZES = {
    "quick": {
        "fig1a_cells": {1e-4: 200, 1e-5: 200, 1e-6: 200},
        "fig1a_years": 10,
        "fig1a_compile": False,
        "lifespans": [3, 10],
        "pilot_cells": 100,
        "fine_cells": 200,
        "cd_lifespans": [3],
        "cd_cells": 100,
    },
    "draft": {
        "fig1a_cells": {1e-4: 3000, 1e-5: 3000, 1e-6: 15000},
        "fig1a_years": 120,
        "fig1a_compile": True,
        "lifespans": LIFESPANS,
        "pilot_cells": 500,
        "fine_cells": 1000,
        "cd_lifespans": LIFESPANS,
        "cd_cells": 200,
    },
    "full": {
        "fig1a_cells": {1e-4: 3000, 1e-5: 3000, 1e-6: 15000},
        "fig1a_years": 120,
        "fig1a_compile": True,
        "lifespans": LIFESPANS,
        "pilot_cells": 500,
        "fine_cells": 3000,
        "cd_lifespans": LIFESPANS,
        "cd_cells": 3000,
    },
}
BATCH = 1000  # cells per crn.ensemble call: bounds the memory of kept paths


def buildModel(pmut) -> crn.Model:
    """Fig 1A/1B: the exact lumped model, two events that change the counts."""
    model = crn.Model("mtdna_random_drift")
    model.add_reaction(
        "delta*wt * ((wt-1) * pmut + mt) / (wt+mt-1)", ["wt"], ["mt"], "wt_to_mt"
    )
    model.add_reaction(
        "delta*mt * wt * (1 - pmut) / (wt+mt-1)", ["mt"], ["wt"], "mt_to_wt"
    )
    model.add_species({"wt": nMolecules, "mt": 0})
    model.add_parameter({"delta": degRate, "pmut": pmut})
    checkReport = crn.check(model)  # validate before simulating
    assert checkReport["valid"], checkReport["diagnostics"]
    return model


def buildCloneModel(pmut) -> crn.Model:
    """Fig 1C/1D: every molecule carries the id of its clone (0 = wild type)."""
    model = crn.Model("mtdna_clones")
    model.add_parameter({"K": nMolecules, "delta": degRate, "pmut": pmut, "r": 10.0})
    model.add_species(
        "mtDNA", attributes=["clone"]
    )  # observable "mtDNA" = all molecules
    model.add_initial("mtDNA", nMolecules, clone=0)
    model.add_observable("mutants", "mtDNA", where="x.clone > 0")
    model.add_observable("types", "mtDNA", reduce="distinct")
    model.add_observable("mutTypes", "mtDNA", where="x.clone > 0", reduce="distinct")
    model.add_observable("largest", "mtDNA", where="x.clone > 0", reduce="max_count")
    refill = (
        "r * max(0, K - mtDNA) / mtDNA"  # per molecule: free places shared by survivors
    )
    model.add_reaction(degRate, ["mtDNA"], [], name="degradation")
    model.add_reaction(
        crn.mass_action(f"{refill} * (1 - pmut)"),
        ["mtDNA"],
        [crn.copy("x"), crn.copy("x")],
        name="replication",
    )
    model.add_reaction(
        crn.mass_action(f"{refill} * pmut"),
        ["mtDNA"],
        [crn.copy("x"), crn.new("mtDNA", clone="new_id()")],
        name="mutation",
    )
    checkReport = crn.check(model)
    assert checkReport["valid"], checkReport["diagnostics"]
    return model


def buildCloneModelJSON(pmut) -> dict:
    """The same clone model as a model document (the tutorial's 0.7.1 form)."""
    x = {"label": "x", "type": "mtDNA"}
    refill = "r * max(0, K - N) / N"
    return {
        "schema": "crn-craft.model/1.0",
        "name": "mtdna_clones",
        "parameters": {"K": nMolecules, "delta": degRate, "pmut": pmut, "r": 10.0},
        "types": [{"name": "mtDNA", "attributes": ["clone"]}],
        "initial": [{"type": "mtDNA", "attributes": {"clone": 0}, "count": nMolecules}],
        "observables": {
            "N": {"type": "mtDNA"},
            "mutants": {"type": "mtDNA", "where": "x.clone > 0"},
            "types": {"type": "mtDNA", "reduce": "distinct"},
            "mutTypes": {"type": "mtDNA", "where": "x.clone > 0", "reduce": "distinct"},
            "largest": {"type": "mtDNA", "where": "x.clone > 0", "reduce": "max_count"},
        },
        "rules": [
            {"name": "degradation", "reactants": [x], "rate": "delta", "products": []},
            {
                "name": "replication",
                "reactants": [x],
                "rate": f"{refill} * (1 - pmut)",
                "products": [
                    {"type": "mtDNA", "copy": "x"},
                    {"type": "mtDNA", "copy": "x"},
                ],
            },
            {
                "name": "mutation",
                "reactants": [x],
                "rate": f"{refill} * pmut",
                "products": [
                    {"type": "mtDNA", "copy": "x"},
                    {"type": "mtDNA", "attributes": {"clone": "new_id()"}},
                ],
            },
        ],
    }


# ---------------------------------------------------------------- helpers


def _runInfo(seconds, cells, seeds, events, **extra):
    return {
        "cells": int(cells),
        "seeds": [int(s) for s in seeds],
        "events": int(events),
        "wall_seconds": round(seconds, 1),
        "threads": nThreads,
        "cpu": platform.processor() or platform.machine(),
        "crn_craft": crn.__version__,
        **extra,
    }


def _batches(nCells):
    return [(start, min(BATCH, nCells - start)) for start in range(0, nCells, BATCH)]


def _write(path, rows, header, info):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)] + [",".join(str(v) for v in row) for row in rows]
    path.with_suffix(".json").write_text(json.dumps(info, indent=1) + "\n")
    path.write_text(
        "\n".join(lines) + "\n"
    )  # last: its presence marks the unit as done


def _readCsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split(",")
    return [dict(zip(header, line.split(","), strict=True)) for line in lines[1:]]


def _fmt(value, digits=6):
    return f"{value:.{digits}g}"


# ---------------------------------------------------------------- Fig 1A


def KowaldKirkwood13Fig1A(size, out):
    """COX-negative cells (%) over 120 years for three mutation rates."""
    cfg = SIZES[size]
    endTime = cfg["fig1a_years"] * hoursPerYear
    reportTimes = np.linspace(0.0, endTime, cfg["fig1a_years"] * 12 + 1)  # monthly
    for k, pmut in enumerate(FIG1A_PMUT):
        path = out / f"fig1a_pmut{pmut:.0e}.csv"
        if path.exists():
            continue
        nCells = cfg["fig1a_cells"][pmut]
        seed = 1000 * (k + 1)
        model = buildModel(pmut)
        modelComp = crn.compile_loop(model) if cfg["fig1a_compile"] else None
        nNegative = np.zeros(len(reportTimes), dtype=np.int64)
        events, seeds, t0 = 0, [], time.perf_counter()
        for batchIdx, (_start, nBatch) in enumerate(_batches(nCells)):
            run = crn.ensemble(
                model,
                end=endTime,
                report_times=reportTimes.tolist(),
                replicates=nBatch,
                seed=seed + batchIdx,  # one seed per batch -> reproducible
                threads=nThreads,
                method="static_direct",
                keep_paths="all",  # every cell's trajectory
                strict=True,  # an incomplete cell is an error, not a silent bias
            )
            result = run.simulate(compiled=modelComp) if modelComp else run.simulate()
            mt = result.paths.sel(variable="mt").to_numpy()  # (cells x times)
            nNegative += (mt > coxThreshold * nMolecules).sum(axis=0)
            events += int(result.records["events"].sum())
            seeds.append(seed + batchIdx)
        rows = [
            (_fmt(t / hoursPerYear), int(n), _fmt(100.0 * n / nCells))
            for t, n in zip(reportTimes, nNegative, strict=True)
        ]
        info = _runInfo(
            time.perf_counter() - t0,
            nCells,
            seeds,
            events,
            pmut=pmut,
            compiled=bool(modelComp),
        )
        _write(path, rows, ["years", "cox_negative", "cox_negative_pct"], info)
        print(
            f"Fig 1A P_mut={pmut:.0e}: {rows[-1][2]} % COX-negative at {cfg['fig1a_years']} y"
        )


# ---------------------------------------------------------------- Fig 1B


def fitLogistic(x, k, n):
    """Binomial GLM logit(f) = b0 + b1*x by iteratively reweighted least squares."""
    X = np.column_stack([np.ones_like(x), x])
    b = np.zeros(2)
    for _ in range(100):
        p = 1.0 / (1.0 + np.exp(-(X @ b)))
        w = n * p * (1.0 - p)
        z = X @ b + (k - n * p) / np.maximum(w, 1e-300)
        bNew = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * z))
        if np.max(np.abs(bNew - b)) < 1e-12:
            b = bNew
            break
        b = bNew
    p = 1.0 / (1.0 + np.exp(-(X @ b)))
    cov = np.linalg.inv(X.T @ ((n * p * (1.0 - p))[:, None] * X))
    chi2 = float(np.sum((k - n * p) ** 2 / (n * p * (1.0 - p))))
    return b, cov, chi2


def invertLogistic(b, cov, target):
    """log10 P_mut for a target fraction, with a 95 % CI by the delta method."""
    logit = math.log(target / (1.0 - target))
    x = (logit - b[0]) / b[1]
    grad = np.array([-1.0 / b[1], -x / b[1]])
    se = math.sqrt(float(grad @ cov @ grad))
    return x, x - 1.96 * se, x + 1.96 * se


def _coxNegativeAtEnd(pmut, endT, nCells, seed):
    """COX-negative cells at the end of life; one seed per grid point."""
    count, events = 0, 0
    for batchIdx, (_start, nBatch) in enumerate(_batches(nCells)):
        run = crn.ensemble(
            buildModel(pmut),
            end=endT,
            report_times=[0.0, endT],
            replicates=nBatch,
            seed=seed + batchIdx,
            threads=nThreads,
            method="static_direct",
            keep_paths="all",
            strict=True,
        )
        result = run.simulate()  # no compile_loop: every grid point has a new parameter
        mtEnd = result.paths.sel(variable="mt").isel(time=-1).to_numpy()
        count += int((mtEnd > coxThreshold * nMolecules).sum())
        events += int(result.records["events"].sum())
    return count, events


def KowaldKirkwood13Fig1B(size, out):
    """Which P_mut gives 5, 10 or 15 % COX-negative cells at the end of life?"""
    cfg = SIZES[size]
    for li, years in enumerate(cfg["lifespans"]):
        path = out / f"fig1b_{years}y.csv"
        if path.exists():
            continue
        endT = years * hoursPerYear
        t0, events, rows, seeds = time.perf_counter(), 0, [], []
        # 1. pilot sweep: 17 points from 1e-6 to 1e-2
        pilot = np.linspace(-6.0, -2.0, 17)
        fractions = []
        for i, x in enumerate(pilot):
            seed = 100_000 * (li + 1) + 1000 * i
            k, ev = _coxNegativeAtEnd(10.0**x, endT, cfg["pilot_cells"], seed)
            rows.append(("pilot", _fmt(x), cfg["pilot_cells"], k))
            fractions.append(k / cfg["pilot_cells"])
            events, seeds = events + ev, seeds + [seed]
        fractions = np.array(fractions)
        # 2. fine sweep: 13 points over the 5-15 % region plus 0.2 decades
        lo = (
            pilot[np.argmax(fractions >= TARGETS[0])]
            if (fractions >= TARGETS[0]).any()
            else pilot[-1]
        )
        hi = (
            pilot[np.argmax(fractions >= TARGETS[-1])]
            if (fractions >= TARGETS[-1]).any()
            else pilot[-1]
        )
        fine = np.linspace(lo - 0.25 - 0.2, hi + 0.2, 13)
        for i, x in enumerate(fine):
            seed = 100_000 * (li + 1) + 50_000 + 1000 * i
            k, ev = _coxNegativeAtEnd(10.0**x, endT, cfg["fine_cells"], seed)
            rows.append(("fine", _fmt(x), cfg["fine_cells"], k))
            events, seeds = events + ev, seeds + [seed]
        # 3. binomial GLM on the fine sweep, inverted for each target
        xs = np.array([float(r[1]) for r in rows if r[0] == "fine"])
        ks = np.array([float(r[3]) for r in rows if r[0] == "fine"])
        ns = np.array([float(r[2]) for r in rows if r[0] == "fine"])
        b, cov, chi2 = fitLogistic(xs, ks, ns)
        estimates = {}
        for p in TARGETS:
            x, xLo, xHi = invertLogistic(b, cov, p)
            estimates[f"{round(100 * p)}"] = {
                "log10_pmut": x,
                "pmut": 10**x,
                "ci95": [10**xLo, 10**xHi],
            }
        info = _runInfo(
            time.perf_counter() - t0,
            len(pilot) * cfg["pilot_cells"] + len(fine) * cfg["fine_cells"],
            seeds,
            events,
            years=years,
            b=b.tolist(),
            cov=cov.tolist(),
            pearson_chi2=chi2,
            df=len(xs) - 2,
            estimates=estimates,
        )
        _write(path, rows, ["sweep", "log10_pmut", "cells", "cox_negative"], info)
        print(
            f"Fig 1B {years} y: "
            + ", ".join(f"{t} % -> {e['pmut']:.3g}" for t, e in estimates.items())
        )
    # summary table, the input of Fig 1C/1D
    summary = []
    for years in cfg["lifespans"]:
        info = json.loads((out / f"fig1b_{years}y.json").read_text())
        for t, e in info["estimates"].items():
            summary.append(
                (years, t, _fmt(e["pmut"]), _fmt(e["ci95"][0]), _fmt(e["ci95"][1]))
            )
    lines = ["years,target_pct,pmut,ci95_low,ci95_high"] + [
        ",".join(map(str, r)) for r in summary
    ]
    (out / "fig1b.csv").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------- Fig 1C/1D


def KowaldKirkwood13Fig1CD(size, out):
    """Number of mtDNA types and the largest clone in COX-negative cells."""
    cfg = SIZES[size]
    pmuts = {
        (int(r["years"]), int(r["target_pct"])): float(r["pmut"])
        for r in _readCsv(out / "fig1b.csv")
    }
    for li, years in enumerate(cfg["cd_lifespans"]):
        for ti, target in enumerate(TARGETS):
            pct = round(100 * target)
            path = out / f"fig1cd_{years}y_{pct}pct.csv"
            if path.exists():
                continue
            pmut, endT = pmuts[(years, pct)], years * hoursPerYear
            model = buildCloneModel(pmut)
            t0, events, rows, seeds, meanN = time.perf_counter(), 0, [], [], []
            for batchIdx, (_start, nBatch) in enumerate(_batches(cfg["cd_cells"])):
                seed = 10_000 * (li + 1) + 1000 * ti + batchIdx
                run = crn.ensemble(
                    model,
                    end=endT,
                    report_times=[0.0, endT],
                    replicates=nBatch,
                    seed=seed,
                    threads=nThreads,
                    keep_paths="all",
                    strict=True,
                )  # no method -> tree method, which creates clones at run time
                result = run.simulate()
                final = result.paths.isel(time=-1)
                total = final.sel(variable="mtDNA").to_numpy()
                mutants = final.sel(variable="mutants").to_numpy()
                negative = mutants > coxThreshold * total
                cols = [
                    final.sel(variable=v).to_numpy()[negative]
                    for v in ["types", "mutTypes", "largest"]
                ]
                for t_, m_, l_, n_ in zip(*cols, total[negative], strict=True):
                    rows.append(
                        (int(t_), int(m_), _fmt(100.0 * l_ / n_), int(t_ - m_ == 0))
                    )
                events += int(result.records["events"].sum())
                meanN.append(float(total.mean()))
                seeds.append(seed)
            info = _runInfo(
                time.perf_counter() - t0,
                cfg["cd_cells"],
                seeds,
                events,
                years=years,
                target_pct=pct,
                pmut=pmut,
                cox_negative=len(rows),
                mean_total_mtdna=float(np.mean(meanN)),
            )
            _write(
                path,
                rows,
                ["types", "mut_types", "largest_pct", "wild_type_lost"],
                info,
            )
            print(
                f"Fig 1C/1D {years} y, {pct} %: {len(rows)} COX-negative cells of {cfg['cd_cells']}"
            )


# ---------------------------------------------------------------- figures


def makeFigures(out, imgDir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    imgDir.mkdir(parents=True, exist_ok=True)
    style = {"dpi": 100, "bbox_inches": "tight"}
    # Fig 1A
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for pmut in FIG1A_PMUT:
        path = out / f"fig1a_pmut{pmut:.0e}.csv"
        if not path.exists():
            continue
        rows = _readCsv(path)
        years = np.array([float(r["years"]) for r in rows])
        pct = np.array([float(r["cox_negative_pct"]) for r in rows])
        keep = pct > 0
        ax.plot(years[keep], pct[keep], ".", markersize=2, color="black")
        ax.annotate(
            f"$P_{{mut}} = 10^{{{round(math.log10(pmut))}}}$",
            (years[-1], pct[-1]),
            xytext=(-60, 8),
            textcoords="offset points",
        )
    ax.set_yscale("log")
    ax.set_xlabel("time (years)")
    ax.set_ylabel("COX-negative cells (%)")
    ax.set_title("(A)", loc="left")
    fig.savefig(imgDir / "fig1a.png", **style)
    plt.close(fig)
    # Fig 1B fits
    lifespans = [y for y in LIFESPANS if (out / f"fig1b_{y}y.csv").exists()]
    if lifespans:
        fig, axes = plt.subplots(
            1, len(lifespans), figsize=(3.0 * len(lifespans), 3.2), squeeze=False
        )
        for ax, years in zip(axes[0], lifespans, strict=True):
            rows = _readCsv(out / f"fig1b_{years}y.csv")
            info = json.loads((out / f"fig1b_{years}y.json").read_text())
            for sweep, marker in [("pilot", "x"), ("fine", "o")]:
                sel = [r for r in rows if r["sweep"] == sweep]
                ax.plot(
                    [float(r["log10_pmut"]) for r in sel],
                    [100 * int(r["cox_negative"]) / int(r["cells"]) for r in sel],
                    marker,
                    markersize=3,
                    color="gray" if sweep == "pilot" else "black",
                    linestyle="none",
                )
            b = info["b"]
            xs = np.linspace(-6.0, -2.0, 400)
            ax.plot(xs, 100 / (1 + np.exp(-(b[0] + b[1] * xs))), color="C0")
            fine = [float(r["log10_pmut"]) for r in rows if r["sweep"] == "fine"]
            ax.set_xlim(min(fine) - 0.5, max(fine) + 0.5)
            ax.set_ylim(0, 40)
            ax.set_title(f"{years} years")
            ax.set_xlabel("log10 P_mut")
        axes[0][0].set_ylabel("COX-negative cells (%)")
        fig.tight_layout()
        fig.savefig(imgDir / "fig1b_fits.png", **style)
        plt.close(fig)
    # Fig 1C and 1D
    stats = cloneStatistics(out)
    for column, name, label, letter in [
        ("types", "fig1c.png", "Number of different mtDNAs per cell", "(C)"),
        ("largest_pct", "fig1d.png", "Most frequent mutant (% of all mtDNA)", "(D)"),
    ]:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for pct, marker, line in [(5, "o", "-"), (10, "s", ":"), (15, "v", "--")]:
            pts = [
                (y, s[column])
                for (y, p), s in sorted(stats.items())
                if p == pct and s["n"] > 0
            ]
            if pts:
                ax.errorbar(
                    [y for y, _ in pts],
                    [m for _, (m, _) in pts],
                    yerr=[e for _, (_, e) in pts],
                    marker=marker,
                    linestyle=line,
                    color="black",
                    capsize=3,
                    label=f"{pct} % COX-negative",
                )
        ax.set_xlabel("years")
        ax.set_ylabel(label)
        ax.set_xlim(0, 130)
        ax.set_ylim(0, 50 if column == "types" else 100)
        ax.set_title(letter, loc="left")
        ax.legend()
        fig.savefig(imgDir / name, **style)
        plt.close(fig)


def cloneStatistics(out):
    """Mean and standard error per (lifespan, target) over COX-negative cells."""
    stats = {}
    for years in LIFESPANS:
        for target in TARGETS:
            pct = round(100 * target)
            path = out / f"fig1cd_{years}y_{pct}pct.csv"
            if not path.exists():
                continue
            rows = _readCsv(path)
            s = {"n": len(rows)}
            for col in ["types", "mut_types", "largest_pct", "wild_type_lost"]:
                v = np.array([float(r[col]) for r in rows])
                s[col] = (
                    (float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))))
                    if len(v) > 1
                    else (math.nan, math.nan)
                )
            stats[(years, pct)] = s
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--size", choices=sorted(SIZES), default="draft")
    parser.add_argument("--part", choices=["1A", "1B", "1CD", "all"], default="all")
    parser.add_argument(
        "--figures", action="store_true", help="only plot from the result files"
    )
    parser.add_argument(
        "--results", type=Path, help="default: examples/kk13/results/<size>"
    )
    parser.add_argument(
        "--images", type=Path, default=ROOT / "docs" / "guide" / "img" / "kk13"
    )
    parser.add_argument(
        "--threads", type=int, help="default: number of CPUs minus 2, at least 1"
    )
    args = parser.parse_args()
    global nThreads
    if args.threads:
        nThreads = args.threads  # results do not depend on the number of threads
    out = args.results or HERE / "results" / args.size
    out.mkdir(parents=True, exist_ok=True)
    if args.figures:
        makeFigures(out, args.images)
        return
    if args.part in ("1A", "all"):
        KowaldKirkwood13Fig1A(args.size, out)
    if args.part in ("1B", "all"):
        KowaldKirkwood13Fig1B(args.size, out)
    if args.part in ("1CD", "all"):
        KowaldKirkwood13Fig1CD(args.size, out)


if __name__ == "__main__":
    main()
