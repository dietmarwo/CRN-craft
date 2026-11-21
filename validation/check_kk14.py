#!/usr/bin/env python3
"""Check the chapter 17 samplers, deterministic targets and stored evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples" / "kk14"))

from campaign import STORED, _selection, endpoint_masks, summarize  # noqa: E402
from native import _end, atp_comparison, drift_pair, marked_hybrid  # noqa: E402
from reference import SAMPLER, Config, Endpoint, simulate, size_capacity  # noqa: E402
from render import comparison_table, read_csv  # noqa: E402
from study import single_cell  # noqa: E402

# Each empirical-CDF check uses its own 0.001 bound. For a comparison of two
# samples, apply the DKW bound to each sample with half the error allowance.
ALPHA = 0.001
N_LAW = 16000


def require(ok, detail):
    if not ok:
        raise AssertionError(detail)


def cdf_bound(n, alpha):
    return math.sqrt(math.log(2 / alpha) / (2 * n))


def empirical_cdf(sample, last):
    return np.bincount(sample.astype(int), minlength=last + 1).cumsum() / len(sample)


def paired_cdf(a, b):
    last = max(int(a.max()), int(b.max()))
    return float(np.max(np.abs(empirical_cdf(a, last) - empirical_cdf(b, last))))


def endpoints():
    data = np.zeros(5, dtype=np.dtype(Endpoint))
    data["wild"] = [400, 399, 700, 0, 0]
    data["mutant"] = [600, 601, 300, 0, 1001]
    data["clones"] = [30, 1, 3, 0, 1]
    data["largest"] = [20, 601, 100, 0, 1001]
    data["status"] = [0, 1, 0, 2, 3]
    known, affected, strict, ties = endpoint_masks(data)
    require(list(known) == [True, True, True, False, False], "known endpoints")
    require(list(affected) == [True, True, False, False, False], "inclusive 60%")
    require(list(strict) == [False, True, False, False, False], "strict 60%")
    require(int(ties.sum()) == 1, "exact tie")
    row = summarize(data, 106)
    require(row["mutant_types"] == 15.5 and row["all_types"] == 16.5, "count wild once")
    require(
        row["survivor_mutant_types"] == 30.0 and row["collapsed"] == 1, "ATP endpoint"
    )
    require(
        row["prevalence_floor"] == 0.5 and row["prevalence_ceiling"] == 0.75,
        "account for empty endpoints without counting computational failures",
    )
    return {
        "known": row["classified"],
        "strict": row["strict_negative"],
        "ties": row["ties"],
    }


def deterministic():
    # The size paper prints these times in years. They match years of 365
    # days (19.33 years for the 10-h case; 19.60 with 360 days), so the
    # printed values are compared in hours through that conversion.
    year = 365 * 24
    with tempfile.TemporaryDirectory() as tmp:
        row = single_cell(Path(tmp))
        hours = row["size_checks_hours"]
        require(abs(row["atp_collapse_days"] - 247.4061) < 0.1, "PNAS Figure 2 zero")
        require(
            abs(row["size_extinction_hours"] - 94.9524 * year) < 0.05 * year,
            "JTB Figure 2",
        )
        require(abs(row["size_wild_busy_equilibrium"] - 0.0058) < 4e-5, "WT busy pool")
        require(
            abs(row["size_final_mutant_busy_fraction"] - 0.0029) < 3e-5,
            "mutant busy pool",
        )
        for mesh in ("mesh_0.5_h", "mesh_0.125_h"):
            require(
                abs(hours[mesh] - row["size_extinction_hours"]) < 1e-4 * year,
                "delay step convergence",
            )
        require(abs(hours["half_life_2_days"] - 3.87 * year) < 0.08 * year, "JTB Fig 3")
        require(
            abs(hours["wild_replication_10_h"] - 19.33 * year) < 0.2 * year,
            "JTB Figure 4",
        )
        require(abs(hours["mutant_size_0.7"] - 160 * year) < 3 * year, "JTB Figure 5")
        for days in (365, 360):
            require(
                row[f"size_extinction_years_{days}d"]
                == row["size_extinction_hours"] / (days * 24),
                "year conversions",
            )
        require(
            abs(row["same_size_final_ratio"] - 999) < 1e-6, "equal-delay neutrality"
        )
        return {
            "atp_zero_days": round(row["atp_collapse_days"], 4),
            "wild_displacement_hours": round(row["size_extinction_hours"], 1),
        }


def drift_law():
    n, m0, hours = 8, 3, 45
    half_life, p = 24, 0.13
    d = math.log(2) / half_life
    rates = np.zeros((n + 1, n + 1))
    for m in range(n):
        w = n - m
        rates[m, m + 1] = d * w * ((w - 1) * p + m) / (n - 1)
        if m:
            rates[m, m - 1] = d * m * w * (1 - p) / (n - 1)
        rates[m, m] = -rates[m].sum()
    target = expm(hours * rates)[m0].cumsum()
    cfg = Config(
        "drift", hours, p, population=n, initial_mutants=m0, half_life=half_life
    )
    bound = cdf_bound(N_LAW, ALPHA)
    drift = {}
    for i, labelled in enumerate((False, True)):
        rows = simulate(
            cfg, np.arange(i * N_LAW, (i + 1) * N_LAW), labelled=labelled, threads=2
        )
        require((rows["wild"] + rows["mutant"] == n).all(), "fixed-N drift")
        gap = float(np.max(np.abs(empirical_cdf(rows["mutant"], n) - target)))
        require(gap < bound, "exact drift count law")
        drift["labelled" if labelled else "aggregate"] = round(gap, 6)
    model = drift_pair(p, n, m0, half_life)
    native = _end(model, "tree", hours, N_LAW, 2, 3298)
    finals = native.paths.isel(time=-1)
    require((native.paths.sel(variable="mtDNA") == n).all(), "native pair conservation")
    gap = float(
        np.max(
            np.abs(empirical_cdf(finals.sel(variable="mutants").to_numpy(), n) - target)
        )
    )
    require(gap < bound, "native pair law")
    return {**drift, "native_pair": round(gap, 6), "cdf_bound": round(bound, 6)}


def literal_queue(cfg, rng):
    """One Bernoulli trial per template. Independent of the batched backend."""
    delta = math.log(2) / cfg.half_life
    c = size_capacity(cfg.population, cfg.half_life, cfg.wild_delay, cfg.initiation)
    free = [0] * (cfg.population - cfg.initial_mutants) + [1] * cfg.initial_mutants
    busy = []
    next_id = int(bool(cfg.initial_mutants)) + 1
    for hour in range(cfg.hours):
        future = []
        for due, clone in busy:
            if due == hour:
                free.append(clone)
                if rng.random() < cfg.mutation:
                    free.append(next_id)
                    next_id += 1
                else:
                    free.append(clone)
            elif rng.random() >= delta:
                future.append((due, clone))
        free = [clone for clone in free if rng.random() >= delta]
        p_start = cfg.initiation * c / (c + len(free) + len(future))
        ready = []
        for clone in free:
            if rng.random() < p_start:
                tau = cfg.mutant_delay if clone else cfg.wild_delay
                future.append((hour + tau, clone))
            else:
                ready.append(clone)
        free, busy = ready, future
    ids = [*free, *(c for _, c in busy)]
    wild = ids.count(0)
    return wild, len(ids) - wild, len(set(ids) - {0})


def size_law():
    cfg = Config("size", 6, 0.23, population=8, initial_mutants=2, half_life=8)
    rng = np.random.default_rng(3101)
    literal = np.asarray([literal_queue(cfg, rng) for _ in range(N_LAW)])
    seeds = np.arange(30100, 30100 + N_LAW)
    batch = simulate(cfg, seeds, labelled=True, threads=2)
    coarse = simulate(cfg, seeds, labelled=False, threads=2)
    bound = 2 * cdf_bound(N_LAW, ALPHA / 2)
    gaps = {}
    for i, key in enumerate(("wild", "mutant", "clones")):
        gap = paired_cdf(batch[key], literal[:, i])
        require(gap < bound, f"hourly queue {key} distribution")
        gaps[key] = round(gap, 6)
    for key in ("wild", "mutant", "status", "time", "steps"):
        require(
            np.array_equal(batch[key], coarse[key]), f"hourly queue marks changed {key}"
        )
    # A one-hour run starts replication but has no completions. Deaths
    # are the only change in genotype counts even though templates move.
    unit = Config("size", 1, 1, population=20, half_life=8)
    first = simulate(unit, seeds, threads=2)
    require((first["mutant"] == 0).all(), "no same-hour replication completion")
    return {**gaps, "cdf_bound": round(bound, 6)}


def atp_and_threads():
    configs = [
        Config("drift", 60, 0.04, population=24, initial_mutants=3),
        Config("size", 60, 0.04, population=24, initial_mutants=3, half_life=8),
        Config(
            "transcription", 60, 0.04, population=24, initial_mutants=3, half_life=20
        ),
    ]
    for cfg in configs:
        seeds = np.arange(101, 164)
        a, b = simulate(cfg, seeds), simulate(cfg, seeds, threads=3)
        require(np.array_equal(a, b), f"{cfg.mechanism} thread replay")
        if cfg.mechanism != "drift":
            c = simulate(cfg, seeds, labelled=False, threads=2)
            for key in ("wild", "mutant", "status", "time", "steps"):
                require(np.array_equal(a[key], c[key]), f"{cfg.mechanism} marked {key}")
    rows = atp_comparison(reps=320)
    a, b = (r["cox_negative"] for r in rows[:2])
    require(abs(a - b) / 320 < 2 * cdf_bound(320, ALPHA / 2), "hybrid vs thinning law")
    marks = marked_hybrid()
    require(marks["matched_counts"], "marked replay changes the native run")
    require(
        any(c > 1 for c in marks["clone_counts"]), "marked run needs multiple lineages"
    )
    return {
        "thread_replay": len(configs),
        "native_events_marked": marks["events"],
        "native_cox_negative": [r["cox_negative"] for r in rows],
    }


def stored():
    rates = read_csv(STORED / "rates.csv")
    results = read_csv(STORED / "comparison.csv")
    refs = read_csv(STORED.parents[1] / "reference.csv")
    deltas = read_csv(STORED / "reference-differences.csv")
    manifest = json.loads((STORED / "manifest.json").read_text())
    code = hashlib.sha256(
        (STORED.parents[1] / "reference.cpp").read_bytes()
    ).hexdigest()
    require(
        manifest["sampler"] == SAMPLER and manifest["reference_cpp_sha256"] == code,
        "stored sampler revision",
    )
    keys = [(r["mechanism"], r["years"]) for r in rates if r["label"] == "primary"]
    require(len(keys) == len(set(keys)) == 12, "all mechanism-age fits")
    # The tissue study counts 360 days per year; 365 days is a sensitivity.
    # Every primary rate is refitted, none imported.
    require(
        manifest["rate_import"] is None
        and all(
            r["days_per_year"] == "360" and r["fit_status"] == "fitted"
            for r in rates
            if r["label"] in ("primary", "size-half-life-2d")
        ),
        "primary fits on the 360-day year",
    )
    clock = [r for r in results if r["label"] == "clock-365d"]
    require(
        len(clock) == 6 and all(r["days_per_year"] == "365" for r in clock),
        "365-day sensitivity",
    )
    validation = [
        r for r in results if r["label"] == "primary" and r["cohort"] == "validation"
    ]
    require(
        len(validation) == 12 and all(int(r["assigned"]) == 20000 for r in validation),
        "full saved validation cohorts",
    )
    require(
        len(refs) == len(deltas) == 42,
        "published reference rows and numeric differences",
    )
    typed = [
        {
            **r,
            "years": int(r["years"]),
            "mutant_types_low": float(r["mutant_types_low"]),
            "mutant_types_high": float(r["mutant_types_high"]),
        }
        for r in validation
    ]
    require(
        _selection(typed) == manifest["three_year_gate"],
        "three-year science gate replay",
    )
    guide = (ROOT / "docs" / "guide" / "17b-mtdna-mechanisms-results.md").read_text()
    part = guide.partition("<!-- kk14-comparison -->\n")[2].partition(
        "<!-- /kk14-comparison -->"
    )[0]
    require(part == comparison_table(STORED), "rendered guide table is out of date")
    # The deterministic size results are hours; the page gives both year
    # conversions next to them.
    cell = json.loads((STORED / "single-cell.json").read_text())
    models = (ROOT / "docs" / "guide" / "17a-mtdna-mechanisms-models.md").read_text()
    hours = {"": cell["size_extinction_hours"], **cell["size_checks_hours"]}
    for name in ("", "half_life_2_days", "wild_replication_10_h", "mutant_size_0.7"):
        text = (
            f"{hours[name]:,.0f} h ({hours[name] / 24:,.0f} days) | "
            f"{hours[name] / 8760:.3f} / {hours[name] / 8640:.3f} |"
        )
        require(text in models, f"deterministic size row {name or 'base'} on 17a")
    return {
        "fits": len(keys),
        "reference_rows": len(refs),
        "cohorts": len(validation),
        "three_year_gate": manifest["three_year_gate"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samplers-only", action="store_true")
    args = parser.parse_args()
    report = {"schema": "crn-craft.validation.kk14/0.1", "checks": {}, "failures": []}
    for check in (
        endpoints,
        deterministic,
        drift_law,
        size_law,
        atp_and_threads,
        stored,
    ):
        if args.samplers_only and check is stored:
            continue
        try:
            report["checks"][check.__name__] = check()
        except Exception as exc:
            report["failures"].append(f"{check.__name__}: {type(exc).__name__}: {exc}")
    report["passed"] = not report["failures"]
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    for failure in report["failures"]:
        print(failure, file=sys.stderr)
    print(f"kk14: {len(report['checks'])} checks, {len(report['failures'])} failures")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
