"""Calibrate mutant fractions and score independently simulated cells."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from reference import HERE, MECHANISMS, SAMPLER, STATUS, Config, simulate

AGES = (3, 10, 40, 80)
YEAR = 360  # days per year of the tissue study
TARGET = 0.10
BAND = (0.095, 0.105)
BASE_SEED = 20261120
BATCH = 256
BOOTSTRAPS = 1999
# Starting brackets only. The comparison uses freshly fitted probabilities.
GUESSES = {
    "drift": (0.0076, 0.00145, 0.00017, 0.00007),
    "size": (0.006, 0.0006, 0.00003, 0.00001),
    "transcription": (5e-6, 1.3e-6, 3.2e-7, 1.6e-7),
}
TIERS = {
    "quick": {"pilot": 96, "refine": 128, "validate": 128, "historical": 0},
    "draft": {"pilot": 240, "refine": 800, "validate": 1000, "historical": 1000},
    "full": {
        "pilot": 400,
        "refine": 1600,
        "local": 1600,
        "confirm": 4800,
        "validate": 20000,
        "historical": 1000,
    },
}
STORED = HERE / "results" / "reconstruction"


def _seed(value):
    data = json.dumps(value, sort_keys=True, allow_nan=False).encode()
    return int.from_bytes(hashlib.blake2b(data, digest_size=8).digest(), "little")


def _code_hash():
    return hashlib.sha256((HERE / "reference.cpp").read_bytes()).hexdigest()


class Samples:
    def __init__(self, path, threads, resume, seed, executor):
        self.path = path
        self.threads = threads
        self.resume = resume
        self.seed = seed
        self.executor = executor
        self.code = _code_hash()
        self.path.mkdir(parents=True, exist_ok=True)
        self.wall = []
        self.identity = {}

    def draw(self, cfg, cells, phase, *, labelled):
        value = [self.code, SAMPLER, self.seed, asdict(cfg), phase, labelled]
        identity = _seed(value)
        self.identity[f"{cfg.mechanism}/{cfg.hours}/{phase}/{identity:016x}"] = {
            "config": asdict(cfg),
            "phase": phase,
            "cells": cells,
            "labelled": labelled,
            "seed_base": f"{identity:016x}",
        }
        rows = []
        start = time.perf_counter()
        batches = 0
        reused = 0
        for offset in range(0, cells, BATCH):
            end = min(offset + BATCH, cells)
            path = self.path / f"{identity:016x}-{offset:06d}-{end:06d}.npy"
            if self.resume and path.is_file():
                batch = np.load(path, allow_pickle=False)
                fields = (
                    "wild",
                    "mutant",
                    "clones",
                    "largest",
                    "status",
                    "time",
                    "steps",
                )
                if batch.dtype.names != fields or batch.shape != (end - offset,):
                    raise ValueError(f"malformed checkpoint {path.name}")
                reused += 1
            else:
                seeds = np.arange(offset, end, dtype=np.uint64) + np.uint64(identity)
                batch = simulate(
                    cfg,
                    seeds,
                    labelled=labelled,
                    threads=self.threads,
                    executor=self.executor,
                )
                with tempfile.NamedTemporaryFile(dir=self.path, delete=False) as tmp:
                    np.save(tmp, batch, allow_pickle=False)
                Path(tmp.name).replace(path)
            rows.append(batch)
            batches += 1
        self.wall.append(
            {
                "case": cfg.mechanism,
                "hours": cfg.hours,
                "phase": phase,
                "cells": cells,
                "batches": batches,
                "reused": reused,
                "seconds": time.perf_counter() - start,
            }
        )
        return np.concatenate(rows)


def endpoint_masks(rows):
    total = rows["wild"] + rows["mutant"]
    known = (total > 0) & (rows["status"] <= 1)
    lower = 5 * rows["mutant"] - 3 * total
    return known, known & (lower >= 0), known & (lower > 0), known & (lower == 0)


def wilson(positive, n):
    if not n:
        return None, None
    z = 1.959963984540054
    p = positive / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    radius = z / denom * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return max(0.0, center - radius), min(1.0, center + radius)


def _diversity(rows, seed):
    if not len(rows) or np.any(rows["clones"] < 0):
        return {}
    values = np.column_stack(
        [
            rows["clones"],
            rows["clones"] + (rows["wild"] > 0),
            rows["largest"] / (rows["wild"] + rows["mutant"]),
        ]
    )
    rng = np.random.default_rng(seed)
    boots = np.empty((BOOTSTRAPS, len(values[0])))
    for start in range(0, BOOTSTRAPS, 128):
        block = boots[start : start + 128]
        indices = rng.integers(0, len(values), (len(block), len(values)))
        block[:] = values[indices].mean(axis=1)
    low, high = np.quantile(boots, (0.025, 0.975), axis=0)
    out = {}
    for j, name in enumerate(("mutant_types", "all_types", "largest_share")):
        out.update(
            {
                name: float(values[:, j].mean()),
                name + "_low": float(low[j]),
                name + "_high": float(high[j]),
            }
        )
    return out


def summarize(rows, seed):
    known, affected, strict, ties = endpoint_masks(rows)
    n = len(rows)
    k = int(affected.sum())
    completed = int((rows["status"] < 3).sum())
    classified = int(known.sum())
    lo, hi = wilson(k, classified)
    band = "unresolved"
    if classified == n and lo is not None:
        band = "compatible"
        if hi < BAND[0] or lo > BAND[1]:
            band = "outside"
        elif lo >= BAND[0] and hi <= BAND[1]:
            band = "within"
    survivors = rows["status"] == 0
    out = {
        "assigned": n,
        "classified": classified,
        "cox_negative": k,
        "strict_negative": int(strict.sum()),
        "ties": int(ties.sum()),
        "collapsed": int((rows["status"] == 1).sum()),
        "empty": int((rows["status"] == 2).sum()),
        "failed": n - completed,
        "prevalence": k / classified if classified else None,
        "wilson_low": lo,
        "wilson_high": hi,
        "target_band": band,
        "prevalence_floor": k / completed if completed else None,
        "prevalence_ceiling": (k + completed - classified) / completed
        if completed
        else None,
        "survivors": int((survivors & known).sum()),
        "survivor_negative": int((survivors & affected).sum()),
        "steps": int(rows["steps"].sum()),
        **_diversity(rows[affected], seed),
    }
    if not survivors.all():
        for key, value in _diversity(rows[affected & survivors], seed ^ 2).items():
            out["survivor_" + key] = value
    if strict.sum() != k:
        for key, value in _diversity(rows[strict], seed ^ 1).items():
            out["strict_" + key] = value
    elif "mutant_types" in out:
        for key in ("mutant_types", "all_types"):
            out["strict_" + key] = out[key]
    return out


def _fit_logistic(points):
    # Reuse the estimator and inverse from the chapter 16 study.
    if str(HERE.parent / "kk13") not in sys.path:
        sys.path.insert(0, str(HERE.parent / "kk13"))
    from kowald_kirkwood_2013 import fitLogistic, invertLogistic

    p, k, n = np.array(
        [(r["mutation"], r["cox_negative"], r["assigned"]) for r in points]
    ).T
    x = np.log10(p)
    offset = float(x.mean())
    # Half a count at separation endpoints keeps the local logistic fit finite.
    b, cov, _ = fitLogistic(x - offset, k + 0.5, n + 1.0)
    if b[1] <= 0 or not np.isfinite(b).all():
        raise ValueError("calibration has no positive local logistic slope")
    return tuple(float(10.0 ** (v + offset)) for v in invertLogistic(b, cov, TARGET))


def calibrate(cfg, samples, tier, guess, label):
    points = []

    def observe(p, cells, phase):
        data = samples.draw(
            replace(cfg, mutation=p), cells, f"{label}-{phase}", labelled=False
        )
        known, affected, _, _ = endpoint_masks(data)
        point = {
            "case": label,
            "mechanism": cfg.mechanism,
            "hours": cfg.hours,
            "phase": phase,
            "mutation": p,
            "cox_negative": int(affected.sum()),
            "classified": int(known.sum()),
            "assigned": cells,
            "failed": int((data["status"] >= len(STATUS) - 1).sum()),
        }
        points.append(point)
        return point

    low, high = guess / 1.7, min(0.9999, guess * 1.7)
    low_point = observe(low, tier["pilot"], "bracket-low")
    high_point = observe(high, tier["pilot"], "bracket-high")

    def prevalence(row):
        return row["cox_negative"] / row["assigned"]

    for i in range(7):
        if prevalence(low_point) <= TARGET <= prevalence(high_point):
            break
        if prevalence(low_point) > TARGET:
            high_point = low_point
            low_point = observe(low_point["mutation"] / 2, tier["pilot"], f"lower-{i}")
        else:
            low_point = high_point
            p = min(0.9999, high_point["mutation"] * 2)
            if p == high_point["mutation"]:
                break
            high_point = observe(p, tier["pilot"], f"raise-{i}")
    a, b = low_point["mutation"], high_point["mutation"]
    if not prevalence(low_point) <= TARGET <= prevalence(high_point) or not a < b:
        return (None, None, None), points, "unbracketed"
    coarse = _fit_logistic([low_point, high_point])[0]
    refined = observe(float(np.clip(coarse, a, b)), tier["refine"], "refine")
    used = [low_point, high_point, refined]
    if any(p["classified"] != p["assigned"] for p in points):
        return (None, None, None), points, "unclassifiable_fit"
    fitted, lo, hi = _fit_logistic(used)
    if tier.get("local"):
        center = float(np.clip(fitted, a, b))
        vicinity = [
            observe(center / 1.20, tier["local"], "local-low"),
            observe(min(0.9999, center * 1.20), tier["local"], "local-high"),
        ]
        if prevalence(vicinity[0]) < prevalence(vicinity[1]):
            local_fit, _, _ = _fit_logistic(vicinity)
            center = float(np.clip(local_fit, a, b))
        confirmation = observe(center, tier["confirm"], "confirm")
        if any(p["classified"] != p["assigned"] for p in points):
            return (None, None, None), points, "unclassifiable_fit"
        if prevalence(vicinity[0]) < prevalence(vicinity[1]):
            fitted, lo, hi = _fit_logistic([*vicinity, confirmation])
        else:
            fitted, lo, hi = _fit_logistic([*used, confirmation])
    return (float(np.clip(fitted, a, b)), lo, hi), points, "fitted"


def save_csv(path, rows):
    cols = list(dict.fromkeys(k for row in rows for k in row))
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", dir=path.parent, delete=False
    ) as file:
        writer = csv.DictWriter(file, fieldnames=cols, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    Path(file.name).replace(path)


def _rates(path):
    if not path or not path.is_file():
        return {}
    with path.open(newline="") as f:
        return {
            (r["mechanism"], int(r["years"])): float(r["mutation"])
            for r in csv.DictReader(f)
            if r["label"] == "primary" and r["mutation"]
        }


def _conditions(ages, mechanisms, sensitivity):
    # The tissue study counts a year as 360 days. The 365-day clock is a
    # sensitivity at the 360-day rates.
    for years in ages:
        for name in mechanisms:
            yield "primary", years, YEAR, Config(name, years * YEAR * 24, 0)
    if sensitivity and "size" in mechanisms and 3 in ages:
        yield (
            "size-half-life-2d",
            3,
            YEAR,
            Config("size", 3 * YEAR * 24, 0, half_life=48),
        )
    if sensitivity:
        for years in (3, 80):
            if years in ages:
                for name in mechanisms:
                    yield "clock-365d", years, 365, Config(name, years * 365 * 24, 0)
    if sensitivity and 3 in ages:
        hours = 3 * YEAR * 24
        if "size" in mechanisms:
            yield "size-init-0.1", 3, YEAR, Config("size", hours, 0, initiation=0.1)
        if "transcription" in mechanisms:
            yield (
                "atp-wild-only",
                3,
                YEAR,
                Config("transcription", hours, 0, wild_only=True),
            )
            yield "atp-half-step", 3, YEAR, Config("transcription", hours, 0, step=0.5)


def _selection(rows):
    data = {
        r["mechanism"]: r
        for r in rows
        if r["label"] == "primary" and r["years"] == 3 and r["cohort"] == "validation"
    }
    if len(data) != 3 or any("mutant_types" not in r for r in data.values()):
        return "insufficient_data"
    if any(r["target_band"] not in ("compatible", "within") for r in data.values()):
        return "prevalence_gate_failed"
    upper = data["transcription"]["mutant_types_high"]
    return (
        "separated"
        if all(upper < data[k]["mutant_types_low"] for k in ("drift", "size"))
        else "overlap"
    )


def compare(
    out,
    size,
    ages,
    threads,
    resume,
    seed,
    sensitivity,
    from_rates=None,
    validation_cells=None,
    mechanisms=MECHANISMS,
):
    out.mkdir(parents=True, exist_ok=True)
    tier = dict(TIERS[size])
    if validation_cells is not None:
        tier["validate"] = validation_cells
    source_path = from_rates or (STORED / "rates.csv" if size == "quick" else None)
    frozen = _rates(source_path)
    rates, fit_points, results, hist = [], [], [], []
    chosen = {}
    tasks = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        samples = Samples(out / "batches", threads, resume, seed, executor)
        for label, years, clock, cfg in _conditions(ages, mechanisms, sensitivity):
            base = {
                "label": label,
                "years": years,
                "days_per_year": clock,
                "mechanism": cfg.mechanism,
                "half_life_h": cfg.half_life,
                "step_h": cfg.step,
                "wild_only": cfg.wild_only,
                "initiation_h": cfg.initiation,
            }
            key = (cfg.mechanism, years)
            guess = GUESSES[cfg.mechanism][AGES.index(years)]
            if label == "primary" and key in frozen:
                fit, status = (frozen[key], None, None), "imported_frozen"
            elif label not in ("primary", "size-half-life-2d"):
                fit, status = (chosen[key], None, None), "primary_rate"
            elif size == "quick" and label == "primary":
                fit, status = (guess, None, None), "smoke_guess"
            else:
                fit, points, status = calibrate(
                    cfg,
                    samples,
                    tier,
                    0.0003 if label == "size-half-life-2d" else guess,
                    label,
                )
                fit_points.extend(points)
            rate, lo, hi = fit
            rates.append(
                {
                    **base,
                    "mutation": rate,
                    "fit_low": lo,
                    "fit_high": hi,
                    "fit_status": status,
                }
            )
            if label == "primary":
                chosen[key] = rate
            if rate is None:
                print(f"{label} {cfg.mechanism} {years}y: {status}", flush=True)
                continue
            tasks.append((base, replace(cfg, mutation=rate), status))
            print(
                f"{label} {cfg.mechanism} {years}y: {status}, p={rate:.7g}", flush=True
            )
        # Persist the complete fit table before opening any independent cohort.
        save_csv(out / "rates.csv", rates)
        save_csv(out / "fit.csv", fit_points)
        for base, tuned, status in tasks:
            label, years = base["label"], base["years"]
            rate = tuned.mutation
            work = [("validation", tier["validate"], True)]
            if label == "primary" and tier["historical"]:
                work.insert(0, ("historical", tier["historical"], False))
            if label not in ("primary", "size-half-life-2d"):
                work = [("validation", min(1000, tier["validate"]), True)]
            for cohort, cells, labelled in work:
                phase = f"{label}-{cohort}"
                batch = samples.draw(tuned, cells, phase, labelled=labelled)
                result = {
                    **base,
                    "cohort": cohort,
                    "mutation": rate,
                    "fit_status": status,
                    **summarize(batch, _seed([seed, base, phase])),
                }
                results.append(result)
                if years == 3 and labelled:
                    _, affected, _, _ = endpoint_masks(batch)
                    counts = Counter(batch["clones"][affected].tolist())
                    hist.extend(
                        {**base, "types": n, "cells": k}
                        for n, k in sorted(counts.items())
                    )
                print(
                    f"{label} {tuned.mechanism} {years}y {cohort}: "
                    f"p={rate:.7g} {result['cox_negative']}/{cells} COX-negative, "
                    f"mutant mean={result.get('mutant_types', float('nan')):.3f}, "
                    f"{result['target_band']}",
                    flush=True,
                )
                save_csv(out / "comparison.csv", results)
    save_csv(out / "diversity-3y.csv", hist)
    manifest = {
        "sampler": SAMPLER,
        "reference_cpp_sha256": samples.code,
        "seed": seed,
        "crn_craft_baseline": "42b27255d228a47a7e2dbb9e882669072276a85c",
        "tier": size,
        "planned_cells": tier,
        "ages": ages,
        "mechanisms": mechanisms,
        "threshold": "5*M >= 3*(W+M)",
        "target": TARGET,
        "target_band": BAND,
        "bootstrap_resamples": BOOTSTRAPS,
        "reference_csv_sha256": hashlib.sha256(
            (HERE / "reference.csv").read_bytes()
        ).hexdigest(),
        "batch_size": BATCH,
        "three_year_gate": _selection(results),
        "batches": samples.identity,
        "rate_import": source_path.name if source_path and frozen else None,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "timings.json").write_text(
        json.dumps({"threads": threads, "draws": samples.wall}, indent=2) + "\n"
    )
    return manifest
