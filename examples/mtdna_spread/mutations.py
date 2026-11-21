"""De novo deletions in the survival-of-the-densest chain, with mutant types.

An extension made in guide chapter 19, not a result of the paper. Every
region starts with wild-type copies only. A replication makes a new
deletion type with probability ``m``, from wild-type and mutant templates
alike; all mutant types share the density parameter ``delta``. The exact
sampler below follows every copy: in region ``i`` every copy, of any type,
is copied at ``lambda_i``, dies at ``mu`` and hops to each neighbor at
``gamma``, so an event picks a region by its total rate and then a copy
uniformly. Labels are kept with mutants first, so the run with labels and
the run without them draw the same random numbers and give the same counts.
"""

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
MANIFEST = json.loads((ROOT / "manifest.json").read_text())
DAYS_PER_YEAR = 360


def _kernel(regions, nss, c, mu, gamma, delta, m, end, seed, labels, capacity, record):
    """Exact path of one fiber; returns counts, labels and event totals.

    ``labels`` keeps one type label per copy (0 = wild type); without it only
    the counts are kept. ``record`` > 0 stores the counts every ``record``
    days. Both choices leave the random stream and the counts unchanged.
    """
    np.random.seed(seed)
    wild = np.full(regions, nss, dtype=np.int64)
    mutant = np.zeros(regions, dtype=np.int64)
    width = capacity if labels else 1
    tag = np.zeros((regions, width), dtype=np.int64)
    rate = np.zeros(regions)
    neighbors = np.full(regions, 2.0)
    if regions == 1:
        neighbors[0] = 0.0
    else:
        neighbors[0] = 1.0
        neighbors[regions - 1] = 1.0
    for i in range(regions):
        n = wild[i] + mutant[i]
        birth = max(0.0, mu + c * (nss - wild[i] - delta * mutant[i]))
        rate[i] = n * (birth + mu + gamma * neighbors[i])
    slots = int(end // record) + 1 if record > 0 else 1
    history = np.zeros((slots, 2, regions), dtype=np.int64)
    events = np.zeros(4, dtype=np.int64)  # births, deaths, hops, mutations
    following = 1
    total = rate.sum()
    t = 0.0
    step = 0
    while total > 0:
        t += -math.log(1.0 - np.random.random()) / total
        if record > 0:
            while following < slots and following * record <= t and t <= end:
                history[following, 0] = wild
                history[following, 1] = mutant
                following += 1
        if t > end:
            break
        target = np.random.random() * total
        i = 0
        while i < regions - 1 and target >= rate[i]:
            target -= rate[i]
            i += 1
        while rate[i] <= 0:
            i -= 1
        n = wild[i] + mutant[i]
        k = int(np.random.random() * n)
        if k >= n:
            k = n - 1
        birth = max(0.0, mu + c * (nss - wild[i] - delta * mutant[i]))
        u = np.random.random() * (birth + mu + gamma * neighbors[i])
        is_mutant = k < mutant[i]
        label = tag[i, k] if labels else 0
        j = -1
        if u < birth:
            events[0] += 1
            if np.random.random() < m:
                events[3] += 1
                child = events[3]
                child_mutant = True
            else:
                child = label
                child_mutant = is_mutant
            j = i
        else:
            # Remove copy k from region i, keeping mutants in [0, M).
            if labels:
                if is_mutant:
                    tag[i, k] = tag[i, mutant[i] - 1]
                    tag[i, mutant[i] - 1] = tag[i, n - 1]
                else:
                    tag[i, k] = tag[i, n - 1]
            if is_mutant:
                mutant[i] -= 1
            else:
                wild[i] -= 1
            if u < birth + mu:
                events[1] += 1
            else:
                events[2] += 1
                if i == 0:
                    j = 1
                elif i == regions - 1:
                    j = regions - 2
                elif u < birth + mu + gamma:
                    j = i - 1
                else:
                    j = i + 1
                child = label
                child_mutant = is_mutant
        if j >= 0:
            size = wild[j] + mutant[j]
            if labels:
                if size >= width:
                    raise ValueError("label capacity exceeded")
                if child_mutant:
                    tag[j, size] = tag[j, mutant[j]]
                    tag[j, mutant[j]] = child
                else:
                    tag[j, size] = 0
            if child_mutant:
                mutant[j] += 1
            else:
                wild[j] += 1
        for r in (i, j):
            if r >= 0:
                size = wild[r] + mutant[r]
                birth = max(0.0, mu + c * (nss - wild[r] - delta * mutant[r]))
                fresh = size * (birth + mu + gamma * neighbors[r])
                total += fresh - rate[r]
                rate[r] = fresh
        step += 1
        if step % 65536 == 0:
            total = rate.sum()
    if record > 0:
        while following < slots:
            history[following, 0] = wild
            history[following, 1] = mutant
            following += 1
        history[0, 0] = nss
    return wild, mutant, tag, events, history


try:  # Numba is optional: the same function, compiled.
    from numba import njit

    _fast = njit(cache=False)(_kernel)
except ImportError:  # pragma: no cover - CI has no Numba
    _fast = None


def simulate(
    regions,
    m,
    days,
    seed,
    *,
    nss=None,
    c=None,
    mu=None,
    gamma=None,
    delta=None,
    labels=True,
    record=0.0,
    compiled=True,
):
    """One fiber of ``regions`` regions over ``days`` days, all wild type at t = 0."""
    p = MANIFEST["fig2_site_parameters"]
    nss = p["copies"] if nss is None else nss
    c = p["c"] if c is None else c
    mu = p["death"] if mu is None else mu
    gamma = p["hop"] if gamma is None else gamma
    delta = p["delta"] if delta is None else delta
    capacity = int(1.6 * nss / delta) + 64
    run = _fast if compiled and _fast is not None else _kernel
    wild, mutant, tag, events, history = run(
        int(regions),
        int(nss),
        float(c),
        float(mu),
        float(gamma),
        float(delta),
        float(m),
        float(days),
        int(seed),
        bool(labels),
        capacity,
        float(record),
    )
    return {
        "wild": wild,
        "mutant": mutant,
        "labels": [tag[i, : mutant[i]].copy() for i in range(regions)]
        if labels
        else None,
        "events": events,
        "history": history,
    }


def crn_model(regions, nss, c, mu, gamma, delta, m):
    """The same reactions as a CRN-craft model with lineage attributes.

    One attributed species, ``mtDNA[clone, region]``; ``clone`` 0 is the wild
    type and ``new_id()`` founds a mutant type. Per region, observables give
    ``W``, ``M`` and the number of mutant types; births have one rule per
    region because their rate reads that region's counts. Hops change the
    ``region`` attribute; a copy at an end hops only inward.
    """
    import crn_craft as crn

    model = crn.Model("spatial_de_novo")
    model.add_parameter(
        {"Nss": nss, "c": c, "mu": mu, "gamma": gamma, "delta": delta, "m": m}
    )
    model.add_species("mtDNA", attributes=["clone", "region"])
    for i in range(regions):
        model.add_initial("mtDNA", nss, clone=0, region=i)
        at = f"x.region == {i}"
        model.add_observable(f"W{i}", "mtDNA", where=f"{at} && x.clone == 0")
        model.add_observable(f"M{i}", "mtDNA", where=f"{at} && x.clone > 0")
        model.add_observable(
            f"T{i}", "mtDNA", where=f"{at} && x.clone > 0", reduce="distinct"
        )
    model.add_observable("types", "mtDNA", where="x.clone > 0", reduce="distinct")
    model.add_reaction(mu, ["mtDNA"], [], name="death")
    for i in range(regions):
        birth = f"max(0, mu + c * (Nss - W{i} - delta * M{i}))"
        here = crn.pattern("mtDNA", where=f"x.region == {i}")
        model.add_reaction(
            crn.mass_action(f"{birth} * (1 - m)"),
            [here],
            [crn.copy("x"), crn.copy("x")],
            name=f"copy{i}",
        )
        model.add_reaction(
            crn.mass_action(f"{birth} * m"),
            [here],
            [crn.copy("x"), crn.copy("x", clone="new_id()")],
            name=f"delete{i}",
        )
    if regions > 1:
        last = regions - 1
        model.add_reaction(
            crn.mass_action("gamma"),
            [crn.pattern("mtDNA", where=f"x.region < {last}")],
            [crn.copy("x", region="x.region + 1")],
            name="hop_right",
        )
        model.add_reaction(
            crn.mass_action("gamma"),
            [crn.pattern("mtDNA", where="x.region > 0")],
            [crn.copy("x", region="x.region - 1")],
            name="hop_left",
        )
    return model


# Cases of the experiment: density parameter and hop rate per directed edge.
CASES = {
    "main": {"delta": 0.2, "gamma": 0.1},
    "neutral": {"delta": 1.0, "gamma": 0.1},
    "isolated": {"delta": 0.2, "gamma": 0.0},
}
# quick: a small model in pure Python, for CI and to read the code by running
# it; draft and full: the Figure 2 parameters over 32 regions.
SIZES = {
    "quick": {
        "model": {"nss": 20, "c": 0.01},
        "regions": 4,
        "days": 120.0,
        "fibers": (2, 2, 8),
        "boot": 199,
        "compiled": False,
    },
    "draft": {
        "model": {},
        "regions": 32,
        "days": None,
        "fibers": (2, 2, 8),
        "boot": 999,
        "compiled": True,
    },
    "full": {
        "model": {},
        "regions": 32,
        "days": None,
        "fibers": (4, 8, 32),
        "boot": 1999,
        "compiled": True,
    },
}
GRID = (0.9, 0.97, 1.04, 1.12)
LOCAL = (1 / 1.02, 1.0, 1.02)
TARGET = 0.10


def seed_for(*parts):
    """A 31-bit seed from the manifest seed and the batch's identity."""
    text = json.dumps([MANIFEST["seed"], *parts], sort_keys=True)
    return int(hashlib.sha256(text.encode()).hexdigest()[:8], 16) % 2**31


def affected(wild, mutant):
    """At least 60% mutants, in integers; an empty unit is not affected."""
    n = np.asarray(wild) + np.asarray(mutant)
    return (5 * np.asarray(mutant) >= 3 * n) & (n > 0)


def summarize(out):
    """Region, zone and segment rows of one simulated fiber."""
    wild, mutant, labels = out["wild"], out["mutant"], out["labels"]
    hit = affected(wild, mutant)
    regions = []
    for i in range(len(wild)):
        values, counts = np.unique(labels[i], return_counts=True)
        regions.append(
            {
                "region": i,
                "wild": int(wild[i]),
                "mutant": int(mutant[i]),
                "affected": bool(hit[i]),
                "types": len(values),
                "largest": int(counts.max()) if len(counts) else 0,
            }
        )
    zones, i = [], 0
    while i < len(wild):
        if not hit[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(wild) and hit[j + 1]:
            j += 1
        pooled = np.concatenate(labels[i : j + 1])
        zones.append(
            {
                "first": i,
                "last": j,
                "wild": int(wild[i : j + 1].sum()),
                "mutant": int(mutant[i : j + 1].sum()),
                "types": len(np.unique(pooled)),
            }
        )
        i = j + 1
    segment = {
        "wild": int(wild.sum()),
        "mutant": int(mutant.sum()),
        "affected": bool(affected(wild.sum(), mutant.sum())),
        "types": len(np.unique(np.concatenate(labels))),
        "mutations": int(out["events"][3]),
    }
    return {"regions": regions, "zones": zones, "segment": segment}


def fiber(case, phase, index, m, size, work, years):
    """One fiber, checkpointed as JSON so an interrupted run resumes."""
    spec = SIZES[size]
    file = work / case / f"{phase}-{m:.8e}-{index:03d}.json"
    if file.is_file():
        return json.loads(file.read_text())
    days = spec["days"] or DAYS_PER_YEAR * years
    seed = seed_for(case, size, years, phase, f"{m:.8e}", index)
    start = time.perf_counter()
    out = simulate(
        spec["regions"],
        m,
        days,
        seed,
        **spec["model"],
        **CASES[case],
        compiled=spec["compiled"],
    )
    result = summarize(out)
    result.update(
        {
            "case": case,
            "phase": phase,
            "index": index,
            "m": m,
            "seed": seed,
            "days": days,
            "events": [int(v) for v in out["events"]],
            "seconds": round(time.perf_counter() - start, 3),
        }
    )
    file.parent.mkdir(parents=True, exist_ok=True)
    pending = file.with_suffix(".tmp")
    pending.write_text(json.dumps(result), encoding="utf-8")
    pending.replace(file)
    return result


def cohort(case, phase, m, fibers, size, work, years):
    return [fiber(case, phase, k, m, size, work, years) for k in range(fibers)]


def counts(group):
    hits = sum(r["affected"] for f in group for r in f["regions"])
    return hits, sum(len(f["regions"]) for f in group)


def logistic(points):
    """Fit P(affected) = 1 / (1 + exp(-(a + b log10 m))) by Newton's method.

    ``points`` holds (m, affected, regions). Returns a, b, the m at TARGET
    and its 95% delta-method interval (regions treated as independent).
    """
    x = np.log10([p[0] for p in points])
    k = np.array([p[1] for p in points], dtype=float)
    n = np.array([p[2] for p in points], dtype=float)
    design = np.column_stack([np.ones_like(x), x])
    coef = np.zeros(2)
    p = k.sum() / n.sum()
    coef[0] = math.log(p / (1 - p)) if 0 < p < 1 else 0.0
    for _ in range(100):
        prob = 1 / (1 + np.exp(-(design @ coef)))
        weight = n * prob * (1 - prob)
        info = design.T @ (design * weight[:, None])
        stepped = np.linalg.solve(info + 1e-12 * np.eye(2), design.T @ (k - n * prob))
        coef = coef + stepped
        if np.max(np.abs(stepped)) < 1e-10:
            break
    prob = 1 / (1 + np.exp(-(design @ coef)))
    info = design.T @ (design * (n * prob * (1 - prob))[:, None])
    cov = np.linalg.inv(info)
    a, b = coef
    logit = math.log(TARGET / (1 - TARGET))
    x0 = (logit - a) / b
    grad = np.array([-1 / b, -(logit - a) / b**2])
    se = math.sqrt(grad @ cov @ grad)
    return a, b, 10**x0, 10 ** (x0 - 1.96 * se), 10 ** (x0 + 1.96 * se)


def calibrate(case, size, work, years, log, attempt=1):
    """Stages 0 to 2 of the predeclared fit; returns the frozen m and fit rows.

    Attempt 2, run only when the first validation lies outside the band,
    repeats stage 2 with new seeds around the first frozen m and fits
    through all fit points.
    """
    spec = SIZES[size]
    days = spec["days"] or DAYS_PER_YEAR * years
    mu = MANIFEST["fig2_site_parameters"]["death"]
    guess = -math.log(0.4) / (mu * days)
    probe = fiber(case, "stage0", 0, guess, size, work, years)
    h = np.array([r["mutant"] / (r["wild"] + r["mutant"]) for r in probe["regions"]])
    wanted = min(0.6 - 1.2816 * h.std(ddof=1), 0.59)
    m0 = guess * math.log(1 - wanted) / math.log(1 - h.mean())
    rows = []

    def add(stage, m, fibers):
        group = cohort(case, stage, m, fibers, size, work, years)
        hits, total = counts(group)
        rows.append((stage, m, fibers, hits, total))
        log(f"{case} {stage} m={m:.6e} affected {hits}/{total}")

    n1, n2, _ = spec["fibers"]
    for factor in GRID:
        add("stage1", m0 * factor, n1)
    for _ in range(4):  # extend the grid until it brackets the target
        fractions = [r[3] / r[4] for r in rows]
        if min(fractions) > TARGET:
            add("stage1", min(r[1] for r in rows) / 1.25, n1)
        elif max(fractions) < TARGET:
            add("stage1", max(r[1] for r in rows) * 1.25, n1)
        else:
            break
    m1 = logistic([(r[1], r[3], r[4]) for r in rows])[2]
    for factor in LOCAL:
        add("stage2", m1 * factor, n2)
    for repeat in range(2, attempt + 1):
        previous = logistic([(r[1], r[3], r[4]) for r in rows])[2]
        for factor in LOCAL:
            add(f"stage2-a{repeat}", previous * factor, n2)
    a, b, m, low, high = logistic([(r[1], r[3], r[4]) for r in rows])
    fit = [
        {
            "case": case,
            "years": years,
            "stage": s,
            "m": mm,
            "fibers": f,
            "affected": k,
            "regions": n,
        }
        for s, mm, f, k, n in rows
    ]
    rate = {
        "case": case,
        "years": years,
        "m": m,
        "m_low95": low,
        "m_high95": high,
        "intercept": a,
        "slope": b,
        "start_m": m0,
        "days": days,
        "attempt": attempt,
    }
    return rate, fit


def wilson(k, n, z=1.959963984540054):
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def band(low, high):
    if 0.095 <= low and high <= 0.105:
        return "within"
    if high < 0.095 or low > 0.105:
        return "outside"
    return "compatible"


def statement(low, high):
    if high < 2:
        return "about one type"
    if low > 10:
        return "many types"
    return "intermediate"


def ratio_mean(values, groups, resamples, rng):
    """Mean over units pooled from fibers; 95% percentile bootstrap over fibers."""
    by_fiber = [np.asarray(v, dtype=float) for v in groups]
    estimate = float(np.mean(values)) if len(values) else float("nan")
    draws = []
    for _ in range(resamples):
        pick = rng.integers(0, len(by_fiber), len(by_fiber))
        pooled = np.concatenate([by_fiber[i] for i in pick])
        if len(pooled):
            draws.append(pooled.mean())
    if not draws:
        return estimate, float("nan"), float("nan"), [], 0
    low, high = np.quantile(draws, [0.025, 0.975])
    return estimate, float(low), float(high), draws, resamples - len(draws)


def validation_tables(case, years, rate, group, resamples):
    """Per-unit rows and the summary row of one validation cohort."""
    regions, zones, segments = [], [], []
    for f in group:
        for r in f["regions"]:
            regions.append({"case": case, "years": years, "fiber": f["index"], **r})
        for z in f["zones"]:
            zones.append({"case": case, "years": years, "fiber": f["index"], **z})
        segments.append(
            {"case": case, "years": years, "fiber": f["index"], **f["segment"]}
        )
    rng = np.random.default_rng(seed_for("bootstrap", case, years, rate["attempt"]))
    hit = [[r["affected"] for r in f["regions"]] for f in group]
    k = sum(map(sum, hit))
    n = sum(map(len, hit))
    _, prev_low, prev_high, _, _ = ratio_mean(
        [v for h in hit for v in h], hit, resamples, rng
    )
    types = [[r["types"] for r in f["regions"] if r["affected"]] for f in group]
    flat = [v for t in types for v in t]
    mean, low, high, draws, empty = ratio_mean(flat, types, resamples, rng)
    largest = [
        [
            r["largest"] / (r["wild"] + r["mutant"])
            for r in f["regions"]
            if r["affected"]
        ]
        for f in group
    ]
    big, big_low, big_high, _, _ = ratio_mean(
        [v for t in largest for v in t], largest, resamples, rng
    )
    per_type = [
        [r["mutant"] / r["types"] for r in f["regions"] if r["affected"]] for f in group
    ]
    copies, copies_low, copies_high, _, _ = ratio_mean(
        [v for t in per_type for v in t], per_type, resamples, rng
    )
    zone = [[z["types"] for z in f["zones"]] for f in group]
    zmean, zlow, zhigh, _, _ = ratio_mean(
        [v for t in zone for v in t], zone, resamples, rng
    )
    seg = [[f["segment"]["types"]] for f in group]
    smean, slow, shigh, _, _ = ratio_mean(
        [v for t in seg for v in t], seg, resamples, rng
    )
    w_low, w_high = wilson(k, n)
    enough = len(flat) >= 20
    summary = {
        "case": case,
        "years": years,
        "days": rate["days"],
        "m": rate["m"],
        "fibers": len(group),
        "regions": n,
        "affected": k,
        "prevalence": k / n,
        "wilson_low": w_low,
        "wilson_high": w_high,
        "fiber_boot_low": prev_low,
        "fiber_boot_high": prev_high,
        "band": band(w_low, w_high),
        "types_mean": mean,
        "types_low": low,
        "types_high": high,
        "types_min": min(flat) if flat else 0,
        "types_max": max(flat) if flat else 0,
        "copies_per_type": copies,
        "copies_per_type_low": copies_low,
        "copies_per_type_high": copies_high,
        "largest_share": big,
        "largest_low": big_low,
        "largest_high": big_high,
        "zones": len([v for t in zone for v in t]),
        "zone_types_mean": zmean,
        "zone_types_low": zlow,
        "zone_types_high": zhigh,
        "segments_affected": sum(s["affected"] for s in segments),
        "segment_types_mean": smean,
        "segment_types_low": slow,
        "segment_types_high": shigh,
        "empty_resamples": empty,
        "statement": statement(low, high) if enough else "too few affected regions",
    }
    return summary, regions, zones, segments, draws


def run(args):
    import pandas as pd

    work = args.work or ROOT / "runs" / "de-novo" / args.size
    output = args.output or work
    output.mkdir(parents=True, exist_ok=True)
    resamples = SIZES[args.size]["boot"]

    def log(text):
        print(text, flush=True)

    rates, fits, summaries, regions, zones, segments, draws = [], [], [], [], [], [], {}
    for case in args.cases:
        for attempt in (1, 2):
            rate, fit = calibrate(case, args.size, work, args.years, log, attempt)
            rates.append(rate)
            fits.extend(
                f for f in fit if attempt == 1 or f["stage"].endswith(f"a{attempt}")
            )
            log(f"{case}: frozen m = {rate['m']:.6e} (attempt {attempt})")
            n3 = SIZES[args.size]["fibers"][2]
            phase = "validation" if attempt == 1 else f"validation-a{attempt}"
            group = cohort(case, phase, rate["m"], n3, args.size, work, args.years)
            summary, r, z, s, d = validation_tables(
                case, args.years, rate, group, resamples
            )
            summary["attempt"] = attempt
            summary["seconds"] = round(sum(f["seconds"] for f in group), 1)
            summaries.append(summary)
            regions += [{**x, "attempt": attempt} for x in r]
            zones += [{**x, "attempt": attempt} for x in z]
            segments += [{**x, "attempt": attempt} for x in s]
            draws[case] = d
            log(json.dumps(summary))
            if summary["band"] != "outside":
                break
    if "main" in draws and "neutral" in draws:
        rng = np.random.default_rng(seed_for("ratio", args.years))
        a, b = np.asarray(draws["main"]), np.asarray(draws["neutral"])
        n = min(len(a), len(b))
        ratios = rng.permutation(a)[:n] / rng.permutation(b)[:n]
        main = [s for s in summaries if s["case"] == "main"][-1]
        neutral = [s for s in summaries if s["case"] == "neutral"][-1]
        low, high = np.quantile(ratios, [0.025, 0.975])
        summaries.append(
            {
                "case": "main/neutral",
                "years": args.years,
                "attempt": max(main["attempt"], neutral["attempt"]),
                "types_mean": main["types_mean"] / neutral["types_mean"],
                "types_low": float(low),
                "types_high": float(high),
            }
        )
    tables = {
        "de-novo-rates.csv": rates,
        "de-novo-fit.csv": fits,
        "de-novo-summary.csv": summaries,
        "de-novo-regions.csv": regions,
        "de-novo-zones.csv": zones,
        "de-novo-segments.csv": segments,
    }
    for name, rows in tables.items():
        pd.DataFrame(rows).to_csv(
            output / name, index=False, lineterminator="\n", float_format="%.6g"
        )
    log(f"tables written to {output}")


SMALL = {"nss": 20, "c": 0.01, "mu": 0.07, "gamma": 0.1, "delta": 0.2, "m": 0.02}


def cdf_gap(a, b):
    """Largest distance between the empirical CDFs of two integer samples."""
    last = int(max(np.max(a), np.max(b)))
    fa = np.bincount(np.asarray(a, dtype=int), minlength=last + 1).cumsum() / len(a)
    fb = np.bincount(np.asarray(b, dtype=int), minlength=last + 1).cumsum() / len(b)
    return float(np.max(np.abs(fa - fb)))


def small_paths(n, seed, regions=3, days=60.0, compiled=True, **model):
    """Endpoints of ``n`` small fibers: W0, M0, M1, types in region 0 and in all."""
    p = {**SMALL, **model}
    rows = []
    for k in range(n):
        out = simulate(
            regions,
            p["m"],
            days,
            seed + k,
            nss=p["nss"],
            c=p["c"],
            mu=p["mu"],
            gamma=p["gamma"],
            delta=p["delta"],
            compiled=compiled,
        )
        rows.append(
            (
                out["wild"][0],
                out["mutant"][0],
                out["mutant"][min(1, regions - 1)],
                len(np.unique(out["labels"][0])),
                len(np.unique(np.concatenate(out["labels"]))),
                out["mutant"].sum(),
                out["events"][3],
            )
        )
    return np.array(rows)


def crn_paths(n, seed, regions=3, days=60.0):
    """The same endpoints from the CRN-craft model (tree method)."""
    import crn_craft as crn

    p = SMALL
    model = crn_model(
        regions, p["nss"], p["c"], p["mu"], p["gamma"], p["delta"], p["m"]
    )
    result = crn.ensemble(
        model,
        end=days,
        report_times=[0.0, days],
        replicates=n,
        seed=seed,
        keep_paths="all",
        keep_instances=n,
        strict=True,
    ).simulate()
    final = result.paths.isel(time=-1)
    pick = ["W0", "M0", f"M{min(1, regions - 1)}", "T0"]
    columns = [final.sel(variable=v).to_numpy() for v in pick]
    # The observable "types" counts instances, that is (clone, region) pairs;
    # distinct clones over the fiber come from the final instances.
    clones = np.zeros(n)
    for cell in result.metadata["final_instances"]:
        found = {i["attributes"][0] for i in cell["instances"] if i["count"] > 0}
        clones[cell["id"]] = len(found - {0})
    return np.column_stack([*columns, clones])


def sampler_checks(args):
    """V1 to V4 of the predeclared sampler validation (see the guide page)."""
    import pandas as pd

    rows = []
    n = args.paths
    alpha = 0.001
    bound = 2 * math.sqrt(math.log(2 / (alpha / 2)) / (2 * n))
    start = time.perf_counter()
    ours = small_paths(n, seed_for("V1", "sampler"), compiled=_fast is not None)
    t_ours = time.perf_counter() - start
    start = time.perf_counter()
    theirs = crn_paths(n, seed_for("V1", "crn") % 2**31)
    t_crn = time.perf_counter() - start
    names = ["W0", "M0", "M1", "types in region 0", "types in fiber"]
    for j, name in enumerate(names):
        gap = cdf_gap(ours[:, j], theirs[:, j])
        rows.append(
            {
                "test": "V1 CDF gap, sampler vs CRN-craft",
                "statistic": name,
                "paths": n,
                "value": gap,
                "bound": bound,
                "mean_sampler": ours[:, j].mean(),
                "mean_reference": theirs[:, j].mean(),
                "passed": gap <= bound,
            }
        )
    # V2: c = 0 holds every birth rate at mu, so E[N] stays at L * Nss.
    L, nss, m, mu, t = 4, 50, 0.05, 0.07, 100.0
    sample = small_paths(
        args.moment_paths,
        seed_for("V2"),
        regions=L,
        days=t,
        nss=nss,
        c=0.0,
        m=m,
        compiled=_fast is not None,
    )
    for name, values, expected in (
        ("mutant copies", sample[:, 5], L * nss * (1 - math.exp(-m * mu * t))),
        ("mutation events", sample[:, 6], m * mu * L * nss * t),
    ):
        z = (values.mean() - expected) / (values.std(ddof=1) / math.sqrt(len(values)))
        rows.append(
            {
                "test": "V2 moment law, c = 0",
                "statistic": name,
                "paths": len(values),
                "value": z,
                "bound": 3.29,
                "mean_sampler": values.mean(),
                "mean_reference": expected,
                "passed": abs(z) < 3.29,
            }
        )
    rows.append(
        {
            "test": "timing",
            "statistic": "seconds sampler / CRN-craft tree",
            "paths": n,
            "value": t_ours,
            "bound": float("nan"),
            "mean_sampler": t_ours,
            "mean_reference": t_crn,
            "passed": True,
        }
    )
    frame = pd.DataFrame(rows)
    frame.to_csv(
        args.output / "de-novo-sampler.csv",
        index=False,
        lineterminator="\n",
        float_format="%.6g",
    )
    print(frame.to_string())


def sensitivity(args):
    """Fibers of another length at a frozen m: affected share and types."""
    import pandas as pd

    rows = []
    for k in range(args.fibers):
        seed = seed_for("length", args.case, args.regions, f"{args.m:.8e}", k)
        out = simulate(
            args.regions, args.m, DAYS_PER_YEAR * args.years, seed, **CASES[args.case]
        )
        s = summarize(out)
        for r in s["regions"]:
            rows.append(
                {"case": args.case, "chain": args.regions, "m": args.m, "fiber": k, **r}
            )
        print(k, sum(r["affected"] for r in s["regions"]), flush=True)
    frame = pd.DataFrame(rows)
    file = args.output / "de-novo-length.csv"
    if file.exists():
        frame = pd.concat([pd.read_csv(file), frame])
    frame.to_csv(file, index=False, lineterminator="\n", float_format="%.6g")


COLORS = {"main": "#2a78d6", "neutral": "#eb6834", "isolated": "#1baf7a"}
LABELS = {
    "main": "Denser mutant, 32 regions",
    "neutral": "Equal density, 32 regions",
    "isolated": "Denser mutant, isolated regions",
}


def figure(args):
    """Draw the calibration and the type counts from the saved tables."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    import pandas as pd

    fit = pd.read_csv(args.results / "de-novo-fit.csv")
    fit = fit[fit.years == 3]
    rates = pd.read_csv(args.results / "de-novo-rates.csv")
    rates = rates[rates.years == 3].sort_values("attempt").groupby("case").last()
    regions = pd.read_csv(args.results / "de-novo-regions.csv")
    regions = regions[regions.years == 3]
    kk14 = pd.read_csv(
        ROOT.parent / "kk14" / "results" / "reconstruction" / "comparison.csv"
    )
    kk14 = kk14[
        (kk14.label == "primary") & (kk14.years == 3) & (kk14.cohort == "validation")
    ].set_index("mechanism")
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(9.2, 4.0), layout="constrained", width_ratios=(1.1, 1)
    )
    for case, color in COLORS.items():
        if case not in rates.index:
            continue
        part = (
            fit[fit.case == case]
            .groupby("m", as_index=False)[["affected", "regions"]]
            .sum()
        )
        left.plot(
            part.m,
            part.affected / part.regions,
            "o",
            color=color,
            ms=6,
            label=LABELS[case],
        )
        r = rates.loc[case]
        grid = np.geomspace(part.m.min() / 1.05, part.m.max() * 1.05, 200)
        curve = 1 / (1 + np.exp(-(r.intercept + r.slope * np.log10(grid))))
        left.plot(grid, curve, color=color, linewidth=2)
    left.axhline(TARGET, color=".55", linewidth=1, linestyle="--")
    left.set(
        xscale="log",
        xlabel="Deletion probability per replication, m",
        ylabel="Share of regions with ≥ 60% mutants",
        title="A  Fit of m (3 years)",
    )
    ticks = np.arange(0.008, 0.0141, 0.001)
    left.set_xticks(ticks, [f"{t:.3f}" for t in ticks])
    left.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    left.legend(frameon=False, fontsize=8, loc="upper left")
    last = regions.groupby("case").attempt.transform("max")
    hit = regions[regions.affected & (regions.attempt == last)]
    order = [c for c in COLORS if c in set(hit.case)]
    data = [hit[hit.case == c].types.to_numpy() for c in order]
    parts = right.boxplot(data, widths=0.5, patch_artist=True, showfliers=True)
    for patch, case in zip(parts["boxes"], order, strict=True):
        patch.set(facecolor=COLORS[case], alpha=0.75, edgecolor=".25")
    for line in parts["medians"]:
        line.set(color=".15", linewidth=1.5)
    for name, mech in (("drift", "Neutral drift"), ("transcription", "Transcription")):
        value = kk14.loc[name, "mutant_types"]
        right.axhline(value, color=".45", linewidth=1, linestyle=":")
        right.text(
            len(order) + 0.45,
            value * 1.08,
            f"{mech}, one cell (ch. 17): {value:.3g}",
            ha="right",
            fontsize=8,
            color=".3",
        )
    right.set_xticks(
        range(1, len(order) + 1),
        [LABELS[c].replace(", ", "\n") for c in order],
        fontsize=8,
    )
    right.set(
        yscale="log",
        ylim=(1, 2.5 * max(max(map(np.max, data)), kk14.mutant_types.max())),
        ylabel="Mutant types per affected region",
        title="B  Types in affected regions",
    )
    args.figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figures / "de-novo-types.png", dpi=155)
    plt.close(fig)


def pilot(args):
    """Time single fibers to size the experiment."""
    for regions, days, labels in args.cases:
        if _fast is not None:
            simulate(2, 0.01, 1.0, 1, labels=labels)  # compile first
        start = time.perf_counter()
        out = simulate(
            regions,
            args.m,
            days,
            MANIFEST["seed"] + 7000,
            labels=labels,
            gamma=args.gamma,
            delta=args.delta,
        )
        seconds = time.perf_counter() - start
        n = out["wild"] + out["mutant"]
        types = [len(np.unique(x)) for x in out["labels"]] if labels else [0]
        events = int(out["events"][:3].sum())
        row = {
            "regions": regions,
            "days": days,
            "labels": labels,
            "compiled": _fast is not None,
            "seconds": round(seconds, 3),
            "events": events,
            "ns_per_event": round(1e9 * seconds / events, 1),
            "mean_h": round(float((out["mutant"] / n).mean()), 4),
            "affected": int(affected(out["wild"], out["mutant"]).sum()),
            "mean_types": round(float(np.mean(types)), 1),
        }
        print(json.dumps(row), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run", help="calibrate m, then simulate the validation cohort")
    p.add_argument("--size", choices=tuple(SIZES), default="quick")
    p.add_argument("--years", type=int, default=3)
    p.add_argument("--case", dest="cases", action="append", choices=tuple(CASES))
    p.add_argument("--work", type=Path, help="checkpoint folder")
    p.add_argument("--output", type=Path, help="folder for the result tables")
    p = sub.add_parser("pilot", help="time single fibers")
    p.add_argument("--m", type=float, default=0.011)
    p.add_argument("--gamma", type=float)
    p.add_argument("--delta", type=float)
    p.add_argument(
        "--case",
        dest="cases",
        nargs=3,
        action="append",
        type=float,
        metavar=("REGIONS", "DAYS", "LABELS"),
    )
    p = sub.add_parser("sampler", help="sampler against CRN-craft and a moment law")
    p.add_argument("--paths", type=int, default=20000)
    p.add_argument("--moment-paths", type=int, default=4000)
    p.add_argument("--output", type=Path, default=ROOT / "runs" / "de-novo")
    p = sub.add_parser("figure", help="draw the figure from saved tables")
    p.add_argument("--results", type=Path, default=ROOT / "results")
    p.add_argument("--figures", type=Path, default=ROOT / "runs" / "de-novo")
    p = sub.add_parser("length", help="fibers of another chain length at fixed m")
    p.add_argument("--case", choices=tuple(CASES), default="main")
    p.add_argument("--m", type=float, required=True)
    p.add_argument("--regions", type=int, default=64)
    p.add_argument("--fibers", type=int, default=8)
    p.add_argument("--years", type=int, default=3)
    p.add_argument("--output", type=Path, default=ROOT / "runs" / "de-novo")
    args = parser.parse_args()
    if args.command == "length":
        args.output.mkdir(parents=True, exist_ok=True)
        sensitivity(args)
    elif args.command == "figure":
        figure(args)
    elif args.command == "sampler":
        args.output.mkdir(parents=True, exist_ok=True)
        sampler_checks(args)
    elif args.command == "run":
        args.cases = args.cases or list(CASES)
        run(args)
    elif args.command == "pilot":
        args.cases = [(int(a), b, bool(c)) for a, b, c in (args.cases or [(8, 20, 1)])]
        pilot(args)


if __name__ == "__main__":
    main()
