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


POOL = None  # process pool of ``run --workers N``; None runs fibers in turn


def prefetch(case, cohorts, size, work, years):
    """Simulate the fibers of several cohorts at once and checkpoint them.

    ``cohorts`` holds (phase, m, fibers). Every fiber has its own seed and
    its own checkpoint file, so the order of simulation does not change any
    result. Without a pool this does nothing; the fibers are then simulated
    when ``cohort`` asks for them.
    """
    if POOL is None:
        return
    jobs = [
        (case, phase, k, m, size, work, years)
        for phase, m, fibers in cohorts
        for k in range(fibers)
    ]
    list(POOL.map(fiber, *zip(*jobs, strict=True)))


def cohort(case, phase, m, fibers, size, work, years):
    prefetch(case, [(phase, m, fibers)], size, work, years)
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
    prefetch(case, [("stage1", m0 * f, n1) for f in GRID], size, work, years)
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
    prefetch(case, [("stage2", m1 * f, n2) for f in LOCAL], size, work, years)
    for factor in LOCAL:
        add("stage2", m1 * factor, n2)
    for repeat in range(2, attempt + 1):
        previous = logistic([(r[1], r[3], r[4]) for r in rows])[2]
        stage = f"stage2-a{repeat}"
        prefetch(case, [(stage, previous * f, n2) for f in LOCAL], size, work, years)
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


def validation_tables(case, years, rate, group, resamples, tag="bootstrap"):
    """Per-unit rows and the summary row of one validation cohort.

    ``tag`` names the bootstrap stream; the fiber design uses its own.
    """
    regions, zones, segments = [], [], []
    for f in group:
        for r in f["regions"]:
            regions.append({"case": case, "years": years, "fiber": f["index"], **r})
        for z in f["zones"]:
            zones.append({"case": case, "years": years, "fiber": f["index"], **z})
        segments.append(
            {"case": case, "years": years, "fiber": f["index"], **f["segment"]}
        )
    rng = np.random.default_rng(seed_for(tag, case, years, rate["attempt"]))
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


def ratio_row(draws, summaries, years, tag):
    """Main over neutral types per affected region, if both cases ran."""
    if "main" not in draws or "neutral" not in draws:
        return
    rng = np.random.default_rng(seed_for(tag, years))
    a, b = np.asarray(draws["main"]), np.asarray(draws["neutral"])
    n = min(len(a), len(b))
    if n == 0:
        return
    ratios = rng.permutation(a)[:n] / rng.permutation(b)[:n]
    main = [s for s in summaries if s["case"] == "main"][-1]
    neutral = [s for s in summaries if s["case"] == "neutral"][-1]
    low, high = np.quantile(ratios, [0.025, 0.975])
    summaries.append(
        {
            "case": "main/neutral",
            "years": years,
            "attempt": max(main["attempt"], neutral["attempt"]),
            "types_mean": main["types_mean"] / neutral["types_mean"],
            "types_low": float(low),
            "types_high": float(high),
        }
    )


def run(args):
    import pandas as pd

    if args.design == "fibers":
        return run_fibers(args)
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
    ratio_row(draws, summaries, args.years, "ratio")
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


# The fiber design (`run --design fibers`, manifest entry `de_novo_fibers`):
# the same target, model and units, but fibers are the sampling unit in the
# fit and in the validation. Fibers per stage-1 value, per stage-2 value and
# in a validation cohort; bootstrap resamples; affected fibers a compatible
# cohort needs.
FIBER_SIZES = {
    "quick": {"fibers": (4, 8, 32), "boot": 199, "min_fibers": 3},
    "draft": {"fibers": (2, 4, 16), "boot": 999, "min_fibers": 3},
    "full": {"fibers": (16, 32, 128), "boot": 1999, "min_fibers": 10},
}
FIBER_GRID = (0.92, 0.97, 1.03, 1.09)
FIBER_LOCAL = (1 / 1.03, 1.0, 1.03)
# CPU-seconds per region and year with Numba, Figure 2 parameters: the
# second machine (Intel Xeon 2.3 GHz) took about 170 s for a 32-region fiber
# at 40 years and 330 s at 80 years. An estimate, not a measurement here.
SECONDS_PER_REGION_YEAR = 0.13
# The rule after a failed fit at attempt 1: "R0" as first declared (the
# calibration ends), "R1" (the refit follows); see `fiber_calibration`.
FIBER_REVISION = "R0"
FULL_PRECISION = ("m", "m_low95", "m_high95", "m_boot_low95", "m_boot_high95")


def logistic_many(x, k, n):
    """Newton fits of the logistic curve of ``logistic`` to many data sets.

    ``x`` holds log10(m) per fit point, ``k`` and ``n`` affected and all
    regions with one row per data set. Returns intercepts, slopes and
    whether each fit converged to a finite curve with a positive slope.
    """
    x = np.asarray(x, dtype=float)
    k = np.atleast_2d(np.asarray(k, dtype=float))
    n = np.atleast_2d(np.asarray(n, dtype=float))
    a = np.zeros(len(k))
    b = np.zeros(len(k))
    p = k.sum(axis=1) / n.sum(axis=1)
    inside = (p > 0) & (p < 1)
    a[inside] = np.log(p[inside] / (1 - p[inside]))
    done = np.zeros(len(k), dtype=bool)
    with np.errstate(all="ignore"):
        for _ in range(100):
            prob = 1 / (1 + np.exp(-(a[:, None] + b[:, None] * x)))
            w = n * prob * (1 - prob)
            r = k - n * prob
            i00 = w.sum(axis=1) + 1e-12
            i01 = (w * x).sum(axis=1)
            i11 = (w * x * x).sum(axis=1) + 1e-12
            g0, g1 = r.sum(axis=1), (r * x).sum(axis=1)
            det = i00 * i11 - i01 * i01
            da = (i11 * g0 - i01 * g1) / det
            db = (i00 * g1 - i01 * g0) / det
            a = np.where(done, a, a + da)
            b = np.where(done, b, b + db)
            done |= np.maximum(np.abs(da), np.abs(db)) < 1e-10
            if done.all():
                break
    ok = done & np.isfinite(a) & np.isfinite(b) & (b > 0)
    return a, b, ok


def fiber_fit(points, resamples, rng):
    """Point estimate and intervals of the m at TARGET from fiber outcomes.

    ``points`` holds (stage, m, hits, sizes): per fit point the affected
    regions and the regions of every fiber. The point estimate is the
    logistic fit through all regions; ``m_low95``/``m_high95`` is its
    delta-method interval (regions independent, as in the first design);
    ``m_boot_low95``/``m_boot_high95`` is the 95% percentile interval of
    refits on resamples of whole fibers within each fit point. A resample
    whose fit fails counts at both ends, so it can only widen the interval.
    A point fit without a positive slope, or whose inverse lies outside the
    range of the fit points, fails (``fit_ok`` False).
    """
    x = np.log10([p[1] for p in points])
    k = np.array([np.sum(p[2]) for p in points], dtype=float)
    n = np.array([np.sum(p[3]) for p in points], dtype=float)
    a, b, ok = logistic_many(x, k, n)
    a, b, ok = float(a[0]), float(b[0]), bool(ok[0])
    logit = math.log(TARGET / (1 - TARGET))
    ok = ok and x.min() <= (logit - a) / b <= x.max()
    out = {"intercept": a, "slope": b, "fit_ok": ok}
    nan = float("nan")
    if not ok:
        out.update(m=nan, m_low95=nan, m_high95=nan)
        out.update(m_boot_low95=nan, m_boot_high95=nan, boot_failed=resamples)
        return out
    prob = 1 / (1 + np.exp(-(a + b * x)))
    w = n * prob * (1 - prob)
    design = np.column_stack([np.ones_like(x), x])
    cov = np.linalg.inv(design.T @ (design * w[:, None]))
    x0 = (logit - a) / b
    grad = np.array([-1 / b, -(logit - a) / b**2])
    se = math.sqrt(grad @ cov @ grad)
    out.update(
        m=10**x0, m_low95=10 ** (x0 - 1.96 * se), m_high95=10 ** (x0 + 1.96 * se)
    )
    kb = np.zeros((resamples, len(points)))
    nb = np.zeros((resamples, len(points)))
    for j, (_, _, hits, sizes) in enumerate(points):
        hits, sizes = np.asarray(hits), np.asarray(sizes)
        pick = rng.integers(0, len(hits), (resamples, len(hits)))
        kb[:, j] = hits[pick].sum(axis=1)
        nb[:, j] = sizes[pick].sum(axis=1)
    ab, bb, okb = logistic_many(x, kb, nb)
    with np.errstate(all="ignore"):
        xb = (logit - ab) / bb
    with np.errstate(all="ignore"):
        low = np.quantile(np.where(okb, xb, -np.inf), 0.025)
        high = np.quantile(np.where(okb, xb, np.inf), 0.975)
    out.update(
        m_boot_low95=float(10**low) if np.isfinite(low) else 0.0,
        m_boot_high95=float(10**high) if np.isfinite(high) else float("inf"),
        boot_failed=int((~okb).sum()),
    )
    return out


def fiber_stages(start, sample, n1, n2, log=print):
    """Stages 1 and 2 of the fiber design; returns the fit points.

    ``sample(stage, ms, fibers)`` returns, per value of m, the affected
    regions and the regions of each fiber. The same code runs on simulated
    fibers and on the synthetic fibers of the coverage experiment.
    """
    points = []

    def add(stage, ms, fibers):
        for m, (hits, sizes) in zip(ms, sample(stage, ms, fibers), strict=True):
            points.append((stage, m, list(hits), list(sizes)))
            log(f"{stage} m={m:.6e} affected {sum(hits)}/{sum(sizes)}")

    add("stage1", [start * f for f in FIBER_GRID], n1)
    for _ in range(4):  # extend the grid until it brackets the target
        shares = [sum(p[2]) / sum(p[3]) for p in points]
        if min(shares) > TARGET:
            add("stage1", [min(p[1] for p in points) / 1.25], n1)
        elif max(shares) < TARGET:
            add("stage1", [max(p[1] for p in points) * 1.25], n1)
        else:
            break
    x = np.log10([p[1] for p in points])
    k = [sum(p[2]) for p in points]
    n = [sum(p[3]) for p in points]
    a, b, ok = logistic_many(x, k, n)
    ok = bool(ok[0])
    if ok:
        x1 = (math.log(TARGET / (1 - TARGET)) - a[0]) / b[0]
        ok = x.min() <= x1 <= x.max()
    if ok:
        m1 = 10**x1
    else:  # no usable inverse: the geometric mean of the stage-1 values
        m1 = 10 ** float(np.mean(x))
        log("stage 1 fit without a usable inverse; stage 2 at the geometric mean")
    add("stage2", [m1 * f for f in FIBER_LOCAL], n2)
    return points, ok


def refit_center(points):
    """Center of the refit: the inverse of the fit through all points.

    Limited to the fit values widened by a factor 1.25 on each side; the
    geometric mean of the fit values if that fit has no positive slope.
    After a successful fit this is the fitted m itself.
    """
    x = np.log10([p[1] for p in points])
    a, b, ok = logistic_many(
        x, [sum(p[2]) for p in points], [sum(p[3]) for p in points]
    )
    if not (np.isfinite(a[0]) and np.isfinite(b[0]) and b[0] > 0):
        return 10 ** float(np.mean(x))
    x0 = (math.log(TARGET / (1 - TARGET)) - a[0]) / b[0]
    wide = math.log10(1.25)
    return 10 ** float(min(max(x0, x.min() - wide), x.max() + wide))


def fiber_refit(points, sample, n2, log=print):
    """Three new stage-2 values around ``refit_center``."""
    center = refit_center(points)
    ms = [center * f for f in FIBER_LOCAL]
    for m, (hits, sizes) in zip(ms, sample("stage2-a2", ms, n2), strict=True):
        points.append(("stage2-a2", m, list(hits), list(sizes)))
        log(f"stage2-a2 m={m:.6e} affected {sum(hits)}/{sum(sizes)}")
    return points


def fallback_point(points):
    """The m of the fallback cohort after a failed fit.

    Among the fit points of the last stage 2 that was run, the one whose
    share of affected regions is closest to TARGET; on a tie, the smaller m.
    """
    last = "stage2-a2" if any(p[0] == "stage2-a2" for p in points) else "stage2"
    pool = [p for p in points if p[0] == last]
    best = min(pool, key=lambda p: (abs(sum(p[2]) / sum(p[3]) - TARGET), p[1]))
    return best[1]


def fiber_calibration(start, sample, judge, spec, rng_for, log, revision):
    """Fit, validation and at most one refit; the same code for simulated
    and synthetic fibers.

    ``judge(attempt, fit)`` validates a frozen m and returns the label of
    the cohort. ``rng_for(attempt)`` gives the bootstrap stream of a fit.
    ``revision`` "R0" is the design as first declared: a failed fit at
    attempt 1 ends the calibration. "R1" refits after it instead. Returns
    the fit points, whether stage 1 gave a usable inverse, and per attempt
    (attempt, fit, label).
    """
    n1, n2, _ = spec["fibers"]
    points, stage1_ok = fiber_stages(start, sample, n1, n2, log)
    attempts = []
    for attempt in (1, 2):
        if attempt == 2:
            points = fiber_refit(points, sample, n2, log)
        fit = fiber_fit(points, spec["boot"], rng_for(attempt))
        log(
            f"frozen m = {fit['m']:.6e}, fiber bootstrap"
            f" {fit['m_boot_low95']:.6e} to {fit['m_boot_high95']:.6e}"
            f" (attempt {attempt})"
        )
        label = judge(attempt, fit) if fit["fit_ok"] else "fit failed"
        attempts.append((attempt, fit, label))
        if label == "compatible":
            break
        if label == "fit failed" and revision == "R0":
            break
    return points, stage1_ok, attempts


def saved_start(case, size, years):
    """The last fitted m of the first design for this case and lifespan.

    Only the draft and full sizes simulate the model of the saved tables.
    """
    if size == "quick":
        return None
    file = ROOT / "results" / "de-novo-rates.csv"
    days = DAYS_PER_YEAR * years
    import csv

    with file.open(newline="") as handle:
        found = [
            r
            for r in csv.DictReader(handle)
            if r["case"] == case and int(r["years"]) == years
            if float(r["days"]) == days
        ]
    if not found:
        return None
    return float(max(found, key=lambda r: int(r["attempt"]))["m"])


def fiber_counts(group):
    """Per fiber: affected regions and all regions."""
    hits = [sum(r["affected"] for r in f["regions"]) for f in group]
    return hits, [len(f["regions"]) for f in group]


def decide(low, high, some, min_fibers):
    """The label of a validation cohort in the fiber design."""
    if not low <= TARGET <= high:
        return "outside"
    if some < min_fibers:
        return "too few affected fibers"
    return "compatible"


def fiber_decision(summary, group, min_fibers):
    """The validation label of the fiber design, added to a summary row."""
    hits, sizes = fiber_counts(group)
    some = sum(h > 0 for h in hits)
    summary["fibers_affected"] = some
    summary["fibers_full"] = sum(h == s for h, s in zip(hits, sizes, strict=True))
    summary["min_fibers"] = min_fibers
    summary["decision"] = decide(
        summary["fiber_boot_low"], summary["fiber_boot_high"], some, min_fibers
    )
    if summary["statement"] != "too few affected regions" and some < min_fibers:
        summary["statement"] = "too few affected fibers"
    return summary


def fiber_plan(args):
    """Print the fibers the run simulates and their expected cost."""
    n1, n2, n3 = FIBER_SIZES[args.size]["fibers"]
    regions = SIZES[args.size]["regions"]
    plain = refit = most = 0
    for case in args.cases:
        probe = saved_start(case, args.size, args.years) is None
        # Each attempt ends with one cohort of n3 fibers: a validation cohort,
        # or after a failed fit the fallback cohort, never both.
        plain += probe + 4 * n1 + 3 * n2 + n3
        refit += probe + 4 * n1 + 6 * n2 + 2 * n3
        most += probe + 8 * n1 + 6 * n2 + 2 * n3
    print(
        f"fiber design, {args.size} size, {args.years} years,"
        f" cases {', '.join(args.cases)}: {plain} fibers without a refit,"
        f" {refit} with it, at most {most} (grid extended four times); a"
        " fallback cohort after a failed fit takes the place of a validation"
        " cohort, so it adds nothing to these counts",
        flush=True,
    )
    if args.size == "quick":
        print("quick model (small regions): seconds in pure Python", flush=True)
        return
    each = SECONDS_PER_REGION_YEAR * regions * args.years
    print(
        f"expected cost with Numba (estimate from a 2.3 GHz Xeon,"
        f" {each:.0f} CPU-seconds per fiber): {plain * each / 3600:.1f}"
        f" CPU-hours without a refit, {refit * each / 3600:.1f} with it,"
        f" at most {most * each / 3600:.1f}",
        flush=True,
    )
    if _fast is None:
        print("Numba is not installed: expect about 100 times longer", flush=True)


def write_tables(tables, output):
    """CSV tables; the fiber design keeps every m with full precision."""
    import pandas as pd

    for name, rows in tables.items():
        # Columns in the order of the fullest row (a "fit failed" row has few).
        order = list(max(rows, key=len)) if rows else []
        order += [k for r in rows for k in r if k not in order]
        frame = pd.DataFrame(rows, columns=list(dict.fromkeys(order)))
        for column in FULL_PRECISION:
            if column in frame and "fibers" in name:
                frame[column] = [
                    repr(float(v)) if v == v else "" for v in frame[column]
                ]
        frame.to_csv(
            output / name, index=False, lineterminator="\n", float_format="%.6g"
        )


def run_fibers(args):
    """The fiber design: fit, validation and at most one refit per case."""
    size = args.size
    work = args.work or ROOT / "runs" / "de-novo-fibers" / size
    output = args.output or work
    output.mkdir(parents=True, exist_ok=True)
    spec = FIBER_SIZES[size]
    n1, n2, n3 = spec["fibers"]
    resamples = spec["boot"]
    years = args.years
    days = SIZES[size]["days"] or DAYS_PER_YEAR * years
    fiber_plan(args)

    def log(text):
        print(text, flush=True)

    rates, fits, summaries, regions, zones, segments, draws = [], [], [], [], [], [], {}
    for case in args.cases:

        def sample(stage, ms, fibers, case=case):
            phase = f"fibers-{stage}"
            prefetch(case, [(phase, m, fibers) for m in ms], size, work, years)
            return [
                fiber_counts(cohort(case, phase, m, fibers, size, work, years))
                for m in ms
            ]

        start = saved_start(case, size, years)
        source = "first design"
        if start is None:
            source = "probe"
            mu = MANIFEST["fig2_site_parameters"]["death"]
            guess = -math.log(0.4) / (mu * days)
            probe = fiber(case, "fibers-stage0", 0, guess, size, work, years)
            h = np.array(
                [r["mutant"] / (r["wild"] + r["mutant"]) for r in probe["regions"]]
            )
            wanted = min(0.6 - 1.2816 * h.std(ddof=1), 0.59)
            start = guess * math.log(1 - wanted) / math.log(1 - h.mean())
        log(f"{case}: start m = {start:.6e} ({source})")

        def judge(attempt, fit, case=case):
            rate = {"m": fit["m"], "days": days, "attempt": attempt}
            phase = "fibers-validation" if attempt == 1 else "fibers-validation-a2"
            group = cohort(case, phase, fit["m"], n3, size, work, years)
            summary, r, z, s, d = validation_tables(
                case, years, rate, group, resamples, tag="fibers-bootstrap"
            )
            summary = fiber_decision(summary, group, spec["min_fibers"])
            summary["attempt"] = attempt
            summary["cohort"] = "validation"
            summary["seconds"] = round(sum(f["seconds"] for f in group), 1)
            summaries.append(summary)
            tag = {"attempt": attempt, "cohort": "validation"}
            regions.extend({**x, **tag} for x in r)
            zones.extend({**x, **tag} for x in z)
            segments.extend({**x, **tag} for x in s)
            draws[case] = d
            log(json.dumps(summary))
            return summary["decision"]

        def rng_for(attempt, case=case):
            return np.random.default_rng(seed_for("fibers-fit", case, years, attempt))

        first = len(summaries)
        points, stage1_ok, attempts = fiber_calibration(
            start, sample, judge, spec, rng_for, log, FIBER_REVISION
        )
        for attempt, fit, label in attempts:
            rates.append(
                {
                    "case": case,
                    "years": years,
                    "m": fit["m"],
                    "m_low95": fit["m_low95"],
                    "m_high95": fit["m_high95"],
                    "intercept": fit["intercept"],
                    "slope": fit["slope"],
                    "start_m": start,
                    "days": days,
                    "attempt": attempt,
                    "m_boot_low95": fit["m_boot_low95"],
                    "m_boot_high95": fit["m_boot_high95"],
                    "boot_failed": fit["boot_failed"],
                    "fit_ok": fit["fit_ok"],
                    "start_source": source,
                    "stage1_fit_ok": stage1_ok,
                    "revision": FIBER_REVISION,
                }
            )
            if label == "fit failed":
                summaries.append(
                    {"case": case, "years": years, "attempt": attempt}
                    | {"decision": label, "cohort": "none"}
                )
        summaries[first:] = sorted(summaries[first:], key=lambda r: r["attempt"])
        done = attempts[-1][2] == "compatible"
        for i, row in enumerate(summaries[first:]):
            last = first + i == len(summaries) - 1
            row["calibration"] = ("done" if done else "failed") if last else "refit"
        if attempts[-1][2] == "fit failed":
            # The declared fallback: one cohort of validation size at the
            # fit point of the last stage 2 closest to 10% affected regions.
            attempt = attempts[-1][0]
            m = fallback_point(points)
            group = cohort(case, "fibers-fallback", m, n3, size, work, years)
            rate = {"m": m, "days": days, "attempt": attempt}
            summary, r, z, s, _ = validation_tables(
                case, years, rate, group, resamples, tag="fibers-fallback"
            )
            summary = fiber_decision(summary, group, spec["min_fibers"])
            summary["decision"] = "calibration failed, fallback cohort"
            summary["attempt"] = attempt
            summary["cohort"] = "fallback"
            summary["calibration"] = "failed"
            summary["seconds"] = round(sum(f["seconds"] for f in group), 1)
            summaries.append(summary)
            tag = {"attempt": attempt, "cohort": "fallback"}
            regions.extend({**x, **tag} for x in r)
            zones.extend({**x, **tag} for x in z)
            segments.extend({**x, **tag} for x in s)
            log(json.dumps(summary))
        for stage, m, hits, sizes in points:
            fits.append(
                {
                    "case": case,
                    "years": years,
                    "stage": stage,
                    "m": m,
                    "fibers": len(hits),
                    "affected": sum(hits),
                    "regions": sum(sizes),
                    "fibers_affected": sum(h > 0 for h in hits),
                    "fibers_full": sum(
                        h == s for h, s in zip(hits, sizes, strict=True)
                    ),
                }
            )
    ratio_row(draws, summaries, years, "fibers-ratio")
    write_tables(
        {
            "de-novo-fibers-rates.csv": rates,
            "de-novo-fibers-fit.csv": fits,
            "de-novo-fibers-summary.csv": summaries,
            "de-novo-fibers-regions.csv": regions,
            "de-novo-fibers-zones.csv": zones,
            "de-novo-fibers-segments.csv": segments,
        },
        output,
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


def length_fiber(case, regions, m, index, size, years, work):
    """One fiber of another chain length, checkpointed like ``fiber``.

    The 3-year full-size seeds are those of the saved 64-region table; other
    lifespans and sizes add both to the seed.
    """
    spec = SIZES[size]
    name = f"{size}-{case}-{regions}-{years}y-{m:.8e}-{index:03d}.json"
    file = work / "length" / name
    if file.is_file():
        return json.loads(file.read_text())
    if size == "full" and years == 3:
        seed = seed_for("length", case, regions, f"{m:.8e}", index)
    else:
        seed = seed_for("length", case, regions, size, years, f"{m:.8e}", index)
    days = spec["days"] or DAYS_PER_YEAR * years
    start = time.perf_counter()
    out = simulate(
        regions,
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
            "phase": "length",
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


def last_row(path, case, years):
    """The row of the last attempt for a case and lifespan in a result table."""
    import csv

    with Path(path).open(newline="") as handle:
        found = [
            r
            for r in csv.DictReader(handle)
            if r["case"] == case and int(r["years"]) == years
        ]
    if not found:
        raise SystemExit(f"no row for {case}, {years} years in {path}")
    return max(found, key=lambda r: int(r["attempt"]))


def sensitivity(args):
    """Fibers of another length at a frozen m: affected share and types.

    Writes the region rows to ``de-novo-length.csv`` (appended, as before)
    and one summary row in the layout of the fiber design to
    ``de-novo-length-summary.csv`` (appended). With ``--reference`` (a
    fiber-design summary) the row states whether chain length matters.
    """
    import pandas as pd

    m = args.m
    if args.m_from is not None:
        m = float(last_row(args.m_from, args.case, args.years)["m"])
    if m is None:
        raise SystemExit("length needs --m or --m-from")
    work = args.work or args.output
    jobs = [
        (args.case, args.regions, m, k, args.size, args.years, work)
        for k in range(args.fibers)
    ]
    if POOL is not None:
        list(POOL.map(length_fiber, *zip(*jobs, strict=True)))
    group = []
    for job in jobs:
        f = length_fiber(*job)
        group.append(f)
        print(job[3], sum(r["affected"] for r in f["regions"]), flush=True)
    rows = [
        {"case": args.case, "chain": args.regions, "m": m, "fiber": f["index"], **r}
        for f in group
        for r in f["regions"]
    ]
    frame = pd.DataFrame(rows)
    file = args.output / "de-novo-length.csv"
    if file.exists():
        frame = pd.concat([pd.read_csv(file), frame])
    frame.to_csv(file, index=False, lineterminator="\n", float_format="%.6g")
    days = SIZES[args.size]["days"] or DAYS_PER_YEAR * args.years
    rate = {"m": m, "days": days, "attempt": 1}
    summary, *_ = validation_tables(
        args.case,
        args.years,
        rate,
        group,
        SIZES[args.size]["boot"],
        tag=f"length-{args.regions}-{args.size}",
    )
    summary = fiber_decision(summary, group, FIBER_SIZES[args.size]["min_fibers"])
    summary = {"chain": args.regions, **summary}
    summary["seconds"] = round(sum(f["seconds"] for f in group), 1)
    summary["reference_types_low"] = summary["reference_types_high"] = float("nan")
    summary["length_statement"] = ""
    if args.reference is not None:
        ref = last_row(args.reference, args.case, args.years)
        low, high = float(ref["types_low"]), float(ref["types_high"])
        summary["reference_types_low"], summary["reference_types_high"] = low, high
        apart = summary["types_high"] < low or summary["types_low"] > high
        summary["length_statement"] = (
            "chain length matters" if apart else "no difference seen"
        )
        bounds = (low, high, summary["types_low"], summary["types_high"])
        if not all(math.isfinite(v) for v in bounds):
            summary["length_statement"] = "no statement (an interval is missing)"
        least = FIBER_SIZES[args.size]["min_fibers"]
        few = [int(ref.get("fibers_affected") or 0), summary["fibers_affected"]]
        if min(few) < least:
            summary["length_statement"] = "no statement (too few affected fibers)"
    print(json.dumps(summary), flush=True)
    frame = pd.DataFrame([summary])
    frame["m"] = [repr(float(m))]
    file = args.output / "de-novo-length-summary.csv"
    if file.exists():
        frame = pd.concat([pd.read_csv(file, dtype={"m": str}), frame])
    frame.to_csv(file, index=False, lineterminator="\n", float_format="%.6g")


def synthetic_sampler(scenario, rng, truth, slope=25.0, regions=32, icc=0.2):
    """Fiber outcomes from a known logistic curve, without simulation.

    ``W``: whole fibers are affected or not. ``C``: the regions of a fiber
    share a Beta-distributed probability with intra-fiber correlation ``icc``.
    """
    logit = math.log(TARGET / (1 - TARGET))

    def sample(stage, ms, fibers):
        out = []
        for m in ms:
            p = 1 / (1 + math.exp(-(logit + slope * math.log10(m / truth))))
            if scenario == "W":
                hits = np.where(rng.random(fibers) < p, regions, 0)
            else:
                s = (1 - icc) / icc
                hits = rng.binomial(regions, rng.beta(p * s, (1 - p) * s, fibers))
            out.append(([int(h) for h in hits], [regions] * fibers))
        return out

    return sample


def end_to_end(scenario, revision, truth, reps, row, fallback=None):
    """T6: whole synthetic calibrations, validation and refit included.

    With a list ``fallback``, every calibration that ends with a failed fit
    adds (true prevalence at the fallback point, affected fibers of a
    fallback cohort drawn there) to it (T7).
    """
    spec = FIBER_SIZES["full"]
    n3 = spec["fibers"][2]
    done = {1: 0, 2: 0}
    failed = covered = 0
    errors = []
    for rep in range(reps):
        rng = np.random.default_rng(seed_for("end-to-end", scenario, rep))
        sample = synthetic_sampler(scenario, rng, truth)
        first = truth * 10 ** rng.uniform(-0.03, 0.03)

        def judge(attempt, fit, sample=sample, rng=rng):
            hits, sizes = sample("validation", [fit["m"]], n3)[0]
            groups = [[1] * h + [0] * (s - h) for h, s in zip(hits, sizes, strict=True)]
            flat = [v for g in groups for v in g]
            _, low, high, _, _ = ratio_mean(flat, groups, spec["boot"], rng)
            some = sum(h > 0 for h in hits)
            return decide(low, high, some, spec["min_fibers"])

        points, _, attempts = fiber_calibration(
            first, sample, judge, spec, lambda _, rng=rng: rng, lambda _: None, revision
        )
        attempt, fit, label = attempts[-1]
        if label == "fit failed" and fallback is not None:
            # T7: the fallback cohort, drawn with its own seed.
            m = fallback_point(points)
            draw = np.random.default_rng(seed_for("fallback", scenario, rep))
            hits, _ = synthetic_sampler(scenario, draw, truth)("fallback", [m], n3)[0]
            slope, logit = 25.0, math.log(TARGET / (1 - TARGET))
            prob = 1 / (1 + math.exp(-(logit + slope * math.log10(m / truth))))
            fallback.append((prob, sum(h > 0 for h in hits)))
        if label == "compatible":
            done[attempt] += 1
            errors.append(abs(math.log10(fit["m"] / truth)))
        else:
            failed += 1
        if fit["fit_ok"]:
            covered += fit["m_boot_low95"] <= truth <= fit["m_boot_high95"]
    name = f"T6 end to end, {revision}"
    share = covered / reps
    row(name, scenario, "done at attempt 1", done[1])
    row(name, scenario, "done at attempt 2", done[2])
    row(name, scenario, "failed calibrations", failed)
    row(name, scenario, "coverage of true m, fiber bootstrap of the last fit", covered)
    median = 100 * (10 ** float(np.median(errors)) - 1) if errors else float("nan")
    row(name, scenario, "median error of m in done calibrations (%)", value=median)
    return share, failed


def coverage(args):
    """Coverage of the fitted-m intervals on synthetic clustered fibers."""
    import pandas as pd

    truth = 1e-4
    n1, n2, n3 = FIBER_SIZES["full"]["fibers"]
    resamples = FIBER_SIZES["full"]["boot"]
    least = FIBER_SIZES["full"]["min_fibers"]
    reps = args.repetitions
    rows = []

    def row(experiment, scenario, quantity, k=None, value=None, label=""):
        low = high = float("nan")
        if k is not None:
            value = k / reps
            low, high = wilson(k, reps)
        rows.append(
            {
                "experiment": experiment,
                "scenario": scenario,
                "quantity": quantity,
                "value": value,
                "low95": low,
                "high95": high,
                "repetitions": reps,
                "resamples": resamples,
                "label": label,
            }
        )

    start = time.perf_counter()
    fallbacks = {"W": [], "C": []}
    for scenario in ("W", "C"):
        boot = delta = failed_fit = failed_resamples = 0
        widths_boot, widths_delta = [], []
        for rep in range(reps):
            rng = np.random.default_rng(seed_for("coverage", scenario, rep))
            sample = synthetic_sampler(scenario, rng, truth)
            first = truth * 10 ** rng.uniform(-0.03, 0.03)
            points, _ = fiber_stages(first, sample, n1, n2, log=lambda _: None)
            fit = fiber_fit(points, resamples, rng)
            if not fit["fit_ok"]:
                failed_fit += 1
                continue
            boot += fit["m_boot_low95"] <= truth <= fit["m_boot_high95"]
            delta += fit["m_low95"] <= truth <= fit["m_high95"]
            failed_resamples += fit["boot_failed"]
            widths_boot.append((fit["m_boot_high95"] - fit["m_boot_low95"]) / truth)
            widths_delta.append((fit["m_high95"] - fit["m_low95"]) / truth)
        share = boot / reps
        label = ""
        if scenario == "W":
            label = "adequate" if 0.90 <= share <= 0.99 else "undercovers"
            if share > 0.99:
                label = "overcovers"
        name = "T4 fit interval"
        row(name, scenario, "coverage of true m, fiber bootstrap", boot, label=label)
        row(name, scenario, "coverage of true m, regions independent", delta)
        row(
            name,
            scenario,
            "median width / true m, fiber bootstrap",
            value=float(np.median(widths_boot)),
        )
        row(
            name,
            scenario,
            "median width / true m, regions independent",
            value=float(np.median(widths_delta)),
        )
        row(
            name,
            scenario,
            "failed point fits (no positive slope or inverse outside the fit points)",
            value=failed_fit,
        )
        row(
            name,
            scenario,
            "share of failed bootstrap refits",
            value=failed_resamples / (reps * resamples),
        )
        contains = enough = compatible = 0
        for rep in range(reps):
            rng = np.random.default_rng(seed_for("validation-rule", scenario, rep))
            sample = synthetic_sampler(scenario, rng, truth)
            hits, sizes = sample("validation", [truth], n3)[0]
            groups = [[1] * h + [0] * (s - h) for h, s in zip(hits, sizes, strict=True)]
            flat = [v for g in groups for v in g]
            _, low, high, _, _ = ratio_mean(flat, groups, resamples, rng)
            inside = low <= TARGET <= high
            many = sum(h > 0 for h in hits) >= least
            contains += inside
            enough += many
            compatible += inside and many
        name = "T5 validation rule at the true m"
        row(name, scenario, "fiber interval of prevalence contains 0.10", contains)
        row(name, scenario, f"at least {least} affected fibers", enough)
        row(name, scenario, "compatible", compatible)
        r0 = end_to_end(scenario, "R0", truth, reps, row, fallbacks[scenario])
        r1 = end_to_end(scenario, "R1", truth, reps, row)
        if scenario == "W":
            adopt = r1[0] >= 0.90 and r1[1] <= r0[1]
            row(
                "T6 decision",
                scenario,
                "R1 replaces R0 (last-fit coverage >= 0.90, no more failures)",
                value=int(adopt),
                label="R1 adopted" if adopt else "R0 kept",
            )
        print(f"{scenario} done after {time.perf_counter() - start:.0f} s", flush=True)
    # T7, appended so that the rows above keep their order: what the fallback
    # cohort of a calibration that ends with a failed fit is worth (R0).
    name = "T7 fallback cohort, R0"
    for scenario, found in fallbacks.items():
        row(
            name,
            scenario,
            "fit-failed calibrations with a fallback cohort",
            value=len(found),
        )
        prob = np.array([p for p, _ in found])
        for q, text in (
            (0.5, "median"),
            (0.025, "2.5% quantile"),
            (0.975, "97.5% quantile"),
        ):
            value = float(np.quantile(prob, q)) if len(prob) else float("nan")
            row(
                name,
                scenario,
                f"true prevalence at the fallback point, {text}",
                value=value,
            )
        enough = sum(k >= least for _, k in found)
        low, high = wilson(enough, len(found)) if found else (float("nan"),) * 2
        rows.append(
            {
                "experiment": name,
                "scenario": scenario,
                "quantity": f"share of fallback cohorts with at least {least} affected fibers",
                "value": enough / len(found) if found else float("nan"),
                "low95": low,
                "high95": high,
                "repetitions": reps,
                "resamples": resamples,
                "label": "",
            }
        )
    frame = pd.DataFrame(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        args.output / "review-fiber-fit-coverage.csv",
        index=False,
        lineterminator="\n",
        float_format="%.6g",
    )
    print(frame.to_string())


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
    kk14 = pd.read_csv(
        ROOT.parent / "kk14" / "results" / "reconstruction" / "comparison.csv"
    )
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(9.2, 4.8), layout="constrained", width_ratios=(1.1, 1)
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
    # B: mean types per affected region against age, last cohort of each
    # case and lifespan, with the 95% bootstrap interval over fibers.
    summary = pd.read_csv(args.results / "de-novo-summary.csv")
    summary = summary[summary.case.isin(COLORS) & summary.types_mean.notna()]
    last = summary.sort_values("attempt").groupby(["case", "years"]).last()
    shift = {"main": 1.0, "neutral": 0.9, "isolated": 1.1}
    for case, color in COLORS.items():
        if case not in last.index.get_level_values(0):
            continue
        part = last.loc[case].sort_index()
        x = part.index.to_numpy() * shift[case]
        whole = part.types_high > part.types_low
        right.plot(x, part.types_mean, color=color, linewidth=1.5, zorder=2)
        right.errorbar(
            x[whole],
            part.types_mean[whole],
            yerr=[
                (part.types_mean - part.types_low)[whole],
                (part.types_high - part.types_mean)[whole],
            ],
            fmt="o",
            color=color,
            ms=6,
            capsize=3,
            label=LABELS[case],
            zorder=3,
        )
        for age, row in part[~whole].iterrows():
            right.plot(
                age * shift[case],
                row.types_mean,
                "o",
                mfc="white",
                mec=color,
                mew=1.5,
                ms=7,
                zorder=3,
            )
            right.annotate(
                "one fiber,\ncalibration failed",
                (age * shift[case], row.types_mean),
                xytext=(-10, -8),
                textcoords="offset points",
                ha="right",
                va="top",
                fontsize=7.5,
                color=".3",
            )
    chapter17 = kk14[(kk14.label == "primary") & (kk14.cohort == "validation")]
    styles = {
        "drift": ("Neutral drift", ":"),
        "size": ("Shorter mtDNA", "-."),
        "transcription": ("Transcription", "--"),
    }
    for name, (mech, style) in styles.items():
        part = chapter17[chapter17.mechanism == name].sort_values("years")
        right.plot(
            part.years,
            part.mutant_types,
            style,
            color=".45",
            linewidth=1.2,
            label=f"{mech}, one cell (ch. 17)",
        )
    ages = sorted(set(last.index.get_level_values(1)) | set(chapter17.years))
    right.set(
        xscale="log",
        yscale="log",
        xlim=(2, 110),
        ylim=(0.8, 5000),
        xlabel="Age (years of 360 days)",
        ylabel="Mutant types per affected region",
        title="B  Types in affected regions",
    )
    right.set_xticks(ages, [str(a) for a in ages])
    right.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    right.legend(
        frameon=False,
        fontsize=7,
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.17),
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


def with_pool(workers, task, args):
    """Run ``task(args)``, with a process pool for the fibers if workers > 1."""
    global POOL
    if workers == 1:
        return task(args)
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=workers) as pool:
        POOL = pool
        try:
            return task(args)
        finally:
            POOL = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run", help="calibrate m, then simulate the validation cohort")
    p.add_argument("--size", choices=tuple(SIZES), default="quick")
    p.add_argument("--years", type=int, default=3)
    p.add_argument("--case", dest="cases", action="append", choices=tuple(CASES))
    p.add_argument("--work", type=Path, help="checkpoint folder")
    p.add_argument("--output", type=Path, help="folder for the result tables")
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        help="processes that simulate fibers at the same time; results do not depend on it",
    )
    p.add_argument(
        "--design",
        choices=("regions", "fibers"),
        default="regions",
        help="regions: the first design (default); fibers: fibers as sampling units",
    )
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
    p.add_argument("--m", type=float)
    p.add_argument("--m-from", type=Path, help="rates table: the m of its last attempt")
    p.add_argument("--regions", type=int, default=64)
    p.add_argument("--fibers", type=int, default=8)
    p.add_argument("--years", type=int, default=3)
    p.add_argument("--size", choices=tuple(SIZES), default="full")
    p.add_argument("--output", type=Path, default=ROOT / "runs" / "de-novo")
    p.add_argument("--work", type=Path, help="checkpoint folder (default: output)")
    p.add_argument(
        "--reference", type=Path, help="fiber-design summary to compare with"
    )
    p.add_argument("--workers", type=int, default=1)
    p = sub.add_parser("coverage", help="fit intervals on synthetic clustered fibers")
    p.add_argument("--repetitions", type=int, default=400)
    p.add_argument("--output", type=Path, default=ROOT / "runs" / "de-novo-fibers")
    args = parser.parse_args()
    if args.command == "length":
        if (args.m is None) == (args.m_from is None):
            parser.error("give one of --m and --m-from")
        if args.workers < 1:
            parser.error("--workers must be at least 1")
        args.output.mkdir(parents=True, exist_ok=True)
        with_pool(args.workers, sensitivity, args)
    elif args.command == "coverage":
        coverage(args)
    elif args.command == "figure":
        figure(args)
    elif args.command == "sampler":
        args.output.mkdir(parents=True, exist_ok=True)
        sampler_checks(args)
    elif args.command == "run":
        args.cases = args.cases or list(CASES)
        if args.workers < 1:
            parser.error("--workers must be at least 1")
        with_pool(args.workers, run, args)
    elif args.command == "pilot":
        args.cases = [(int(a), b, bool(c)) for a, b, c in (args.cases or [(8, 20, 1)])]
        pilot(args)


if __name__ == "__main__":
    main()
