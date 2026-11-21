"""Check guide chapter 19: exact migration laws, the archive reader, the saved
results and the tables printed on the guide pages, and the de novo extension:
exact bookkeeping laws of the lineage sampler, a moment law, its quick
calibration, and the saved tables of pages 19b and 19c.

Run with --write to print the generated tables of the pages.
"""

import csv
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "mtdna_spread"
RESULTS = EXAMPLE / "results"
PAGE = REPO / "docs" / "guide" / "19-mtdna-spread.md"
DE_NOVO_PAGE = REPO / "docs" / "guide" / "19c-mtdna-spread-de-novo-results.md"
KK14 = REPO / "examples" / "kk14" / "results" / "reconstruction" / "comparison.csv"
sys.path.insert(0, str(EXAMPLE))
import mutations  # noqa: E402
from study import build, paths  # noqa: E402


def require(condition, message):
    if not condition:
        raise SystemExit(f"check_mtdna_spread: {message}")


def rows(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file))


def crossing(x, h, level):
    i = [k for k in range(len(h) - 1) if h[k] >= level > h[k + 1]]
    require(len(i) == 1, f"one descending crossing at {level}")
    k = i[0]
    return x[k] + (x[k + 1] - x[k]) * (level - h[k]) / (h[k + 1] - h[k])


def archived_speeds():
    data = rows(EXAMPLE / "reference-2e.csv")
    require(len(data) == 279, "Figure 2E reference has 139 + 140 rows")
    out = {}
    for level in (0.4, 0.5, 0.6):
        ends = []
        for day in ("0", "1480"):
            part = [r for r in data if r["day"] == day]
            x = [float(r["position_um"]) for r in part]
            h = [float(r["mean"]) for r in part]
            ends.append(crossing(x, h, level))
        out[level] = (ends[0], ends[1], (ends[1] - ends[0]) / 1480)
    return out


def figure1_table():
    data = {
        (r["condition"], r["degradation"], float(r["day"])): r
        for r in rows(RESULTS / "two-regions.csv")
    }
    tests = {r["condition"]: r for r in rows(RESULTS / "terminal-checks.csv")}
    lines = [
        "| Condition at day 700 | Archived mean | Fresh mean | Fresh minus archive | Simultaneous 99% interval |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    names = {
        "isolated": "Isolated",
        "equal": "Equal density",
        "dense": "Higher mutant density",
    }
    for key, name in names.items():
        t = tests[key]
        lines.append(
            f"| {name}, neutral | {float(t['mean_reference']):.5f} | {float(t['mean']):.5f}"
            f" | {float(t['difference']):+.5f} | {float(t['low']):+.5f} to {float(t['high']):+.5f} |"
        )
    for key, name in {"deterministic": "Deterministic", **names}.items():
        r = data[(key, "extra", 700.0)]
        lines.append(
            f"| {name}, mutants degraded faster | {float(r['mean_reference']):.5f}"
            f" | {float(r['mean']):.5f} | {float(r['mean']) - float(r['mean_reference']):+.5f} | not tested |"
        )
    return "\n".join(lines) + "\n"


def figure2_table():
    lines = [
        "| Source of the profiles | Fibers | Half-height speed (µm/day) | 99% interval | Label |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    archive = archived_speeds()
    lines.append(
        f"| Authors' archive, days 0 and 1,480 | 800 | {archive[0.5][2]:.5f} | not available | reference |"
    )
    det = {
        r["variant"]: float(r["value"])
        for r in rows(RESULTS / "review-experiments.csv")
        if r["experiment"] == "T1-deterministic" and r["threshold"] == "0.5"
    }
    lines.append(
        f"| Same chain without noise (Nss = 3,500) | none | {det['ode Nss=3500']:.5f} | none | deterministic part |"
    )
    for r in rows(RESULTS / "approx-wave-speeds.csv"):
        if r["threshold"] != "0.5":
            continue
        low, high = float(r["low_99"]), float(r["high_99"])
        label = (
            "unresolved discrepancy"
            if high < float(r["archived_speed"])
            else "inconclusive"
        )
        lines.append(
            f"| Approximate, {float(r['step_day']):g}-day step, Nss = 3,500 | {r['fibers']}"
            f" | {float(r['speed_um_day']):.5f} | {low:.5f}–{high:.5f} | {label} |"
        )
    names = {
        "T3-step-0.0625": "Approximate, 0.0625-day step, Nss = 3,500",
        "T3b-step-0.5": "Approximate, 0.5-day step, Nss = 3,500, second cohort",
        "T2-copies-3000": "Approximate, 0.125-day step, Nss = 3,000",
    }
    review = rows(RESULTS / "review-experiments.csv")
    for key, name in names.items():
        r = next(
            x for x in review if x["experiment"] == key and x["threshold"] == "0.5"
        )
        lines.append(
            f"| {name} | {r['fibers']} | {float(r['value']):.5f}"
            f" | {float(r['low_99']):.5f}–{float(r['high_99']):.5f} | {r['label']} |"
        )
    return "\n".join(lines) + "\n"


CASE_NAMES = {
    "main": "Survival of the densest, 32 regions (this guide)",
    "neutral": "Equal density, 32 regions (this guide)",
    "isolated": "Survival of the densest, isolated regions (this guide)",
}
KK14_NAMES = {
    "drift": "Neutral drift, one cell (chapter 17)",
    "size": "Faster copying of shorter mtDNA, one cell (chapter 17)",
    "transcription": "Transcription-linked replication, one cell (chapter 17)",
}


def interval(r, mean, low, high, spec, scale=1, unit=""):
    """A mean with its fiber-bootstrap interval; a dash when no unit is
    affected, and no interval when all affected units lie in one fiber."""
    if r[mean] == "":
        return "–"
    m, lo, hi = (scale * float(r[k]) for k in (mean, low, high))
    if lo == hi:
        return f"{m:{spec}}{unit} (one fiber, no interval)"
    return f"{m:{spec}}{unit} ({lo:{spec}}–{hi:{spec}}{unit})"


def fibers_hit():
    """Per cohort: fibers with an affected region, and fibers affected whole."""
    hits = {}
    for x in rows(RESULTS / "de-novo-regions.csv"):
        key = (x["case"], x["years"], x["attempt"])
        per = hits.setdefault(key, {})
        per[x["fiber"]] = per.get(x["fiber"], 0) + (x["affected"] == "True")
    return {
        key: (sum(v > 0 for v in per.values()), sum(v == 32 for v in per.values()))
        for key, per in hits.items()
    }


def de_novo_table():
    """The comparison table of page 19c, in the layout of chapter 17b."""
    lines = [
        "| Age | Mechanism and unit | m per replication | Affected / units | Prevalence (95% Wilson) | Mutant types (95% bootstrap) | Target band |",
        "| ---: | --- | --- | ---: | --- | --- | --- |",
    ]
    summary = rows(RESULTS / "de-novo-summary.csv")
    ours = {r["years"] for r in summary}
    for r in summary:
        if r["case"] not in CASE_NAMES:
            continue
        cohort = "first fit" if r["attempt"] == "1" else "refit"
        lines.append(
            f"| {r['years']} y | {CASE_NAMES[r['case']]}, {cohort} | {float(r['m']):.4g}"
            f" | {r['affected']} / {r['regions']} regions"
            f" | {100 * float(r['prevalence']):.2f}% ({100 * float(r['wilson_low']):.2f}–{100 * float(r['wilson_high']):.2f}%)"
            f" | {interval(r, 'types_mean', 'types_low', 'types_high', '.0f')}"
            f" | {r['band']} |"
        )
    for r in rows(KK14):
        if (r["label"], r["cohort"]) != ("primary", "validation"):
            continue
        if r["years"] not in ours:
            continue
        lines.append(
            f"| {r['years']} y | {KK14_NAMES[r['mechanism']]} | {float(r['mutation']):.3g}"
            f" | {r['cox_negative']} / {r['classified']} cells"
            f" | {100 * float(r['prevalence']):.2f}% ({100 * float(r['wilson_low']):.2f}–{100 * float(r['wilson_high']):.2f}%)"
            f" | {float(r['mutant_types']):.2f} ({float(r['mutant_types_low']):.2f}–{float(r['mutant_types_high']):.2f})"
            f" | {r['target_band']} |"
        )
    return "\n".join(lines) + "\n"


def de_novo_units():
    """Second table of page 19c: zones, segments and the size of the types."""
    lines = [
        "| Case, cohort | Fiber 95% interval of prevalence | Fibers with affected regions (all 32 affected) | Mutant copies per type | Largest type, share of copies | Zones | Types per zone | Types per 32-region segment |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    names = {
        "main": "Denser mutant",
        "neutral": "Equal density",
        "isolated": "Isolated",
    }
    hit = fibers_hit()
    for r in rows(RESULTS / "de-novo-summary.csv"):
        if r["case"] not in names:
            continue
        cohort = "first fit" if r["attempt"] == "1" else "refit"
        some, whole = hit[(r["case"], r["years"], r["attempt"])]
        lines.append(
            f"| {names[r['case']]}, {r['years']} y, {cohort}"
            f" | {100 * float(r['fiber_boot_low']):.1f}–{100 * float(r['fiber_boot_high']):.1f}%"
            f" | {some} ({whole})"
            f" | {interval(r, 'copies_per_type', 'copies_per_type_low', 'copies_per_type_high', '.2f')}"
            f" | {interval(r, 'largest_share', 'largest_low', 'largest_high', '.2f', 100, '%')}"
            f" | {r['zones']}"
            f" | {interval(r, 'zone_types_mean', 'zone_types_low', 'zone_types_high', '.0f')}"
            f" | {interval(r, 'segment_types_mean', 'segment_types_low', 'segment_types_high', '.0f')} |"
        )
    return "\n".join(lines) + "\n"


def de_novo_laws():
    """Exact bookkeeping of the lineage sampler on small pure-Python runs."""
    small = {"nss": 20, "c": 0.01, "mu": 0.07, "gamma": 0.1, "delta": 0.2}
    for seed in (11, 12, 13):
        a = mutations.simulate(4, 0.05, 60.0, seed, compiled=False, **small)
        b = mutations.simulate(
            4, 0.05, 60.0, seed, compiled=False, labels=False, **small
        )
        # Labels do not change the copy-number path.
        require(np.array_equal(a["wild"], b["wild"]), "replay: wild type")
        require(np.array_equal(a["mutant"], b["mutant"]), "replay: mutants")
        require(np.array_equal(a["events"], b["events"]), "replay: events")
        births, deaths, _, made = a["events"]
        total = int(a["wild"].sum() + a["mutant"].sum())
        require(births - deaths == total - 4 * 20, "births minus deaths")
        for i, labels in enumerate(a["labels"]):
            require(len(labels) == a["mutant"][i], "one label per mutant copy")
            require(np.all((labels >= 1) & (labels <= made)), "labels are types")
        require(made >= 1, "the small run makes deletions")
    # Hops alone conserve every type: no birth or death when c = mu = 0.
    out = mutations.simulate(
        3, 0.5, 30.0, 5, nss=10, c=0.0, mu=0.0, gamma=0.3, delta=0.2, compiled=False
    )
    require(out["events"][0] == out["events"][1] == 0, "no birth or death")
    require(int(out["wild"].sum()) == 30, "hops conserve copies")
    # Moment law with c = 0: E[mutants] = L Nss (1 - exp(-m mu t)).
    sample = mutations.small_paths(
        600, 4242, regions=2, days=50.0, nss=10, c=0.0, m=0.1, compiled=False
    )
    expected = 20 * (1 - math.exp(-0.1 * 0.07 * 50))
    z = (sample[:, 5].mean() - expected) / (
        sample[:, 5].std(ddof=1) / math.sqrt(len(sample))
    )
    require(abs(z) < 3.29, f"mutant moment law, z = {z:.2f}")
    # The quick calibration runs end to end and brackets 10%.
    with tempfile.TemporaryDirectory() as work:
        rate, fit = mutations.calibrate("main", "quick", Path(work), 3, lambda _: None)
        share = [f["affected"] / f["regions"] for f in fit]
        require(min(share) <= 0.10 <= max(share), "quick fit brackets 10%")
        require(rate["m_low95"] < rate["m"] < rate["m_high95"], "quick fit interval")
    # Fibers simulated by a process pool equal the fibers simulated in turn.
    from concurrent.futures import ProcessPoolExecutor

    def cohort(work):
        group = mutations.cohort("main", "pool", 0.02, 4, "quick", Path(work), 3)
        return [{k: v for k, v in f.items() if k != "seconds"} for f in group]

    with tempfile.TemporaryDirectory() as work:
        serial = cohort(work)
    with tempfile.TemporaryDirectory() as work, ProcessPoolExecutor(2) as pool:
        mutations.POOL = pool
        try:
            pooled = cohort(work)
        finally:
            mutations.POOL = None
    require(serial == pooled, "fibers from a process pool equal fibers in turn")


def de_novo_saved():
    summary = rows(RESULTS / "de-novo-summary.csv")
    cohorts = [(r["years"], r["case"]) for r in summary]
    expected = [("3", c) for c in ["main", "main", "neutral", "neutral", "isolated"]]
    expected += [("3", "main/neutral")]
    expected += [(y, "main") for y in ("10", "40", "80") for _ in (1, 2)]
    require(cohorts == expected, "de novo cohorts")
    require(len(rows(RESULTS / "de-novo-rates.csv")) == 11, "de novo rates")
    regions = rows(RESULTS / "de-novo-regions.csv")
    require(len(regions) == 11 * 32 * 32, "de novo region rows")
    require(len(rows(RESULTS / "de-novo-segments.csv")) == 11 * 32, "segments")
    for r in summary:
        if r["case"] not in CASE_NAMES:
            continue
        # The refit rule: a second cohort only after one outside the band.
        if r["attempt"] == "2":
            first = next(
                x for x in summary if (x["case"], x["years"]) == (r["case"], r["years"])
            )
            require(first["band"] == "outside", "refit only after a miss")
        part = [
            x
            for x in regions
            if (x["case"], x["attempt"], x["years"])
            == (r["case"], r["attempt"], r["years"])
        ]
        hits = sum(x["affected"] == "True" for x in part)
        require(hits == int(r["affected"]), f"{r['case']} affected count")
        for x in part:
            w, m = int(x["wild"]), int(x["mutant"])
            require((5 * m >= 3 * (w + m)) == (x["affected"] == "True"), "60% rule")
        if hits == 0:
            require(r["types_mean"] == "", "no types mean without affected regions")
            continue
        mean = np.mean([int(x["types"]) for x in part if x["affected"] == "True"])
        require(abs(mean - float(r["types_mean"])) < 1e-3 * mean, "types mean")
    length = rows(RESULTS / "de-novo-length.csv")
    require(len(length) == 8 * 64, "chain-length rows")
    require(sum(r["affected"] == "True" for r in length) == 79, "64-region affected")
    checks = rows(RESULTS / "de-novo-sampler.csv")
    require(len(checks) == 8, "sampler rows")
    require(all(r["passed"] == "True" for r in checks), "sampler checks passed")


def region(page, name):
    start, end = f"<!-- {name} -->\n", f"<!-- /{name} -->"
    require(start in page and end in page, f"page region {name}")
    return page.partition(start)[2].partition(end)[0]


def main():
    if "--write" in sys.argv:
        print(figure1_table())
        print(figure2_table())
        print(de_novo_table())
        print(de_novo_units())
        return
    # Reference reader: the three day-700 neutral means of Figure 1C.
    reference = rows(EXAMPLE / "reference-1c.csv")
    require(len(reference) == 48, "Figure 1C reference has 48 rows")
    values = {
        r["condition"]: float(r["mean"])
        for r in reference
        if r["day"] == "700.0" and r["degradation"] == "neutral"
    }
    targets = {"isolated": 0.503004, "equal": 0.499436, "dense": 0.540258}
    for key, expected in targets.items():
        require(abs(values[key] - expected) <= 5e-7, f"Figure 1C {key} mean")
    # Archived wave: crossings and speed to the printed precision.
    speed = archived_speeds()
    require(round(speed[0.5][0], 2) == 1695.52, "day-0 half-height crossing")
    require(round(speed[0.5][1], 2) == 1976.33, "day-1,480 half-height crossing")
    require(abs(speed[0.5][2] - 0.18973) < 5e-6, "archived speed")
    require(round(speed[0.4][2], 4) == 0.2035, "archived speed at 0.4")
    require(round(speed[0.6][2], 4) == 0.1747, "archived speed at 0.6")
    # Exact law: one copy hopping between two sites at rate 0.2 each way.
    case = build([(0, 1), (0, 0)], 0, 1, 0.2, death=0, c=0)
    result = case.fsp(2.0, [2.0], {"W0": 0, "W1": 0, "M0": 1, "M1": 1}).metadata
    i = list(result["species"]).index("M1")
    final = result["snapshots"][0]["probabilities"]
    right = sum(p for s, p in zip(result["states"], final, strict=True) if s[i])
    require(abs(right - (1 - math.exp(-0.8)) / 2) < 1e-10, "two-state migration law")
    x = paths([(7, 3), (3, 7)], (0, 1, 0.2), [0.0, 2.0], 64, 1936, 1, death=0, c=0)
    require(np.all(x[..., 0].sum(axis=2) == 10), "migration conserves wild types")
    require(np.all(x[..., 1].sum(axis=2) == 10), "migration conserves mutants")
    # Saved results: structure.
    require(len(rows(RESULTS / "two-regions.csv")) == 48, "two-regions rows")
    require(len(rows(RESULTS / "terminal-checks.csv")) == 3, "terminal-check rows")
    require(
        len(rows(RESULTS / "approx-wave-speeds.csv")) == 9, "approximate speed rows"
    )
    require(len(rows(RESULTS / "approx-wave-profiles.csv")) == 2000, "profile rows")
    require(len(rows(RESULTS / "review-experiments.csv")) == 46, "review rows")
    # The guide page prints the saved results.
    page = PAGE.read_text(encoding="utf-8")
    require(
        region(page, "mtdna-spread-1c") == figure1_table(),
        "Figure 1C table out of date",
    )
    require(
        region(page, "mtdna-spread-2e") == figure2_table(),
        "Figure 2E table out of date",
    )
    de_novo_laws()
    de_novo_saved()
    page = DE_NOVO_PAGE.read_text(encoding="utf-8")
    require(
        region(page, "mtdna-spread-de-novo") == de_novo_table(),
        "de novo table out of date",
    )
    require(
        region(page, "mtdna-spread-de-novo-units") == de_novo_units(),
        "de novo unit table out of date",
    )
    print("spatial archive, migration laws, saved results and guide tables verified")
    print("de novo sampler laws, saved tables and page 19c verified")


if __name__ == "__main__":
    main()
