"""Check guide chapter 19: exact migration laws, the archive reader, the saved
results and the tables printed on the guide page.

Run with --write to print the two generated tables of the page.
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "mtdna_spread"
RESULTS = EXAMPLE / "results"
PAGE = REPO / "docs" / "guide" / "19-mtdna-spread.md"
sys.path.insert(0, str(EXAMPLE))
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


def region(page, name):
    start, end = f"<!-- {name} -->\n", f"<!-- /{name} -->"
    require(start in page and end in page, f"page region {name}")
    return page.partition(start)[2].partition(end)[0]


def main():
    if "--write" in sys.argv:
        print(figure1_table())
        print(figure2_table())
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
    print("spatial archive, migration laws, saved results and guide tables verified")


if __name__ == "__main__":
    main()
