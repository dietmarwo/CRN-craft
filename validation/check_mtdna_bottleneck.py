"""Check the bottleneck laws, the reference data and the guide's result tables.

With --write the generated table regions of the guide pages are rewritten
from the stored result files; without it they must already be equal.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "mtdna_bottleneck"
RESULTS = EXAMPLE / "results"
PAGE = ROOT / "docs" / "tutorials" / "simulation" / "S6b-mtdna-bottleneck-results.md"
sys.path.insert(0, str(EXAMPLE))
from model import (  # noqa: E402
    birth_death_draw,
    build,
    forecast_sample,
    native_lineage,
    predictions,
    protocol_lineages,
)

THRESHOLD = 40.0


def require(condition, message):
    if not condition:
        raise SystemExit(f"check failed: {message}")


def rows(name):
    with (RESULTS / name).open(newline="") as f:
        return list(csv.DictReader(f))


def laws():
    rng = np.random.default_rng(2167)
    part = rng.binomial(2, 0.5, size=20000)
    empirical = np.bincount(part) / len(part)
    require(np.max(np.abs(empirical - [0.25, 0.5, 0.25])) <= 0.014, "binomial law")
    model = build((4, 2), 0.15, 0.10)
    fsp = model.fsp(2.0, [2.0], {"W": 25, "M": 25}, max_states=700).metadata
    i = list(fsp["species"]).index("M")
    result = np.bincount(
        np.asarray(fsp["states"])[:, i], weights=fsp["snapshots"][0]["probabilities"]
    )
    samples = birth_death_draw(np.full(10000, 2), 0.15, 0.1, 2.0, rng)
    gap = np.cumsum(np.bincount(samples, minlength=26) / len(samples))
    require(np.max(np.abs(gap - np.cumsum(result))) <= 0.022, "descendant law")
    paths = native_lineage(
        (4, 2), [(2.0, 0.15, 0.1, True), (2.0, 0.15, 0.1, True)], 9217
    )
    require(paths.shape == (2, 2) and np.all(paths >= 0), "native lineage driver")
    protocol_recurrence()
    h = forecast_sample(80, 0.2, np.log(80), 0.0, 100, rng)
    require(np.max(abs(h - 0.2)) <= 1e-12, "deterministic ancestral count")
    # Mean recurrence: growth per phase, halving per division (29 + 7).
    theta = np.zeros(16)
    theta[0] = 1e5
    theta[1:6] = [0.07, 0.05, 0.002, -0.001, 0.0]
    theta[7:12] = 0.03
    theta[13:16] = [10.0, 10.0, 10.0]
    base = math.log(1e5) + 0.07 * 203 + 0.05 * 112 - 36 * math.log(2)
    later = base + 24 * 10 * (0.002 - 0.001 + 0.0)
    got = predictions(theta, np.array([13.125, 43.125]))[0]
    require(abs(got[0, 0] - base) < 1e-9, "mean after the last division")
    require(abs(got[1, 0] - later) < 1e-9, "mean at the end of phase five")
    # First-order variance with no turnover is the sum of 1/(2 E) per division.
    flat = theta.copy()
    flat[1:6], flat[7:12] = 0, 0
    z = predictions(flat, np.array([7 / 24, 14 / 24]))[0][:, 1]
    require(abs(z[0] - 1 / 1e5) < 1e-15, "one division adds 1/(2 E)")
    require(abs(z[1] - 3 / 1e5) < 1e-15, "two divisions")


def protocol_recurrence():
    """The native scheduled lineage against the phase and division recurrences.

    Within a phase each copy founds an independent linear birth-death
    family; a division keeps each copy with probability 1/2, so
    E[X+] = E[X-] / 2 and Var(X+) = Var(X-) / 4 + E[X-] / 4. Every report
    mean and variance must lie within four standard errors (declared
    before the first run; 2,000 lineages, 8 comparisons per species).
    """
    phases = [
        (2.0, 0.30, 0.10, True),
        (1.5, 0.20, 0.20, False),
        (2.0, 0.05, 0.25, True),
        (1.0, 0.40, 0.10, True),
    ]
    initial = (40, 20)
    n = 2000
    result = protocol_lineages(initial, phases, 4051, replicates=n, threads=2)
    again = protocol_lineages(initial, phases, 4051, replicates=n, threads=1)
    paths = result.paths.sel(variable=["W", "M"]).values
    require(
        np.array_equal(paths, again.paths.sel(variable=["W", "M"]).values),
        "scheduled lineage depends on the worker count",
    )
    mean = np.array(initial, dtype=float)
    var = np.zeros(2)
    for k, (hours, birth, death, divide) in enumerate(phases):
        net = birth - death
        g = math.exp(net * hours)
        v = (birth + death) * hours if net == 0 else (birth + death) / net * g * (g - 1)
        mean, var = mean * g, var * g * g + mean * v
        if divide:
            mean, var = mean / 2, var / 4 + mean / 4
        x = paths[:, k, :]
        mu, s2 = x.mean(axis=0), x.var(axis=0, ddof=1)
        m4 = ((x - mu) ** 4).mean(axis=0)
        z_mean = (mu - mean) / np.sqrt(s2 / n)
        z_var = (s2 - var) / np.sqrt((m4 - s2 * s2) / n)
        require(
            np.all(np.abs(z_mean) < 4) and np.all(np.abs(z_var) < 4),
            f"scheduled lineage moments after phase {k + 1}",
        )


def references():
    t1 = pd.read_csv(EXAMPLE / "reference-table1.csv")
    require(len(t1) == 74, "Table 1 rows")
    counts = t1.groupby(["observable", "study"]).size().to_dict()
    require(
        counts
        == {
            ("copies", "Cao"): 14,
            ("copies", "Cree"): 10,
            ("copies", "Wai"): 10,
            ("normalized_variance", "Jenuth"): 16,
            ("normalized_variance", "Wai"): 24,
        },
        "Table 1 rows per study",
    )
    require((t1.day == 100).sum() == 5, "five mature-oocyte rows at 100 dpc")
    moved = sorted(set(t1.day[np.isclose(t1.day * 24 % 7, 0) & (t1.day < 1)]))
    require(moved == [0.0, 7 / 24, 14 / 24, 21 / 24], "division-timed ages")
    raw = pd.read_csv(EXAMPLE / "reference-hb.csv")
    source = pd.read_csv(EXAMPLE / "reference-table2.csv").set_index("group")
    require(len(raw) == 339 and len(source) == 15, "HB cells and groups")
    stats = raw.groupby("group").heteroplasmy.agg(["mean", "var", "count"])
    require((stats["count"] == source.n).all(), "HB group sizes")
    require((abs(stats["mean"] - source["mean"]) <= 5e-5).all(), "HB group means")
    ratio = stats["var"] / (stats["mean"] * (1 - stats["mean"]))
    require((abs(ratio - source.normalized_variance) <= 5e-5).all(), "HB V'")


def stored():
    abc = rows("abc.csv")
    require(len(abc) == 4, "four ensembles")
    for r in abc:
        require(int(r["sweeps"]) == 360000 and int(r["walkers"]) == 16, "run size")
        require(float(r["max_distance"]) <= THRESHOLD, "stored distances")
    fig = rows("figure3-comparison.csv")
    require(len(fig) == 12, "12 Figure 3 anchors")
    hb = rows("hb-comparison.csv")
    require(len(hb) == 15, "15 HB groups")
    traj = pd.read_csv(RESULTS / "trajectories.csv")
    require(len(traj) == 284 and traj.isna().sum().sum() == 0, "trajectory bands")
    require((traj.q025 <= traj.q975).all(), "band order")


def fmt(x, digits=3):
    return f"{float(x):.{digits}g}"


def headline():
    abc = pd.read_csv(RESULTS / "abc.csv")
    fig = pd.read_csv(RESULTS / "figure3-comparison.csv")
    hb = pd.read_csv(RESULTS / "hb-comparison.csv")
    vp = pd.read_csv(RESULTS / "review-vprime.csv")
    shift = pd.read_csv(RESULTS / "review-age-shift.csv").iloc[0]
    ens = pd.read_csv(RESULTS / "review-ensembles.csv")
    c = fig[fig.observable == "copies"]
    v = fig[fig.observable != "copies"]
    late = c[c.day == 100].iloc[0]
    outside = int(
        ((fig["median"] < fig.outer_low) | (fig["median"] > fig.outer_high)).sum()
    )
    draws = vp.groupby(["set", "draw"]).first()
    random = vp[(vp.set == "random") & vp.day.isin([10, 20, 35, 50, 70, 100])]
    hold = hb[hb.holdout]
    params = ens[ens.quantity.str.startswith("theta")]
    obs = ens[~ens.quantity.str.startswith("theta")]
    lines = [
        f"- Run: {len(abc)} ensembles x {int(abc.walkers[0])} walkers x "
        f"{int(abc.sweeps[0]):,} sweeps = {int((abc.sweeps * abc.walkers).sum()):,} "
        f"proposals; acceptance {100 * abc.acceptance.min():.2f}% to "
        f"{100 * abc.acceptance.max():.2f}%; largest stored distance "
        f"{abc.max_distance.max():.5f}.",
        f"- Effective size: sum over ensembles of the smallest of eight monitored "
        f"observables, {abc.ess_min_observable.sum():,.0f} (target 1,000).",
        f"- Figure 3: all {int(fig.center_in_95.sum())} of {len(fig)} digitized central "
        f"values inside the reconstructed 95% bands; median absolute difference "
        f"{fmt(c.difference.abs().median())} in natural-log copies and "
        f"{fmt(v.difference.abs().median())} in V'(h).",
        f"- Figure 3 discrepancy: at 100 dpc the reconstructed median is "
        f"{late['median']:,.0f} copies against {late.center:,.0f} digitized "
        f"(log ratio {late.difference:.2f}); reconstructed medians outside the "
        f"digitized outer envelope: {outside} of {len(fig)}.",
        f"- HB groups: {int(hb.in_sample_band.sum())} of {len(hb)} inside their "
        f"group-sized 95% bands; {int(hold.in_sample_band.sum())} of {len(hold)} "
        f"for the two flagged litters; {int(((hb.prediction_low <= hb.normalized_variance) & (hb.normalized_variance <= hb.prediction_high)).sum())} "
        f"of {len(hb)} inside the population band.",
        f"- First-order V'(h) against exact lineages: median relative error "
        f"{100 * random.relative_error.abs().median():.2f}%, largest "
        f"{100 * vp.relative_error.abs().max():.1f}%; "
        f"{int((draws.distance_exact <= THRESHOLD).sum())} of {len(draws)} draws "
        f"stay below 40, all {int((draws.distance_exact <= THRESHOLD + 2 * draws.distance_se).sum())} "
        f"within two Monte Carlo errors.",
        f"- Division timing: with the rounded labels 0.29, 0.58 and 0.88 dpc only "
        f"{100 * shift.accepted_with_rounded_labels:.1f}% of {int(shift.states):,} "
        f"stored states stay below 40.",
        f"- Ensemble agreement: largest R-hat {params.rhat_walkers.max():.2f} over "
        f"parameters and {obs.rhat_walkers.max():.2f} over the monitored observables.",
    ]
    return "\n".join(lines) + "\n"


def figure3_table():
    fig = pd.read_csv(RESULTS / "figure3-comparison.csv")
    out = [
        "| Day | Observable | Digitized center | Digitized outer envelope "
        "| Reconstructed median | Reconstructed 95% band |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for r in fig.itertuples():
        d = 0 if r.observable == "copies" else 4
        f = (lambda x: f"{x:,.0f}") if d == 0 else (lambda x: f"{x:.4f}")
        name = "copies" if d == 0 else "V'(h)"
        out.append(
            f"| {r.day:.0f} | {name} | {f(r.center)} | {f(r.outer_low)} to "
            f"{f(r.outer_high)} | {f(r.median)} | {f(r.q025)} to {f(r.q975)} |"
        )
    return "\n".join(out) + "\n"


def hb_table():
    hb = pd.read_csv(RESULTS / "hb-comparison.csv")
    out = [
        "| Group | Litter | Day | Cells | Observed V'(h) | Population band "
        "| Group-sized band | Inside |",
        "|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in hb.itertuples():
        flag = "*" if r.holdout else ""
        out.append(
            f"| {r.group} | {r.litter}{flag} | {r.day} | {r.n} | "
            f"{r.normalized_variance:.4f} | {r.prediction_low:.4f} to "
            f"{r.prediction_high:.4f} | {r.expected_sample_low:.4f} to "
            f"{r.expected_sample_high:.4f} | {'yes' if r.in_sample_band else 'no'} |"
        )
    return "\n".join(out) + "\n"


def phase6_table():
    d = pd.read_csv(RESULTS / "review-final-degradation.csv")
    out = [
        "| Phase-six degradation | Acceptance | Effective size | MAD log copies "
        "| MAD V'(h) | Copies at 100 dpc, median (95%) | HB inside |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in d.itertuples():
        label = (
            "fixed at 0"
            if r.phase6_degradation == "fixed_zero"
            else "free on [0, 1] per hour"
        )
        out.append(
            f"| {label} | {100 * r.acceptance_min:.1f}% to {100 * r.acceptance_max:.1f}% "
            f"| {r.ess_sum_of_minima:.0f} | {r.mad_log_copies:.3f} | {r.mad_vprime:.5f} "
            f"| {r.day100_copies_median:,.0f} ({r.day100_copies_q025:,.0f} to "
            f"{r.day100_copies_q975:,.0f}) | {r.hb_in_sample_band} of 15 |"
        )
    return "\n".join(out) + "\n"


def prose():
    """Numbers quoted in the running text of the chapter's pages."""
    fig = pd.read_csv(RESULTS / "figure3-comparison.csv")
    c = fig[fig.observable == "copies"].set_index("day")
    ratio = np.log(c.q975 / c.q025) / np.log(c.outer_high / c.outer_low)
    vp = pd.read_csv(RESULTS / "review-vprime.csv")
    low = vp[vp.copy_minimum].mean_copies
    draws = vp.groupby(["set", "draw"]).first()
    cross = draws[draws.distance_exact > THRESHOLD].iloc[0]
    ens = pd.read_csv(RESULTS / "review-ensembles.csv")
    sizes = ens.filter(like="ess_").dropna()
    hb = pd.read_csv(RESULTS / "hb-comparison.csv")
    late = c.loc[100.0]
    expected = {
        "S6-mtdna-bottleneck.md": [
            f"{math.exp(-late.difference):.0f} times below",
            f"{int(hb.in_sample_band.sum())} of their {len(hb)}",
        ],
        "S6b-mtdna-bottleneck-results.md": [
            f"{ratio.loc[10.0:50.0].min():.2f} to {ratio.loc[10.0:50.0].max():.2f}",
            f"about {ratio.loc[70.0:].mean():.0f} times",
            f"{low.min():,.0f} to {round(low.max(), -2):,.0f} copies",
            f"distance {cross.distance_first_order:.3f}",
            f"by {cross.distance_exact - THRESHOLD:.3f}",
            f"error ({cross.distance_se:.2f})",
            f"{100 * vp.fixed_or_lost_fraction.max():.0f}% of lineages",
            f"it is {sizes.sum(axis=1).min():,.0f}",
        ],
    }
    for page, items in expected.items():
        text = (PAGE.parent / page).read_text(encoding="utf-8")
        text = " ".join(text.split())
        for item in items:
            require(item in text, f"{page} should quote {item!r}")


REGIONS = {
    "mtdna-bottleneck-headline": headline,
    "mtdna-bottleneck-figure3": figure3_table,
    "mtdna-bottleneck-hb": hb_table,
    "mtdna-bottleneck-phase6": phase6_table,
}


def pages(write):
    text = PAGE.read_text(encoding="utf-8")
    for name, render in REGIONS.items():
        start, end = f"<!-- {name} -->\n", f"<!-- /{name} -->"
        head, found, rest = text.partition(start)
        require(found, f"region {name} missing")
        old, found, tail = rest.partition(end)
        require(found, f"region {name} not closed")
        new = render()
        if write:
            text = head + start + new + end + tail
        else:
            require(old == new, f"guide region {name} is out of date")
    if write:
        PAGE.write_text(text, encoding="utf-8", newline="\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true", help="rewrite the page regions")
    a = p.parse_args()
    laws()
    references()
    stored()
    pages(a.write)
    prose()
    print("bottleneck laws, references, stored results and guide tables verified")


if __name__ == "__main__":
    main()
