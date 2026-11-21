"""Check chapter 20: network hazards, the reference reader and the saved results.

Every expected hazard below is a number worked out by hand from the
reaction list of the paper (Equations 1-10) and the authors' simulator
(hazards 0-16 of closed_loop_deg_single.c, xi = 0), not from the model
code under test.
"""

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "mtdna_network"
sys.path.insert(0, str(EXAMPLE))
from study import NAMES, build, predicted_slope, render  # noqa: E402

DT = 1e-7
GUIDE = ROOT / "docs" / "guide" / "20a-mtdna-network-results.md"


def require(ok, reason):
    if not ok:
        raise SystemExit(reason)


def one_event_rates(counts, params, bound=8):
    """Hazard of each one-event destination: P(state at DT) / DT."""
    meta = build(counts, **params).fsp(
        DT, [DT], dict.fromkeys(NAMES, bound), tolerance=1e-14
    )
    meta = meta.metadata
    ids = [list(meta["species"]).index(k) for k in NAMES]
    law = {
        tuple(row[i] for i in ids): p
        for row, p in zip(
            meta["states"], meta["snapshots"][0]["probabilities"], strict=True
        )
    }
    law.pop(tuple(counts), None)
    return {k: v / DT for k, v in law.items() if v > 1e-15}


def close(rates, expected, what):
    keys = set(rates) | set(expected)
    worst = max(
        (abs(rates.get(k, 0.0) - expected.get(k, 0.0)) for k in keys), default=0.0
    )
    require(worst < 1e-5, f"{what}: hazards differ by {worst:.3g}")


def fusion_identical():
    # Two identical singletons fuse at gamma per unordered pair:
    # gamma * k * (k - 1) / 2 with gamma = 0.4 gives 0, 0, 0.4, 1.2 for k = 0..3.
    # A factor-two error (gamma*k*(k-1) or gamma*k**2/2) gives 0.8 or 0.2 at
    # k = 1 or 2 and is caught.
    off = {"mu": 0.0, "beta": 0.0, "gamma": 0.4, "kappa": 0.0, "b": 0.0, "delta": 1.0}
    for k, rate in ((0, 0.0), (1, 0.0), (2, 0.4), (3, 1.2)):
        expected = {(k - 2, 0, 2, 0): rate} if rate else {}
        close(one_event_rates((k, 0, 0, 0), off), expected, f"WS+WS at {k}")
        mirrored = {(0, k - 2, 0, 2): rate} if rate else {}
        close(one_event_rates((0, k, 0, 0), off), mirrored, f"MS+MS at {k}")


def full_network():
    # State (WS, MS, WF, MF) = (2, 2, 1, 1); beta 0.3, gamma 0.2, mu 0.1,
    # b 0.01, kappa 10, delta 0.5. Replication rate per copy:
    # 0.1 + 0.01 * (10 - 3 - 0.5 * 3) = 0.155.
    params = {
        "mu": 0.1,
        "beta": 0.3,
        "gamma": 0.2,
        "kappa": 10.0,
        "b": 0.01,
        "delta": 0.5,
    }
    expected = {
        (0, 2, 3, 1): 0.2,  # WS+WS: 1 pair
        (2, 0, 1, 3): 0.2,  # MS+MS: 1 pair
        (1, 2, 2, 1): 0.8,  # WF+WS (1*2) and MF+WS (1*2), each 0.2 per pair
        (2, 1, 1, 2): 0.8,  # MF+MS and WF+MS
        (1, 1, 2, 2): 0.8,  # WS+MS: 4 pairs
        (3, 2, 0, 1): 0.3,  # WF fission
        (2, 3, 1, 0): 0.3,  # MF fission
        (1, 2, 1, 1): 0.2,  # WS degraded: 2 * 0.1
        (2, 1, 1, 1): 0.2,  # MS degraded
        (1, 2, 3, 1): 0.31,  # WS replicates into two WF: 2 * 0.155
        (2, 1, 1, 3): 0.31,  # MS replicates into two MF
        (2, 2, 2, 1): 0.155,  # WF replicates
        (2, 2, 1, 2): 0.155,  # MF replicates
    }
    close(one_event_rates((2, 2, 1, 1), params), expected, "full network")
    # Replication is clamped at zero: kappa = -100 makes the feedback negative.
    clamped = {**params, "kappa": -100.0}
    rates = one_event_rates((2, 2, 1, 1), clamped)
    births = [(1, 2, 3, 1), (2, 1, 1, 3), (2, 2, 2, 1), (2, 2, 1, 2)]
    require(
        all(k not in rates for k in births), "negative replication rate not clamped"
    )
    # Fused copies are never degraded (xi = 0): no state loses a fused copy alone.
    require(
        (2, 2, 0, 1) not in rates and (2, 2, 1, 0) not in rates, "fused copy degraded"
    )


def single_genotype():
    params = {
        "mu": 0.1,
        "beta": 0.3,
        "gamma": 0.2,
        "kappa": 6.0,
        "b": 0.01,
        "delta": 1.0,
    }
    pure = (
        build((3, 0, 1, 0), **params).fsp(0.5, [0.5], dict.fromkeys(NAMES, 8)).metadata
    )
    mutant = [i for i, label in enumerate(pure["species"]) if label.startswith("M")]
    require(
        not any(any(row[i] for i in mutant) for row in pure["states"]),
        "single-genotype model acquired a mutant",
    )


def equation13():
    # n = 999, h = 300/999, f_s = 466/999, mu = 0.023 per day:
    # 2 * 0.023 * (466/999) * (300/999) * (699/999) / 999 = 4.51315475815e-6.
    value = predicted_slope((326, 140, 373, 160), 0.023)
    require(abs(value / 4.51315475815e-6 - 1) < 1e-10, "wrong Equation 13 slope")


def saved_results():
    conditions = pd.read_csv(EXAMPLE / "conditions.csv")
    reference = pd.read_csv(EXAMPLE / "reference-variance.csv")
    slopes = pd.read_csv(EXAMPLE / "results" / "slopes.csv")
    manifest = json.loads((EXAMPLE / "manifest.json").read_text())
    require(len(conditions) == 39, "conditions.csv must have 39 rows")
    series = reference.groupby(["panel", "index"])
    require(series.ngroups == 39 and len(reference) == 1989, "archived series count")
    require(
        series.day.max().eq(500).all() and series.size().eq(51).all(), "archive grid"
    )
    first = reference[(reference.panel == "E") & (reference["index"] == 0)]
    require(
        abs(first.variance.iloc[1] - 3.36930317045e-05) < 1e-16, "archive value E:0"
    )
    require(
        list(slopes.condition)
        == list(conditions.panel + ":" + conditions["index"].astype(str)),
        "rows",
    )
    require(set(slopes.n) == {1024, 2000}, "path counts")
    require(set(manifest["primary"]) <= set(slopes.condition), "primary rows")
    for name in (
        "review-primary-intervals",
        "review-h65",
        "review-archive",
        "review-archive-traces",
    ):
        require((EXAMPLE / "results" / f"{name}.csv").exists(), f"missing {name}.csv")
    text = GUIDE.read_text()
    for region, table in render(EXAMPLE / "results").items():
        start, end = f"<!-- {region} -->", f"<!-- /{region} -->"
        require(start in text and end in text, f"missing region {region}")
        shown = text.split(start, 1)[1].split(end, 1)[0].strip()
        require(
            shown == table.strip(),
            f"guide table {region} differs from the saved results",
        )


def main():
    fusion_identical()
    full_network()
    single_genotype()
    equation13()
    saved_results()
    print(
        "network hazards, feedback, no-mutation boundary, Equation 13 and saved tables verified"
    )


if __name__ == "__main__":
    main()
