#!/usr/bin/env python3
"""Validate the cross-sectional inference tutorial against independent laws."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import Bounds

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples" / "kk26"))

import crn_craft as crn  # noqa: E402
from evidence import PAPER_MODES  # noqa: E402
from fitting import (  # noqa: E402
    STORED,
    TOLERANCE,
    code_hash,
    fit,
    fit_summary,
    read_csv,
)
from model import (  # noqa: E402
    BUILD_FLAGS,
    Model,
    cohort,
    forward,
    native_model,
    passage,
    tails,
)
from render import recovery_table  # noqa: E402

GUIDE = ROOT / "docs" / "guide" / "18b-mtdna-inference-results.md"

N_LAW = 16000
ALPHA = 0.001


def require(ok, reason):
    if not ok:
        raise AssertionError(reason)


def small_chain(model):
    """Construct Eq. 7 independently from per-parent mutation and replacement."""
    states = [
        (model.copies - u - a, u, a)
        for u in range(model.copies + 1)
        for a in range(model.copies - u + 1)
    ]
    index = {state: i for i, state in enumerate(states)}
    matrix = np.zeros((len(states), len(states)))
    mutate = np.array(
        [
            [
                1 - model.mutation,
                model.mutation * (1 - model.advantageous),
                model.mutation * model.advantageous,
            ],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    for i, state in enumerate(states):
        f = np.array(state) * [1, 1, 1 + model.advantage]
        daughter = (f / f.sum()) @ mutate
        for born in range(3):
            for dead in range(3):
                p = daughter[born] * state[dead] / model.copies
                if p:
                    after = list(state)
                    after[born] += 1
                    after[dead] -= 1
                    matrix[i, index[tuple(after)]] += p
    return np.asarray(states), matrix


def exact_laws():
    model = Model(0.06, 0.30, 0.60, copies=7, half_days=4)
    ages = np.array([0.0, 6.0, 17.0, 36.0])
    states, jump = small_chain(model)
    require(np.max(np.abs(jump.sum(axis=1) - 1)) < 1e-15, "Eq. 7 rows")
    initial = np.zeros(len(states))
    initial[0] = 1
    ticks = np.floor(model.rate * ages).astype(int)
    expected = np.array([initial @ np.linalg.matrix_power(jump, t) for t in ticks])
    saved, final = forward(model, ages, states=True)
    gap = float(np.max(np.abs(expected[-1] - final)))
    require(gap <= 1e-11, "full joint discrete law")
    require(np.max(np.abs(saved[:, 5:])) <= 1e-11, "complete lattice mass")
    total_mutant = states[:, 1] + states[:, 2]
    exact_tails = np.stack(
        [
            expected[:, 100 * total_mutant > threshold * model.copies].sum(axis=1)
            for threshold in (0, 5, 20, 50)
        ],
        axis=1,
    )
    require(
        np.max(np.abs(tails(saved[:, :5]) - exact_tails)) <= 1e-11,
        "strict threshold law",
    )

    capped = forward(model, ages, cap=1)
    require(capped[-1, 5] > 0, "test must produce leakage")
    lo = tails(capped[:, :5])
    require(np.all(lo <= exact_tails + 1e-11), "projection lower bound")
    require(
        np.all(lo + capped[:, 5, None] >= exact_tails - 1e-11), "projection upper bound"
    )

    native = native_model(model)
    generator = model.rate * (jump - np.eye(len(jump)))
    continuous = np.asarray([initial @ expm(t * generator) for t in ages])
    fsp = native.fsp(
        end=ages[-1],
        report_times=ages.tolist(),
        bounds={"W": 7, "U": 7, "A": 7},
        max_states=len(states),
        tolerance=1e-11,
    )
    native_order = {tuple(s): i for i, s in enumerate(fsp.metadata["states"])}
    indices = [native_order[tuple(s)] for s in states]
    held = np.array([s["probabilities"] for s in fsp.metadata["snapshots"]])[:, indices]
    fsp_gap = float(np.max(np.abs(held - continuous)))
    require(fsp_gap < 1e-10, "native FSP continuous law")

    run = crn.ensemble(
        native,
        end=ages[-1],
        report_times=[0.0, ages[-1]],
        replicates=N_LAW,
        threads=2,
        seed=718,
        method="static_direct",
        keep_paths="all",
        strict=True,
    ).simulate()
    end = run.paths.isel(time=-1)
    mutants = end.sel(variable=["U", "A"]).sum("variable").to_numpy().astype(int)
    require(
        np.all(end.sum("variable").to_numpy() == model.copies),
        "native fixed copy count",
    )
    empirical = np.bincount(mutants, minlength=model.copies + 1).cumsum() / N_LAW
    target = np.bincount(
        total_mutant, weights=continuous[-1], minlength=model.copies + 1
    ).cumsum()
    native_gap = float(np.max(np.abs(empirical - target)))
    bound = math.sqrt(math.log(2 / ALPHA) / (2 * N_LAW))
    require(native_gap < bound, "native Static Direct empirical CDF")

    counts, _ = cohort(model, N_LAW, 129903, ages, threads=3)
    geometric_gap = float(np.max(np.abs(tails(counts[-1]) / N_LAW - exact_tails[-1])))
    require(geometric_gap < bound, "geometric-skip cohort law")
    return {
        "discrete_joint_gap": round(gap, 13),
        "fsp_max_gap": round(fsp_gap, 13),
        "native_static_direct_cdf_gap": round(native_gap, 6),
        "discrete_geometric_cdf_gap": round(geometric_gap, 6),
        "law_bound": round(bound, 6),
    }


def clock_and_replay():
    model = Model(0.03, 0.2, 0.8, copies=22, half_days=8)
    ages = [0, 2, 5]
    seed = 2**64 - 3
    exact = cohort(model, 257, seed, ages, threads=1)[0]
    threads = cohort(model, 257, seed, ages, threads=4)[0]
    require(np.array_equal(exact, threads), "64-bit cohort seeds changed on replay")
    require(np.all(exact.sum(axis=1) == 257), "five bins do not partition cells")
    require(
        np.array_equal(tails(np.array([6, 5, 4, 3, 2])), [14, 9, 5, 2]), "tail indexing"
    )
    neutral = passage(
        Model(0, 0, 1, copies=10, half_days=3), max_days=1000, interval=10
    )
    require(abs(neutral["success"] - 1 / 9) < 1e-14, "neutral gambler's ruin")
    selected = passage(
        Model(0, 0.4, 1, copies=10, half_days=3), max_days=1000, interval=10
    )
    require(
        abs(selected["absorption"][-1, 1] - selected["success"]) < 1e-10,
        "computed success differs from geometric result",
    )
    require(np.all(np.diff(selected["cdf"]) >= -1e-11), "conditional CDF ordering")
    return {
        "thread_rows": len(ages),
        "neutral_target": neutral["target"],
        "selected_success": round(selected["success"], 8),
    }


def optimizer():
    # Tiny deterministic inverse check with known counts. The actual
    # bootstrap uses the grid fitted against the study-specific forward law.
    class Proxy:
        days = np.array([0, 30])
        model = Model(1e-5, 0.1, 0.8)
        bounds = Bounds([-6, 0.01, 0.1], [-4, 0.45, 1])

        def predict(self, x):
            p, s, f = x
            return np.array(
                [
                    [0, 0, 0, 0],
                    [
                        0.2 + (p + 5) * 0.08,
                        0.12 + (s - 0.1) * 0.10,
                        0.07 + (f - 0.8) * 0.07,
                        0.025,
                    ],
                ]
            )

    obs = np.array([[1000, 0, 0, 0, 0], [800, 80, 50, 45, 25]])
    r = fit(Proxy(), obs, 418405, 1100)
    require(r["ssq"] < 1e-10, "DE did not recover the target tails")
    json.dumps(r, allow_nan=False)
    return {"final_objective_below_1e_10": True}


def saved_evidence():
    data = read_csv(STORED / "fits.csv")
    counts = Counter(r["regime"] for r in data)
    require(counts == {"baseline": 40, "strong": 40}, "80 fits")
    manifest = json.loads((STORED / "manifest.json").read_text())
    require(
        manifest["moran_cpp_sha256"] == code_hash(), "saved solver source fingerprint"
    )
    require(
        len(read_csv(STORED / "fig6-differences.csv")) == 296, "two digitized panels"
    )
    grids = read_csv(STORED / "grid-checks.csv")
    kinds = Counter(r["kind"] for r in grids)
    require(
        kinds == {"off-grid": 6, "fitted": 6, "box": 32},
        "off-grid, whole-box and fitted interpolation checks",
    )
    require(
        all(float(r["max_curve_difference"]) <= TOLERANCE for r in grids),
        "grid accuracy",
    )
    check = [r for r in read_csv(STORED / "forward-checks.csv") if r["accepted"] != ""]
    require(
        len(check) == 9 and all(r["accepted"] == "True" for r in check),
        "full-size forward accuracy",
    )
    paper = read_csv(ROOT / "examples" / "kk26" / "paper-figure6.csv")
    compared = read_csv(STORED / "fig6-differences.csv")
    require(
        [r["value"] for r in paper] == [r["value"] for r in compared],
        "Figure 6 differences use the digitized values",
    )
    # The summary replays from the fit rows; numpy may sum in another order.
    summary = read_csv(STORED / "bootstrap.csv")
    replay = []
    for regime, grid in manifest["grids"].items():
        require(grid["build"] == " ".join(BUILD_FLAGS), f"{regime} grid build flags")
        axes = grid["axes"]
        bounds = Bounds([a[0] for a in axes], [a[-1] for a in axes])
        replay += fit_summary(
            regime, [r for r in data if r["regime"] == regime], bounds
        )
    require(len(replay) == len(summary) == 6, "six summary rows")
    for new, old in zip(replay, summary, strict=True):
        for key, value in new.items():
            if isinstance(value, float):
                ok = math.isclose(value, float(old[key]), rel_tol=1e-12, abs_tol=0)
            else:
                ok = str(value) == old[key]
            require(ok, f"bootstrap.csv {new['regime']} {new['parameter']} {key}")
    near_edge = sum(int(r["fits_within_1pct_of_box_edge"]) for r in summary)
    page = GUIDE.read_text()
    part = page.partition("<!-- kk26-recovery -->\n")[2].partition(
        "<!-- /kk26-recovery -->"
    )[0]
    require(part == recovery_table(STORED), "rendered guide table is out of date")
    rows = [
        r
        for r in read_csv(STORED / "takeover-summary.csv")
        if r["label"].startswith("s=")
    ]
    require(len(rows) == 4, "all figure 5 modes")
    for row in rows:
        require(
            int(row["mode_bin_lower_days"]) == PAPER_MODES[float(row["advantage"])],
            "figure 5 modal-bin agreement",
        )
    require(
        all(
            float(a["median_days"]) > float(b["median_days"])
            for a, b in zip(rows[:-1], rows[1:], strict=True)
        ),
        "takeover medians",
    )
    for name, grid in manifest["grids"].items():
        require(
            (STORED / "grids" / grid["file"]).is_file(),
            f"{name} stored forward cache missing",
        )
    return {
        "fits": len(data),
        "fits_within_1pct_of_box_edge": near_edge,
        "validated_interpolations": len(grids),
        "accepted_forward_checks": len(check),
        "paper_modes": len(rows),
    }


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--output", required=True, type=Path)
    args.add_argument("--laws-only", action="store_true")
    settings = args.parse_args()
    report = {"schema": "crn-craft.validation.kk26/0.1", "checks": {}, "failures": []}
    for test in (exact_laws, clock_and_replay, optimizer, saved_evidence):
        if settings.laws_only and test is saved_evidence:
            continue
        try:
            report["checks"][test.__name__] = test()
        except Exception as exc:
            report["failures"].append(f"{test.__name__}: {type(exc).__name__}: {exc}")
    settings.output.write_text(json.dumps(report, indent=2) + "\n")
    for msg in report["failures"]:
        print(msg, file=sys.stderr)
    print(f"kk26: {len(report['checks'])} checks, {len(report['failures'])} failures")
    return bool(report["failures"])


if __name__ == "__main__":
    sys.exit(main())
