"""Check the constant-signal limit, both costs and processed observation counts."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.linalg import expm_multiply

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "examples" / "yeast_rna_fit")
)
from fitting import Data, observations  # noqa: E402
from model import ROOT, exact_means, matrix, native, nonspatial, quick_means  # noqa: E402
from render import TABLES, region  # noqa: E402
from sites import intensity_map, site_data  # noqa: E402


def main():
    x = np.array(
        [
            0.02,
            0.05,
            0.0,
            0.03,
            0.02,
            0.01,
            0.015,
            0.01,
            0.02,
            0.05,
            0.03,
            0.02,
            0.0,
            0.025,
            0.01,
        ]
    )
    a, b = matrix(x, 5, with_sink=True)
    init = np.zeros(a.shape[0])
    init[0] = 1
    law = expm_multiply((a + x[1] * b) * 60, init)
    model = native(x)
    result = model.fsp(
        60.0, [60.0], {"N": 5, "S1": 1, "S2": 1, "S3": 1, "S4": 1}
    ).metadata
    pos = {label: k for k, label in enumerate(result["species"])}
    ids = [
        4 * row[pos["N"]] + np.argmax([row[pos[f"S{k}"]] for k in range(1, 5)])
        for row in result["states"]
    ]
    if np.max(np.abs(law[ids] - result["snapshots"][0]["probabilities"])) > 5e-10:
        raise SystemExit("Native FSP and explicit forced model disagree")
    try:
        p = nonspatial(x, 0.2, 5, 0.05, minutes=1)[1].ravel()
    except RuntimeError as exc:
        # The fits need the example's C++ solver; the data, likelihood and
        # result checks below do not.
        print(f"C++ solver not checked: {exc}")
    else:
        if np.max(np.abs(p - law[:-1])) > 2e-5:
            raise SystemExit("Banded BDF2 does not converge to the native limit")
        predicted = quick_means(x, 0.2, 0.1, minutes=1)
        exact = exact_means(x, 0.2, [0, 1])[1, 0]
        if abs(predicted[1, 0] - exact) > 1e-4 * exact:
            raise SystemExit("Closed moment calculation failed")
    frame = pd.DataFrame(
        {
            "dose": [0.2] * 3,
            "replicate": ["a"] * 3,
            "minute": [1] * 3,
            "total": [0, 1, 2],
            "cells": [7, 2, 1],
        }
    )
    probabilities = [0.6, 0.3, 0.1]
    exact = -np.dot([7, 2, 1], np.log(probabilities))
    score = Data(frame).counts_cost({0.2: {1: np.array(probabilities)}})
    if abs(score - exact) > 1e-12:
        raise SystemExit("Likelihood does not preserve counts")
    rows = observations()
    n = rows.groupby("dose").cells.sum().tolist()
    if n != [25511, 40132]:
        raise SystemExit("Mature cell totals changed from the pinned fit archive")
    if Data(rows).n() != 59557 or rows[rows.minute > 0].cells.sum() != 59557:
        raise SystemExit("Usable cells after time zero changed")
    differ = rows[rows.total != rows.nuclear + rows.cytoplasm].cells.sum()
    audit = pd.read_csv(ROOT / "source-audit.csv", dtype={"replicate": str})
    if differ != 235 or audit.total_disagrees.sum() != 235:
        raise SystemExit("Compartment accounting changed")
    grouped = rows.groupby(["dose", "replicate", "minute"]).cells.sum()
    if (
        audit.cells.tolist()
        != grouped.loc[
            list(zip(audit.dose, audit.replicate, audit.minute, strict=True))
        ].tolist()
    ):
        raise SystemExit("Count rows and audit disagree")
    n_site = site_data().groupby("dose").cells.sum().tolist()
    if n_site != [25511, 39943]:
        raise SystemExit("Site archive no longer matches the reported cell totals")
    source = json.loads((ROOT / "source-parameters.json").read_text())
    if (
        abs(source["parameters"]["Spatial-Distributions"]["delay"] - 226.8641829723936)
        > 1e-8
    ):
        raise SystemExit("The spatial archive's forcing delay was lost")
    a, b = matrix(x, 3, 4, with_sink=True)
    init = np.zeros(a.shape[0])
    init[0] = 1.0
    law = expm_multiply((a + x[1] * b) * 20, init)
    native_joint = (
        native(x, spatial=True)
        .fsp(20.0, [20.0], {"N": 3, "C": 4, "S1": 1, "S2": 1, "S3": 1, "S4": 1})
        .metadata
    )
    pos = {label: k for k, label in enumerate(native_joint["species"])}
    ids = [
        4 * (4 * row[pos["C"]] + row[pos["N"]])
        + np.argmax([row[pos[f"S{k}"]] for k in range(1, 5)])
        for row in native_joint["states"]
    ]
    snap = native_joint["snapshots"][0]
    if (
        np.max(np.abs(law[ids] - snap["probabilities"])) > 5e-10
        or abs(law[-1] - snap["physical_leakage_estimate"]) > 5e-10
    ):
        raise SystemExit("Native joint FSP and the nuclear export generator disagree")
    intensities = intensity_map(8, 10)
    if np.any(intensities[:3]) or abs(intensities[3].sum() - 1 / 6) > 1e-12:
        raise SystemExit("Transcription-site active threshold mismatch")
    results()
    print("yeast CME, count likelihood, cell totals and site intensity verified")


def results():
    """Structure of the saved results and the tables of guide page 23b."""
    scores = pd.read_csv(ROOT / "results" / "scores.csv")
    if sorted(
        zip(scores.preset, scores.partition, scores.objective, strict=True)
    ) != sorted(
        [("full", "all", k) for k in ("mean", "histogram")]
        + [
            ("draft", p, k)
            for p in ("all", "0.2", "0.4")
            for k in ("mean", "histogram")
        ]
    ):
        raise SystemExit("scores.csv rows changed")
    if scores.coarse_refined_nats_per_cell.max() >= 0.01:
        raise SystemExit("A reported fit misses the refinement target")
    work = pd.read_csv(ROOT / "results" / "optimizer-work.csv")
    full = work[work.preset == "full"].groupby("objective").evaluations.sum()
    if full.to_dict() != {"histogram": 24175, "mean": 24174}:
        raise SystemExit("Full-run evaluation counts changed")
    start = pd.read_csv(ROOT / "results" / "review-optimizer-start.csv")
    best = start.groupby("box").change_vs_archived.min()
    if len(start) != 8 or not best["wide"] <= 0.01:
        raise SystemExit("Search-start review result changed")
    grid = pd.read_csv(ROOT / "results" / "spatial-fits.csv").set_index("objective")
    if grid.loc["histogram", ["log_i3_factor", "log_decay_factor"]].any() or set(
        grid.verified_selection
    ) != {"grid"}:
        raise SystemExit("The spatial histogram result is no longer the archive")
    page = (ROOT.parents[1] / "docs/guide/23b-yeast-rna-fit-results.md").read_text()
    for name, make in TABLES.items():
        if region(page, name) != make():
            raise SystemExit(f"Guide table {name} is out of date")


if __name__ == "__main__":
    main()
