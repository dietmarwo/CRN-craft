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
from model import ROOT, matrix, native, nonspatial, quick_means  # noqa: E402
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
    p = nonspatial(x, 0.2, 5, 0.05, minutes=1)[1].ravel()
    if np.max(np.abs(p - law[:-1])) > 2e-5:
        raise SystemExit("Banded BDF2 does not converge to the native limit")
    predicted = quick_means(x, 0.2, 0.1, minutes=1)
    if not 0.1 < predicted[1, 0] < 1.0:
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
    n = observations().groupby("dose").cells.sum().tolist()
    if n != [25511, 40132]:
        raise SystemExit("Mature cell totals changed from the pinned fit archive")
    audit = pd.read_csv(ROOT / "source-audit.csv")
    if audit.total_disagrees.sum() != 235:
        raise SystemExit("Compartment accounting changed")
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
    print("yeast CME, count likelihood, cell totals and site intensity verified")


if __name__ == "__main__":
    main()
