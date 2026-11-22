"""Review checks of the nonspatial fits: search start and data choices.

Each command writes one ``review-*.csv`` file to ``results/``. The question,
runs and decision rule of ``start`` were fixed before it ran; page O3b
gives them. ``choices`` rescores saved parameter vectors and fits nothing.
"""

import argparse
import json
from time import perf_counter

import numpy as np
import pandas as pd
from fcmaes import decpp
from scipy.optimize import Bounds

from fitting import Data, nonspatial_laws, observations
from model import KEYS, ROOT, nonspatial, parameters

RESULTS = ROOT / "results"
ANCHORS = json.loads((ROOT / "source-parameters.json").read_text())["parameters"]
NARROW = np.log([3, 4, 4, 3, 3, 3, 3, 4, 3, 4, 4, 3, 1.6])
# Two decades each way for the 12 rates, a factor of two for the delay.
WIDE = np.r_[np.full(12, np.log(100.0)), np.log(2.0)]


def start(budget=2000, seeds=(7001, 7002, 7003, 7004)):
    """Histogram objective on 0.2 M from a box that does not favor the archive."""
    anchor = parameters(ANCHORS["NonSpatial-Distributions"])
    frame = observations()
    data = Data(frame[frame.dose == 0.2])
    n = data.n()

    def score(x):
        return data.counts_cost(nonspatial_laws(x, data, 250, 12))

    a0 = score(anchor)
    # The wide box is centered one random offset (up to one decade per
    # rate) away from the archived vector.
    offset = np.r_[np.random.default_rng(4242).uniform(-1, 1, 12) * np.log(10), 0]
    rows = []
    for name, center, half in (
        ("narrow", np.zeros(13), NARROW),
        ("wide", offset, WIDE),
    ):
        for seed in seeds:
            calls = invalid = 0
            t0 = perf_counter()

            def cost(z, center=center):
                nonlocal calls, invalid
                calls += 1
                x = anchor.copy()
                x[:13] *= np.exp(center + z)
                value = score(x)
                if not np.isfinite(value):
                    invalid += 1
                    return 1e18
                return value

            found = decpp.minimize(
                cost,
                bounds=Bounds(-half, half),
                popsize=26,
                max_evaluations=budget,
                rg=np.random.default_rng(seed),
            )
            log10 = (center + found.x) / np.log(10)
            rows.append(
                {
                    "box": name,
                    "seed": seed,
                    "requested": budget,
                    "evaluations": calls,
                    "invalid": invalid,
                    "cells": n,
                    "nll_per_cell": found.fun / n,
                    "archived_nll_per_cell": a0 / n,
                    "change_vs_archived": (found.fun - a0) / n,
                    "max_abs_log10_from_archived": float(np.max(np.abs(log10))),
                    "seconds": perf_counter() - t0,
                    **{
                        f"log10_{k}": float(v)
                        for k, v in zip(KEYS[:13], log10, strict=True)
                    },
                }
            )
            print(rows[-1]["box"], seed, rows[-1]["change_vs_archived"], flush=True)
    pd.DataFrame(rows).to_csv(RESULTS / "review-optimizer-start.csv", index=False)


def choices():
    """Score three vectors on the paper's and on this example's cell choices."""
    counts = observations()
    counts = counts[counts.minute > 0].assign(
        compartments=lambda f: f.nuclear + f.cytoplasm
    )
    fits = json.loads((RESULTS / "parameters.json").read_text())
    vectors = {
        "archived histogram fit": parameters(ANCHORS["NonSpatial-Distributions"]),
        "histogram fit": parameters(fits["histogram"]),
        "mean fit": parameters(fits["mean"]),
    }
    sets = {
        "all replicates, recorded total": (counts, "total"),
        "all replicates, nuclear plus cytoplasmic": (counts, "compartments"),
        "authors' fit replicates, nuclear plus cytoplasmic": (
            counts[counts.replicate != "20131022"],
            "compartments",
        ),
    }
    sets.update(
        {
            f"replicate {rep}, recorded total": (part, "total")
            for rep, part in counts.groupby("replicate")
        }
    )
    rows = []
    for name, x in vectors.items():
        laws = {dose: nonspatial(x, dose, 720, 1.5).sum(-1) for dose in (0.2, 0.4)}
        for label, (frame, column) in sets.items():
            loss = 0.0
            for dose, part in frame.groupby("dose"):
                p = laws[dose][part.minute.to_numpy(), part[column].to_numpy()]
                loss -= float(np.dot(part.cells, np.log(p)))
            rows.append(
                {
                    "vector": name,
                    "cells": label,
                    "n": int(frame.cells.sum()),
                    "nll_per_cell": loss / frame.cells.sum(),
                }
            )
    pd.DataFrame(rows).to_csv(RESULTS / "review-data-choices.csv", index=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("start", "choices"))
    {"start": start, "choices": choices}[p.parse_args().command]()


if __name__ == "__main__":
    main()
