"""Check binomial partitioning, descendant counts and the reference data."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "examples" / "mtdna_bottleneck")
)
from model import birth_death_draw, build, forecast_sample, native_lineage  # noqa: E402
from study import ROOT  # noqa: E402


def main():
    rng = np.random.default_rng(2167)
    part = rng.binomial(2, 0.5, size=20000)
    empirical = np.bincount(part) / len(part)
    if np.max(np.abs(empirical - [0.25, 0.5, 0.25])) > 0.014:
        raise SystemExit("Wrong binomial law")
    model = build((4, 2), 0.15, 0.10)
    fsp = model.fsp(2.0, [2.0], {"W": 25, "M": 25}, max_states=700).metadata
    i = list(fsp["species"]).index("M")
    snapshot = fsp["snapshots"][0]
    result = np.bincount(
        np.asarray(fsp["states"])[:, i], weights=snapshot["probabilities"]
    )
    samples = birth_death_draw(np.full(10000, 2), 0.15, 0.1, 2.0, rng)
    if (
        np.max(
            np.abs(
                np.cumsum(np.bincount(samples, minlength=26) / len(samples))
                - np.cumsum(result)
            )
        )
        > 0.022
    ):
        raise SystemExit("Descendant sampler differs from the matrix law")
    paths = native_lineage(
        (4, 2), [(2.0, 0.15, 0.1, True), (2.0, 0.15, 0.1, True)], 9217
    )
    if paths.shape != (2, 2) or np.any(paths < 0):
        raise SystemExit("The native lineage driver failed")
    h = forecast_sample(80, 0.2, np.log(80), 0.0, 100, rng)
    if np.max(abs(h - 0.2)) > 1e-12:
        raise SystemExit("Deterministic ancestral count check failed")
    raw = pd.read_csv(ROOT / "reference-hb.csv")
    source = pd.read_csv(ROOT / "reference-table2.csv").set_index("group")
    stats = raw.groupby("group").heteroplasmy.agg(["mean", "var", "count"])
    if (stats["count"] != source.n).any() or (
        abs(stats["mean"] - source["mean"]) > 5e-5
    ).any():
        raise SystemExit("Re-extracted HB groups differ from the source table")
    print("division, descendant laws and HB reference groups verified")


if __name__ == "__main__":
    main()
