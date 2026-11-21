"""Check spatial copy conservation, the exactly solvable migration law and the archive."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "examples" / "mtdna_spread")
)
from study import ROOT, build, paths, speeds  # noqa: E402


def main():
    archive = pd.read_csv(ROOT / "reference-1c.csv").query(
        "day == 700 and degradation == 'neutral'"
    )
    values = archive.set_index("condition")["mean"]
    targets = {"isolated": 0.503004, "equal": 0.499436, "dense": 0.540258}
    if any(abs(values[key] - expected) > 5e-7 for key, expected in targets.items()):
        raise SystemExit("Incorrect Figure 1C means")
    wave = speeds(pd.read_csv(ROOT / "reference-2e.csv"))
    if abs(wave.iloc[1].speed_um_day - 0.18973014373038566) > 1e-10:
        raise SystemExit("Incorrect Figure 2E reference")
    initial = [(0, 1), (0, 0)]
    case = build(initial, 0, 1, 0.2, death=0, c=0)
    result = case.fsp(2.0, [2.0], {"W0": 0, "W1": 0, "M0": 1, "M1": 1}).metadata
    i = list(result["species"]).index("M1")
    final = result["snapshots"][0]["probabilities"]
    p_right = sum(
        prob for st, prob in zip(result["states"], final, strict=True) if st[i]
    )
    if abs(p_right - (1 - np.exp(-0.8)) / 2) > 1e-10:
        raise SystemExit("Incorrect two-state migration law")
    x = paths([(7, 3), (3, 7)], (0, 1, 0.2), [0.0, 2.0], 64, 1936, 1, death=0, c=0)
    if np.any(x.sum(axis=2) != 10):
        raise SystemExit("Migration changed the combined genotype counts")
    print("spatial archive, migration law and genotype conservation verified")


if __name__ == "__main__":
    main()
