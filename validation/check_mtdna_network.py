"""Check all enabled network hazards and the initially single-genotype law."""

import sys
from pathlib import Path

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "examples" / "mtdna_network")
)
from study import NAMES, build, predicted_slope  # noqa: E402


def main():
    ws, ms, wf, mf = initial = (2, 2, 1, 1)
    beta, gamma, mu = 0.3, 0.2, 0.1
    cases = {
        (2, 0, 1, 3): gamma,
        (0, 2, 3, 1): gamma,
        (1, 2, 2, 1): gamma * (wf * ws + mf * ws),
        (2, 1, 1, 2): gamma * (mf * ms + wf * ms),
        (1, 1, 2, 2): gamma * ws * ms,
        (3, 2, 0, 1): beta * wf,
        (2, 3, 1, 0): beta * mf,
        (1, 2, 1, 1): mu * ws,
        (2, 1, 1, 1): mu * ms,
        (1, 2, 3, 1): mu * ws,
        (2, 1, 1, 3): mu * ms,
        (2, 2, 2, 1): mu * wf,
        (2, 2, 1, 2): mu * mf,
    }
    params = {
        "mu": mu,
        "beta": beta,
        "gamma": gamma,
        "kappa": 6.0,
        "b": 0.0,
        "delta": 1.0,
    }
    dt = 1e-7
    result = (
        build(initial, **params)
        .fsp(dt, [dt], {k: 7 for k in NAMES}, tolerance=1e-14)
        .metadata
    )
    ids = [list(result["species"]).index(k) for k in NAMES]
    law = {
        tuple(row[i] for i in ids): p
        for row, p in zip(
            result["states"], result["snapshots"][0]["probabilities"], strict=True
        )
    }
    if max(abs(law[k] / dt - v) for k, v in cases.items()) > 2e-6:
        raise SystemExit("Mismatch in one-event fusion/fission/turnover rates")
    pure = build((3, 0, 1, 0), **params).fsp(0.5, [0.5], {k: 8 for k in NAMES}).metadata
    mutant = [i for i, label in enumerate(pure["species"]) if label.startswith("M")]
    if any(any(row[i] for i in mutant) for row in pure["states"]):
        raise SystemExit("Single-genotype model acquired a mutant")
    if abs(predicted_slope((326, 140, 373, 160), 0.023) - 4.51315475815e-6) > 1e-8:
        raise SystemExit("Wrong Equation 13 normalization")
    print("network hazards, no-mutation boundary and Equation 13 verified")


if __name__ == "__main__":
    main()
