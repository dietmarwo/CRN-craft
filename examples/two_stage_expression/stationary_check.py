"""Check the public stationary solver on the archived expression sweep."""

from argparse import ArgumentParser
import json
from pathlib import Path

import crn_craft as crn
import numpy as np
import pandas as pd

from model import build
from review import binned, kl, negative_binomial, stationary

ROOT = Path(__file__).resolve().parent


def sweep():
    inputs = pd.read_csv(ROOT / "results/review-stationary-check.csv")
    results = []
    for row in inputs.itertuples():
        m, p = int(row.m_bound), int(row.p_bound)
        model = build(row.a, row.b, row.gamma)
        anchor = (min(m, round(row.a / row.gamma)), min(p, round(row.a * row.b)))
        s = crn.stationary_fsp(
            model, {"M": m, "P": p}, return_state=anchor, max_states=120_000
        )
        pmf = s.marginal("P")
        oracle, _ = stationary(row.a, row.b, row.gamma, m, p)
        gap = float(np.max(abs(np.cumsum(pmf - oracle))))
        mean = float(pmf @ np.arange(len(pmf)))
        fano = float(pmf @ (np.arange(len(pmf)) - mean) ** 2 / mean)
        relative_mean = abs(mean / (row.a * row.b) - 1)
        relative_fano = abs(fano / (1 + row.b * row.gamma / (1 + row.gamma)) - 1)
        divergence = kl(binned(pmf), binned(negative_binomial(row.a, row.b, len(pmf))))
        results.append(
            {
                "a": row.a,
                "b": row.b,
                "gamma": row.gamma,
                "states": len(s.states),
                "cdf": gap,
                "mean_relative_error": relative_mean,
                "fano_relative_error": relative_fano,
                "kl_bits": divergence / np.log(2),
                "kl_gap_bits": abs(divergence - row.kl_nats) / np.log(2),
                "flux": s.outward_flux,
                "residual": s.residual_l1,
                "pass": gap < 0.001 and relative_mean < 0.01 and relative_fano < 0.01,
            }
        )
    return results


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = sweep()
    args.output.write_text(json.dumps(rows, indent=2) + "\n")
    if not all(r["pass"] for r in rows):
        raise SystemExit("stationary sweep exceeds a declared tolerance")


if __name__ == "__main__":
    main()
