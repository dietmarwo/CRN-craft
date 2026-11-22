"""Check tutorial S7: exact small laws, stored results and the page.

Offline and fast. It verifies the native FSP against a matrix exponential,
the generating-function marginal against the joint count law, the
digitized reference, the structure and the solver identities of the
stored sweep, the DKW bounds and gaps of the four ensembles, and that the
tables of the guide page equal those rendered from the stored files.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.linalg import expm_multiply
from scipy.stats import nbinom

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "two_stage_expression"
RESULTS = EXAMPLE / "results"
sys.path.insert(0, str(EXAMPLE))
from model import joint_generator, negative_binomial, protein_pmf  # noqa: E402
from study import native_check  # noqa: E402
from tables import TABLES  # noqa: E402


def require(condition, message):
    if not condition:
        raise SystemExit(f"check_two_stage_expression: {message}")


def exact_laws():
    check = native_check()
    require(max(check.values()) <= 1e-9, "independent and native FSP disagree")
    law = joint_generator(0.7, 1.2, 2.0, 12, 30)
    init = np.zeros(law.shape[0])
    init[0] = 1.0
    p = expm_multiply(law, init)
    approximate = protein_pmf(0.7, 1.2, 2.0, 1.0, 512)
    gap = np.max(np.abs(p[:-1].reshape(13, 31).sum(axis=0) - approximate[:31]))
    require(gap <= 1e-9, "protein PMF disagrees with joint count law")
    # Shape a, success probability 1/(1+b): mean a*b, variance a*b*(1+b).
    q = negative_binomial(0.5, 100, 20000)
    k = np.arange(q.size)
    require(abs(k @ q - 50) < 1e-6, "negative-binomial mean")
    require(abs((k - 50) ** 2 @ q / 50 - 101) < 1e-4, "negative-binomial Fano factor")
    require(np.allclose(q[:5], nbinom.pmf(np.arange(5), 0.5, 1 / 101)), "nbinom")


def stored_results():
    manifest = json.loads((EXAMPLE / "manifest.json").read_text())
    reference = pd.read_csv(EXAMPLE / "reference-figure1d.csv")
    require(len(reference) == 18, "Figure 1D reference has 18 markers")
    top = reference.sort_values("published_kl_display").iloc[-1]
    require(
        (top.a, top.b, round(top.published_kl_display, 5)) == (0.5, 100.0, 1.41713),
        "largest digitized Figure 1D value",
    )
    div = pd.read_csv(RESULTS / "divergence.csv")
    require(len(div) == 20 and bool(div.numerical_checks.all()), "sweep rows")
    fano = 1 + div.b * div.gamma / (1 + div.gamma)
    require(np.abs(div["mean"] / (div.a * div.b) - 1).max() < 1e-9, "stored means")
    require(np.abs(div.fano / fano - 1).max() < 1e-8, "stored Fano factors")
    require(div.tail_beyond_2047.max() < manifest["numerical_tail_target"], "tail")
    require(
        div.horizon_cdf_gap.max() < manifest["cdf_double_horizon_tolerance"],
        "horizon doubling",
    )
    review = pd.read_csv(RESULTS / "review-stationary-check.csv")
    require(len(review) == 20, "review rows")
    require(np.abs(review.kl_minus_shipped).max() < 1e-9, "two stationary solvers")
    hist = pd.read_csv(RESULTS / "histograms.csv")
    require(len(hist) == 4 * 152, "four histograms of 152 bins")
    simulated = div[div.n > 0]
    require(len(simulated) == 4 and (simulated.n == 10000).all(), "four ensembles")
    bound = np.sqrt(np.log(2 * 4 / 0.01) / (2 * 10000))
    require(np.allclose(simulated.dkw_bound, bound), "DKW bound")
    for (a, gamma), s in hist.groupby(["a", "gamma"]):
        s = s.sort_values("count")
        counts = s.ssa.to_numpy() * 10000
        require(np.abs(counts - np.round(counts)).max() < 1e-6, "integer counts")
        gap = np.abs(np.cumsum(s.ssa) - np.cumsum(s.reaction_pmf)).max()
        row = simulated[(simulated.a == a) & np.isclose(simulated.gamma, gamma)]
        require(abs(gap - row.cdf_gap.item()) < 1e-12 and gap <= bound, "CDF gap")


def page_tables():
    page = (
        ROOT / "docs" / "tutorials" / "simulation" / "S7-two-stage-expression.md"
    ).read_text(encoding="utf-8")
    for marker, render in TABLES.items():
        part = page.partition(f"<!-- {marker} -->\n")[2]
        part = part.partition(f"<!-- /{marker} -->")[0]
        require(part == render(), f"guide table {marker} is out of date")


def main():
    exact_laws()
    stored_results()
    page_tables()
    print("native FSP, matrix exponential and generating function verified")
    print("stored sweep, ensembles and guide tables verified")


if __name__ == "__main__":
    main()
