"""Review experiments: an independent stationary solver and Figure 1D conventions.

The stationary law of the four reactions is computed here a second way, as
the null vector of a truncated joint master equation (reflecting box,
sparse LU with one state pinned). It shares no code with ``model.py``
except the negative-binomial reference, which is written out again.
The second table recomputes the Figure 1D divergence under the conventions
the paper leaves open (direction, logarithm base, binning and tail,
estimation from a finite simulated sample) and sets each next to the
digitized published values. Targets and decision rules were declared
before the run; the guide page states them.

    python examples/two_stage_expression/review.py   # about two minutes
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parent
M = json.loads((ROOT / "manifest.json").read_text())
LN2 = np.log(2)


def stationary(a, b, gamma, mmax, pmax):
    """Protein and mRNA marginals of the stationary vector on a box."""
    nm, npp = mmax + 1, pmax + 1
    n = nm * npp
    # Number the smaller species fastest: the matrix is then banded.
    if nm < npp:
        p, m = np.divmod(np.arange(n), nm)

        def index(mm, pp):
            return pp * nm + mm
    else:
        m, p = np.divmod(np.arange(n), npp)

        def index(mm, pp):
            return mm * npp + pp

    rows, cols, vals = [], [], []
    for allowed, dm, dp, rate in (
        (m < mmax, 1, 0, np.full(n, float(a))),
        (p < pmax, 0, 1, b * gamma * m),
        (m > 0, -1, 0, gamma * m),
        (p > 0, 0, -1, p.astype(float)),
    ):
        src = np.nonzero(allowed & (rate > 0))[0]
        rows += [index(m[src] + dm, p[src] + dp), src]
        cols += [src, src]
        vals += [rate[src], -rate[src]]
    q = coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsc()
    pinned = index(min(round(a / gamma), mmax), min(int(a * b), pmax))
    keep = np.r_[0:pinned, pinned + 1 : n]
    x = np.empty(n)
    x[pinned] = 1.0
    x[keep] = spsolve(
        q[keep][:, keep].tocsc(),
        -q[keep][:, [pinned]].toarray().ravel(),
        permc_spec="NATURAL",
    )
    x = np.maximum(x, 0.0) / np.maximum(x, 0.0).sum()
    joint = x.reshape(npp, nm).T if nm < npp else x.reshape(nm, npp)
    return joint.sum(axis=0), joint.sum(axis=1)


def box(a, b, gamma):
    mrna = a / gamma
    mean, fano = a * b, 1 + b * gamma / (1 + gamma)
    mmax = int(mrna + 9 * np.sqrt(mrna) + 10)
    pmax = int(mean + 12 * np.sqrt(mean * fano) + 40)
    return mmax, max(pmax, 2600) if a < 1 else pmax


def negative_binomial(a, b, n):
    k = np.arange(n)
    p = 1 / (1 + b)
    return np.exp(
        gammaln(a + k) - gammaln(k + 1) - gammaln(a) + a * np.log(p) + k * np.log1p(-p)
    )


def binned(pmf, cap=150, width=1):
    pmf = np.pad(pmf, (0, max(0, cap + 2 - pmf.size)))
    body = np.add.reduceat(pmf[: cap + 1], np.arange(0, cap + 1, width))
    # The tail bin sums the computed tail and the mass beyond the support.
    return np.r_[body, pmf[cap + 1 :].sum() + max(0.0, 1 - pmf.sum())]


def kl(p, q):
    nz = p > 0
    with np.errstate(divide="ignore"):
        return float(p[nz] @ np.log(p[nz] / q[nz]))


def main():
    shipped = pd.read_csv(ROOT / "results" / "divergence.csv")
    reference = pd.read_csv(ROOT / "reference-figure1d.csv")
    rng = np.random.default_rng(M["seed"] + 7)
    checks, conventions = [], []
    for _, row in shipped.iterrows():
        a, b, gamma = row.a, row.b, row.gamma
        mmax, pmax = box(a, b, gamma)
        base, mrna = stationary(a, b, gamma, mmax, pmax)
        m2 = 2 * mmax if (2 * mmax + 1) * (2 * pmax + 1) <= 120_000 else mmax
        p2 = 2 * pmax if (m2 + 1) * (2 * pmax + 1) <= 150_000 else int(1.5 * pmax)
        pmf, _ = stationary(a, b, gamma, m2, p2)
        x = np.arange(pmf.size)
        mean = x @ pmf
        fano = (x - mean) ** 2 @ pmf / mean
        q = negative_binomial(a, b, pmf.size)
        pb, qb = binned(pmf), binned(q)
        checks.append(
            {
                "a": a,
                "b": b,
                "gamma": gamma,
                "m_bound": mmax,
                "p_bound": pmax,
                "edge_mass_m": float(mrna[-1]),
                "edge_mass_p": float(base[-1]),
                "larger_box_cdf_gap": float(
                    np.abs(np.cumsum(base) - np.cumsum(pmf)[: base.size]).max()
                ),
                "mean_relative_error": mean / (a * b) - 1,
                "fano_relative_error": fano / (1 + b * gamma / (1 + gamma)) - 1,
                "kl_nats": kl(pb, qb),
                "kl_minus_shipped": kl(pb, qb) - row.kl,
                "tail_beyond_150": float(pmf[151:].sum()),
                "tail_beyond_2047": float(pmf[2048:].sum()),
            }
        )
        published = reference[(reference.a == a) & np.isclose(reference.gamma, gamma)]
        if published.empty:
            continue
        candidate = {
            "a": a,
            "b": b,
            "gamma": gamma,
            "published": published.published_kl_display.item(),
            "nats_reaction_vs_nb": kl(pb, qb),
            "nats_nb_vs_reaction": kl(qb, pb),
            "bits_reaction_vs_nb": kl(pb, qb) / LN2,
            "nats_no_tail_bin": kl(pmf[:151], q[:151]),
            "nats_renormalized_0_150": kl(
                pmf[:151] / pmf[:151].sum(), q[:151] / q[:151].sum()
            ),
            "nats_full_support": kl(pmf, q),
            "nats_tail_above_100": kl(binned(pmf, 100), binned(q, 100)),
            "nats_tail_above_300": kl(binned(pmf, 300), binned(q, 300)),
            "nats_width5_bins": kl(binned(pmf, 149, 5), binned(q, 149, 5)),
        }
        for n in (1000, 10000):
            draws = [
                kl(rng.multinomial(n, pb / pb.sum()) / n, qb) / LN2 for _ in range(20)
            ]
            candidate[f"bits_sampled_{n}_mean"] = float(np.mean(draws))
            candidate[f"bits_sampled_{n}_sd"] = float(np.std(draws, ddof=1))
        conventions.append(candidate)
        print(f"{a:g},{b:g},{gamma:.4g}", flush=True)
    out = ROOT / "results"
    pd.DataFrame(checks).to_csv(out / "review-stationary-check.csv", index=False)
    pd.DataFrame(conventions).to_csv(
        out / "review-figure1d-conventions.csv", index=False
    )


if __name__ == "__main__":
    main()
