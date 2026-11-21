"""Native model and an independent protein count law from immigration of mRNAs."""

import crn_craft as crn
import numpy as np
from scipy.special import roots_legendre
from scipy.stats import nbinom


def build(a, b, gamma, initial=(0, 0)):
    model = crn.Model("two_stage_gene_expression")
    model.add_species({"M": int(initial[0]), "P": int(initial[1])})
    model.add_reaction(a, [], ["M"], name="transcription")
    model.add_reaction(b * gamma, ["M"], ["M", "P"], name="translation")
    model.add_reaction(gamma, ["M"], [], name="mrna_decay")
    model.add_reaction(1.0, ["P"], [], name="protein_decay")
    return model


def rules(a, b, gamma, m, p):
    return (
        (m + 1, p, a),
        (m, p + 1, b * gamma * m),
        (m - 1, p, gamma * m),
        (m, p - 1, p),
    )


def joint_generator(a, b, gamma, mmax, pmax):
    """Small independent finite state projection with a final absorbing sink."""
    from scipy.sparse import coo_matrix

    count = (mmax + 1) * (pmax + 1)
    rows, columns, entries = [], [], []
    for m in range(mmax + 1):
        for p in range(pmax + 1):
            origin = m * (pmax + 1) + p
            total = 0.0
            for mm, pp, rate in rules(a, b, gamma, m, p):
                if rate == 0:
                    continue
                total += rate
                dest = mm * (pmax + 1) + pp if mm <= mmax and pp <= pmax else count
                rows.append(dest)
                columns.append(origin)
                entries.append(rate)
            rows.append(origin)
            columns.append(origin)
            entries.append(-total)
    return coo_matrix((entries, (rows, columns)), shape=(count + 1,) * 2).tocsc()


def grid(limit, breaks, order):
    roots, weights = roots_legendre(order)
    bounds = np.array([0, *(x for x in breaks if x < limit), limit], dtype=float)
    mids = (bounds[:-1] + bounds[1:]) / 2
    scales = np.diff(bounds) / 2
    return (mids[:, None] + scales[:, None] * roots).ravel(), (
        scales[:, None] * weights
    ).ravel()


def protein_pmf(a, b, gamma, horizon, support=2048, order=10):
    """Evaluate the protein marginal via compound Poisson generating functions.

    RNA arrivals before the report form a Poisson process. An RNA of age s
    and lifetime l contributes a Poisson number with mean
    b*gamma*(exp(-(s-min(l,s)))-exp(-s)). Integrate over still living and
    dead RNAs separately, then invert on roots of unity. Support and
    quadrature convergence must be checked by the caller.
    """
    age, age_weights = grid(
        min(gamma * horizon, 36), (0.2, 0.75, 2, 4, 8, 15, 24), order
    )
    lifetime = age / gamma
    living_mean = -b * gamma * np.expm1(-lifetime)
    live_weight = a / gamma * age_weights * np.exp(-age)
    means = [living_mean]
    weights = [live_weight]
    for t, intensity, weight in zip(
        lifetime, living_mean, a * age_weights * np.exp(-age), strict=True
    ):
        end = min(36.0, horizon - t)
        lag, lag_w = grid(end, (0.2, 0.75, 2, 4, 8, 15, 24), order)
        means.append(intensity * np.exp(-lag))
        weights.append(weight * lag_w)
    means = np.concatenate(means)
    weights = np.concatenate(weights)
    z = np.expm1(-2j * np.pi * np.arange(support // 2 + 1) / support)
    logpgf = np.zeros_like(z)
    for start in range(0, means.size, 256):
        sl = slice(start, start + 256)
        logpgf += weights[sl] @ np.expm1(means[sl, None] * z)
    values = np.fft.irfft(np.exp(logpgf), n=support)
    if values.min() < -1e-9 or abs(values.sum() - 1) > 1e-10:
        raise RuntimeError(
            "Generating-function inversion failed probability accounting"
        )
    return np.maximum(0.0, values)


def negative_binomial(a, b, support):
    return nbinom.pmf(np.arange(support), a, 1 / (1 + b))


def moments(a, b, gamma):
    return a * b, 1 + b * gamma / (1 + gamma)


def comparison(p, q, cap=150):
    # Keep all unseen high counts in the tail, using the integrated support.
    pp = np.r_[p[: cap + 1], p[cap + 1 :].sum()]
    qq = np.r_[q[: cap + 1], max(0.0, 1 - q[: cap + 1].sum())]
    nz = pp > 0
    kl = np.dot(pp[nz], np.log(pp[nz] / qq[nz]))
    return float(kl), float(np.sum(np.abs(pp - qq)) / 2)
