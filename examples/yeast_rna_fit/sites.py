"""Predict independent transcription-site intensities from frozen parameters."""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import irwinhall

from model import ROOT, hog1, matrix, switching

ELONGATION_TIME = 1630 / 63


def propagator(x, dose, times, occupancy, cap):
    if occupancy:
        a = np.zeros((4, 4))
        b = np.zeros_like(a)
        b[0, 1], b[1, 1] = 1.0, -1.0
        for src, dest, k in switching(x, 0):
            a[dest, src] += k
            a[src, src] -= k
    else:
        clone = x.copy()
        clone[11] = 0.0
        a, b = matrix(clone, cap)

    def rate(t):
        return max(0.0, x[1] - x[2] * float(hog1(t, dose, x[12])))

    def rhs(t, p):
        return a @ p + rate(t) * (b @ p)

    if occupancy:
        y = np.array([1.0, 0.0, 0.0, 0.0])
        p = solve_ivp(
            rhs,
            (0.0, max(times)),
            y,
            method="BDF",
            t_eval=times,
            rtol=1e-8,
            atol=1e-11,
            jac=lambda t, _, a=a, b=b, rate=rate: a + rate(t) * b,
            max_step=8,
        )
        if not p.success:
            raise RuntimeError(p.message)
        return p.y.T
    return a, b, rate, rhs


def intensity_map(cap, bins):
    """Each active polymerase adds an independent Uniform(0,1) intensity."""
    index = np.arange(2, bins + 1)
    out = np.zeros((cap + 1, len(index)))
    for k in range(1, cap + 1):
        cdf = irwinhall.cdf(index, k)
        out[k] = np.r_[np.diff(cdf), 1 - cdf[-1]]
    return out


def prediction(x, site_data):
    cap, bins = 80, 30
    mixing = intensity_map(cap, bins)
    rows, leakage = [], []
    for dose, values in site_data.groupby("dose"):
        minute = np.sort(values.minute.unique())
        start = np.maximum(0, 60 * minute - ELONGATION_TIME)
        state = propagator(x, dose, start, True, cap)
        a, b, rate, rhs = propagator(x, dose, [], False, cap)
        for t, left, promoter in zip(minute, start, state, strict=True):
            initial = np.zeros(a.shape[0])
            initial[:4] = promoter
            if t:
                result = solve_ivp(
                    rhs,
                    (left, 60 * t),
                    initial,
                    method="BDF",
                    rtol=2e-8,
                    atol=1e-11,
                    jac=lambda v, _, a=a, b=b, rate=rate: a + rate(v) * b,
                )
                if not result.success:
                    raise RuntimeError(result.message)
                initial = result.y[:, -1]
            leakage.append(max(0.0, 1 - initial.sum()))
            p = initial.reshape(-1, 4).sum(axis=1) @ mixing
            rows.extend((dose, t, k, prob) for k, prob in enumerate(p, 2))
            rows.append((dose, t, 0, 1 - p.sum()))
    return pd.DataFrame(
        rows, columns=["dose", "minute", "intensity_bin", "model_probability"]
    ), max(leakage)


def site_data():
    return pd.read_csv(ROOT / "source-sites.csv", dtype={"replicate": str})
