"""Trajectory-wise diagnostics; exploratory oscillation evidence, not a proof."""

import numpy as np


def oscillation(times, values, *, burn_fraction=0.2):
    t = np.asarray(times, dtype=float)
    x = np.asarray(values, dtype=float)
    if (
        t.ndim != 1
        or x.shape != t.shape
        or len(t) < 16
        or not np.isfinite(t).all()
        or not np.isfinite(x).all()
    ):
        raise ValueError("finite one-dimensional samples, at least 16, required")
    dt = np.diff(t)
    if (dt <= 0).any() or not np.allclose(dt, dt[0], rtol=1e-10, atol=0):
        raise ValueError("spectrum requires an increasing uniform grid")
    if not 0 <= burn_fraction < 0.8:
        raise ValueError("invalid burn fraction")
    x = x[int(len(x) * burn_fraction) :]
    if len(x) < 16:
        raise ValueError("at least 16 samples must remain after the burn-in")
    scale = float(np.abs(x).max())
    x = x - x.mean()
    energy = float(x @ x)
    # A constant series leaves only rounding residue after removing the mean.
    if energy <= len(x) * (64 * np.finfo(float).eps * scale) ** 2:
        return {
            "classification": "nonoscillatory",
            "peak_fraction": 0.0,
            "autocorrelation_peak": 0.0,
            "period": None,
            "cycles": 0.0,
        }
    fft = np.fft.rfft(x)
    power = np.abs(fft) ** 2
    power[0] = 0
    index = int(power.argmax())
    if index == 0:
        return {
            "classification": "nonoscillatory",
            "peak_fraction": 0.0,
            "autocorrelation_peak": 0.0,
            "period": None,
            "cycles": 0.0,
        }
    fraction = float(power[index] / power.sum())
    period = len(x) * float(dt[0]) / index
    cycles = float(index)
    lag = int(round(period / dt[0]))
    ac = float(x[:-lag] @ x[lag:] / energy) if 0 < lag < len(x) else 0.0
    classification = (
        "oscillatory"
        if cycles >= 3 and fraction >= 0.3 and ac >= 0.2
        else ("ambiguous" if cycles < 3 else "nonoscillatory")
    )
    return {
        "classification": classification,
        "peak_fraction": fraction,
        "autocorrelation_peak": ac,
        "period": period,
        "cycles": cycles,
    }


def sampled_first_passage(times, paths, threshold):
    """First saved time at or above a threshold, with only the bounds it supports.

    A path can cross the threshold and fall back between saved times, so the
    saved grid can miss a passage and its earlier sample says nothing about
    when the passage happened. What the grid does give:

    - a detected path passed at or before its first saved hit, so that time is
      an upper bound on its passage time;
    - min(passage time, horizon) is at most that time, or the horizon when
      nothing was detected, so ``restricted_mean_upper`` bounds the restricted
      mean from above;
    - ``detected_fraction`` is a lower bound on P(passage <= horizon) and
      ``not_detected_fraction`` an upper bound on P(passage > horizon).

    There is no lower bound on the passage time from saved samples alone.
    """
    t = np.asarray(times, dtype=float)
    x = np.asarray(paths, dtype=float)
    if (
        t.ndim != 1
        or x.ndim != 2
        or x.shape[1] != len(t)
        or len(t) < 2
        or len(x) == 0
        or not np.isfinite(x).all()
        or not np.isfinite(t).all()
        or not np.isfinite(threshold)
        or t[0] != 0
        or (np.diff(t) <= 0).any()
    ):
        raise ValueError("finite paths on an increasing grid starting at zero required")
    hit = x >= threshold
    detected = hit.any(axis=1)
    upper = np.where(detected, t[hit.argmax(axis=1)], t[-1])
    return {
        "samples": len(x),
        "detected": int(detected.sum()),
        "not_detected": int((~detected).sum()),
        "detected_fraction": float(detected.mean()),
        "not_detected_fraction": float((~detected).mean()),
        "restricted_mean_upper": float(upper.mean()),
        "horizon": float(t[-1]),
        "resolution": "saved_grid_upper_bounds",
    }
