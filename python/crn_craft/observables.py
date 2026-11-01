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
    x = x - x.mean()
    energy = float(x @ x)
    if energy == 0:
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
    """Saved-grid hitting times: interval-censored, not exact event-time FPT."""
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
    observed = hit.any(axis=1)
    index = hit.argmax(axis=1)
    upper = np.where(observed, t[index], t[-1])
    lower = np.where(observed, t[np.maximum(index - 1, 0)], t[-1])
    return {
        "samples": len(x),
        "censored": int((~observed).sum()),
        "censoring_probability": float((~observed).mean()),
        "restricted_mean_lower": float(lower.mean()),
        "restricted_mean_upper": float(upper.mean()),
        "horizon": float(t[-1]),
        "resolution": "interval_censored_saved_grid",
    }
