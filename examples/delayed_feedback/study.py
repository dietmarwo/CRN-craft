"""Delayed feedback benchmark using equations 15–17 of the cited paper.

The native run (a run/1.1 protocol with a fixed production delay) and an
independent delayed direct method written here in Python give the first
power-spectrum peak in the band 0.09–0.22 rad per time unit and the
autocorrelation at twice the delay. With the default 12 paths per sampler
the peak estimate varies between seeds by about 0.005–0.01 rad per time
unit, so the 0.03 tolerance (about ten frequency bins) detects only large
shifts of the oscillation period. The autocorrelation at twice the delay is
near 0.01 for both samplers, so its 0.10 tolerance cannot detect a
difference of that size. pi / delay is the angular frequency of a period of
twice the delay; it is a reference value, not a number read from the paper.
"""

from argparse import ArgumentParser
import heapq
import json
from pathlib import Path

import crn_craft as crn
import numpy as np
from scipy.signal import welch

# Figure 6A, main text. The forward dimerization constant multiplies
# the ordered-pair count, so a binomial reaction uses twice that constant.
A, B = 20.0, 4.0
K1, K1_REVERSE = 100.0, 1000.0
K2, K2_REVERSE = 200.0, 1000.0
DELAY = 20.0
DT, BURN, HORIZON = 0.5, 200.0, 2200.0
SEARCH = (0.09, 0.22)


def build():
    """Return the model of immediate reactions and fixed production delay."""
    model = crn.Model("delayed_feedback")
    model.add_species({"X": 0, "X2": 0, "D0": 1, "D1": 0})
    model.add_reaction(2 * K2, ["X", "X"], ["X2"], name="dimerize")
    model.add_reaction(K2_REVERSE, ["X2"], ["X", "X"], name="dissociate")
    model.add_reaction(K1, ["D0", "X2"], ["D1"], name="bind")
    model.add_reaction(K1_REVERSE, ["D1"], ["D0", "X2"], name="unbind")
    model.add_reaction(B, ["X"], [], name="decay")
    model.add_reaction(A, ["D0"], ["D0"], name="initiate")
    protocol = crn.Protocol().delay("initiate", DELAY, products={"X": 1})
    return model, protocol


def reference(times, seed):
    """Independent direct method for the same delayed-reaction process."""
    np.random.seed(seed)
    sites = np.array([0, 0, 1, 0], dtype=np.int64)
    delta = np.array(
        [
            [-2, 1, 0, 0],
            [2, -1, 0, 0],
            [0, -1, -1, 1],
            [0, 1, 1, -1],
            [-1, 0, 0, 0],
            [0, 0, 0, 0],
        ],
        dtype=np.int64,
    )
    pending = [0.0]
    pending.pop()
    samples = np.zeros(len(times))
    tick = 0
    at = 0.0
    while tick < len(times):
        x, d, o, c = sites
        hazards = (
            K2 * x * (x - 1),
            K2_REVERSE * d,
            K1 * o * d,
            K1_REVERSE * c,
            B * x,
            A * o,
        )
        total = sum(hazards)
        jump = at + np.random.exponential(1 / total) if total else np.inf
        due = pending[0] if pending else np.inf
        change = min(jump, due)
        while tick < len(times) and times[tick] < change:
            samples[tick] = x
            tick += 1
        if tick == len(times):
            break
        if due <= jump:
            heapq.heappop(pending)
            sites[0] += 1
        else:
            target = np.random.random() * total
            channel = 0
            while channel < 5 and target >= hazards[channel]:
                target -= hazards[channel]
                channel += 1
            sites += delta[channel]
            if channel == 5:
                heapq.heappush(pending, change + DELAY)
        at = change
    return samples


def spectrum(paths):
    """Average the replicate spectra before choosing the first band peak."""
    rates, psd = welch(paths, fs=1 / DT, nperseg=min(4096, paths.shape[1]), axis=1)
    angular = 2 * np.pi * rates
    band = (angular >= SEARCH[0]) & (angular <= SEARCH[1])
    aggregate = psd.mean(axis=0)
    peak = float(angular[band][np.argmax(aggregate[band])])
    centered = paths - paths.mean(axis=1, keepdims=True)
    variance = (centered**2).mean()
    lag = round(2 * DELAY / DT)
    corr = float((centered[:, :-lag] * centered[:, lag:]).mean() / variance)
    return {
        "peak_angular": peak,
        "acf_at_twice_delay": corr,
        "frequency_bin_width": float(angular[1]),
        "mean": float(paths.mean()),
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--paths", type=int, default=12)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    from numba import njit

    times = np.arange(BURN, HORIZON + DT / 2, DT)
    model, protocol = build()
    run = crn.ensemble(
        model,
        end=HORIZON,
        report_times=times,
        protocol=protocol,
        quantiles=(),
        replicates=args.paths,
        keep_paths=args.paths,
        max_output_bytes=256 << 20,
        seed=640110,
        threads=args.threads,
        strict=True,
    )
    actual = spectrum(run.simulate().paths.sel(variable="X").values)
    oracle = njit(cache=False)(reference)
    expected = spectrum(
        np.stack([oracle(times, 641110 + i) for i in range(args.paths)])
    )
    # The angular frequency of a period of twice the delay.
    reference_peak = np.pi / DELAY
    result = {
        "native": actual,
        "independent": expected,
        "pi_over_delay": reference_peak,
        "phase_tolerance": 0.03,
        "acf_tolerance": 0.10,
        "paths": args.paths,
        "native_vs_reference_peak": abs(
            actual["peak_angular"] - expected["peak_angular"]
        ),
        "native_vs_pi_over_delay": abs(actual["peak_angular"] - reference_peak),
        "acf_difference": abs(
            actual["acf_at_twice_delay"] - expected["acf_at_twice_delay"]
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    if (
        result["native_vs_reference_peak"] > 0.03
        or result["native_vs_pi_over_delay"] > 0.03
        or result["acf_difference"] > 0.10
    ):
        raise SystemExit("spectral comparison failed its declared tolerance")


if __name__ == "__main__":
    main()
