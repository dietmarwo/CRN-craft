"""Python interface to the fixed-work, study-specific reference samplers."""

from __future__ import annotations

import ctypes as ct
import math
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

HERE = Path(__file__).resolve().parent
MECHANISMS = ("drift", "size", "transcription")
SAMPLER = "kk14-splitmix64/endpoints/1"
STATUS = ("horizon", "atp_collapse", "empty", "population_limit")


def size_capacity(population, half_life, wild_delay, initiation):
    delta = math.log(2) / half_life
    denominator = initiation * (2 * math.exp(-delta * wild_delay) - 1) / delta - 1
    if denominator <= 0:
        raise ValueError("size model requires positive equilibrium capacity")
    return population / denominator


def atp_equilibrium(population=1000, half_life=240):
    delta = math.log(2) / half_life
    atp = population * (0.1 - 0.01) / 0.2
    return delta, atp, atp * delta / (1 - delta)


@dataclass(frozen=True)
class Config:
    mechanism: str
    hours: int
    mutation: float
    population: int = 1000
    initial_mutants: int = 0
    half_life: float = 240.0
    wild_delay: int = 2
    mutant_delay: int = 1
    initiation: float = 1.0
    wild_only: bool = False
    step: float = 1.0

    def encode(self):
        if self.mechanism not in MECHANISMS:
            raise ValueError("unknown mechanism")
        if not 0 <= self.hours <= 5_000_000:
            raise ValueError("the horizon must be between 0 and 5,000,000 hours")
        if not 2 <= self.population <= 10_000:
            raise ValueError("the initial population must be between 2 and 10,000")
        if not 0 <= self.initial_mutants <= self.population:
            raise ValueError("initial mutants must fit in the initial population")
        if not 0 <= self.mutation <= 1 or not math.isfinite(self.half_life):
            raise ValueError("invalid mutation probability or half-life")
        if self.half_life <= math.log(2):
            raise ValueError("half-life must exceed ln 2 hours")
        if self.step not in (0.25, 0.5, 1.0):
            raise ValueError("step must be 1, 0.5 or 0.25 hour")
        if self.mechanism != "transcription" and self.step != 1.0:
            raise ValueError("sub-hour sensitivity is defined for transcription")
        if self.mechanism == "drift" and self.wild_only:
            raise ValueError("the drift comparison mutates daughters of all templates")
        if self.mechanism == "size":
            if not 1 <= self.mutant_delay <= self.wild_delay <= 48:
                raise ValueError("size delays must be ordered whole hours in [1, 48]")
            if not 0 < self.initiation <= 1:
                raise ValueError("hourly initiation multiplier must be in (0, 1]")
            size_capacity(
                self.population, self.half_life, self.wild_delay, self.initiation
            )
        return _Config(
            MECHANISMS.index(self.mechanism),
            self.population,
            self.initial_mutants,
            self.hours,
            self.wild_delay,
            self.mutant_delay,
            self.wild_only,
            self.half_life,
            self.mutation,
            self.initiation,
            self.step,
        )


class _Config(ct.Structure):
    _fields_ = [
        *(
            (name, ct.c_int)
            for name in (
                "mechanism",
                "population",
                "initial_mutants",
                "hours",
                "wild_delay",
                "mutant_delay",
                "wild_only",
            )
        ),
        *(
            (name, ct.c_double)
            for name in (
                "half_life",
                "mutation",
                "initiation",
                "step",
            )
        ),
    ]


class Endpoint(ct.Structure):
    _fields_ = [
        ("wild", ct.c_int),
        ("mutant", ct.c_int),
        ("clones", ct.c_int),
        ("largest", ct.c_int),
        ("status", ct.c_int),
        ("time", ct.c_double),
        ("steps", ct.c_uint64),
    ]


_U64 = np.ctypeslib.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS")
_POINTS = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="C_CONTIGUOUS")
_BASENAME = "kk14-reference"


@lru_cache(maxsize=1)
def library():
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    flags = ["-dynamiclib"] if sys.platform == "darwin" else ["-shared"]
    with TemporaryDirectory(prefix="kk14-cxx-") as tmp:
        output = Path(tmp) / (
            _BASENAME + (".dylib" if sys.platform == "darwin" else ".so")
        )
        try:
            subprocess.run(
                [
                    *compiler,
                    "-std=c++17",
                    "-O3",
                    # No fused multiply-add: compilers contract by default on
                    # targets with FMA (arm64), which changes the rounding.
                    "-ffp-contract=off",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-fPIC",
                    *flags,
                    str(HERE / "reference.cpp"),
                    "-o",
                    str(output),
                ],
                check=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError("install a C++17 compiler or set CXX") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"{' '.join(compiler)} could not build reference.cpp as C++17; "
                "set CXX to another compiler"
            ) from exc
        lib = ct.CDLL(str(output))
    lib.kk14_batch.argtypes = [
        ct.POINTER(_Config),
        ct.c_int,
        ct.c_int,
        _U64,
        ct.POINTER(Endpoint),
    ]
    lib.kk14_batch.restype = ct.c_int
    lib.kk14_delay.argtypes = [ct.c_double] * 5 + [ct.c_int, ct.c_int, _POINTS]
    lib.kk14_delay.restype = ct.c_int
    return lib


def simulate(cfg, seeds, *, labelled=True, threads=1, executor=None):
    """Return one structured row per seed in the input order."""
    encoded = cfg.encode()
    seeds = np.ascontiguousarray(seeds, dtype=np.uint64)
    if seeds.ndim != 1 or len(seeds) == 0 or threads < 1:
        raise ValueError("need a vector of seeds and at least one worker")
    lib = library()

    def batch(part):
        result = (Endpoint * len(part))()
        if lib.kk14_batch(ct.byref(encoded), labelled, len(part), part, result):
            raise RuntimeError(f"{cfg.mechanism} sampler failed inside the batch")
        return np.ctypeslib.as_array(result).copy()

    partitions = np.array_split(seeds, min(len(seeds), threads * 4))
    if threads == 1:
        return batch(seeds)
    if executor is not None:
        return np.concatenate(list(executor.map(batch, partitions)))
    with ThreadPoolExecutor(max_workers=threads) as pool:
        return np.concatenate(list(pool.map(batch, partitions)))


def solve_delay(
    step=0.25,
    half_life=240.0,
    wild_delay=2.0,
    mutant_size=0.5,
    initiation=1.0,
    days_per_year=365,
    max_years=200,
):
    """Method of steps with aligned delays and Heun integration."""
    tm = wild_delay * mutant_size
    size_capacity(1000, half_life, wild_delay, initiation)
    if not 0 < step <= tm <= wild_delay or step * math.log(2) / half_life > 0.1:
        raise ValueError("invalid method-of-steps mesh")
    if any(
        not math.isclose(round(t / step) * step, t, abs_tol=1e-12)
        for t in (tm, wild_delay)
    ):
        raise ValueError("both delays must be positive integer numbers of steps")
    n_steps = round(max_years * days_per_year * 24 / step)
    if not 1 <= n_steps <= 100_000_000:
        raise ValueError("method-of-steps work limit exceeded")
    stride = max(1, n_steps // 1600)
    points = np.zeros((2 + n_steps // stride, 5))
    written = library().kk14_delay(
        wild_delay,
        tm,
        half_life,
        initiation,
        step,
        n_steps,
        stride,
        points,
    )
    if written < 0 or not np.isfinite(points[:written]).all():
        raise RuntimeError("delay integration failed")
    return points[:written]


def solve_atp(rtol=1e-9, atol=1e-11):
    from scipy.integrate import solve_ivp

    delta, atp0, c = atp_equilibrium()

    def rhs(_hour, y):
        w, m, atp = y
        rate = c / (c + atp)
        return [
            (rate - delta) * w,
            (1.5 * rate - delta) * m,
            0.1 * w - 0.2 * atp - 0.01 * (w + m),
        ]

    def depleted(_hour, y):
        return y[2]

    depleted.terminal = True
    depleted.direction = -1
    # Clamp the trial-step size near zero to keep implicit trial stages in
    # the domain; the terminal event itself is located at zero.
    result = solve_ivp(
        rhs,
        (0, 24 * 400),
        [1000, 1, atp0],
        events=depleted,
        method="Radau",
        rtol=rtol,
        atol=atol,
        dense_output=True,
        max_step=12,
    )
    if not result.success or not len(result.t_events[0]):
        raise RuntimeError(
            f"ATP integral did not reach the depletion event: {result.message}"
        )
    times = np.linspace(0, result.t_events[0][0], 501)
    return np.column_stack([times, result.sol(times).T])
