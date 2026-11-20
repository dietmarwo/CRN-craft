"""Fixed-population Moran law and the synthetic ATP comparison."""

from __future__ import annotations

import ctypes as ct
import math
import os
import platform
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
AGES = np.arange(0, 1081, 30, dtype=float)
THRESHOLDS = (0, 5, 20, 50)
DISCRETE = "kk26-discrete-moran/1"
SPLIT = "kk26-poisson-strang-cn/1"
SAMPLER = "kk26-discrete-geometric-splitmix64/1"
ATP = "kk26-hourly-atp-reconstruction/1"


@dataclass(frozen=True)
class Model:
    mutation: float
    advantage: float
    advantageous: float
    copies: int = 1000
    half_days: float = 10.0

    @property
    def rate(self):
        return self.copies * math.log(2) / self.half_days

    def encode(self, cap=None):
        if not 2 <= self.copies <= 1500 or not 0.1 <= self.half_days <= 1000:
            raise ValueError(
                "copies must be in [2, 1500] and half-life in [0.1, 1000] days"
            )
        if not 0 <= self.mutation <= 1 or not 0 <= self.advantageous <= 1:
            raise ValueError("mutation and advantageous fraction must be probabilities")
        if not 0 <= self.advantage <= 10:
            raise ValueError("selection advantage must be in [0, 10]")
        cap = self.copies if cap is None else cap
        if not 0 <= cap <= self.copies:
            raise ValueError("the neutral-count cap must be in [0, copies]")
        return _Model(
            self.copies,
            cap,
            self.mutation,
            self.advantage,
            self.advantageous,
            self.half_days,
        )


BASE = Model(6.8e-6, 0.1, 0.8)
SECOND = Model(3e-5, 0.5, 0.2)
REGIMES = {"baseline": BASE, "strong": SECOND}


class _Model(ct.Structure):
    _fields_ = [
        ("copies", ct.c_int),
        ("cap", ct.c_int),
        ("mutation", ct.c_double),
        ("advantage", ct.c_double),
        ("advantageous", ct.c_double),
        ("half_days", ct.c_double),
    ]


def _vec(dtype, dims):
    return np.ctypeslib.ndpointer(dtype=dtype, ndim=dims, flags="C_CONTIGUOUS")


@lru_cache(maxsize=1)
def library():
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    link = "-dynamiclib" if sys.platform == "darwin" else "-shared"
    cpu = ["-march=native"] if platform.machine() in ("x86_64", "AMD64") else []
    with TemporaryDirectory(prefix="kk26-") as temp:
        output = Path(temp) / "moran.so"
        try:
            subprocess.run(
                [
                    *compiler,
                    "-std=c++17",
                    "-O3",
                    *cpu,
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-fPIC",
                    link,
                    str(HERE / "moran.cpp"),
                    "-o",
                    str(output),
                ],
                check=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError("install a C++17 compiler or set CXX") from exc
        lib = ct.CDLL(str(output))
    common = [ct.POINTER(_Model), _vec(np.float64, 1), _vec(np.uint64, 1), ct.c_int]
    lib.kk26_forward.argtypes = [
        *common,
        ct.c_double,
        _vec(np.float64, 2),
        ct.POINTER(ct.c_double),
    ]
    lib.kk26_forward.restype = ct.c_int
    lib.kk26_sample.argtypes = [
        *common,
        ct.c_double,
        ct.c_int,
        ct.c_uint64,
        _vec(np.int32, 2),
        _vec(np.int32, 1),
    ]
    lib.kk26_sample.restype = ct.c_int
    lib.kk26_passage.argtypes = [
        ct.POINTER(_Model),
        ct.c_int,
        _vec(np.uint64, 1),
        ct.c_int,
        _vec(np.float64, 2),
    ]
    lib.kk26_passage.restype = None
    return lib


def _ages(model, ages):
    days = np.ascontiguousarray(ages, dtype=np.float64)
    if days.ndim != 1 or not days.size or not (0 <= days[0] <= days[-1] <= 20000):
        raise ValueError("ages must be a nonempty vector in [0, 20000] days")
    if not np.isfinite(days).all() or np.any(np.diff(days) <= 0):
        raise ValueError("ages must be finite and strictly increasing")
    return days, np.floor(days * model.rate).astype(np.uint64)


def tails(bins):
    """Four strict exceedance probabilities or counts from five exclusive bins."""
    return np.cumsum(bins[..., :0:-1], axis=-1)[..., ::-1]


def forward(model, ages=AGES, *, cap=None, step=None, states=False):
    """Return bin masses, leakage and negative mass at each report age."""
    config = model.encode(cap)
    days, ticks = _ages(model, ages)
    if step is not None and not 0 < step <= 30:
        raise ValueError("split step must be in (0, 30] days")
    if states and step is not None:
        raise ValueError("joint state output is available for the discrete backend")
    out = np.zeros((len(days), 7))
    vector = None
    if states:
        n = config.cap + 1
        vector = np.zeros(n * (config.copies + 1) - n * (n - 1) // 2)
    ptr = None if vector is None else vector.ctypes.data_as(ct.POINTER(ct.c_double))
    res = library().kk26_forward(
        ct.byref(config),
        days,
        ticks,
        len(days),
        0.0 if step is None else step,
        out,
        ptr,
    )
    if res or not np.isfinite(out).all():
        raise RuntimeError("forward solver failed")
    return (out, vector) if states else out


def cohort(model, cells, seed, ages=AGES, *, threads=1, atp_step=None):
    """Simulate the same cells at all ages. Return five bin counts and ATP stops."""
    config = model.encode()
    days, ticks = _ages(model, ages)
    if not 1 <= cells <= 10_000_000 or not 0 <= seed < 2**64 or threads < 1:
        raise ValueError("invalid cohort size, seed or worker count")
    if atp_step is not None and (
        atp_step not in (1.0, 0.5)
        or np.any(np.abs(days * 24 / atp_step - np.round(days * 24 / atp_step)) > 1e-8)
    ):
        raise ValueError("ATP ages must lie on the hourly or half-hour mesh")
    lib = library()

    def draw(bounds):
        start, end = bounds
        counts = np.zeros((len(days), 5), dtype=np.int32)
        stops = np.zeros(len(days), dtype=np.int32)
        status = lib.kk26_sample(
            ct.byref(config),
            days,
            ticks,
            len(days),
            atp_step or 0,
            end - start,
            (seed + int(start)) % 2**64,
            counts,
            stops,
        )
        if status:
            raise RuntimeError("synthetic cohort failed")
        return counts, stops

    if threads == 1:
        return draw((0, cells))
    edges = np.linspace(0, cells, min(cells, threads * 4) + 1, dtype=int)
    with ThreadPoolExecutor(max_workers=threads) as pool:
        results = list(pool.map(draw, zip(edges[:-1], edges[1:], strict=True)))
    return np.sum([c for c, _ in results], axis=0), np.sum(
        [s for _, s in results], axis=0
    )


def native_model(model, initial=None):
    """The Poisson-clock counterpart of Eq. 7 in CRN-craft."""
    import crn_craft as crn

    model.encode()
    doc = crn.Model("mtdna_moran_inference")
    n = model.copies
    w, u, a = (n, 0, 0) if initial is None else initial
    if min(w, u, a) < 0 or w + u + a != n:
        raise ValueError("the initial state must sum to the population")
    doc.add_species({"W": w, "U": u, "A": a})
    doc.add_parameter(
        {
            "rb": model.rate / n,
            "mut": model.mutation,
            "fitness": 1 + model.advantage,
            "frac": model.advantageous,
        }
    )
    weights = {"W": "W*(1-mut)", "U": "U+W*mut*(1-frac)", "A": "fitness*A+W*mut*frac"}
    for died in weights:
        for daughter, weight in weights.items():
            if died != daughter:
                doc.add_reaction(
                    f"rb*{died}*({weight})/(W+U+fitness*A)",
                    [died],
                    [daughter],
                    name=f"{died}_to_{daughter}",
                )
    return doc


def passage(model, max_days=6000, interval=30):
    """Conditional time from one advantageous copy to at least 90% abundance."""
    if not 0 < interval <= max_days <= 20000 or model.copies < 2:
        raise ValueError("invalid first-passage reporting grid")
    days = np.append(np.arange(0, max_days, interval, dtype=float), float(max_days))
    days, ticks = _ages(model, days)
    config = model.encode()
    config.mutation = 0
    target = (9 * model.copies + 9) // 10
    log_fitness = math.log1p(model.advantage)
    success = (
        math.expm1(-log_fitness) / math.expm1(-target * log_fitness)
        if log_fitness
        else 1 / target
    )
    out = np.zeros((len(days), 3))
    library().kk26_passage(ct.byref(config), target, ticks, len(days), out)
    cdf = out[:, 1] / success
    known = {
        f"p{int(p * 100):02d}"
        if p != 0.025 and p != 0.975
        else f"p{100 * p:.1f}": float(np.interp(p, cdf, days)) if p <= cdf[-1] else None
        for p in (0.025, 0.5, 0.975)
    }
    return {
        "days": days,
        "steps": ticks,
        "absorption": out,
        "cdf": cdf,
        "target": target,
        "success": success,
        "right_tail": max(0.0, 1 - cdf[-1]),
        "quantiles_days": known,
    }
