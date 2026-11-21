"""Four promoter states, osmotic forcing and the mature RNA master equation."""

import ctypes
import functools
import hashlib
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import crn_craft as crn
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
# Fixed flags: no -march=native and no fused multiply-add, which compilers
# use by default on targets with FMA and which change the last bits.
BUILD_FLAGS = ("-std=c++17", "-O3", "-DNDEBUG", "-ffp-contract=off")
KEYS = (
    "k12",
    "k21",
    "k21_signal",
    "k23",
    "k32",
    "k34",
    "k43",
    "i1",
    "i2",
    "i3",
    "i4",
    "decay",
    "delay",
    "export",
    "nuclear_decay",
)


def parameters(record):
    return np.array([record[k] for k in KEYS])


def hog1(time, dose, delay):
    t = np.maximum(0.0, np.asarray(time, dtype=float) - delay)
    rise = -np.expm1(-0.006094509227679889 * t)
    response = rise * np.exp(-0.001679176848801883 / (dose + 0.04096887973322170) * t)
    return (
        response / (1 + 46.51137270724188 * response)
    ) ** 5.890086837674720 / 1.0727183838607e-10


def switching(x, k21):
    k12, _, _, k23, k32, k34, k43 = x[:7]
    return (
        (0, 1, k12),
        (1, 0, k21),
        (1, 2, k23),
        (2, 1, k32),
        (2, 3, k34),
        (3, 2, k43),
    )


def native(x, signal=0.0, spatial=False):
    model = crn.Model("yeast_mature_rna")
    s = [f"S{k + 1}" for k in range(4)]
    model.add_species({p: int(k == 0) for k, p in enumerate(s)})
    model.add_species("N", 0)
    for src, dst, rate in switching(x, max(0.0, x[1] - x[2] * signal)):
        model.add_reaction(float(rate), [s[src]], [s[dst]], name=f"switch{src}{dst}")
    for k, p in enumerate(s):
        model.add_reaction(float(x[7 + k]), [p], [p, "N"], name=f"initiate{k}")
    model.add_reaction(
        float(x[14] if spatial else x[11]), ["N"], [], name="nuclear_decay"
    )
    if spatial:
        model.add_species("C", 0)
        model.add_reaction(float(x[13]), ["N"], ["C"], name="export")
        model.add_reaction(float(x[11]), ["C"], [], name="cytoplasmic_decay")
    return model


@functools.cache
def compiled():
    """Build one shared library per source, compiler and flags outside the repository."""
    source = ROOT / "cme.cpp"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    link = "-dynamiclib" if sys.platform == "darwin" else "-shared"
    command = [*compiler, *BUILD_FLAGS, "-fPIC", link]
    digest = hashlib.sha256(
        source.read_bytes() + "\0".join(command).encode()
    ).hexdigest()[:20]
    user = getattr(os, "getuid", lambda: "user")()
    directory = Path(tempfile.gettempdir()) / f"crn-craft-yeast-{user}"
    directory.mkdir(mode=0o700, exist_ok=True)
    lib = directory / f"{digest}.so"
    if not lib.is_file():
        pending = directory / f"{digest}-{os.getpid()}.tmp.so"
        try:
            subprocess.run([*command, str(source), "-o", str(pending)], check=True)
        except FileNotFoundError as exc:
            raise RuntimeError("install a C++17 compiler or set CXX") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"{' '.join(compiler)} could not build cme.cpp as C++17; "
                "set CXX to another compiler"
            ) from exc
        os.replace(pending, lib)
    dll = ctypes.CDLL(str(lib))
    method = dll.advance
    real = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    method.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        real,
        real,
        real,
    ]
    method.restype = None
    mean = dll.advance_mean
    mean.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        real,
        real,
        real,
    ]
    mean.restype = None
    return method, mean


def nonspatial(x, dose, cap, step, minutes=55):
    stride = round(60 / step)
    if not np.isclose(stride * step, 60):
        raise ValueError("Step must divide one minute")
    steps = stride * minutes
    forcing = hog1(step * np.arange(steps + 1), dose, x[12])
    out = np.empty((minutes + 1, cap + 1, 4))
    compiled()[0](cap, steps, stride, step, np.ascontiguousarray(x), forcing, out)
    return out


def quick_means(x, dose, step, spatial_model=False, minutes=55):
    stride = round(60 / step)
    steps = stride * minutes
    forcing = hog1(step * np.arange(steps + 1), dose, x[12])
    out = np.empty((minutes + 1, 2))
    compiled()[1](
        steps, stride, step, int(spatial_model), np.ascontiguousarray(x), forcing, out
    )
    return out if spatial_model else out[:, :1]


def matrix(x, cap, cyt_cap=None, with_sink=False):
    """Return the constant part and the coefficient of the S2→S1 switch."""
    c_count = 1 if cyt_cap is None else cyt_cap + 1
    n_count = cap + 1
    n = 4 * c_count * n_count
    shape = (n + int(with_sink),) * 2
    states = np.arange(n).reshape(c_count, n_count, 4)
    rows, cols, rates = [], [], []

    def events(source, dest, rate):
        source = np.ravel(source)
        dest = np.broadcast_to(dest, source.shape)
        rate = np.broadcast_to(rate, source.shape)
        rows.append(source)
        cols.append(source)
        rates.append(-rate)
        if with_sink:
            rows.append(dest)
            cols.append(source)
            rates.append(rate)
        else:
            present = dest < n
            rows.append(dest[present])
            cols.append(source[present])
            rates.append(rate[present])

    for src, dst, rate in switching(x, 0.0):
        if rate:
            events(states[..., src], states[..., dst].ravel(), rate)
    for k in range(4):
        events(states[:, :-1, k], states[:, 1:, k].ravel(), x[7 + k])
        events(states[:, -1, k], np.int64(n), x[7 + k])
    n_decay = x[11] if cyt_cap is None else x[14]
    counts_n = np.broadcast_to(np.arange(1, n_count)[None, :, None], (c_count, cap, 4))
    events(states[:, 1:, :], states[:, :-1, :].ravel(), n_decay * counts_n.ravel())
    if cyt_cap is not None:
        counts_c = np.broadcast_to(
            np.arange(1, c_count)[:, None, None], (cyt_cap, n_count, 4)
        )
        events(states[1:], states[:-1].ravel(), x[11] * counts_c.ravel())
        events(states[:-1, 1:], states[1:, :-1].ravel(), x[13] * counts_n[:-1].ravel())
        events(states[-1, 1:], np.int64(n), x[13] * counts_n[-1].ravel())
    a = coo_matrix(
        (np.concatenate(rates), (np.concatenate(rows), np.concatenate(cols))),
        shape=shape,
    ).tocsc()
    originate = states[..., 1].ravel()
    arrive = states[..., 0].ravel()
    b = coo_matrix(
        (
            np.r_[np.ones(len(arrive)), -np.ones(len(arrive))],
            (np.r_[arrive, originate], np.tile(originate, 2)),
        ),
        shape=shape,
    ).tocsc()
    return a, b


def spatial(x, dose, minutes, ncap, ccap, tolerance=2e-6):
    """Sparse BDF on a joint projection, with the sink from missing mass."""
    a, b = matrix(x, ncap, ccap)
    initial = np.zeros(a.shape[0])
    initial[0] = 1.0

    def k(t):
        return max(0, x[1] - x[2] * float(hog1(t, dose, x[12])))

    def rhs(t, state):
        return a @ state + k(t) * (b @ state)

    analysis = solve_ivp(
        rhs,
        (0, 60 * max(minutes)),
        initial,
        method="BDF",
        t_eval=np.asarray(minutes) * 60,
        jac=lambda t, p: a + k(t) * b,
        rtol=tolerance,
        atol=tolerance * 1e-4,
        max_step=20.0,
    )
    if not analysis.success:
        raise RuntimeError(analysis.message)
    return analysis.y.T.reshape(len(minutes), ccap + 1, ncap + 1, 4)


def exact_means(x, dose, minutes, spatial_model=False):
    """Linear closed equations for promoter occupancy and each first moment."""
    minutes = np.asarray(minutes)
    initial = np.zeros(6 if spatial_model else 5)
    initial[0] = 1
    if max(minutes) == 0:
        return initial[None, 4:]

    def rhs(t, y):
        out = np.zeros_like(y)
        back = max(0.0, x[1] - x[2] * float(hog1(t, dose, x[12])))
        for src, dst, rate in switching(x, back):
            flow = y[src] * rate
            out[src] -= flow
            out[dst] += flow
        init = np.dot(x[7:11], y[:4])
        if spatial_model:
            out[4] = init - (x[13] + x[14]) * y[4]
            out[5] = x[13] * y[4] - x[11] * y[5]
        else:
            out[4] = init - x[11] * y[4]
        return out

    solution = solve_ivp(
        rhs,
        (0.0, 60.0 * max(minutes)),
        initial,
        method="BDF",
        rtol=1e-7,
        atol=1e-9,
        max_step=20,
        t_eval=60 * minutes,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[4:].T
