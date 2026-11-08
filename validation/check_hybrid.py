"""Exact-law gates of the hybrid method and its ATP equivalence with thinning.

The protocol is declared in validation/hybrid/RESULTS.md: five DKW gates
against exact laws (20,000 replicates each, family alpha 0.001),
deterministic checks of bounds and of the integrator, and four predeclared
two-sample tests of hybrid on the continuous ATP model against exact thinning
on the legacy ATP model (400 replicates per arm, or 4,000 with --large;
family alpha 0.01, Bonferroni). For rules-model/0.4 independent variables:
four DKW gates (C, 20,000 replicates each, family alpha 0.001) and four
two-sample tests of independent against dependent integration (D, family
alpha 0.01, Bonferroni). Uses only the installed wheel, NumPy and the
standard library. Every row is reported, including failures; the report holds
fixed-work evidence only and is identical under python -O. Run times are
printed to standard error.
"""

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

from crn_craft import rules

ROOT = Path(__file__).resolve().parents[1]
GATE_N = 20_000
GATE_ALPHA = 0.001
GATE_K = 5
EPSILON = math.sqrt(math.log(2 * GATE_K / GATE_ALPHA) / (2 * GATE_N))
RTOL, ATOL = 1e-6, 1e-9
HYBRID = {"kind": "hybrid", "rtol": RTOL, "atol": ATOL}
ATP_ALPHA = 0.01
ATP_TESTS = 4
ATP_LEVEL = ATP_ALPHA / ATP_TESTS
SIZES = {"default": (400, 7101, 7102), "large": (4000, 7201, 7202)}
TAKEOVER = "N > 0 && mutants >= 0.9 * N"
ATP_END = 25_920.0
ATP_TIMES = [0.0, 8640.0, 17280.0, ATP_END]
INDEPENDENT_K = 4
INDEPENDENT_EPSILON = math.sqrt(math.log(2 * INDEPENDENT_K / GATE_ALPHA) / (2 * GATE_N))
EQUIVALENCE_ALPHA = 0.01
EQUIVALENCE_TESTS = 4
EQUIVALENCE_LEVEL = EQUIVALENCE_ALPHA / EQUIVALENCE_TESTS


def model(variables, rules_, *, initial=None, parameters=None, schema="0.3"):
    return {
        "schema": f"crn-craft.rules-model/{schema}",
        "name": "hybrid_gate",
        "parameters": parameters or {},
        "types": [{"name": "A"}],
        "initial": [{"type": "A", "count": n} for n in [initial] if n],
        "observables": {"A": {"type": "A"}},
        "rules": rules_,
        "continuous": variables,
    }


def death(rate):
    return {
        "name": "death",
        "reactants": [{"label": "a", "type": "A"}],
        "rate": rate,
        "products": [],
    }


def arrival(propensity):
    return {
        "name": "arrival",
        "reactants": [],
        "propensity": propensity,
        "products": [{"type": "A"}],
    }


def simulate(document, *, end, times, replicates, seed, method, trackers=()):
    """The result and its run time in seconds."""
    prepared = rules.ensemble(
        document,
        end=end,
        report_times=times,
        replicates=replicates,
        seed=seed,
        threads=2,
        keep_paths=replicates,
        trackers=list(trackers),
        method=method,
        max_events=100_000_000,
        max_sampling_steps=1_000_000_000,
        max_output_bytes=1 << 30,
    )
    started = time.perf_counter()
    result = prepared.simulate()
    return result, time.perf_counter() - started


def column(meta, name):
    """Replicate × report time values of one observable."""
    j = meta["observables"].index(name)
    return np.array([[row[j] for row in p["values"]] for p in meta["paths"]])


def ecdf_distance_continuous(samples, cdf, horizon):
    """sup over [0, horizon] of |F_n - F| for times censored at the horizon."""
    n = len(samples)
    seen = np.sort(np.asarray([s for s in samples if s <= horizon], dtype=float))
    distance = 0.0
    for i, x in enumerate(seen, start=1):
        f = cdf(float(x))
        distance = max(distance, i / n - f, f - (i - 1) / n)
    return max(distance, cdf(horizon) - len(seen) / n)


def ecdf_distance_discrete(samples, pmf):
    """max over integers k of |F_n(k) - F(k)|, with pmf(k) for k >= 0."""
    values = np.asarray(samples, dtype=float)
    if not (np.all(values == np.round(values)) and np.all(values >= 0)):
        raise RuntimeError("counts must be nonnegative integers")
    values = values.astype(np.int64)
    counts = np.bincount(values)
    empirical, exact, distance, k = 0, 0.0, 0.0, 0
    while k < len(counts) or exact < 1.0 - 1e-12:
        empirical += int(counts[k]) if k < len(counts) else 0
        exact += pmf(k)
        distance = max(distance, abs(empirical / len(values) - exact))
        k += 1
    return distance


def poisson(mean):
    return lambda k: math.exp(k * math.log(mean) - mean - math.lgamma(k + 1))


def binomial(n, p):
    return lambda k: math.comb(n, k) * p**k * (1 - p) ** (n - k) if k <= n else 0.0


def ks_two_sample(first, second):
    """Two-sample Kolmogorov-Smirnov distance and asymptotic p-value.

    The p-value is the Kolmogorov law at (sqrt(m) + 0.12 + 0.11 / sqrt(m)) D
    with m = n1 n2 / (n1 + n2) (Stephens 1970); ties, including +inf for
    values that never occur, make it conservative.
    """
    a = np.sort(np.asarray(first, dtype=float))
    b = np.sort(np.asarray(second, dtype=float))
    pooled = np.unique(np.concatenate([a, b]))
    fa = np.searchsorted(a, pooled, side="right") / len(a)
    fb = np.searchsorted(b, pooled, side="right") / len(b)
    distance = float(np.max(np.abs(fa - fb)))
    m = len(a) * len(b) / (len(a) + len(b))
    lam = (math.sqrt(m) + 0.12 + 0.11 / math.sqrt(m)) * distance
    if lam < 0.3:
        return distance, 1.0
    total, sign = 0.0, 1.0
    for j in range(1, 101):
        term = sign * 2.0 * math.exp(-2.0 * j * j * lam * lam)
        total += term
        if abs(term) <= 1e-12 * abs(total):
            break
        sign = -sign
    return distance, min(1.0, max(0.0, total))


def fisher_exact(a, n1, c, n2):
    """Two-sided Fisher exact test of a events in n1 against c in n2."""
    k, n = a + c, n1 + n2

    def probability(x):
        return math.comb(k, x) * math.comb(n - k, n1 - x) / math.comb(n, n1)

    observed = probability(a)
    support = range(max(0, k - n2), min(k, n1) + 1)
    return min(
        1.0, sum(p for x in support if (p := probability(x)) <= observed * (1 + 1e-7))
    )


def exact_gates(timings):
    gates, checks, power = [], [], []

    def gate(name, distance):
        gates.append(
            {
                "gate": name,
                "replicates": GATE_N,
                "distance": distance,
                "epsilon": EPSILON,
                "passed": distance <= EPSILON,
            }
        )

    def check(name, condition, detail=None):
        row = {"check": name, "passed": bool(condition)}
        if detail is not None:
            row["detail"] = detail
        checks.append(row)

    def complete(result):
        meta = result.metadata
        return meta["completed_replicates"] == meta["requested_replicates"]

    # A1: x' = -k x drives the death of one molecule at rate c x.
    k, x0, c, horizon = 0.5, 2.0, 1.0, 10.0
    doc = model(
        [{"name": "x", "initial": x0, "derivative": "-k * x"}],
        [death("c * x")],
        initial=1,
        parameters={"k": k, "c": c},
    )
    result, timings["A1"] = simulate(
        doc,
        end=horizon,
        times=[horizon],
        replicates=GATE_N,
        seed=6101,
        method=HYBRID,
        trackers=[{"observable": "A"}],
    )
    check("A1_complete", complete(result))
    times = [
        math.inf if t is None else t
        for t in (o["trackers"][0]["last_extinction"] for o in result.observations)
    ]
    check("A1_observed", len(times) == GATE_N)

    def death_cdf(rate):
        return lambda t: 1.0 - math.exp(-rate * x0 * (1.0 - math.exp(-k * t)) / k)

    gate(
        "A1_decay_driven_death", ecdf_distance_continuous(times, death_cdf(c), horizon)
    )
    power.append(
        {
            "gate": "A1",
            "wrong_law": "c = 1.1",
            "distance": ecdf_distance_continuous(times, death_cdf(1.1), horizon),
        }
    )

    # A2: x' = b and immigration with propensity x: Poisson(b T^2 / 2).
    b, horizon = 2.0, 3.0
    doc = model(
        [{"name": "x", "initial": 0, "derivative": "b"}],
        [arrival("x")],
        parameters={"b": b},
    )
    result, timings["A2"] = simulate(
        doc, end=horizon, times=[horizon], replicates=GATE_N, seed=6102, method=HYBRID
    )
    check("A2_complete", complete(result))
    count = column(result.metadata, "A")[:, 0]
    gate("A2_linear_rate_poisson", ecdf_distance_discrete(count, poisson(b * 4.5)))
    power.append(
        {
            "gate": "A2",
            "wrong_law": "b = 2.2",
            "distance": ecdf_distance_discrete(count, poisson(2.2 * 4.5)),
        }
    )

    # A3: x' = -1 from 3 with a stop at 0; deaths at rate 0.2 meanwhile.
    doc = model(
        [
            {
                "name": "x",
                "initial": 3,
                "derivative": "-1",
                "lower_bound": 0,
                "on_bound": "stop",
            }
        ],
        [death("0.2")],
        initial=20,
    )
    result, timings["A3"] = simulate(
        doc, end=5.0, times=[1.5, 4.0, 5.0], replicates=GATE_N, seed=6103, method=HYBRID
    )
    meta = result.metadata
    stops = np.asarray(result.records["stop_time"], dtype=float)
    check("A3_all_stopped", meta["stopped_replicates"] == GATE_N)
    check(
        "A3_stop_time",
        bool(np.all(np.abs(stops - 3.0) <= 1e-9)),
        float(np.max(np.abs(stops - 3.0))),
    )
    x = column(meta, "x")
    check("A3_x_before_stop", bool(np.all(np.abs(x[:, 0] - 1.5) <= 1e-9)))
    check("A3_x_after_stop", bool(np.all(x[:, 1:] == 0.0)))
    gate(
        "A3_count_at_stop_binomial",
        ecdf_distance_discrete(column(meta, "A")[:, 1], binomial(20, math.exp(-0.6))),
    )

    # A4: x' = -1 from 2, clamped at 0; immigration with propensity x.
    doc = model(
        [
            {
                "name": "x",
                "initial": 2,
                "derivative": "-1",
                "lower_bound": 0,
                "on_bound": "clamp",
            }
        ],
        [arrival("x")],
    )
    result, timings["A4"] = simulate(
        doc, end=5.0, times=[3.0, 5.0], replicates=GATE_N, seed=6104, method=HYBRID
    )
    meta = result.metadata
    check(
        "A4_all_completed",
        complete(result) and meta["stopped_replicates"] == 0,
    )
    check("A4_x_clamped", bool(np.all(column(meta, "x") == 0.0)))
    gate(
        "A4_clamp_poisson",
        ecdf_distance_discrete(column(meta, "A")[:, 1], poisson(2.0)),
    )

    # A5: immigration-death without continuous variables.
    doc = model(
        [],
        [arrival("10"), death("1")],
    )
    del doc["continuous"]
    result, timings["A5"] = simulate(
        doc, end=3.0, times=[3.0], replicates=GATE_N, seed=6105, method=HYBRID
    )
    check("A5_complete", complete(result))
    count = column(result.metadata, "A")[:, 0]
    scale = 1.0 - math.exp(-3.0)
    gate(
        "A5_jump_process_poisson", ecdf_distance_discrete(count, poisson(10.0 * scale))
    )
    power.append(
        {
            "gate": "A5",
            "wrong_law": "immigration 10.5",
            "distance": ecdf_distance_discrete(count, poisson(10.5 * scale)),
        }
    )

    # A6: no event; the reports follow the ODE solutions.
    doc = model(
        [
            {"name": "x", "initial": 1, "derivative": "-x"},
            {"name": "z", "initial": 0.1, "derivative": "z * (1 - z)"},
        ],
        [death("1")],
    )
    grid = [float(t) for t in range(11)]
    result, _ = simulate(
        doc, end=10.0, times=grid, replicates=2, seed=6106, method=HYBRID
    )
    meta = result.metadata
    worst = 0.0
    for name, exact in (
        ("x", lambda t: math.exp(-t)),
        ("z", lambda t: 1.0 / (1.0 + 9.0 * math.exp(-t))),
    ):
        values = column(meta, name)
        for j, t in enumerate(grid):
            allowed = 100.0 * (ATOL + RTOL * abs(exact(t)))
            worst = max(worst, float(np.max(np.abs(values[:, j] - exact(t)))) / allowed)
    check("A6_ode_accuracy", worst <= 1.0, worst)
    for row in power:
        row["rejected"] = row["distance"] > EPSILON
    return gates, checks, power


def gompertz(independent, integrated=False):
    """C1 (m drives the deaths) or C2 (x' = m drives them)."""
    variables = [
        {
            "name": "m",
            "initial": "m0",
            "derivative": "g * m",
            "independent": independent,
        }
    ]
    if integrated:
        variables.append(
            {"name": "x", "initial": 0, "derivative": "m", "independent": False}
        )
    return model(
        variables,
        [death("x" if integrated else "m")],
        initial=20,
        parameters={"g": 0.5, "m0": 0.1},
        schema="0.4",
    )


def stop_bound(independent):
    """C3: x' = -1 from 3 with a stop at 0; deaths at rate 0.2."""
    return model(
        [
            {
                "name": "x",
                "initial": 3,
                "derivative": "-1",
                "lower_bound": 0,
                "on_bound": "stop",
                "independent": independent,
            }
        ],
        [death("0.2")],
        initial=20,
        schema="0.4",
    )


def pulse(independent):
    """C4: an eight-compartment linear chain drives arrivals."""
    variables = [
        {
            "name": "u1",
            "initial": 1,
            "derivative": "-k * u1",
            "independent": independent,
        }
    ] + [
        {
            "name": f"u{i}",
            "initial": 0,
            "derivative": f"k * (u{i - 1} - u{i})",
            "independent": independent,
        }
        for i in range(2, 9)
    ]
    return model(
        variables,
        # The integrated u8 can come out just below 0 near t = 0.
        [arrival("c * max(0, u8)"), death("d")],
        parameters={"k": 2.0, "c": 50.0, "d": 0.1},
        schema="0.4",
    )


def pulse_mean():
    """Poisson mean of C4: c e^{-dT} (k/a)^7 P(8, aT) / a with a = k - d."""
    k, c, d, horizon = 2.0, 50.0, 0.1, 12.0
    a = k - d
    x = a * horizon
    lower = 1.0 - math.exp(-x) * sum(x**i / math.factorial(i) for i in range(8))
    return c * math.exp(-d * horizon) * (k / a) ** 7 * lower / a


# label, model builder, horizon, report times, seed of the true arm and of
# the false arm
INDEPENDENT_CASES = (
    ("C1", lambda flag: gompertz(flag), 4.0, [0.0, 1.0, 2.0, 3.0, 4.0], 6201, 6211),
    ("C2", lambda flag: gompertz(flag, True), 4.0, [4.0], 6202, 6212),
    ("C3", stop_bound, 5.0, [1.5, 4.0, 5.0], 6203, 6213),
    ("C4", pulse, 12.0, [12.0], 6204, 6214),
)


def independent_gates(timings):
    """C: exact laws with independent variables; returns the true-arm counts."""
    gates, checks, counts = [], [], {}

    def gate(name, distance):
        gates.append(
            {
                "gate": name,
                "replicates": GATE_N,
                "distance": distance,
                "epsilon": INDEPENDENT_EPSILON,
                "passed": distance <= INDEPENDENT_EPSILON,
            }
        )

    def check(name, condition, detail=None):
        row = {"check": name, "passed": bool(condition)}
        if detail is not None:
            row["detail"] = detail
        checks.append(row)

    metas = {}
    for label, build, end, times, seed, _ in INDEPENDENT_CASES:
        result, timings[label] = simulate(
            build(True),
            end=end,
            times=times,
            replicates=GATE_N,
            seed=seed,
            method=HYBRID,
        )
        meta = result.metadata
        metas[label] = (result, meta)
        counts[label] = column(meta, "A")[:, -1]
        if label != "C3":
            check(
                f"{label}_complete",
                meta["completed_replicates"] == GATE_N
                and meta["stopped_replicates"] == 0,
            )
    survival = math.exp(-0.1 * (math.exp(2.0) - 1.0) / 0.5)
    gate(
        "C1_gompertz_binomial",
        ecdf_distance_discrete(counts["C1"], binomial(20, survival)),
    )
    m = column(metas["C1"][1], "m")
    check("C1_m_shared", bool(np.all(m == m[0])))
    worst = max(
        abs(float(m[0, j]) - 0.1 * math.exp(0.5 * t))
        / (100.0 * (ATOL + RTOL * 0.1 * math.exp(0.5 * t)))
        for j, t in enumerate([0.0, 1.0, 2.0, 3.0, 4.0])
    )
    check("C1_m_accuracy", worst <= 1.0, worst)
    integrated = (0.1 / 0.5) * ((math.exp(2.0) - 1.0) / 0.5 - 4.0)
    gate(
        "C2_dependent_reads_independent_binomial",
        ecdf_distance_discrete(counts["C2"], binomial(20, math.exp(-integrated))),
    )
    result, meta = metas["C3"]
    stops = np.asarray(result.records["stop_time"], dtype=float)
    check("C3_all_stopped", meta["stopped_replicates"] == GATE_N)
    check(
        "C3_stop_time",
        bool(np.all(np.abs(stops - 3.0) <= 1e-9)),
        float(np.max(np.abs(stops - 3.0))),
    )
    x = column(meta, "x")
    check("C3_x_before_stop", bool(np.all(np.abs(x[:, 0] - 1.5) <= 1e-9)))
    check("C3_x_after_stop", bool(np.all(x[:, 1:] == 0.0)))
    gate(
        "C3_independent_stop_binomial",
        ecdf_distance_discrete(counts["C3"], binomial(20, math.exp(-0.6))),
    )
    gate(
        "C4_pulse_poisson", ecdf_distance_discrete(counts["C4"], poisson(pulse_mean()))
    )
    by_threads = []
    for threads in (1, 2):
        run = rules.ensemble(
            pulse(True),
            end=12.0,
            report_times=[6.0, 12.0],
            replicates=200,
            seed=6205,
            threads=threads,
            keep_paths=200,
            method=HYBRID,
        )
        meta = dict(run.simulate().metadata)
        meta.pop("run_sha256")
        by_threads.append(meta)
    check("C4_threads", by_threads[0] == by_threads[1])
    return gates, checks, counts


def independent_equivalence(counts, timings):
    """D: the C models with every variable dependent, against the C samples."""
    rows = []
    for (label, build, end, times, _, seed), test in zip(
        INDEPENDENT_CASES, ("D1", "D2", "D3", "D4"), strict=True
    ):
        result, timings[test] = simulate(
            build(False),
            end=end,
            times=times,
            replicates=GATE_N,
            seed=seed,
            method=HYBRID,
        )
        meta = result.metadata
        finished = (
            meta["stopped_replicates"] == GATE_N
            if label == "C3"
            else meta["completed_replicates"] == GATE_N
            and meta["stopped_replicates"] == 0
        )
        other = column(meta, "A")[:, -1]
        distance, p = ks_two_sample(counts[label], other)
        rows.append(
            {
                "test": f"{test}_{label}_independent_vs_dependent",
                "kind": "two-sample Kolmogorov-Smirnov",
                "replicates_per_arm": GATE_N,
                "seed_dependent": seed,
                "mean_count": [float(np.mean(counts[label])), float(np.mean(other))],
                "statistic": distance,
                "p_value": p,
                "passed": p >= EQUIVALENCE_LEVEL and finished,
            }
        )
    return rows


def atp_equivalence(size, timings):
    replicates, seed_thinning, seed_hybrid = SIZES[size]
    tracker = [{"observable": "mutants", "first_passage": TAKEOVER}]
    arms = {}
    for name, path, method, seed in (
        ("thinning", "mtdna-three-class-atp.rules.json", "thinning", seed_thinning),
        ("hybrid", "mtdna-three-class-atp-continuous.rules.json", HYBRID, seed_hybrid),
    ):
        document = rules.load(ROOT / "models/rules" / path)
        result, timings[f"B_{name}"] = simulate(
            document,
            end=ATP_END,
            times=ATP_TIMES,
            replicates=replicates,
            seed=seed,
            method=method,
            trackers=tracker,
        )
        meta = result.metadata
        passage = [
            math.inf if t is None else t
            for t in (o["trackers"][0]["first_passage"] for o in result.observations)
        ]
        arms[name] = {
            "sampler": meta["sampler"],
            "complete": meta["completed_replicates"] == replicates,
            "stopped": meta["stopped_replicates"],
            "events": meta["events"],
            "takeover": passage,
            "atp": column(meta, "ATP")[:, -1],
            "n": column(meta, "N")[:, -1],
        }
    one, two = arms["thinning"], arms["hybrid"]
    rows = []
    for test, key in (
        ("B1_takeover_time", "takeover"),
        ("B2_final_atp", "atp"),
        ("B3_final_population", "n"),
    ):
        distance, p = ks_two_sample(one[key], two[key])
        rows.append(
            {
                "test": test,
                "kind": "two-sample Kolmogorov-Smirnov",
                "statistic": distance,
                "p_value": p,
                "passed": p >= ATP_LEVEL and one["complete"] and two["complete"],
            }
        )
    p = fisher_exact(one["stopped"], replicates, two["stopped"], replicates)
    rows.append(
        {
            "test": "B4_atp_exhaustion",
            "kind": "Fisher exact, two-sided",
            "statistic": [one["stopped"], two["stopped"]],
            "p_value": p,
            "passed": p >= ATP_LEVEL and one["complete"] and two["complete"],
        }
    )

    def summary(arm):
        takeover = np.asarray(arm["takeover"])
        return {
            "sampler": arm["sampler"],
            "complete": arm["complete"],
            "events": arm["events"],
            "stopped": arm["stopped"],
            "taken_over": int(np.sum(np.isfinite(takeover))),
            "median_takeover_hours": (
                float(np.median(takeover[np.isfinite(takeover)]))
                if np.any(np.isfinite(takeover))
                else None
            ),
            "mean_final_atp": float(np.mean(arm["atp"])),
            "mean_final_population": float(np.mean(arm["n"])),
        }

    return {
        "size": size,
        "replicates_per_arm": replicates,
        "seeds": {"thinning": seed_thinning, "hybrid": seed_hybrid},
        "family_alpha": ATP_ALPHA,
        "per_test_level": ATP_LEVEL,
        "arms": {name: summary(arm) for name, arm in arms.items()},
        "tests": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--large", action="store_true", help="4,000 replicates per ATP arm"
    )
    args = parser.parse_args()
    timings = {}
    gates, checks, power = exact_gates(timings)
    equivalence = atp_equivalence("large" if args.large else "default", timings)
    independent, independent_checks, counts = independent_gates(timings)
    against_dependent = independent_equivalence(counts, timings)
    failed = (
        [g["gate"] for g in gates if not g["passed"]]
        + [c["check"] for c in checks if not c["passed"]]
        + [t["test"] for t in equivalence["tests"] if not t["passed"]]
        + [g["gate"] for g in independent if not g["passed"]]
        + [c["check"] for c in independent_checks if not c["passed"]]
        + [t["test"] for t in against_dependent if not t["passed"]]
    )
    report = {
        "schema": "crn-craft.hybrid-validation/0.1",
        "protocol": "validation/hybrid/RESULTS.md",
        "tolerances": {"rtol": RTOL, "atol": ATOL},
        "gates": {
            "family_alpha": GATE_ALPHA,
            "comparisons": GATE_K,
            "epsilon": EPSILON,
            "rows": gates,
        },
        "checks": checks,
        "power": power,
        "atp_equivalence": equivalence,
        "independent_gates": {
            "family_alpha": GATE_ALPHA,
            "comparisons": INDEPENDENT_K,
            "epsilon": INDEPENDENT_EPSILON,
            "rows": independent,
        },
        "independent_checks": independent_checks,
        "independent_vs_dependent": {
            "family_alpha": EQUIVALENCE_ALPHA,
            "per_test_level": EQUIVALENCE_LEVEL,
            "tests": against_dependent,
        },
        "failed": failed,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(text)
    print(text, end="")
    print(
        "run times (s): "
        + ", ".join(f"{key} {value:.1f}" for key, value in timings.items()),
        file=sys.stderr,
    )
    if failed:
        raise SystemExit("hybrid validation failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()
