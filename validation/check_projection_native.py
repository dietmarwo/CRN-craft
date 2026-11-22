"""Native prepared projections against the pure-Python/SciPy reference.

Compares ``crn_craft.prepare_fsp`` (native enumeration, generator assembly,
interventions, uniformization, exponential midpoint rule and the
variable-order NDF integrator of ``solver="bdf"``) and
``crn_craft.stationary_fsp`` with ``reference_projection.py`` (a pure
Python and SciPy implementation with SciPy's BDF and sparse LU) on a
declared case set, and the native binomial and Poisson-tail kernels with
SciPy. Gates, declared before the first run:

- identical ``states`` arrays and prepared ids;
- per snapshot, the L1 distance of (probabilities, leakage) at most
  1e-12 + B * tolerance, B = 1 for uniformization and 2 for the midpoint
  rule and BDF; identical times, actions applied and solver;
- series-tail bounds at most the tolerance in both;
- work within 0.1 % + 2 units for uniformization;
- identical exception types and messages on input and budget failures;
- closed-form laws (``analytic``), independent of the reference: the
  largest probability gap at most 1e-11 for uniformization (binomial
  survival) or 1e-9 (Poisson thinning by a partition) and at most 2e-6 for the
  midpoint rule at tolerance 1e-7 and BDF at 1e-9 (Poisson laws under a
  sinusoidal callable and under a step table read at its jump), with a
  ``max_step`` whose multiples do not sum to the segment ends in floating
  point, and at most 1e-8 for the native BDF there; for BDF under
  ``2 + sin(3 t)`` (case C7) an L1 gap to the exact Poisson law of at most
  10 * tolerance at tolerances 1e-6 and 1e-8; stationary laws within 1e-12
  of a binomial and of a Poisson law;
- stationary laws: identical success or exception type and message; on
  success the largest probability difference at most 1e-10, both
  residuals within the solver's own acceptance, identical histories
  (bounds, convergence); cases with two closed classes and transient
  states, an absorbing state, a nearly decomposable stiff chain and a
  redirected birth-death chain with expansions; a BDF solve and a
  stationary law above ``max_factor_bytes`` are refused with a message
  that names the budget;
- binomial pmf within a relative 1e-12 of SciPy for n <= 1e5 and, where
  the two differ by more than that (SciPy loses accuracy for large n),
  within 1e-13 of a 60-digit Decimal evaluation; Poisson tail within a
  relative 1e-12 of SciPy; both where SciPy's value is at least 1e-280;
- forward sensitivities (``crn_craft.sensitivities``): for an
  immigration-death law at tolerance 1e-10 the derivative of the Poisson
  law in the birth rate within 1e-7; for a Hill-repressed birth through a
  parameter switch at tolerance 1e-10 the derivatives in the Hill
  coefficient and in the switched rate within 1e-6 of centered differences
  of uniformization solves at relative steps 1e-3 and 5e-4, with the error
  ratio of the two steps between 3 and 5 (second order) unless both gaps
  are below 5e-8; total probability within 1e-12 of one and each
  derivative sum, sink included, within 1e-9 of zero;
- retained state sets: retaining every state of the box gives the
  ordinary projection's states, prepared id and law (identical arrays);
  a cut set keeps total probability within 2e-12 of one (solve tolerance
  1e-12) with the cut mass in the sink.

``--quick`` skips the total-RNA yeast cases (midpoint and BDF, minutes 0
to 10, about a minute each in the reference); ``--full-yeast`` solves its 55 minutes (several
minutes in the reference). Timings go to standard error; the JSON on standard output is
the same under ``python`` and ``python -O``.
"""

import argparse
import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
from time import perf_counter

import crn_craft as crn
import numpy as np
from crn_craft import _native
from scipy.special import gammainc
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_projection as reference  # noqa: E402

FAILURES = []


def require(condition, message):
    if not condition:
        FAILURES.append(message)
    return bool(condition)


getcontext().prec = 70
PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640"
)
BERNOULLI = [Fraction(1, 6), Fraction(-1, 30), Fraction(1, 42), Fraction(-1, 30)]
BERNOULLI += [Fraction(5, 66), Fraction(-691, 2730), Fraction(7, 6)]


def log_factorial(n):
    """ln(n!) to about 60 digits: exact below 2000, Stirling's series above."""
    if n < 2000:
        return Decimal(math.factorial(n)).ln() if n > 1 else Decimal(0)
    x = Decimal(n)
    value = (x + Decimal("0.5")) * x.ln() - x + (2 * PI).ln() / 2
    for m, b in enumerate(BERNOULLI, 1):
        value += (
            Decimal(b.numerator)
            / Decimal(b.denominator)
            / (2 * m * (2 * m - 1) * x ** (2 * m - 1))
        )
    return value


def binomial_reference(k, n, p):
    """The binomial pmf of the double p (failure 1 - p exactly), 60 digits."""
    p = Decimal(p)
    log = (
        log_factorial(n)
        - log_factorial(k)
        - log_factorial(n - k)
        + k * p.ln()
        + (n - k) * (1 - p).ln()
    )
    return float(log.exp())


def special_functions():
    """Binomial pmf against SciPy (n <= 1e5) and, where SciPy and the native
    kernel disagree, against a 60-digit Decimal value; Poisson tail against
    SciPy."""
    worst_scipy, worst_exact, disputed = 0.0, 0.0, 0
    scipy_above = {}
    for p in (1e-9, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.999, 1 - 1e-9):
        for n in [*range(1, 61), 100, 1000, 10**4, 10**5, 10**6, 2**40]:
            if n <= 1000:
                ks = np.arange(n + 1)
            else:
                sd = np.sqrt(n * p * (1 - p))
                ks = np.unique(
                    np.clip(
                        np.round(n * p + np.linspace(-10, 10, 201) * sd), 0, n
                    ).astype(np.int64)
                )
                ks = np.unique(np.r_[0, 1, ks, n - 1, n])
            want = binom.pmf(ks, n, p)
            got = np.array([_native._binomial_pmf(float(k), float(n), p) for k in ks])
            ok = want >= 1e-280
            error = np.abs(got[ok] - want[ok]) / want[ok]
            if n <= 10**5 and ok.any():
                worst_scipy = max(worst_scipy, float(error.max()))
            for i in np.flatnonzero(error > 1e-12):
                k = int(ks[ok][i])
                exact = binomial_reference(k, n, p)
                worst_exact = max(worst_exact, float(abs(got[ok][i] - exact) / exact))
                disputed += 1
                scipy_above[n] = max(scipy_above.get(n, 0.0), float(error[i]))
    require(worst_scipy <= 1e-12, f"binomial pmf against SciPy {worst_scipy}")
    require(worst_exact <= 1e-13, f"binomial pmf against Decimal {worst_exact}")
    tail_worst = 0.0
    for mu in (1e-12, 1e-3, 0.5, 1, 4, 8, 15.9, 16, 16.0001, 40, 100, 1000):
        for k in range(int(mu) + 201):
            want = float(gammainc(k + 1, mu))
            if want >= 1e-280:
                got = _native._poisson_tail(k + 1, float(mu))
                tail_worst = max(tail_worst, abs(got - want) / want)
    require(tail_worst <= 1e-12, f"Poisson tail relative error {tail_worst}")
    return {
        "binomial_against_scipy_n_up_to_1e5_within_1e-12": worst_scipy <= 1e-12,
        "binomial_rows_where_scipy_differs_by_more_than_1e-12": disputed,
        "binomial_against_decimal_on_those_rows_within_1e-13": bool(
            worst_exact <= 1e-13
        ),
        "scipy_disagreement_by_n": {str(n): v > 1e-12 for n, v in scipy_above.items()},
        "poisson_tail_against_scipy_within_1e-12": tail_worst <= 1e-12,
    }


def model(name, species, parameters, reactions, stop=None):
    m = crn.Model(name)
    m.add_species(species).add_parameter(parameters)
    for rate, reactants, products, label in reactions:
        m.add_reaction(rate, reactants, products, name=label)
    if stop is not None:
        m.set_stop(stop)
    return m


def outcome(call):
    try:
        call()
    except Exception as error:  # noqa: BLE001 - the type is the result
        return type(error).__name__, str(error)
    return None


def compare(case, prepare, solves, bound_kind):
    """Prepare with both implementations and compare every solve."""
    start = perf_counter()
    native = prepare(crn.prepare_fsp)
    native_seconds = perf_counter() - start
    start = perf_counter()
    ref = prepare(reference.prepare_fsp)
    reference_seconds = perf_counter() - start
    row = {"states": len(native.states)}
    require(
        native.states.dtype == ref.states.dtype
        and np.array_equal(native.states, ref.states),
        f"{case}: states differ",
    )
    require(native._id == ref._id, f"{case}: prepared ids differ")
    for label, kwargs in solves:
        tol = kwargs.get("tolerance", 1e-8)
        kind = bound_kind if "signals" in kwargs else "uniform"
        bound = 1e-12 + (1 if kind == "uniform" else 2) * tol
        start = perf_counter()
        a = native.solve(**kwargs)
        native_seconds += perf_counter() - start
        start = perf_counter()
        b = ref.solve(**kwargs)
        reference_seconds += perf_counter() - start
        require(
            len(a.snapshots) == len(b.snapshots)
            and a.solver == b.solver
            and all(
                (x.time, x.actions_applied) == (y.time, y.actions_applied)
                for x, y in zip(a.snapshots, b.snapshots, strict=True)
            ),
            f"{case}/{label}: snapshot grid or solver differs",
        )
        distance = max(
            float(
                np.abs(x.probabilities - y.probabilities).sum()
                + abs(x.physical_leakage - y.physical_leakage)
            )
            for x, y in zip(a.snapshots, b.snapshots, strict=True)
        )
        require(distance <= bound, f"{case}/{label}: L1 distance {distance} > {bound}")
        tails = (a.snapshots[-1].series_tail_bound, b.snapshots[-1].series_tail_bound)
        require(max(tails) <= tol, f"{case}/{label}: tail bound above tolerance")
        if kind == "uniform":
            require(
                abs(a.work - b.work) <= 0.001 * b.work + 2,
                f"{case}/{label}: work {a.work} against {b.work}",
            )
        row[label] = {
            "solver": a.solver,
            "within_bound": distance <= bound,
            "leakage": float(a.snapshots[-1].physical_leakage),
            "work_equal": a.work == b.work,
        }
    print(
        f"{case}: native {native_seconds:.2f} s, reference {reference_seconds:.2f} s",
        file=sys.stderr,
    )
    return row


def same_error(case, call):
    """Both implementations fail with the same type and message."""
    a, b = (
        outcome(lambda: call(crn.prepare_fsp)),
        outcome(lambda: call(reference.prepare_fsp)),
    )
    require(a is not None and a == b, f"{case}: {a} against {b}")
    return a[0] if a else None


def birth_death(birth=2.0, death=1.0, start=0):
    return model(
        "birth_death",
        {"X": start},
        {"b": birth, "d": death},
        [
            (crn.mass_action("b"), [], ["X"], "birth"),
            (crn.mass_action("d"), ["X"], [], "death"),
        ],
    )


def cases():
    rows = {}
    # C1 state-dependent (non-mass-action) rate and a stop condition.
    feedback = model(
        "feedback",
        {"X": 0, "Y": 0},
        {"k": 4.0, "d": 0.5, "e": 1.5},
        [
            ("k / (1 + X)", [], ["X"], "make"),
            (crn.mass_action("d"), ["X"], [], "loss"),
            (crn.mass_action("e"), ["X"], ["X", "Y"], "mark"),
        ],
        stop="Y >= 5",
    )
    rows["C1_state_dependent_stop"] = compare(
        "C1",
        lambda prep: prep(feedback, {"X": 30, "Y": 6}),
        [("uniform", {"end": 3, "report_times": [0.5, 1, 3], "tolerance": 1e-10})],
        "uniform",
    )
    # C2 an initial law with several states.
    rows["C2_initial_mixture"] = compare(
        "C2",
        lambda prep: prep(
            birth_death(), {"X": 40}, initial={(0,): 0.2, (3,): 0.5, (7,): 0.3}
        ),
        [("uniform", {"end": 2, "report_times": [0, 1, 2], "tolerance": 1e-12})],
        "uniform",
    )
    # C3 interventions of every kind at 0, inside and at the horizon; C4
    # escape through the bounds.
    pair = model(
        "pair",
        {"X": 2, "Y": 0},
        {"b": 3.0, "d": 0.7, "c": 0.4},
        [
            (crn.mass_action("b"), [], ["X"], "birth"),
            (crn.mass_action("d"), ["X"], [], "death"),
            (crn.mass_action("c"), ["X"], ["X", "Y"], "copy"),
        ],
    )
    policy = crn.Protocol()
    policy.at(0, set_counts={"X": 5})
    policy.at(0.5, change={"X": 3}, partition={"Y": 0.3})
    policy.at(1.0, partition={"X": 0.0}, change={"Y": 2})
    policy.at(1.5, partition={"X": 1.0}, set_counts={"Y": 30})
    policy.at(2.0, partition={"X": 0.6}, parameters={"b": 1.0})
    rows["C3_C4_interventions_escape"] = compare(
        "C3",
        lambda prep: prep(pair, {"X": 14, "Y": 12}, protocol=policy),
        [
            (
                "uniform",
                {"end": 2, "report_times": [0, 0.5, 1, 1.5, 2], "tolerance": 1e-12},
            )
        ],
        "uniform",
    )
    leak = crn.prepare_fsp(pair, {"X": 14, "Y": 12}, protocol=policy).solve(
        2, [2], tolerance=1e-12
    )
    require(leak.snapshots[-1].physical_leakage > 1e-3, "C4: no escape")
    removal = crn.Protocol().at(0.5, change={"X": -6})
    rows["C3_removal_error"] = same_error(
        "C3 removal",
        lambda prep: prep(pair, {"X": 14, "Y": 12}, protocol=removal).solve(
            1, [1], tolerance=1e-12
        ),
    )
    # C5-C7 signals: a step table with a discontinuity inside the horizon,
    # a linear table and a callable, by the midpoint rule and BDF.
    decaying = birth_death()
    step = crn.Signal.table([0, 1.25], [1.0, 6.0], interpolation="step")
    linear = crn.Signal.table([0, 0.7, 2], [1.0, 5.0, 2.0])
    wave = crn.Signal(lambda t: 2 + np.sin(3 * t))
    grid = {"end": 2, "report_times": [0.5, 1.25, 1.5, 2]}
    for solver in ("midpoint", "bdf"):
        rows[f"C5_C7_signals_{solver}"] = compare(
            f"C5-C7 {solver}",
            lambda prep: prep(decaying, {"X": 40}),
            [
                (
                    name,
                    {
                        **grid,
                        "signals": {"birth": signal},
                        "solver": solver,
                        "tolerance": 1e-8,
                    },
                )
                for name, signal in (("step", step), ("linear", linear), ("call", wave))
                if solver == "midpoint" or name != "linear"
            ],
            solver,
        )
    # C8 parameter rebinding; C8b a binding that changes the initial counts.
    seeded = {
        "schema": "crn-craft.model/1.0",
        "name": "seeded",
        "parameters": {"b": 2.0, "d": 1.0, "x0": 3},
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": "x0"}],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {
                "name": "birth",
                "reactants": [],
                "rate": "b",
                "products": [{"type": "X"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "X"}],
                "rate": "d",
                "products": [],
            },
        ],
    }
    rows["C8_rebinding"] = compare(
        "C8",
        lambda prep: prep(seeded, {"X": 40}),
        [
            (
                "bound",
                {
                    "end": 2,
                    "report_times": [1, 2],
                    "parameters": {"b": 3.5, "d": 0.5},
                    "tolerance": 1e-12,
                },
            )
        ],
        "uniform",
    )
    rows["C8b_topology_refused"] = same_error(
        "C8b",
        lambda prep: prep(seeded, {"X": 40}).solve(1, [1], parameters={"x0": 4}),
    )
    # C9 stiff: rate * duration / 16 above 1000 pieces.
    stiff = birth_death(200.0, 1000.0, 0)
    rows["C9_stiff"] = compare(
        "C9",
        lambda prep: prep(stiff, {"X": 20}),
        [
            (
                "uniform",
                {"end": 1, "report_times": [1], "tolerance": 1e-9, "max_work": 10**7},
            )
        ],
        "uniform",
    )
    # C11 budgets and inputs.
    errors = {
        "state_budget": lambda prep: prep(birth_death(), {"X": 100}, max_states=5),
        "edge_budget": lambda prep: prep(birth_death(), {"X": 100}, max_edges=3),
        "work_budget": lambda prep: prep(birth_death(), {"X": 40}).solve(
            2, [2], max_work=10
        ),
        "pieces_budget": lambda prep: prep(stiff, {"X": 20}).solve(
            1, [1], max_work=100
        ),
        "expansion_budget": lambda prep: prep(
            birth_death(start=10**6),
            {"X": 2 * 10**6},
            protocol=crn.Protocol().at(1, partition={"X": 0.5}),
            max_edges=1000,
        ),
        "outside_initial": lambda prep: prep(birth_death(start=50), {"X": 40}),
        "unprepared_initial": lambda prep: prep(birth_death(), {"X": 40}).solve(
            1, [1], initial={(41,): 1.0}
        ),
        "bad_law": lambda prep: prep(birth_death(), {"X": 40}, initial={(1,): 0.5}),
        "state_dependent_signal": lambda prep: prep(feedback, {"X": 10, "Y": 6}).solve(
            1, [1], signals={"make": wave}, solver="midpoint"
        ),
        "unknown_signal": lambda prep: prep(birth_death(), {"X": 10}).solve(
            1, [1], signals={"nothing": wave}, solver="midpoint"
        ),
        "negative_signal": lambda prep: prep(birth_death(), {"X": 10}).solve(
            1, [1], signals={"birth": crn.Signal(lambda t: -1.0)}, solver="midpoint"
        ),
        "raising_signal": lambda prep: prep(birth_death(), {"X": 10}).solve(
            1, [1], signals={"birth": crn.Signal(lambda t: 1 / 0)}, solver="midpoint"
        ),
        "midpoint_rounding": lambda prep: prep(birth_death(), {"X": 48}).solve(
            2, [2], signals={"birth": wave}, solver="midpoint", tolerance=1e-14
        ),
        "midpoint_work": lambda prep: prep(birth_death(), {"X": 48}).solve(
            2, [2], signals={"birth": wave}, solver="midpoint", max_work=5
        ),
        "kick_work": lambda prep: prep(pair, {"X": 14, "Y": 12}, protocol=policy).solve(
            2, [2], max_work=3
        ),
    }
    rows["C11_errors"] = {name: same_error(name, call) for name, call in errors.items()}
    # C12 a second-order reaction.
    dimer = model(
        "dimer",
        {"A": 30, "B": 0},
        {"k": 0.05, "u": 0.3},
        [
            (crn.mass_action("k"), ["A", "A"], ["B"], "bind"),
            (crn.mass_action("u"), ["B"], ["A", "A"], "unbind"),
        ],
    )
    rows["C12_second_order"] = compare(
        "C12",
        lambda prep: prep(dimer, {"A": 30, "B": 15}),
        [("uniform", {"end": 3, "report_times": [1, 3], "tolerance": 1e-12})],
        "uniform",
    )
    return rows


def analytic():
    """Laws with a closed form, compared with nothing but that form."""
    rows = {}

    def poisson(mean, top):
        return np.array(
            [
                math.exp(i * math.log(mean) - mean - math.lgamma(i + 1))
                for i in range(top + 1)
            ]
        )

    def gate(name, graph, result, expected, bound):
        worst = 0.0
        for snapshot, law in zip(result.snapshots, expected, strict=True):
            found = np.bincount(
                graph.states[:, 0], snapshot.probabilities, minlength=len(law)
            )
            worst = max(worst, float(np.abs(found - law).max()))
        require(worst <= bound, f"analytic {name}: gap {worst:.3e} above {bound:.0e}")
        rows[name] = {"gap": worst, "bound": bound}

    # Thirty molecules dying at rate k: binomial survival.
    death = model(
        "death", {"X": 30}, {"k": 0.3}, [(crn.mass_action("k"), ["X"], [], "die")]
    )
    graph = crn.prepare_fsp(death, {"X": 30})
    times = [0.5, 2.0]
    keep = [math.exp(-0.3 * t) for t in times]
    gate(
        "binomial_survival",
        graph,
        graph.solve(2.0, times, tolerance=1e-12),
        [
            np.array([math.comb(30, i) * q**i * (1 - q) ** (30 - i) for i in range(31)])
            for q in keep
        ],
        1e-11,
    )

    # Immigration and death from zero: Poisson with the mean of m' = b(t) - k m.
    a, w, k, top = 4.0, 1.3, 0.7, 80
    queue = model(
        "queue",
        {"X": 0},
        {"b": a, "k": k},
        [
            (crn.mass_action("b"), [], ["X"], "in"),
            (crn.mass_action("k"), ["X"], [], "out"),
        ],
    )
    graph = crn.prepare_fsp(queue, {"X": top})

    def wave_mean(t):
        decay = math.exp(-k * t)
        swing = k * math.sin(w * t) - w * math.cos(w * t) + w * decay
        return a / k * (1 - decay) + a * swing / (k * k + w * w)

    def step_mean(t):
        if t <= 1:
            return 5 / k * (1 - math.exp(-k * t))
        at_jump = 5 / k * (1 - math.exp(-k))
        return at_jump * math.exp(-k * (t - 1)) + 0.5 / k * (1 - math.exp(-k * (t - 1)))

    wave = crn.Signal(lambda t: a * (1 + math.sin(w * t)))
    step = crn.Signal.table([0.0, 1.0], [5.0, 0.5], interpolation="step")
    for solver in ("midpoint", "bdf"):
        for name, signal, end, times, mean in (
            ("wave", wave, 5.0, [0.7, 2.0, 5.0], wave_mean),
            ("step", step, 3.0, [1.0, 1.5, 3.0], step_mean),
        ):
            result = graph.solve(
                end,
                times,
                signals={"in": signal},
                # The second-order midpoint rule refuses goals near 1e-9
                # on these spans (its step error goal reaches rounding).
                tolerance=1e-9 if solver == "bdf" else 1e-7,
                max_step=0.05,
                solver=solver,
                max_work=10**8,
            )
            gate(
                f"poisson_{name}_{solver}",
                graph,
                result,
                [poisson(mean(t), top) for t in times],
                2e-6,
            )
            if solver == "bdf":
                gate(
                    f"poisson_{name}_bdf_tight",
                    graph,
                    result,
                    [poisson(mean(t), top) for t in times],
                    1e-8,
                )

    # Birth 2 + sin(3 t), death 1 from 0 (case C7): an L1 gap of the native
    # BDF to the exact Poisson law of at most 10 * tolerance.
    decaying = birth_death()
    graph = crn.prepare_fsp(decaying, {"X": 40})
    times = [0.5, 1.25, 1.5, 2.0]

    def call_mean(t):
        return (
            2 * (1 - math.exp(-t))
            + (math.sin(3 * t) - 3 * math.cos(3 * t) + 3 * math.exp(-t)) / 10
        )

    for tol in (1e-6, 1e-8):
        result = graph.solve(
            2.0,
            times,
            signals={"birth": crn.Signal(lambda t: 2 + math.sin(3 * t))},
            tolerance=tol,
            max_work=10**8,
        )
        worst = max(
            float(
                np.abs(
                    np.bincount(graph.states[:, 0], s.probabilities, minlength=41)
                    - poisson(call_mean(t), 40)
                ).sum()
            )
            for s, t in zip(result.snapshots, times, strict=True)
        )
        require(
            worst <= 10 * tol,
            f"analytic poisson_call_bdf {tol:.0e}: L1 gap {worst:.3e}",
        )
        rows[f"poisson_call_bdf_{tol:.0e}"] = {"l1_gap": worst, "bound": 10 * tol}

    # Stationary laws: 400 molecules exchanging between two forms are
    # binomial; immigration and death with a return state are Poisson.
    flip = model(
        "flip",
        {"A": 400, "B": 0},
        {"f": 0.7, "r": 0.2},
        [
            (crn.mass_action("f"), ["A"], ["B"], "ab"),
            (crn.mass_action("r"), ["B"], ["A"], "ba"),
        ],
    )
    law = crn.stationary_fsp(flip, {"A": 400, "B": 400})
    share = 0.2 / 0.9
    binomial = np.exp(
        [
            math.lgamma(401)
            - math.lgamma(i + 1)
            - math.lgamma(401 - i)
            + i * math.log(share)
            + (400 - i) * math.log(1 - share)
            for i in range(401)
        ]
    )
    still = crn.stationary_fsp(queue, {"X": top}, return_state=(0,))
    for name, found, expected in (
        (
            "stationary_binomial",
            np.bincount(law.states[:, 0], law.probabilities, minlength=401),
            binomial,
        ),
        (
            "stationary_poisson",
            still.probabilities[np.argsort(still.states[:, 0])],
            poisson(a / k, top),
        ),
    ):
        worst = float(np.abs(found - expected).max())
        require(worst <= 1e-12, f"analytic {name}: gap {worst:.3e} above 1e-12")
        rows[name] = {"gap": worst, "bound": 1e-12}

    # Keeping each molecule with probability 0.3 at t = 1 thins the Poisson law.
    thinned = crn.prepare_fsp(
        queue, {"X": top}, protocol=crn.Protocol().at(1.0, partition={"X": 0.3})
    )
    after = 0.3 * a / k * (1 - math.exp(-k))
    gate(
        "poisson_thinning",
        thinned,
        thinned.solve(2.0, [1.0, 2.0], tolerance=1e-12, max_work=10**8),
        [
            poisson(after, top),
            poisson(after * math.exp(-k) + a / k * (1 - math.exp(-k)), top),
        ],
        1e-9,
    )
    return rows


def stationary_cases():
    """Stationary laws: native closed classes and sparse LU against SciPy."""
    rows = {}

    def switch(name, rates, start):
        """One molecule moving among slots S0..S5 with first-order rates."""
        slots = [f"S{i}" for i in range(6)]
        m = crn.Model(name)
        m.add_species({s: int(i == start) for i, s in enumerate(slots)})
        for (a, b), rate in rates.items():
            m.add_parameter(f"k{a}{b}", rate)
            m.add_reaction(
                crn.mass_action(f"k{a}{b}"), [slots[a]], [slots[b]], name=f"k{a}{b}"
            )
        return m, dict.fromkeys(slots, 1)

    def state(i):
        return tuple(int(j == i) for j in range(6))

    def run(case, call):
        a = outcome_value(lambda: call(crn.stationary_fsp))
        b = outcome_value(lambda: call(reference.stationary_fsp))
        if isinstance(a, tuple) or isinstance(b, tuple):
            brief = [x if isinstance(x, tuple) else "a result" for x in (a, b)]
            require(a == b, f"{case}: {brief[0]} against {brief[1]}")
            rows[case] = a[0] if isinstance(a, tuple) else "no error"
            return
        same = (
            np.array_equal(a.states, b.states)
            and a.method == b.method
            and a.converged == b.converged
            and len(a.history) == len(b.history)
            and all(
                x["bounds"] == y["bounds"] and x["states"] == y["states"]
                for x, y in zip(a.history, b.history, strict=True)
            )
        )
        gap = float(np.abs(a.probabilities - b.probabilities).max())
        require(same, f"{case}: results differ in shape or history")
        require(gap <= 1e-10, f"{case}: probability gap {gap} above 1e-10")
        rows[case] = {
            "states": len(a.states),
            "within_bound": gap <= 1e-10,
            "converged": a.converged,
            "residual_below_1e-12": a.residual_l1 <= 1e-12,
        }

    # S2 two closed classes ({S0, S1}, {S3, S4}) and transient S2, S5.
    reducible, slots = switch(
        "reducible",
        {
            (0, 1): 2.0,
            (1, 0): 1.0,
            (5, 2): 0.5,
            (2, 0): 1.0,
            (2, 3): 1.0,
            (3, 4): 1.0,
            (4, 3): 3.0,
        },
        5,
    )
    run("S2_two_classes_refused", lambda f: f(reducible, slots))
    run("S2_first_class", lambda f: f(reducible, slots, class_state=state(1)))
    run("S2_second_class", lambda f: f(reducible, slots, class_state=state(3)))
    run("S2_transient_refused", lambda f: f(reducible, slots, class_state=state(2)))
    # S3 an absorbing state.
    absorbing, slots = switch(
        "absorbing", {(0, 1): 1.0, (1, 0): 2.0, (1, 2): 0.5, (2, 3): 4.0}, 0
    )
    run("S3_absorbing", lambda f: f(absorbing, slots))
    # S4 nearly decomposable and stiff: fast blocks {S0, S1}, {S2, S3, S4}.
    stiff, slots = switch(
        "nearly_decomposable",
        {
            (0, 1): 1e3,
            (1, 0): 2e3,
            (1, 2): 1e-7,
            (2, 1): 3e-7,
            (2, 3): 3e3,
            (3, 2): 1e3,
            (3, 4): 5e2,
            (4, 3): 2e3,
        },
        0,
    )
    run("S4_nearly_decomposable", lambda f: f(stiff, slots))
    # S5 a redirected birth-death chain with face expansions.
    chain = birth_death(2.0, 1.0, 0)
    run(
        "S5_redirect_expansions",
        lambda f: f(chain, {"X": 3}, return_state=(0,), expansions=3, tolerance=1e-4),
    )
    run(
        "S5_return_outside",
        lambda f: f(chain, {"X": 3}, return_state=(9,)),
    )
    run("S6_escapes_refused", lambda f: f(chain, {"X": 3}))
    stopper = model(
        "stopper",
        {"X": 0},
        {"b": 1.0},
        [(crn.mass_action("b"), [], ["X"], "birth")],
        stop="X >= 2",
    )
    run("S6_stop_refused", lambda f: f(stopper, {"X": 4}, return_state=(0,)))
    # The factorization budget (native only): a refusal that names it.
    graph = crn.prepare_fsp(chain, {"X": 30})
    refusals = (
        outcome(
            lambda: graph.solve(
                1,
                [1],
                signals={"birth": crn.Signal(lambda t: 1 + t)},
                max_factor_bytes=999,
            )
        ),
        outcome(
            lambda: crn.stationary_fsp(
                chain, {"X": 30}, return_state=(0,), max_factor_bytes=999
            )
        ),
    )
    named = all(
        r is not None and r[0] == "RuntimeError" and "max_factor_bytes=999" in r[1]
        for r in refusals
    )
    require(named, f"factorization budget refusals: {refusals}")
    rows["S7_factorization_budget_refused"] = named
    return rows


def outcome_value(call):
    """The result of a call, or the type and message of its exception."""
    try:
        return call()
    except Exception as error:  # noqa: BLE001 - the type is the result
        return type(error).__name__, str(error)


def yeast(full):
    """C10: the total-RNA model of tutorial O3 at one dose, to minute
    10 (the full 55 minutes with ``full``)."""
    folder = ROOT / "examples" / "yeast_rna_fit"
    sys.path.insert(0, str(folder))
    from fitting import Data, observations
    from model import KEYS, hog1, parameters
    from native_fit import Objective

    anchors = json.loads((folder / "source-parameters.json").read_text())
    x = parameters(anchors["parameters"]["NonSpatial-Distributions"])
    data = Data(observations(), spatial_model=False)
    objective = Objective(x, data, (360,), 1e-6, 60)
    graph = objective.graph
    document = graph._document
    bounds = dict(zip(graph.species, graph._bounds, strict=True))
    dose = data.dose[0]

    def rate(t):
        return max(0.0, x[1] - x[2] * float(hog1(t, dose, x[12])))

    minutes = objective.minutes if full else objective.minutes[objective.minutes <= 10]
    solve = {
        "end": float(minutes[-1] * 60),
        "report_times": minutes * 60,
        "parameters": dict(zip(KEYS, map(float, x), strict=True)),
        "signals": {"k21": crn.Signal(rate, change_times=(x[12],))},
        "tolerance": 1e-6,
        "max_step": 60,
        "solver": "midpoint",
        "max_work": 10**9,
    }
    rows = compare(
        "C10",
        lambda prep: prep(document, bounds),
        [("midpoint", solve)],
        "midpoint",
    )
    # C10b: the same segment by the native NDF and SciPy's BDF.
    rows.update(
        compare(
            "C10b",
            lambda prep: prep(document, bounds),
            [("bdf", {**solve, "solver": "bdf"})],
            "bdf",
        )
    )
    return rows


def sensitivities():
    """Forward sensitivities against closed forms and centered differences."""
    from crn_craft.sensitivities import forward_sensitivities

    rows = {}
    graph = crn.prepare_fsp(birth_death(2.0, 0.3), {"X": 45})
    result = forward_sensitivities(
        graph, [(4.0, {"b": 2.0, "d": 0.3})], [1.0, 4.0], ["b"], tolerance=1e-10
    )
    worst = 0.0
    for row, at in enumerate((1.0, 4.0)):
        slope = -math.expm1(-0.3 * at) / 0.3
        mean = 2.0 * slope
        for i, (x,) in enumerate(graph.states):
            p = math.exp(x * math.log(mean) - mean - math.lgamma(x + 1))
            exact = p * (x / mean - 1.0) * slope
            worst = max(worst, abs(result.derivatives[row, 0, i] - exact))
    require(worst <= 1e-7, f"sensitivity birth-death: gap {worst:.3e}")
    rows["birth_death_poisson"] = {"gap": worst, "bound": 1e-7}

    hill = crn.Model("hill_switch")
    hill.add_species({"x": 0, "y": 0}).add_parameter(
        {"b": 0.5, "k": 6.0, "a": 0.3, "c": 2.0, "d": 1.0}
    )
    hill.add_observable("cy", "y")
    hill.add_reaction("b+k/(1+a*cy^2)", [], ["x"], name="repressed")
    hill.add_reaction(crn.mass_action("c"), [], ["y"], name="switched")
    hill.add_reaction(crn.mass_action("d"), ["x"], [], name="loss_x")
    hill.add_reaction(crn.mass_action("d"), ["y"], [], name="loss_y")
    graph = crn.prepare_fsp(hill, {"x": 30, "y": 25})
    base = {"b": 0.5, "k": 6.0, "a": 0.3, "c": 2.0, "d": 1.0}
    segments = [(1.0, base), (2.0, base | {"c": 6.0})]
    names = ("a", "c")
    result = forward_sensitivities(graph, segments, [1.0, 2.0], names, tolerance=1e-10)
    mass = float(np.abs(result.probabilities.sum(axis=1) - 1).max())
    drift = float(np.abs(result.derivatives.sum(axis=2)).max())
    require(mass <= 1e-12 and drift <= 1e-9, f"sensitivity mass {mass} {drift}")

    def laws(name, delta):
        initial, out = None, []
        for (end, binding), start in zip(segments, (0.0, 1.0), strict=True):
            binding = binding | {name: binding[name] + delta}
            solved = graph.solve(
                end - start,
                [end - start],
                initial=initial,
                parameters=binding,
                tolerance=1e-13,
            )
            p = np.maximum(solved.snapshots[-1].probabilities, 0.0)
            out.append(p)
            initial = {
                tuple(int(n) for n in s): q
                for s, q in zip(graph.states, p / p.sum(), strict=True)
                if q > 0
            }
        return np.array(out)

    for j, name in enumerate(names):
        gaps = []
        for step in (1e-3, 5e-4):
            h = step * base[name]
            difference = (laws(name, h) - laws(name, -h)) / (2 * h)
            gaps.append(float(np.abs(result.derivatives[:, j, :-1] - difference).max()))
        ratio = gaps[0] / gaps[1]
        require(
            max(gaps) <= 1e-6 and (3 <= ratio <= 5 or max(gaps) < 5e-8),
            f"sensitivity {name}: gaps {gaps}",
        )
        rows[f"hill_switch_{name}"] = {"gaps": gaps, "bound": 1e-6, "ratio": ratio}
    rows["hill_switch_mass"] = {"mass": mass, "derivative_sum": drift}

    full = crn.prepare_fsp(birth_death(1.0, 0.5, start=1), {"X": 12})
    box = crn.prepare_fsp(
        birth_death(1.0, 0.5, start=1),
        {"X": 12},
        retained=[(x,) for x in range(13)],
    )
    same_states = bool(np.array_equal(full.states, box.states))
    same_id = full._id == box._id
    one, two = (g.solve(3.0, [3.0]).snapshots[0] for g in (full, box))
    same_law = bool(np.array_equal(one.probabilities, two.probabilities))
    cut = (
        crn.prepare_fsp(
            birth_death(1.0, 0.5, start=1),
            {"X": 12},
            retained=[(x,) for x in range(4)],
        )
        .solve(3.0, [3.0], tolerance=1e-12)
        .snapshots[0]
    )
    total = float(cut.probabilities.sum() + cut.physical_leakage)
    require(
        same_states
        and same_id
        and same_law
        and abs(total - 1) <= 2e-12
        and cut.physical_leakage > 1e-3,
        "retained set: full box differs or the cut loses mass",
    )
    rows["retained"] = {
        "full_box_identical": same_states and same_id and same_law,
        "cut_sink": float(cut.physical_leakage),
        "cut_total": total,
    }
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--full-yeast", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "special": special_functions(),
        "cases": cases(),
        "analytic": analytic(),
        "stationary": stationary_cases(),
        "sensitivities": sensitivities(),
    }
    if not args.quick:
        result["cases"]["C10_yeast_total_rna"] = yeast(args.full_yeast)
    result["failures"] = FAILURES
    text = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    if FAILURES:
        raise SystemExit("check_projection_native: " + "; ".join(FAILURES))


if __name__ == "__main__":
    main()
