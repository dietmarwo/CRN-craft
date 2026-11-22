"""Fit 36 parameters to 16 reconstructed three-step pathway trajectories."""

from itertools import product

import numpy as np
from scipy.integrate import solve_ivp

from common import digest, require, write
from .concentrations import Solver, document
from .search import PACKAGES as PACKAGES
from .search import PENALTY, check, options_resume, save_plot, search

DOI = "10.1101/gr.1262503"
STATES = ("G1", "G2", "G3", "E1", "E2", "E3", "M1", "M2")
INITIAL = (0.66667, 0.57254, 0.41758, 0.4, 0.36409, 0.29457, 1.419, 0.93464)
NAMES = (
    *(f"{name}{i}" for i in (1, 2, 3) for name in ("V", "Ki", "ni", "Ka", "na", "k")),
    *(f"{name}{i}" for i in (4, 5, 6) for name in ("V", "K", "k")),
    *(name for i in (1, 2, 3) for name in (f"kcat{i}", f"Km{2 * i - 1}", f"Km{2 * i}")),
)
TRUTH = np.array([1.0, 1.0, 2.0, 1.0, 2.0, 1.0] * 3 + [0.1, 1.0, 0.1] * 3 + [1.0] * 9)
CONDITIONS = tuple(product((0.1, 0.46416, 2.1544, 10.0), (0.05, 0.13572, 0.36840, 1.0)))
BROAD = np.array(
    [(-1.0, 1.0) if n.startswith(("ni", "na")) else (-12.0, 12.0) for n in NAMES]
)
# A free-start, local-scale comparison distinct from the original large box.
FOCUSED = np.column_stack((np.log10(TRUTH) - 1, np.log10(TRUTH) + 1))
FOCUSED = np.clip(FOCUSED, BROAD[:, 0, None], BROAD[:, 1, None])
PUBLISHED = 0.0013


def model(literal):
    regulators = ("S", "M1", "M2")
    # After multiplying numerator and denominator by (regulator/Ka)^na,
    # the activation term is u / (1 + (1 + inhibition)*u).
    defs = [(f"sat{i}", f"({r}/Ka{i})^na{i}") for i, r in enumerate(regulators, 1)]
    genes = tuple(
        f"V{i} * sat{i} / (1 + (1 + (P/Ki{i})^ni{i})*sat{i}) - k{i}*G{i}"
        for i in (1, 2, 3)
    )
    enzymes = tuple(
        f"V{i + 3}*G{2 if literal and i == 3 else i}/(K{i + 3}+G{i})-k{i + 3}*E{i}"
        for i in (1, 2, 3)
    )
    left, right = ("S", "M1", "M2"), ("M1", "M2", "P")
    defs.extend(
        (
            f"flux{i}",
            f"kcat{i}*E{i}*({a}-{b})/(Km{2 * i - 1}+{a}+{b}*Km{2 * i - 1}/Km{2 * i})",
        )
        for i, (a, b) in enumerate(zip(left, right, strict=True), 1)
    )
    return document(
        STATES,
        (*NAMES, "S", "P"),
        tuple(map(str, INITIAL)),
        (*genes, *enzymes, "flux1-flux2", "flux2-flux3"),
        assignments=defs,
    )


def reference(parameters, times, literal=False):
    theta = np.asarray(parameters)
    tx = theta[:18].reshape(3, 6)
    tl = theta[18:27].reshape(3, 3)
    cat = theta[27:].reshape(3, 3)
    results = []
    for a, b in CONDITIONS:

        def dynamics(t, concentrations, a=a, b=b):
            x = np.maximum(concentrations, 0)
            g, e = x[:3], x[3:6]
            inputs = np.array([a, *x[-2:]])
            with np.errstate(divide="ignore", over="ignore"):
                source = tx[:, 0] / (
                    1 + (b / tx[:, 1]) ** tx[:, 2] + (tx[:, 3] / inputs) ** tx[:, 4]
                )
            transcribed = source - tx[:, 5] * g
            translated = (
                tl[:, 0] * (g[[0, 1, 1]] if literal else g) / (tl[:, 1] + g)
                - tl[:, 2] * e
            )
            products = np.array([*x[-2:], b])
            flux = (
                cat[:, 0]
                * e
                * (inputs - products)
                / cat[:, 1]
                / (1 + inputs / cat[:, 1] + products / cat[:, 2])
            )
            return np.r_[transcribed, translated, flux[:-1] - flux[1:]]

        solved = solve_ivp(
            dynamics,
            [0.0, times[-1]],
            INITIAL,
            t_eval=times,
            method="LSODA",
            atol=1e-12,
            rtol=2e-11,
        )
        require(solved.success, solved.message)
        results.append(solved.y.T)
    return np.asarray(results)


def weighted_error(observations, prediction, scales):
    return float(np.square((np.asarray(prediction) - observations) / scales).sum())


def options(parser):
    parser.add_argument(
        "--evaluations", type=int, help="Final nested evaluation checkpoint."
    )
    parser.add_argument("--retries", type=int)
    parser.add_argument("--broad-only", action="store_true")
    parser.add_argument("--printed-equation-12", action="store_true")
    options_resume(parser)


def run(args, report):
    """Reconstruct the weighted least-squares problem and compare DE with CMA-ES."""
    times = np.arange(0.0, 121.0, 6.0)
    literal = args.printed_equation_12
    observed = reference(TRUTH, times, literal)
    scales = np.max(np.abs(observed[:, 1:]), axis=1, keepdims=True)
    require(np.min(scales) > 0, "all observed states need a positive weight scale")
    write(
        args.output / "reconstruction.json",
        dict(
            doi=DOI,
            times=times,
            observations=observed,
            scales=scales,
            parameters=dict(zip(NAMES, TRUTH, strict=True)),
            conditions=CONDITIONS,
            initial=INITIAL,
            initial_and_grid_source="Müller et al. (2008), Table 1 and Data generation",
            variant="printed-G2" if literal else "G3 translation",
        ),
    )
    context = digest(args.output / "reconstruction.json")
    fast = args.preset == "smoke"
    budget = args.evaluations or (24 if fast else 200_000)
    retries = args.retries or (1 if fast else 4)
    checkpoints = (2000, 20_000, 200_000)
    results = []
    with Solver(args.ode_binary, args.output / "model.json", model(literal)) as solver:

        def solve(p, t=times, tight=False):
            y, _ = solver.batch(
                [(*p, a, b) for a, b in CONDITIONS],
                t,
                **(dict(rtol=3e-10, atol=1e-12) if tight else {}),
            )
            return y

        native = solve(TRUTH, tight=True)
        off = TRUTH * np.exp(np.linspace(-0.1, 0.1, 36))
        native_off = solve(off, tight=True)
        off_gap = np.max(np.abs(native_off - reference(off, times, literal)))
        nominal_gap = np.max(np.abs(native - observed))
        nominal_cost = weighted_error(observed[:, 1:], native[:, 1:], scales)
        probe_cost = weighted_error(observed[:, 1:], native_off[:, 1:], scales)
        check(
            report,
            "sixteen-conditions-independent-trajectory",
            nominal_gap < 2e-6 and off_gap < 2e-6 and probe_cost > nominal_cost,
            nominal_state_gap=nominal_gap,
            off_nominal_state_gap=off_gap,
            nominal_cost=nominal_cost,
            perturbed_cost=probe_cost,
            conditions=len(CONDITIONS),
        )
        check(
            report,
            "weighted-loss-fixture",
            weighted_error(
                np.array([1.0, 2.0]), np.array([2.0, 4.0]), np.array([2.0, 4.0])
            )
            == 0.5,
        )

        def score(logp):
            prediction = solve(10**logp)
            return weighted_error(observed[:, 1:], prediction[:, 1:], scales)

        regions = [("broad", BROAD, budget)]
        if not args.broad_only:
            regions.append(("near-table-scale", FOCUSED, min(budget, 2000)))
        for label, bounds, calls in regions:
            for algorithm in ("de", "cmaes"):
                for retry in range(retries):
                    fit = search(
                        args,
                        f"{label}-{algorithm}-{retry}",
                        algorithm,
                        score,
                        bounds,
                        51941 + retry + args.seed_offset,
                        calls,
                        pop=24,
                        checkpoints=checkpoints,
                        context=context,
                    )
                    row = dict(
                        region=label,
                        optimizer=algorithm,
                        retry=retry,
                        cost=fit.value,
                        calls=len(fit.scores),
                        failures=fit.failures,
                        parameters=10**fit.best,
                        checkpoints={
                            k: float(np.min(fit.scores[:k, 0]))
                            for k in (*checkpoints, calls)
                            if k <= calls
                        },
                    )
                    rejected = []
                    seen = set()
                    ordering = np.argsort(fit.scores[:, 0], kind="stable")
                    # Try the best distinct final points, then check the
                    # earlier checkpoint incumbents if those solves fail.
                    probes = [int(j) for j in ordering if fit.scores[j, 0] < PENALTY]
                    provision = [
                        int(np.argmin(fit.scores[:k, 0]))
                        for k in reversed(checkpoints)
                        if k < calls
                    ]
                    candidates = []
                    for j in probes:
                        token = fit.candidates[j].tobytes()
                        if token not in seen:
                            candidates.append(j)
                            seen.add(token)
                        if len(candidates) >= 24:
                            break
                    for j in dict.fromkeys(candidates + provision):
                        if fit.scores[j, 0] >= PENALTY:
                            continue
                        logp = fit.candidates[j]
                        try:
                            prediction = solve(10**logp, tight=True)
                        except ValueError as error:
                            rejected.append(dict(call=int(j + 1), reason=str(error)))
                            continue
                        row["parameters"] = 10**logp
                        row["rescored_call"] = int(j + 1)
                        row["verified"] = weighted_error(
                            observed[:, 1:], prediction[:, 1:], scales
                        )
                        row["state_rmse"] = np.sqrt(
                            ((prediction - observed) ** 2).mean(axis=(0, 1))
                        )
                        row["log_parameter_error"] = logp - np.log10(TRUTH)
                        break
                    row["verification_failures"] = rejected
                    results.append(row)
        eligible = [row for row in results if "verified" in row]
        require(bool(eligible), "optimizer finished without a feasible trajectory")
        best = min(eligible, key=lambda row: row["verified"])
        output = solve(best["parameters"], tight=True)
        dense = np.arange(0.0, 121.0, 3.0)
        dense_obs = reference(TRUTH, dense, literal)
        dense_fit = solve(best["parameters"], dense, tight=True)
        dense_score = weighted_error(
            dense_obs[:, 1:],
            dense_fit[:, 1:],
            np.max(np.abs(dense_obs[:, 1:]), axis=1, keepdims=True),
        )
    report.update(
        doi=DOI,
        results=results,
        best=best,
        published_objective=PUBLISHED,
        dense_grid_score=dense_score,
        time_grid="6 to 120 s by 6 s",
        original_raw_data="unavailable",
        condition_solves_attempted=16 * sum(r["calls"] for r in results),
        printed_equation_12=literal,
    )

    def draw(fig):
        for ax, i in zip(fig.subplots(1, 2), [7, 3], strict=True):
            for j in range(16):
                color = f"C{j % 10}"
                ax.plot(times, observed[j, :, i], c=color, lw=1)
                ax.plot(times, output[j, :, i], c=color, lw=1, ls="--")
            ax.set(
                xlabel="Time (s)",
                ylabel=STATES[i],
                title="Reference and fitted trajectories",
            )

    save_plot(args, "pathway-fits.png", draw)
