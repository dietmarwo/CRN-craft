"""Compare MODE with random designs on the Summit SnAr virtual reactor."""

import csv
from itertools import product
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import mannwhitneyu

from .concentrations import Solver, document
from .search import PACKAGES as PACKAGES
from .search import check, search
from .sources import data_option, fixtures

DOI = "10.1002/cmtd.202000051"
# tau (min), pyrrolidine equivalents, inlet 2,4-dinitrofluorobenzene (M), T (C)
BOUNDS = ((0.5, 2.0), (1.0, 5.0), (0.1, 0.5), (30.0, 120.0))
HV_POINT = np.array([-2957.0, 10.7])
MOLAR = np.array([159.09, 71.12, 210.21, 210.21, 261.33])
# The rate constant of Q formation. Summit releases 0.5.0 to 0.8.3, which
# produced the archived benchmark, use the ortho constant ka there; the pinned
# source revision uses kb. `archive` reproduces the archived outputs.
KINETICS = {"archive": 0, "pinned": 1}
MODEL = document(
    ("A", "B", "P", "Q", "R"),
    ("a0", "b0", "ka", "kb", "kc", "kd", "kq"),
    ("a0", "b0", "0", "0", "0"),
    ("-v1-v2", "-v1-v2-v3-v4", "v1-v3", "vq-v4", "v3+v4"),
    assignments=(
        ("v1", "ka*A*B"),
        ("v2", "kb*A*B"),
        ("vq", "kq*A*B"),
        ("v3", "kc*P*B"),
        ("v4", "kd*Q*B"),
    ),
)
BUDGET = 50
STARTS = 20
LONG_BUDGET = 500


def inputs(design, kinetics="archive"):
    """Initial concentrations (M) and rate constants (1/(M min)) of a design."""
    _, equivalents, a0, temp = design
    delta = 1 / (temp + 273.71) - 1 / (90 + 273.71)
    rates = (
        0.6
        * np.array([57.9, 2.70, 0.865, 1.63])
        * np.exp(-np.array([33.3, 35.3, 38.9, 44.8]) * delta / 0.008314)
    )
    return (a0, equivalents * a0, *rates, rates[KINETICS[kinetics]])


def outlets_reference(design, *, loose=False, kinetics="archive"):
    """Outlet concentrations by SciPy: tight, or Summit's RK45 defaults."""
    a0, b0, ka, kb, kc, kd, kq = inputs(design, kinetics)

    def rhs(_, state):
        a, b, p, q, _ = state
        if loose and a < 1e-6 * a0:
            a = 0.0
        if loose and b < 1e-6 * b0:
            b = 0.0
        consumption = (ka + kb) * a * b
        return [
            -consumption,
            -consumption - kc * b * p - kd * b * q,
            ka * a * b - kc * b * p,
            kq * a * b - kd * b * q,
            kc * b * p + kd * b * q,
        ]

    options = {} if loose else dict(rtol=2e-11, atol=1e-13)
    solution = solve_ivp(rhs, [0.0, design[0]], [a0, b0, 0.0, 0.0, 0.0], **options)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[:, -1]


def metrics(residence, outlets):
    """(-STY in kg/(m3 h), E-factor) from the outlet concentrations in M."""
    outlet = np.maximum(outlets, 0.0)
    product_mass = MOLAR[2] * outlet[2]
    other_mass = MOLAR @ outlet - product_mass
    sty = max(60 * product_mass / residence, 1e-6)
    waste = (
        1000.0
        if np.isclose(outlet[2], 0.0)
        else min(1000.0, (789.0 + other_mass) / product_mass)
    )
    return np.array([-sty, waste])


def hypervolume(points, point=HV_POINT):
    """Area of the union of minimizing rectangles anchored at point."""
    ordered = points[np.argsort(points[:, 0], kind="stable")]
    previous = point[1]
    area = 0.0
    for x, y in ordered:
        if x < point[0] and y < previous:
            area += (point[0] - x) * (previous - y)
            previous = y
    return float(area)


def nondominated(points):
    """Indices of the points that no other point weakly dominates, by first objective."""
    order = np.lexsort((points[:, 1], points[:, 0]))
    kept, best = [], np.inf
    for i in order:
        if points[i, 1] < best:
            kept.append(int(i))
            best = points[i, 1]
    return kept


def options(parser):
    data_option(parser)
    parser.add_argument(
        "--kinetics",
        choices=sorted(KINETICS),
        default="archive",
        help="Q-formation constant: archive (ka, as in the archived runs) or pinned (kb).",
    )


def table(path, header, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def archive_groups(archived):
    groups = {}
    for row in archived:
        key = (
            row["strategy_name"],
            row["transform_name"],
            row["sty_tolerance"],
            row["e_factor_tolerance"],
            row["expression"],
        )
        groups.setdefault(key, []).append(float(row["terminal_hypervolume"]))
    rows = []
    for key, values in groups.items():
        q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
        rows.append([*key, len(values), median, q25, q75])
    return sorted(rows, key=lambda r: -r[-3])


def run(args, report):
    """Check the emulator against the archive, then run 20 starts of 50 evaluations."""
    sources = fixtures(args.data_dir, "C7")
    archived = read(sources["snar_results.csv"])
    noise = sorted(set(row["noise_level"] for row in archived))
    check(
        report,
        "archived-reference-noise",
        noise == ["0"],
        rows=len(archived),
        noise=noise,
    )
    sample = read(sources["pareto_front_snar.csv"])[1:]  # row 1 holds Summit's types
    kinetics = args.kinetics
    smoke = args.preset == "smoke"
    starts = 2 if smoke else STARTS
    out = args.output

    # 1. The emulator against the archived Summit sample, at Summit's settings.
    picks = np.random.default_rng(4983).choice(
        len(sample), 8 if smoke else 200, replace=False
    )
    archive_gap = []
    for i in picks:
        row = sample[i]
        design = [
            float(row[k]) for k in ("tau", "equiv_pldn", "conc_dfnb", "temperature")
        ]
        ours = metrics(
            design[0], outlets_reference(design, loose=True, kinetics=kinetics)
        )
        stored = np.array([-float(row["sty"]), float(row["e_factor"])])
        archive_gap.append(np.abs(ours - stored) / np.abs(stored))
    archive_gap = np.max(archive_gap, axis=0)
    check(
        report,
        "summit-archive-outputs",
        kinetics != "archive" or bool(np.all(archive_gap < 1e-6)),
        kinetics=kinetics,
        designs=len(picks),
        max_relative_sty=archive_gap[0],
        max_relative_e_factor=archive_gap[1],
    )

    lows, highs = np.asarray(BOUNDS).T
    designs = np.array(
        [*product(*BOUNDS), *np.random.default_rng(4984).uniform(lows, highs, (8, 4))]
    )
    with Solver(args.ode_binary, out / "model.json", MODEL) as engine:

        def outlets(design, *, tight=False):
            state, _ = engine.batch(
                [inputs(design, kinetics)],
                [0.0, design[0]],
                **(dict(rtol=2e-10, atol=1e-12) if tight else {}),
            )
            return state[0, -1]

        # 2. Native solve against a tight SciPy solve and Summit's RK45 setting.
        parity = []
        for design in designs:
            native = outlets(design)
            fine = outlets_reference(design, kinetics=kinetics)
            summit = outlets_reference(design, loose=True, kinetics=kinetics)
            parity.append(
                [
                    *design,
                    np.max(np.abs(native - fine)),
                    *np.abs(metrics(design[0], native) - metrics(design[0], fine)),
                    *np.abs(metrics(design[0], summit) - metrics(design[0], fine)),
                ]
            )
        parity = np.array(parity)
        state_err = float(np.max(parity[:, 4]))
        score_err = float(np.max(parity[:, 5:7]))
        check(
            report,
            "reactor-outlet-parity",
            state_err < 1e-6 and score_err < 0.02,
            designs=len(designs),
            state_max_abs=state_err,
            objective_max_abs=score_err,
            summit_rk45_objective_max_abs=float(np.max(parity[:, 7:9])),
        )
        check(
            report,
            "hypervolume-hand-fixture",
            hypervolume(
                np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 4.0], [6.0, 0.0]]),
                np.array([5.0, 5.0]),
            )
            == 11.0,
        )

        # 3. The searches: every start has exactly its budget of evaluations.
        arms = [("mode", "mode", BUDGET), ("random", "random", BUDGET)]
        if not smoke:
            arms.append(("mode-500", "mode", LONG_BUDGET))
        terminal, curves, fronts = [], [], []
        for label, algorithm, budget in arms:
            for seed in range(starts):
                fit = search(
                    args,
                    f"{label}-{seed}",
                    algorithm,
                    lambda x: metrics(x[0], outlets(x)),
                    BOUNDS,
                    9971 + seed + args.seed_offset,
                    budget,
                    pop=10,
                    observables=2,
                )
                hv = [hypervolume(fit.scores[:n]) for n in range(1, budget + 1)]
                curves.append((label, hv))
                terminal.append(
                    [label, seed, budget, len(fit.scores), fit.failures, hv[-1]]
                )
                for j in nondominated(fit.scores):
                    fronts.append(
                        [
                            label,
                            seed,
                            j + 1,
                            *fit.candidates[j],
                            -fit.scores[j, 0],
                            fit.scores[j, 1],
                        ]
                    )
        check(
            report,
            "exact-evaluations",
            all(r[3] == r[2] and r[4] == 0 for r in terminal),
            runs=len(terminal),
        )
        rescore = max(
            float(
                np.max(
                    np.abs(
                        metrics(d[0], outlets(d, tight=True))
                        - metrics(d[0], outlets(d))
                    )
                )
            )
            for d in designs
        )
        report["tight_objective_difference"] = rescore

    # 4. Saved tables; the figures are drawn from them alone.
    table(
        out / "terminal.csv",
        ["arm", "seed", "budget", "evaluations", "failures", "hypervolume"],
        terminal,
    )
    rows = []
    for label in dict.fromkeys(name for name, _ in curves):
        values = np.array([hv for name, hv in curves if name == label])
        q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75], axis=0)
        rows.extend(
            [label, n + 1, q25[n], median[n], q75[n]] for n in range(values.shape[1])
        )
    table(out / "trajectory.csv", ["arm", "evaluations", "q25", "median", "q75"], rows)
    table(
        out / "fronts.csv",
        [
            "arm",
            "seed",
            "call",
            "tau",
            "equiv_pldn",
            "conc_dfnb",
            "temperature",
            "sty",
            "e_factor",
        ],
        fronts,
    )
    table(
        out / "archive-groups.csv",
        [
            "strategy",
            "transform",
            "sty_tolerance",
            "e_factor_tolerance",
            "expression",
            "runs",
            "median",
            "q25",
            "q75",
        ],
        archive_groups(archived),
    )
    points = np.array([[-float(r["sty"]), float(r["e_factor"])] for r in sample])
    table(
        out / "archive-front.csv",
        ["sty", "e_factor"],
        [[-points[i, 0], points[i, 1]] for i in nondominated(points)],
    )
    table(
        out / "parity.csv",
        [
            "tau",
            "equiv_pldn",
            "conc_dfnb",
            "temperature",
            "native_state_abs",
            "native_sty_abs",
            "native_e_factor_abs",
            "summit_rk45_sty_abs",
            "summit_rk45_e_factor_abs",
        ],
        parity,
    )

    by_arm = {
        label: np.array([r[5] for r in terminal if r[0] == label])
        for label in dict.fromkeys(r[0] for r in terminal)
    }
    summary = {
        label: dict(
            runs=len(v),
            median=float(np.median(v)),
            q25_q75=np.quantile(v, [0.25, 0.75]),
        )
        for label, v in by_arm.items()
    }
    test = mannwhitneyu(by_arm["mode"], by_arm["random"], alternative="two-sided")
    rng = np.random.default_rng(7)
    boot = [
        np.median(rng.choice(by_arm["mode"], len(by_arm["mode"])))
        - np.median(rng.choice(by_arm["random"], len(by_arm["random"])))
        for _ in range(2000)
    ]
    report.update(
        doi=DOI,
        kinetics=kinetics,
        reference_point=HV_POINT,
        terminal=summary,
        mode_vs_random=dict(
            mann_whitney_u=float(test.statistic),
            p_two_sided=float(test.pvalue),
            median_difference=float(
                np.median(by_arm["mode"]) - np.median(by_arm["random"])
            ),
            bootstrap_95=np.quantile(boot, [0.025, 0.975]),
        ),
    )
    draw(out, out)


ARCHIVE_STRATEGIES = (
    "TSEMO",
    "SNOBFIT",
    "SOBO",
    "GRYFFIN",
    "Random",
    "DRO",
    "NelderMead",
)


def draw(tables: Path, folder: Path):
    """C7-pareto.png and C7-hypervolume.png from the saved tables."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    tables, folder = Path(tables), Path(folder)
    fronts = read(tables / "fronts.csv")
    archive = np.array(
        [
            [float(r["sty"]), float(r["e_factor"])]
            for r in read(tables / "archive-front.csv")
        ]
    )

    figure = Figure(figsize=(6.4, 4.4), layout="constrained")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    axis.plot(
        archive[:, 0],
        archive[:, 1],
        c="0.55",
        lw=1.5,
        label="Summit archive: nondominated designs of its NSGA-II sample",
    )
    for label, marker, color in (
        ("random", "x", "C1"),
        ("mode", "o", "C0"),
        ("mode-500", ".", "C2"),
    ):
        rows = [r for r in fronts if r["arm"] == label]
        if rows:
            axis.scatter(
                [float(r["sty"]) for r in rows],
                [float(r["e_factor"]) for r in rows],
                s=12,
                marker=marker,
                c=color,
                alpha=0.7,
                label=f"ours, {label}: nondominated designs of each start",
            )
    axis.axvline(-HV_POINT[0], c="k", ls=":", lw=0.8)
    axis.axhline(
        HV_POINT[1], c="k", ls=":", lw=0.8, label="hypervolume reference (2957, 10.7)"
    )
    axis.set(
        xlabel="Space–time yield (kg m⁻³ h⁻¹)",
        ylabel="E-factor",
        ylim=(8, 20),
        xlim=(0, 13000),
    )
    axis.legend(fontsize=7, loc="upper right")
    figure.savefig(folder / "C7-pareto.png", dpi=155)

    figure = Figure(figsize=(9.0, 4.2), layout="constrained")
    FigureCanvasAgg(figure)
    axis, longer = figure.subplots(1, 2, width_ratios=(3, 2))
    trajectory = read(tables / "trajectory.csv")

    def band(ax, label, color, text):
        rows = [r for r in trajectory if r["arm"] == label]
        n = [int(r["evaluations"]) for r in rows]
        ax.plot(n, [float(r["median"]) for r in rows], c=color, label=text)
        ax.fill_between(
            n,
            [float(r["q25"]) for r in rows],
            [float(r["q75"]) for r in rows],
            color=color,
            alpha=0.15,
        )

    band(axis, "mode", "C0", "ours, MODE: median and quartiles of 20 starts")
    band(axis, "random", "C1", "ours, random: median and quartiles of 20 starts")
    groups = read(tables / "archive-groups.csv")
    shown = []
    for strategy in ARCHIVE_STRATEGIES:
        best = max(
            (g for g in groups if g["strategy"] == strategy),
            key=lambda g: float(g["median"]),
            default=None,
        )
        if best is not None:
            shown.append(best)
    for k, g in enumerate(shown):
        x = 51 + 0.9 * k
        axis.errorbar(
            x,
            float(g["median"]),
            yerr=[
                [float(g["median"]) - float(g["q25"])],
                [float(g["q75"]) - float(g["median"])],
            ],
            fmt="s",
            ms=4,
            c=f"C{k + 3}",
            capsize=2,
            label=f"archive {g['strategy']} ({g['transform']}), n={g['runs']}",
        )
    axis.set(
        xlabel="Evaluations",
        ylabel="Hypervolume, reference (−2957, 10.7)",
        xlim=(0, 58),
    )
    axis.legend(fontsize=6.5, loc="upper left")
    axis.set_title(
        "At 50 evaluations; squares: archived runs (median, quartiles)", fontsize=8
    )
    if any(r["arm"] == "mode-500" for r in trajectory):
        band(longer, "mode-500", "C2", "ours, MODE, 500 evaluations")
        tsemo = next((g for g in shown if g["strategy"] == "TSEMO"), None)
        if tsemo is not None:
            longer.axhline(
                float(tsemo["median"]),
                c="C3",
                ls="--",
                lw=1,
                label="archive TSEMO median at 50",
            )
        longer.axvline(50, c="k", ls=":", lw=0.8)
        longer.set(xlabel="Evaluations", xlim=(0, 500))
        longer.legend(fontsize=6.5, loc="lower right")
        longer.set_title("Larger budget (not the paper's protocol)", fontsize=8)
    figure.savefig(folder / "C7-hypervolume.png", dpi=155)
