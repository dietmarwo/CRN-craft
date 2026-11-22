"""Reproduce the STAT5 fit using the archived PEtab measurements."""

import csv
from functools import reduce
import operator
from pathlib import Path

from crn_craft._xml import _read_xml
import libsbml as sbml
import numpy as np
from scipy.integrate import solve_ivp

from common import require, write
from . import petab
from .concentrations import Solver, document
from .search import PACKAGES, PENALTY, check, options_resume, search
from .sources import data_option, fixtures

PACKAGES = (*PACKAGES, "PyYAML", "python-libsbml")
DOI = "10.1021/pr5006923"
REFERENCE_NLL = 138.22199970618203
STATES = ("STAT5A", "STAT5B", "pApA", "pApB", "pBpB", "nucpApA", "nucpApB", "nucpBpB")
KINETIC = (
    "Epo_degradation_BaF3",
    "k_exp_hetero",
    "k_exp_homo",
    "k_imp_hetero",
    "k_imp_homo",
    "k_phos",
    "ratio",
    "specC17",
)


def model(problem):
    assignments = (
        ("light", "1.25e-7 * exp(-Epo_degradation_BaF3 * time) * k_phos"),
        ("aa", "light * STAT5A^2"),
        ("ab", "light * STAT5A * STAT5B"),
        ("bb", "light * STAT5B^2"),
        ("iaq", "k_imp_homo * pApA"),
        ("ihq", "k_imp_hetero * pApB"),
        ("ibq", "k_imp_homo * pBpB"),
        ("eaq", "(0.45/1.4) * k_exp_homo * nucpApA"),
        ("ehq", "(0.45/1.4) * k_exp_hetero * nucpApB"),
        ("ebq", "(0.45/1.4) * k_exp_homo * nucpBpB"),
    )
    rhs = (
        "-2*aa-ab+2*eaq+ehq",
        "-2*bb-ab+2*ebq+ehq",
        "aa-iaq",
        "ab-ihq",
        "bb-ibq",
        "(1.4/0.45)*(iaq-eaq)",
        "(1.4/0.45)*(ihq-ehq)",
        "(1.4/0.45)*(ibq-ebq)",
    )
    initial = ("207.6*ratio", "207.6*(1-ratio)", *["0"] * 6)
    return document(
        STATES,
        KINETIC,
        initial,
        rhs,
        assignments=assignments,
        observations=[r["observableFormula"] for r in problem.observations.values()],
    )


OPERATIONS = {
    sbml.AST_PLUS: operator.add,
    sbml.AST_MINUS: operator.sub,
    sbml.AST_TIMES: operator.mul,
    sbml.AST_DIVIDE: operator.truediv,
    sbml.AST_POWER: operator.pow,
    sbml.AST_FUNCTION_POWER: operator.pow,
}


def arithmetic(node, env):
    """Evaluate a libSBML formula tree of names, numbers, + - * / ^ and exp."""
    kind = node.getType()
    if kind == sbml.AST_NAME_TIME:
        return env["time"]
    if node.isName():
        name = node.getName()
        require(name in env, f"SBML formula uses an unknown name {name!r}")
        return env[name]
    if node.isInteger():
        return node.getInteger()
    if node.isReal():
        return node.getReal()
    require(
        kind in OPERATIONS or kind == sbml.AST_FUNCTION_EXP,
        f"unsupported SBML math {node.getName() or kind}",
    )
    children = [arithmetic(node.getChild(i), env) for i in range(node.getNumChildren())]
    if kind == sbml.AST_FUNCTION_EXP:
        return np.exp(children[0])
    if kind == sbml.AST_MINUS and len(children) == 1:
        return -children[0]
    return reduce(OPERATIONS[kind], children)


def sbml_model(problem):
    """Read the SBML file through the package's XML guard, then with libSBML."""
    raw = problem.model.read_bytes()
    _read_xml(raw, petab.LIMIT)  # refuses DTDs, entities, non-UTF-8, deep nesting
    source = sbml.readSBMLFromString(raw.decode("utf-8"))
    sbmodel = source.getModel()
    require(
        source.getLevel() == 2
        and sbmodel is not None
        and source.getNumErrors(sbml.LIBSBML_SEV_ERROR) == 0
        and source.getNumErrors(sbml.LIBSBML_SEV_FATAL) == 0,
        "expected a valid SBML Level 2 model",
    )
    require(
        sbmodel.getNumEvents() == 0
        and sbmodel.getNumFunctionDefinitions() == 0
        and sbmodel.getNumConstraints() == 0,
        "events, function definitions and constraints are not supported",
    )
    require(
        all(r.isAssignment() for r in sbmodel.getListOfRules()),
        "only SBML assignment rules are supported",
    )
    require(
        all(
            not s.getHasOnlySubstanceUnits()
            and not s.getBoundaryCondition()
            and not s.getConstant()
            for s in sbmodel.getListOfSpecies()
        ),
        "only variable concentration species are supported",
    )
    require(
        all(c.getConstant() for c in sbmodel.getListOfCompartments()),
        "compartment sizes must be constant",
    )
    require(
        all(
            r.getKineticLaw() is not None
            and r.getKineticLaw().getNumParameters() == 0
            and all(
                ref.getStoichiometryMath() is None
                for side in (r.getListOfReactants(), r.getListOfProducts())
                for ref in side
            )
            for r in sbmodel.getListOfReactions()
        ),
        "reactions need a kinetic law without local parameters"
        " and constant stoichiometry",
    )
    require(
        sorted(s.getId() for s in sbmodel.getListOfSpecies()) == sorted(STATES),
        "the SBML species differ from the hand-mapped STAT5 states",
    )
    return sbmodel


def reference(problem, parameters, *, rtol=2e-10, atol=1e-12):
    """Solve the SBML reactions in amount space from libSBML formula trees.

    This route does not use the hand-written equations of `model`: it builds
    the stoichiometry and the rates from the SBML file, integrates molecule
    amounts (concentration times compartment size) with SciPy's Radau method
    and divides by the sizes afterwards.
    """
    sbmodel = sbml_model(problem)
    compartments = {c.getId(): c.getSize() for c in sbmodel.getListOfCompartments()}
    species = {s.getId(): s for s in sbmodel.getListOfSpecies()}
    volumes = np.array([compartments[species[s].getCompartment()] for s in STATES])
    env = compartments | parameters | {"time": 0.0}
    init = {
        a.getSymbol(): arithmetic(a.getMath(), env)
        for a in sbmodel.getListOfInitialAssignments()
    }
    y0 = [init.get(s, species[s].getInitialConcentration()) for s in STATES] * volumes
    reactions = list(sbmodel.getListOfReactions())
    rules = [(r.getVariable(), r.getMath()) for r in sbmodel.getListOfRules()]
    stoich = np.zeros((len(STATES), len(reactions)))
    for j, reaction in enumerate(reactions):
        for sign, side in [
            (-1, reaction.getListOfReactants()),
            (1, reaction.getListOfProducts()),
        ]:
            for species_ref in side:
                stoich[STATES.index(species_ref.getSpecies()), j] += (
                    sign * species_ref.getStoichiometry()
                )
    kinetics = [reaction.getKineticLaw().getMath() for reaction in reactions]

    def rhs(t, amount):
        state = env | dict(zip(STATES, amount / volumes, strict=True)) | {"time": t}
        for name, expression in rules:
            state[name] = arithmetic(expression, state)
        return stoich @ [arithmetic(k, state) for k in kinetics]

    result = solve_ivp(
        rhs,
        [0.0, problem.times[-1]],
        y0,
        method="Radau",
        t_eval=problem.times,
        rtol=rtol,
        atol=atol,
    )
    require(result.success, result.message)
    states = result.y.T / volumes
    signals = [
        sbml.parseL3Formula(r["observableFormula"])
        for r in problem.observations.values()
    ]
    require(all(s is not None for s in signals), "unparsable observable formula")
    observations = [
        [arithmetic(s, env | dict(zip(STATES, y, strict=True))) for s in signals]
        for y in states
    ]
    return states, np.asarray(observations)


# Log10 shifts of the six kinetic rates for the two off-nominal checks; the
# second also moves the two fixed parameters so that their mapping is tested.
OFF_NOMINAL = {
    "off-nominal-a": (
        dict(
            Epo_degradation_BaF3=0.3,
            k_exp_hetero=1,
            k_exp_homo=-0.5,
            k_imp_hetero=0.5,
            k_imp_homo=-3,
            k_phos=-0.5,
        ),
        {},
    ),
    "off-nominal-b": (
        dict(
            Epo_degradation_BaF3=-0.5,
            k_exp_hetero=3,
            k_exp_homo=1,
            k_imp_hetero=-1,
            k_imp_homo=-5,
            k_phos=0.5,
        ),
        dict(ratio=0.5, specC17=0.2),
    ),
}


def vectors(problem):
    """The nominal vector and two declared off-nominal vectors inside the bounds."""
    result = {"nominal": problem.nominal}
    for name, (shifts, fixed) in OFF_NOMINAL.items():
        p = problem.nominal | fixed
        for key, shift in shifts.items():
            p[key] *= 10**shift
        for key, value in p.items():
            meta = problem.parameters[key]
            require(
                float(meta["lowerBound"]) <= value <= float(meta["upperBound"]),
                f"{name} leaves the bounds of {key}",
            )
        result[name] = p
    return result


def options(parser):
    data_option(parser)
    parser.add_argument("--retries", type=int)
    parser.add_argument(
        "--profile-evaluations",
        type=int,
        default=0,
        help="Optional fcmaes calls with Gaussian noise scales solved conditionally.",
    )
    parser.add_argument("--joint-refine-evaluations", type=int, default=0)
    options_resume(parser)
    parser.add_argument(
        "--evaluations", type=int, help="Total per retry across DE and CMA-ES."
    )


def per_signal(problem, costs):
    return [float(costs[problem.o_index == i].sum()) for i in range(3)]


def write_csv(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(
            [f"{v:.12g}" if isinstance(v, float) else v for v in row] for row in rows
        )


def run(args, report):
    """Verify the archived likelihood and fit all nine PEtab coordinates."""
    files = fixtures(args.data_dir, "C9")
    problem = petab.read(files["problem.yaml"])
    costs = problem.score_archive(files["simulations.tsv"])
    nominal = problem.nominal
    archived_nll = float(costs.sum())
    check(
        report,
        "archived-gaussian-likelihood",
        len(costs) == 48
        and len(problem.fits) == 9
        and abs(archived_nll - REFERENCE_NLL) < 1e-6,
        nll=archived_nll,
        per_signal=per_signal(problem, costs),
    )
    fast = args.preset == "smoke"
    retries = args.retries or (1 if fast else 8)
    budget = args.evaluations or (50 if fast else 5000)
    tables = args.output / "tables"
    tables.mkdir()
    with Solver(args.ode_binary, args.output / "model.json", model(problem)) as solver:

        def simulate(parameters, *, tight=False, times=problem.times):
            y, outputs = solver.batch(
                [[parameters[key] for key in KINETIC]],
                times,
                **(dict(rtol=2e-10, atol=1e-12) if tight else {}),
            )
            return y[0], outputs[0]

        def score(x):
            p = problem.bind(x)
            _, ys = simulate(p)
            return float(problem.score(ys, p).sum())

        # The hand-mapped native model against the libSBML route at the nominal
        # vector and two off-nominal vectors; the search tolerance against the
        # tight one. Targets are the plan's: 0.01 signal units and 0.05 nats.
        verification = []
        for name, parameters in vectors(problem).items():
            ref_state, ref_signal = reference(problem, parameters)
            state, signal = simulate(parameters, tight=True)
            _, loose = simulate(parameters)
            row = dict(
                vector=name,
                signal_max_abs=float(np.max(np.abs(signal - ref_signal))),
                state_max_abs=float(np.max(np.abs(state - ref_state))),
                nll_native=float(problem.score(signal, parameters).sum()),
                nll_reference=float(problem.score(ref_signal, parameters).sum()),
                nll_search_tolerance=float(problem.score(loose, parameters).sum()),
            )
            verification.append(row)
            if name == "nominal":
                native, readouts, nominal_reference = state, signal, ref_signal
        write_csv(
            tables / "verification.csv",
            list(verification[0]),
            [list(row.values()) for row in verification],
        )
        nll = verification[0]["nll_native"]
        scale = np.array([1.4] * 5 + [0.45] * 3)
        incidence = np.array([[1, 0, 2, 1, 0, 2, 1, 0], [0, 1, 0, 1, 2, 0, 1, 2]])
        mass = (native * scale) @ incidence.T
        mass_error = np.max(np.abs(mass - mass[0]))
        check(
            report,
            "native-nominal-signal-and-nll",
            all(
                row["signal_max_abs"] < 0.01
                and abs(row["nll_native"] - row["nll_reference"]) < 0.05
                and abs(row["nll_search_tolerance"] - row["nll_native"]) < 0.05
                for row in verification
            )
            and abs(nll - archived_nll) < 0.05
            and mass_error < 1e-4,
            nll=nll,
            signal_max_abs=verification[0]["signal_max_abs"],
            state_max_abs=verification[0]["state_max_abs"],
            archive_nll_gap=nll - archived_nll,
            dimer_amount_drift=mass_error,
            vectors=verification,
        )
        results = []
        for i in range(retries):
            seed = 47425 + i + args.seed_offset
            de = search(
                args,
                f"repeat-{i}-de",
                "de",
                score,
                problem.bounds,
                seed,
                budget // 2,
                pop=12,
            )
            cma = search(
                args,
                f"repeat-{i}-cma",
                "cmaes",
                score,
                problem.bounds,
                seed + 1000,
                budget - budget // 2,
                pop=12,
                x0=de.best,
            )
            best = min((de, cma), key=lambda fit: fit.value)
            require(best.value < PENALTY, "no finite STAT5 fit")
            params = problem.bind(best.best)
            _, fitted_signal = simulate(params, tight=True)
            tight_costs = problem.score(fitted_signal, params)
            results.append(
                (
                    dict(
                        retry=i,
                        seed=seed,
                        parameters=params,
                        proposed=best.value,
                        verified=float(tight_costs.sum()),
                        per_signal=per_signal(problem, tight_costs),
                        de_best=de.value,
                        calls=len(de.scores) + len(cma.scores),
                        failures=de.failures + cma.failures,
                    ),
                    fitted_signal,
                )
            )
        best, signal = min(results, key=lambda pair: pair[0]["verified"])
        fits = {"best plain retry": (best["parameters"], best["calls"])}
        if args.profile_evaluations:
            # Kinetic rates only; each noise scale is replaced by its conditional
            # maximum-likelihood value, the root mean square residual of its
            # signal (clipped to the declared bounds). `nominal` supplies only the
            # two fixed parameters; its noise scales are overwritten.
            kinetics = [i for i, name in enumerate(problem.fits) if name in KINETIC]
            labels = [problem.fits[i] for i in kinetics]

            def profiled(coords, *, tight=False):
                p = nominal | dict(zip(labels, 10**coords, strict=True))
                _, ys = simulate(p, tight=tight)
                residual = problem.y - ys[problem.t_index, problem.o_index]
                for label in sorted(set(problem.noise)):
                    rows = np.array([x == label for x in problem.noise])
                    meta = problem.parameters[label]
                    p[label] = float(
                        np.clip(
                            np.sqrt(np.mean(residual[rows] ** 2)),
                            float(meta["lowerBound"]),
                            float(meta["upperBound"]),
                        )
                    )
                return float(problem.score(ys, p).sum()), p, ys

            fit = search(
                args,
                "profiled-noise-de",
                "de",
                lambda x: profiled(x)[0],
                problem.bounds[kinetics],
                7458 + args.seed_offset,
                args.profile_evaluations // 2,
                pop=24,
            )
            second = search(
                args,
                "profiled-noise-cma",
                "cmaes",
                lambda x: profiled(x)[0],
                problem.bounds[kinetics],
                7459 + args.seed_offset,
                args.profile_evaluations - args.profile_evaluations // 2,
                pop=16,
                x0=fit.best,
                sigma=0.06,
            )
            candidate = min((fit, second), key=lambda arm: arm.value)
            value, pp, vv = profiled(candidate.best, tight=True)
            calls = len(fit.scores) + len(second.scores)
            report["profiled_noise"] = dict(
                verified=value,
                per_signal=per_signal(problem, problem.score(vv, pp)),
                parameters=pp,
                calls=calls,
            )
            fits["profiled noise scales"] = (pp, calls)
            write(
                args.output / "profiled-signals.json",
                dict(times=problem.times, signals=vv),
            )
            if args.joint_refine_evaluations:
                # CMA-ES on all nine coordinates, started at the profiled optimum.
                joint = search(
                    args,
                    "joint-refinement-cma",
                    "cmaes",
                    score,
                    problem.bounds,
                    7460 + args.seed_offset,
                    args.joint_refine_evaluations,
                    x0=problem.scaled(pp),
                    sigma=0.015,
                    pop=14,
                )
                joint_parameters = problem.bind(joint.best)
                _, joint_signals = simulate(joint_parameters, tight=True)
                joint_costs = problem.score(joint_signals, joint_parameters)
                report["joint_refinement"] = dict(
                    verified=float(joint_costs.sum()),
                    per_signal=per_signal(problem, joint_costs),
                    parameters=joint_parameters,
                    calls=len(joint.scores),
                    start="profiled-noise optimum",
                )
                fits["profiled, then joint refinement"] = (
                    joint_parameters,
                    calls + len(joint.scores),
                )
                write(
                    args.output / "joint-refined-signals.json",
                    dict(times=problem.times, signals=joint_signals),
                )
        # Saved tables: the inputs of `draw`.
        archived = {}
        for row in petab.table(files["simulations.tsv"], "simulations"):
            archived[petab.key(row)] = float(row["simulation"])
        dense = np.linspace(0.0, float(problem.times[-1]), 161)
        curves = {"nominal (native)": simulate(nominal, tight=True, times=dense)[1]}
        measured = {"nominal (native)": readouts}
        nll_rows = [
            ["archive", "archived simulations", *per_signal(problem, costs), 0],
            [
                "nominal (native)",
                "PEtab nominal vector, native solve",
                *per_signal(problem, problem.score(readouts, nominal)),
                0,
            ],
            [
                "nominal (libSBML route)",
                "PEtab nominal vector, SciPy solve",
                *per_signal(problem, problem.score(nominal_reference, nominal)),
                0,
            ],
        ]
        for label, (parameters, calls) in fits.items():
            _, at_data = simulate(parameters, tight=True)
            measured[label] = at_data
            curves[label] = simulate(parameters, tight=True, times=dense)[1]
            nll_rows.append(
                [
                    label,
                    "fcmaes search",
                    *per_signal(problem, problem.score(at_data, parameters)),
                    calls,
                ]
            )
    ids = list(problem.observations)
    write_csv(
        tables / "nll.csv",
        ["source", "route", *ids, "total", "calls"],
        [[*r[:5], sum(r[2:5]), r[5]] for r in nll_rows],
    )
    write_csv(
        tables / "measurements.csv",
        ["observable", "time", "measurement", "archive", *measured],
        [
            [
                row["observableId"],
                float(row["time"]),
                float(row["measurement"]),
                archived[petab.key(row)],
                *[
                    float(v[problem.t_index[i], problem.o_index[i]])
                    for v in measured.values()
                ],
            ]
            for i, row in enumerate(problem.measurements)
        ],
    )
    write_csv(
        tables / "curves.csv",
        ["time", *[f"{k}|{o}" for k in curves for o in ids]],
        [
            [float(t), *[float(v[j, i]) for v in curves.values() for i in range(3)]]
            for j, t in enumerate(dense)
        ],
    )
    fitted = [*[(f"retry {r['retry']}", r) for r, _ in results]]
    if "profiled_noise" in report:
        fitted.append(("profiled", report["profiled_noise"]))
    if "joint_refinement" in report:
        fitted.append(("joint", report["joint_refinement"]))
    write_csv(
        tables / "fits.csv",
        ["fit", "nll", "calls", *[f"log10 {k}" for k in problem.fits]],
        [["nominal", nll, 0, *problem.scaled(nominal)]]
        + [
            [name, r["verified"], r["calls"], *problem.scaled(r["parameters"])]
            for name, r in fitted
        ],
    )
    report.update(
        doi=DOI,
        nominal_nll=nll,
        archive_nll=archived_nll,
        best=best,
        fcs=[row for row, _ in results],
        nominal_params=nominal,
    )
    write(
        args.output / "trajectories.json",
        dict(
            times=problem.times,
            nominal=readouts,
            archive_model=nominal_reference,
            fitted=signal,
        ),
    )
    draw(tables, args.output)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def draw(tables: Path, folder: Path):
    """Draw the C9 figures from the saved tables only."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    tables, folder = Path(tables), Path(folder)
    measurements = read_csv(tables / "measurements.csv")
    curves = read_csv(tables / "curves.csv")
    ids = list(dict.fromkeys(r["observable"] for r in measurements))
    labels = list(dict.fromkeys(k.split("|")[0] for k in curves[0] if "|" in k))
    figure = Figure(figsize=(10, 3.6), layout="constrained")
    styles = {
        "nominal (native)": dict(color="0.35", lw=2.2),
        "best plain retry": dict(color="tab:red", lw=1.2, ls="--"),
        "profiled noise scales": dict(color="tab:orange", lw=1.0, ls=":"),
        "profiled, then joint refinement": dict(color="tab:blue", lw=1.0),
    }
    for i, (ax, name) in enumerate(zip(figure.subplots(1, 3), ids, strict=True)):
        rows = [r for r in measurements if r["observable"] == name]
        t = [float(r["time"]) for r in rows]
        ax.scatter(
            t,
            [float(r["measurement"]) for r in rows],
            s=14,
            color="k",
            label="measurement (archive)",
            zorder=3,
        )
        ax.scatter(
            t,
            [float(r["archive"]) for r in rows],
            s=30,
            marker="x",
            color="tab:green",
            label="archived simulation",
            zorder=4,
        )
        for label in labels:
            ax.plot(
                [float(r["time"]) for r in curves],
                [float(r[f"{label}|{name}"]) for r in curves],
                label=f"ours: {label}",
                **styles.get(label, {}),
            )
        ax.set(xlabel="time (min)", ylabel="percent" if i == 0 else None, title=name)
    figure.axes[0].legend(fontsize=6.5, loc="lower right")
    figure.savefig(folder / "C9-signals.png", dpi=155)

    fits = read_csv(tables / "fits.csv")
    names = [k for k in fits[0] if k.startswith("log10 ")]
    figure = Figure(figsize=(10, 4), layout="constrained")
    left, right = figure.subplots(1, 2, width_ratios=[2.2, 1])
    colors = {"nominal": "k", "profiled": "tab:orange", "joint": "tab:blue"}
    for row in fits:
        retry = row["fit"].startswith("retry")
        left.plot(
            range(len(names)),
            [float(row[k]) for k in names],
            marker="o",
            ms=3,
            lw=0.8,
            color="tab:red" if retry else colors[row["fit"]],
            alpha=0.45 if retry else 1.0,
            label=None
            if retry and row["fit"] != "retry 0"
            else ("plain retries (8)" if retry else row["fit"]),
        )
    left.set_xticks(
        range(len(names)), [k[6:] for k in names], rotation=40, ha="right", fontsize=7
    )
    left.set(
        ylabel="log10 value (bounds -5 to 5)",
        ylim=(-5.3, 5.3),
        title="fitted coordinates",
    )
    left.legend(fontsize=7)
    order = sorted(fits, key=lambda r: float(r["nll"]))
    right.barh(
        [r["fit"] for r in order],
        [float(r["nll"]) - float(order[0]["nll"]) + 0.01 for r in order],
        color=[colors.get(r["fit"], "tab:red") for r in order],
    )
    right.set_xscale("log")
    right.invert_yaxis()
    right.set(
        xlabel=f"NLL above {float(order[0]['nll']):.4f} (+0.01)", title="verified NLL"
    )
    right.tick_params(labelsize=7)
    figure.savefig(folder / "C9-fits.png", dpi=155)
