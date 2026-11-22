"""Reproduce the STAT5 fit using the archived PEtab measurements."""

from functools import reduce
import operator

from crn_craft._xml import _read_xml
import libsbml as sbml
import numpy as np
from scipy.integrate import solve_ivp

from common import require, write
from . import petab
from .concentrations import Solver, document
from .search import PACKAGES, PENALTY, check, options_resume, save_plot, search
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
    kind = node.getType()
    if kind == sbml.AST_NAME_TIME:
        return env["time"]
    if node.isName():
        return env[node.getName()]
    if node.isInteger():
        return node.getInteger()
    if node.isReal():
        return node.getReal()
    children = [arithmetic(node.getChild(i), env) for i in range(node.getNumChildren())]
    if kind == sbml.AST_FUNCTION_EXP:
        return np.exp(children[0])
    if kind == sbml.AST_MINUS and len(children) == 1:
        return -children[0]
    return reduce(OPERATIONS[kind], children)


def reference(problem, parameters):
    # Independent amount-space interpretation of the pinned SBML source.
    raw = problem.model.read_bytes()
    _read_xml(raw, 4 << 20)
    source = sbml.readSBMLFromString(raw.decode("utf-8"))
    sbmodel = source.getModel()
    require(
        source.getLevel() == 2 and sbmodel is not None,
        "expected the archived SBML model",
    )
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
        rtol=2e-10,
        atol=1e-12,
    )
    require(result.success, result.message)
    states = result.y.T / volumes
    signals = [
        sbml.parseL3Formula(r["observableFormula"])
        for r in problem.observations.values()
    ]
    observations = [
        [arithmetic(s, env | dict(zip(STATES, y, strict=True))) for s in signals]
        for y in states
    ]
    return states, np.asarray(observations)


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
        per_signal=[float(costs[problem.o_index == i].sum()) for i in range(3)],
    )
    ref_state, ref_signal = reference(problem, nominal)
    fast = args.preset == "smoke"
    retries = args.retries or (1 if fast else 8)
    budget = args.evaluations or (50 if fast else 5000)
    with Solver(args.ode_binary, args.output / "model.json", model(problem)) as solver:

        def simulate(parameters, *, tight=False):
            y, outputs = solver.batch(
                [[parameters[key] for key in KINETIC]],
                problem.times,
                **(dict(rtol=2e-10, atol=1e-12) if tight else {}),
            )
            return y[0], outputs[0]

        def score(x):
            p = problem.bind(x)
            _, ys = simulate(p)
            return float(problem.score(ys, p).sum())

        native, readouts = simulate(nominal, tight=True)
        discrepancy = float(np.max(np.abs(readouts - ref_signal)))
        nll = float(problem.score(readouts, nominal).sum())
        scale = np.array([1.4] * 5 + [0.45] * 3)
        incidence = np.array([[1, 0, 2, 1, 0, 2, 1, 0], [0, 1, 0, 1, 2, 0, 1, 2]])
        mass = (native * scale) @ incidence.T
        mass_error = np.max(np.abs(mass - mass[0]))
        check(
            report,
            "native-nominal-signal-and-nll",
            discrepancy < 0.01 and abs(nll - archived_nll) < 0.05 and mass_error < 1e-4,
            nll=nll,
            signal_max_abs=discrepancy,
            state_max_abs=np.max(np.abs(native - ref_state)),
            archive_nll_gap=nll - archived_nll,
            dimer_amount_drift=mass_error,
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
            tight_nll = float(problem.score(fitted_signal, params).sum())
            results.append(
                (
                    dict(
                        retry=i,
                        seed=seed,
                        parameters=params,
                        proposed=best.value,
                        verified=tight_nll,
                        calls=len(de.scores) + len(cma.scores),
                        failures=de.failures + cma.failures,
                    ),
                    fitted_signal,
                )
            )
        best, signal = min(results, key=lambda pair: pair[0]["verified"])
        if args.profile_evaluations:
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
            report["profiled_noise"] = dict(
                verified=value,
                parameters=pp,
                calls=len(fit.scores) + len(second.scores),
            )
            write(
                args.output / "profiled-signals.json",
                dict(times=problem.times, signals=vv),
            )
            if args.joint_refine_evaluations:
                seed = [
                    np.log10(pp[k])
                    if problem.parameters[k]["parameterScale"] == "log10"
                    else pp[k]
                    for k in problem.fits
                ]
                joint = search(
                    args,
                    "joint-refinement-cma",
                    "cmaes",
                    score,
                    problem.bounds,
                    7460 + args.seed_offset,
                    args.joint_refine_evaluations,
                    x0=seed,
                    sigma=0.015,
                    pop=14,
                )
                joint_parameters = problem.bind(joint.best)
                _, joint_signals = simulate(joint_parameters, tight=True)
                joint_nll = float(problem.score(joint_signals, joint_parameters).sum())
                report["joint_refinement"] = dict(
                    verified=joint_nll,
                    parameters=joint_parameters,
                    calls=len(joint.scores),
                )
                write(
                    args.output / "joint-refined-signals.json",
                    dict(times=problem.times, signals=joint_signals),
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
            archive_model=ref_signal,
            fitted=signal,
        ),
    )

    def draw(fig):
        for i, ax in enumerate(fig.subplots(1, 3)):
            choose = problem.o_index == i
            ax.scatter(
                problem.times[problem.t_index[choose]],
                problem.y[choose],
                s=12,
                label="Measurements",
            )
            ax.plot(problem.times, readouts[:, i], label="Nominal")
            ax.plot(problem.times, signal[:, i], label="fcmaes")
            ax.set(xlabel="Time (min)", title=list(problem.observations)[i])
        fig.axes[0].legend(fontsize=7)

    save_plot(args, "signal-fits.png", draw)
