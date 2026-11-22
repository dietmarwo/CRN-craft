"""Optimize single-cell sampling schedules with probability sensitivities."""

from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
from itertools import product
import json
import re
from pathlib import Path

import crn_craft as crn
from crn_craft.sensitivities import fisher_information, forward_sensitivities
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import poisson

from common import require, write
from .search import PACKAGES as PACKAGES
from .search import check, save_plot, search
from .sources import MANIFEST, data_option, fixtures, read_tar

DOI = "10.1371/journal.pcbi.1006365"
A = (2.0, 3.0, 9.0)
B = (3.0, 4.0, 6.0)
AXES = (np.arange(1.0, 5.0), np.arange(1.0, 7.0, 0.5), np.arange(1.0, 15.0))
TOGGLE_BOUNDS = [(g[0], g[-1]) for g in AXES]
SWEEP = np.geomspace(0.02, 700.0, 60)
FILE_PATTERN = re.compile(
    r"/all_fims_sampled_fim/exp_design_FIM_alpha_(.+)_delta_(.+)_uv_(.+)_sampled.txt$"
)
FILES = {
    "bursting_gene/out/E_opt_fsamp_fsp1.0_100.txt": "burst",
    "toggle/parameter_samples.txt": "samples",
    "toggle/ref_pars.txt": "ref",
    "toggle/out/initial_condition.txt": "initial",
}
DIRECTIONS = ("bx", "by", "kx", "ky", "ayx", "axy", "dy")
TOGGLE_NAMES = ("bx", "by", "kx", "ky", "ayx", "axy", "nyx", "nxy", "dx", "dy")


def eigenvalue(probabilities, derivatives, cells, floor=1e-8):
    fim, clipped = fisher_information(
        probabilities, derivatives, cells, probability_floor=floor
    )
    return float(np.linalg.eigvalsh(fim)[0]), float(np.max(clipped))


def _table(data):
    """Parse a whitespace-separated numeric text member; nothing else."""
    values = np.loadtxt(io.StringIO(data.decode("ascii")), dtype=float, ndmin=1)
    require(np.isfinite(values).all(), "nonfinite value in the design archive")
    return values


def archive(path, spec):
    """Read the archived Fisher matrices and reference tables.

    Returns the smallest eigenvalue per design and sample, the 100 matrices
    of each design (archive parameter order, natural-log coordinates) and
    the reference tables. The content is verified while it is read.
    """
    designs, matrices, references = {}, {}, {}

    def handle(name, data):
        match = FILE_PATTERN.search("/" + name)
        if match:
            values = _table(data)
            require(values.size == 4900, f"unexpected matrix size in {name}")
            stack = values.reshape(7, 7, 100).transpose(2, 0, 1)
            key = tuple(map(float, match.groups()))
            matrices[key] = stack
            designs[key] = np.linalg.eigvalsh(stack)[:, 0]
        elif name in FILES:
            references[FILES[name]] = _table(data)

    read_tar(path, spec, handle)
    require(
        len(designs) == 672 and set(references) == set(FILES.values()),
        "incomplete design archive",
    )
    return designs, matrices, references


def analytic_check(report):
    net = crn.Model("sensitivity_birth_death")
    net.add_species("rna", 0)
    net.add_parameter("birth", 2.0)
    net.add_parameter("loss", 0.3)
    net.add_reaction(crn.mass_action("birth"), [], ["rna"])
    net.add_reaction(crn.mass_action("loss"), ["rna"], [])
    g = crn.prepare_fsp(net, {"rna": 45})
    segments = [(1.0, {"birth": 2.0, "loss": 0.3}), (3.0, {"birth": 3.0, "loss": 0.3})]
    times = (1.0, 3.0)
    result = forward_sensitivities(g, segments, times, ("birth",), tolerance=2e-8)
    bins, p, tangent = result.observed(("rna",))
    a, b = -np.expm1(-0.3) / 0.3, -np.expm1(-0.6) / 0.3
    mu = np.array([2 * a, 2 * a * np.exp(-0.6) + 3 * b])[:, None]
    exact = poisson.pmf(bins[:, 0], mu)
    derivative = (
        exact * (bins[:, 0] / mu - 1.0) * np.array([a, a * np.exp(-0.6) + b])[:, None]
    )
    law_gap = np.max(np.abs(exact - p))
    analytic_gap = np.max(np.abs(tangent[:, 0] - derivative))

    def finite_difference(h):
        history = []
        for sign in (1, -1):
            p0 = None
            series = []
            last = 0.0
            for end, binding in segments:
                solved = g.solve(
                    end - last,
                    [end - last],
                    initial=p0,
                    parameters=binding | {"birth": binding["birth"] + sign * h},
                    tolerance=1e-11,
                )
                series.append(solved.marginal("rna"))
                snap = solved.snapshots[-1]
                require(snap.physical_leakage < 3e-13, "test law truncated")
                law = np.maximum(snap.probabilities, 0)
                p0 = dict(zip(map(tuple, g.states), law / law.sum(), strict=True))
                last = end
            history.append(np.asarray(series))
        return np.max(np.abs(tangent[:, 0] - (history[0] - history[1]) / (2 * h)))

    gaps = [finite_difference(h) for h in (1e-3, 5e-4)]
    fim, _ = fisher_information(p, tangent, 1000.0, probability_floor=1e-12)
    check(
        report,
        "probability-sensitivities-and-Fisher-fixture",
        max(law_gap, analytic_gap, *gaps) < 2e-7 and fim[0, 0] > 0,
        law_gap=law_gap,
        analytic_gap=analytic_gap,
        full_law_difference=gaps,
        mass_error=np.max(np.abs(result.probabilities.sum(axis=1) - 1.0)),
        derivative_mass_error=np.max(np.abs(result.derivatives.sum(axis=2))),
        fim=fim,
    )


def burst_sweep(report, archive_curve):
    net = crn.Model("burst_sampling")
    net.add_species({"off": 1, "on": 0, "rna": 0})
    binding = dict(kon=0.05, koff=0.15, synth=5.0, removal=0.05)
    for name, value in binding.items():
        net.add_parameter(name, value)
    for name, lhs, rhs in (
        ("kon", ["off"], ["on"]),
        ("koff", ["on"], ["off"]),
        ("synth", ["on"], ["on", "rna"]),
        ("removal", ["rna"], []),
    ):
        net.add_reaction(crn.mass_action(name), lhs, rhs, name=name)
    graph = crn.prepare_fsp(net, {"on": 1, "off": 1, "rna": 215})
    script_grid = 25 * np.logspace(-3, 1.5, 50)
    schedules = np.r_[SWEEP, script_grid][:, None] * np.arange(1, 5)[None, :]
    times = np.unique(schedules)
    res = forward_sensitivities(
        graph,
        [(times[-1], binding)],
        times,
        ("koff", "synth"),
        tolerance=1e-8,
        max_work=300_000,
    )
    _, p, s = res.observed(("rna",))
    table = []
    for schedule in schedules:
        idx = np.searchsorted(times, schedule)
        table.append(eigenvalue(p[idx], s[idx], 1000.0)[0])
    values = np.asarray(table)
    native_peak = float(SWEEP[np.argmax(values[: len(SWEEP)])])
    archive_peak = float(script_grid[np.argmax(archive_curve)])
    script_peak = float(script_grid[np.argmax(values[len(SWEEP) :])])
    near = np.abs(np.log(SWEEP / 6.1)).argmin()
    check(
        report,
        "bursting-interval-near-6.1-minutes",
        SWEEP[near - 1] <= native_peak <= SWEEP[near + 1]
        and abs(archive_peak - 6.1) < 0.1,
        native_peak=native_peak,
        archive_peak=archive_peak,
        native_on_archive_grid=script_peak,
        leakage=res.probabilities[-1, -1],
        archive_grid_points=len(script_grid),
    )
    return net, values[: len(SWEEP)], script_grid


def toggle_network(reference):
    net = crn.Model("toggle_design")
    net.add_species({"x": 0, "y": 0})
    net.add_observable("cx", "x")
    net.add_observable("cy", "y")
    for name, value in zip(TOGGLE_NAMES, reference, strict=True):
        net.add_parameter(name, float(value))
    net.add_reaction("bx+kx/(1+ayx*cy^nyx)", [], ["x"], name="birth_x")
    net.add_reaction("by+ky/(1+axy*cx^nxy)", [], ["y"], name="birth_y")
    net.add_reaction(crn.mass_action("dx"), ["x"], [], name="death_x")
    net.add_reaction(crn.mass_action("dy"), ["y"], [], name="death_y")
    return net


# UV raises the decay of the archive's x (lambda cI, the paper's y):
# basal rate, amplitude and shape of rate = basal + a UV^2 / (shape + UV^3).
# The archived matrices follow the constants in the archive code; the
# article's Eq. 28 prints other constants (checked in the lesson).
PULSE_LAWS = {"archive-code": (3.84e-4, 0.02, 2500.0), "paper": (3.8e-4, 0.002, 1250.0)}
# The archived matrices use three snapshots at dt, 2 dt and 3 dt (time
# zero carries no information: the initial law is fixed) and 100 cells.
OBSERVATIONS = 3
CELLS = 100.0


class ToggleDesign:
    def __init__(self, reference, bounds=(70, 100, 2500), *, pulse_law="archive-code"):
        self.indices = np.array([0, 1, 2, 3, 4, 5, 9])
        self.reference = reference
        self.net = toggle_network(reference["ref"])
        self.bounds = bounds
        self.pulse_law = pulse_law
        mx, my, shape = bounds
        retained = [
            (x, y)
            for x, y in product(range(mx + 1), range(my + 1))
            if (x - 4) * (y - 4) < shape
        ]
        original = [
            (x, y)
            for x, y in product(range(71), range(101))
            if (x - 4) * (y - 4) < 2500
        ]
        self.initial = dict(zip(original, reference["initial"], strict=True))
        self.graph = crn.prepare_fsp(
            self.net, {"x": mx, "y": my}, initial=self.initial, retained=retained
        )
        initial_hash = hashlib.sha256(reference["initial"].astype("<f8").tobytes())
        self.version = (
            f"toggle-{OBSERVATIONS}-snapshots-{CELLS:g}-cells-{pulse_law}-"
            + initial_hash.hexdigest()
        )

    def solve(self, design, sample, *, tolerance=3e-6):
        """Probabilities and log-parameter sensitivities at the snapshots."""
        duration, spacing, magnitude = design
        par = self.reference["ref"].copy()
        par[self.indices] = self.reference["samples"][sample]
        binding = dict(zip(TOGGLE_NAMES, map(float, par), strict=True))
        basal, amplitude, shape = PULSE_LAWS[self.pulse_law]
        binding["dx"] = basal
        pulse = binding | {
            "dx": float(basal + amplitude * magnitude**2 / (shape + magnitude**3))
        }
        times = tuple(3600 * spacing * i for i in range(1, OBSERVATIONS + 1))
        end = times[-1]
        knots = [(3600.0, binding), (3600 * (1 + duration), pulse), (end, binding)]
        segments = [(stop, value) for stop, value in knots if stop < end]
        active = binding if end <= 3600.0 or end > 3600 * (1 + duration) else pulse
        segments.append((end, active))
        return forward_sensitivities(
            self.graph,
            segments,
            times,
            DIRECTIONS,
            logarithmic=True,
            tolerance=tolerance,
            max_work=150_000,
        )

    def matrix(self, design, sample, *, tolerance=3e-6):
        """The 7 x 7 Fisher matrix in the archive's parameter order."""
        r = self.solve(design, sample, tolerance=tolerance)
        # Both counts are measured, so every retained state is a bin.
        fim, _ = fisher_information(
            r.probabilities[:, :-1],
            r.derivatives[:, :, :-1],
            CELLS,
            probability_floor=1e-8,
        )
        return fim, r

    def one(self, design, sample, *, tolerance=3e-6):
        r = self.solve(design, sample, tolerance=tolerance)
        value, dropped = eigenvalue(
            r.probabilities[:, :-1], r.derivatives[:, :, :-1], CELLS
        )
        return value, max(r.probabilities[:, -1]), dropped, r.work

    def scores(self, design, samples, cache=None, tolerance=3e-6):
        rows = []
        for sample in samples:
            inputs = json.dumps(
                [
                    self.version,
                    self.bounds,
                    list(map(float, design)),
                    int(sample),
                    tolerance,
                    self.reference["samples"][sample].tolist(),
                    self.reference["ref"].tolist(),
                ],
                allow_nan=False,
            )
            key = hashlib.sha256(inputs.encode()).hexdigest()
            path = None if cache is None else cache / key[:2] / f"{key}.json"
            if path is not None and path.is_file():
                rows.append(json.loads(path.read_text())["value"])
                continue
            value = self.one(design, sample, tolerance=tolerance)
            if path is not None:
                from tempfile import NamedTemporaryFile

                path.parent.mkdir(parents=True, exist_ok=True)
                with NamedTemporaryFile(
                    mode="w", dir=path.parent, delete=False
                ) as handle:
                    json.dump(dict(input=json.loads(inputs), value=value), handle)
                (path.parent / handle.name).replace(path)
            rows.append(value)
        return np.asarray(rows)


def options(parser):
    data_option(parser)
    parser.add_argument(
        "--pulse-law", choices=tuple(PULSE_LAWS), default="archive-code"
    )
    parser.add_argument(
        "--native-samples",
        type=int,
        help="Number of parameter sets per native candidate.",
    )
    parser.add_argument(
        "--native-evaluations", type=int, help="Native candidate evaluations per start."
    )
    parser.add_argument("--native-retries", type=int)
    parser.add_argument(
        "--finalist-samples",
        type=int,
        help="Parameter samples for the independent final comparison.",
    )
    parser.add_argument(
        "--expanded-bounds",
        type=int,
        nargs=3,
        default=(80, 110, 3300),
        metavar=("X", "Y", "SHAPE"),
    )

    parser.add_argument("--cache", type=Path)
    parser.add_argument(
        "--native-pair-samples",
        type=int,
        default=0,
        help="Independently compare A and B on this many native parameter samples.",
    )
    parser.add_argument(
        "--native-grid",
        action="store_true",
        help="Replay the full native grid with the same cache.",
    )


def run(args, report):
    """Check the 6.1-minute design and optimize the robust toggle experiment."""
    analytic_check(report)
    file = fixtures(args.data_dir, "C8")["fsp-fim.tar.gz"]
    spec = json.loads(MANIFEST.read_text())["C8"]["fsp-fim.tar.gz"]
    designs, matrices, references = archive(file, spec)
    burst, values, script_grid = burst_sweep(report, references["burst"])
    write(args.output / "bursting-model.json", burst.to_dict())
    archive_winner = max(designs, key=lambda k: designs[k].mean())
    ranking = int(np.sum(designs[A] > designs[B]))
    check(
        report,
        "archived-grid-and-97-of-100-ranking",
        len(designs) == 672 and archive_winner == A and ranking == 97,
        best=archive_winner,
        comparison=ranking,
        a_score=np.mean(designs[A]),
        b_score=np.mean(designs[B]),
    )
    grid = np.array([designs[k].mean() for k in product(*AXES)]).reshape(
        tuple(map(len, AXES))
    )
    interpolate = RegularGridInterpolator(AXES, grid)
    fit_archive = search(
        args,
        "archive-interpolation-de",
        "de",
        lambda x: -float(interpolate(x)[0]),
        TOGGLE_BOUNDS,
        72611 + args.seed_offset,
        160,
    )
    toggle = ToggleDesign(references, pulse_law=args.pulse_law)
    write(args.output / "toggle-model.json", toggle.net.to_dict())
    fast = args.preset == "smoke"
    n = args.native_samples or (1 if fast else 100)
    require(1 <= n <= 100, "native samples must be in 1..100")
    samples = np.arange(n)

    def at(guess):
        return toggle.scores(guess, samples, args.cache)

    def cost(guess):
        return -float(np.mean(at(guess)[:, 0]))

    tested = {label: at(d) for label, d in (("A", A), ("B", B))}
    searches = []
    candidates = [A, B]
    for retry in range(args.native_retries or (1 if fast else 4)):
        fit = search(
            args,
            f"native-de-{retry}",
            "de",
            cost,
            TOGGLE_BOUNDS,
            72580 + retry + args.seed_offset,
            args.native_evaluations or (4 if fast else 100),
        )
        searches.append(
            dict(
                candidate=fit.best,
                score=-fit.value,
                calls=len(fit.scores),
                failures=fit.failures,
            )
        )
        if fit.failures < len(fit.scores):
            candidates.append(tuple(fit.best))
    final = max(candidates, key=lambda d: -cost(d))
    bounds = tuple(args.expanded_bounds)
    require(
        all(a <= b for a, b in zip(toggle.bounds, bounds, strict=True)),
        "verification bounds must include the original retained set",
    )
    expanded = ToggleDesign(references, bounds, pulse_law=args.pulse_law)
    final_count = args.finalist_samples or n
    require(1 <= final_count <= 100, "finalist sample count must be in 1..100")

    def verify(i):
        return (
            expanded.scores(final, [i], args.cache, tolerance=1e-6)[0],
            expanded.scores(A, [i], args.cache, tolerance=1e-6)[0],
        )

    with ThreadPoolExecutor(args.threads) as pool:
        rows = np.asarray(list(pool.map(verify, range(final_count))))
    verified, verified_a = rows[:, 0, :], rows[:, 1, :]
    if args.native_pair_samples:
        require(
            1 <= args.native_pair_samples <= 100,
            "native pair samples must be in 1..100",
        )
        ids = range(args.native_pair_samples)

        def pair(i):
            return (
                i,
                toggle.scores(A, [i], args.cache)[0],
                toggle.scores(B, [i], args.cache)[0],
            )

        with ThreadPoolExecutor(args.threads) as pool:
            comparison = list(pool.map(pair, ids))
        write(args.output / "native-pair.json", comparison)
        candidates_a = np.array([a for _, a, _ in comparison])
        candidates_b = np.array([b for _, _, b in comparison])
        report["native_pair"] = dict(
            samples=len(comparison),
            wins=int(sum(candidates_a[:, 0] > candidates_b[:, 0])),
            archive_wins=int(sum(designs[A][: len(ids)] > designs[B][: len(ids)])),
            average_a=np.mean(candidates_a[:, 0]),
            average_b=np.mean(candidates_b[:, 0]),
            max_leakage=np.max([candidates_a[:, 1], candidates_b[:, 1]]),
        )
    if args.native_grid:
        require(args.cache is not None, "native grid needs --cache for restart")
        write(
            args.output / "native-grid.json",
            [dict(design=d, score=at(d)) for d in product(*AXES)],
        )
    rows = [dict(design=d, score=raw.mean()) for d, raw in designs.items()]
    write(args.output / "archive-grid.json", rows)
    report.update(
        doi=DOI,
        pulse_law=args.pulse_law,
        native_samples=n,
        native_designs=tested,
        searches=searches,
        native_winner=final,
        native_tight_final=verified,
        native_finalist_sample_count=final_count,
        expanded_bounds=bounds,
        native_tight_A=verified_a,
        archive_de=dict(candidate=fit_archive.best, interpolated=-fit_archive.value),
        experiment_A=A,
        experiment_B=B,
        native_pair_wins=int(np.sum(tested["A"][:, 0] > tested["B"][:, 0])),
    )

    def draw(fig):
        burst_axis, pair_axis = fig.subplots(1, 2)
        burst_axis.semilogx(SWEEP, values / np.max(values), label="CRN-craft")
        burst_axis.semilogx(
            script_grid,
            references["burst"] / np.max(references["burst"]),
            ls="--",
            label="Archived",
        )
        burst_axis.set(
            xlabel="Sampling interval (min)", ylabel="Normalized information"
        )
        burst_axis.legend()
        pair_axis.loglog(designs[B], designs[A], ".")
        extent = [
            min(designs[A].min(), designs[B].min()),
            max(designs[A].max(), designs[B].max()),
        ]
        pair_axis.plot(extent, extent, color="grey")
        pair_axis.set(xlabel="Design B information", ylabel="Design A information")

    save_plot(args, "information-designs.png", draw)
