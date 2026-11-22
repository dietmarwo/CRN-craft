"""Fit the seven two-state stimulus hypotheses to the archived RNA/site counts."""

from itertools import combinations

import crn_craft as crn
from crn_craft.likelihood import HistogramGroup, histogram_score, measured_counts
from crn_craft.stationary import _solve
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply, spsolve

from common import require, write
from .search import PACKAGES as PACKAGES
from .search import (
    PENALTY,
    check,
    objective_context,
    options_resume,
    save_plot,
    search,
)
from .sources import data_option, fixtures

DOI = "10.1186/s13059-017-1297-9"
TIMES = (0, 5, 15, 25)
TOTALS = (137, 135, 38, 66)
FIXED_DECAY = 0.0559
NAMES = ("kon", "koff", "basal", "active")
TABLE_BASAL = (0.0093, 0.1108, 0.0078, 8.14)
TABLE_STIMULATED = (0.0839, 0.0189, 0.0078, 8.14)
HYPOTHESES = tuple(h for n in (1, 2, 3) for h in combinations((0, 1, 3), n))
BOUNDS = ((1e-5, 3.0), (1e-5, 3.0), (1e-6, 3.0), (0.01, 25.0))


MAX_BYTES = 1 << 20
MAX_RNA = 10_000


def parse_counts(raw, name):
    """One cell per nonempty row: two site intensities and an RNA count.

    A site is active when its intensity is positive, the two-state rule of
    the BayFish data reader. UTF-16 needs a byte-order mark; other input
    must be UTF-8 without NUL bytes.
    """
    require(len(raw) <= MAX_BYTES, f"{name}: more than {MAX_BYTES} bytes")
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        text = raw.decode("utf-16")
    else:
        require(b"\0" not in raw, f"{name}: NUL byte without a UTF-16 mark")
        text = raw.decode("utf-8-sig")
    pairs = []
    for number, line in enumerate(text.splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        where = f"{name}:{number}"
        require(len(fields) == 3, f"{where}: expected 3 fields, found {len(fields)}")
        try:
            a, b = float(fields[0]), float(fields[1])
            copies = int(fields[2])
        except ValueError:
            raise ValueError(f"{where}: non-numeric field") from None
        require(
            np.isfinite(a) and np.isfinite(b) and a >= 0 and b >= 0,
            f"{where}: site intensity must be finite and nonnegative",
        )
        require(0 <= copies <= MAX_RNA, f"{where}: RNA count out of range")
        pairs.append((int(a > 0) + int(b > 0), copies))
    require(pairs, f"{name}: no cells")
    return np.array(pairs, dtype=int)


def read_counts(files):
    return [
        parse_counts(files[name].read_bytes(), name)
        for name in (f"myData_Npas4_t{t}_List.txt" for t in TIMES)
    ]


def model(limit, *, reflect=False):
    net = crn.Model("joint_promoter_rna")
    net.add_species({"off": 2, "on": 0, "rna": 0})
    for name, value in zip(NAMES, TABLE_BASAL, strict=True):
        net.add_parameter(name, value)
    net.add_parameter("decay", FIXED_DECAY)
    for name in ("off", "on", "rna"):
        net.add_observable(f"{name}_count", name)
    net.add_reaction(crn.mass_action("kon"), ["off"], ["on"], name="switch_on")
    net.add_reaction(crn.mass_action("koff"), ["on"], ["off"], name="switch_off")
    net.add_reaction(crn.mass_action("decay"), ["rna"], [], name="rna_loss")
    for name, rate in (("off", "basal"), ("on", "active")):
        propensity = f"{rate} * {name}_count"
        if reflect:
            propensity += f" * (rna_count < {limit})"
        net.add_reaction(propensity, [name], [name, "rna"], name=f"transcribe_{name}")
    return net


def cme(rates, limit, *, reflect):
    """Independent generator in lexicographic (sites, RNA) coordinates."""
    k_on, k_off, basal, active = rates
    size = 3 * (limit + 1)
    rows, cols, vals = [], [], []
    for sites in range(3):
        for m in range(limit + 1):
            j = sites * (limit + 1) + m
            changes = []
            if sites < 2:
                changes.append((j + limit + 1, (2 - sites) * k_on))
            if sites:
                changes.append((j - limit - 1, sites * k_off))
            if m:
                changes.append((j - 1, FIXED_DECAY * m))
            rate = (2 - sites) * basal + sites * active
            if m < limit:
                changes.append((j + 1, rate))
            elif not reflect:
                changes.append((size, rate))
            changes.append((j, -sum(x[1] for x in changes)))
            for i, value in changes:
                rows.append(i)
                cols.append(j)
                vals.append(value)
    return coo_matrix((vals, (rows, cols)), shape=(size + 1, size + 1)).tocsc()


def normalized_law(rates, limit):
    q = cme(rates, limit, reflect=True)[:-1, :-1].tolil()
    q[0, :] = 1.0
    b = np.zeros(q.shape[0])
    b[0] = 1
    p = spsolve(q.tocsc(), b)
    require(
        np.min(p) >= -1e-10 and abs(np.sum(p) - 1) < 1e-10,
        "normalized stationary solve returned an invalid law",
    )
    p = np.maximum(p, 0.0)
    return p / np.sum(p)


def reference(before, after, limit):
    q = cme(after, limit, reflect=False)
    p = np.r_[normalized_law(before, limit), 0.0]
    series = [p[:-1].reshape(3, limit + 1)]
    for dt in np.diff(TIMES):
        p = expm_multiply(q * dt, p)
        series.append(p[:-1].reshape(3, limit + 1))
    return np.array(series)


class StimulusFit:
    def __init__(self, groups, limit):
        self.limit = limit
        # A count above the cutoff would wrap into the next site row.
        require(
            all(int(group[:, 1].max()) <= limit for group in groups),
            f"an observed RNA count exceeds the cutoff {limit}",
        )
        bounds = {"off": 2, "on": 2, "rna": limit}
        self.basal = crn.prepare_fsp(model(limit, reflect=True), bounds)
        self.post = crn.prepare_fsp(model(limit), bounds)
        self.obs = tuple(
            np.bincount(
                group[:, 0] * (limit + 1) + group[:, 1], minlength=3 * (limit + 1)
            )
            for group in groups
        )
        self.initial_states = tuple(tuple(row) for row in self.basal.states)
        i, j = measured_counts(self.post.states, self.post.species, ("on", "rna")).T
        self.readout = i * (limit + 1) + j
        self.probability_shape = (len(TIMES), 3, limit + 1)

    def probabilities(self, before, after):
        parameters = dict(zip(NAMES, before, strict=True))
        kon, koff, basal, active = before
        f = kon / (kon + koff)
        sites = min(2, int(3 * f))
        mean = 2 * ((1 - f) * basal + f * active) / FIXED_DECAY
        by_name = dict(off=2 - sites, on=sites, rna=min(self.limit, int(mean)))
        anchor = tuple(by_name[key] for key in self.basal.species)
        try:
            law, residual, _, _ = _solve(
                self.basal, parameters, None, anchor, 256 << 20
            )
        except RuntimeError as error:
            if "invalid probabilities" not in str(error):
                raise
            # The normalized balance equations avoid a rare state's unit anchor.
            p = normalized_law(before, self.limit)
            layout = self.basal.species
            linear = self.basal.states[:, layout.index("on")] * (self.limit + 1)
            linear = linear + self.basal.states[:, layout.index("rna")]
            law = p[linear]
            residual = np.linalg.norm(
                cme(before, self.limit, reflect=True) @ np.r_[p, 0.0], ord=1
            )
            require(residual < 1e-8, "normalized stationary solve has a large residual")
        start = dict(zip(self.initial_states, law, strict=True))
        result = self.post.solve(
            TIMES[-1],
            TIMES,
            initial=start,
            parameters=dict(zip(NAMES, after, strict=True)),
            tolerance=2e-9,
            max_work=100_000,
        )
        pmf = np.array(
            [
                np.bincount(
                    self.readout, snap.probabilities, minlength=3 * (self.limit + 1)
                )
                for snap in result.snapshots
            ]
        ).reshape(self.probability_shape)
        boundary_flux = float(
            pmf[0, :, -1]
            @ np.array([2 * before[2], before[2] + before[3], 2 * before[3]])
        )
        return pmf, boundary_flux, result.snapshots[-1].physical_leakage, residual

    def score(self, pmf):
        # Restore the original probability floor for fixed-vector parity.
        groups = [
            HistogramGroup(
                "stimulus", "archive", float(t), count, np.maximum(p.ravel(), 1e-100)
            )
            for t, count, p in zip(TIMES, self.obs, pmf, strict=True)
        ]
        return histogram_score(groups)

    def loss(self, vector, changed):
        rates = 10**vector
        before = rates[:4]
        after = before.copy()
        after[list(changed)] = rates[4:]
        pmf, _, _, _ = self.probabilities(before, after)
        s = self.score(pmf)
        return s.cells * s.nats_per_cell


def information(nll, order, cells):
    """AIC and BIC from the total NLL, as in the BayFish figure code."""
    return 2 * nll + 2 * order, 2 * nll + order * np.log(cells)


def corrected_aic(nll, order, cells):
    """The small-sample AIC that the paper's Methods section writes."""
    return 2 * nll + 2 * order + 2 * order * (order + 1) / (cells - order - 1)


def nested_pairs(models, tolerance=0.01):
    """Every (superset, subset) pair and whether the superset fits no worse."""
    by_set = {frozenset(m["changed_indices"]): m for m in models}
    return [
        dict(
            superset=big["changed"],
            subset=small["changed"],
            gap=big["nll"] - small["nll"],
            passed=big["nll"] <= small["nll"] + tolerance,
        )
        for key, big in by_set.items()
        for other, small in by_set.items()
        if other < key
    ]


def options(parser):
    data_option(parser)
    parser.add_argument(
        "--evaluations", type=int, help="Evaluations per model and start."
    )
    parser.add_argument("--retries", type=int)
    parser.add_argument("--nested-evaluations", type=int)
    options_resume(parser)


def run(args, report):
    """Check a joint likelihood and select a stimulus model from observed cells."""
    files = fixtures(args.data_dir, "C6")
    groups = read_counts(files)
    cells = tuple(map(len, groups))
    check(report, "376-archived-cells", cells == TOTALS, by_time=cells)
    limit = 300
    case = StimulusFit(groups, limit)
    large = StimulusFit(groups, 450)
    before, after = TABLE_BASAL, TABLE_STIMULATED
    expected = reference(before, after, limit)
    prediction, flux, leakage, residual = case.probabilities(before, after)
    independent = case.score(expected)
    native = case.score(prediction)
    error = abs(
        sum(c[-2] for c in native.groups) - sum(c[-2] for c in independent.groups)
    )
    expanded, _, _, _ = large.probabilities(before, after)
    cutoff_error = (
        abs(large.score(expanded).nats_per_cell - native.nats_per_cell) * native.cells
    )
    check(
        report,
        "fixed-vector-joint-likelihood",
        error <= 0.05 and cutoff_error <= 0.05,
        native=376 * native.nats_per_cell,
        reference=376 * independent.nats_per_cell,
        discrepancy=error,
        cutoff_difference=cutoff_error,
        leakage=leakage,
        stationary_boundary_flux=flux,
        stationary_residual=residual,
    )
    fast = args.preset == "smoke"
    attempts = args.retries or (1 if fast else 8)
    calls = args.evaluations or (30 if fast else 2000)
    all_models = []
    warm_calls = (
        args.nested_evaluations
        if args.nested_evaluations is not None
        else (0 if fast else 300)
    )
    for j, changed in enumerate(HYPOTHESES):
        label = "-".join(NAMES[i] for i in changed)
        bounds = np.log10([*BOUNDS, *[BOUNDS[i] for i in changed]])
        # A resumed search replays saved scores: tie them to this objective.
        context = objective_context(
            "C6",
            objective="joint count NLL, version 1",
            hypothesis=label,
            rna_limit=limit,
            decay=FIXED_DECAY,
        )
        fits = [
            search(
                args,
                f"{label}-{i}",
                "de",
                lambda x, changed=changed: case.loss(x, changed),
                bounds,
                4497 + i + args.seed_offset,
                calls,
                pop=16 if not fast else 10,
                context=context,
            )
            for i in range(attempts)
        ]
        best = min(fits, key=lambda row: row.value)
        require(best.value < PENALTY, f"no valid fit for {label}")
        competing = [(best.value, best.best)]
        warm = []
        if warm_calls:
            for previous in all_models:
                if set(previous["changed_indices"]) < set(changed):
                    rates = np.r_[
                        previous["before"], np.asarray(previous["after"])[list(changed)]
                    ]
                    x = np.log10(rates)
                    cost = case.loss(x, changed)
                    warm.append(dict(parent=previous["changed"], nll=cost))
                    competing.append((cost, x))
            start = min(competing, key=lambda row: row[0])[1]
            polish = search(
                args,
                f"{label}-cma",
                "cmaes",
                lambda x, changed=changed: case.loss(x, changed),
                bounds,
                44497 + args.seed_offset,
                warm_calls,
                pop=12,
                x0=start,
                sigma=0.04,
                context=context,
            )
            fits.append(polish)
            competing.append((polish.value, polish.best))
        nll_search, chosen = min(competing, key=lambda row: row[0])
        rates = 10**chosen
        before = rates[:4]
        after = before.copy()
        after[list(changed)] = rates[4:]
        pmf, flux, leakage, residual = large.probabilities(before, after)
        score = large.score(pmf)
        nll = score.nats_per_cell * score.cells
        aic, bic = information(nll, len(rates), score.cells)
        retries = [f.value for f in fits[:attempts]]
        all_models.append(
            dict(
                changed=label,
                before=before,
                after=after,
                size=len(rates),
                aic=aic,
                bic=bic,
                aicc=corrected_aic(nll, len(rates), score.cells),
                nll=nll,
                retry_nll=retries,
                retry_spread=max(retries) - min(retries),
                search_nll=nll_search,
                cutoff_difference=nll - nll_search,
                by_time=score.groups,
                leakage=leakage,
                flux=flux,
                residual=residual,
                failures=sum(f.failures for f in fits),
                calls=sum(len(f.scores) for f in fits) + len(warm),
                nested_probes=warm,
                changed_indices=changed,
                index=j,
            )
        )
    all_models.sort(key=lambda m: m["bic"])
    for m in all_models:
        m.update(
            delta_bic=m["bic"] - all_models[0]["bic"],
            delta_aic=m["aic"] - min(row["aic"] for row in all_models),
        )
    pairs = nested_pairs(all_models)
    # A superset worse than its subset means an unconverged search.
    # The smoke budget is too small for this gate; it only records the pairs.
    if fast:
        report["nested_pairs"] = pairs
    else:
        check(report, "nested-hypotheses", all(p["passed"] for p in pairs), pairs=pairs)
    winner = all_models[0]
    pmf, *_ = large.probabilities(winner["before"], winner["after"])
    report.update(
        doi=DOI,
        models=all_models,
        selected=winner["changed"],
        paper_winner="kon-koff",
        max_counts=(limit, 450),
        observed_cells=cells,
        aic_winner=min(all_models, key=lambda m: m["aic"])["changed"],
        degradation=FIXED_DECAY,
    )
    write(args.output / "model.json", model(limit).to_dict())
    write(
        args.output / "distributions.json",
        dict(times=TIMES, predicted=pmf, samples=groups),
    )

    def draw(fig):
        for ax, frame, counts, t in zip(
            fig.subplots(2, 2).flat, pmf, groups, TIMES, strict=True
        ):
            ymax = min(450, max(60, int(np.max(counts[:, 1])) + 12))
            for a, row in enumerate(frame):
                x = counts[counts[:, 0] == a, 1]
                observed = np.bincount(x, minlength=ymax + 1)[: ymax + 1] / len(counts)
                ax.plot(row[: ymax + 1], label=f"{a} active sites")
                ax.scatter(
                    np.arange(len(observed))[observed > 0], observed[observed > 0], s=6
                )
            ax.set(title=f"{t} min", xlabel="RNA copies", ylabel="Joint probability")
        fig.axes[0].legend(fontsize=7)

    save_plot(args, "joint-distributions.png", draw)


SITE_COLORS = ("#2a78d6", "#eb6834", "#1baf7a")


def draw(tables, folder):
    """Draw C6-distributions.png and C6-criteria.png from results/C6 tables."""
    import csv
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tables, folder = Path(tables), Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with (tables / "distributions.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    with (tables / "models.csv").open(newline="") as stream:
        models = [m for m in csv.DictReader(stream) if m["run"] == "delivered"]
    width = 10
    edges = np.arange(0, 301 + width, width)
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True, layout="constrained")
    for ax, t in zip(axes.flat, TIMES, strict=True):
        sub = [row for row in rows if int(row["time"]) == t]
        cells = sum(int(row["observed"]) for row in sub)
        for sites, color in enumerate(SITE_COLORS):
            part = [row for row in sub if int(row["sites"]) == sites]
            rna = np.array([int(row["rna"]) for row in part])
            bins = np.digitize(rna, edges) - 1
            size = len(edges) - 1
            predicted = np.bincount(
                bins, [float(row["predicted"]) for row in part], size
            )
            observed = np.bincount(bins, [int(row["observed"]) for row in part], size)
            centers = edges[:-1] + width / 2
            ax.plot(centers, predicted, color=color, lw=2, label=f"fit, {sites} sites")
            seen = observed > 0
            ax.plot(
                centers[seen],
                observed[seen] / cells,
                "o",
                ms=5,
                mfc="none",
                color=color,
                label=f"cells, {sites} sites",
            )
        ax.set_title(f"{t} min, {cells} cells", fontsize=10)
        ax.set_yscale("log")
        ax.set_ylim(5e-4, 1)
        ax.grid(alpha=0.3, lw=0.5)
    for ax in axes[1]:
        ax.set_xlabel("mature RNA copies per cell (bins of 10)")
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of cells")
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle(
        "Npas4 joint counts: archive cells (circles) and our maximum-likelihood"
        " kon-koff fit (lines)",
        fontsize=10,
    )
    fig.savefig(folder / "C6-distributions.png", dpi=155)
    plt.close(fig)

    order = sorted(models, key=lambda m: float(m["delta_bic"]))
    fig, ax = plt.subplots(figsize=(7, 3.6), layout="constrained")
    y = np.arange(len(order))[::-1]
    for key, marker, color, label in (
        ("delta_aic", "s", SITE_COLORS[0], "ΔAIC (ours)"),
        ("delta_bic", "o", SITE_COLORS[1], "ΔBIC (ours)"),
    ):
        ax.plot(
            [float(m[key]) for m in order],
            y,
            marker,
            ms=8,
            color=color,
            mfc="none" if marker == "s" else color,
            label=label,
            ls="none",
        )
    ax.set_yticks(y, [f"{m['changed']} (k={m['k']})" for m in order])
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlim(-0.3, 2000)
    ax.set_xlabel("difference to the best model (symmetric log scale)")
    ax.axvline(0, color="#52514e", lw=0.8)
    ax.annotate(
        "paper's selection: (k1, k0)",
        (0, y[0]),
        xytext=(3, y[0] - 0.4),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", lw=0.8),
    )
    ax.grid(axis="x", alpha=0.3, lw=0.5)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title(
        "Seven stimulus hypotheses, 376 cells, delivered study fits", fontsize=10
    )
    fig.savefig(folder / "C6-criteria.png", dpi=155)
    plt.close(fig)
