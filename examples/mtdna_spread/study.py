"""Simulate spatial mtDNA birth, death and migration with CRN-craft."""

import argparse
import json
from hashlib import sha256
from tempfile import NamedTemporaryFile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import crn_craft as crn
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent
MANIFEST = json.loads((ROOT / "manifest.json").read_text())
# Regions, Nss, delta and hop rate per directed edge (per copy and day). The
# archive headers give delta and the hop rate; Nss = 150 in every panel is a
# reconstruction choice (the per-region spread at day 700 matches the archive).
# Each region starts at h = 0.5 on its steady state W + delta * M = Nss.
CONDITIONS = {
    "deterministic": (2, 150, 0.5, 0.0),
    "isolated": (1, 150, 0.5, 0.0),
    "equal": (2, 150, 1.0, 0.035),
    "dense": (2, 150, 0.5, 0.035),
}


def build(
    initial, copies, delta, hop, *, death=0.07, c=0.001, extra=0.0, fixed_left=False
):
    model = crn.Model("spatial_mtdna")
    model.add_parameter(
        {
            "target": copies,
            "delta": delta,
            "migration": hop,
            "feedback": c,
            "death": death,
        }
    )
    for i, (wild, mutant) in enumerate(initial):
        model.add_species({f"W{i}": int(wild), f"M{i}": int(mutant)})
    for i in range(len(initial)):
        for label, removal in (("W", death), ("M", death + extra)):
            species = f"{label}{i}"
            if not (fixed_left and i == 0):
                birth = f"max(0, death + feedback * (target - W{i} - delta * M{i}))"
                model.add_reaction(
                    crn.mass_action(birth),
                    [species],
                    [species, species],
                    name=f"birth{species}",
                )
                model.add_reaction(removal, [species], [], name=f"death{species}")
            for j in (i - 1, i + 1):
                if 0 <= j < len(initial):
                    old = [] if fixed_left and i == 0 else [species]
                    new = [] if fixed_left and j == 0 else [f"{label}{j}"]
                    rate = (
                        f"migration * {int(initial[i][int(label == 'M')])}"
                        if not old
                        else crn.mass_action("migration")
                    )
                    model.add_reaction(rate, old, new, name=f"move{species}to{j}")
    return model


def build_graph(
    initial, copies, delta, hop, *, death=0.07, c=0.001, extra=0.0, fixed_left=False
):
    """The same two-genotype reactions built from a compartment graph."""
    regions = [str(i) for i in range(len(initial))]
    graph = crn.SpatialModel("spatial_graph_mtdna", dict.fromkeys(regions, 1.0))
    graph.model.add_parameter(
        {"target": copies, "delta": delta, "feedback": c, "death": death}
    )
    active = regions[1:] if fixed_left else regions
    for label, counts in zip(("W", "M"), np.asarray(initial).T, strict=True):
        graph.add_species(
            label,
            dict(zip(regions, map(int, counts), strict=True)),
            reservoirs=regions[:1] if fixed_left else (),
        )
    for label, removal in (("W", death), ("M", death + extra)):
        graph.local(
            f"birth_{label}",
            "max(0, death + feedback * (target - {W} - delta * {M}))",
            [label],
            [label, label],
            combinatorial=True,
            regions=active,
        )
        graph.local(f"death_{label}", removal, [label], [], regions=active)
        for a, b in zip(regions[:-1], regions[1:], strict=True):
            graph.transport(label, a, b, hop).transport(label, b, a, hop)
    return graph


def crossing(x, h, threshold=0.5):
    x, h = np.asarray(x), np.asarray(h)
    if x.shape != h.shape or np.any(np.diff(x) <= 0):
        raise ValueError("A front requires matching, increasing position samples")
    candidates = np.flatnonzero((h[:-1] >= threshold) & (h[1:] < threshold))
    if len(candidates) != 1:
        raise ValueError(f"Expected one descending crossing, got {len(candidates)}")
    i = candidates[0]
    return float(x[i] + (x[i + 1] - x[i]) * (threshold - h[i]) / (h[i + 1] - h[i]))


def speeds(profiles, strict=True):
    ages = sorted(profiles.day.unique())
    out = []
    for threshold in MANIFEST["fig2_crossings"]:
        positions = []
        for _, group in profiles.sort_values("position_um").groupby("day"):
            try:
                positions.append(crossing(group.position_um, group["mean"], threshold))
            except ValueError:
                if strict:
                    raise
                positions.append(float("nan"))
        out.append(
            (threshold, *positions, (positions[1] - positions[0]) / (ages[1] - ages[0]))
        )
    return pd.DataFrame(out, columns=["threshold", "x0_um", "x1_um", "speed_um_day"])


def fractions(counts):
    total = counts.sum(axis=-1)
    return np.divide(
        counts[..., 1], total, out=np.zeros_like(total, dtype=float), where=total > 0
    )


def paths(initial, params, times, n, seed, threads, **kwargs):
    nsites = len(initial)
    model = build(initial, *params, **kwargs)
    result = crn.ensemble(
        model,
        end=float(times[-1]),
        report_times=list(times),
        replicates=n,
        seed=seed,
        threads=threads,
        method="static_direct",
        keep_paths="all",
        strict=True,
        max_output_bytes=256 * 2**20,
        max_events=3_000_000_000 if nsites >= 100 else 100_000_000,
    ).simulate()
    labels = [f"{allele}{i}" for i in range(nsites) for allele in ("W", "M")]
    return result.paths.sel(variable=labels).values.reshape(n, len(times), nsites, 2)


def figure1(size, threads):
    n = 128 if size == "quick" else 2000
    times = np.array(
        [0, 5, 20.0] if size == "quick" else MANIFEST["fig1_times"], dtype=float
    )
    rows = []
    for index, (name, (regions, copies, delta, hop)) in enumerate(CONDITIONS.items()):
        for treatment, extra in (("neutral", 0), ("extra", 0.00007)):
            if name == "deterministic":
                h = 1 / (1 + np.exp(extra * times))
                err, empty, actual = np.zeros_like(h), np.zeros_like(h), 0
            else:
                half = round(copies / (1 + delta))
                p = paths(
                    [(half, half)] * regions,
                    (copies, delta, hop),
                    times,
                    n,
                    MANIFEST["seed"] + index * 100 + int(extra > 0),
                    threads,
                    extra=extra,
                )
                f = fractions(p).mean(axis=-1)
                h, err, empty, actual = (
                    f.mean(axis=0),
                    f.std(axis=0, ddof=1) / np.sqrt(n),
                    (p.sum(-1) == 0).sum(axis=(0, 2)),
                    n,
                )
            rows.extend(
                (name, treatment, t, v, e, z, actual)
                for t, v, e, z in zip(times, h, err, empty, strict=True)
            )
    return pd.DataFrame(
        rows,
        columns=[
            "condition",
            "degradation",
            "day",
            "mean",
            "sem",
            "empty_regions",
            "replicates",
        ],
    )


def terminal_checks(comparison):
    rows = comparison.query(
        "day == 700 and degradation == 'neutral' and replicates > 0"
    ).copy()
    if rows.empty:
        return rows, {}
    z = norm.ppf(1 - 0.01 / (2 * 4))
    error = np.hypot(rows["sem"], rows["sem_upper"])
    rows["difference"] = rows["mean"] - rows["mean_reference"]
    rows["low"] = rows["difference"] - z * error
    rows["high"] = rows["difference"] + z * error
    rows["consistent"] = (rows.low <= 0) & (rows.high >= 0)
    contrast = rows.set_index("condition").loc[["dense", "equal"]]
    gap = contrast["mean"].iloc[0] - contrast["mean"].iloc[1]
    lower = gap - z * np.sqrt(np.sum(contrast["sem"] ** 2))
    return rows, {"dense_minus_equal": float(gap), "contrast_lower": float(lower)}


def chain_initial(p):
    x = p["spacing"] * np.arange(p["sites"])
    h = 1 / (1 + np.exp(p["initial_slope"] * x - p["initial_intercept"]))
    h[0] = 1
    h[h < p["right_cutoff"]] = 0
    target = p["copies"] / (1 - h + p["delta"] * h)
    return x, np.column_stack([target * (1 - h), target * h]).round().astype(int)


def small_chain(size, threads):
    sites = 12 if size == "quick" else 24
    params = {
        "sites": sites,
        "spacing": 30,
        "initial_slope": 0.025,
        "initial_intercept": 0.375 * sites,
        "right_cutoff": 0.005,
        "copies": 80,
        "delta": 0.5,
    }
    x, counts = chain_initial(params)
    times = [0.0, 25.0, 50.0] if size == "quick" else [0.0, 50.0, 100.0, 200.0]
    sample = paths(
        counts,
        (80, 0.5, 0.035),
        times,
        24 if size == "quick" else 128,
        MANIFEST["seed"] + 500,
        threads,
        fixed_left=True,
    )
    values = fractions(sample).mean(axis=0)
    return pd.DataFrame(
        [
            (t, d, f)
            for t, row in zip(times, values, strict=True)
            for d, f in zip(x, row, strict=True)
        ],
        columns=["day", "position_um", "mean"],
    )


def full_chain(output, threads, fibers=800, sampler="static_direct", start=0):
    """Exact source-scale route with independent, restartable single-fiber jobs."""
    p = MANIFEST["fig2_site_parameters"]
    x, initial = chain_initial(p)
    fingerprint = sha256(
        json.dumps(
            {
                "parameters": p,
                "horizon": 1480,
                "seed": MANIFEST["seed"],
                "method": "chacha8-u53-ziggurat/spatial-turnover/2",
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()[:16]
    checkpoint = output / (
        "fiber-checkpoints"
        if sampler == "static_direct"
        else f"local-fiber-checkpoints-{fingerprint}"
    )
    checkpoint.mkdir(parents=True, exist_ok=True)
    if fibers < 1 or start < 0:
        raise ValueError("At least one fiber and a nonnegative start are required")

    def one(i):
        file = checkpoint / f"{i:04d}.npy"
        if not file.exists():
            if sampler == "local":
                result = crn.spatial_turnover(
                    initial,
                    end=1480,
                    target=p["copies"],
                    feedback=p["c"],
                    mutant_weight=p["delta"],
                    death=p["death"],
                    hop=p["hop"],
                    fixed_left=True,
                    seed=MANIFEST["seed"] + 1000 + i,
                )
                counts = np.stack(
                    (initial, np.stack((result["wild"], result["mutant"]), axis=1))
                )
                sample = fractions(counts)
            else:
                counts = paths(
                    initial,
                    (p["copies"], p["delta"], p["hop"]),
                    [0.0, 1480.0],
                    1,
                    MANIFEST["seed"] + 1000 + i,
                    1,
                    death=p["death"],
                    c=p["c"],
                    fixed_left=True,
                )
                sample = fractions(counts)[0]
            with NamedTemporaryFile(
                dir=checkpoint, suffix=".npy", delete=False
            ) as temporary:
                tmp = Path(temporary.name)
                try:
                    np.save(temporary, sample)
                except BaseException:
                    tmp.unlink(missing_ok=True)
                    raise
            tmp.replace(file)
        saved = np.load(file, allow_pickle=False)
        if saved.shape != (2, len(initial)) or not np.isfinite(saved).all():
            raise ValueError(f"invalid fiber checkpoint {file}")
        return saved

    with ThreadPoolExecutor(max_workers=threads) as pool:
        sample = np.asarray(list(pool.map(one, range(start, start + fibers))))
    mean = sample.mean(axis=0)
    frame = pd.DataFrame(
        [
            (day, pos, h)
            for day, row in zip((0, 1480), mean, strict=True)
            for pos, h in zip(x, row, strict=True)
        ],
        columns=["day", "position_um", "mean"],
    )
    rng = np.random.default_rng(MANIFEST["seed"] + 2000)
    estimates, ambiguous = [], 0
    for weights in rng.multinomial(fibers, np.full(fibers, 1 / fibers), size=500):
        profile = np.einsum("i,ijk->jk", weights / fibers, sample)
        try:
            estimates.append((crossing(x, profile[1]) - crossing(x, profile[0])) / 1480)
        except ValueError:
            ambiguous += 1
    interval = (
        np.quantile(estimates, [0.005, 0.995]).tolist() if estimates else [None, None]
    )
    reference_speed = (
        speeds(pd.read_csv(ROOT / "reference-2e.csv")).iloc[1].speed_um_day
    )
    margin = MANIFEST["fig2_half_height_margin_um_day"]
    status = "inconclusive"
    if not ambiguous and interval[0] is not None:
        if (
            reference_speed - margin
            <= interval[0]
            <= interval[1]
            <= reference_speed + margin
        ):
            status = "equivalent"
        elif (
            interval[1] < reference_speed - margin
            or interval[0] > reference_speed + margin
        ):
            status = "discrepancy"
    return frame, {
        "interval_99": interval,
        "ambiguous_bootstraps": ambiguous,
        "fibers": fibers,
        "first_fiber": start,
        "sampler": sampler,
        "checkpoint_directory": checkpoint.name,
        "status": status,
    }


def plot(results, folder, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        2, 2, figsize=(8, 5.8), sharex=True, sharey=True, layout="constrained"
    )
    for (name, group), ax in zip(
        results.groupby("condition", sort=False), axes.flat, strict=True
    ):
        for kind, color in (("neutral", "#234a74"), ("extra", "#bc4b34")):
            data = group[group.degradation == kind]
            ax.plot(
                data.day,
                data["mean"],
                color=color,
                label="neutral" if kind == "neutral" else "mutants degraded faster",
            )
            ax.errorbar(
                data.day,
                data.mean_reference,
                yerr=data.sem_upper,
                color=color,
                fmt="o",
                ms=3,
                alpha=0.6,
            )
        ax.set(
            title=name.capitalize(), xlabel="Time (days)", ylabel="Mean mutant fraction"
        )
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(folder / "two-regions.png", dpi=155)
    plt.close(fig)
    archive = pd.read_csv(ROOT / "reference-2e.csv")
    fig, ax = plt.subplots(figsize=(7.6, 4.2), layout="constrained")
    for day, group in archive.groupby("day"):
        ax.plot(
            group.position_um, group["mean"], label=f"Authors’ profile, day {day:,}"
        )
    ax.axhline(0.5, color="0.45", linestyle="--", linewidth=1)
    ax.set(xlabel="Position in archived window (µm)", ylabel="Mean mutant fraction")
    ax.legend(frameon=False)
    fig.savefig(folder / "wave-reference.png", dpi=155)
    plt.close(fig)

    candidate = output / "approx-wave-profiles.csv"
    if candidate.exists():
        fresh = pd.read_csv(candidate)
        fig, ax = plt.subplots(figsize=(7.6, 4.2), layout="constrained")
        archive_initial = archive[archive.day == 0]
        archived_origin = crossing(archive_initial.position_um, archive_initial["mean"])
        archived_end = archive[archive.day == 1480]
        ax.plot(
            archived_end.position_um - archived_origin,
            archived_end["mean"],
            color="#2f4c6b",
            linewidth=2,
            label="Archived day 1,480",
        )
        initial = fresh[fresh.day == 0]
        origin = crossing(initial.position_um, initial["mean"])
        ax.plot(
            initial.position_um - origin,
            initial["mean"],
            color=".52",
            linestyle="--",
            label="Fresh day 0",
        )
        for dt, final in fresh[fresh.day > 0].groupby("step_day"):
            ax.plot(
                final.position_um - origin,
                final["mean"],
                label=f"Approximate, {dt:g} day step",
            )
        ax.axhline(0.5, color=".4", linewidth=0.7)
        ax.set(
            xlim=(-1800, 2400),
            xlabel="Distance from each day-zero front (µm)",
            ylabel="Mean mutant fraction",
        )
        ax.legend(frameon=False, fontsize=8)
        fig.savefig(folder / "wave-approximation.png", dpi=155)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", choices=("quick", "draft", "full"), default="quick")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--figures", type=Path)
    parser.add_argument(
        "--wave-sampler", choices=("static_direct", "local"), default="static_direct"
    )
    parser.add_argument("--fiber-start", type=int, default=0)
    parser.add_argument(
        "--source-chain",
        action="store_true",
        help="Run the exact 800-fiber source-scale cohort",
    )
    parser.add_argument(
        "--fibers",
        type=int,
        default=800,
        help="Count for a checkpointed source-scale cohort",
    )
    args = parser.parse_args()
    args.output = args.output or ROOT / (
        "results" if args.figures else f"runs/{args.size}"
    )
    args.output.mkdir(parents=True, exist_ok=True)
    if args.figures:
        plot(pd.read_csv(args.output / "two-regions.csv"), args.figures, args.output)
        return
    reference = pd.read_csv(ROOT / "reference-1c.csv").rename(
        columns={"mean": "mean_reference"}
    )
    data = figure1(args.size, args.threads).merge(reference, how="left")
    data.to_csv(args.output / "two-regions.csv", index=False)
    archive = pd.read_csv(ROOT / "reference-2e.csv")
    speeds(archive).to_csv(args.output / "wave-reference-speeds.csv", index=False)
    small_chain(args.size, args.threads).to_csv(
        args.output / "chain-demo.csv", index=False
    )
    table, summary = terminal_checks(data)
    if not table.empty:
        table.to_csv(args.output / "terminal-checks.csv", index=False)
    if args.source_chain or args.size == "full":
        frame, review = full_chain(
            args.output, args.threads, args.fibers, args.wave_sampler, args.fiber_start
        )
        frame.to_csv(args.output / "wave-fresh.csv", index=False)
        speeds(frame, strict=False).to_csv(
            args.output / "wave-fresh-speeds.csv", index=False
        )
        summary["wave"] = review
    print(json.dumps({"size": args.size, **summary}, sort_keys=True))


if __name__ == "__main__":
    main()
