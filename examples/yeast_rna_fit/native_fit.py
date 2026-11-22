"""Reuse prepared count projections for the chapter's forced RNA likelihood."""

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import json
from pathlib import Path
from time import perf_counter

import crn_craft as crn
from crn_craft.likelihood import HistogramGroup, histogram_score
import numpy as np

from fitting import Data, observations, one_fit, spatial_laws
from model import KEYS, ROOT, exact_means, hog1, nonspatial, parameters


class Objective:
    """Prepared topology and data for both doses and one RNA readout."""

    def __init__(self, anchor, data, limits, tolerance, max_step):
        self.joint = data.spatial
        self.data = data
        self.limits = limits
        self.tolerance = tolerance
        self.max_step = max_step
        self.minutes = np.unique([0, *[t for _, _, t, *_ in data.groups]])
        model = crn.Model("forced_rna")
        promoter = [f"S{k + 1}" for k in range(4)]
        model.add_species({p: int(i == 0) for i, p in enumerate(promoter)})
        model.add_species("N", 0)
        if self.joint:
            model.add_species("C", 0)
        for name, value in zip(KEYS, anchor, strict=True):
            model.add_parameter(name, float(value))
        for rate, source, target in (
            ("k12", 0, 1),
            ("k21", 1, 0),
            ("k23", 1, 2),
            ("k32", 2, 1),
            ("k34", 2, 3),
            ("k43", 3, 2),
        ):
            model.add_reaction(
                crn.mass_action(rate), [promoter[source]], [promoter[target]], name=rate
            )
        for i, p in enumerate(promoter, 1):
            model.add_reaction(crn.mass_action(f"i{i}"), [p], [p, "N"], name=f"rna{i}")
        model.add_reaction(
            crn.mass_action("nuclear_decay" if self.joint else "decay"),
            ["N"],
            [],
            name="loss_n",
        )
        bounds = {**dict.fromkeys(promoter, 1), "N": limits[0]}
        if self.joint:
            bounds["C"] = limits[1]
            model.add_reaction(
                crn.mass_action("export"), ["N"], ["C"], name="transport"
            )
            model.add_reaction(crn.mass_action("decay"), ["C"], [], name="loss_c")
        self.graph = crn.prepare_fsp(model, bounds)
        n = self.graph.states[:, self.graph.species.index("N")]
        if self.joint:
            c = self.graph.states[:, self.graph.species.index("C")]
            self.keys = c * (limits[0] + 1) + n
            self.shape = (limits[1] + 1, limits[0] + 1)
        else:
            self.keys = n
            self.shape = (limits[0] + 1,)
        self.all_bins = int(np.prod(self.shape))
        self.histograms = {}
        for dose, rep, minute, _, observed, weights, _, _ in data.groups:
            unique, indices = np.unique(observed, axis=0, return_inverse=True)
            total = np.zeros(len(unique), dtype=np.int64)
            np.add.at(total, indices, weights)
            self.histograms[dose, rep, minute] = (unique, total)

    def laws(self, theta, *, solver="bdf"):
        laws, stats = {}, {}
        for dose in self.data.dose:

            def rate(t, dose=dose):
                return max(0.0, theta[1] - theta[2] * float(hog1(t, dose, theta[12])))

            result = self.graph.solve(
                float(self.minutes[-1] * 60),
                self.minutes * 60,
                parameters=dict(zip(KEYS, map(float, theta), strict=True)),
                signals={"k21": crn.Signal(rate, change_times=(theta[12],))},
                tolerance=self.tolerance,
                max_step=self.max_step,
                solver=solver,
            )
            laws[dose] = {
                minute: np.bincount(
                    self.keys, snapshot.probabilities, minlength=self.all_bins
                ).reshape(self.shape)
                for minute, snapshot in zip(self.minutes, result.snapshots, strict=True)
            }
            stats[dose] = {
                "work": result.work,
                "leakage": max(s.physical_leakage for s in result.snapshots),
                "mass_defect": max(abs(s.mass_defect) for s in result.snapshots),
                "min_probability": min(s.minimum_probability for s in result.snapshots),
                "sink": dict(
                    zip(
                        map(int, self.minutes),
                        (s.physical_leakage for s in result.snapshots),
                        strict=True,
                    )
                ),
            }
        return laws, stats

    def score(self, theta, *, groups=None, legacy_tail=True):
        return self.score_laws(
            *self.laws(theta), groups=groups, legacy_tail=legacy_tail
        )

    def score_laws(self, laws, diagnostics, *, groups=None, legacy_tail=True):
        samples = []
        for dose, rep, minute, _, xs, weights, _, _ in (
            self.data.groups if groups is None else groups
        ):
            xs, weights = self.histograms[dose, rep, minute]
            law = laws[dose][minute]
            if self.joint:
                inside = (xs[:, 0] < law.shape[1]) & (xs[:, 1] < law.shape[0])
                observed = law[xs[inside, 1], xs[inside, 0]]
            else:
                inside = xs[:, 0] < len(law)
                observed = law[xs[inside, 0]]
            missing = min(1.0, max(0.0, diagnostics[dose]["sink"][minute]))
            count = int(weights[~inside].sum())
            probs = np.r_[observed, missing if legacy_tail else 0.0]
            samples.append(
                HistogramGroup(
                    str(dose),
                    rep,
                    minute,
                    np.r_[weights[inside], count],
                    probs,
                    missing,
                    "minute",
                    "joint" if self.joint else "total",
                )
            )
        return histogram_score(samples)

    def means(self, theta, groups=None):
        means = {
            d: exact_means(theta, d, self.minutes, self.joint) for d in self.data.dose
        }
        out = 0.0
        for dose, _, minute, _, _, _, mean, inverse in (
            self.data.groups if groups is None else groups
        ):
            error = mean - means[dose][np.searchsorted(self.minutes, minute)]
            out += float(error @ inverse @ error / 2)
        return out


def gap(p, q, spatial):
    if not spatial:
        return float(np.max(abs(np.cumsum(p.ravel() - q.ravel()))))
    return max(float(np.max(abs(np.cumsum((p - q).sum(axis=i))))) for i in (0, 1))


def check(folder):
    anchors = json.loads((ROOT / "source-parameters.json").read_text())["parameters"]
    output = {}
    for joint, name, limits in (
        (False, "NonSpatial-Distributions", (360,)),
        (True, "Spatial-Distributions", (28, 305)),
    ):
        x = parameters(anchors[name])
        data = Data(observations(), spatial_model=joint)
        objective = Objective(x, data, limits, 4e-7 if joint else 1e-8, 20)
        start = perf_counter()
        actual, diagnostics = objective.laws(x)
        elapsed = perf_counter() - start
        if joint:
            expected = spatial_laws(x, data, *limits, 4e-7)
        else:
            expected = {
                d: dict(enumerate(nonspatial(x, d, limits[0], 1.5).sum(-1)))
                for d in data.dose
            }
        difference = {
            str(d): max(
                gap(actual[d][t], expected[d][t], joint) for t in objective.minutes
            )
            for d in data.dose
        }
        joint_tv = None
        if joint:
            joint_tv = {
                str(d): max(
                    float(
                        (
                            np.abs(actual[d][t] - expected[d][t]).sum()
                            + abs(actual[d][t].sum() - expected[d][t].sum())
                        )
                        / 2
                    )
                    for t in objective.minutes
                )
                for d in data.dose
            }
        mean_error = {}
        for d in data.dose:
            reference = exact_means(x, d, objective.minutes, joint)
            if joint:
                n = np.arange(limits[0] + 1)
                c = np.arange(limits[1] + 1)
                expected_means = [
                    [p.sum(axis=0) @ n, p.sum(axis=1) @ c] for p in actual[d].values()
                ]
            else:
                expected_means = [[p @ np.arange(len(p))] for p in actual[d].values()]
            mean_error[str(d)] = float(
                np.max(np.abs(reference - np.asarray(expected_means)))
            )
        score = objective.score_laws(actual, diagnostics)
        estimate = 2.65027 if joint else 2.04126
        if (
            max(difference.values()) > 0.001
            or abs(score.nats_per_cell - estimate) > 0.01
        ):
            raise RuntimeError(
                f"paper replay exceeds targets: {name} {difference} {score.nats_per_cell}"
            )
        if not joint and max(v["leakage"] for v in diagnostics.values()) > 1e-5:
            raise RuntimeError("total RNA projection leakage exceeds target")
        output[name] = {
            "states": len(objective.graph.states),
            "cdf_error": difference,
            "joint_total_variation": joint_tv,
            "mean_absolute_error": mean_error,
            "nll": score.nats_per_cell,
            "cells": score.cells,
            "groups": len(score.groups),
            "diagnostics": diagnostics,
            "seconds": elapsed,
        }
        (folder / "paper-check.json").write_text(json.dumps(output, indent=2) + "\n")
    return output


def fit(
    folder,
    budget,
    restarts,
    workers,
    joint=True,
    holdout=False,
    *,
    draft_joint=False,
    training_only=False,
):
    anchors = json.loads((ROOT / "source-parameters.json").read_text())["parameters"]
    name = "Spatial-Distributions" if joint else "NonSpatial-Distributions"
    data = Data(observations(), spatial_model=joint)
    x = parameters(anchors[name])
    if training_only:
        # Preset in rate units, independent of the authors' fit or any held-out
        # measurements. The model topology and input waveform remain the same.
        x = np.array(
            [
                0.015,
                0.05,
                0.5,
                0.025,
                0.008,
                0.013,
                0.01,
                0.001,
                0.03,
                0.45,
                1.5,
                0.01,
                120,
                0.04,
                0.025,
            ]
        )
    limits = ((11, 120) if draft_joint else (21, 240)) if joint else (280,)
    refined = (28, 305) if joint else (360,)
    objective = Objective(
        x, data, limits, 2e-5 if draft_joint else 2e-6, 60 if draft_joint else 30
    )
    refined_obj = Objective(x, data, refined, 4e-7, 20)
    widths = np.full(15 if joint else 13, np.log(1.35 if joint else 3.0))
    if np.any(x[: len(widths)] <= 0):
        raise ValueError("the log fit requires a positive anchor")
    rng = np.random.default_rng(230521)

    def transform(z):
        values = x.copy()
        values[: len(z)] *= np.exp(z)
        return values

    if joint and not draft_joint:
        probe = [np.zeros(15)]
        for a, b in product((-1, 0, 1), repeat=2):
            z = np.zeros(15)
            z[[9, 11]] = np.log(1.35) * np.array([a, b])
            probe.append(z)
        probe.append(rng.uniform(-widths, widths))
        rows = []
        for z in probe:
            start = perf_counter()
            trial = transform(z)
            score = objective.score(trial)
            rows.append(
                {
                    "z": z.tolist(),
                    "seconds": perf_counter() - start,
                    "nll": score.nats_per_cell
                    if np.isfinite(score.nats_per_cell)
                    else None,
                }
            )
            (folder / "joint-pilot.json").write_text(json.dumps(rows, indent=2) + "\n")
        if not budget:
            return rows

    splits = [(None, data.groups, [])]
    if holdout:
        for label in ("20140429", "20160509"):
            splits.append(
                (
                    label,
                    [row for row in data.groups if row[1] != label],
                    [row for row in data.groups if row[1] == label],
                )
            )

    jobs = [
        (label, train, test, method, 230531 + i)
        for label, train, test in splits
        for method in ("histogram", "mean")
        for i in range(restarts)
    ]

    def optimize(task):
        label, train, test, method, seed = task

        def cost(z):
            values = transform(z)
            try:
                return (
                    objective.score(values, groups=train).nats_per_cell
                    if method == "histogram"
                    else objective.means(values, train)
                )
            except (ValueError, RuntimeError, FloatingPointError):
                return np.inf

        z, work = one_fit(cost, widths, budget, seed)
        values = transform(z)
        predictions = refined_obj.laws(values)
        train_score = refined_obj.score_laws(*predictions, groups=train).nats_per_cell
        test_score = (
            refined_obj.score_laws(*predictions, groups=test).nats_per_cell
            if test
            else None
        )
        return {
            "held_out": label,
            "objective": method,
            "seed": seed,
            "work": work,
            "parameters": dict(zip(KEYS, values.tolist(), strict=True)),
            "training_nats": train_score if np.isfinite(train_score) else None,
            "held_out_nats": test_score
            if test_score is not None and np.isfinite(test_score)
            else None,
            "free": len(widths),
            "coarse_bounds": limits,
            "refined_bounds": refined,
            "training_only_seed_center": training_only,
        }

    with ThreadPoolExecutor(max_workers=workers) as pool:
        from concurrent.futures import as_completed

        results = []
        for task in as_completed([pool.submit(optimize, job) for job in jobs]):
            results.append(task.result())
            results.sort(
                key=lambda row: (row["held_out"] or "", row["objective"], row["seed"])
            )
            (folder / "refits.json").write_text(json.dumps(results, indent=2) + "\n")
    return results


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("check", "fit"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--budget", type=int, default=0)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--total", action="store_true")
    parser.add_argument("--holdout", action="store_true")
    parser.add_argument(
        "--draft-joint",
        action="store_true",
        help="fit a coarse joint projection, then rescore with the refined box",
    )
    parser.add_argument(
        "--training-only",
        action="store_true",
        help="use a preset independent of the authors' fit",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.operation == "check":
        check(args.output)
    else:
        fit(
            args.output,
            args.budget,
            args.restarts,
            args.workers,
            not args.total,
            args.holdout,
            draft_joint=args.draft_joint,
            training_only=args.training_only,
        )


if __name__ == "__main__":
    main()
