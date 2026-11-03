"""C3: a common CRN objective adapter for native MODE and grid MAP-Elites."""

import copy
import json
import time
import numpy as np
import crn_craft as crn
from crn_craft.observables import oscillation
import fcmaes_rust
from common import ROOT, arm_manifest, call_json, checked, digest, require, write

# Seeds: training paths select designs; validation, the cross-backend control
# and the reference backend each use their own paths.
TRAINING_SEED, VALIDATION_SEED, CONTROL_SEED, REFERENCE_SEED = 8101, 9101, 9201, 9301
# 201 samples over 125 time units: after the 20% burn-in the retained window is
# about 100.6, so the FFT resolves periods 100.6/k, including about 20.1 (k=5).
HORIZON, SAMPLES, BURN = 125.0, 201, 0.2


def network(edges=((2, 0, -1), (0, 1, -1), (1, 2, -1))):
    """One basal channel plus one Hill channel per edge preserves additive rates."""
    model = dict(
        schema="crn-craft.model/0.1",
        name="three_species_additive_hill",
        convention="stochastic_counts",
        omega=1.0,
        parameters={"basal": 0.2, "vmax": 8.0, "k": 8.0, "n": 3.0, "decay": 0.3},
        species=[
            dict(name=f"X{i}", initial_count=6 if i == 0 else 0) for i in range(3)
        ],
        reactions=[],
    )
    for i in range(3):
        model["reactions"] += [
            dict(
                name=f"basal{i}",
                reactants={},
                products={f"X{i}": 1},
                rate=dict(kind="mass_action", coefficient="basal"),
            ),
            dict(
                name=f"decay{i}",
                reactants={f"X{i}": 1},
                products={},
                rate=dict(kind="mass_action", coefficient="decay"),
            ),
        ]
    for j, (source, target, sign) in enumerate(edges):
        require(
            sign in (-1, 1) and source in range(3) and target in range(3),
            "invalid signed edge",
        )
        model["reactions"].append(
            dict(
                name=f"edge{j}",
                reactants={},
                products={f"X{target}": 1},
                rate=dict(
                    kind="hill_activation" if sign == 1 else "hill_repression",
                    regulator=f"X{source}",
                    vmax="vmax",
                    basal=0.0,
                    k="k",
                    n="n",
                ),
            )
        )
    return model


def generator_check():
    # Multiple incoming edges, both signs and zero populations exercise the sum.
    edges = [(2, 0, -1), (0, 1, 1), (2, 1, -1), (1, 2, -1)]
    rng = np.random.default_rng(51)
    largest = 0.0
    for counts in [[0, 0, 0], [8, 8, 8], *rng.integers(0, 60, (20, 3)).tolist()]:
        model = network(edges)
        for s, n in zip(model["species"], counts, strict=True):
            s["initial_count"] = n
        analysis = crn.analyze(model)["analysis"]
        actual_birth = np.zeros(3)
        actual_death = np.zeros(3)
        for reaction, rate in zip(
            analysis["reactions"], analysis["propensities"], strict=True
        ):
            change = reaction["net_changes"][0]
            (actual_birth if change["coefficient"] > 0 else actual_death)[
                change["species"]
            ] += rate
        expected = np.full(3, 0.2)
        for source, target, sign in edges:
            h = (counts[source] / 8.0) ** 3
            expected[target] += 8 * (h / (1 + h) if sign == 1 else 1 / (1 + h))
        largest = max(largest, float(np.max(np.abs(actual_birth - expected))))
        require(
            np.allclose(actual_birth, expected, atol=1e-12, rtol=1e-12),
            "additive Hill generator differs",
        )
        require(
            np.allclose(actual_death, 0.3 * np.asarray(counts)),
            "degradation generator differs",
        )
    return largest


def score(times, paths):
    metrics = []
    for values in paths:
        diagnostic = oscillation(times, values[:, 0], burn_fraction=BURN)
        retained = values[int(len(times) * BURN) :]
        amplitude = float(
            np.quantile(retained[:, 0], 0.9) - np.quantile(retained[:, 0], 0.1)
        )
        period = diagnostic["period"] if diagnostic["period"] is not None else 160.0
        metrics.append(
            dict(
                period=period,
                amplitude=amplitude,
                quality=abs(period - 20) / 20 + 1 - diagnostic["peak_fraction"],
                failed=diagnostic["classification"] != "oscillatory",
                cost=float(retained.sum(axis=1).mean()) / 100,
            )
        )
    periods = np.array([m["period"] for m in metrics])
    objectives = [
        float(np.mean([m["quality"] for m in metrics])),
        float(np.mean([m["cost"] for m in metrics])),
        float(
            np.mean([m["failed"] for m in metrics])
            + np.std(periods) / max(periods.mean(), 1e-12)
        ),
    ]
    descriptors = [
        float(periods.mean()),
        float(np.mean([m["amplitude"] for m in metrics])),
    ]
    return objectives, descriptors, metrics


def oscillatory_fraction(metrics):
    """Share of paths classified oscillatory; descriptors average all paths."""
    return float(np.mean([not m["failed"] for m in metrics])) if metrics else 0.0


def nondominated(rows):
    return [
        r
        for r in rows
        if r["valid"]
        and not any(
            s["valid"]
            and np.all(np.asarray(s["objectives"]) <= r["objectives"])
            and np.any(np.asarray(s["objectives"]) < r["objectives"])
            for s in rows
        )
    ]


def run(args, report):
    """Search a Pareto set and a behavior archive, then measure held-out drift."""
    budget = 32 if args.preset == "smoke" else 256
    training = 4 if args.preset == "smoke" else 16
    validation = 16 if args.preset == "smoke" else 128
    report["budgets"] = dict(
        requested_per_search=budget,
        training_paths=training,
        validation_paths=validation,
        training_seed=TRAINING_SEED,
        validation_seed=VALIDATION_SEED,
        control_seed=CONTROL_SEED,
        reference_seed=REFERENCE_SEED,
        max_events_per_path=200000,
        archive_descriptor_upper=[80, 80],
        horizon=HORIZON,
        sample_count=SAMPLES,
        burn_fraction=BURN,
        resolvable_periods="retained window / k for k = 1, 2, ...; oscillatory needs k >= 3",
        threads=args.threads,
    )
    report["generator_max_error"] = generator_check()
    report["checks"].append(dict(name="additive-generator-equivalence", passed=True))
    model = network()
    write(args.output / "hill.model.json", model)
    times = np.linspace(0, HORIZON, SAMPLES).tolist()
    request = crn.ssa_spec(
        model,
        aggregator="cached_direct",
        end=HORIZON,
        save_times=times,
        replicates=training,
        keep=training,
        seed=TRAINING_SEED,
        threads=args.threads,
        max_events=200000,
    )
    report["budgets"].update(
        request_schema=request["schema"], fast=request["method"]["fast"]
    )
    template = crn.prepare_parameters(model, request)
    test_request = copy.deepcopy(request)
    test_request["ensemble"].update(replicates=validation, root_seed=VALIDATION_SEED)
    test_request["outputs"]["keep_trajectories"] = validation
    test_template = crn.prepare_parameters(model, test_request)
    records = []

    def evaluate(x, holdout=False):
        parameters = dict(model["parameters"])
        parameters.update(
            dict(zip(["vmax", "k", "decay"], np.exp(x).tolist(), strict=True))
        )
        result = (test_template if holdout else template).simulate(parameters)
        valid = (
            result["termination"] == "completed"
            and result["completed_replicates"] == result["requested_replicates"]
            and len(result["trajectories"]) == result["requested_replicates"]
        )
        if valid:
            paths = np.asarray([t["counts"] for t in result["trajectories"]])
            objectives, descriptors, metrics = score(times, paths)
        else:
            objectives, descriptors, metrics = [1e6, 1e6, 1e6], [0, 0], []
        row = dict(
            x=list(x),
            objectives=objectives,
            descriptors=descriptors,
            oscillatory_fraction=oscillatory_fraction(metrics),
            metrics=metrics,
            valid=valid,
            sampler=result.get("sampler"),
            termination=result["termination"],
            events=result["events"],
            requested_paths=result["requested_replicates"],
            completed_paths=result["completed_replicates"],
            holdout=holdout,
        )
        records.append(row)
        return row

    lower, upper = np.log([4, 4, 0.15]), np.log([16, 16, 0.6])
    mode = fcmaes_rust.MODE(3, 3, 0, lower, upper, [False] * 3, popsize=8, seed=718)
    start = time.perf_counter()
    mode_rows = []
    for _ in range(budget // 8):
        xs = mode.ask()
        rows = [evaluate(x) for x in xs]
        mode_rows.extend(rows)
        mode.tell([r["objectives"] for r in rows])
    front = nondominated(mode_rows)
    require(front, "MODE found no complete numerical evaluation")
    arm_manifest(
        args.output,
        "oscillator-mode",
        budget,
        len(mode_rows),
        time.perf_counter() - start,
        718,
        args.threads,
        formulation="mo",
        model="three-species-additive-Hill",
    )
    archive = fcmaes_rust.Archive(
        3,
        lower,
        upper,
        [0, 0],
        [80, 80],
        16,
        samples_per_niche=0,
        seed=719,
        seed_parents=False,
    )
    before = len(records)
    start = time.perf_counter()

    def quality(x):
        row = evaluate(x)
        return (
            sum(row["objectives"]) if row["valid"] else float("inf"),
            row["descriptors"],
        )

    archive.optimize_map_elites(quality, generations=budget // 8, chunk_size=8)
    archive_calls = len(records) - before
    require(archive.occupied > 0, "empty behavior archive")
    arm_manifest(
        args.output,
        "oscillator-qd",
        budget,
        archive_calls,
        time.perf_counter() - start,
        719,
        args.threads,
        formulation="qd",
        model="three-species-additive-Hill",
    )
    selected = [min(front, key=lambda r: r["objectives"][i]) for i in range(3)]
    qd = []
    for i, (x, y, descriptor) in enumerate(
        zip(archive.xs(), archive.ys(), archive.descriptors(), strict=True)
    ):
        if np.isfinite(y):
            require(
                grid_niche(descriptor) == i,
                "Python grid audit and native archive index differ",
            )
            qd.append(
                dict(
                    niche=i,
                    x=x.tolist(),
                    quality=float(y),
                    descriptors=descriptor.tolist(),
                )
            )
    # The archive keeps x and descriptors; look up each elite's training
    # evaluation to report how many of its paths oscillate.
    for elite in qd:
        elite["oscillatory_fraction"] = next(
            (
                r["oscillatory_fraction"]
                for r in records
                if not r["holdout"] and np.array_equal(r["x"], elite["x"])
            ),
            None,
        )
    validation_rows = []
    for origin, row in [("pareto", r) for r in selected] + [("archive", r) for r in qd]:
        held = evaluate(row["x"], True)
        validation_rows.append(
            dict(
                origin=origin,
                training=row,
                validation=held,
                descriptor_distance=float(
                    np.linalg.norm(np.asarray(held["descriptors"]) - row["descriptors"])
                ),
                training_niche=grid_niche(row["descriptors"]),
                heldout_niche=grid_niche(held["descriptors"]),
                training_oscillatory_fraction=row["oscillatory_fraction"],
                heldout_oscillatory_fraction=held["oscillatory_fraction"],
                niche_changed=grid_niche(row["descriptors"])
                != grid_niche(held["descriptors"]),
            )
        )
    report["backend_control"] = backend_control(args, times, model, selected[0]["x"])
    report["checks"].append(
        dict(name="common-scoring-and-independent-backend-law-control", passed=True)
    )
    # The canonical nine-species Vilar model is a separate short reference control.
    vilar = json.loads((ROOT / "models/examples/vilar-noise.model.json").read_text())
    require(
        len(vilar["species"]) == 9 and len(vilar["reactions"]) == 16,
        "Vilar fixture identity changed",
    )
    vrequest = crn.ssa_spec(
        vilar,  # fast policy, the default
        aggregator="cached_direct",
        end=2,
        save_times=[0, 1, 2],
        replicates=2,
        keep=2,
        seed=831,
        max_events=2000000,
    )
    control = checked(crn.simulate_ssa(vilar, vrequest))
    names = control["species"]
    for trajectory in control["trajectories"]:
        values = np.asarray(trajectory["counts"])
        require(
            np.all(values[:, names.index("DA")] + values[:, names.index("DAp")] == 1),
            "Vilar activator promoter conservation",
        )
        require(
            np.all(values[:, names.index("DR")] + values[:, names.index("DRp")] == 1),
            "Vilar repressor promoter conservation",
        )
    write(args.output / "vilar-control.json", control)
    write(args.output / "evaluations.json", records)
    write(args.output / "pareto.json", front)
    write(args.output / "archive.json", qd)
    report.update(
        pareto_points=len(front),
        occupied_niches=archive.occupied,
        validation=validation_rows,
        training_calls=len([r for r in records if not r["holdout"]]),
        validation_calls=len(validation_rows),
        total_events=sum(r["events"] for r in records)
        + control["events"]
        + sum(r["crn_events"] for r in report["backend_control"]["rows"]),
        failed_evaluations=sum(not r["valid"] for r in records),
        clipped_descriptors=sum(
            any(d < 0 or d > 80 for d in r["descriptors"]) for r in records
        ),
        vilar_control="Distinct canonical 9-species/16-reaction model; two-time-unit conservation control, not an oscillation study.",
        scoring="Common score for CRN and the topology tutorial local ReBop compatibility backend; distinct from the older tutorials' score and parameterization. No speed claim from smoke timings.",
    )
    report["checks"].append(
        dict(name="mode-archive-holdout-and-distinct-vilar-control", passed=True)
    )


def grid_niche(descriptor):
    # Fixed 4 by 4 grid, same declared descriptor bounds as the native archive.
    scaled = np.clip(np.asarray(descriptor) / [80, 80], 0, np.nextafter(1.0, 0.0))
    return int(scaled[1] * 4) * 4 + int(scaled[0] * 4)


def backend_control(args, times, model, selected):
    count = 64 if args.preset == "smoke" else 512
    rows = []
    for index, point in enumerate([np.log([8, 8, 0.3]), selected]):
        parameters = dict(model["parameters"])
        parameters.update(
            dict(zip(["vmax", "k", "decay"], np.exp(point).tolist(), strict=True))
        )
        request = dict(
            **parameters,
            edges=[(2, 0, -1), (0, 1, -1), (1, 2, -1)],
            times=times,
            paths=count,
            seed=REFERENCE_SEED,
            initial=[s["initial_count"] for s in model["species"]],
            states=[[0, 0, 0], [8, 8, 8], [1, 40, 3]],
        )
        path = args.output / f"backend-{index}.request.json"
        write(path, request)
        reference = call_json([args.reference_binary.resolve(), path], 30)
        write(args.output / f"backend-{index}.reference.json", reference)
        m = copy.deepcopy(model)
        m["parameters"] = parameters
        for state, births in zip(request["states"], reference["births"], strict=True):
            for sp, n in zip(m["species"], state, strict=True):
                sp["initial_count"] = n
            analysis = crn.analyze(m)["analysis"]
            actual = np.zeros(3)
            for reaction, propensity in zip(
                analysis["reactions"], analysis["propensities"], strict=True
            ):
                change = reaction["net_changes"][0]
                if change["coefficient"] > 0:
                    actual[change["species"]] += propensity
            require(
                np.allclose(actual, births, rtol=1e-12, atol=1e-12),
                "cross-backend summed propensities differ",
            )
        m = copy.deepcopy(model)
        m["parameters"] = parameters
        spec = crn.ssa_spec(
            m,
            end=times[-1],
            save_times=times,
            replicates=count,
            keep=count,
            seed=CONTROL_SEED,
            threads=args.threads,
            aggregator="cached_direct",
            max_events=200000,
        )
        result = checked(crn.simulate_ssa(m, spec))
        paths = np.asarray([p["counts"] for p in result["trajectories"]])
        other = np.asarray(reference["paths"])
        # Path-averaged counts, not correlated samples treated as independent.
        burn = int(len(times) * BURN)
        left, right = paths[:, burn:, :].mean(axis=1), other[:, burn:, :].mean(axis=1)
        se = np.sqrt(
            left.var(axis=0, ddof=1) / count + right.var(axis=0, ddof=1) / count
        )
        z = np.abs(left.mean(axis=0) - right.mean(axis=0)) / np.maximum(se, 1 / count)
        require(np.all(z < 6), "cross-backend path-average law check failed")
        rows.append(
            dict(
                parameters=parameters,
                paths_per_backend=count,
                max_standard_errors=float(max(z)),
                crn_score=score(times, paths)[:2],
                reference_score=score(times, other)[:2],
                crn_events=result["events"],
                crn_sampler=result["sampler"],
                reference_event_count=None,
                reference_backend=reference["backend"],
            )
        )
    return dict(
        binary_sha256=digest(args.reference_binary),
        rows=rows,
        scope="Same six-channel generator aggregated from CRN's nine channels; shared scoring. CRN runs its default fast policy (Xoshiro256++ per replicate, seeded from a keyed stream) on its own seed; the reference seeds SplitMix64 per path from another. Neither side reuses training or validation paths; no path equality claim.",
    )
