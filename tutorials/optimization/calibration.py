"""C1: complete fresh-versus-bound parameter fits with independent moments."""

import copy
import json
import time
import numpy as np
import crn_craft._legacy as crn
import fcmaes_rust
from common import ROOT, arm_manifest, checked, require, write


START, LOWER, UPPER = 1.25, 0.35, 2.5  # multiples of the true rates


def options(parser):
    parser.add_argument(
        "--include-checked",
        action="store_true",
        help="add a separate arm with the checked exact policy (fast=False)",
    )


def binding_model():
    return dict(
        schema="crn-craft.model/0.1",
        name="reversible_binding",
        convention="stochastic_counts",
        omega=1.0,
        parameters={"forward": 0.4, "reverse": 0.2},
        species=[dict(name="A", initial_count=30), dict(name="B", initial_count=0)],
        reactions=[
            dict(
                name=name,
                reactants={a: 1},
                products={b: 1},
                rate=dict(kind="mass_action", coefficient=parameter),
            )
            for name, a, b, parameter in [
                ("bind", "A", "B", "forward"),
                ("release", "B", "A", "reverse"),
            ]
        ],
    )


def moments(kind, values, times):
    a, b = values
    t = np.asarray(times)
    if kind == "immigration":
        mean = (a / b * -np.expm1(-b * t))[:, None]
        return mean, mean.copy()
    p = a / (a + b) * -np.expm1(-(a + b) * t)
    return np.column_stack((30 * (1 - p), 30 * p)), np.repeat(
        (30 * p * (1 - p))[:, None], 2, axis=1
    )


def loss(result, target):
    statistics = result["statistics"]
    mean, variance = target
    return float(
        np.mean(((np.asarray(statistics["mean"]) - mean) / (1 + mean)) ** 2)
        + np.mean(
            ((np.asarray(statistics["sample_variance"]) - variance) / (1 + variance))
            ** 2
        )
    )


def run(args, report):
    """Fit two analytic fixtures; report estimator error separately from parameter error."""
    budget = 32 if args.preset == "smoke" else 128
    training = 16 if args.preset == "smoke" else 32
    holdout = 1024 if args.preset == "smoke" else 4096
    seeds = [7001] if args.preset == "smoke" else [7001, 7002, 7003]
    seeds = [seed + args.seed_offset for seed in seeds]
    times = np.linspace(0, 6, 7).tolist()
    report["budgets"] = dict(
        requested_calls=budget,
        training_paths=training,
        holdout_paths=holdout,
        optimizer_seeds=seeds,
        training_seed=77201,
        validation_seed=99173,
        # One training sample serves every optimizer seed: the seeds replicate
        # the search, not the data.
        training_sample="one ensemble, seed 77201, shared by all optimizer seeds",
        start_multiple_of_truth=START,
        bounds_multiple_of_truth=[LOWER, UPPER],
        policies=["fast", "checked"] if args.include_checked else ["fast"],
        event_cap_per_path=100000,
        threads=args.threads,
    )
    for kind, model in [
        (
            "immigration",
            json.loads(
                (ROOT / "models/examples/immigration-death.model.json").read_text()
            ),
        ),
        ("binding", binding_model()),
    ]:
        names = sorted(model["parameters"])
        truth = np.array([model["parameters"][name] for name in names])
        target = moments(kind, truth, times)
        write(args.output / f"{kind}.model.json", model)
        for fast in [True, False] if args.include_checked else [True]:
            policy = "fast" if fast else "checked"
            setup = time.perf_counter()
            request = crn.ssa_spec(
                model,
                fast=fast,
                aggregator="cached_direct",
                end=6,
                save_times=times,
                replicates=training,
                seed=77201,
                threads=args.threads,
                max_events=100000,
            )
            template = crn.prepare_parameters(model, request)
            require(
                list(template.parameter_names) == names,
                "parameter coordinate order changed",
            )
            setup_seconds = time.perf_counter() - setup
            # The training loss of the true rates on the same sample tells
            # overfitting apart from search failure.
            at_truth = checked(template.simulate(truth))
            require(
                at_truth["method"]["fast"] is fast
                and at_truth["schema"] == "crn-craft.ssa-result/0.3",
                f"{kind}-{policy}: result does not state the requested policy",
            )
            truth_training_loss = loss(at_truth, target)
            for seed in seeds:
                comparison = []
                for route in ["fresh", "reusable"]:
                    trace = []
                    work = dict(events=0, paths=0)

                    def objective(
                        x,
                        route=route,
                        template=template,
                        model=model,
                        names=names,
                        fast=fast,
                        target=target,
                        work=work,
                        trace=trace,
                    ):
                        rates = np.exp(x)
                        if route == "reusable":
                            result = checked(template.simulate(rates))
                        else:
                            candidate = copy.deepcopy(model)
                            candidate["parameters"] = dict(
                                zip(names, rates.tolist(), strict=True)
                            )
                            spec = crn.ssa_spec(
                                candidate,
                                fast=fast,
                                aggregator="cached_direct",
                                end=6,
                                save_times=times,
                                replicates=training,
                                seed=77201,
                                threads=args.threads,
                                max_events=100000,
                            )
                            result = checked(crn.simulate_ssa(candidate, spec))
                        value = loss(result, target)
                        work["events"] += result["events"]
                        work["paths"] += result["completed_replicates"]
                        trace.append(dict(x=list(x), value=value))
                        return value

                    start = time.perf_counter()
                    x0 = np.log(truth * START)
                    answer = fcmaes_rust.optimize_bite(
                        objective,
                        x0,
                        np.log(truth * LOWER),
                        np.log(truth * UPPER),
                        seed=seed,
                        max_evaluations=budget,
                    )
                    search_seconds = time.perf_counter() - start
                    x, value, actual = answer[:3]
                    require(actual == len(trace), "optimizer call accounting differs")
                    validation_start = time.perf_counter()
                    validation_request = crn.ssa_spec(
                        model,
                        fast=fast,
                        aggregator="cached_direct",
                        end=6,
                        save_times=times,
                        replicates=holdout,
                        seed=99173,
                        threads=args.threads,
                        max_events=100000,
                    )
                    validation_template = crn.prepare_parameters(
                        model, validation_request
                    )
                    fitted = checked(validation_template.simulate(np.exp(x)))
                    reference = checked(validation_template.simulate(truth))
                    exact_fitted = moments(kind, np.exp(x), times)[0]
                    empirical = np.asarray(fitted["statistics"]["mean"])
                    se = np.asarray(fitted["statistics"]["mean_standard_error"])
                    error = float(
                        np.max(
                            np.abs(empirical - exact_fitted)
                            / np.maximum(se, 1 / holdout)
                        )
                    )
                    require(
                        error < 6,
                        "held-out means disagree with independent exact moments",
                    )
                    bound = template.bind(np.exp(x))
                    replay = checked(
                        crn.simulate_ssa(bound.model_bytes, bound.request_bytes)
                    )
                    repeat = checked(bound.simulate())
                    require(
                        replay == repeat, "ordinary replay and reusable binding differ"
                    )
                    tag = f"{kind}-{policy}-{seed}-{route}"
                    write(args.output / f"{tag}.trace.json", trace)
                    validation_work = dict(
                        events=fitted["events"]
                        + reference["events"]
                        + replay["events"]
                        + repeat["events"],
                        paths=2 * holdout + 2 * training,
                    )
                    row = dict(
                        model=kind,
                        fast=fast,
                        policy=policy,
                        sampler=fitted["sampler"],
                        request_schema=validation_request["schema"],
                        route=route,
                        seed=seed,
                        rates=np.exp(x),
                        truth=truth,
                        parameter_names=names,
                        training_loss=value,
                        truth_training_loss=truth_training_loss,
                        truth_training_work=dict(
                            events=at_truth["events"],
                            paths=at_truth["completed_replicates"],
                        ),
                        best_is_start=bool(np.allclose(x, x0, rtol=0, atol=1e-12)),
                        holdout_loss=loss(fitted, target),
                        truth_holdout_loss=loss(reference, target),
                        parameter_relative_error=float(
                            np.linalg.norm(np.exp(x) - truth) / np.linalg.norm(truth)
                        ),
                        independent_moment_max_standard_errors=error,
                        requested_evaluations=budget,
                        actual_evaluations=actual,
                        setup_seconds=setup_seconds,
                        search_seconds=search_seconds,
                        validation_seconds=time.perf_counter() - validation_start,
                        search_work=work,
                        validation_work=validation_work,
                    )
                    report["arms"].append(row)
                    arm_manifest(
                        args.output,
                        tag,
                        budget,
                        actual,
                        search_seconds,
                        seed,
                        args.threads,
                    )
                    comparison.append(trace)
                    write(args.output / "summary.json", report)
                same = comparison[0] == comparison[1]
                require(
                    same,
                    "fresh/reusable optimizer candidate streams or objective values differ",
                )
                report["checks"].append(
                    dict(
                        name=f"{kind}-{policy}-{seed}-fresh-bound-replay-and-moments",
                        passed=True,
                    )
                )
