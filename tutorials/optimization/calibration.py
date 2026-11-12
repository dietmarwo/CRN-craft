"""C1: complete fresh-versus-bound parameter fits with independent moments."""

import copy
import time
import numpy as np
import crn_craft as crn
import fcmaes_rust
from common import arm_manifest, checked, require, write


START, LOWER, UPPER = 1.25, 0.35, 2.5  # multiples of the true rates
TIMES = np.linspace(0, 6, 7).tolist()
METHODS = {"direct": "static_direct", "next_reaction": "next_reaction"}


def options(parser):
    parser.add_argument(
        "--include-next-reaction",
        action="store_true",
        help="add a separate arm with the next-reaction method (its own random stream)",
    )


def binding_model():
    """The reversible binding fixture as a crn-craft.model/1.0 document."""
    model = crn.Model("reversible_binding")
    model.add_species({"A": 30, "B": 0})
    model.add_parameter({"forward": 0.4, "reverse": 0.2})
    model.add_reaction("forward * A", ["A"], ["B"], name="bind")
    model.add_reaction("reverse * B", ["B"], ["A"], name="release")
    return model.to_dict()


def without_hashes(result):
    """Result metadata without the hashes of the serialized model and request.

    A binding serializes the bound model itself, so its hashes need not equal
    those of a model prepared afresh from the same values; every other field,
    including the sampler, events and statistics, must.
    """
    return {
        k: v
        for k, v in result.metadata.items()
        if k not in ("model_sha256", "run_sha256")
    }


def request(model, method, **settings):
    """One Static Direct (or next-reaction) ensemble of the model document."""
    return crn.ensemble(
        crn.Model.from_dict(model),
        end=6,
        report_times=TIMES,
        method=method,
        quantiles=(),
        max_events=100000,
        **settings,
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
    statistics = result.statistics
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
    times = TIMES
    policies = ["direct", "next_reaction"] if args.include_next_reaction else ["direct"]
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
        policies=policies,
        event_cap_per_path=100000,
        threads=args.threads,
    )
    for kind, model in [
        ("immigration", crn.load_example("immigration-death").to_dict()),
        ("binding", binding_model()),
    ]:
        names = sorted(model["parameters"])
        truth = np.array([model["parameters"][name] for name in names])
        truth_values = dict(zip(names, truth.tolist(), strict=True))
        target = moments(kind, truth, times)
        write(args.output / f"{kind}.model.json", model)
        for policy in policies:
            method = METHODS[policy]
            setup = time.perf_counter()
            template = request(
                model, method, replicates=training, seed=77201, threads=args.threads
            )
            setup_seconds = time.perf_counter() - setup
            # The training loss of the true rates on the same sample tells
            # overfitting apart from search failure.
            at_truth = checked(template.bind_parameters(truth_values).simulate())
            require(
                at_truth.metadata["method"]["kind"] == method
                and at_truth.metadata["schema"] == "crn-craft.result/1.0",
                f"{kind}-{policy}: result does not state the requested method",
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
                        method=method,
                        target=target,
                        work=work,
                        trace=trace,
                    ):
                        rates = dict(zip(names, np.exp(x).tolist(), strict=True))
                        if route == "reusable":
                            result = checked(template.bind_parameters(rates).simulate())
                        else:
                            candidate = copy.deepcopy(model)
                            candidate["parameters"] = rates
                            result = checked(
                                request(
                                    candidate,
                                    method,
                                    replicates=training,
                                    seed=77201,
                                    threads=args.threads,
                                ).simulate()
                            )
                        value = loss(result, target)
                        work["events"] += result.metadata["events"]
                        work["paths"] += result.metadata["completed_replicates"]
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
                    validation_template = request(
                        model,
                        method,
                        replicates=holdout,
                        seed=99173,
                        threads=args.threads,
                    )
                    fitted_rates = dict(zip(names, np.exp(x).tolist(), strict=True))
                    fitted = checked(
                        validation_template.bind_parameters(fitted_rates).simulate()
                    )
                    reference = checked(
                        validation_template.bind_parameters(truth_values).simulate()
                    )
                    exact_fitted = moments(kind, np.exp(x), times)[0]
                    empirical = np.asarray(fitted.statistics["mean"])
                    se = np.asarray(fitted.statistics["mean_standard_error"])
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
                    # The ordinary replay prepares the fitted model afresh; the
                    # reusable binding must give the same result.
                    ordinary = copy.deepcopy(model)
                    ordinary["parameters"] = fitted_rates
                    replay = checked(
                        request(
                            ordinary,
                            method,
                            replicates=holdout,
                            seed=99173,
                            threads=args.threads,
                        ).simulate()
                    )
                    repeat = checked(
                        validation_template.bind_parameters(fitted_rates).simulate()
                    )
                    require(
                        without_hashes(replay) == without_hashes(repeat)
                        and replay.statistics["mean"].equals(repeat.statistics["mean"])
                        and replay.statistics["sample_variance"].equals(
                            repeat.statistics["sample_variance"]
                        ),
                        "ordinary replay and reusable binding differ",
                    )
                    tag = f"{kind}-{policy}-{seed}-{route}"
                    write(args.output / f"{tag}.trace.json", trace)
                    validation_work = dict(
                        events=fitted.metadata["events"]
                        + reference.metadata["events"]
                        + replay.metadata["events"]
                        + repeat.metadata["events"],
                        paths=2 * holdout + 2 * training,
                    )
                    row = dict(
                        model=kind,
                        policy=policy,
                        sampler=fitted.metadata["sampler"],
                        request_schema=fitted.metadata["schema"],
                        route=route,
                        seed=seed,
                        rates=np.exp(x),
                        truth=truth,
                        parameter_names=names,
                        training_loss=value,
                        truth_training_loss=truth_training_loss,
                        truth_training_work=dict(
                            events=at_truth.metadata["events"],
                            paths=at_truth.metadata["completed_replicates"],
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
