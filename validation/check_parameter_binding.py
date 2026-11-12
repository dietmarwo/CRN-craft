"""Rules parameter binding and the benchmark worker before optimization timing.

The migrated counterpart of check_parameter_interface.py (which keeps testing
the network engine's ParameterTemplate until M8): the frozen fixtures are
converted to crn-craft.model/1.0, every candidate binds its values to a
prepared Static Direct run (PreparedEnsemble.bind_parameters, the worker's
`reuse` route) and must equal a fresh preparation of the model with those
values (the worker's `fresh` route), in Python and in the worker.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import tempfile

import crn_craft._legacy as crn
from crn_craft import rules
from parameter_optimization import Worker, save
from parameter_optimization_models import fixtures

PATHS = 4
SEED = 7301


def require(value, message):
    if not value:
        raise RuntimeError(message)


def request(model, times, wait, **budget):
    """The worker's training request for `model` (bytes)."""
    return {
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": hashlib.sha256(model).hexdigest()},
        "method": {"kind": "static_direct", "waiting_time": wait},
        "time": {"end": times[-1], "report_times": times},
        "ensemble": {"replicates": PATHS, "root_seed": SEED},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": PATHS,
            "keep_instances": 0,
        },
        "budget": {
            "threads": 1,
            "max_events_per_replicate": 1000000,
            "max_sampling_steps_per_replicate": 1000000000,
            "max_instances": 1000000,
            "max_output_bytes": 64000000,
        }
        | budget,
    }


def fresh(converted, values, times, wait, **budget):
    """A fresh preparation of the converted model with `values`."""
    model = json.loads(converted)
    model["parameters"] = {**model["parameters"], **values}
    data = crn.model_bytes(model)
    return rules.prepare(data, request(data, times, wait, **budget)).simulate()


def same_result(a, b, *, hashes=True):
    """Equal records, rule counts, statistics and kept paths (and hashes)."""
    keys = ["statistics", "events", "termination", "sampler", "paths"]
    if hashes:
        keys += ["model_sha256", "run_sha256"]
    return (
        a.records.equals(b.records)
        and a.rule_events.equals(b.rule_events)
        and all(a.metadata.get(k) == b.metadata.get(k) for k in keys)
    )


def untimed(reply):
    """A worker reply without its timings."""
    return {k: v for k, v in reply.items() if not k.endswith("_seconds")}


def refused(action, code):
    try:
        action()
    except (crn.CrnError, ValueError, TypeError) as error:
        diagnostics = getattr(error, "diagnostics", None)
        return code is None or (diagnostics and diagnostics[0]["code"] == code)
    return False


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worker", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = dict(
        schema="crn-craft.parameter-binding-validation/0.1",
        passed=False,
        bindings=0,
        worker_evaluations=0,
        invalid_inputs=0,
    )
    worker = Worker(args.worker)
    try:
        for invalid, reason in [
            (dict(command="eval", values=[1.0], unknown=True), "unknown field"),
            (dict(command="init", paths="4"), "unsigned integer"),
            # The network engine's options and the generated route are gone.
            (dict(command="init", algorithm="direct"), "unknown field"),
            (dict(command="init", fast=True), "unknown field"),
            (dict(command="validate_kernel"), "unknown command"),
        ]:
            try:
                worker.call(**invalid)
            except RuntimeError as error:
                require(reason in str(error), f"wrong refusal: {error}")
                report["invalid_inputs"] += 1
            else:
                raise RuntimeError(f"malformed worker request accepted: {invalid}")
    finally:
        worker.close()
    for fixture in fixtures():
        converted = crn.convert(fixture["model"])
        names = sorted(fixture["model"]["parameters"])
        truth = [fixture["model"]["parameters"][name] for name in names]
        times = fixture["times"]
        types = [t["name"] for t in json.loads(converted)["types"]]
        for wait in ["inverse", "ziggurat"]:
            template = rules.prepare(converted, request(converted, times, wait))
            baseline = template.simulate()
            workers = {}
            for engine in ["fresh", "reuse"]:
                workers[engine] = Worker(args.worker)
                workers[engine].call(
                    command="init",
                    engine=engine,
                    case=fixture["case"],
                    model=fixture["model"],
                    times=times,
                    paths=PATHS,
                    heldout_paths=PATHS,
                    seed=SEED,
                    keep=PATHS,
                    wait=wait,
                )
            try:
                for multiplier in [[1.0, 1.0], [0.0, 1.7], [1.3, 0.6]]:
                    values = [a * b for a, b in zip(truth, multiplier, strict=True)]
                    named = dict(zip(names, values, strict=True))
                    bound = template.bind_parameters(named).simulate()
                    # The bound model bytes are the engine's re-serialization,
                    # so only the hashes differ from the fresh preparation.
                    require(
                        same_result(
                            bound, fresh(converted, named, times, wait), hashes=False
                        ),
                        "binding differs from a fresh preparation",
                    )
                    documents = {
                        engine: w.call(command="eval", values=values, full=True)
                        for engine, w in workers.items()
                    }
                    require(
                        untimed(documents["fresh"]) == untimed(documents["reuse"]),
                        "worker routes differ",
                    )
                    document = documents["reuse"]["result"]
                    require(
                        all(
                            bound.metadata[k] == document[k]
                            for k in ["model_sha256", "run_sha256", "statistics"]
                        ),
                        "worker result differs from the Python binding",
                    )
                    mean = bound.statistics["mean"][types].to_numpy().tolist()
                    require(
                        documents["reuse"]["moments"]["mean"] == mean,
                        "worker moments are not in the model's type order",
                    )
                    report["bindings"] += 1
                    report["worker_evaluations"] += 2
            finally:
                for w in workers.values():
                    w.close()
            require(
                same_result(template.simulate(), baseline),
                "binding mutated the template",
            )
        save(args.output, report)
    f = fixtures()[0]
    converted = crn.convert(f["model"])
    times = f["times"]
    template = rules.prepare(converted, request(converted, times, "ziggurat"))
    for overrides in [
        {"unknown": 1.0},
        {"birth": math.nan},
        {"birth": math.inf},
        {"birth": "1"},
    ]:
        require(
            refused(lambda o=overrides: template.bind_parameters(o), None),
            f"invalid parameter values accepted: {overrides}",
        )
        report["invalid_inputs"] += 1
    worker = Worker(args.worker)
    try:
        worker.call(
            command="init",
            engine="reuse",
            case=f["case"],
            model=f["model"],
            times=times,
            paths=PATHS,
            heldout_paths=PATHS,
        )
        for values in [[1.0], [-1.0, 1.0], [1.0, 2.0, 3.0]]:
            try:
                worker.call(command="eval", values=values)
            except RuntimeError as error:
                require("invalid benchmark values" in str(error), str(error))
                report["invalid_inputs"] += 1
            else:
                raise RuntimeError(f"worker accepted {values}")
    finally:
        worker.close()

    def run(values):
        return template.bind_parameters(dict(birth=values[0], death=values[1]))

    values = [[10.0, 0.4], [20.0, 0.6], [30.0, 0.9], [0.0, 0.1]] * 3
    with ThreadPoolExecutor(max_workers=4) as pool:
        concurrent = list(pool.map(lambda v: run(v).simulate(), values))
    require(
        all(
            same_result(a, run(v).simulate())
            for a, v in zip(concurrent, values, strict=True)
        ),
        "concurrent binding mutated another candidate",
    )
    report["concurrent_bindings"] = len(values)
    with tempfile.TemporaryDirectory() as folder:
        cancel = Path(folder) / "cancel"
        cancel.write_text("cancel")
        stopped = run([20.0, 0.6]).simulate(cancel_file=cancel)
        require(
            stopped.metadata["termination"] == "cancelled"
            and stopped.statistics is None,
            "binding lost cancellation",
        )
    for budget, limit in [
        ("max_events_per_replicate", 1),
        ("max_sampling_steps_per_replicate", 1),
    ]:
        limited = rules.prepare(
            converted, request(converted, times, "ziggurat", **{budget: limit})
        )
        bound = limited.bind_parameters(dict(birth=20.0, death=0.6)).simulate()
        require(
            same_result(
                bound,
                fresh(
                    converted,
                    dict(birth=20.0, death=0.6),
                    times,
                    "ziggurat",
                    **{budget: limit},
                ),
                hashes=False,
            )
            and bound.metadata["termination"] != "completed",
            f"binding changed the {budget} stop",
        )
    many = rules.prepare(converted, request(converted, times, "ziggurat", threads=4))
    one = run([20.0, 0.6]).simulate()
    other = many.bind_parameters(dict(birth=20.0, death=0.6)).simulate()
    require(
        same_result(one, other, hashes=False)
        and one.metadata["run_sha256"] != other.metadata["run_sha256"],
        "parameter binding changed with the thread count",
    )
    # Why there is no generated route: a generated loop compiles the
    # parameter values in, so each candidate has its own loop source.
    source = template._native.loop_source()
    candidate = run([20.0, 0.6])._native.loop_source()
    require(
        source.startswith("//") and source != candidate,
        "generated loop source does not depend on the parameter values",
    )
    report["passed"] = True
    save(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
