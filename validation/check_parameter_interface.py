"""The network engine's ParameterTemplate: binding and replay checks.

This keeps testing the 0.6.0 Python path (`prepare_parameters`,
`ParameterTemplate`, `simulate_ssa`, `mass_action_kernel_source`) until M6a
and M8 remove it. The benchmark worker runs the rules engine since M5g; its
checks moved to check_parameter_binding.py.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import itertools
import json
from pathlib import Path
import tempfile
import crn_craft._legacy as crn
from parameter_optimization_models import FOLDER, fixtures


def save(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def require(value, message):
    if not value:
        raise RuntimeError(message)


def spec(model, times, algorithm, wait, fast, **kwargs):
    return crn.ssa_spec(
        model,
        aggregator=algorithm,
        waiting_time=wait,
        fast=fast,
        end=times[-1],
        save_times=times,
        replicates=4,
        keep=4,
        seed=7301,
        max_events=1000000,
        max_output_bytes=64000000,
        **kwargs,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--output", type=Path, default=FOLDER / "evidence/parameter-validation.json"
    )
    args = p.parse_args()
    report = dict(
        passed=False,
        bindings=0,
        invalid_inputs=0,
    )
    for fixture in fixtures():
        model = crn.model_bytes(fixture["model"])
        truth = [
            fixture["model"]["parameters"][name]
            for name in sorted(fixture["model"]["parameters"])
        ]
        for algorithm, wait, fast in itertools.product(
            ["direct", "cached_direct", "direct_cr", "sorting_direct", "rssa_cr"],
            ["inverse", "ziggurat"],
            [False, True],
        ):
            request = spec(model, fixture["times"], algorithm, wait, fast)
            template = crn.prepare_parameters(model, request)
            for multiplier in [[1.0, 1.0], [0.0, 1.7], [1.3, 0.6]]:
                values = [a * b for a, b in zip(truth, multiplier, strict=True)]
                bound = template.bind(values)
                fresh = crn.simulate_ssa(bound.model_bytes, bound.request_bytes)
                actual = bound.simulate()
                require(actual == fresh, "Binding differs from fresh preparation")
                require(
                    template.simulate(
                        dict(zip(template.parameter_names, values, strict=True))
                    )
                    == actual,
                    "Named binding differs",
                )
                require(template.default_values == tuple(truth), "Template mutated")
                report["bindings"] += 1
                save(args.output, report)
    f = fixtures()[0]
    model = crn.model_bytes(f["model"])
    template = crn.prepare_parameters(
        model, spec(model, f["times"], "rssa_cr", "ziggurat", True)
    )
    for values in [
        [],
        [1.0],
        [float("nan"), 1.0],
        [1.0, float("inf")],
        [-1.0, 1.0],
        [True, 1.0],
        ["1", 1.0],
        {"birth": 1.0},
        {"birth": 1.0, "death": 1.0, "unknown": 1.0},
    ]:
        try:
            template.bind(values)
        except (crn.CrnError, ValueError, TypeError):
            report["invalid_inputs"] += 1
        else:
            raise RuntimeError("Invalid parameter input accepted")
    with ThreadPoolExecutor(max_workers=4) as pool:
        values = [[10.0, 0.4], [20.0, 0.6], [30.0, 0.9], [0.0, 0.1]] * 3
        concurrent = list(pool.map(template.simulate, values))
    require(
        concurrent == [template.simulate(v) for v in values],
        "Concurrent binding mutated another candidate",
    )
    report["concurrent_bindings"] = len(values)
    with tempfile.TemporaryDirectory() as folder:
        cancel = Path(folder) / "cancel"
        cancel.write_text("cancel")
        stopped = template.simulate([20.0, 0.6], cancel_file=cancel)
        require(
            stopped["termination"] == "cancelled" and stopped["statistics"] is None,
            "Binding lost cancellation",
        )
    for budget, value in [
        ("max_events_per_replicate", 1),
        ("max_sampling_steps_per_replicate", 1),
    ]:
        request = spec(model, f["times"], "direct_cr", "ziggurat", True)
        request["budget"][budget] = value
        t = crn.prepare_parameters(model, request)
        bound = t.bind([20.0, 0.6])
        require(
            bound.simulate()
            == crn.simulate_ssa(bound.model_bytes, bound.request_bytes),
            "Binding changed budget stop",
        )
    one = template.bind([20.0, 0.6], seed=774).simulate()
    many_request = spec(model, f["times"], "rssa_cr", "ziggurat", True, threads=4)
    many = (
        crn.prepare_parameters(model, many_request)
        .bind([20.0, 0.6], seed=774)
        .simulate()
    )
    one.pop("run_sha256")
    many.pop("run_sha256")
    require(one == many, "Parameter binding changed worker replay")
    require(
        crn.mass_action_kernel_source(model).startswith("// Generated"),
        "Generator not exposed",
    )
    custom = (
        Path(__file__).resolve().parents[1] / "models/examples/repressilator.model.json"
    ).read_bytes()
    try:
        crn.mass_action_kernel_source(custom)
    except crn.CrnError:
        report["invalid_inputs"] += 1
    else:
        raise RuntimeError("Generator accepted unsupported custom laws")
    report["passed"] = True
    save(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
