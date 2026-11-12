"""Native/Python contracts, controls and prepared reuse of every exact selector.

The 0.6.0 checks of the `ssa-run` contracts on the unified engine: for every
Static Direct selector and waiting time (and for `next_reaction`) a
`crn-craft.run/1.0` request built by `crn_craft.ensemble` validates against
the `run-1.0` schema, its result against `result-1.0`, the Python result
equals the result `crn-craft run` and `python -m crn_craft run` print, a
prepared run replays identically from three threads, every budget ends a run
with its termination and no statistics, a cancel file set before the start
ends it with no random words, and unknown fields, `null` options and duplicate
keys are refused with diagnostics.

Not carried over from 0.6.0 (docs/campaign/merge-0.7/M8a-log.md): the `fast`
policy and the `ssa-run` request versions 0.1 to 0.4 (xoshiro256++ and the
request versions are removed), and the cached-direct variant. Replay of the
rules request versions is checked by check_rules_python.py and the stored
requests of the Rust tests.
"""

import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import crn_craft as crn
from crn_craft.rules import model_bytes
from jsonschema import Draft7Validator
from check_alternative_ssa import VARIANTS

ROOT = Path(__file__).resolve().parents[1]
# Every variant of the law gate, plus the exact next-reaction method.
CASES = [*VARIANTS, ("next_reaction", None)]


def require(value, message):
    if not value:
        raise RuntimeError(message)


def method_for(algorithm, wait):
    if algorithm == "next_reaction":
        return "next_reaction", None
    return {"kind": "static_direct", "waiting_time": wait}, algorithm


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    binary = str(args.binary.resolve())
    rows = []
    validators = {
        k: Draft7Validator(crn.describe(k)) for k in ["run-1.0", "result-1.0"]
    }
    network = (
        ROOT / "validation/ssa-fixtures/immigration-nonzero.model.json"
    ).read_bytes()
    model = crn.convert(network)
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        mp = folder / "model.rules.json"
        rp = folder / "request.json"
        mp.write_bytes(model)

        def options(**more):
            base = {
                "end": 2.0,
                "report_times": [0.0, 0.25, 1.0, 2.0],
                "replicates": 16,
                "keep_paths": 16,
                "threads": 2,
                "seed": 800,
                "max_events": 10000,
                "max_sampling_steps": 100000,
                "max_output_bytes": 1000000,
                "quantiles": (),
            }
            base.update(more)
            return base

        # A bad selector or wait is refused before anything runs.
        for bad in ["fast", "", "RSSA_CR", 7]:
            try:
                crn.ensemble(model, **options(selector=bad, method="static_direct"))
            except ValueError:
                pass
            else:
                raise RuntimeError(f"ensemble accepted selector={bad!r}")
        try:
            crn.ensemble(
                model,
                **options(
                    method={"kind": "static_direct", "waiting_time": "fast"},
                    selector="linear",
                ),
            )
        except crn.CrnError:
            pass
        else:
            raise RuntimeError("ensemble accepted waiting_time='fast'")
        # Absent options mean the defaults of run/1.0: auto and Ziggurat waits,
        # which for two channels is the linear scan.
        default = crn.ensemble(model, **options(method="static_direct")).simulate()
        explicit = crn.ensemble(
            model,
            **options(
                method={"kind": "static_direct", "waiting_time": "ziggurat"},
                selector="linear",
            ),
        ).simulate()
        require(
            default.metadata["sampler"] == "chacha8-u53-ziggurat/rules-static-direct/1"
            and default.paths.equals(explicit.paths),
            "run/1.0 default is the linear scan with Ziggurat waits",
        )
        for algorithm, wait in CASES:
            method, selector = method_for(algorithm, wait)
            prepared = crn.ensemble(model, **options(method=method, selector=selector))
            expected = prepared.simulate()
            request = json.loads(json.dumps(expected.request))
            request["model"]["path"] = mp.name
            validators["run-1.0"].validate(request)
            validators["result-1.0"].validate(_full_result(binary, mp, request, folder))
            require(
                request["schema"] == "crn-craft.run/1.0"
                and expected.metadata["schema"] == "crn-craft.result/1.0"
                and expected.metadata["sampler"].startswith("chacha8-u53")
                and expected.metadata["termination"] == "completed",
                "contract identity",
            )
            if algorithm != "next_reaction":
                require(
                    expected.metadata["method"]["selector"] == algorithm
                    and expected.metadata["method"]["waiting_time"] == wait,
                    "method identity",
                )
            rp.write_bytes(model_bytes(request))
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                concurrent_results = list(
                    pool.map(lambda _, run=prepared: run.simulate(), range(3))
                )
            require(
                all(
                    v.metadata == expected.metadata and v.paths.equals(expected.paths)
                    for v in concurrent_results
                ),
                "concurrent replay",
            )
            for command in [[binary], [sys.executable, "-m", "crn_craft"]]:
                result = subprocess.run(
                    [*command, "run", str(mp), str(rp)],
                    capture_output=True,
                    check=False,
                    timeout=30,
                )
                require(result.returncode == 0 and not result.stderr, "CLI failure")
                native = json.loads(result.stdout)
                require(
                    {k: v for k, v in native.items() if k != "replicates"}
                    == {k: v for k, v in expected.metadata.items() if k != "replicates"}
                    and [r["events"] for r in native["replicates"]]
                    == expected.records["events"].tolist(),
                    "CLI identity",
                )
                # Incomplete runs and rejected requests exit 1 with JSON on stdout.
                limited = json.loads(json.dumps(request))
                limited["budget"]["max_events_per_replicate"] = 1
                unknown = dict(request, unknown=True)
                for changed, check in [
                    (limited, lambda r: r["termination"] == "event_limit"),
                    (
                        unknown,
                        lambda r: (
                            r["success"] is False
                            and r["diagnostics"]
                            and all("code" in d for d in r["diagnostics"])
                        ),
                    ),
                ]:
                    bad = folder / "changed.json"
                    bad.write_bytes(model_bytes(changed))
                    result = subprocess.run(
                        [*command, "run", str(mp), str(bad)],
                        capture_output=True,
                        check=False,
                        timeout=30,
                    )
                    require(
                        result.returncode == 1 and check(json.loads(result.stdout)),
                        f"CLI exit status: {command[-1]}",
                    )
            for field, value, stop in [
                ("max_sampling_steps_per_replicate", 1, "sampling_limit"),
                ("max_events_per_replicate", 1, "event_limit"),
                ("wall_time_s", 1e-12, "wall_time_limit"),
            ]:
                changed = json.loads(json.dumps(request))
                changed["budget"][field] = value
                incomplete = crn.prepare(model, changed, quantiles=()).simulate()
                require(
                    incomplete.metadata["termination"] == stop
                    and incomplete.statistics is None,
                    stop,
                )
            cancel = folder / "CANCEL"
            cancel.touch()
            incomplete = prepared.simulate(cancel_file=cancel)
            cancel.unlink()
            require(
                incomplete.metadata["termination"] == "cancelled"
                and int(incomplete.records["rng_words"].sum()) == 0,
                "pre-cancel",
            )
            changes = [dict(request, unknown=True)]
            if algorithm != "next_reaction":
                changes += [
                    dict(request, method=dict(request["method"], unknown=True)),
                    dict(request, method=dict(request["method"], selector=None)),
                ]
            for changed in changes:
                try:
                    crn.prepare(model, changed, quantiles=())
                except crn.CrnError:
                    pass
                else:
                    raise RuntimeError("unknown field or null option accepted")
            duplicate = model_bytes(request).replace(
                b"{", b'{"schema":"crn-craft.run/1.0",', 1
            )
            try:
                crn.prepare(model, duplicate, quantiles=())
            except crn.CrnError:
                pass
            else:
                raise RuntimeError("duplicate key accepted")
            rows.append(
                dict(
                    algorithm=algorithm,
                    waiting_time=wait,
                    passed=True,
                    result_sha256=hashlib.sha256(
                        json.dumps(expected.metadata, sort_keys=True).encode()
                    ).hexdigest(),
                )
            )
    report = dict(passed=True, checks=rows)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def _full_result(binary, model_path, request, folder):
    """The CLI's result document of `request` (for the schema check)."""
    path = folder / "schema-request.json"
    path.write_bytes(model_bytes(request))
    done = subprocess.run(
        [binary, "run", str(model_path), str(path)],
        capture_output=True,
        check=False,
        timeout=30,
    )
    require(done.returncode == 0 and not done.stderr, "CLI failure")
    return json.loads(done.stdout)


if __name__ == "__main__":
    main()
