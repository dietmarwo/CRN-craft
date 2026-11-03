"""Native/Python alternative SSA contracts, controls and prepared reuse."""

import argparse
import concurrent.futures
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import crn_craft as crn
from jsonschema import Draft7Validator
from check_alternative_ssa import VARIANTS

ROOT = Path(__file__).resolve().parents[1]


def require(value, message):
    if not value:
        raise RuntimeError(message)


def equal_numbers(a, b):
    """Results equal apart from the request identity and the echoed `fast`."""

    def strip(r):
        out = {k: v for k, v in r.items() if k not in ("schema", "run_sha256")}
        out["method"] = {k: v for k, v in r["method"].items() if k != "fast"}
        return out

    return strip(a) == strip(b)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = []
    validators = {
        k: Draft7Validator(json.loads((ROOT / f"schemas/{k}.schema.json").read_text()))
        for k in ["ssa-run", "ssa-result"]
    }
    model = (
        ROOT / "validation/ssa-fixtures/immigration-nonzero.model.json"
    ).read_bytes()
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        mp = folder / "model.json"
        rp = folder / "request.json"
        mp.write_bytes(model)
        default = crn.ssa_spec(model, end=1.0, save_times=[1.0])
        require(
            default["schema"] == "crn-craft.ssa-run/0.3"
            and default["method"]["fast"] is True,
            "ssa_spec defaults to fast",
        )
        for invalid in [None, 1, "true"]:
            try:
                crn.ssa_spec(model, end=1.0, save_times=[1.0], fast=invalid)
            except TypeError:
                pass
            else:
                raise RuntimeError(f"ssa_spec accepted fast={invalid!r}")
        # Next reaction is checked only, so its default must be runnable.
        nr = crn.ssa_spec(
            model, aggregator="next_reaction", end=1.0, save_times=[0.0, 1.0]
        )
        require(
            nr["schema"] == "crn-craft.ssa-run/0.4" and nr["method"]["fast"] is False,
            "ssa_spec next_reaction defaults to the checked policy",
        )
        require(
            crn.simulate_ssa(model, nr)["termination"] == "completed",
            "default next_reaction request runs",
        )
        for (algorithm, wait), fast in itertools.product(VARIANTS, (False, True)):
            request = crn.ssa_spec(
                model,
                fast=fast,
                aggregator=algorithm,
                waiting_time=wait,
                end=2.0,
                save_times=[0.0, 0.25, 1.0, 2.0],
                replicates=16,
                keep=16,
                threads=2,
                seed=800,
                max_events=10000,
                max_sampling_steps=100000,
                max_output_bytes=1000000,
            )
            validators["ssa-run"].validate(request)
            expected = crn.simulate_ssa(model, request)
            validators["ssa-result"].validate(expected)
            require(
                request["schema"] == "crn-craft.ssa-run/0.3"
                and request["method"]["fast"] is fast
                and expected["schema"] == "crn-craft.ssa-result/0.3"
                and expected["sampler"].startswith(
                    "xoshiro256pp" if fast else "chacha8"
                )
                and expected["method"]["fast"] is fast,
                "fast identity",
            )

            # Without `fast`, ssa-run/0.3 is fast. The same policy under 0.2
            # (fast) or 0.1 (checked) gives the same numbers.
            implicit = json.loads(json.dumps(request))
            del implicit["method"]["fast"]
            implicit_result = crn.simulate_ssa(model, implicit)
            require(
                implicit_result["method"]["fast"] is True
                and equal_numbers(implicit_result, expected) == fast,
                "ssa-run/0.3 default is fast",
            )
            earlier = json.loads(json.dumps(request))
            if fast:
                earlier["schema"] = "crn-craft.ssa-run/0.2"
            else:
                earlier["schema"] = "crn-craft.ssa-run/0.1"
                del earlier["method"]["fast"]
            earlier_result = crn.simulate_ssa(model, earlier)
            validators["ssa-result"].validate(earlier_result)
            require(
                equal_numbers(earlier_result, expected)
                and earlier_result["schema"]
                == ("crn-craft.ssa-result/0.2" if fast else "crn-craft.ssa-result/0.1")
                and earlier_result["method"].get("fast", False) is fast,
                "earlier request versions",
            )
            for invalid in [None, 1, "true"]:
                changed = json.loads(json.dumps(request))
                changed["method"]["fast"] = invalid
                try:
                    crn.simulate_ssa(model, changed)
                except crn.CrnError:
                    pass
                else:
                    raise RuntimeError(f"fast={invalid!r} accepted")
            if fast:
                old = dict(request, schema="crn-craft.ssa-run/0.1")
                try:
                    crn.simulate_ssa(model, old)
                except crn.CrnError:
                    pass
                else:
                    raise RuntimeError("fast accepted under ssa-run/0.1")
            prepared = crn.prepare_ssa(model, request)
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                require(
                    all(
                        v == expected
                        for v in pool.map(
                            lambda _, prepared=prepared: prepared.simulate(), range(3)
                        )
                    ),
                    "concurrent replay",
                )
            rp.write_bytes(crn.model_bytes(request))
            for command in [
                [str(args.binary.resolve())],
                [sys.executable, "-m", "crn_craft"],
            ]:
                result = subprocess.run(
                    [*command, "ssa", str(mp), str(rp)],
                    capture_output=True,
                    check=False,
                    timeout=30,
                )
                require(result.returncode == 0 and not result.stderr, "CLI failure")
                require(json.loads(result.stdout) == expected, "CLI identity")
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
                    bad.write_bytes(crn.model_bytes(changed))
                    result = subprocess.run(
                        [*command, "ssa", str(mp), str(bad)],
                        capture_output=True,
                        check=False,
                        timeout=30,
                    )
                    require(
                        result.returncode == 1 and check(json.loads(result.stdout)),
                        f"CLI exit status: {command[-1]}",
                    )
            for field, value, stop in [
                ("max_sampling_steps_per_replicate", 0, "sampling_limit"),
                ("max_events_per_replicate", 1, "event_limit"),
                ("wall_time_s", 1e-12, "wall_time_limit"),
            ]:
                changed = json.loads(json.dumps(request))
                changed["budget"][field] = value
                incomplete = crn.simulate_ssa(model, changed)
                validators["ssa-result"].validate(incomplete)
                require(
                    incomplete["termination"] == stop
                    and incomplete["statistics"] is None,
                    stop,
                )
            cancel = folder / "CANCEL"
            cancel.touch()
            incomplete = prepared.simulate(cancel_file=cancel)
            cancel.unlink()
            require(
                incomplete["termination"] == "cancelled"
                and incomplete["work"]["rng_words"] == 0,
                "pre-cancel",
            )
            for changed in [
                dict(request, unknown=True),
                dict(request, method=dict(request["method"], unknown=True)),
            ]:
                try:
                    crn.simulate_ssa(model, changed)
                except crn.CrnError:
                    pass
                else:
                    raise RuntimeError("unknown field accepted")
            duplicate = crn.model_bytes(request).replace(
                b"{", b'{"schema":"crn-craft.ssa-run/0.1",', 1
            )
            try:
                crn.simulate_ssa(model, duplicate)
            except crn.CrnError:
                pass
            else:
                raise RuntimeError("duplicate key accepted")
            rows.append(
                dict(
                    algorithm=algorithm,
                    waiting_time=wait,
                    fast=fast,
                    passed=True,
                    result_sha256=hashlib.sha256(
                        json.dumps(expected, sort_keys=True).encode()
                    ).hexdigest(),
                )
            )
    report = dict(passed=True, checks=rows)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
