"""Actual pinned fcmaes worker/optimization-run end-to-end acceptance campaign.

The main report holds only reproducible content: requested points, objectives,
work counts and asserted checks. Runner responses and native manifests, which
carry paths, timings and attempt identifiers, go to a separate raw evidence file.
"""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import crn_craft as crn


def require(value, message):
    if not value:
        raise RuntimeError(message)


def flatten(rows):
    return [v for row in rows for v in row]


def recompute(model, data, parameters, x):
    """Independent Python evaluation of the worker objective through the binding.

    Same engine (cached Direct SSA), seeds and replicate streams; the loss, seed
    mean and circuit amplitude penalty are recomputed here from saved statistics.
    """
    # The worker's documented default: every model parameter, in sorted order.
    names = parameters.get("parameter_names", sorted(model["parameters"]))
    replicates = parameters["replicates"]
    circuit = "amplitude_species" in parameters
    times = data["times"]
    loss = 0.0
    for seed in parameters["seed_set"]:
        candidate = copy.deepcopy(model)
        for name, log_rate in zip(names, x, strict=True):
            candidate["parameters"][name] = math.exp(log_rate)
        run = crn.run_spec(
            candidate,
            end=times[-1],
            save_times=times,
            replicates=replicates,
            seed=seed,
            keep=replicates if circuit else 0,
            max_events=parameters.get("max_events", 1000000),
            max_output_bytes=268435456,
        )
        result = crn.simulate(candidate, run, cached=True)
        require(result["termination"] == "completed", "reference run incomplete")
        stats = result["statistics"]
        for actual, target, variance, target_v in zip(
            flatten(stats["mean"]),
            flatten(data["mean"]),
            flatten(stats["sample_variance"]),
            flatten(data["variance"]),
            strict=True,
        ):
            loss += ((actual - target) / max(abs(target), 1.0)) ** 2 + (
                (variance - target_v) / max(abs(target_v), 1.0)
            ) ** 2
        if circuit:
            species = [s["name"] for s in model["species"]].index(
                parameters["amplitude_species"]
            )
            burn_in = parameters.get("burn_in", 0.0)
            keep = [i for i, t in enumerate(times) if t > burn_in]
            penalty = 0.0
            for path in result["trajectories"]:
                values = [float(path["counts"][i][species]) for i in keep]
                mean = sum(values) / len(values)
                variance = sum((v - mean) ** 2 for v in values) / len(values)
                penalty += 1.0 / (1.0 + variance)
            require(len(result["trajectories"]) == replicates, "missing paths")
            loss += penalty / replicates
    return loss / len(parameters["seed_set"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--worker", type=Path, required=True)
    p.add_argument("--runner", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--raw-output",
        type=Path,
        help="runner responses and native manifests (default: OUTPUT with .raw.json)",
    )
    args = p.parse_args()
    raw_output = args.raw_output or args.output.with_suffix(".raw.json")
    worker = args.worker.resolve()
    runner = args.runner.resolve()
    report = {
        "schema": "crn-craft.jobs-validation/0.2",
        "rows": [],
        "jobs": [],
        "protocol_pins": {
            "fcmaes": "583f33a9a815ede0e1115c8ddd69029a874593f9",
            "worker_SDK": "602d31afacf8bd2ecb0b0e941d18bd8e2821c1e3",
            "runner": "cc6b03d594f73d1f3b38f088e7e3bd8663814d75",
        },
    }
    raw = {"schema": "crn-craft.jobs-validation-raw/0.1", "jobs": [], "other": []}

    def save():
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        raw_output.write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n")

    def row(check, passed, **details):
        report["rows"].append({"check": check, "passed": bool(passed), **details})
        save()
        require(passed, f"{check} failed; retained report")

    def invoke(*argv, expected=0):
        proc = subprocess.run(
            list(map(str, argv)), capture_output=True, check=False, timeout=180
        )
        require(
            proc.returncode == expected,
            f"exit {proc.returncode}: {proc.stdout.decode()} {proc.stderr.decode()}",
        )
        return json.loads(proc.stdout)

    discovery = invoke(worker, "describe")
    with tempfile.TemporaryDirectory(prefix="crn-jobs-") as temp:
        base = Path(temp)
        node = base / "node"
        prefix = [runner, "--root", node]
        invoke(*prefix, "node", "init", "--capacity", "8", "--allow-commands")
        invoke(*prefix, "workers", "add", "--id", "crn", "--executable", worker)

        def spec(example, key):
            folder = base / key
            folder.mkdir()
            inputs = []
            for item in example.get("inputs", []):
                content = item.get("content_text")
                content = (
                    content.encode()
                    if content is not None
                    else (json.dumps(item["content"]) + "\n").encode()
                )
                (folder / item["file"]).write_bytes(content)
                inputs.append(
                    {
                        "name": item["name"],
                        "file": item["file"],
                        "sha256": hashlib.sha256(content).hexdigest(),
                    }
                )
            s = {
                "schema_version": "optimization.job/0.1",
                "request_key": key,
                "worker": "crn",
                "payload_schema": example["profile"],
                "resources": example["resources"],
                "limits": {"wall_time_seconds": 60.0, "grace_seconds": 1.0},
                "inputs": inputs,
                "payload": copy.deepcopy(example["payload"]),
            }
            return folder, s

        def execute(folder, s):
            path = folder / "job.json"
            path.write_text(json.dumps(s))
            invoke(*prefix, "validate", "--spec", path)
            result = invoke(*prefix, "run", "--spec", path)
            # Preserve the runner's own response rather than assuming a manifest path.
            candidates = [
                json.loads(path.read_text())
                for path in Path(result["directory"]).rglob("run.json")
            ]
            manifest = next(
                (
                    m
                    for m in candidates
                    if m.get("schema_version") == "optimization.result/0.2"
                ),
                None,
            )
            require(manifest is not None, "missing native result manifest")
            raw["jobs"].append(
                {"key": s["request_key"], "response": result, "manifest": manifest}
            )
            report["jobs"].append(
                {
                    "key": s["request_key"],
                    "state": result["state"],
                    "problem": s["payload"]["problem"]["id"],
                    "parameters": s["payload"]["problem"]["parameters"],
                    "points": s["payload"].get("points"),
                    "results": [
                        {
                            "x": r["outcome"]["x"],
                            "objective": r["outcome"]["objective"],
                            "feasible": r["outcome"]["feasible"],
                        }
                        for r in manifest["results"]
                    ],
                    "work": manifest["work"],
                }
            )
            save()
            return manifest

        examples = copy.deepcopy(discovery["examples"])
        fit_template = next(
            x
            for x in examples
            if x["payload"]["problem"]["id"] == "app.crn.fit-moments"
        )
        for item in examples:
            item["inputs"] = copy.deepcopy(fit_template["inputs"])
            item["payload"]["problem"] = copy.deepcopy(
                fit_template["payload"]["problem"]
            )
        evaluate = next(x for x in examples if x["profile"] == "fcmaes.evaluate/1")
        truth = [1.0986122886681098, -0.2231435513142097]  # log 3, log 0.8
        off = [0.5, -0.5]
        for mode in ["base", "base_repeat", "base_reversed", "held_out", "event_limit"]:
            folder, s = spec(evaluate, mode)
            s["payload"]["points"] = (
                [off, truth] if mode == "base_reversed" else [truth, off]
            )
            if mode == "held_out":
                s["payload"]["problem"]["parameters"]["seed_set"] = [20261201, 20261202]
            if mode == "event_limit":
                s["payload"]["problem"]["parameters"]["max_events"] = 1
            execute(folder, s)
        for replicas in [8, 64]:
            for profile in ["fcmaes.retry/1", "fcmaes.surrogate/1"]:
                if profile.endswith("retry/1"):
                    template = next(
                        x
                        for x in examples
                        if x["profile"] == profile
                        and x["payload"]["optimizer"]["kind"] == "biteopt"
                    )
                else:
                    template = next(x for x in examples if x["profile"] == profile)
                folder, s = spec(
                    template, f"{replicas}_{profile.split('.')[1].split('/')[0]}"
                )
                s["payload"]["problem"]["parameters"]["replicates"] = replicas
                if profile.endswith("retry/1"):
                    s["payload"]["retries"] = 1
                    s["payload"]["evaluations_per_retry"] = 32
                else:
                    s["payload"]["evaluations"] = 32
                execute(folder, s)
        # Each optimized point receives a separately counted disjoint-seed validation job.
        search_keys = ["8_retry", "8_surrogate", "64_retry", "64_surrogate"]
        for key in search_keys:
            search = next(j for j in report["jobs"] if j["key"] == key)
            point = search["results"][0]["x"]
            folder, s = spec(evaluate, "validation_" + key)
            s["payload"]["points"] = [point]
            s["payload"]["problem"]["parameters"].update(
                {"seed_set": [20261250, 20261251], "replicates": 256}
            )
            execute(folder, s)
        # Same fixed seeds and numerical objective with one and eight allowed threads.
        template = next(x for x in examples if x["profile"] == "fcmaes.surrogate/1")
        for threads in [1, 8]:
            folder, s = spec(template, f"threads{threads}")
            s["resources"] = {
                "cpu_slots": threads,
                "outer_workers": 1,
                "inner_threads": threads,
            }
            s["payload"]["evaluations"] = 8
            s["payload"]["evaluation_threads"] = threads
            execute(folder, s)
        # Two explicit topology fixtures; held-out validation is an independently budgeted job.
        circuit = next(
            x
            for x in discovery["examples"]
            if x["profile"] == "fcmaes.evaluate/1"
            and x["payload"]["problem"]["id"] == "app.crn.circuit-score"
        )
        circuit_inputs = {}
        for topology in ["ring", "broken_ring"]:
            folder, s = spec(circuit, topology)
            s["payload"]["points"] = [[5.375278407684165, 0.0]]
            if topology == "broken_ring":
                mp = folder / "model.json"
                model = json.loads(mp.read_text())
                for reaction in model["reactions"]:
                    if reaction["name"] == "tx3":
                        reaction["rate"]["vmax"] = 0.0
                content = json.dumps(model).encode()
                mp.write_bytes(content)
                s["inputs"][0]["sha256"] = hashlib.sha256(content).hexdigest()
            # m2 is not species 0, so the amplitude species must be resolved by name.
            s["payload"]["problem"]["parameters"].update(
                {"seed_set": [20261210, 20261211], "amplitude_species": "m2"}
            )
            circuit_inputs[topology] = (
                json.loads((folder / "model.json").read_text()),
                json.loads((folder / "data.json").read_text()),
            )
            execute(folder, s)

        # Contract: seed_set and amplitude_species have no defaults.
        rejected = {}
        for key, example, field in [
            ("no_seed_set", evaluate, "seed_set"),
            ("no_circuit_seed_set", circuit, "seed_set"),
            ("no_amplitude_species", circuit, "amplitude_species"),
        ]:
            folder, s = spec(example, key)
            s["payload"]["points"] = [list(example["payload"]["points"][0])]
            del s["payload"]["problem"]["parameters"][field]
            path = folder / "job.json"
            path.write_text(json.dumps(s))
            proc = subprocess.run(
                list(map(str, [*prefix, "validate", "--spec", path])),
                capture_output=True,
                check=False,
                timeout=180,
            )
            text = (proc.stdout + proc.stderr).decode()
            raw["other"].append({"key": key, "exit": proc.returncode, "output": text})
            rejected[key] = (
                proc.returncode != 0 and f"requires parameter {field}" in text
            )
        row(
            "required_parameters_rejected_when_missing",
            all(rejected.values()),
            jobs=rejected,
        )

        folder, s = spec(evaluate, "active_deadline")
        s["limits"]["wall_time_seconds"] = 0.01
        s["payload"]["points"] = [truth]
        s["payload"]["problem"]["parameters"]["replicates"] = 10000
        path = folder / "job.json"
        path.write_text(json.dumps(s))
        timed = invoke(*prefix, "run", "--spec", path)
        raw["other"].append({"key": "active_deadline", "response": timed})
        row("active_deadline", timed["state"] in ["timed_out", "cancelled"])
        # The same native CLI can run through the existing shell-free command protocol.
        command = {
            "schema_version": "execution.command/0.1",
            "request_key": "crn-command-discovery",
            "argv": [str(worker), "describe"],
            "env": {},
            "resources": {"cpu_slots": 1},
            "limits": {"wall_seconds": 10.0, "grace_seconds": 1.0},
            "replay": "never",
            "inputs": [],
            "outputs": [],
        }
        path = base / "command.json"
        path.write_text(json.dumps(command))
        invoke(*prefix, "validate", "--spec", path)
        cmd = invoke(*prefix, "run", "--spec", path)
        raw["other"].append({"key": "command", "response": cmd})
        row("existing_command_protocol", cmd.get("state") == "completed")
        manifests = []
        for path in node.rglob("run.json"):
            try:
                manifest = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            manifests.append(manifest)
        raw["manifests"] = manifests
        row(
            "worker_protocol_and_profiles",
            all(j["state"] == "completed" for j in report["jobs"]),
            runs=len(report["jobs"]),
            manifests=len(manifests),
        )

        by_key = {j["key"]: j for j in report["jobs"]}

        def objectives(key):
            return [r["objective"] for r in by_key[key]["results"]]

        row(
            "failed_objective_not_accepted",
            all(
                r["objective"] is None and not r["feasible"]
                for r in by_key["event_limit"]["results"]
            ),
        )
        for key in search_keys:
            work = by_key[key]["work"]
            row(
                f"search_and_validation_accounting_{key}",
                work["evaluations"] == 32 and work["validation_evaluations"] >= 1,
            )
        row(
            "thread_objective_identity",
            objectives("threads1") == objectives("threads8"),
            objectives=objectives("threads1"),
        )
        # Every evaluate job reports one result per requested point, in request order.
        evaluated = [j for j in report["jobs"] if j["points"] is not None]
        row(
            "results_follow_requested_point_order",
            all(
                [r["x"] for r in j["results"]] == j["points"]
                and j["work"]["evaluations"] == len(j["points"])
                for j in evaluated
            ),
            jobs=len(evaluated),
        )
        row(
            "repeat_and_permutation_identity",
            objectives("base") == objectives("base_repeat")
            and objectives("base_reversed") == objectives("base")[::-1],
            objectives=objectives("base"),
        )
        row(
            "true_rates_score_below_offset_rates_on_both_seed_sets",
            all(
                0.0 <= objectives(k)[0] < objectives(k)[1] for k in ["base", "held_out"]
            )
            and objectives("base") != objectives("held_out"),
            base=objectives("base"),
            held_out=objectives("held_out"),
        )
        row(
            "ring_scores_below_broken_ring",
            objectives("ring")[0] < objectives("broken_ring")[0],
            ring=objectives("ring")[0],
            broken_ring=objectives("broken_ring")[0],
        )
        # Values: every finite evaluate objective equals an independent recomputation.
        fit_inputs = {
            i["name"]: json.loads(i["content_text"]) for i in fit_template["inputs"]
        }
        comparisons = []
        for job in evaluated:
            if job["key"] == "event_limit":
                continue
            model, data = circuit_inputs.get(
                job["key"], (fit_inputs["model"], fit_inputs["data"])
            )
            for r in job["results"]:
                expected = recompute(model, data, job["parameters"], r["x"])
                comparisons.append(
                    {
                        "job": job["key"],
                        "x": r["x"],
                        "worker": r["objective"],
                        "recomputed": expected,
                        "relative_difference": abs(r["objective"] - expected)
                        / max(abs(expected), 1e-300),
                    }
                )
        row(
            "objective_values_match_independent_recomputation",
            all(c["relative_difference"] <= 1e-12 for c in comparisons),
            tolerance="relative 1e-12",
            comparisons=comparisons,
        )
        report["passed"] = all(r["passed"] for r in report["rows"])
        save()
        print(
            json.dumps(
                {
                    "passed": report["passed"],
                    "jobs": len(report["jobs"]),
                    "manifests": len(manifests),
                    "rows": len(report["rows"]),
                }
            )
        )


if __name__ == "__main__":
    main()
