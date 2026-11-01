"""Actual pinned fcmaes worker/optimization-run end-to-end acceptance campaign."""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def require(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--worker", type=Path, required=True)
    p.add_argument("--runner", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    worker = args.worker.resolve()
    runner = args.runner.resolve()
    report = {
        "schema": "crn-craft.jobs-validation/0.1",
        "rows": [],
        "jobs": [],
        "protocol_pins": {
            "fcmaes": "583f33a9a815ede0e1115c8ddd69029a874593f9",
            "worker_SDK": "602d31afacf8bd2ecb0b0e941d18bd8e2821c1e3",
            "runner": "cc6b03d594f73d1f3b38f088e7e3bd8663814d75",
        },
    }

    def save():
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

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
            report["jobs"].append(
                {
                    "key": s["request_key"],
                    "response": result,
                    "manifest": manifest,
                    "parameters": s["payload"]["problem"]["parameters"],
                }
            )
            save()
            return result

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
        for mode in ["base", "held_out", "event_limit"]:
            folder, s = spec(evaluate, mode)
            s["payload"]["points"] = [
                [1.0986122886681098, -0.2231435513142097],
                [0.5, -0.5],
            ]
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
        for key in ["8_retry", "8_surrogate", "64_retry", "64_surrogate"]:
            search = next(j for j in report["jobs"] if j["key"] == key)
            point = search["manifest"]["results"][0]["outcome"]["x"]
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
            s["payload"]["problem"]["parameters"]["seed_set"] = [20261210, 20261211]
            execute(folder, s)
        folder, s = spec(evaluate, "active_deadline")
        s["limits"]["wall_time_seconds"] = 0.01
        s["payload"]["points"] = [[1.0986122886681098, -0.2231435513142097]]
        s["payload"]["problem"]["parameters"]["replicates"] = 10000
        path = folder / "job.json"
        path.write_text(json.dumps(s))
        timed = invoke(*prefix, "run", "--spec", path)
        require(
            timed["state"] in ["timed_out", "cancelled"], "deadline did not stop work"
        )
        report["rows"].append(
            {"check": "active_deadline", "passed": True, "state": timed["state"]}
        )
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
        report["rows"].append(
            {"check": "existing_command_protocol", "passed": True, "response": cmd}
        )
        manifests = []
        for path in node.rglob("run.json"):
            try:
                manifest = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            manifests.append(manifest)
        report["manifests"] = manifests
        report["rows"].append(
            {
                "check": "worker_protocol_and_profiles",
                "passed": True,
                "runs": len(report["jobs"]),
            }
        )
        by_key = {j["key"]: j["manifest"] for j in report["jobs"]}
        require(
            all(
                r["outcome"]["objective"] is None and not r["outcome"]["feasible"]
                for r in by_key["event_limit"]["results"]
            ),
            "failed objective accepted",
        )
        for key in ["8_retry", "8_surrogate", "64_retry", "64_surrogate"]:
            require(
                by_key[key]["work"]["evaluations"] == 32
                and by_key[key]["work"]["validation_evaluations"] >= 1,
                "search/validation accounting",
            )
        require(
            by_key["threads1"]["results"][0]["outcome"]["objective"]
            == by_key["threads8"]["results"][0]["outcome"]["objective"],
            "thread objective identity",
        )
        report["rows"].append(
            {
                "check": "incomplete_rejected_exact_budget_rechecks_thread_identity",
                "passed": True,
            }
        )
        report["passed"] = True
        save()
        print(
            json.dumps(
                {
                    "passed": True,
                    "jobs": len(report["jobs"]),
                    "manifests": len(manifests),
                }
            )
        )


if __name__ == "__main__":
    main()
