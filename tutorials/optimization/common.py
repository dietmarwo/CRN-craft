"""Bounded tutorial execution and evidence, independent of job infrastructure."""

import argparse
from dataclasses import asdict, is_dataclass
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def plain(value):
    if is_dataclass(value):
        return plain(asdict(value))
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, tuple | list):
        return [plain(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plain(value), indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(result):
    require(
        result["termination"] == "completed", f"incomplete SSA: {result['termination']}"
    )
    require(
        result["completed_replicates"] == result["requested_replicates"],
        "missing paths",
    )
    require(result["statistics"] is not None, "missing ensemble statistics")
    return result


def call_json(command, timeout=30):
    result = subprocess.run(
        list(map(str, command)), capture_output=True, text=True, timeout=timeout
    )
    require(
        result.returncode == 0,
        f"native command failed: {result.stderr[:1500]} {result.stdout[:1500]}",
    )
    return json.loads(result.stdout)


def arm_manifest(
    out,
    name,
    requested,
    actual,
    seconds,
    seed,
    workers=1,
    formulation="scalar",
    **extra,
):
    # Existing fcmaes tutorial run.json v1 convention; all objectives minimized.
    write(
        Path(out) / f"{name}.run.json",
        dict(
            schema_version=1,
            tutorial="crn-" + name,
            formulation=formulation,
            command="CRN-craft tutorials/optimization/run.py",
            seed=seed,
            workers=workers,
            requested_evaluations=requested,
            actual_evaluations=actual,
            elapsed_seconds=seconds,
            status="completed",
            objectives=[],
            descriptors=[],
            artifacts={"summary": "summary.json"},
            **extra,
        ),
    )


def main(chapter, function, options=None):
    """Run one chapter; `options` adds the chapter's own command-line flags."""
    parser = argparse.ArgumentParser(description=function.__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--preset", choices=["smoke", "study"], default="smoke")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed-offset", type=int, default=0)
    if options is not None:
        options(parser)
    parser.add_argument(
        "--reference-binary",
        type=Path,
        default=ROOT.parent
        / "fcmaes-rust/tutorials/oscillator-topology-search/target/release/examples/crn_reference",
    )
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument(
        "--ode-binary",
        type=Path,
        default=ROOT / "integrations/diffsol/target/release/crn-diffsol",
    )
    args = parser.parse_args()
    require(1 <= args.threads <= 8, "threads must be in 1..8")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(
        schema_version="crn-craft.optimization-tutorial/0.1",
        chapter=chapter,
        preset=args.preset,
        state="running",
        passed=False,
        checks=[],
        arms=[],
    )
    start = time.perf_counter()
    try:
        report["versions"] = {
            name: importlib.metadata.version(name)
            for name in ["crn-craft", "fcmaes-rust", "numpy", "scipy", "sqp-craft"]
        }
        report["source_sha256"] = {
            p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))
        }
        report["setup_scope"] = (
            "Installed native libraries; build/import time is excluded and must be accounted separately in deployment studies."
        )
        function(args, report)
        require(all(c["passed"] for c in report["checks"]), "tutorial gate failed")
        report["passed"] = True
    except Exception as error:
        report["error"] = f"{type(error).__name__}: {error}"
    finally:
        report["state"] = "complete"
        report["total_seconds"] = time.perf_counter() - start
        write(args.output / "summary.json", report)
    print(
        json.dumps(
            {k: report[k] for k in ["chapter", "passed", "total_seconds"]}
            | {"error": report.get("error")}
        )
    )
    return 0 if report["passed"] else 1
