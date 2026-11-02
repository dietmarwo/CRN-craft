"""Run final independent numerical/contract gates and keep every output."""

import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation/julia-gap10/final"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    python = sys.executable
    binary = "/tmp/crn-gap10-target/release/crn-craft"
    rows = []

    def run(name, command, cwd=ROOT, env=None):
        with (OUT / (name + ".log")).open("w") as log:
            p = subprocess.run(
                command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, env=env
            )
        rows.append(dict(check=name, passed=p.returncode == 0, returncode=p.returncode))
        (OUT / "summary.json").write_text(
            json.dumps(
                dict(passed=all(r["passed"] for r in rows), checks=rows), indent=2
            )
            + "\n"
        )
        print(name, p.returncode, flush=True)
        if p.returncode:
            raise RuntimeError(name + " failed; output retained")

    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = "/tmp/crn-gap10-target"
    run("fmt", ["cargo", "fmt", "--all", "--check"])
    run(
        "binding-fmt",
        [
            "cargo",
            "fmt",
            "--manifest-path",
            "bindings/python/Cargo.toml",
            "--all",
            "--check",
        ],
    )
    run(
        "clippy",
        [
            "cargo",
            "clippy",
            "--workspace",
            "--all-targets",
            "--locked",
            "--",
            "-D",
            "warnings",
        ],
        env=env,
    )
    run(
        "binding-clippy",
        [
            "cargo",
            "clippy",
            "--manifest-path",
            "bindings/python/Cargo.toml",
            "--all-targets",
            "--locked",
            "--",
            "-D",
            "warnings",
        ],
        env=env,
    )
    run(
        "release-tests",
        ["cargo", "test", "--release", "--workspace", "--locked"],
        env=env,
    )
    msrv = os.environ.copy()
    msrv["CARGO_TARGET_DIR"] = "/tmp/crn-gap10-msrv"
    msrv["CARGO_BUILD_JOBS"] = "4"
    run(
        "msrv-tests",
        ["cargo", "+1.85.0", "test", "--release", "--workspace", "--locked"],
        env=msrv,
    )
    for script in [
        "check_contracts",
        "check_alternative_contracts",
        "check_python",
        "check_alternative_ssa",
        "check_ssa",
        "check_cme",
    ]:
        for opt in [False, True]:
            name = script + ("-optimized" if opt else "")
            command = (
                [python]
                + (["-O"] if opt else [])
                + [f"validation/{script}.py", "--output", str(OUT / (name + ".json"))]
            )
            if script not in ["check_alternative_ssa"]:
                command += ["--binary", binary]
            run(name, command)
        if (OUT / (script + ".json")).read_bytes() != (
            OUT / (script + "-optimized.json")
        ).read_bytes():
            raise RuntimeError(script + " optimized report differs")
    for script in ["check_ssa", "check_cme"]:
        name = script + "-cached"
        run(
            name,
            [
                python,
                f"validation/{script}.py",
                "--binary",
                binary,
                "--command",
                "simulate-cached",
                "--output",
                str(OUT / (name + ".json")),
            ],
        )
        if (OUT / (name + ".json")).read_bytes() != (
            OUT / (script + ".json")
        ).read_bytes():
            raise RuntimeError("legacy cached identity differs")
    for opt in [False, True]:
        name = "dependencies" + ("-optimized" if opt else "")
        run(
            name,
            [python] + (["-O"] if opt else []) + ["validation/check_dependencies.py"],
        )
    if (OUT / "dependencies.log").read_bytes() != (
        OUT / "dependencies-optimized.log"
    ).read_bytes():
        raise RuntimeError("dependency reports differ")
    for opt in [False, True]:
        name = "installed" + ("-optimized" if opt else "")
        run(
            name,
            [python, "-I"]
            + (["-O"] if opt else [])
            + [
                str(ROOT / "validation/check_installed.py"),
                "--output",
                str(OUT / (name + ".json")),
            ],
            cwd=Path("/tmp"),
        )
    if (OUT / "installed.json").read_bytes() != (
        OUT / "installed-optimized.json"
    ).read_bytes():
        raise RuntimeError("installed reports differ")
    print("All final gates passed.", flush=True)


if __name__ == "__main__":
    main()
