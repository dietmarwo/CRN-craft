"""Execute the final parameter-interface numerical/contract gates, retaining logs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from parameter_optimization_models import ROOT, FOLDER


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    args = parser.parse_args()
    target = args.target.resolve()
    out = FOLDER / "evidence/final"
    out.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = str(target)
    env["CARGO_BUILD_JOBS"] = "4"
    binary = str(target / "release/crn-craft")
    rows = []

    def run(name, command, *, cwd=ROOT, environment=env):
        with (out / f"{name}.log").open("w") as log:
            result = subprocess.run(
                command, cwd=cwd, env=environment, stdout=log, stderr=subprocess.STDOUT
            )
        rows.append(
            dict(
                check=name, passed=result.returncode == 0, returncode=result.returncode
            )
        )
        (out / "summary.json").write_text(
            json.dumps(
                dict(
                    complete=False, passed=all(r["passed"] for r in rows), checks=rows
                ),
                indent=2,
            )
            + "\n"
        )
        print(name, result.returncode, flush=True)
        if result.returncode:
            raise RuntimeError(name + " failed; evidence retained")

    def same(a, b):
        if a.read_bytes() != b.read_bytes():
            raise RuntimeError(f"Fixed-work reports differ: {a.name}, {b.name}")

    run("fmt", ["cargo", "fmt", "--all", "--check"])
    run(
        "binding-fmt",
        ["cargo", "fmt", "--manifest-path", "bindings/python/Cargo.toml", "--check"],
    )
    run(
        "worker-fmt",
        [
            "cargo",
            "fmt",
            "--manifest-path",
            "validation/parameter-optimization/worker/Cargo.toml",
            "--check",
        ],
    )
    run(
        "python-lint",
        [
            sys.executable,
            "-m",
            "ruff",
            "check",
            "validation",
            "python/crn_craft/__init__.py",
        ],
    )
    run(
        "python-fmt",
        [
            sys.executable,
            "-m",
            "ruff",
            "format",
            "--check",
            "validation",
            "python/crn_craft/__init__.py",
        ],
    )
    run(
        "build",
        ["cargo", "build", "--release", "--workspace", "--examples", "--locked"],
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
    )
    run(
        "worker-clippy",
        [
            "cargo",
            "clippy",
            "--manifest-path",
            "validation/parameter-optimization/worker/Cargo.toml",
            "--all-features",
            "--all-targets",
            "--locked",
            "--",
            "-D",
            "warnings",
        ],
    )
    run("release-tests", ["cargo", "test", "--release", "--workspace", "--locked"])
    msrv = dict(env, CARGO_TARGET_DIR=str(target.parent / "crn-parameters-msrv"))
    run(
        "msrv-tests",
        ["cargo", "+1.85.0", "test", "--release", "--workspace", "--locked"],
        environment=msrv,
    )
    wheel_dir = target / "final-wheel"
    run(
        "build-wheel",
        [
            sys.executable,
            "-m",
            "maturin",
            "build",
            "--release",
            "--locked",
            "--out",
            str(wheel_dir),
        ],
    )
    run(
        "install-wheel",
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--force-reinstall",
            "--no-deps",
            *map(str, sorted(wheel_dir.glob("*.whl"))),
        ],
    )
    for script in [
        "check_contracts",
        "check_alternative_contracts",
        "check_python",
        "check_ssa",
        "check_cme",
    ]:
        for optimized in [False, True]:
            name = script + ("-optimized" if optimized else "")
            run(
                name,
                [sys.executable]
                + (["-O"] if optimized else [])
                + [
                    f"validation/{script}.py",
                    "--binary",
                    binary,
                    "--output",
                    str(out / f"{name}.json"),
                ],
            )
        same(out / f"{script}.json", out / f"{script}-optimized.json")
    for fast in [False, True]:
        for optimized in [False, True]:
            name = ("fast-laws" if fast else "default-laws") + (
                "-optimized" if optimized else ""
            )
            run(
                name,
                [sys.executable]
                + (["-O"] if optimized else [])
                + ["validation/check_alternative_ssa.py"]
                + (["--fast"] if fast else [])
                + ["--output", str(out / f"{name}.json")],
            )
        name = "fast-laws" if fast else "default-laws"
        same(out / f"{name}.json", out / f"{name}-optimized.json")
    for script in ["check_ssa", "check_cme"]:
        name = script + "-cached"
        run(
            name,
            [
                sys.executable,
                f"validation/{script}.py",
                "--binary",
                binary,
                "--command",
                "simulate-cached",
                "--output",
                str(out / f"{name}.json"),
            ],
        )
        same(out / f"{script}.json", out / f"{name}.json")
    for optimized in [False, True]:
        name = "dependencies" + ("-optimized" if optimized else "")
        run(
            name,
            [sys.executable]
            + (["-O"] if optimized else [])
            + ["validation/check_dependencies.py"],
        )
    same(out / "dependencies.log", out / "dependencies-optimized.log")
    for optimized in [False, True]:
        name = "installed" + ("-optimized" if optimized else "")
        run(
            name,
            [sys.executable, "-I"]
            + (["-O"] if optimized else [])
            + [
                str(ROOT / "validation/check_installed.py"),
                "--output",
                str(out / f"{name}.json"),
            ],
            cwd=Path("/tmp"),
        )
    same(out / "installed.json", out / "installed-optimized.json")
    for optimized in [False, True]:
        name = "parameter-validation" + ("-optimized" if optimized else "")
        run(
            name,
            [sys.executable]
            + (["-O"] if optimized else [])
            + [
                "validation/check_parameter_interface.py",
                "--output",
                str(FOLDER / f"evidence/{name}.json"),
            ],
        )
    same(
        FOLDER / "evidence/parameter-validation.json",
        FOLDER / "evidence/parameter-validation-optimized.json",
    )
    (out / "summary.json").write_text(
        json.dumps(
            dict(complete=True, passed=all(r["passed"] for r in rows), checks=rows),
            indent=2,
        )
        + "\n"
    )
    print("All final gates passed", flush=True)


if __name__ == "__main__":
    main()
