"""Build/test/time one declared iteration and retain its exact decision evidence."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "validation/julia-gap10"
BINARIES = Path("/tmp/crn-gap10-binaries")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("iteration")
    p.add_argument("--parent", type=Path, required=True)
    args = p.parse_args()
    label = args.iteration
    EVIDENCE.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = "/tmp/crn-gap10-target"
    prefix = EVIDENCE / label
    source = subprocess.check_output(
        ["git", "diff", "--name-only", "--", "crates", "bindings", "Cargo.toml"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    if not source:
        raise RuntimeError("no candidate source edits")
    parent_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    subprocess.run(["cargo", "fmt", "--all"], cwd=ROOT, check=True)
    patch = subprocess.check_output(
        ["git", "diff", "--binary", "--", *source], cwd=ROOT
    )
    prefix.with_suffix(".patch").write_bytes(patch)
    buildstart = time.monotonic()
    for suffix, command in [
        (
            "build",
            ["cargo", "build", "--release", "--workspace", "--examples", "--locked"],
        ),
        ("tests", ["cargo", "test", "--release", "--workspace", "--locked"]),
    ]:
        # A retry gets a new log; failed development runs remain reviewable.
        log = prefix.with_suffix("." + suffix + ".log")
        index = 1
        while log.exists():
            log = prefix.with_suffix(f".{suffix}-{index}.log")
            index += 1
        with log.open("w") as f:
            subprocess.run(
                command,
                cwd=ROOT,
                env=env,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
    build_seconds = time.monotonic() - buildstart
    binary = BINARIES / label
    shutil.copyfile("/tmp/crn-gap10-target/release/examples/bench_ssa", binary)
    binary.chmod(0o755)
    reportpath = prefix.with_suffix(".json")
    with prefix.with_suffix(".timing.log").open("w") as f:
        subprocess.run(
            [
                os.sys.executable,
                str(ROOT / "validation/run_gap10.py"),
                "--parent",
                str(args.parent),
                "--candidate",
                str(binary),
                "--iteration",
                label,
                "--output",
                str(reportpath),
            ],
            cwd=ROOT,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )
    report = json.loads(reportpath.read_text())
    report.update(
        parent_commit=parent_commit,
        build_and_test_seconds=build_seconds,
        parent_binary_sha256=hashlib.sha256(args.parent.read_bytes()).hexdigest(),
        candidate_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        patch_sha256=hashlib.sha256(patch).hexdigest(),
    )
    reportpath.write_text(json.dumps(report, indent=2) + "\n")
    if not report["accepted"]:
        subprocess.run(
            ["git", "restore", "--source=HEAD", "--", *source], cwd=ROOT, check=True
        )
    log = ROOT / "docs/campaign/julia-gap10-implementation-log-2026-11-02.md"
    with log.open("a") as f:
        f.write(
            f"\n## {label}\n\nDecision: **{'kept' if report['accepted'] else 'reverted'}**. Target geometric speedup {report['target_speedup']:.4f}×; all-row geometric speedup {report['guard_speedup']:.4f}×. Target aggregate rounds: "
            + ", ".join(f"{v:.4f}×" for v in report["target_rounds"])
            + f". All 35 numerical identities and the full release suite passed. Build/test time {build_seconds:.1f} s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/{label}.*`.\n\n"
        )
        f.write("| Row | Target | Paired median speedup |\n|---|---|---:|\n")
        for name, row in report["rows"].items():
            f.write(
                f"| {name} | {'yes' if row['target'] else 'guard'} | {row['speedup']:.4f}× |\n"
            )
    subprocess.run(
        [
            "git",
            "add",
            "--",
            *source,
            str(log.relative_to(ROOT)),
            "validation/julia-gap10",
        ],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [
            "git",
            "commit",
            "-m",
            f"{label}: "
            + (
                "keep tested Julia-gap optimization"
                if report["accepted"]
                else "retain rejected optimization evidence and restore parent"
            ),
        ],
        cwd=ROOT,
        check=True,
    )
    print(
        json.dumps(
            dict(
                iteration=label,
                accepted=report["accepted"],
                target_speedup=report["target_speedup"],
                guard_speedup=report["guard_speedup"],
                next_parent=str(binary if report["accepted"] else args.parent),
            ),
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
