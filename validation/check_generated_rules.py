"""Build and verify a generated mtDNA rules kernel against the interpreter."""

import argparse
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter

from benchmark_generic_rules import FAST, TRACKERS, mt
from benchmark_mtdna import GRID
from check_mtdna import ROOT, require, rules_files

RUNNER = """
use crn_craft::rules::{Control, Termination, execute_with_rules_kernel, prepare};
include!("kernel.rs");
fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let model = std::fs::read(&args[0]).map_err(|e| e.to_string())?;
    let run = std::fs::read(&args[1]).map_err(|e| e.to_string())?;
    let prepared = prepare(&model, &run).map_err(|e| format!("{e:?}"))?;
    let result = execute_with_rules_kernel(&prepared, &Control::default(), &GeneratedRulesKernel)
        .map_err(|e| format!("{e:?}"))?;
    println!("{}", serde_json::to_value(&result).map_err(|e| e.to_string())?);
    if result.termination != Termination::Completed { std::process::exit(1); }
    Ok(())
}
"""


def run(args):
    start = perf_counter()
    result = subprocess.run(args, capture_output=True, check=False, timeout=300)
    return perf_counter() - start, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("generator", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        source = folder / "src"
        source.mkdir()
        model = mt()
        files = rules_files(folder, model, n=1, method=FAST)
        _, generated = run([args.generator, *files])
        require(generated.returncode == 0, generated.stderr)
        (source / "kernel.rs").write_bytes(generated.stdout)
        (source / "main.rs").write_text(RUNNER)
        crate = ROOT / "crates/crn-craft"
        manifest = folder / "Cargo.toml"
        manifest.write_text(
            '[package]\nname = "rules-kernel-probe"\nversion = "0.1.0"\nedition = "2024"\n'
            f'[dependencies]\ncrn-craft = {{path = "{crate}"}}\nserde_json = "1"\n'
        )
        target = Path(os.environ.get("CARGO_TARGET_DIR", folder / "target")).resolve()
        build, compiled = run(
            [
                "cargo",
                "build",
                "--offline",
                "--release",
                "--target-dir",
                target,
                "--manifest-path",
                manifest,
            ]
        )
        require(compiled.returncode == 0, compiled.stderr.decode())
        exe = target / "release/rules-kernel-probe"
        for method in ("static_direct", FAST):
            for workers in (1, 8):
                files = rules_files(
                    folder,
                    model,
                    n=4,
                    times=GRID,
                    threads=workers,
                    method=method,
                    trackers=TRACKERS,
                )
                _, old = run([args.binary, "rules", *files])
                _, new = run([exe, *files])
                require(
                    old.returncode == new.returncode == 0 and old.stdout == new.stdout,
                    f"generated replay failed: {method}/{workers}: {new.stderr}",
                )
        other = mt()
        other["parameters"]["mu"] *= 2
        _, invalid = run([exe, *rules_files(folder, other, n=1, method=FAST)])
        require(
            invalid.returncode != 0 and b"fingerprint mismatch" in invalid.stderr,
            "mismatched kernel accepted",
        )

        for workers in (1, 8):
            files = rules_files(
                folder,
                model,
                n=512,
                times=GRID,
                threads=workers,
                seed=17,
                method=FAST,
                trackers=TRACKERS,
            )
            order = [
                ("generic", [args.binary, "rules", *files]),
                ("generated", [exe, *files]),
            ]
            for repeat in range(3):
                payloads = []
                for label, command in order[:: (-1 if repeat % 2 else 1)]:
                    seconds, result = run(command)
                    require(result.returncode == 0, result.stderr)
                    payloads.append(result.stdout)
                    row = {
                        "workers": workers,
                        "repeat": repeat,
                        "backend": label,
                        "seconds": seconds,
                    }
                    print(json.dumps(row), flush=True)
                    rows.append(row)
                require(payloads[0] == payloads[1], "timed outputs differ")
    args.output.write_text(
        json.dumps(
            {
                "build_seconds": build,
                "rows": rows,
                "replay": True,
                "fingerprint_rejected": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
