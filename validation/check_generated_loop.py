"""Build generated rules loops, replay them against the CLI and time them.

Each replay case generates a loop for its model, compiles all loops into one
runner crate and compares the runner's output with `crn-craft rules` byte for
byte, with one and eight workers, under default Static Direct and the
Ziggurat and sparse-tracker options. The timing section runs the 128-cell
mtDNA probe with one worker under the interpreter, `mtdna_direct`, the
rates-only kernel and the generated loop, in alternating order. Timing is
host-specific.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter

from benchmark_generic_rules import FAST, TRACKERS, mt
from benchmark_mtdna import GRID, sparse
from check_mtdna import ROOT, read_model, require, rules_files

RUNNER = """
use crn_craft_rules::{Control, execute, execute_with_rules_kernel,
                      execute_with_rules_loop, prepare};
mod kernel { include!("kernel.rs"); }
@MODULES@
fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let model = std::fs::read(&args[2]).map_err(|e| e.to_string())?;
    let run = std::fs::read(&args[3]).map_err(|e| e.to_string())?;
    let prepared = prepare(&model, &run).map_err(|e| format!("{e:?}"))?;
    let control = Control::default();
    let result = match (args[0].as_str(), args[1].as_str()) {
        ("interpreter", _) => execute(&prepared, &control),
        ("kernel", _) => execute_with_rules_kernel(
            &prepared, &control, &kernel::GeneratedRulesKernel),
@ARMS@
        other => return Err(format!("unknown backend {other:?}")),
    }
    .map_err(|e| format!("{e:?}"))?;
    println!("{}", serde_json::to_value(&result).map_err(|e| e.to_string())?);
    Ok(())
}
"""


def cases():
    model = read_model("immigration-death.rules.json")
    stiff = read_model("immigration-death.rules.json")
    stiff["parameters"]["k_in"] = 1e-30
    stiff["initial"] = [{"type": "X", "count": 10}]
    aggregate = read_model("aggregation-constant.rules.json")
    aggregate["rules"][0]["products"][0]["attributes"]["size"] = 1
    overflow = read_model("immigration-death.rules.json")
    overflow["observables"]["cold"] = {"type": "X", "weight": "1e308"}
    undefined = read_model("immigration-death.rules.json")
    undefined["observables"]["cold"] = {"type": "X", "weight": "2"}
    undefined["stop"] = "if(X > 0, 0/0, 0)"
    small = read_model("mtdna-three-class-cme.rules.json")
    small["stop"] = "N > 0 && mutants >= 0.5 * N"
    # label, model, request settings, expected termination
    return [
        ("mt", mt(False), {}, "completed"),
        ("stop", mt(), {}, "completed"),
        ("trackers", mt(), {"trackers": TRACKERS}, "completed"),
        ("log", mt(), {"log": {"max_entries": 4}}, "completed"),
        ("sparse", sparse(), {}, "completed"),
        # Max and distinct projections are not generated.
        ("aggregate", aggregate, {}, "unsupported"),
        ("immigration", model, {}, "completed"),
        ("small-stop", small, {"trackers": TRACKERS}, "completed"),
        ("cold-overflow", overflow, {}, "failed"),
        ("stop-diagnostic", undefined, {}, "failed"),
        ("resolution", stiff, {}, "failed"),
        ("instance", sparse(), {"max_instances": 1}, "instance_limit"),
    ]


def run(args):
    start = perf_counter()
    result = subprocess.run(args, capture_output=True, check=False, timeout=600)
    return perf_counter() - start, result


def generate(generator, files, supported=True):
    _, result = run([generator, *files])
    if not supported:
        require(
            result.returncode != 0
            and b"generated rules loops support" in result.stderr,
            "unsupported model accepted",
        )
        return None
    require(result.returncode == 0, result.stderr.decode())
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path, help="crn-craft CLI")
    parser.add_argument("examples", type=Path, help="directory with the generators")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    loop_generator = args.examples / "generate_rules_loop"
    kernel_generator = args.examples / "generate_rules_kernel"
    report = {"replay": [], "timing": []}
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        source = folder / "src"
        source.mkdir()
        selected = cases()
        modules, arms = [], []
        for index, (label, model, _, expected) in enumerate(selected):
            case = folder / f"case-{index}"
            case.mkdir()
            files = rules_files(case, model, method="static_direct", n=1)
            code = generate(loop_generator, files, expected != "unsupported")
            if code is None:
                report["replay"].append({"case": label, "generated": False})
                continue
            (source / f"case_{index}.rs").write_bytes(code)
            modules.append(f'mod case_{index} {{ include!("case_{index}.rs"); }}')
            arms.append(
                f'        ("loop", "{label}") => execute_with_rules_loop('
                f"&prepared, &control, &case_{index}::GeneratedRulesLoop),"
            )
        probe = folder / "probe"
        probe.mkdir()
        files = rules_files(probe, mt(), method="static_direct", n=1)
        (source / "kernel.rs").write_bytes(generate(kernel_generator, files))
        main_rs = RUNNER.replace("@MODULES@", "\n".join(modules))
        (source / "main.rs").write_text(main_rs.replace("@ARMS@", "\n".join(arms)))
        crate = ROOT / "crates/crn-craft-rules"
        (folder / "Cargo.toml").write_text(
            '[package]\nname = "rules-loop-probe"\nversion = "0.1.0"\n'
            'edition = "2024"\n[workspace]\n'
            f'[dependencies]\ncrn-craft-rules = {{path = "{crate}"}}\n'
            'serde_json = "1"\n'
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
                folder / "Cargo.toml",
            ]
        )
        require(compiled.returncode == 0, compiled.stderr.decode())
        report["build_seconds"] = build
        exe = target / "release/rules-loop-probe"

        for index, (label, model, extra, expected) in enumerate(selected):
            case = folder / f"case-{index}"
            if expected == "unsupported":
                continue
            for method in ("static_direct", FAST):
                for workers in (1, 8):
                    settings = {
                        "method": method,
                        "n": 8,
                        "times": [0, 10, 100],
                        "threads": workers,
                    } | extra
                    files = rules_files(case, model, **settings)
                    _, old = run([args.binary, "rules", *files])
                    _, new = run([exe, "loop", label, *files])
                    # The CLI exits with 1 for failed ensembles; the runner
                    # always prints its result.
                    require(
                        new.returncode == 0 and old.stdout == new.stdout,
                        f"replay {label}/{method}/{workers}: {new.stderr.decode()}",
                    )
                    termination = json.loads(new.stdout)["termination"]
                    require(termination == expected, f"{label}: {termination}")
                    report["replay"].append(
                        {
                            "case": label,
                            "method": method,
                            "workers": workers,
                            "termination": termination,
                            "identical": True,
                        }
                    )
        other = mt()
        other["parameters"]["mu"] *= 2
        files = rules_files(probe, other, method="static_direct", n=1)
        _, invalid = run([exe, "loop", "stop", *files])
        require(
            invalid.returncode != 0 and b"fingerprint mismatch" in invalid.stderr,
            "mismatched loop accepted",
        )
        report["fingerprint_rejected"] = True

        def probe_files(method):
            return rules_files(
                probe,
                mt(),
                method=method,
                n=128,
                times=GRID,
                threads=1,
                seed=17,
                trackers=TRACKERS,
            )

        backends = [
            ("interpreter", lambda files: [args.binary, "rules", *files]),
            ("kernel", lambda files: [exe, "kernel", "-", *files]),
            ("loop", lambda files: [exe, "loop", "stop", *files]),
        ]
        medians = {}
        for request, method in (("generic", FAST), ("default", "static_direct")):
            payloads = {}
            for repeat in range(args.repeats):
                order = backends if repeat % 2 == 0 else backends[::-1]
                for label, command in order:
                    files = probe_files(method)
                    seconds, result = run(command(files))
                    require(result.returncode == 0, result.stderr.decode())
                    payloads.setdefault(label, result.stdout)
                    require(payloads[label] == result.stdout, "nondeterministic")
                    events = json.loads(result.stdout)["events"]
                    report["timing"].append(
                        {
                            "request": request,
                            "backend": label,
                            "repeat": repeat,
                            "seconds": seconds,
                            "events": events,
                        }
                    )
            require(len(set(payloads.values())) == 1, f"{request} outputs differ")
        for repeat in range(args.repeats):
            files = probe_files("mtdna_direct")
            seconds, result = run([args.binary, "rules", *files])
            require(result.returncode == 0, result.stderr.decode())
            report["timing"].append(
                {
                    "request": "special",
                    "backend": "mtdna_direct",
                    "repeat": repeat,
                    "seconds": seconds,
                    "events": json.loads(result.stdout)["events"],
                }
            )
        for row in report["timing"]:
            key = f"{row['request']}/{row['backend']}"
            medians.setdefault(key, []).append(row)
        summary = {}
        for key, rows in medians.items():
            seconds = sorted(r["seconds"] for r in rows)[len(rows) // 2]
            summary[key] = {
                "median_seconds": seconds,
                "events": rows[0]["events"],
                "ns_per_event": 1e9 * seconds / rows[0]["events"],
            }
            print(
                f"{key:28s} {seconds:.3f} s {summary[key]['ns_per_event']:.1f} ns/event"
            )
        report["summary"] = summary
        ns = {k: v["ns_per_event"] for k, v in summary.items()}
        report["loop_speedup_generic"] = ns["generic/interpreter"] / ns["generic/loop"]
        report["loop_speedup_default"] = ns["default/interpreter"] / ns["default/loop"]
        report["loop_gap_to_mtdna_direct"] = (
            ns["generic/loop"] / ns["special/mtdna_direct"]
        )
        print(json.dumps({k: report[k] for k in report if k.startswith("loop_")}))
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
