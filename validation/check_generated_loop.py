"""Build generated rules loops, replay them against the CLI and time them.

Each replay case generates a loop for its model, compiles all loops into
one runner crate and compares the runner's output with `crn-craft rules`
byte for byte, with one and eight workers, under default Static Direct and
the Ziggurat and sparse-tracker options. A dense network whose loop would
exceed 1 MiB must be refused. The timing section runs, with one worker and
in alternating order:

- the 128-cell mtDNA probe under the interpreter, the rates-only kernel,
  the generated loop and `mtdna_direct`;
- the same probe with a stop predicate outside the fused fraction form;
- the 96-channel sparse and 512-channel wide networks.

Timing is host-specific.
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
import probe_structural_rules as probe_models

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

LAYOUTS = ("unrolled",)


def dense():
    """128 instances sharing one vacancy base: 256 densely coupled channels."""
    return {
        "schema": "crn-craft.rules-model/0.1",
        "name": "dense_128",
        "parameters": {"K": 2000, "r": 1.0, "d": 0.5},
        "types": [{"name": "A", "attributes": ["s"]}],
        "initial": [
            {"type": "A", "attributes": {"s": i}, "count": 8} for i in range(128)
        ],
        "observables": {
            "N": {"type": "A"},
            "D": {"type": "A", "weight": "1 + x.s / 64"},
        },
        "rules": [
            {
                "name": "rep",
                "reactants": [{"label": "x", "type": "A"}],
                "rate": "max(0, K - N) / D * r",
                "products": [{"type": "A", "copy": "x"}, {"type": "A", "copy": "x"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "A"}],
                "rate": "d",
                "products": [],
            },
        ],
    }


def general_stop():
    """The mtDNA probe with a stop the fused fraction form does not cover."""
    model = mt()
    model["stop"] = "N > 0 && mutants / N >= 0.9"
    return model


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
        # A closed aggregation model with a max projection.
        ("aggregate", aggregate, {}, "completed"),
        ("immigration", model, {}, "completed"),
        ("small-stop", small, {"trackers": TRACKERS}, "completed"),
        ("cold-overflow", overflow, {}, "failed"),
        ("stop-diagnostic", undefined, {}, "failed"),
        ("resolution", stiff, {}, "failed"),
        ("instance", sparse(), {"max_instances": 1}, "instance_limit"),
        ("general-stop", general_stop(), {"trackers": TRACKERS}, "completed"),
        # Refused: the loop would exceed 1 MiB of source.
        ("dense", dense(), {}, "completed"),
    ]


def run(args):
    start = perf_counter()
    result = subprocess.run(args, capture_output=True, check=False, timeout=600)
    return perf_counter() - start, result


def generate(generator, files, supported=True, layout=None):
    _, result = run([generator, *files, *([layout] if layout else [])])
    if not supported:
        require(
            result.returncode != 0
            and b"generated rules loops support" in result.stderr,
            "unsupported model accepted",
        )
        return None
    require(result.returncode == 0, result.stderr.decode())
    return result.stdout


def build(folder, loops, kernel_source):
    """Compile a runner with the rates kernel and every (layout, label) loop."""
    source = folder / "src"
    source.mkdir(parents=True)
    modules, arms = [], []
    for index, ((layout, label), code) in enumerate(loops.items()):
        (source / f"loop_{index}.rs").write_bytes(code)
        modules.append(f'mod loop_{index} {{ include!("loop_{index}.rs"); }}')
        arms.append(
            f'        ("{layout}", "{label}") => execute_with_rules_loop('
            f"&prepared, &control, &loop_{index}::GeneratedRulesLoop),"
        )
    (source / "kernel.rs").write_bytes(kernel_source)
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
    seconds, compiled = run(
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
    return seconds, target / "release/rules-loop-probe"


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
        loops = {}
        selected = cases()
        for index, (label, model, _, expected) in enumerate(selected):
            case = folder / f"case-{index}"
            case.mkdir()
            files = rules_files(case, model, method="static_direct", n=1)
            if expected == "unsupported":
                generate(loop_generator, files, False)
                report["replay"].append({"case": label, "generated": False})
                continue
            for layout in LAYOUTS:
                _, result = run([loop_generator, *files])
                if result.returncode == 0:
                    loops[(layout, label)] = result.stdout
                else:
                    # Only oversized loops may be refused.
                    require(
                        label == "dense" and b"1 MiB" in result.stderr,
                        result.stderr.decode(),
                    )
                    report["replay"].append(
                        {"case": label, "layout": layout, "generated": False}
                    )
        timing_models = {
            "mtdna": mt(),
            "general-stop": general_stop(),
            "sparse": probe_models.sparse(),
            "wide": probe_models.wide(),
        }
        for label, model in timing_models.items():
            work = folder / f"timing-{label}"
            work.mkdir()
            files = rules_files(work, model, method="static_direct", n=1)
            for layout in LAYOUTS:
                loops[(layout, f"timing-{label}")] = generate(loop_generator, files)
        probe = folder / "probe"
        probe.mkdir()
        files = rules_files(probe, mt(), method="static_direct", n=1)
        build_seconds, exe = build(
            folder / "runner", loops, generate(kernel_generator, files)
        )
        report["build_seconds"] = build_seconds

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
                    termination = json.loads(old.stdout)["termination"]
                    require(termination == expected, f"{label}: {termination}")
                    for layout in LAYOUTS:
                        if (layout, label) not in loops:
                            continue
                        _, new = run([exe, layout, label, *files])
                        # The CLI exits with 1 for failed ensembles; the
                        # runner always prints its result.
                        require(
                            new.returncode == 0 and old.stdout == new.stdout,
                            f"replay {label}/{layout}/{method}/{workers}: "
                            f"{new.stderr.decode()}",
                        )
                        report["replay"].append(
                            {
                                "case": label,
                                "layout": layout,
                                "method": method,
                                "workers": workers,
                                "termination": termination,
                                "identical": True,
                            }
                        )
        other = mt()
        other["parameters"]["mu"] *= 2
        files = rules_files(probe, other, method="static_direct", n=1)
        _, invalid = run([exe, "unrolled", "stop", *files])
        require(
            invalid.returncode != 0 and b"fingerprint mismatch" in invalid.stderr,
            "mismatched loop accepted",
        )
        report["fingerprint_rejected"] = True

        # workload: (model label, method, cells, report times, trackers)
        workloads = {
            "generic": ("mtdna", FAST, 128, GRID, TRACKERS),
            "default": ("mtdna", "static_direct", 128, GRID, TRACKERS),
            "general-stop": ("general-stop", FAST, 128, GRID, TRACKERS),
            "sparse": ("sparse", "static_direct", 64, [0, 50, 100], ()),
            "wide": ("wide", "static_direct", 8, [0, 50, 100], ()),
        }
        for name, (label, method, cells, times, trackers) in workloads.items():
            work = folder / f"timing-{label}"
            files = rules_files(
                work,
                timing_models[label],
                method=method,
                n=cells,
                times=times,
                threads=1,
                seed=17,
                trackers=trackers,
            )
            backends = [("interpreter", [args.binary, "rules", *files])]
            if label == "mtdna":
                backends.append(("kernel", [exe, "kernel", "-", *files]))
            for layout in LAYOUTS:
                if (layout, f"timing-{label}") in loops:
                    backends.append((layout, [exe, layout, f"timing-{label}", *files]))
            payloads = {}
            for repeat in range(args.repeats):
                order = backends if repeat % 2 == 0 else backends[::-1]
                for backend, command in order:
                    seconds, result = run(command)
                    require(result.returncode == 0, result.stderr.decode())
                    payloads.setdefault(backend, result.stdout)
                    require(payloads[backend] == result.stdout, "nondeterministic")
                    report["timing"].append(
                        {
                            "request": name,
                            "backend": backend,
                            "repeat": repeat,
                            "seconds": seconds,
                            "events": json.loads(result.stdout)["events"],
                        }
                    )
            require(len(set(payloads.values())) == 1, f"{name} outputs differ")
        files = rules_files(
            folder / "timing-mtdna",
            mt(),
            method="mtdna_direct",
            n=128,
            times=GRID,
            threads=1,
            seed=17,
            trackers=TRACKERS,
        )
        for repeat in range(args.repeats):
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
        grouped = {}
        for row in report["timing"]:
            grouped.setdefault(f"{row['request']}/{row['backend']}", []).append(row)
        summary = {}
        for key, rows in grouped.items():
            seconds = sorted(r["seconds"] for r in rows)[len(rows) // 2]
            summary[key] = {
                "median_seconds": seconds,
                "events": rows[0]["events"],
                "ns_per_event": 1e9 * seconds / rows[0]["events"],
            }
            ns = summary[key]["ns_per_event"]
            print(f"{key:28s} {seconds:.3f} s {ns:.1f} ns/event", flush=True)
        report["summary"] = summary
        ns = {k: v["ns_per_event"] for k, v in summary.items()}
        report["loop_speedup_generic"] = (
            ns["generic/interpreter"] / ns["generic/unrolled"]
        )
        report["loop_speedup_default"] = (
            ns["default/interpreter"] / ns["default/unrolled"]
        )
        report["loop_gap_to_mtdna_direct"] = (
            ns["generic/unrolled"] / ns["special/mtdna_direct"]
        )
        print(
            json.dumps({k: report[k] for k in report if "speedup" in k or "gap" in k})
        )
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
