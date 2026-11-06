"""Compiled Static Direct event loops from Python, against the interpreter.

Builds the event loop of the three-class mtDNA model with the takeover stop
through `rules.compile_loop`, against this checkout (`--source`, default
`CRN_CRAFT_SOURCE` or the repository root), and requires the metadata and
every record column to be byte-identical to the interpreter's:

- 1 and 2 threads, inverse and Ziggurat waits, dense and sparse trackers,
  each with three trackers (one with a first-passage predicate), an event
  log, kept paths and kept instances; at least one cell takes over;
- the plain `"static_direct"` method, with the defaults of the request
  version `rules.ensemble` writes, trackers and an event log;
- event, sampling, instance and output budgets, a wall-time budget that is
  not reached, and a cancel file that exists before the run.

It checks the refusals: other methods, models outside the generator's forms
(expression rates, `propensity` rules, continuous variables, ATP blocks,
open models), callbacks, a missing cargo, a failed build, a loop for other
model bytes or a bound candidate, a request the program rejects, and a
program whose result differs. A second `compile_loop` call must reuse the
cached build.

Without cargo the check is skipped with a message, unless `--require-cargo`
is given. The report holds no timings or paths, so runs under `-O` compare
equal. Builds go to the user cache directory unless `--cache-dir` is given;
set `CARGO_TARGET_DIR` to reuse an existing release build of the checkout.
"""

import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import stat
import sys
import tempfile
import time

from crn_craft import CrnError, _rules_loop, rules

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "models" / "rules"
STOP = "N > 0 && mutants >= 0.9 * N"
GRID = [24.0 * day for day in range(0, 1081, 30)]
TRACKERS = [
    {"observable": "neutral_mutants", "keep_episodes": 2},
    {"observable": "advantaged_mutants"},
    {"observable": "mutants", "first_passage": "N > 0 && mutants >= 0.5 * N"},
]
LOG = {"max_entries": 8, "rules": ["mutation_mt", "mutation_mtAdv"]}


def require(condition, message):
    if not condition:
        raise SystemExit(f"check_compiled_loop: {message}")


def read(name):
    return json.loads((MODELS / name).read_text(encoding="utf-8"))


def takeover():
    model = read("mtdna-three-class.rules.json")
    model["stop"] = STOP
    return model


def identical(label, run, loop, rows, cancel_file=None):
    cancel = None if cancel_file is None else str(cancel_file)
    interpreted = run._native.simulate(cancel, None)
    request = run._compiled_request(loop)
    compiled = loop.execute(request, cancel_file=cancel_file)
    require(
        interpreted[0] == compiled[0] and interpreted[1] == compiled[1],
        f"{label}: compiled result differs from the interpreter's",
    )
    # The public path gives the same result.
    result = run.simulate(compiled=loop, cancel_file=cancel_file)
    require(
        result.metadata == json.loads(compiled[0])["result"], f"{label}: simulate()"
    )
    meta = result.metadata
    row = {
        "case": label,
        "identical": True,
        "termination": meta["termination"],
        "events": meta["events"],
        "stopped": meta["stopped_replicates"],
        "sampler": meta["sampler"],
    }
    rows.append(row)
    print(json.dumps(row), file=sys.stderr, flush=True)
    return result


def refused(label, action, rows, expected):
    """Run `action`; it must raise `expected` (a code and path, or a type)."""
    try:
        action()
    except CrnError as error:
        found = (error.diagnostics[0]["code"], error.diagnostics[0]["path"])
    except rules.CompiledLoopError as error:
        found = ("CompiledLoopError", error.reason)
        if error.reason == "build_failed":
            require(error.log is not None and error.log.is_file(), f"{label}: log")
            require(str(error.log) in str(error), f"{label}: log path in message")
    except (TypeError, ValueError) as error:
        found = (type(error).__name__, "")
    else:
        found = None
    require(found == expected, f"{label}: expected {expected}, got {found}")
    rows.append({"case": label, "refused": list(expected)})


def executable(folder, name, body):
    path = folder / name
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def ensemble(model, loop=None, **settings):
    options = {
        "end": 25_920.0,
        "report_times": GRID,
        "replicates": 24,
        "seed": 11,
        "threads": 2,
        "method": "static_direct",
    }
    return rules.ensemble(model, **(options | settings), compiled=loop)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--require-cargo", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if shutil.which("cargo") is None:
        message = "cargo not found; install Rust (https://rustup.rs) to run it"
        require(not args.require_cargo, message)
        print(f"check_compiled_loop: skipped, {message}", file=sys.stderr)
        report = {"schema": "crn-craft.compiled-loop-check/0.1", "skipped": message}
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        return
    source = args.source or Path(os.environ.get(_rules_loop.SOURCE_ENV, ROOT))
    build = {"source": source, "cache_dir": args.cache_dir}
    model = takeover()

    start = time.perf_counter()
    loop = rules.compile_loop(model, **build)
    print(
        f"check_compiled_loop: compile_loop took {time.perf_counter() - start:.1f} s"
        f" (cargo {loop.build_seconds}, cached {loop.cached})",
        file=sys.stderr,
    )
    again = rules.compile_loop(model, progress=False, **build)
    require(again.cached and again.path == loop.path, "second build not cached")
    require(
        _rules_loop.fingerprint_of(rules.loop_source(model)) == loop.fingerprint,
        "loop_source fingerprint",
    )
    rows = []

    # Threads, both waiting-time policies and both tracker paths.
    for threads in (1, 2):
        for wait in ("inverse", "ziggurat"):
            for sparse in (False, True):
                method = {
                    "kind": "static_direct",
                    "waiting_time": wait,
                    "sparse_trackers": sparse,
                }
                run = ensemble(
                    model,
                    loop,
                    threads=threads,
                    method=method,
                    trackers=TRACKERS,
                    event_log=LOG,
                    keep_paths=3,
                    keep_instances=2,
                )
                identical(
                    f"threads {threads}, {wait}, sparse {sparse}", run, loop, rows
                )
    require(any(r["stopped"] > 0 for r in rows), "no cell took over")
    # The plain method name, with the defaults of the request version that
    # rules.ensemble writes for it.
    default = ensemble(model, loop, trackers=TRACKERS, event_log=LOG)
    identical("default static_direct", default, loop, rows)

    # Budgets and cancellation.
    for label, settings in (
        ("event budget", {"max_events": 20_000}),
        ("sampling budget", {"max_sampling_steps": 30_000}),
        ("instance budget", {"max_instances": 1}),
        ("output budget", {"keep_paths": "all", "max_output_bytes": 24 * 400}),
        ("wall time not reached", {"wall_time_s": 600.0}),
    ):
        identical(label, ensemble(model, loop, **settings), loop, rows)
    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        cancel = scratch / "cancel"
        cancel.write_text("stop\n")
        identical("cancel file", ensemble(model, loop), loop, rows, cancel_file=cancel)
        late = ensemble(model, loop, wall_time_s=1e-6).simulate()
        require(late.metadata["termination"] == "wall_time_limit", "wall time reached")
        rows.append({"case": "wall time reached", "termination": "wall_time_limit"})

        # Refusals with codes.
        run = ensemble(model)
        for method in ("next_reaction", "mtdna_direct", None):
            refused(
                f"method {method}",
                lambda m=method: ensemble(model, loop, method=m),
                rows,
                ("CRN-R020", "/method"),
            )
        expression = takeover()
        expression["rules"][0]["rate"] = "delta * (1 + mutants / K)"
        propensity = rules.Model("birth_death")
        propensity.add_parameter("k", 10.0).add_species("A")
        propensity.add_reaction("k", [], ["A"], name="birth")
        propensity.add_reaction(1.0, ["A"], [], name="death")
        for label, other in (
            ("expression rate", expression),
            ("propensity rule", propensity),
            (
                "continuous variables",
                read("mtdna-three-class-atp-continuous.rules.json"),
            ),
            ("ATP block", read("mtdna-three-class-atp.rules.json")),
            ("open model", read("mtdna-fixed.rules.json")),
        ):
            refused(
                f"loop_source: {label}",
                lambda m=other: rules.loop_source(m),
                rows,
                ("CRN-R020", "/method"),
            )
            refused(
                f"compile_loop: {label}",
                lambda m=other: rules.compile_loop(m, cargo="cargo-missing", **build),
                rows,
                ("CRN-R020", "/method"),
            )
        refused(
            "callback",
            lambda: run.simulate(compiled=loop, callback=lambda i, e: None),
            rows,
            ("ValueError", ""),
        )
        refused(
            "batch_callback",
            lambda: run.simulate(compiled=loop, batch_callback=lambda i, e: None),
            rows,
            ("ValueError", ""),
        )
        refused(
            "not a compiled loop",
            lambda: run.simulate(compiled="loop"),
            rows,
            ("TypeError", ""),
        )
        other = takeover()
        other["stop"] = "N > 0 && mutants >= 0.8 * N"
        refused(
            "loop for other model bytes",
            lambda: ensemble(other, loop),
            rows,
            ("CRN-R020", "/model"),
        )
        refused(
            "bound candidate",
            lambda: run.bind_parameters({"mu": 3e-5}).simulate(compiled=loop),
            rows,
            ("CRN-R020", "/model"),
        )
        refused(
            "request for other model bytes",
            lambda: loop.execute(ensemble(other)._documents[1]),
            rows,
            ("CRN-R020", "/model/sha256"),
        )
        refused(
            "cargo missing",
            lambda: rules.compile_loop(model, cargo=scratch / "no-cargo", **build),
            rows,
            ("CompiledLoopError", "cargo_missing"),
        )
        if os.name == "posix":
            fake = executable(
                scratch,
                "cargo",
                "import sys\nprint('error: simulated build failure', file=sys.stderr)\n"
                "sys.exit(101)\n",
            )
            # A model of its own, so the failure cannot touch the cached build.
            failing = takeover()
            failing["name"] = "build_failure"
            refused(
                "build failure",
                lambda: rules.compile_loop(
                    failing, cargo=fake, cache_dir=scratch / "cache", progress=False
                ),
                rows,
                ("CompiledLoopError", "build_failed"),
            )
            # A program whose result differs in one byte is refused.
            tampered = copy.copy(loop)
            tampered.directory = scratch / "tampered"
            tampered.path = executable(
                scratch,
                "tampered-loop",
                "import subprocess, sys\n"
                f"done = subprocess.run([{str(loop.path)!r}, *sys.argv[1:]],\n"
                "    stdin=sys.stdin, capture_output=True)\n"
                "out = bytearray(done.stdout)\nout[-1] ^= 1\n"
                "sys.stdout.buffer.write(out)\nsys.exit(done.returncode)\n",
            )
            refused(
                "differing program",
                lambda: tampered.verify(run),
                rows,
                ("CompiledLoopError", "mismatch"),
            )
            require(
                (scratch / "tampered" / "mismatch" / "compiled.json").is_file(),
                "mismatch payloads saved",
            )
    require(loop.verify(run).metadata["events"] > 0, "verify")
    report = {
        "schema": "crn-craft.compiled-loop-check/0.1",
        "model": "mtdna-three-class.rules.json with the takeover stop",
        "cached_second_call": True,
        "rows": rows,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"passed": True, "rows": len(rows)}))


if __name__ == "__main__":
    main()
