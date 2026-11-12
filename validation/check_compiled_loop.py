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
  not reached, and a cancel file that exists before a one-worker run.

String rates: four more models, written as JSON, are compiled and must be
byte-identical to the interpreter with 1 and 2 threads, inverse and
Ziggurat waits, trackers (one with a first-passage predicate), an event log
and the model's stop:

- the SIR model as `rules.Model` writes it, with the string rates
  `"beta * S * I / N"` and `"gamma * I"` (whole-rule propensities), and a
  stop at extinction or a final size;
- Michaelis-Menten conversion `V * A / (Km + A)` beside constant rates;
- the three-class model with degradation rates `delta * (1 + mutants / K)`,
  which read observables, and the takeover stop;
- rates with `max`, `min`, `if` and the other functions over `distinct`,
  `max` and `min` observables;
- higher-order rules of a `crn-craft.model/1.0` model (`3A -> 0`,
  `2A + B -> C` with an expression rate, a whole propensity over `3A + B`
  and a rule of three patterns), run with `crn-craft.run/1.0` requests.

It checks the refusals: other methods, models that only another method runs
(continuous variables, ATP blocks, open models), callbacks, a missing cargo,
a failed build, a loop for other model bytes or a bound candidate, a
request the program rejects, and a program whose result differs. A second
`compile_loop` call must reuse the cached build. The generated Cargo
project must depend on the `crn-craft` crate alone, at exactly the
installed version: from the checkout, and from crates.io when no source is
given (that manifest is checked as written, not built).

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
import tomllib

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


def sir():
    """The SIR model as `rules.Model` writes it, plus a stop."""
    return {
        "schema": "crn-craft.rules-model/0.3",
        "name": "sir",
        "parameters": {"beta": 0.3, "gamma": 0.1, "N": 1000.0},
        "types": [{"name": "S"}, {"name": "I"}, {"name": "R"}],
        "initial": [{"type": "S", "count": 990}, {"type": "I", "count": 10}],
        "observables": {s: {"type": s} for s in ("S", "I", "R")},
        "rules": [
            {
                "name": "infection",
                "reactants": [
                    {"label": "x", "type": "S"},
                    {"label": "y", "type": "I"},
                ],
                "propensity": "beta * S * I / N",
                "products": [{"type": "I"}, {"type": "I"}],
            },
            {
                "name": "recovery",
                "reactants": [{"label": "x", "type": "I"}],
                "propensity": "gamma * I",
                "products": [{"type": "R"}],
            },
        ],
        "stop": "I == 0 || R >= 600",
    }


def michaelis_menten():
    return {
        "schema": "crn-craft.rules-model/0.3",
        "name": "michaelis_menten",
        "parameters": {"V": 50.0, "Km": 20.0, "k_in": 40.0, "d": 0.5},
        "types": [{"name": "A"}, {"name": "P"}],
        "initial": [{"type": "A", "count": 30}],
        "observables": {"A": {"type": "A"}, "P": {"type": "P"}},
        "rules": [
            {
                "name": "influx",
                "reactants": [],
                "rate": "k_in",
                "products": [{"type": "A"}],
            },
            {
                "name": "conversion",
                "reactants": [{"label": "x", "type": "A"}],
                "propensity": "V * A / (Km + A)",
                "products": [{"type": "P"}],
            },
            {
                "name": "decay",
                "reactants": [{"label": "x", "type": "P"}],
                "rate": "d",
                "products": [],
            },
        ],
        "stop": "A >= 140",
    }


def three_class_expression():
    model = takeover()
    for rule in model["rules"]:
        if rule["name"].startswith("degradation"):
            rule["rate"] = "delta * (1 + mutants / K)"
    return model


def functions():
    """Rates with max, min, if and other functions over distinct, max and min
    observables, and a whole-rule identical pair."""
    return {
        "schema": "crn-craft.rules-model/0.3",
        "name": "functions",
        "parameters": {"d": 0.3, "p": 0.02, "q": 0.01, "e": 0.5, "c": 0.2},
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C", "attributes": ["s"]}],
        "initial": [
            {"type": "A", "count": 40},
            {"type": "C", "attributes": {"s": 1}, "count": 5},
            {"type": "C", "attributes": {"s": 2}, "count": 3},
        ],
        "observables": {
            "nA": {"type": "A"},
            "nB": {"type": "B"},
            "w": {"type": "B", "weight": "2.5"},
            "kinds": {"type": "C", "reduce": "distinct"},
            "spread": {"type": "C", "weight": "x.s * x.s", "reduce": "distinct"},
            "top": {"type": "C", "weight": "x.s", "reduce": "max"},
            "bottom": {"type": "C", "weight": "-x.s", "reduce": "min"},
        },
        "rules": [
            {
                "name": "death_a",
                "reactants": [{"label": "x", "type": "A"}],
                "rate": "d * (1 + floor(nA / 10) / 4) * if(top > 2, 1.25, 1)",
                "products": [],
            },
            {
                "name": "pair_aa",
                "reactants": [
                    {"label": "x", "type": "A"},
                    {"label": "y", "type": "A"},
                ],
                "propensity": "p * nA * (nA - 1) / 2 * (nB < 9 || !(w > 4))",
                "products": [{"type": "B"}],
            },
            {
                "name": "pair_ab",
                "reactants": [
                    {"label": "x", "type": "A"},
                    {"label": "y", "type": "B"},
                ],
                "rate": "q * exp(-nB / 20) * sqrt(1 + nA) / ceil(sqrt(1 + nA))",
                "products": [{"type": "A"}, {"type": "A"}],
            },
            {
                "name": "death_b",
                "reactants": [{"label": "x", "type": "B"}],
                "rate": "e * max(1, log(1 + nB)) ^ 0.5 - 0 * abs(bottom)",
                "products": [],
            },
            {
                "name": "create_c3",
                "reactants": [],
                "propensity": "if(kinds >= 3 && bottom != -3, 0.05,"
                " 0.5 * abs(top - 3) + min(0.1, spread / 100))",
                "products": [{"type": "C", "attributes": {"s": 3}}],
            },
            {
                "name": "death_c",
                "reactants": [{"label": "x", "type": "C"}],
                "rate": "c + (nA == 0) * 0.3 + (top <= 2) * 0.1"
                " + (spread != 5) * 0.05 + 2 ^ -3 * (kinds >= 2)",
                "products": [],
            },
        ],
        "stop": "nB >= 9 || (nA == 0 && w > 4) || (kinds == 0 && nA < 5)",
    }


def higher_order():
    """Higher-order rules (crn-craft.model/1.0): pattern counts and three
    patterns, a rate, an expression rate and a whole propensity."""
    return {
        "schema": "crn-craft.model/1.0",
        "name": "higher_order",
        "parameters": {
            "kin": 30.0,
            "kb": 4.0,
            "c3": 0.0005,
            "c2": 0.001,
            "w": 2.0,
            "dc": 0.05,
        },
        "types": [{"name": "A"}, {"name": "B"}, {"name": "C"}],
        "initial": [{"type": "A", "count": 20}, {"type": "B", "count": 5}],
        "observables": {"nA": {"type": "A"}, "nB": {"type": "B"}, "nC": {"type": "C"}},
        "rules": [
            {
                "name": "source_a",
                "reactants": [],
                "rate": "kin",
                "products": [{"type": "A"}],
            },
            {
                "name": "source_b",
                "reactants": [],
                "rate": "kb",
                "products": [{"type": "B"}],
            },
            {
                "name": "third_order",
                "reactants": [{"label": "a", "type": "A", "count": 3}],
                "rate": "c3",
                "products": [],
            },
            {
                "name": "bind",
                "reactants": [
                    {"label": "a", "type": "A", "count": 2},
                    {"label": "b", "type": "B"},
                ],
                "rate": "c2 * (1 + nC / 50)",
                "products": [{"type": "C"}],
            },
            {
                "name": "whole",
                "reactants": [
                    {"label": "a", "type": "A", "count": 3},
                    {"label": "b", "type": "B"},
                ],
                "propensity": "w * nB / (1 + nB)",
                "products": [{"type": "A"}, {"type": "C"}],
            },
            {
                "name": "three_patterns",
                "reactants": [
                    {"label": "x", "type": "A"},
                    {"label": "y", "type": "B"},
                    {"label": "z", "type": "C"},
                ],
                "rate": "c2",
                "products": [{"type": "C"}, {"type": "C"}],
            },
            {
                "name": "decay_c",
                "reactants": [{"label": "c", "type": "C"}],
                "rate": "dc",
                "products": [],
            },
        ],
        "stop": "nC >= 60",
    }


# label, model, run settings, trackers, event log
STRING_RATES = (
    (
        "string-rate SIR",
        sir,
        {"end": 300.0, "report_times": [30.0 * i for i in range(11)]},
        [
            {"observable": "I", "keep_episodes": 2},
            {"observable": "R", "first_passage": "R >= 300"},
        ],
        {"max_entries": 6, "rules": ["recovery"]},
    ),
    (
        "Michaelis-Menten",
        michaelis_menten,
        {"end": 40.0, "report_times": [0.0, 0.5, 1.0, 10.0, 20.0, 40.0]},
        [{"observable": "A"}, {"observable": "P", "first_passage": "P > 60"}],
        {"max_entries": 5, "rules": ["conversion"]},
    ),
    (
        "three-class expression rates",
        three_class_expression,
        {"end": 25_920.0, "report_times": GRID},
        TRACKERS,
        LOG,
    ),
    (
        "max, min and if",
        functions,
        {"end": 40.0, "report_times": [0.0, 1e-6, 1.0, 10.0, 40.0]},
        [
            {"observable": "top", "keep_episodes": 2},
            {"observable": "spread"},
            {"observable": "nB", "first_passage": "nB >= 3"},
        ],
        {"max_entries": 6, "rules": ["create_c3", "pair_aa"]},
    ),
    (
        "higher-order rules",
        higher_order,
        {"end": 60.0, "report_times": [0.0, 1.0, 15.0, 30.0, 60.0]},
        [
            {"observable": "nA", "keep_episodes": 2},
            {"observable": "nC", "first_passage": "nC >= 10"},
        ],
        {"max_entries": 6, "rules": ["third_order", "bind", "three_patterns"]},
    ),
)


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

    # The loop project depends on the crn-craft crate alone, at exactly the
    # installed version: from the checkout here, and from crates.io when no
    # source is given (the manifest of a published install).
    project = tomllib.loads((loop.directory / "Cargo.toml").read_text("utf-8"))
    checkout = (source.expanduser().resolve() / "crates" / "crn-craft").as_posix()
    require(
        project["dependencies"]
        == {
            "crn-craft": {"version": f"={loop.version}", "path": checkout},
            "serde_json": "1",
        },
        f"loop project dependencies: {project['dependencies']}",
    )
    rust_version, _ = _rules_loop._native.rules_loop_driver()
    published = tomllib.loads(
        _rules_loop._manifest("crn-craft-loop-check", loop.version, rust_version, None)
    )
    require(
        published["dependencies"]
        == {"crn-craft": {"version": f"={loop.version}"}, "serde_json": "1"},
        f"crates.io loop project dependencies: {published['dependencies']}",
    )
    rows.append(
        {
            "case": "loop project dependencies",
            "checkout": ["crn-craft (path, =installed version)", "serde_json"],
            "crates.io": ["crn-craft (=installed version)", "serde_json"],
        }
    )

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

    # Quantiles need crn-craft.run/1.0 and run the model as
    # crn-craft.model/1.0: a loop compiled for that document serves them,
    # byte-identically, also through verify().
    quantiles = (0.025, 0.5, 0.975)
    unified = {**model, "schema": rules.MODEL_1_0}
    unified_loop = rules.compile_loop(unified, progress=False, **build)
    for threads in (1, 2):
        quantile_run = ensemble(
            model,
            unified_loop,
            threads=threads,
            quantiles=quantiles,
            keep_paths="all",
            trackers=TRACKERS,
        )
        result = identical(
            f"quantiles, threads {threads}", quantile_run, unified_loop, rows
        )
        require(
            result.metadata["schema"] == "crn-craft.result/1.0"
            and list(result.statistics["quantiles"]) == list(quantiles)
            and result.statistics["median"] is result.statistics["quantiles"][0.5],
            f"quantiles, threads {threads}: quantile statistics",
        )
    verified = unified_loop.verify(quantile_run)
    require(
        verified.metadata == result.metadata
        and list(verified.statistics["quantiles"]) == list(quantiles),
        "verify() of a quantile run",
    )

    # String rates and whole-rule propensities.
    for label, build_model, settings, trackers, log in STRING_RATES:
        other = build_model()
        other_loop = rules.compile_loop(other, progress=False, **build)
        cells = []
        for threads in (1, 2):
            for wait in ("inverse", "ziggurat"):
                method = {"kind": "static_direct", "waiting_time": wait}
                run = rules.ensemble(
                    other,
                    replicates=24,
                    seed=13,
                    threads=threads,
                    method=method,
                    trackers=trackers,
                    event_log=log,
                    keep_paths=3,
                    keep_instances=2,
                    compiled=other_loop,
                    **settings,
                )
                cells.append(
                    identical(
                        f"{label}: threads {threads}, {wait}", run, other_loop, rows
                    )
                )
        require(
            any(c.metadata["stopped_replicates"] > 0 for c in cells),
            f"{label}: the stop never fired",
        )

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
        # One worker: with more, which replicates start before the first poll
        # sees the file depends on timing, in the interpreter too.
        cancelled = ensemble(model, loop, threads=1)
        identical("cancel file", cancelled, loop, rows, cancel_file=cancel)
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
        for label, other in (
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
            "quantiles with the rules-model loop",
            lambda: ensemble(model, loop, quantiles=(0.5,)),
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
