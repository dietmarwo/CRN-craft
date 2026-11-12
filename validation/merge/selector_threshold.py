"""Declared benchmark for the `auto` selector of Static Direct (M4c).

Synthetic `crn-craft.model/1.0` families over a grid of channel counts, each
run with every selector through the CLI under test (`crn-craft rules`, one
process per run, one replicate, one thread, root seed 4242, Ziggurat waits).
The horizon of a run is chosen from the initial total propensity a0 so that
about `--events` events (default 1,000,000) fire (end = events / a0); the time of a point is the
median wall time of three runs after a warm-up, divided by the events of the
run (identical in every repetition).

Families (declared in `validation/merge/RESULTS.md`, "M4c ... declaration"):

- `ring`: n species X_i, n channels X_i -> X_(i+1 mod n) at rate 1, 100
  molecules each (sparse coupling, constant rates, wide intervals);
- `ring_small`: the same with 4 molecules each (narrow intervals);
- `dense`: m species, every ordered pair X_i -> X_j at rate 1 (n = m(m - 1)
  channels), 32 molecules each (dense coupling);
- `michaelis`: the ring with the whole propensity vmax X_i / (km + X_i),
  vmax 1, km 10, 100 molecules each (expression rates of one species);
- `two_reads`: the ring with the whole propensity X_i X_(i+1) / (1 + X_(i+1)),
  100 molecules each (rates that read two species: RSSA-CR pins them).

The rule fixed from the result (also declared):

- N_cr is the smallest grid count g such that composition-rejection is
  faster than the linear scan at every point with n >= g of every family;
- N_rssa is the smallest grid count g such that RSSA-CR is faster than both
  the linear scan and composition-rejection at every point with n >= g of
  every family except `two_reads`;
- `auto` picks RSSA-CR when every rule's rate reads at most one species and
  n >= N_rssa, else composition-rejection when n >= N_cr, else the linear
  scan; sorting is never picked. A count no grid point satisfies is infinite.

Only the standard library is needed. Run it under the CPU lock:

    flock /home/claude/craft5/cpu.lock python3 validation/merge/selector_threshold.py \\
        --binary target/release/crn-craft --output validation/merge/selector-threshold.json
"""

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

GRID = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
DENSE = [3, 4, 6, 8, 12, 16, 23, 32, 46, 64]
SELECTORS = ["linear", "composition_rejection", "sorting", "rssa_cr"]
FAMILIES = ["ring", "ring_small", "dense", "michaelis", "two_reads"]
INTERVAL_FAMILIES = ["ring", "ring_small", "dense", "michaelis"]
SEED = 4242
REPEATS = 3


def dump(value):
    return (json.dumps(value, separators=(",", ":")) + "\n").encode()


def ring_model(family, n):
    """The model of a ring family with n channels, and its total a0."""
    count = 4 if family == "ring_small" else 100
    names = [f"X{i}" for i in range(n)]
    rules = []
    for i, name in enumerate(names):
        after = names[(i + 1) % n]
        rule = {
            "name": f"r{i}",
            "reactants": [{"label": "x", "type": name}],
            "products": [{"type": after}],
        }
        if family == "michaelis":
            rule["propensity"] = f"vmax * {name} / (km + {name})"
        elif family == "two_reads":
            rule["propensity"] = f"{name} * {after} / (1 + {after})"
        else:
            rule["rate"] = "1"
        rules.append(rule)
    a0 = {
        "ring": n * count,
        "ring_small": n * count,
        "michaelis": n * count / (10 + count),
        "two_reads": n * count * count / (1 + count),
    }[family]
    model = {
        "schema": "crn-craft.model/1.0",
        "name": f"{family}_{n}",
        "parameters": {"vmax": 1.0, "km": 10.0},
        "types": [{"name": x} for x in names],
        "initial": [{"type": x, "count": count} for x in names],
        "observables": {x: {"type": x} for x in names},
        "rules": rules,
    }
    return model, a0


def dense_model(m):
    names = [f"X{i}" for i in range(m)]
    rules = [
        {
            "name": f"r{i}_{j}",
            "reactants": [{"label": "x", "type": names[i]}],
            "rate": "1",
            "products": [{"type": names[j]}],
        }
        for i in range(m)
        for j in range(m)
        if i != j
    ]
    model = {
        "schema": "crn-craft.model/1.0",
        "name": f"dense_{m}",
        "parameters": {},
        "types": [{"name": x} for x in names],
        "initial": [{"type": x, "count": 32} for x in names],
        "observables": {x: {"type": x} for x in names},
        "rules": rules,
    }
    return model, (m - 1) * 32 * m


def points(family):
    if family == "dense":
        for m in DENSE:
            yield m * (m - 1), *dense_model(m)
    else:
        for n in GRID:
            yield n, *ring_model(family, n)


def request(model_bytes, end, selector, types):
    return {
        "schema": "crn-craft.run/1.0",
        "model": {
            "path": "model.json",
            "sha256": hashlib.sha256(model_bytes).hexdigest(),
        },
        "time": {"end": end, "report_times": [0.0, end]},
        "ensemble": {"replicates": 1, "root_seed": SEED},
        "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
        "budget": {
            "max_events_per_replicate": 10**9,
            "max_sampling_steps_per_replicate": 10**12,
            "max_instances": types,
            "threads": 1,
            "max_output_bytes": 1 << 30,
        },
        "method": {"kind": "static_direct", "selector": selector},
    }


def timed(binary, folder):
    runs = []
    output = None
    for i in range(REPEATS + 1):
        start = time.perf_counter()
        done = subprocess.run(
            [str(binary), "rules", "model.json", "run.json"],
            cwd=folder,
            capture_output=True,
            check=False,
        )
        seconds = time.perf_counter() - start
        if done.returncode != 0:
            raise RuntimeError(done.stderr.decode()[:400])
        if output is not None and done.stdout != output:
            raise RuntimeError("runs printed different bytes")
        output = done.stdout
        if i:
            runs.append(seconds)
    return statistics.median(runs), runs, json.loads(output)


def threshold(rows, families, faster):
    """The smallest grid count from which `faster(row)` holds at every point
    with at least that many channels of every family; None if none."""
    for g in GRID:
        if all(
            faster(r) for r in rows if r["family"] in families and r["channels"] >= g
        ):
            return g
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--events", type=int, default=1_000_000)
    parser.add_argument("--families", default=",".join(FAMILIES))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    report = {
        "schema": "crn-craft.selector-threshold/0.1",
        "declared_in": "validation/merge/RESULTS.md",
        "host": {"machine": platform.machine(), "python": platform.python_version()},
        "events_per_run": args.events,
        "repeats": REPEATS,
        "boundary": "wall time of one crn-craft rules process, median of three after a warm-up",
        "rows": [],
    }
    for family in args.families.split(","):
        for channels, model, a0 in points(family):
            model_bytes = dump(model)
            end = args.events / a0
            row = {"family": family, "channels": channels, "end": end, "selectors": {}}
            with tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                (folder / "model.json").write_bytes(model_bytes)
                for selector in SELECTORS:
                    (folder / "run.json").write_bytes(
                        dump(request(model_bytes, end, selector, len(model["types"])))
                    )
                    seconds, runs, out = timed(binary, folder)
                    events = out["events"]
                    row["selectors"][selector] = {
                        "seconds": seconds,
                        "runs": runs,
                        "events": events,
                        "ns_per_event": 1e9 * seconds / max(events, 1),
                        "sampler": out["sampler"],
                    }
            report["rows"].append(row)
            line = "  ".join(
                f"{s[:6]} {row['selectors'][s]['ns_per_event']:8.1f}" for s in SELECTORS
            )
            print(f"{family:10s} {channels:5d}  {line}  ns/event", flush=True)
            args.output.write_text(json.dumps(report, indent=1) + "\n")

    def ns(row, selector):
        return row["selectors"][selector]["ns_per_event"]

    rows = report["rows"]
    n_cr = threshold(
        rows, FAMILIES, lambda r: ns(r, "composition_rejection") < ns(r, "linear")
    )
    n_rssa = threshold(
        rows,
        INTERVAL_FAMILIES,
        lambda r: ns(r, "rssa_cr")
        < min(ns(r, "linear"), ns(r, "composition_rejection")),
    )
    report["rule"] = {
        "N_cr": n_cr,
        "N_rssa": n_rssa,
        "auto": "rssa_cr if every rule reads at most one species and channels >= N_rssa; "
        "else composition_rejection if channels >= N_cr; else linear "
        "(null: no grid count qualifies, the branch never applies)",
    }
    print(json.dumps(report["rule"]), flush=True)
    args.output.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
