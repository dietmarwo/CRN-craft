"""Declared benchmark and fit of the `auto` selector of Static Direct (M4f).

Declared in `validation/merge/RESULTS.md` ("M4f ... declaration", item 4).
Synthetic `crn-craft.model/1.0` families over the M4c grid of channel
counts, each point run with `linear`, `composition_rejection` and `rssa_cr`
through the CLI under test as in M4c (`selector_threshold.py`): one
replicate, one thread, root seed 4242, Ziggurat waits, about `--events`
events, the median wall time of three processes after a warm-up, divided by
the events.

Families: `ring4`, `ring32`, `ring256` (X_i -> X_(i+1 mod n) at rate 1 with
4, 32 or 256 molecules per species), `dense4`, `dense32`, `dense256` (every
X_i -> X_j, i != j, at rate 1), `michaelis` (the ring with the whole
propensity vmax X_i / (km + X_i), 100 molecules) and `two_reads` (the ring
with X_i X_(i+1) / (1 + X_(i+1)), 100 molecules; its rules read two
species).

`--fit SPEED_JSON` (a `compare_speed.py` report with the selectors linear,
composition_rejection, sorting and rssa_cr for W1-W7 and the examples) then
chooses the rule

    auto = rssa_cr                if not pinned and n >= (N_small if P < P_cut else N_large)
           composition_rejection  else if n >= N_cr
           linear                 otherwise

over the declared parameter grid by the declared criteria (a)-(e), where n
is the lowered channel count, `pinned` whether some rule's rate reads two or
more lowered species, and P the mean initial count of the lowered species
with a positive initial count (infinite if there is none).

Only the standard library is needed. Run it under the CPU lock:

    flock /home/claude/craft5/cpu.lock python3 validation/merge/selector_rule.py \\
        --binary target/release/crn-craft --output validation/merge/selector-rule.json
    python3 validation/merge/selector_rule.py --output validation/merge/selector-rule.json \\
        --fit validation/merge/speed-m4f-selectors.json
"""

import argparse
import itertools
import json
import math
import platform
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from selector_threshold import DENSE, GRID, SEED, dump, timed  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SELECTORS = ["linear", "composition_rejection", "rssa_cr"]
FAMILIES = [
    "ring4",
    "ring32",
    "ring256",
    "dense4",
    "dense32",
    "dense256",
    "michaelis",
    "two_reads",
]
THRESHOLDS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, math.inf]
CUTS = [16, 8, 32]
NEAR = 1.10
TARGETS = [
    "W1",
    "W2",
    "W3",
    "W4",
    "W5",
    "W6",
    "W7",
    "enzyme-mechanism",
    "immigration-death",
    "michaelis-menten",
    "repressilator",
    "vilar-noise",
]


def model_of(names, rules, count, name, parameters=None):
    return {
        "schema": "crn-craft.model/1.0",
        "name": name,
        "parameters": parameters or {},
        "types": [{"name": x} for x in names],
        "initial": [{"type": x, "count": count} for x in names],
        "observables": {x: {"type": x} for x in names},
        "rules": rules,
    }


def ring(family, n):
    """A ring family's model with n channels, its total a0 and molecules."""
    count = {"ring4": 4, "ring32": 32, "ring256": 256}.get(family, 100)
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
        "michaelis": n * count / (10 + count),
        "two_reads": n * count * count / (1 + count),
    }.get(family, n * count)
    parameters = {"vmax": 1.0, "km": 10.0} if family == "michaelis" else {}
    return model_of(names, rules, count, f"{family}_{n}", parameters), a0, count


def dense(family, m):
    count = {"dense4": 4, "dense32": 32, "dense256": 256}[family]
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
    return model_of(names, rules, count, f"{family}_{m}"), (m - 1) * count * m, count


def points(family):
    if family.startswith("dense"):
        for m in DENSE:
            yield m * (m - 1), *dense(family, m)
    else:
        for n in GRID:
            yield n, *ring(family, n)


def request(model_bytes, end, selector, types):
    import hashlib

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


def bench(args, report):
    done = {(r["family"], r["channels"]) for r in report["rows"]}
    for family in args.families.split(","):
        for channels, model, a0, count in points(family):
            if (family, channels) in done:
                continue
            model_bytes = dump(model)
            end = args.events / a0
            row = {
                "family": family,
                "channels": channels,
                "population": count,
                "pinned": family == "two_reads",
                "end": end,
                "selectors": {},
            }
            with tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                (folder / "model.json").write_bytes(model_bytes)
                for selector in SELECTORS:
                    (folder / "run.json").write_bytes(
                        dump(request(model_bytes, end, selector, len(model["types"])))
                    )
                    seconds, runs, out = timed(args.binary, folder)
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


def features(converted):
    """n, pinned and P of an attribute-free converted model, as the engine
    lowers it: a type is a lowered species if it starts with molecules or is
    a product; a rule is a channel if every reactant type is lowered."""
    initial = {e["type"]: e["count"] for e in converted["initial"]}
    lowered = {t for t, c in initial.items() if c > 0}
    for rule in converted["rules"]:
        for product in rule["products"]:
            lowered.add(product["type"])
    observables = converted["observables"]
    channels = 0
    pinned = False
    for rule in converted["rules"]:
        if all(p["type"] in lowered for p in rule["reactants"]):
            channels += 1
        text = rule.get("propensity", rule.get("rate", ""))
        names = set()
        for token in "".join(
            c if c.isalnum() or c == "_" else " " for c in text
        ).split():
            if token in observables and observables[token]["type"] in lowered:
                names.add(observables[token]["type"])
        pinned |= len(names) > 1
    positive = [c for c in initial.values() if c > 0]
    population = sum(positive) / len(positive) if positive else math.inf
    return channels, pinned, population


def pick(rule, channels, pinned, population):
    n_small, n_large, n_cr, cut = rule
    if not pinned and channels >= (n_small if population < cut else n_large):
        return "rssa_cr"
    if channels >= n_cr:
        return "composition_rejection"
    return "linear"


def fit(report, speed):
    targets = []
    for entry in speed["models"]:
        if entry["model"] not in TARGETS:
            continue
        times = {s: entry["selectors"][s]["seconds"] for s in entry["selectors"]}
        converted_path = ROOT / next(
            m["path"]
            for m in json.loads((ROOT / "validation/merge/models.json").read_text())[
                "rows"
            ]
            if m["name"] == entry["model"]
        )
        targets.append((entry["model"], times, converted_path))
    missing = set(TARGETS) - {t[0] for t in targets}
    if missing:
        raise SystemExit(f"target timings missing: {sorted(missing)}")
    converted = {}
    for name, _, path in targets:
        converted[name] = json.loads(
            __import__("subprocess")
            .run(
                [str(report["binary"]), "convert", str(path)],
                capture_output=True,
                check=True,
            )
            .stdout
        )
    target_rows = []
    for name, times, _ in targets:
        n, pinned, population = features(converted[name])
        fastest = min(
            times[s] for s in ["linear", "composition_rejection", "sorting", "rssa_cr"]
        )
        target_rows.append((name, n, pinned, population, times, fastest))
    grid_rows = [
        (
            r["channels"],
            r["pinned"],
            r["population"],
            {s: r["selectors"][s]["ns_per_event"] for s in SELECTORS},
        )
        for r in report["rows"]
    ]

    def score(rule):
        misses = 0
        worst = 0.0
        for _, n, pinned, population, times, fastest in target_rows:
            ratio = times[pick(rule, n, pinned, population)] / fastest
            misses += ratio > NEAR
            worst = max(worst, ratio)
        grid_worst = max(
            ns[pick(rule, n, pinned, population)] / min(ns.values())
            for n, pinned, population, ns in grid_rows
        )
        n_small, n_large, n_cr, cut = rule
        return (misses, worst, grid_worst, -n_cr, -n_large, -n_small, CUTS.index(cut))

    rules = itertools.product(THRESHOLDS, THRESHOLDS, THRESHOLDS, CUTS)
    best = min(rules, key=score)
    misses, worst, grid_worst, *_ = score(best)

    def number(x):
        return None if x == math.inf else x

    return {
        "N_small": number(best[0]),
        "N_large": number(best[1]),
        "N_cr": number(best[2]),
        "P_cut": best[3],
        "targets_more_than_1_10_slower": misses,
        "largest_target_ratio": worst,
        "largest_grid_ratio": grid_worst,
        "targets": [
            {
                "model": name,
                "channels": n,
                "pinned": pinned,
                "population": number(population),
                "pick": pick(best, n, pinned, population),
                "ratio": times[pick(best, n, pinned, population)] / fastest,
                "fastest": min(
                    ["linear", "composition_rejection", "sorting", "rssa_cr"],
                    key=lambda s, times=times: times[s],
                ),
            }
            for name, n, pinned, population, times, fastest in target_rows
        ],
        "auto": "rssa_cr if not pinned and channels >= (N_small if population < P_cut "
        "else N_large); else composition_rejection if channels >= N_cr; else linear "
        "(null: infinite)",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--events", type=int, default=1_000_000)
    parser.add_argument("--families", default=",".join(FAMILIES))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--fit", type=Path, help="compare_speed.py report of the targets"
    )
    args = parser.parse_args()
    if args.output.exists():
        report = json.loads(args.output.read_text())
    else:
        report = {
            "schema": "crn-craft.selector-rule/0.1",
            "declared_in": "validation/merge/RESULTS.md",
            "host": {
                "machine": platform.machine(),
                "python": platform.python_version(),
            },
            "events_per_run": args.events,
            "repeats": 3,
            "boundary": "wall time of one crn-craft rules process, median of three after a warm-up",
            "rows": [],
        }
    if args.binary is not None:
        args.binary = args.binary.resolve()
        report["binary"] = str(args.binary)
        if report["events_per_run"] != args.events:
            raise SystemExit("the existing output has another event count")
        bench(args, report)
    if args.fit is not None:
        measured = {(r["family"], r["channels"]) for r in report["rows"]}
        expected = {(f, n) for f in FAMILIES for n, *_ in points(f)}
        if measured != expected:
            raise SystemExit(f"{len(expected - measured)} grid points missing")
        report["rows"].sort(key=lambda r: (FAMILIES.index(r["family"]), r["channels"]))
        report["fit_from"] = str(args.fit)
        report["rule"] = fit(report, json.loads(args.fit.read_text()))
        print(json.dumps(report["rule"], indent=1), flush=True)
    args.output.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
