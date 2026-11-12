"""Exact-law tests of higher-order reactions against FSP (task M4b).

For every row of `validation/merge/higher-order/models.json` (`3A -> 0` and
`2A + B -> C`, network models) this converts the model with the engine under
test (`crn_craft.convert`, which writes reactant multiplicities as `count`),
runs the converted `crn-craft.model/1.0` model with each of the row's exact
methods through `crn_craft.rules` (`crn-craft.run/1.0`, every path kept),
chunk k with root seed `gate_seed + k`, and compares the sample with the
0.6.0 FSP law stored by `higher_order_references.py`:

1. mean and variance of every species at every report time, with the law
   gate's formulas 3 and 4 (`merge_common.compare_to_law`);
2. the probability of every joint state whose FSP probability (normalized by
   the retained mass) is at least the manifest's `state_probability_floor`:
   |p_hat - p| <= 5 sqrt(p (1 - p) / n);
3. the total probability of the other states, with the same bound.

A bound of zero demands equality (report time 0). The declarations are in
`validation/merge/RESULTS.md` (M4b). The report holds every row, failed rows
included, and no timings, so runs under `-O` compare byte for byte. Run it
with the Python that has the wheel under test, under the CPU lock:

    flock /home/claude/craft5/cpu.lock python validation/merge/check_higher_order.py

Exit status: 0 when every comparison passes and every model converts.
"""

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_laws import checked_paths, dump, request_for  # noqa: E402
from merge_common import (  # noqa: E402
    HERE,
    ROOT,
    SIGMAS,
    compare_to_law,
    load_manifest,
    sha256,
    statistics_from_power_sums,
    summarize,
)

MANIFEST = HERE / "higher-order" / "models.json"
REFERENCE = HERE / "reference" / "higher-order"


def probability_row(kind, time, state, observed, p, n):
    """One binomial comparison of a sample frequency with an exact probability."""
    scale = math.sqrt(max(0.0, p * (1.0 - p)) / n)
    value = observed / n
    difference = abs(value - p)
    bound = SIGMAS * scale
    return {
        "comparison": kind,
        "time": time,
        "state": state,
        "count": observed,
        "value": value,
        "reference": p,
        "bound": bound,
        "z": difference / scale if scale > 0 else (0.0 if difference == 0 else None),
        "passed": difference <= bound,
    }


def state_rows(fsp, histograms, times, n, floor):
    """Comparisons 2 and 3 for every report time."""
    states = [tuple(s) for s in fsp["states"]]
    rows = []
    for t, (time, snapshot) in enumerate(zip(times, fsp["snapshots"], strict=True)):
        if snapshot["time"] != time:
            raise RuntimeError(f"FSP snapshot at {snapshot['time']}, expected {time}")
        p = snapshot["probabilities"]
        mass = math.fsum(p)
        seen = 0
        for state, q in zip(states, p, strict=True):
            q /= mass
            if q >= floor:
                observed = histograms[t].get(state, 0)
                seen += observed
                rows.append(probability_row("state", time, list(state), observed, q, n))
        tail = math.fsum(q for q in p if q / mass < floor) / mass
        rows.append(probability_row("other_states", time, None, n - seen, tail, n))
    return rows


def run_method(rules, np, power_sums, converted, document, row, method, species):
    """Power sums and joint-state histograms of one method's chunks."""
    times = row["report_times"]
    chunk = row["chunk"]
    histograms = [{} for _ in times]
    totals = None
    sampler = None
    events = 0
    columns = None
    for k in range(row["replicates"] // chunk):
        request = request_for(
            document,
            row,
            chunk=chunk,
            seed=row["gate_seed"] + k,
            threads=2,
            method=method,
        )
        request["model"]["sha256"] = sha256(converted)
        result = rules.prepare(converted, request).simulate().metadata
        values = checked_paths(result, chunk, times, row["name"])
        if columns is None:
            sampler = result["sampler"]
            columns = [result["observables"].index(s) for s in species]
        if result["sampler"] != sampler:
            raise RuntimeError(f"{row['name']}: sampler changed")
        events += result["events"]
        x = np.asarray(values, dtype=np.float64)[:, :, columns]
        counts = x.astype(np.int64)
        if not np.array_equal(counts, x):
            raise RuntimeError(f"{row['name']}: non-integral counts")
        sums = power_sums(counts, np)
        totals = (
            sums
            if totals is None
            else [a + b for a, b in zip(totals, sums, strict=True)]
        )
        for t in range(len(times)):
            unique, number = np.unique(counts[:, t, :], axis=0, return_counts=True)
            for state, m in zip(unique.tolist(), number.tolist(), strict=True):
                key = tuple(state)
                histograms[t][key] = histograms[t].get(key, 0) + m
    return {
        "sampler": sampler,
        "events": events,
        "request_schema": request["schema"],
        "power_sums": {
            "n": row["replicates"],
            **{f"s{p + 1}": totals[p].tolist() for p in range(4)},
        },
        "histograms": histograms,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--output", type=Path, default=HERE / "laws" / "higher-order.json"
    )
    parser.add_argument("--replicates", type=int, help="smoke run with N replicates")
    args = parser.parse_args()

    import crn_craft as crn
    import numpy as np
    from crn_craft import _native, rules
    from references import power_sums

    manifest = json.loads(MANIFEST.read_text())
    floor = manifest["state_probability_floor"]
    report = {
        "schema": "crn-craft.merge-higher-order-gate/0.1",
        "declared_in": "validation/merge/RESULTS.md",
        "manifest_sha256": sha256(MANIFEST.read_bytes()),
        "smoke": args.replicates is not None,
        "engine": {
            "crn_craft": crn.__version__,
            "native_sha256": sha256(Path(_native.__file__).read_bytes()),
            "numpy": np.__version__,
        },
        "models": [],
    }
    passed = True
    every_row = []
    for row in load_manifest(path=MANIFEST):
        if args.replicates:
            row = {
                **row,
                "replicates": args.replicates,
                "chunk": min(row["chunk"], args.replicates),
            }
        name = row["name"]
        entry = {"name": name, "path": row["path"], "sha256": row["sha256"]}
        network = (ROOT / row["path"]).read_bytes()
        try:
            converted = crn.convert(network)
        except crn.CrnError as error:
            entry["status"] = "refused"
            entry["diagnostics"] = error.diagnostics
            report["models"].append(entry)
            passed = False
            print(f"{name}: refused {error.diagnostics}", file=sys.stderr, flush=True)
            continue
        document = json.loads(converted)
        law = json.loads((REFERENCE / f"{name}.fsp-moments.json").read_text())
        fsp = json.loads(
            gzip.decompress((REFERENCE / f"{name}.fsp.json.gz").read_bytes())
        )
        species = law["species"]
        times = row["report_times"]
        n = row["replicates"]
        entry.update(
            {
                "status": "compared",
                "converted_sha256": sha256(converted),
                "reactants": [rule["reactants"] for rule in document["rules"]],
                "methods": [],
            }
        )
        for method in row["methods"]:
            answer = run_method(
                rules, np, power_sums, converted, document, row, method, species
            )
            stats = statistics_from_power_sums(answer["power_sums"])
            rows = compare_to_law(stats, law, n, times, species)
            rows += state_rows(fsp, answer["histograms"], times, n, floor)
            summary = summarize(rows)
            passed &= summary["failures"] == 0
            every_row += rows
            entry["methods"].append(
                {
                    "method": method,
                    "request_schema": answer["request_schema"],
                    "sampler": answer["sampler"],
                    "replicates": n,
                    "events": answer["events"],
                    "summary": summary,
                    "rows": rows,
                }
            )
            print(f"{name} {method['kind']}: {summary}", file=sys.stderr, flush=True)
        report["models"].append(entry)
    report["summary"] = {**summarize(every_row), "passed": passed}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(dump(report))
    print("| Model | Method | Comparisons | Failures | max \\|z\\| |")
    print("| --- | --- | ---: | ---: | ---: |")
    for model in report["models"]:
        if model["status"] == "refused":
            codes = sorted({d["code"] for d in model["diagnostics"]})
            print(f"| {model['name']} | refused ({', '.join(codes)}) | | | |")
            continue
        for m in model["methods"]:
            s = m["summary"]
            z = "–" if s["max_z"] is None else f"{s['max_z']:.2f}"
            print(
                f"| {model['name']} | {m['method']['kind']} | {s['comparisons']} "
                f"| {s['failures']} | {z} |"
            )
    print(json.dumps(report["summary"], indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
