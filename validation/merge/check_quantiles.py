"""Quantile statistics of the 0.7.0 engine merge (task M4e): exact-law gate.

Declared in `validation/merge/RESULTS.md` ("M4e (quantile statistics):
declaration") before the first run. The immigration-death rules model
(k_in = 10, k_out = 1, X(0) = 0), with only its schema string rewritten to
`crn-craft.model/1.0`, runs 100,000 replicates to t = 10 with a
`crn-craft.run/1.0` request asking for the quantiles 0.025, 0.5 and 0.975,
under `static_direct`, `tree` and `next_reaction`, through the native CLI
(`crn-craft rules MODEL REQUEST`). The exact law at t is Poisson with mean
10 (1 - e^-t).

Gate: at t = 10 each reported quantile lies within one count of the Poisson
quantile Q(p), the smallest k with F(k) >= p (9 comparisons), and at t = 0
every quantile is exactly 0. Contract checks: the request and the results
validate against the committed run-1.0 and result-1.0 schemas, the result
lists the quantiles in request order shaped like `statistics.mean`, and the
same request as `crn-craft.rules-run/0.6` is refused with `CRN-R020` at
`/schema`. Every row is written to the report, failed rows included; the
report holds no timings.

    flock /home/claude/craft5/cpu.lock python validation/merge/check_quantiles.py \\
        --binary target/release/crn-craft --output quantiles.json

Needs Python with `jsonschema` (no NumPy). Exit status 0 when every check
passes, 1 otherwise.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "models/rules/immigration-death.rules.json"
K_IN, K_OUT = 10.0, 1.0
END = 10.0
REPLICATES = 100_000
SEED = 20261111
PROBABILITIES = [0.025, 0.5, 0.975]
METHODS = [{"kind": "static_direct"}, {"kind": "tree"}, {"kind": "next_reaction"}]


def model_1_0():
    """The model bytes with only the schema string rewritten."""
    data = MODEL.read_bytes()
    old = b'"crn-craft.rules-model/0.1"'
    if data.count(old) != 1:
        raise SystemExit("check_quantiles: the model's schema string is not unique")
    return data.replace(old, b'"crn-craft.model/1.0"')


def poisson_quantile(mean, p):
    """The smallest k with P(X <= k) >= p for X ~ Poisson(mean)."""
    k, pmf = 0, math.exp(-mean)
    total = pmf
    while total < p:
        k += 1
        pmf *= mean / k
        total += pmf
    return k


def request(model, schema, method):
    return {
        "schema": schema,
        "model": {
            "path": "model.json",
            "sha256": hashlib.sha256(model).hexdigest(),
        },
        "time": {"end": END, "report_times": [0.0, END]},
        "ensemble": {"replicates": REPLICATES, "root_seed": SEED},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 0,
            "keep_instances": 0,
            "quantiles": PROBABILITIES,
        },
        "budget": {
            "max_events_per_replicate": 1_000_000,
            "max_sampling_steps_per_replicate": 10_000_000,
            "max_instances": 4,
            "threads": 2,
            "max_output_bytes": 1 << 30,
        },
        "method": method,
    }


def run(binary, folder, model, document):
    (folder / "model.json").write_bytes(model)
    path = folder / "request.json"
    path.write_text(json.dumps(document))
    done = subprocess.run(
        [str(binary), "rules", str(folder / "model.json"), str(path)],
        capture_output=True,
        check=False,
    )
    return done.returncode, json.loads(done.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    schemas = {
        kind: Draft7Validator(
            json.loads((ROOT / f"schemas/{kind}.schema.json").read_text())
        )
        for kind in ("run-1.0", "result-1.0")
    }
    model = model_1_0()
    mean = K_IN / K_OUT * (1.0 - math.exp(-K_OUT * END))
    expected = {p: poisson_quantile(mean, p) for p in PROBABILITIES}
    rows, checks = [], []

    def check(name, passed, **detail):
        checks.append({"check": name, "passed": bool(passed), **detail})

    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        for method in METHODS:
            label = method["kind"]
            document = request(model, "crn-craft.run/1.0", method)
            check(f"{label}: request valid", schemas["run-1.0"].is_valid(document))
            code, result = run(args.binary, folder, model, document)
            check(f"{label}: exit status 0", code == 0, status=code)
            check(f"{label}: result valid", schemas["result-1.0"].is_valid(result))
            stats = result["statistics"]
            quantiles = stats["quantiles"]
            check(
                f"{label}: quantiles in request order, shaped like mean",
                [q["p"] for q in quantiles] == PROBABILITIES
                and all(
                    len(q["values"]) == len(stats["mean"])
                    and all(
                        len(a) == len(b)
                        for a, b in zip(q["values"], stats["mean"], strict=True)
                    )
                    for q in quantiles
                ),
            )
            column = result["observables"].index("X")
            check(
                f"{label}: quantiles at t = 0 are 0",
                all(q["values"][0][column] == 0.0 for q in quantiles),
            )
            for q in quantiles:
                value = q["values"][1][column]
                target = expected[q["p"]]
                row = {
                    "method": label,
                    "p": q["p"],
                    "time": END,
                    "quantile": value,
                    "poisson_quantile": target,
                    "difference": value - target,
                    "exact": value == target,
                    "passed": abs(value - target) <= 1.0,
                }
                rows.append(row)
            x_mean = stats["mean"][1][column]
            se = stats["mean_standard_error"][1][column]
            check(
                f"{label}: completed, mean informative",
                result["termination"] == "completed",
                mean=x_mean,
                law_mean=mean,
                z=(x_mean - mean) / se,
            )
        # The field is refused before run/1.0.
        old = request(model, "crn-craft.rules-run/0.6", METHODS[0])
        old_model = MODEL.read_bytes()
        old["model"]["sha256"] = hashlib.sha256(old_model).hexdigest()
        code, refusal = run(args.binary, folder, old_model, old)
        # The CLI prints {"success": false, "diagnostics": [...]}.
        first = (refusal.get("diagnostics") or [{}])[0]
        check(
            "rules-run/0.6 refuses outputs.quantiles",
            code != 0
            and refusal.get("success") is False
            and (first.get("code"), first.get("path")) == ("CRN-R020", "/schema"),
            status=code,
        )
    failures = [r for r in rows if not r["passed"]]
    report = {
        "schema": "crn-craft.merge-quantile-check/0.1",
        "model": "models/rules/immigration-death.rules.json as crn-craft.model/1.0",
        "model_sha256": hashlib.sha256(model).hexdigest(),
        "replicates": REPLICATES,
        "root_seed": SEED,
        "law_mean": mean,
        "poisson_quantiles": {str(p): q for p, q in expected.items()},
        "comparisons": len(rows),
        "failures": len(failures),
        "rows": rows,
        "checks": checks,
        "passed": not failures and all(c["passed"] for c in checks),
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "comparisons": len(rows),
                "failures": len(failures),
                "checks_failed": sum(not c["passed"] for c in checks),
            }
        )
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
