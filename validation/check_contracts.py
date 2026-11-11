"""Independent JSON Schema and subprocess acceptance. Safe under python -O."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    rows = []

    def invoke(*command):
        result = subprocess.run(
            [str(binary), *map(str, command)],
            cwd=ROOT,
            capture_output=True,
            timeout=15,
            check=False,
        )
        require(not result.stderr, f"unexpected stderr: {result.stderr!r}")
        return result.returncode, json.loads(result.stdout)

    def record(name, condition):
        rows.append({"check": name, "passed": bool(condition)})

    validators = {}
    for kind in (
        "model",
        "run",
        "result",
        "result-v1",
        "archive",
        "archive-check",
        "check",
        "state",
        "analysis",
        "analysis-v1",
        "export",
        "fsp-result",
        "tau-request",
        "tau-result",
        "ssa-run",
        "ssa-result",
        "rules-model",
        "rules-check",
        "rules-run",
        "rules-result",
    ):
        code, schema = invoke("describe", kind)
        Draft7Validator.check_schema(schema)
        validators[kind] = Draft7Validator(schema)
        snapshot = json.loads((ROOT / f"schemas/{kind}.schema.json").read_text())
        record(f"schema_snapshot_{kind}", code == 0 and snapshot == schema)

    code, capabilities = invoke("describe")
    record(
        "honest_capabilities",
        code == 0
        and capabilities["simulation_available"] is True
        and capabilities["rate_evaluation_available"] is True
        and capabilities["commands"]
        == [
            "describe",
            "check",
            "analyze",
            "analyze-v1",
            "export",
            "simulate",
            "verify",
            "simulate-cached",
            "fsp",
            "tau",
            "ssa",
            "rules-check",
            "rules",
            "convert",
        ]
        and capabilities["dynamic_species"]["model"] == "crn-craft.rules-model/0.1"
        and capabilities["schema"] == "crn-craft.describe/0.14"
        and capabilities["dynamic_species"]["accepted_requests"][-1]
        == "crn-craft.rules-run/0.6",
    )

    examples = ROOT / "models/examples"
    for path in sorted(examples.glob("*.model.json")):
        validators["model"].validate(json.loads(path.read_text()))
        code, report = invoke("check", path)
        validators["check"].validate(report)
        record(f"model_{path.name}", code == 0 and report["valid"])
        code, analysis = invoke("analyze", path)
        validators["analysis"].validate(analysis)
        record(
            f"analysis_{path.name}",
            code == 0 and analysis["success"] and not analysis["simulation_available"],
        )
        code, legacy = invoke("analyze-v1", path)
        validators["analysis-v1"].validate(legacy)
        record(
            f"legacy_analysis_{path.name}",
            code == 0
            and legacy["success"]
            and legacy["schema"] == "crn-craft.analysis/0.1",
        )
        code, exported = invoke("export", path, "gillespy2")
        validators["export"].validate(exported)
        supported = path.name in (
            "enzyme-mechanism.model.json",
            "immigration-death.model.json",
            "vilar-noise.model.json",
        )
        record(
            f"export_{path.name}",
            (code == 0 and exported["success"] and exported["artifact"] is not None)
            if supported
            else (
                code == 1
                and not exported["success"]
                and exported["artifact"] is None
                and any(d["code"] == "CRN-E050" for d in exported["diagnostics"])
            ),
        )
    for path in sorted(examples.glob("*.run.json")):
        run = json.loads(path.read_text())
        validators["run"].validate(run)
        model = path.parent / run["model"]["path"]
        expected = hashlib.sha256(model.read_bytes()).hexdigest()
        code, report = invoke("check", model, path)
        validators["check"].validate(report)
        record(
            f"run_{path.name}",
            code == 0
            and report["valid"]
            and not report["execution_available"]
            and run["model"]["sha256"] == expected == report["model_sha256"]
            and report["run_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(),
        )

    for path in sorted((ROOT / "models/rules").glob("*.rules.json")):
        validators["rules-model"].validate(json.loads(path.read_text()))
        code, report = invoke("rules-check", path)
        validators["rules-check"].validate(report)
        record(
            f"rules_{path.name}",
            code == 0
            and report["valid"]
            and report["model_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(),
        )

    # rules-model/0.3 through the CLI: a hybrid run of the continuous ATP
    # model (rules-run/0.5, rules-result/0.3), its refusal by Static Direct,
    # and a whole-rule propensity under Static Direct (rules-result/0.2, and
    # rules-result/0.3 with the Ziggurat default of rules-run/0.6).
    with tempfile.TemporaryDirectory(prefix="crn-rules-03-") as folder:
        temp = Path(folder)
        atp = temp / "atp.rules.json"
        atp.write_bytes(
            (
                ROOT / "models/rules/mtdna-three-class-atp-continuous.rules.json"
            ).read_bytes()
        )
        propensity = temp / "propensity.rules.json"
        propensity.write_text(
            json.dumps(
                {
                    "schema": "crn-craft.rules-model/0.3",
                    "name": "birth_death",
                    "parameters": {"k": 5.0, "g": 1.0},
                    "types": [{"name": "A"}],
                    "initial": [],
                    "observables": {"A": {"type": "A"}},
                    "rules": [
                        {
                            "name": "birth",
                            "reactants": [],
                            "propensity": "k",
                            "products": [{"type": "A"}],
                        },
                        {
                            "name": "death",
                            "reactants": [{"label": "x", "type": "A"}],
                            "propensity": "g * A",
                            "products": [],
                        },
                    ],
                }
            )
        )

        def request(model, schema, method):
            return {
                "schema": schema,
                "model": {
                    "path": model.name,
                    "sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
                },
                "time": {"end": 50.0, "report_times": [0.0, 50.0]},
                "ensemble": {"replicates": 3, "root_seed": 9},
                "outputs": {
                    "statistics": ["mean"],
                    "keep_paths": 1,
                    "keep_instances": 0,
                },
                "budget": {
                    "max_events_per_replicate": 1000000,
                    "max_sampling_steps_per_replicate": 10000000,
                    "max_instances": 10,
                    "threads": 1,
                    "max_output_bytes": 1048576,
                },
                "method": method,
            }

        cases = [
            (
                "rules_run_hybrid_continuous",
                atp,
                request(
                    atp,
                    "crn-craft.rules-run/0.5",
                    {"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9},
                ),
            ),
            (
                "rules_run_static_refuses_continuous",
                atp,
                request(atp, "crn-craft.rules-run/0.2", {"kind": "static_direct"}),
            ),
            (
                "rules_run_static_propensity",
                propensity,
                request(
                    propensity, "crn-craft.rules-run/0.2", {"kind": "static_direct"}
                ),
            ),
            (
                "rules_run_static_defaults",
                propensity,
                request(
                    propensity, "crn-craft.rules-run/0.6", {"kind": "static_direct"}
                ),
            ),
        ]
        for name, model, document in cases:
            validators["rules-model"].validate(json.loads(model.read_text()))
            validators["rules-run"].validate(document)
            path = temp / f"{name}.run.json"
            path.write_text(json.dumps(document))
            code, report = invoke("rules", model, path)
            if name == "rules_run_static_refuses_continuous":
                record(
                    name,
                    code == 1
                    and report["diagnostics"][0]["code"] == "CRN-R020"
                    and report["diagnostics"][0]["path"] == "/method",
                )
                continue
            validators["rules-result"].validate(report)
            hybrid = name == "rules_run_hybrid_continuous"
            defaults = name == "rules_run_static_defaults"
            if hybrid:
                shape = (
                    report["method"] == document["method"]
                    and report["continuous"] == ["ATP"]
                    and report["observables"][-1] == "ATP"
                )
            elif defaults:
                # The result echoes the request; the identity names the waits.
                shape = (
                    report["method"] == document["method"]
                    and "continuous" not in report
                    and report["sampler"]
                    == "chacha8-u53-ziggurat/rules-static-direct/1"
                )
            else:
                shape = "method" not in report and "continuous" not in report
            record(
                name,
                code == 0
                and report["termination"] == "completed"
                and report["schema"]
                == (
                    "crn-craft.rules-result/0.3"
                    if hybrid or defaults
                    else "crn-craft.rules-result/0.2"
                )
                and shape,
            )

    invalid = ROOT / "models/invalid"
    for row in json.loads((invalid / "expected.json").read_text()):
        command = ["check", invalid / row["model"]]
        if row["run"] is not None:
            command.append(invalid / row["run"])
        code, report = invoke(*command)
        validators["check"].validate(report)
        record(
            row["name"],
            code == 1
            and not report["valid"]
            and any(
                d["code"] == row["code"] and d["path"] == row["path"]
                for d in report["diagnostics"]
            ),
        )

    # Schema validation is structural; semantic cross-references are Rust checks.
    original = json.loads((examples / "immigration-death.model.json").read_text())
    for name, value in [
        ("unknown_top_level", {**original, "typo": 1}),
        ("wrong_version", {**original, "schema": "crn-craft.model/99"}),
        ("zero_omega", {**original, "omega": 0}),
    ]:
        record(f"schema_rejects_{name}", not validators["model"].is_valid(value))
    for path in ("species", "reactions", "rate"):
        value = json.loads(json.dumps(original))
        target = value["reactions"][0]["rate"] if path == "rate" else value[path][0]
        target["typo"] = 1
        record(
            f"schema_rejects_unknown_{path}", not validators["model"].is_valid(value)
        )

    with tempfile.TemporaryDirectory(prefix="crn-schema-check-") as folder:
        temp = Path(folder)
        mutated = temp / "model.json"
        for name, raw, expected_code in [
            ("trailing_document", b"{}{}", "CRN-E001"),
            ("non_utf8", b"\xff", "CRN-E001"),
            ("deep_input", b"[" * 1000 + b"0" + b"]" * 1000, "CRN-E004"),
            ("oversized_input", b" " * (16 * 1024 * 1024 + 1), "CRN-E004"),
        ]:
            mutated.write_bytes(raw)
            code, report = invoke("check", mutated)
            record(
                name,
                code == 1
                and any(d["code"] == expected_code for d in report["diagnostics"]),
            )

    for command in ("run", "report"):
        code, report = invoke(command)
        record(
            f"unimplemented_{command}_fails",
            code == 1 and report["diagnostics"][0]["code"] == "CRN-E031",
        )
    code, report = invoke("export")
    record(
        "export_missing_arguments_fails",
        code == 1 and report["diagnostics"][0]["code"] == "CRN-E031",
    )
    first = invoke("check", examples / "immigration-death.model.json")
    second = invoke("check", examples / "immigration-death.model.json")
    record("deterministic_check_report", first == second)

    summary = {
        "schema": "crn-craft.contract-validation/0.1",
        "scope": "contracts_only_no_simulation",
        "checks": rows,
        "passed": sum(row["passed"] for row in rows),
        "failed": sum(not row["passed"] for row in rows),
    }
    text = json.dumps(summary, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(text)
    print(text, end="")
    require(summary["failed"] == 0, "contract acceptance failed; inspect the report")


if __name__ == "__main__":
    main()
