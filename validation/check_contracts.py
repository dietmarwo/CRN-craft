"""Independent JSON Schema and subprocess acceptance of the unified command
line (describe, check, run, verify, analyze, convert, export; rules and
rules-check are aliases). Safe under python -O.

The network engine's commands and contracts were removed in task M8b; rows
named `removed_*` check that they are gone. Network models
(`crn-craft.model/0.1`, the originals in models/legacy) are still read through
the converter; models/examples holds their converted `crn-craft.model/1.0`
files.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = ["describe", "check", "run", "verify", "analyze", "convert", "export"]
# The network engine's commands, removed in M8b.
REMOVED = ["simulate", "simulate-cached", "ssa", "tau", "fsp"]
# The 1.0 contracts and current reports, then the older contracts still read.
SCHEMAS = [
    "model-1.0",
    "run-1.0",
    "result-1.0",
    "check-1.0",
    "fsp-result-1.0",
    "analysis-1.0",
    "archive-1.0",
    "archive-check-1.0",
    "export",
    "state",
]
STILL_READ = [
    "rules-model",
    "rules-run",
    "rules-result",
    "rules-check",
    "model",
]
# The schemas of the network engine's contracts, removed in M8b: no
# `describe KIND` kind and no file under schemas/.
REMOVED_SCHEMAS = [
    "run",
    "result",
    "result-v1",
    "check",
    "fsp-result",
    "tau-request",
    "tau-result",
    "ssa-run",
    "ssa-result",
    "archive",
    "archive-check",
    "analysis",
    "analysis-v1",
]


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

    def call(*command):
        result = subprocess.run(
            [str(binary), *map(str, command)],
            cwd=ROOT,
            capture_output=True,
            timeout=30,
            check=False,
        )
        require(not result.stderr, f"unexpected stderr: {result.stderr!r}")
        return result.returncode, result.stdout

    def invoke(*command):
        code, stdout = call(*command)
        return code, json.loads(stdout)

    def record(name, condition):
        rows.append({"check": name, "passed": bool(condition)})

    def coded(report, code, path=None):
        return any(
            d["code"] == code and (path is None or d["path"] == path)
            for d in report.get("diagnostics", [])
        )

    validators = {}
    for kind in SCHEMAS + STILL_READ:
        code, schema = invoke("describe", kind)
        Draft7Validator.check_schema(schema)
        validators[kind] = Draft7Validator(schema)
        snapshot = json.loads((ROOT / f"schemas/{kind}.schema.json").read_text())
        record(f"schema_snapshot_{kind}", code == 0 and snapshot == schema)
    for kind in REMOVED_SCHEMAS:
        code, report = invoke("describe", kind)
        record(
            f"removed_schema_{kind}",
            code == 1
            and coded(report, "CRN-E031")
            and not (ROOT / f"schemas/{kind}.schema.json").exists(),
        )

    # Discovery, help and version list the unified command set only.
    code, capabilities = invoke("describe")
    record(
        "describe_1_0",
        code == 0
        and capabilities["schema"] == "crn-craft.describe/1.0"
        and capabilities["crate"] == "crn-craft"
        and capabilities["commands"] == COMMANDS
        and capabilities["aliases"] == {"rules": "run", "rules-check": "check"}
        and capabilities["schemas"] == SCHEMAS
        and capabilities["legacy"]["schemas"] == STILL_READ
        and capabilities["contracts"]["model"] == "crn-craft.model/1.0"
        and capabilities["contracts"]["request"] == "crn-craft.run/1.0"
        and capabilities["contracts"]["result"] == "crn-craft.result/1.0"
        and capabilities["legacy"]["requests"][-1] == "crn-craft.rules-run/0.6"
        and "crn-craft.model/0.1" in capabilities["legacy"]["models"]
        and capabilities["legacy"]["rules"]["crate"] == "crn-craft",
    )
    record(
        "describe_lists_no_removed_command_or_contract",
        not set(REMOVED) & set(capabilities["commands"])
        and not set(REMOVED_SCHEMAS)
        & set(capabilities["schemas"] + capabilities["legacy"]["schemas"])
        and "archives" not in capabilities["legacy"]
        and "crn-craft.archive/0.1" not in json.dumps(capabilities),
    )
    code, help_text = call("--help")
    listed = [
        line.split()[1]
        for line in help_text.decode().splitlines()
        if line.startswith("  crn-craft ")
    ]
    record(
        "help_lists_the_unified_commands",
        code == 0
        and listed == COMMANDS + ["--version"]
        and not any(f"crn-craft {name} " in help_text.decode() for name in REMOVED)
        and b"archive/1.0 and 0.1" not in help_text,
    )
    code, version = invoke("--version")
    record(
        "version_and_milestone",
        code == 0 and set(version) == {"version", "milestone"},
    )

    # Network examples: check (of the conversion), convert, analyze, export.
    examples = ROOT / "models/legacy"
    converted_examples = ROOT / "models/examples"
    for path in sorted(examples.glob("*.model.json")):
        validators["model"].validate(json.loads(path.read_text()))
        code, report = invoke("check", path)
        validators["check-1.0"].validate(report)
        record(
            f"model_{path.name}",
            code == 0
            and report["valid"]
            and report["schema"] == "crn-craft.check/1.0"
            and report["model_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        code, converted = call("convert", path)
        document = json.loads(converted)
        validators["model-1.0"].validate(document)
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "converted.json"
            target.write_bytes(converted)
            check_code, checked = invoke("check", target)
            alias_code, alias = invoke("rules-check", target)
        record(
            f"convert_{path.name}",
            code == 0
            and document["schema"] == "crn-craft.model/1.0"
            and (check_code, alias_code) == (0, 0)
            and checked == alias
            and checked["schema"] == "crn-craft.check/1.0",
        )
        # models/examples holds exactly the converted bytes.
        record(
            f"converted_example_{path.name}",
            (converted_examples / path.name).read_bytes() == converted,
        )
        code, analysis = invoke("analyze", path)
        validators["analysis-1.0"].validate(analysis)
        record(
            f"analysis_{path.name}",
            code == 0 and analysis["success"] and not analysis["simulation_available"],
        )
        # analysis/1.0 replaces analysis/0.2 and the legacy analysis/0.1 (M5b).
        code, legacy = invoke("analyze-v1", path)
        record(
            f"legacy_analysis_removed_{path.name}",
            code == 1
            and legacy["diagnostics"][0]["code"] == "CRN-E031"
            and analysis["schema"] == "crn-craft.analysis/1.0",
        )
        code, exported = invoke("export", path, "gillespy2")
        validators["export"].validate(exported)
        # Every example exports since M5f: the Hill and Michaelis-Menten laws
        # of michaelis-menten and repressilator as custom propensities.
        custom = path.name in (
            "michaelis-menten.model.json",
            "repressilator.model.json",
        )
        record(
            f"export_{path.name}",
            code == 0
            and exported["success"]
            and exported["artifact"] is not None
            and ("propensity_function=" in exported["artifact"]["source"]) == custom,
        )
    record(
        "examples_are_the_converted_legacy_models",
        sorted(p.name for p in converted_examples.glob("*.model.json"))
        == sorted(p.name for p in examples.glob("*.model.json")),
    )
    # The network run/0.1 requests of the examples are no request any more:
    # `check` refuses a request for a network model (CRN-E033), and the
    # removed commands are usage errors.
    for path in sorted(examples.glob("*.run.json")):
        run = json.loads(path.read_text())
        model = path.parent / run["model"]["path"]
        code, report = invoke("check", model, path)
        validators["check-1.0"].validate(report)
        record(
            f"removed_network_run_{path.name}",
            code == 1
            and not report["valid"]
            and coded(report, "CRN-E033", "/schema")
            and run["model"]["sha256"] == hashlib.sha256(model.read_bytes()).hexdigest(),
        )
        with tempfile.TemporaryDirectory() as folder:
            outcomes = [
                invoke(command, model, path, Path(folder) / "out")
                if command.startswith("simulate")
                else invoke(command, model, path)
                for command in REMOVED
            ]
            created = (Path(folder) / "out").exists()
        record(
            f"removed_commands_{path.name}",
            not created
            and all(
                code == 1 and coded(report, "CRN-E031") for code, report in outcomes
            ),
        )

    for path in sorted((ROOT / "models/rules").glob("*.rules.json")):
        validators["rules-model"].validate(json.loads(path.read_text()))
        code, report = invoke("check", path)
        validators["rules-check"].validate(report)
        alias = invoke("rules-check", path)
        record(
            f"rules_{path.name}",
            code == 0
            and report["valid"]
            and report["schema"] == "crn-craft.rules-check/0.1"
            and alias == (code, report)
            and report["model_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        # Analysis reads rules models. Models with ATP, continuous variables,
        # random product attributes or products that create attribute tuples
        # during a run have no flat view (CRN-R017).
        code, analysis = invoke("analyze", path)
        validators["analysis-1.0"].validate(analysis)
        closed = path.name in (
            "immigration-death.rules.json",
            "immigration-dimerization.rules.json",
            "mtdna-three-class-cme.rules.json",
            "mtdna-three-class.rules.json",
        )
        record(
            f"rules_analysis_{path.name}",
            (code == 0 and analysis["success"])
            if closed
            else (
                code == 1
                and analysis["analysis"] is None
                and analysis["diagnostics"][0]["code"] == "CRN-R017"
            ),
        )

    # rules-model/0.3 through `run` (and its alias `rules`): a hybrid run of
    # the continuous ATP model (rules-run/0.5, rules-result/0.3), its refusal
    # by Static Direct, and a whole-rule propensity under Static Direct
    # (rules-result/0.2, and rules-result/0.3 with the Ziggurat default of
    # rules-run/0.6). `check MODEL REQUEST` reports what `run` refuses.
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
            code, text = call("run", model, path)
            alias = call("rules", model, path)
            report = json.loads(text)
            check_code, checked = invoke("check", model, path)
            validators["check-1.0"].validate(checked)
            if name == "rules_run_static_refuses_continuous":
                record(
                    name,
                    code == 1
                    and alias == (code, text)
                    and coded(report, "CRN-R020", "/method")
                    and check_code == 1
                    and checked["diagnostics"] == report["diagnostics"],
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
                and alias == (code, text)
                and check_code == 0
                and checked["valid"]
                and checked["run_sha256"]
                == hashlib.sha256(path.read_bytes()).hexdigest()
                and report["termination"] == "completed"
                and report["schema"]
                == (
                    "crn-craft.rules-result/0.3"
                    if hybrid or defaults
                    else "crn-craft.rules-result/0.2"
                )
                and shape,
            )

        # crn-craft.run/1.0 with method fsp: crn-craft.fsp-result/1.0 from
        # `run`, and a missing bound refused at its label by run and check.
        fsp = {
            "schema": "crn-craft.run/1.0",
            "model": {
                "path": propensity.name,
                "sha256": hashlib.sha256(propensity.read_bytes()).hexdigest(),
            },
            "time": {"end": 5.0, "report_times": [0.0, 5.0]},
            "method": {
                "kind": "fsp",
                "bounds": {"A": 40},
                "max_states": 100,
                "tolerance": 1e-12,
            },
            "budget": {"max_output_bytes": 1048576},
        }
        validators["run-1.0"].validate(fsp)
        path = temp / "fsp.run.json"
        path.write_text(json.dumps(fsp))
        code, report = invoke("run", propensity, path)
        validators["fsp-result-1.0"].validate(report)
        check_code, checked = invoke("check", propensity, path)
        record(
            "run_1_0_fsp",
            code == 0
            and check_code == 0
            and checked["valid"]
            and report["schema"] == "crn-craft.fsp-result/1.0"
            and report["species"] == ["A"]
            and len(report["snapshots"]) == 2,
        )
        fsp["method"]["bounds"] = {}
        path.write_text(json.dumps(fsp))
        code, report = invoke("run", propensity, path)
        check_code, checked = invoke("check", propensity, path)
        record(
            "run_1_0_fsp_missing_bound",
            code == 1
            and coded(report, "CRN-R020", "/method/bounds/A")
            and check_code == 1
            and checked["diagnostics"] == report["diagnostics"],
        )

        # Strict request handling of `run` and `check`: unknown fields,
        # duplicate keys, size limits; every refusal exits 1 with JSON.
        base = request(propensity, "crn-craft.rules-run/0.6", {"kind": "static_direct"})
        text = json.dumps(base)
        for name, raw, expected in [
            ("unknown_field", json.dumps({**base, "typo": 1}).encode(), "CRN-E003"),
            (
                "duplicate_key",
                text.replace('"time":', '"ensemble": {}, "time":', 1).encode(),
                "CRN-E002",
            ),
            ("trailing_document", (text + "{}").encode(), "CRN-E001"),
            ("oversized", b" " * (16 * 1024 * 1024 + 1), "CRN-E004"),
            ("deep", b"[" * 1000 + b"0" + b"]" * 1000, "CRN-E004"),
        ]:
            path = temp / f"{name}.run.json"
            path.write_bytes(raw)
            code, report = invoke("run", propensity, path)
            check_code, checked = invoke("check", propensity, path)
            record(
                f"run_rejects_{name}",
                code == 1
                and coded(report, expected)
                and check_code == 1
                and coded(checked, expected),
            )
        # A network model where a rules or 1.0 model is read: CRN-E033.
        network = temp / "network.json"
        network.write_bytes((examples / "immigration-death.model.json").read_bytes())
        document = request(network, "crn-craft.run/1.0", {"kind": "static_direct"})
        path = temp / "network.run.json"
        path.write_text(json.dumps(document))
        code, report = invoke("run", network, path)
        check_code, checked = invoke("check", network, path)
        record(
            "run_refuses_network_models",
            code == 1
            and coded(report, "CRN-E033", "/schema")
            and check_code == 1
            and coded(checked, "CRN-E033", "/schema"),
        )

    invalid = ROOT / "models/invalid"
    for row in json.loads((invalid / "expected.json").read_text()):
        # A network model alone is checked by converting it (check/1.0).
        code, report = invoke("check", invalid / row["model"])
        validators["check-1.0"].validate(report)
        record(
            row["name"],
            code == 1
            and not report["valid"]
            and coded(report, row["code"], row["path"]),
        )

    # The four invalid network run/0.1 requests of 0.6.0's table (removed with
    # that contract in M8b), written as crn-craft.run/1.0 requests for the
    # converted model: refused by `check` and `run` with the unified code, at
    # the same place.
    with tempfile.TemporaryDirectory(prefix="crn-invalid-run-") as folder:
        temp = Path(folder)
        code, converted = call("convert", invalid / "model.json")
        require(code == 0, "convert of the invalid-table model")
        (temp / "model.json").write_bytes(converted)
        digest = hashlib.sha256(converted).hexdigest()
        base = {
            "schema": "crn-craft.run/1.0",
            "model": {"path": "model.json", "sha256": digest},
            "method": {"kind": "static_direct"},
            "time": {"end": 10.0, "report_times": [0.0, 1.0, 2.0, 5.0, 10.0]},
            "ensemble": {"replicates": 64, "root_seed": 42},
            "outputs": {
                "statistics": ["mean", "variance"],
                "keep_paths": 2,
                "keep_instances": 0,
            },
            "budget": {
                "max_events_per_replicate": 1000000,
                "max_sampling_steps_per_replicate": 1000000,
                "max_instances": 100,
                "wall_time_s": 60.0,
                "threads": 1,
                "max_output_bytes": 1048576,
            },
        }
        fsp = {
            "schema": "crn-craft.run/1.0",
            "model": {"path": "model.json", "sha256": digest},
            "time": {"end": 10.0, "report_times": [0.0, 1.0, 2.0, 5.0, 10.0]},
            "method": {
                "kind": "fsp",
                "bounds": {},
                "max_states": 102,
                "tolerance": 1e-10,
            },
            "budget": {"wall_time_s": 60.0, "max_output_bytes": 1048576},
        }
        for name, document, expected_path in [
            (
                "hash-mismatch",
                {**base, "model": {"path": "model.json", "sha256": "0" * 64}},
                "/model/sha256",
            ),
            ("missing-fsp-bound", fsp, "/method/bounds/A"),
            (
                "duplicate-report-time",
                {
                    **base,
                    "time": {"end": 10.0, "report_times": [0.0, 1.0, 1.0, 5.0, 10.0]},
                },
                "/time/report_times/2",
            ),
            (
                "parent-path",
                {**base, "model": {"path": "../model.json", "sha256": digest}},
                "/model/path",
            ),
        ]:
            path = temp / f"{name}.run.json"
            path.write_text(json.dumps(document))
            check_code, checked = invoke("check", temp / "model.json", path)
            run_code, report = invoke("run", temp / "model.json", path)
            validators["check-1.0"].validate(checked)
            record(
                f"run_1_0_{name}",
                check_code == 1
                and run_code == 1
                and not checked["valid"]
                and coded(checked, "CRN-R020", expected_path)
                and coded(report, "CRN-R020", expected_path),
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
    unified = json.loads(
        subprocess.run(
            [str(binary), "convert", examples / "immigration-death.model.json"],
            capture_output=True,
            check=True,
        ).stdout
    )
    for name, value in [
        ("unknown_top_level", {**unified, "typo": 1}),
        ("wrong_version", {**unified, "schema": "crn-craft.model/99"}),
        ("rules_version", {**unified, "schema": "crn-craft.rules-model/0.4"}),
    ]:
        record(
            f"schema_1_0_rejects_{name}", not validators["model-1.0"].is_valid(value)
        )
    pattern = json.loads(json.dumps(unified))
    pattern["rules"][1]["reactants"][0]["typo"] = 1
    record(
        "schema_1_0_rejects_unknown_pattern_field",
        not validators["model-1.0"].is_valid(pattern)
        and validators["model-1.0"].is_valid(unified),
    )

    with tempfile.TemporaryDirectory(prefix="crn-schema-check-") as folder:
        temp = Path(folder)
        mutated = temp / "model.json"
        for name, raw, expected_code in [
            ("trailing_document", b"{}{}", "CRN-E001"),
            ("non_utf8", b"\xff", "CRN-E001"),
            ("deep_input", b"[" * 1000 + b"0" + b"]" * 1000, "CRN-E004"),
            ("oversized_input", b" " * (16 * 1024 * 1024 + 1), "CRN-E004"),
            (
                "duplicate_key",
                b'{"schema": "crn-craft.model/1.0", "name": "a", "name": "b"}',
                "CRN-E002",
            ),
            (
                "unknown_field",
                json.dumps({**unified, "typo": 1}).encode(),
                "CRN-E003",
            ),
        ]:
            mutated.write_bytes(raw)
            code, report = invoke("check", mutated)
            record(name, code == 1 and coded(report, expected_code))

    for command in (["run"], ["rules"], ["report"], ["export"], ["check"], []):
        code, report = invoke(*command)
        hint = report["diagnostics"][0]["hint"] if report.get("diagnostics") else ""
        record(
            f"usage_error_{'_'.join(command) or 'none'}",
            code == 1
            and report["diagnostics"][0]["code"] == "CRN-E031"
            and not any(name in hint for name in ("simulate", "ssa", "tau")),
        )
    first = invoke("check", examples / "immigration-death.model.json")
    second = invoke("check", examples / "immigration-death.model.json")
    record("deterministic_check_report", first == second)

    summary = {
        "schema": "crn-craft.contract-validation/0.2",
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
