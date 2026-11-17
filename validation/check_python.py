"""The unified Python API: binding/CLI equivalence, strict inputs, the 1.0
documents of the top-level names, the Model methods and held-out observable
checks.

Task M6a rewrote this script for `crn_craft.Model`, `ensemble`, `prepare`
and `check`. Task M8a deleted validation/legacy/check_python_network.py, which
tested the removed network functions of 0.6.0 (`run_spec`, `simulate`,
`prepare_parameters` and the rest of the old API): this script covers the
same behavior on the unified API (binding and CLI equivalence for 1, 2 and 8
threads, reusable parameters and generated code, strict inputs, the
held-out observable checks, cancel files and keyboard interrupt). Since task
M8b it also checks that the old API is gone: the private module
`crn_craft._legacy`, the native functions behind it and the network
commands of the Python command line.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np
import crn_craft as crn
from crn_craft import rules
from crn_craft.observables import oscillation, sampled_first_passage

# The names of the unified top level (plan section 3, task M6a).
TOP_LEVEL = {
    "CompiledLoopError",
    "CrnError",
    "IncompleteEnsembleError",
    "IncompleteEnsembleWarning",
    "Model",
    "archive",
    "check",
    "compile_loop",
    "convert",
    "copy",
    "describe",
    "ensemble",
    "load",
    "load_example",
    "load_report",
    "loop_source",
    "mass_action",
    "new",
    "pattern",
    "prepare",
    "read_bngl",
    "read_sbml",
    "verify_archive",
}
# The network functions of 0.6.0 that left the top level.
REMOVED = {
    "run_spec",
    "simulate",
    "ssa_spec",
    "simulate_ssa",
    "prepare_ssa",
    "prepare_parameters",
    "ParameterTemplate",
    "PreparedSimulation",
    "BoundParameterSimulation",
    "mass_action_kernel_source",
    "tau",
    "fsp",
    "analyze",
    "model_bytes",
}
# The native functions and classes only the network API used, removed in M8b.
NATIVE_REMOVED = {
    "BoundParameterRun",
    "ParameterTemplateRun",
    "PreparedRun",
    "PreparedSsaRun",
    "archive_json",
    "check_files_json",
    "check_json",
    "fsp_json",
    "legacy_describe_json",
    "mass_action_kernel_source",
    "prepare_json",
    "prepare_parameters_json",
    "prepare_ssa_json",
    "simulate_cached_json",
    "simulate_json",
    "ssa_files_json",
    "ssa_json",
    "tau_files_json",
    "tau_json",
}
# The network engine's commands, removed from both command lines in M8b.
REMOVED_COMMANDS = ("simulate", "simulate-cached", "ssa", "tau", "fsp")
# Every public name python/crn_craft/rules.py defined or imported in 0.6.0.
RULES_0_6_0 = {
    "CompiledLoop",
    "CompiledLoopError",
    "CrnError",
    "HYBRID_TOLERANCES",
    "Model",
    "Path",
    "PreparedEnsemble",
    "Result",
    "TERMINATIONS",
    "cached_property",
    "compile_loop",
    "describe",
    "ensemble",
    "fingerprint_of",
    "hashlib",
    "json",
    "load",
    "loop_source",
    "model_bytes",
    "np",
    "operator",
    "pd",
    "prepare",
    "xarray",
}
# rules.Model("sir") of the manual as release 0.6.0 wrote it (rules-model/0.3).
SIR_0_6_0 = (
    b'{"schema":"crn-craft.rules-model/0.3","name":"sir","parameters":{"beta":0.0002,'
    b'"gamma":0.01},"types":[{"name":"S"},{"name":"I"},{"name":"R"}],"initial":'
    b'[{"type":"S","count":999},{"type":"I","count":1}],"observables":{"S":{"type":'
    b'"S"},"I":{"type":"I"},"R":{"type":"R"}},"rules":[{"name":"infection",'
    b'"reactants":[{"label":"x","type":"S"},{"label":"y","type":"I"}],"propensity":'
    b'"beta * S * I","products":[{"type":"I"},{"type":"I"}]},{"name":"recovery",'
    b'"reactants":[{"label":"x","type":"I"}],"propensity":"gamma * I","products":'
    b'[{"type":"R"}]},{"name":"dimer","reactants":[{"label":"x","type":"I"},'
    b'{"label":"y","type":"I"}],"rate":"0.5","products":[]}],"stop":"I == 0"}\n'
)


def require(value, message):
    if not value:
        raise RuntimeError(message)


def raises(error, call, *args, **kwargs):
    try:
        call(*args, **kwargs)
    except error as caught:
        return caught
    raise RuntimeError(f"{call.__name__} did not raise {error.__name__}")


def sir(model_class):
    m = model_class("sir")
    m.add_species({"S": 999, "I": 1, "R": 0}).add_parameter(
        {"beta": 2e-4, "gamma": 0.01}
    )
    m.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
    m.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    m.add_reaction(0.5, ["I", "I"], [], name="dimer")
    return m.set_stop("I == 0")


def cli_identity(cli, result):
    """A result printed by `crn-craft run` equals a Python Result."""
    cli = dict(cli)
    # Result.metadata holds the records as columns, not as a list.
    records, cli["replicates"] = cli["replicates"], []
    frame = result.records
    return (
        cli == result.metadata
        and [r["id"] for r in records] == frame.index.tolist()
        and [r["events"] for r in records] == frame["events"].tolist()
        and [r["termination"] for r in records] == frame["termination"].tolist()
        and [r["rng_words"] for r in records] == frame["rng_words"].tolist()
        and [r["rule_events"] for r in records] == result.rule_events.values.tolist()
        and {r["id"]: r["diagnostic"] for r in records if r.get("diagnostic")}
        == result.record_diagnostics
    )


def check_names(rows):
    """The top level exports the unified names; crn_craft.rules its 0.6.0 ones."""
    public = {n for n in dir(crn) if not n.startswith("_")}
    modules = {"builder", "importers", "observables", "report", "rules"}
    require(set(crn.__all__) == TOP_LEVEL, f"__all__: {sorted(crn.__all__)}")
    require(public - modules == TOP_LEVEL, f"top level: {sorted(public - modules)}")
    require(not REMOVED & public, f"removed names exported: {REMOVED & public}")
    native = set(dir(crn._native))
    require(
        not NATIVE_REMOVED & native,
        f"removed native names present: {sorted(NATIVE_REMOVED & native)}",
    )
    try:
        import crn_craft._legacy  # noqa: F401
    except ImportError:
        pass
    else:
        raise RuntimeError("crn_craft._legacy is still importable")
    missing = RULES_0_6_0 - set(dir(rules))
    require(not missing, f"crn_craft.rules lost 0.6.0 names: {sorted(missing)}")
    require(
        rules.CrnError is crn.CrnError
        and rules.load_report is crn.load_report
        and rules.archive is crn.archive
        and rules.verify_archive is crn.verify_archive
        and issubclass(crn.Model, rules.Model)
        and crn.Model is not rules.Model,
        "shared objects of the two namespaces",
    )
    rows.append({"check": "names", "passed": True, "top_level": sorted(TOP_LEVEL)})


def check_documents(rows):
    """The top-level names write the 1.0 contracts; crn_craft.rules the 0.6.0 ones."""
    old, new = sir(rules.Model), sir(crn.Model)
    require(
        rules.model_bytes(old.to_dict()) == SIR_0_6_0,
        "rules.Model writes the 0.6.0 document byte for byte",
    )
    require(
        new.to_dict() == {**old.to_dict(), "schema": "crn-craft.model/1.0"}
        and list(new.to_dict()) == list(old.to_dict()),
        "crn.Model writes the same document as model/1.0",
    )
    settings = {
        "end": 50.0,
        "report_times": [0.0, 25.0, 50.0],
        "replicates": 20,
        "seed": 3,
        "method": "static_direct",
    }
    plain = rules.ensemble(old, **settings)
    model, request = (json.loads(d) for d in plain._documents)
    plain_request = request
    require(
        model["schema"] == "crn-craft.rules-model/0.3"
        and request["schema"] == "crn-craft.rules-run/0.6"
        and "quantiles" not in request["outputs"]
        and "selector" not in request["method"]
        and not plain._unified,
        "crn_craft.rules.ensemble writes the 0.6.0 documents",
    )
    unified = crn.ensemble(old, **settings)
    model_1_0, request_1_0 = unified._documents
    expected = dict(request, schema="crn-craft.run/1.0")
    expected["model"] = {
        "path": request["model"]["path"],
        "sha256": hashlib.sha256(model_1_0).hexdigest(),
    }
    expected["outputs"] = {**request["outputs"], "quantiles": [0.025, 0.5, 0.975]}
    require(
        model_1_0 == plain._documents[0].replace(b"rules-model/0.3", b"model/1.0", 1)
        and json.loads(request_1_0) == expected
        and unified._unified,
        "crn.ensemble writes model/1.0 and run/1.0 with the default quantiles",
    )
    require(
        crn.ensemble(new, **settings)._documents == unified._documents,
        "a crn.Model and a rules.Model give the same top-level documents",
    )
    a, b = plain.simulate(), unified.simulate()
    require(
        b.metadata["schema"] == "crn-craft.result/1.0"
        and list(b.statistics["quantiles"]) == [0.025, 0.5, 0.975]
        and b.statistics["median"] is b.statistics["quantiles"][0.5]
        and b.metadata["sampler"] == a.metadata["sampler"]
        and b.records.equals(a.records)
        and b.statistics["mean"].equals(a.statistics["mean"]),
        "the 1.0 documents change the contracts, not the numbers",
    )
    for none in (None, ()):
        off = crn.ensemble(old, quantiles=none, **settings)
        require(
            "quantiles" not in json.loads(off._documents[1])["outputs"],
            f"quantiles={none!r} switches them off",
        )
    # selector= on both names; through crn_craft.rules it needs the 1.0
    # documents, as quantiles= does.
    for name, call in (("crn", crn.ensemble), ("rules", rules.ensemble)):
        chosen = call(old, selector="composition_rejection", **settings)
        model, request = (json.loads(d) for d in chosen._documents)
        result = chosen.simulate()
        require(
            model["schema"] == "crn-craft.model/1.0"
            and request["schema"] == "crn-craft.run/1.0"
            and request["method"]
            == {"kind": "static_direct", "selector": "composition_rejection"}
            and "rules-static-direct-cr" in result.metadata["sampler"],
            f"{name}.ensemble(selector=)",
        )
    raises(ValueError, crn.ensemble, old, selector="fastest", **settings)
    raises(
        ValueError,
        crn.ensemble,
        old,
        selector="linear",
        **{**settings, "method": "next_reaction"},
    )
    # prepare(): a rules-run request keeps its meaning as run/1.0, the
    # default quantiles are added unless the request states its own.
    for given in (
        plain_request,
        json.dumps(plain_request),
        json.dumps(plain_request).encode(),
    ):
        prepared = crn.prepare(plain._documents[0], given)
        document = json.loads(prepared._documents[1])
        require(
            document["schema"] == "crn-craft.run/1.0"
            and document["method"] == {"kind": "static_direct", "selector": "linear"}
            and document["outputs"]["quantiles"] == [0.025, 0.5, 0.975]
            and json.loads(prepared._documents[0])["schema"] == "crn-craft.model/1.0",
            "crn.prepare rewrites a rules-run request as run/1.0",
        )
    exact = crn.prepare(*unified._documents)
    require(
        exact._documents == unified._documents,
        "crn.prepare keeps the bytes of run/1.0 documents with quantiles",
    )
    carried = crn.prepare(model_1_0, dict(json.loads(request_1_0)), selector="sorting")
    require(
        json.loads(carried._documents[1])["method"]["selector"] == "sorting",
        "crn.prepare(selector=)",
    )
    chosen = rules.prepare(*plain._documents, selector="linear")
    require(
        json.loads(chosen._documents[1])["method"]
        == {"kind": "static_direct", "selector": "linear"}
        and json.loads(chosen._documents[1])["schema"] == "crn-craft.run/1.0",
        "rules.prepare(selector=) writes run/1.0",
    )
    # A network model through the top level is converted.
    network = (ROOT / "models/legacy/immigration-death.model.json").read_bytes()
    converted = crn.ensemble(network, end=1.0, report_times=[0.0, 1.0])
    require(
        converted._documents[0] == crn.convert(network),
        "crn.ensemble converts a network model",
    )
    refused = raises(crn.CrnError, rules.ensemble, network, end=1.0, report_times=[1.0])
    require(
        refused.diagnostics[0]["code"] == "CRN-E033",
        "crn_craft.rules refuses a network model with CRN-E033",
    )
    rows.append({"check": "documents_1_0", "passed": True})


def check_models(rows):
    """add_reaction with repeated names, from_dict, load_example, importers."""
    pair = rules.Model("p").add_species("A", 10).add_reaction(0.1, ["A", "A"], [])
    require(
        pair.to_dict()["rules"][0]["reactants"]
        == [{"label": "x", "type": "A"}, {"label": "y", "type": "A"}]
        and pair.to_dict()["schema"] == "crn-craft.rules-model/0.3",
        "two identical names keep the 0.6.0 pair",
    )
    third = rules.Model("t").add_species("A", 10)
    third.add_reaction(0.1, ["A", "A", "A"], ["B"])
    third.add_reaction("k * A", ["A", "A", "B", "C"], [], name="mixed")
    third.add_parameter("k", 2.0)
    document = third.to_dict()
    require(
        document["schema"] == "crn-craft.model/1.0"
        and document["rules"][0]["reactants"]
        == [{"label": "x", "type": "A", "count": 3}]
        and document["rules"][1]["reactants"]
        == [
            {"label": "x1", "type": "A", "count": 2},
            {"label": "x2", "type": "B"},
            {"label": "x3", "type": "C"},
        ]
        and crn.check(third)["valid"],
        "three or more molecules use count and model/1.0",
    )
    # from_dict keeps a document exactly; a converted model is editable.
    death = crn.load_example("immigration-death")
    packaged = crn.convert(
        (ROOT / "models/legacy/immigration-death.model.json").read_bytes()
    )
    require(
        rules.model_bytes(death.to_dict()) == rules.model_bytes(json.loads(packaged))
        and list(death.to_dict()) == list(json.loads(packaged))
        and death.species == ("A",)
        and death.reaction_names == ("birth", "death"),
        "load_example wraps the converted document",
    )
    for name in (
        "enzyme-mechanism",
        "michaelis-menten",
        "repressilator",
        "vilar-noise",
    ):
        model = crn.load_example(name)
        path = ROOT / f"models/legacy/{name}.model.json"
        require(
            model.to_dict() == json.loads(crn.convert(path.read_bytes())),
            f"load_example({name!r})",
        )
    raises(ValueError, crn.load_example, "missing")
    edited = crn.load_example("immigration-death").set_initial(A=7)
    require(
        edited.to_dict()["initial"] == [{"type": "A", "count": 7}]
        and edited.to_dict()["schema"] == "crn-craft.model/1.0"
        and crn.check(edited)["valid"],
        "an edited converted model writes the builder's document",
    )
    three = json.loads((ROOT / "models/rules/mtdna-three-class.rules.json").read_text())
    fixed = rules.Model.from_dict(three)
    require(
        fixed.to_dict() == three
        and crn.Model.from_dict(three).to_dict()
        == {**three, "schema": "crn-craft.model/1.0"}
        and fixed.species == ("mtDNA",),
        "a document the builder cannot write is kept whole",
    )
    raises(ValueError, fixed.add_species, "B")
    raises(ValueError, fixed.set_initial, mtDNA=3)
    error = raises(crn.CrnError, crn.Model.from_dict, {**three, "rules": []})
    require(error.diagnostics, "from_dict checks the document")
    require(
        crn.Model.from_dict(
            (ROOT / "models/legacy/repressilator.model.json").read_bytes()
        ).to_dict()
        == crn.load_example("repressilator").to_dict(),
        "from_dict converts a network model",
    )
    # Importers return Models.
    with tempfile.TemporaryDirectory() as folder:
        net = Path(folder) / "decay.net"
        net.write_bytes(
            b"begin parameters\n1 k 0.5\nend parameters\nbegin species\n1 A() 10\n"
            b"end species\nbegin reactions\n1 1 0 k\nend reactions\n"
        )
        imported = crn.read_bngl(net, convention="binomial")
        require(
            isinstance(imported, crn.Model)
            and imported.to_dict()["schema"] == "crn-craft.model/1.0"
            and imported.import_report["schema"] == "crn-craft.import/1.0"
            and "model" not in imported.import_report,
            "read_bngl returns a Model with the import provenance",
        )
    rows.append({"check": "models", "passed": True})


def check_methods(rows):
    """Model.fsp, Model.analyze and Model.to_gillespy2."""
    death = crn.load_example("immigration-death")
    result = death.fsp(end=10.0, report_times=[0.0, 5.0, 10.0], bounds={"A": 120})
    times = np.array([0.0, 5.0, 10.0])
    mean = 20.0 * (1.0 - np.exp(-0.2 * times))
    require(
        isinstance(result, rules.FspResult)
        and result.states.shape == (121, 1)
        and list(result.states.columns) == ["A"]
        and result.probabilities.dims == ("time", "state")
        and result.probabilities.shape == (3, 121)
        and np.allclose(result.mean["A"].to_numpy(), mean, rtol=0, atol=1e-9)
        and list(result.truncation.columns)
        == [
            "physical_leakage_estimate",
            "series_tail_upper_bound",
            "discarded_weight_estimate",
            "mass_defect",
            "roundoff_mass_residual",
        ]
        and result.metadata["schema"] == "crn-craft.fsp-result/1.0",
        "Model.fsp: the Poisson mean of immigration-death",
    )
    report = crn.load_example("enzyme-mechanism").analyze()
    require(
        report["schema"] == "crn-craft.analysis/1.0"
        and report["success"]
        and report["analysis"]["species"] == ["S", "E", "ES", "P"],
        "Model.analyze",
    )
    try:
        import gillespy2  # noqa: F401
    except ImportError:
        error = raises(ImportError, crn.load_example("michaelis-menten").to_gillespy2)
        gillespy = f"skipped: {type(error).__name__}"
    else:
        target = crn.load_example("michaelis-menten").to_gillespy2()
        require(
            sorted(target.listOfSpecies) == ["S0", "S1"]
            and len(target.listOfReactions) == 1,
            "Model.to_gillespy2",
        )
        gillespy = "built"
    refused = raises(crn.CrnError, sir(crn.Model).to_gillespy2)
    require(refused.diagnostics[0]["code"] == "CRN-E050", "a stop is not exported")
    rows.append({"check": "model_methods", "passed": True, "to_gillespy2": gillespy})


def check_cli(rows, binary):
    """The Python command line is the native command set."""
    native = [str(binary)]
    python = [sys.executable, "-m", "crn_craft"]

    def run(command, *args):
        done = subprocess.run(
            [*command, *map(str, args)], capture_output=True, check=False, timeout=60
        )
        require(not done.stderr, f"stderr: {done.stderr!r}")
        return done.returncode, done.stdout

    help_lines = run(python, "--help")[1].decode().splitlines()
    commands = [
        line.split()[1] for line in help_lines if line.startswith("  crn-craft ")
    ]
    require(
        commands
        == ["describe", "check", "run", "verify", "analyze", "convert", "export"]
        + ["--version"],
        f"python --help commands: {commands}",
    )
    require(
        json.loads(run(python, "describe")[1])
        == json.loads(run(native, "describe")[1])
        == crn.describe(),
        "describe is the same document",
    )
    death = crn.load_example("immigration-death")
    prepared = crn.ensemble(
        death,
        end=2.0,
        report_times=[0.0, 1.0, 2.0],
        replicates=16,
        seed=4,
        keep_paths=2,
        method="static_direct",
    )
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        model, request = folder / "model.rules.json", folder / "request.json"
        model.write_bytes(prepared._documents[0])
        request.write_bytes(prepared._documents[1])
        for alias in ("run", "rules"):
            code, by_python = run(python, alias, model, request)
            require(code == 0, by_python.decode())
            require(
                by_python == run(native, alias, model, request)[1],
                f"python -m crn_craft {alias} prints the native bytes",
            )
        for alias in ("check", "rules-check"):
            for args in ((model,), (model, request)):
                require(
                    json.loads(run(python, alias, *args)[1])
                    == json.loads(run(native, alias, *args)[1]),
                    f"{alias} reports agree",
                )
        code, converted = run(
            python, "convert", ROOT / "models/legacy/vilar-noise.model.json"
        )
        require(
            code == 0
            and converted
            == run(native, "convert", ROOT / "models/legacy/vilar-noise.model.json")[1],
            "convert prints the native bytes",
        )
        network = ROOT / "models/legacy/immigration-death.model.json"
        legacy_run = ROOT / "models/legacy/immigration-death.run.json"
        for command in REMOVED_COMMANDS:
            for args in ((), (network, legacy_run), (model, request, folder / "out")):
                code, usage = run(python, command, *args)
                require(
                    code == 1
                    and json.loads(usage)["schema"] == "crn-craft.check/1.0"
                    and json.loads(usage)["diagnostics"][0]["code"] == "CRN-E031"
                    and b"simulate" not in usage
                    and json.loads(usage) == json.loads(run(native, command, *args)[1]),
                    f"the removed command {command} is a usage error",
                )
        require(not (folder / "out").exists(), "no archive from a removed command")
    version = json.loads(run(python, "--version")[1])
    require(
        version
        == {
            "version": crn.__version__,
            "interface": "python",
            "milestone": "unified-engine",
        },
        "python --version",
    )
    rows.append({"check": "command_line", "passed": True})


ROOT = Path(__file__).resolve().parents[1]


# Release 0.8.0: the limits of ensemble, the incomplete-ensemble warning and
# the builder for species with attributes.
U64_MAX = 2**64 - 1
# The clone model of a published mtDNA drift study, as its author wrote it by
# hand (model/1.0): one attribute `clone`, new_id() for a new mutant clone.
DEG_RATE = math.log(2) / 240.0


def clone_document(pmut):
    x = {"label": "x", "type": "mtDNA"}
    refill = "r * max(0, K - N) / N"
    return {
        "schema": "crn-craft.model/1.0",
        "name": "mtdna_clones",
        "parameters": {"K": 1000, "delta": DEG_RATE, "pmut": pmut, "r": 10.0},
        "types": [{"name": "mtDNA", "attributes": ["clone"]}],
        "initial": [{"type": "mtDNA", "attributes": {"clone": 0}, "count": 1000}],
        "observables": {
            "N": {"type": "mtDNA"},
            "mutants": {"type": "mtDNA", "where": "x.clone > 0"},
            "types": {"type": "mtDNA", "reduce": "distinct"},
            "mutTypes": {"type": "mtDNA", "where": "x.clone > 0", "reduce": "distinct"},
        },
        "rules": [
            {"name": "degradation", "reactants": [x], "rate": "delta", "products": []},
            {
                "name": "replication",
                "reactants": [x],
                "rate": f"{refill} * (1 - pmut)",
                "products": [
                    {"type": "mtDNA", "copy": "x"},
                    {"type": "mtDNA", "copy": "x"},
                ],
            },
            {
                "name": "mutation",
                "reactants": [x],
                "rate": f"{refill} * pmut",
                "products": [
                    {"type": "mtDNA", "copy": "x"},
                    {"type": "mtDNA", "attributes": {"clone": "new_id()"}},
                ],
            },
        ],
    }


def clone_model(model_class, pmut, degradation):
    """The clone model with the builder; `degradation` is its rate."""
    model = model_class("mtdna_clones")
    model.add_parameter({"K": 1000, "delta": DEG_RATE, "pmut": pmut, "r": 10.0})
    model.add_species("mtDNA", attributes=["clone"])
    model.add_initial("mtDNA", 1000, clone=0)
    model.add_observable("N", "mtDNA")
    model.add_observable("mutants", "mtDNA", where="x.clone > 0")
    model.add_observable("types", "mtDNA", reduce="distinct")
    model.add_observable("mutTypes", "mtDNA", where="x.clone > 0", reduce="distinct")
    refill = "r * max(0, K - N) / N"
    model.add_reaction(degradation, ["mtDNA"], [], name="degradation")
    model.add_reaction(
        crn.mass_action(f"{refill} * (1 - pmut)"),
        ["mtDNA"],
        [crn.copy("x"), crn.copy("x")],
        name="replication",
    )
    model.add_reaction(
        crn.mass_action(f"{refill} * pmut"),
        ["mtDNA"],
        [crn.copy("x"), crn.new("mtDNA", clone="new_id()")],
        name="mutation",
    )
    return model


def budget_of(prepared):
    return json.loads(prepared._documents[1])["budget"]


def check_limits(rows, binary):
    """L1-L3: max_events 10**8 by default, sampling follows it, None is u64::MAX."""
    model = crn.load_example("immigration-death")
    short = {"end": 30.0, "report_times": [0.0, 30.0]}
    for name, call in (("crn", crn.ensemble), ("rules", rules.ensemble)):
        budget = budget_of(call(model, **short))
        require(
            budget["max_events_per_replicate"] == 10**8
            and budget["max_sampling_steps_per_replicate"] == 16 * 10**8,
            f"{name}.ensemble writes the default budgets",
        )
        budget = budget_of(call(model, max_events=1000, **short))
        require(
            budget["max_sampling_steps_per_replicate"] == 16_000,
            f"{name}.ensemble: the sampling budget follows max_events",
        )
        budget = budget_of(call(model, max_events=1000, max_sampling_steps=5, **short))
        require(
            budget["max_events_per_replicate"] == 1000
            and budget["max_sampling_steps_per_replicate"] == 5,
            f"{name}.ensemble: an explicit sampling budget wins",
        )
        budget = budget_of(call(model, max_events=2**63, **short))
        require(
            budget["max_sampling_steps_per_replicate"] == U64_MAX,
            f"{name}.ensemble: the derived sampling budget is capped at 2**64 - 1",
        )
        unlimited = call(model, max_events=None, method="static_direct", **short)
        budget = budget_of(unlimited)
        require(
            budget["max_events_per_replicate"] == U64_MAX
            and budget["max_sampling_steps_per_replicate"] == U64_MAX
            and b'"max_events_per_replicate":18446744073709551615'
            in unlimited._documents[1],
            f"{name}.ensemble(max_events=None) writes 2**64 - 1",
        )
        result = unlimited.simulate()
        require(
            result.metadata["termination"] == "completed"
            and result.request["budget"]["max_events_per_replicate"] == U64_MAX
            and crn.check(*unlimited._documents)["valid"],
            f"{name}: the engine reads the u64 budgets",
        )
        budget = budget_of(
            call(model, max_events=None, max_sampling_steps=10**6, **short)
        )
        require(
            budget["max_events_per_replicate"] == U64_MAX
            and budget["max_sampling_steps_per_replicate"] == 10**6,
            f"{name}.ensemble: no event limit with an explicit sampling budget",
        )
        for keyword in ("max_events", "max_sampling_steps"):
            for bad in (-1, 2**64):
                raises(ValueError, call, model, **{keyword: bad}, **short)
            raises(TypeError, call, model, **{keyword: 1.5}, **short)
        raises(TypeError, call, model, strict="yes", **short)
        for keyword in (
            "max_events",
            "max_sampling_steps",
            "max_instances",
            "max_output_bytes",
            "threads",
            "wall_time_s",
        ):
            for flag in (True, np.bool_(False)):
                raises(TypeError, call, model, **{keyword: flag}, **short)
    # The command line runs the u64 request too.
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        unlimited = crn.ensemble(
            model, max_events=None, method="static_direct", replicates=2, **short
        )
        (folder / "model.rules.json").write_bytes(unlimited._documents[0])
        (folder / "run.json").write_bytes(unlimited._documents[1])
        done = subprocess.run(
            [str(binary), "run", "model.rules.json", "run.json"],
            cwd=folder,
            capture_output=True,
            check=False,
        )
        require(
            done.returncode == 0
            and json.loads(done.stdout)["termination"] == "completed",
            "the command line runs a request with max_events=None",
        )
    rows.append({"check": "limits_0_8", "passed": True})


def caught_warnings(call):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = call()
    return result, [w for w in caught if w.category is crn.IncompleteEnsembleWarning]


def check_incomplete(rows):
    """L4: simulate() warns when replicates end early; strict=True raises."""
    model = crn.load_example("immigration-death")
    short = {
        "end": 30.0,
        "report_times": [0.0, 30.0],
        "replicates": 8,
        "seed": 1,
        "method": "static_direct",
    }
    require(
        issubclass(crn.IncompleteEnsembleWarning, UserWarning)
        and issubclass(crn.IncompleteEnsembleError, RuntimeError)
        and rules.IncompleteEnsembleWarning is crn.IncompleteEnsembleWarning
        and rules.IncompleteEnsembleError is crn.IncompleteEnsembleError,
        "the warning and error classes",
    )
    for name, call in (("crn", crn.ensemble), ("rules", rules.ensemble)):
        result, caught = caught_warnings(call(model, max_events=50, **short).simulate)
        message = str(caught[0].message) if caught else ""
        require(
            len(caught) == 1
            and caught[0].filename == __file__
            and message.startswith(
                "8 of 8 replicates ended early: 8 event_limit (raise max_events)"
            )
            and caught[0].message.counts == {"event_limit": 8}
            and result.statistics is None,
            f"{name}: event_limit warns at the caller with the counts: {message}",
        )
        error = raises(
            crn.IncompleteEnsembleError,
            call(model, max_events=50, strict=True, **short).simulate,
        )
        require(
            str(error).startswith(message.split(" Pass strict=True")[0])
            and error.counts == {"event_limit": 8}
            and isinstance(error.result, rules.Result),
            f"{name}: strict=True raises with the same message",
        )
        complete, caught = caught_warnings(call(model, **short).simulate)
        require(
            not caught and complete.metadata["termination"] == "completed",
            f"{name}: a complete run does not warn",
        )
    # Every reason names its keyword.
    _, caught = caught_warnings(
        crn.ensemble(model, max_sampling_steps=40, **short).simulate
    )
    require(
        "8 sampling_limit (raise max_sampling_steps)" in str(caught[0].message),
        "sampling_limit names max_sampling_steps",
    )
    failing = crn.Model("failing").add_species("A", 0)
    failing.add_reaction(1.0, [], ["B"], name="clock")
    failing.add_reaction("5 - B", [], ["A"], name="birth")
    _, caught = caught_warnings(crn.ensemble(failing, **short).simulate)
    require(
        "8 failed (see result.record_diagnostics)" in str(caught[0].message),
        "failed replicates point to record_diagnostics",
    )
    wide = dict(short, report_times=list(range(31)), keep_paths="all")
    dropped, caught = caught_warnings(
        crn.ensemble(model, max_output_bytes=3000, **wide).simulate
    )
    require(
        dropped.metadata["termination"] == "output_limit"
        and str(caught[0].message).startswith("The result exceeded max_output_bytes"),
        "an output limit names max_output_bytes",
    )
    # Stopped and cancelled replicates are not reported.
    stopping = crn.load_example("immigration-death").set_stop("A >= 3")
    stopped, caught = caught_warnings(crn.ensemble(stopping, **short).simulate)
    require(
        not caught and set(stopped.records["termination"]) == {"stopped"},
        "stopped replicates do not warn",
    )
    with tempfile.TemporaryDirectory() as folder:
        flag = Path(folder) / "stop"
        flag.touch()
        run = crn.ensemble(model, strict=True, **short)
        cancelled = run.simulate(cancel_file=flag)
        require(
            cancelled.metadata["termination"] == "cancelled",
            "a cancelled run is not an incomplete-ensemble error",
        )
    # strict survives bind_parameters, and prepare takes it.
    template = crn.ensemble(model, max_events=50, strict=True, **short)
    raises(
        crn.IncompleteEnsembleError,
        template.bind_parameters({"birth": 5.0}).simulate,
    )
    _, caught = caught_warnings(
        crn.ensemble(model, max_events=50, **short)
        .bind_parameters({"birth": 5.0})
        .simulate
    )
    require(len(caught) == 1, "a bound candidate warns too")
    batches = crn.ensemble(model, max_events=50, **short)
    with warnings.catch_warnings(record=True) as caught:
        for _ in range(2):
            batches.simulate()  # one call site, as in a loop over batches
    require(
        sum(w.category is crn.IncompleteEnsembleWarning for w in caught) == 2,
        "two batches from one call site warn twice under the default filters",
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("ignore", crn.IncompleteEnsembleWarning)
        batches.simulate()
    require(not caught, "a filter of the user wins over the package's filter")
    for name, call in (("crn", crn.prepare), ("rules", rules.prepare)):
        documents = rules.ensemble(model, max_events=50, **short)._documents
        raises(crn.IncompleteEnsembleError, call(*documents, strict=True).simulate)
        _, caught = caught_warnings(call(*documents).simulate)
        require(len(caught) == 1, f"{name}.prepare warns")
    rows.append({"check": "incomplete_ensemble_warning", "passed": True})


def check_attributed_builder(rows):
    """B1-B4: species with attributes, patterns, copy/new, mass_action and
    observables in the builder."""
    hand = clone_document(1e-3)
    built = clone_model(crn.Model, 1e-3, crn.mass_action("delta"))
    document = built.to_dict()
    require(
        crn.check(built)["valid"]
        and list(document["observables"]) == ["mtDNA", *hand["observables"]]
        and document["observables"].pop("mtDNA") == {"type": "mtDNA"}
        and document == hand
        and json.loads(built.to_json()) == built.to_dict(),
        "the clone model with the builder is the hand-written document plus "
        "the species count",
    )
    numeric = clone_model(crn.Model, 1e-3, DEG_RATE).to_dict()
    require(
        numeric["rules"][0]["rate"] == repr(DEG_RATE),
        "a numeric rate on a species with attributes is its rate per molecule",
    )
    old = clone_model(rules.Model, 1e-3, crn.mass_action("delta"))
    old.add_observable("largest", "mtDNA", where="x.clone > 0", reduce="max_count")
    require(
        old.to_dict()
        == {**built.to_dict(), "observables": old.to_dict()["observables"]}
        and old.to_dict()["observables"]["largest"]
        == {"type": "mtDNA", "where": "x.clone > 0", "reduce": "max_count"}
        and crn.check(old)["valid"],
        "rules.Model writes an attributed model as model/1.0",
    )
    # Each call of release 0.8.0 makes rules.Model write model/1.0; plain
    # calls keep rules-model/0.3.
    for call in (
        lambda m: m.add_species("B", attributes=["a"]),
        lambda m: m.add_observable("twice", "A", weight=2),
        lambda m: m.add_reaction(0.1, [crn.pattern("A")], []),
        lambda m: m.add_reaction(0.1, ["A"], [crn.copy("x")]),
        lambda m: m.add_reaction(0.1, ["A"], [crn.new("A")]),
        lambda m: m.add_reaction(0.1, [{"label": "x", "type": "A"}], []),
    ):
        plain = rules.Model("plain").add_species("A", 3)
        plain.add_reaction(0.1, ["A"], [], name="decay")
        require(
            plain.to_dict()["schema"] == "crn-craft.rules-model/0.3",
            "plain rules.Model calls write rules-model/0.3",
        )
        call(plain)
        require(
            plain.to_dict()["schema"] == "crn-craft.model/1.0"
            and crn.check(plain)["valid"],
            "a call of 0.8.0 makes rules.Model write model/1.0",
        )
    # A legacy document of from_dict keeps its version until such a call.
    legacy = dict(clone_document(1e-3), schema="crn-craft.rules-model/0.3")
    loaded = rules.Model.from_dict(legacy).add_parameter("unused", 1.0)
    require(
        loaded.to_dict()["schema"] == "crn-craft.rules-model/0.3"
        and crn.check(loaded)["valid"],
        "a loaded rules-model/0.3 document keeps its version",
    )
    loaded.add_observable("largest", "mtDNA", reduce="max_count")
    require(
        loaded.to_dict()["schema"] == "crn-craft.model/1.0"
        and crn.check(loaded)["valid"],
        "add_observable on a loaded legacy document writes model/1.0",
    )
    # add_initial: every attribute, counts add up per attribute vector.
    m = rules.Model("init").add_parameter("n0", 4.0)
    m.add_species("C", attributes=["strain", "w"])
    m.add_initial("C", 10, strain=0, w=1)
    m.add_initial("C", 5, w=1, strain=0)
    m.add_initial("C", "n0", strain=1, w=0.5)
    m.add_initial("C", 0, strain=2, w=1)
    m.add_species("A", 3).add_initial("A", 4)
    m.add_reaction(0.1, ["A"], [], name="decay")
    require(
        m.to_dict()["types"]
        == [{"name": "C", "attributes": ["strain", "w"]}, {"name": "A"}]
        and m.to_dict()["initial"]
        == [
            {"type": "C", "attributes": {"strain": 0, "w": 1}, "count": 15},
            {"type": "C", "attributes": {"strain": 1, "w": 0.5}, "count": "n0"},
            {"type": "A", "count": 7},
        ]
        and m.to_dict()["observables"] == {"C": {"type": "C"}, "A": {"type": "A"}}
        and crn.check(m)["valid"],
        "add_initial",
    )
    # Patterns, copy and new, and their dictionary forms.
    a = crn.Model("forms").add_species("C", attributes=["s"]).add_initial("C", 9, s=1)
    a.add_observable("heavy", "C", weight="x.s", label="x")
    a.add_observable("biggest", "C", weight=2, reduce="max")
    a.add_reaction(
        crn.mass_action("0.1"),
        [crn.pattern("C", where="x.s < 5", factor="x.s")],
        [crn.copy("x"), crn.new("C", s="x.s + 1")],
        name="helpers",
    )
    a.add_reaction(
        crn.mass_action("0.1"),
        [{"label": "x", "type": "C", "where": "x.s < 5", "factor": "x.s"}],
        [{"type": "C", "copy": "x"}, {"type": "C", "attributes": {"s": "x.s + 1"}}],
        name="dictionaries",
    )
    a.add_reaction(0.01, ["C", crn.pattern("C")], [crn.copy("y", s="x.s + y.s")])
    rules_of = {r["name"]: r for r in a.to_dict()["rules"]}
    require(
        rules_of["helpers"]["reactants"] == rules_of["dictionaries"]["reactants"]
        and rules_of["helpers"]["products"] == rules_of["dictionaries"]["products"]
        and rules_of["helpers"]["reactants"]
        == [{"label": "x", "type": "C", "where": "x.s < 5", "factor": "x.s"}]
        and rules_of["r2"]["reactants"]
        == [{"label": "x", "type": "C"}, {"label": "y", "type": "C"}]
        and rules_of["r2"]["products"]
        == [{"type": "C", "copy": "y", "attributes": {"s": "x.s + y.s"}}]
        and a.to_dict()["observables"]["heavy"]
        == {"type": "C", "label": "x", "weight": "x.s"}
        and a.to_dict()["observables"]["biggest"]
        == {"type": "C", "weight": "2.0", "reduce": "max"}
        and crn.check(a)["valid"],
        "patterns, copy and new, as helpers and as dictionaries",
    )
    # Early errors, and nothing changes when a call fails.
    before = a.to_dict()
    for error, call in (
        (ValueError, lambda: a.add_reaction(0.1, ["C"], ["C"])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [crn.copy("z")])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [crn.copy("x", size=1)])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [crn.new("C")])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [crn.new("C", s=1, t=2)])),
        (KeyError, lambda: a.add_reaction(0.1, ["C"], [crn.new("D", s=1)])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [{"type": "B", "copy": "x"}])),
        (ValueError, lambda: a.add_reaction(0.1, ["C"], [{"type": "C", "cpy": "x"}])),
        (ValueError, lambda: a.add_reaction(0.1, [{"type": "C", "lable": "x"}], [])),
        (ValueError, lambda: a.add_reaction(0.1, ["C", "C", "C"], [])),
        (
            ValueError,
            lambda: a.add_reaction(0.1, [crn.pattern("C", label="y"), "C"], []),
        ),
        (
            ValueError,
            lambda: a.add_reaction("k", [crn.pattern("C", where="x.s > 0")], []),
        ),
        (TypeError, lambda: a.add_reaction(0.1, [3], [])),
        (TypeError, lambda: a.add_reaction(0.1, ["C"], [3])),
        (ValueError, lambda: a.add_initial("C", 3)),
        (ValueError, lambda: a.add_initial("C", 3, s=1, t=2)),
        (KeyError, lambda: a.add_initial("D", 3)),
        (ValueError, lambda: a.set_initial({"C": 3})),
        (ValueError, lambda: a.add_species("D", 5, ["x1"])),
        (TypeError, lambda: a.add_species({"D": 0, "E": 0}, attributes=["x1"])),
        (TypeError, lambda: a.add_species({"D": 3}, attributes=["x1"])),
        (TypeError, lambda: a.add_species("D", 0, "x1")),
        (
            ValueError,
            lambda: a.add_reaction(0.1, [crn.pattern("C", where="x.s > 0"), "C"], []),
        ),
        (
            ValueError,
            lambda: a.add_reaction(
                0.1,
                [crn.pattern("C", factor="x.s"), crn.pattern("C", factor="y.w")],
                [],
            ),
        ),
        (ValueError, lambda: a.add_species("D", 0, ["a", "a"])),
        (ValueError, lambda: a.add_observable("C", "C")),
        (ValueError, lambda: a.add_observable("heavy", "C")),
        (ValueError, lambda: a.add_observable("q", "C", None, None, "mean")),
        (KeyError, lambda: a.add_observable("q", "D")),
    ):
        raises(error, call)
    require(a.to_dict() == before, "failed calls leave the model unchanged")
    message = str(raises(ValueError, a.add_species, "D", 5, ["x1"]))
    require(
        "add_species takes no initial count" in message
        and "initial=" not in message
        and "add_initial('D', count, x1=...)" in message,
        "the message for an initial count of a species with attributes",
    )
    message = str(raises(TypeError, a.add_species, {"D": 0}, attributes=["x1"]))
    require(
        "in its own call" in message,
        "a mapping with attributes asks for one call per species",
    )
    message = str(
        raises(
            ValueError,
            a.add_reaction,
            0.1,
            [crn.pattern("C", where="x.s > 0"), "C"],
            [],
        )
    )
    require(
        "two different patterns over the same species 'C'" in message
        and "identical" in message,
        "two different patterns over one species are refused",
    )
    # Identical patterns, whatever their labels and spacing, form a pair.
    pairs = crn.Model("pairs").add_species("C", attributes=["s"])
    pairs.add_initial("C", 9, s=1)
    pairs.add_reaction(
        0.1,
        [crn.pattern("C", where="x.s > 0"), crn.pattern("C", where="(y.s>0)")],
        [crn.copy("x")],
        name="pair",
    )
    require(crn.check(pairs)["valid"], "identical patterns with other labels pair")
    message = str(raises(ValueError, a.add_reaction, 0.1, ["C"], ["C"]))
    require(
        "crn_craft.copy" in message and "crn_craft.new" in message,
        "a bare product of a species with attributes names copy and new",
    )
    implicit = crn.Model("implicit").add_reaction(0.1, ["A"], ["B"])
    raises(ValueError, implicit.add_species, "A", attributes=["a"])
    raises(ValueError, implicit.add_parameter, "A", 1.0)
    raises(ValueError, implicit.add_observable, "B", "A")
    # Observable names share the namespace.
    named = (
        crn.Model("names").add_species("A", 1).add_observable("twice", "A", weight=2)
    )
    raises(ValueError, named.add_parameter, "twice", 1.0)
    raises(ValueError, named.add_species, "twice", 1)
    raises(ValueError, named.add_ode, "twice", 0.0, "1")
    # Attribute-free models keep their 0.7.1 documents.
    plain = sir(rules.Model)
    plain.add_initial("S", 0)
    require(
        rules.model_bytes(plain.to_dict()) == SIR_0_6_0,
        "attribute-free calls write the documents of 0.7.1",
    )
    # to_dict / from_dict round trips stay editable.
    for model in (built, old, m, a):
        document = model.to_dict()
        back = type(model).from_dict(document)
        require(
            not back._fixed and back.to_dict() == document,
            f"from_dict of {model.name}: editable and unchanged",
        )
        back.add_parameter("extra_parameter", 1.0)
        edited = back.to_dict()
        edited["parameters"].pop("extra_parameter")
        require(edited == document, f"from_dict of {model.name}: the edit only")
        require(
            type(model).from_dict(model.to_json()).to_dict() == document,
            f"from_dict(to_json()) of {model.name}",
        )
    # The hand-written document is editable and keeps its observables.
    loaded = crn.Model.from_dict(hand)
    require(
        loaded.to_dict() == hand and loaded.species == ("mtDNA",),
        "from_dict keeps the hand-written document",
    )
    loaded.add_observable("largest", "mtDNA", where="x.clone > 0", reduce="max_count")
    require(
        list(loaded.to_dict()["observables"])
        == ["N", "mutants", "types", "mutTypes", "largest"]
        and loaded.to_dict()["rules"] == hand["rules"]
        and crn.check(loaded)["valid"],
        "an edited hand-written document gets no extra count observable",
    )
    # crn.Model also edits rules-model/0.1 documents; rules.Model keeps them.
    three = json.loads((ROOT / "models/rules/mtdna-three-class.rules.json").read_text())
    unified = crn.Model.from_dict(three)
    unified.add_observable("largest", "mtDNA", reduce="max_count")
    require(
        crn.check(unified)["valid"]
        and unified.to_dict()["schema"] == "crn-craft.model/1.0"
        and unified.to_dict()["rules"] == three["rules"]
        and rules.Model.from_dict(three)._fixed,
        "crn.Model edits a rules-model/0.1 document",
    )
    # Exported helpers and their representation.
    require(
        crn.copy is rules.copy
        and crn.new is rules.new
        and crn.pattern is rules.pattern
        and crn.mass_action is rules.mass_action
        and repr(crn.pattern("C", where="x.s > 0")) == "pattern('C', where='x.s > 0')"
        and repr(crn.copy("x", s=1)) == "copy('x', s=1)"
        and repr(crn.new("C", s="new_id()")) == "new('C', s='new_id()')"
        and repr(crn.mass_action("k")) == "mass_action('k')",
        "builder helpers",
    )
    raises(TypeError, crn.copy, "x", s=True)
    raises(ValueError, crn.mass_action, "")
    raises(ValueError, crn.pattern, "C", where="")
    rows.append({"check": "builder_attributes", "passed": True})


def check_clone_model(rows):
    """The builder's clone model runs the hand-written model's paths, and
    max_count equals the largest clone of the final instances."""
    settings = {
        "end": 8640.0,
        "report_times": [0.0, 4320.0, 8640.0],
        "replicates": 8,
        "seed": 3,
        "threads": 2,
        "keep_paths": "all",
    }
    hand = crn.ensemble(clone_document(1e-3), **settings).simulate()
    for degradation in (DEG_RATE, crn.mass_action("delta")):
        model = clone_model(crn.Model, 1e-3, degradation)
        model.add_observable(
            "largest", "mtDNA", where="x.clone > 0", reduce="max_count"
        )
        result = crn.ensemble(model, keep_instances=8, **settings).simulate()
        require(
            result.metadata["sampler"] == hand.metadata["sampler"]
            and result.records["events"].equals(hand.records["events"])
            and all(
                np.array_equal(
                    result.paths.sel(variable=name).values,
                    hand.paths.sel(variable=name).values,
                )
                for name in ("N", "mutants", "types", "mutTypes")
            ),
            f"builder and hand-written clone models: identical paths ({degradation})",
        )
        final = result.paths.sel(variable="largest").isel(time=-1).values
        largest = []
        for cell in result.metadata["final_instances"]:
            counts = [i["count"] for i in cell["instances"] if i["attributes"][0] > 0]
            largest.append(max(counts, default=0))
        require(
            final.tolist() == largest and max(largest) > 1,
            "max_count is the largest mutant clone of the final instances",
        )

    # A whole propensity on a species with attributes: the tree method fires
    # at its value and takes a molecule uniformly (the process of the
    # per-molecule rate); the static methods refuse it over two clones.
    def two_clones(rate):
        model = crn.Model("two").add_parameter({"delta": 0.01, "K": 100})
        model.add_species("M", attributes=["clone"])
        model.add_initial("M", 50, clone=0).add_initial("M", 50, clone=1)
        model.add_observable("c1", "M", where="x.clone == 1")
        model.add_reaction(rate, ["M"], [], name="degradation")
        model.add_reaction(
            crn.mass_action("max(0, K - M) / M"),
            ["M"],
            [crn.copy("x"), crn.copy("x")],
            name="replication",
        )
        return model

    short = {"end": 200.0, "report_times": [0.0, 200.0], "replicates": 200, "seed": 1}
    per_molecule = crn.ensemble(two_clones(crn.mass_action("delta")), **short)
    whole = crn.ensemble(two_clones("delta * M"), **short)
    a, b = per_molecule.simulate(), whole.simulate()
    difference = abs(
        a.statistics["mean"]["c1"].iloc[-1] - b.statistics["mean"]["c1"].iloc[-1]
    )
    error = math.sqrt(
        (
            a.statistics["sample_variance"]["c1"].iloc[-1]
            + b.statistics["sample_variance"]["c1"].iloc[-1]
        )
        / 200
    )
    require(difference <= 4 * error, "whole propensity and mass_action agree")
    refused = raises(
        crn.CrnError,
        crn.ensemble,
        two_clones("delta * M"),
        method="static_direct",
        **short,
    )
    require(
        refused.diagnostics[0]["path"] == "/method"
        and crn.ensemble(
            two_clones(crn.mass_action("delta")), method="static_direct", **short
        )
        .simulate()
        .metadata["termination"]
        == "completed",
        "static methods take mass_action, not a whole propensity, over two clones",
    )
    # A one-rule Moran step "x replaces y" through crn.rules: written as
    # model/1.0 and run/1.0, so the pair's labels are random (rules-direct/3)
    # and clone 1 keeps its mean of 20 (neutral drift).
    moran = rules.Model("moran").add_parameter("k", 0.01)
    moran.add_species("M", attributes=["clone"])
    moran.add_initial("M", 20, clone=0).add_initial("M", 20, clone=1)
    moran.add_reaction(
        crn.mass_action("k"),
        ["M", "M"],
        [crn.copy("x"), crn.copy("x")],
        name="x_replaces_y",
    )
    run = rules.ensemble(
        moran,
        end=200.0,
        report_times=[0.0, 200.0],
        replicates=400,
        seed=1,
        threads=2,
        keep_instances=400,
    )
    result = run.simulate()
    clone_1 = np.array(
        [
            sum(i["count"] for i in cell["instances"] if i["attributes"][0] == 1)
            for cell in result.metadata["final_instances"]
        ],
        dtype=float,
    )
    require(
        json.loads(run._documents[0])["schema"] == "crn-craft.model/1.0"
        and json.loads(run._documents[1])["schema"] == "crn-craft.run/1.0"
        and result.metadata["sampler"] == "chacha8-u53/rules-direct/3"
        and abs(clone_1.mean() - 20) <= 4 * clone_1.std(ddof=1) / np.sqrt(400),
        "crn.rules runs an attributed pair rule with random labels "
        f"({result.metadata['sampler']}, clone 1 mean {clone_1.mean()})",
    )
    rows.append({"check": "clone_model_builder", "passed": True})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = []
    check_names(rows)
    model = crn.load_example("immigration-death")

    def ensemble(threads, **extra):
        return crn.ensemble(
            model,
            end=2.0,
            report_times=[0.0, 1.0, 2.0],
            replicates=64,
            seed=5001,
            threads=threads,
            keep_paths=64,
            method="static_direct",
            **extra,
        )

    expected = ensemble(1).simulate()
    for threads in [1, 2, 8]:
        prepared = ensemble(threads)
        actual = prepared.simulate()
        normalized = dict(actual.metadata)
        normalized.pop("run_sha256")
        exp = dict(expected.metadata)
        exp.pop("run_sha256")
        require(
            normalized == exp
            and actual.records.equals(expected.records)
            and actual.rule_events.equals(expected.rule_events),
            "binding thread identity",
        )
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            (folder / "model.rules.json").write_bytes(prepared._documents[0])
            (folder / "run.json").write_bytes(prepared._documents[1])
            proc = subprocess.run(
                [
                    str(args.binary.resolve()),
                    "run",
                    str(folder / "model.rules.json"),
                    str(folder / "run.json"),
                    "--archive",
                    str(folder / "out"),
                ],
                capture_output=True,
                check=False,
            )
            require(proc.returncode == 0, proc.stdout.decode())
            require(
                cli_identity(
                    json.loads((folder / "out/result.json").read_bytes()), actual
                ),
                "CLI binding identity",
            )
        rows.append(
            {
                "check": f"CLI_binding_threads_{threads}",
                "passed": True,
                "numerical_sha256": hashlib.sha256(
                    json.dumps(normalized, sort_keys=True).encode()
                ).hexdigest(),
            }
        )
    # Parameter search reuses the prepared run; every candidate replays from
    # the documents it reports, and the loop source needs no toolchain.
    run = ensemble(1)
    values = {name: value * 1.2 for name, value in model.parameters.items()}
    bound = run.bind_parameters(values)
    first = bound.simulate()
    # A candidate reports its model with the bound values and the request
    # bound to the engine's re-serialized bytes (M6b); prepared again with
    # the hash of these bytes, they replay the candidate's numbers.
    document, request = first.model, first.request
    request["model"]["sha256"] = hashlib.sha256(rules.model_bytes(document)).hexdigest()
    replay = crn.prepare(document, request).simulate()
    without = ("model_sha256", "run_sha256")
    require(
        {k: v for k, v in replay.metadata.items() if k not in without}
        == {k: v for k, v in first.metadata.items() if k not in without}
        and document["parameters"] == {k: float(v) for k, v in values.items()}
        and replay.records.equals(first.records)
        and replay.rule_events.equals(first.rule_events)
        and replay.paths.equals(first.paths),
        "parameter binding replay",
    )
    require(
        run.bind_parameters(values).simulate().metadata == first.metadata,
        "repeated parameter binding",
    )
    error = raises(crn.CrnError, run.bind_parameters, {"birth": True})
    require(
        error.diagnostics[0]["path"] == "/parameters/birth",
        "boolean parameter was refused",
    )
    require(
        crn.loop_source(model).startswith("// Generated"),
        "generated loop source interface",
    )
    rows.append({"check": "reusable_parameters_and_codegen", "passed": True})
    for invalid in [b'{"name":"a","name":"b"}', b'{"schema":"crn-craft.model/99"}']:
        require(not crn.check(invalid)["valid"], "strict binding parser")
    rows.append({"check": "strict_inputs", "passed": True})
    check_documents(rows)
    check_models(rows)
    check_methods(rows)
    check_cli(rows, args.binary.resolve())
    check_limits(rows, args.binary.resolve())
    check_incomplete(rows)
    check_attributed_builder(rows)
    check_clone_model(rows)
    t = np.linspace(0, 100, 501)
    for split, seeds in [
        ("development", [101, 102]),
        ("held_out", [201, 202, 203, 204]),
    ]:
        correct = 0
        for seed in seeds:
            rng = np.random.default_rng(seed)
            positive = np.sin(2 * np.pi * t / 10) + 0.1 * rng.normal(size=len(t))
            negative = rng.normal(size=len(t))
            correct += oscillation(t, positive)["classification"] == "oscillatory"
            correct += oscillation(t, negative)["classification"] == "nonoscillatory"
        require(correct == 2 * len(seeds), "synthetic held-out detector checks")
        rows.append(
            {
                "check": "synthetic_" + split,
                "passed": True,
                "correct": correct,
                "total": 2 * len(seeds),
                "scope": "labelled synthetic controls; no biological sensitivity/specificity claim",
            }
        )
    first = sampled_first_passage([0, 1, 2], [[0, 2, 3], [0, 0, 0], [3, 3, 3]], 2)
    require(
        first["not_detected"] == 1
        and first["detected_fraction"] == 2 / 3
        and first["restricted_mean_upper"] == 1.0
        and "restricted_mean_lower" not in first,
        "saved-grid passage accounting",
    )
    # A crossing between saved times: the coarse grid misses it, so it must
    # neither report a passage nor a lower bound; the finer grid sees it.
    fine_t = [0.0, 0.5, 1.0, 1.5, 2.0]
    fine_x = [[0, 3, 1, 1, 0]]
    coarse = sampled_first_passage(fine_t[::2], [row[::2] for row in fine_x], 2)
    fine = sampled_first_passage(fine_t, fine_x, 2)
    require(
        coarse["detected"] == 0
        and fine["detected"] == 1
        and fine["restricted_mean_upper"] <= coarse["restricted_mean_upper"],
        "nested saved grids",
    )
    rows.append({"check": "saved_grid_first_passage", "passed": True})
    # A constant non-integer series is nonoscillatory, not an error.
    flat = oscillation(np.arange(100.0), np.full(100, 0.1))
    require(flat["classification"] == "nonoscillatory", "constant series")
    try:
        oscillation(np.arange(20.0), np.sin(np.arange(20.0)), burn_fraction=0.79)
    except ValueError:
        pass
    else:
        raise RuntimeError("too few samples after the burn-in were accepted")
    rows.append({"check": "oscillation_edge_cases", "passed": True})
    # The detector's three documented gates, each checked on its own:
    # removing the mean (an offset changes nothing), the autocorrelation gate
    # (a spectral peak whose lag-period autocorrelation is below 0.2 is not
    # called oscillatory) and the burn-in (samples before it are ignored).
    rng = np.random.default_rng(5)
    base = np.sin(2 * np.pi * t / 10) + 0.1 * rng.normal(size=len(t))
    plain, shifted = oscillation(t, base), oscillation(t, base + 1000.0)
    require(
        plain["classification"] == shifted["classification"] == "oscillatory"
        and abs(plain["autocorrelation_peak"] - shifted["autocorrelation_peak"]) < 1e-9,
        "oscillation is offset invariant",
    )
    noise = np.random.default_rng(12).normal(scale=np.sqrt(0.8), size=500)
    n = np.arange(500.0)
    aliased = oscillation(n, np.sin(2 * np.pi * 0.4 * n) + noise, burn_fraction=0)
    require(
        aliased["peak_fraction"] >= 0.3
        and aliased["cycles"] >= 3
        and aliased["autocorrelation_peak"] < 0.2
        and aliased["classification"] == "nonoscillatory",
        "autocorrelation gate",
    )
    transient = base.copy()
    transient[: int(len(t) * 0.2)] += 1000.0
    require(
        oscillation(t, transient) == plain
        and oscillation(t, transient, burn_fraction=0)["classification"]
        != "oscillatory",
        "burn-in excludes the transient",
    )
    rows.append({"check": "oscillation_gates", "passed": True})
    # Rejections carry their diagnostics; a cancel file stops a run; Ctrl-C
    # interrupts a long native call promptly instead of being ignored.
    request = json.loads(ensemble(1)._documents[1])
    try:
        crn.prepare(model, dict(request, schema="crn-craft.run/9"))
    except crn.CrnError as error:
        require(
            isinstance(error, ValueError)
            and error.diagnostics[0]["code"] == "CRN-E003"
            and str(error).startswith("CRN-E003"),
            "CrnError diagnostics",
        )
    else:
        raise RuntimeError("an invalid run was accepted")
    with tempfile.TemporaryDirectory() as tmp:
        stop = Path(tmp) / "stop"
        stop.write_text("")
        long = crn.ensemble(
            model,
            end=1e6,
            report_times=[0.0, 1e6],
            replicates=8,
            max_events=10**9,
            max_sampling_steps=10**10,
            method="static_direct",
        )
        cancelled = long.simulate(cancel_file=stop)
        require(cancelled.metadata["termination"] == "cancelled", "binding cancel file")
    rows.append({"check": "binding_errors_and_cancel_file", "passed": True})
    if os.name == "posix":
        script = (
            "import crn_craft as c, time\n"
            "m = c.load_example('immigration-death')\n"
            "r = c.ensemble(m, end=1e9, report_times=[0.0, 1e9], replicates=4,\n"
            "               max_events=10**12, max_sampling_steps=10**13,\n"
            "               method='static_direct')\n"
            "print('started', flush=True)\n"
            "try:\n"
            "    r.simulate()\n"
            "except KeyboardInterrupt:\n"
            "    print('interrupted', flush=True)\n"
        )
        child = subprocess.Popen(
            [sys.executable, "-I", "-c", script],
            stdout=subprocess.PIPE,
            text=True,
        )
        require(child.stdout.readline().strip() == "started", "child started")
        time.sleep(0.5)
        sent = time.monotonic()
        child.send_signal(signal.SIGINT)
        out, _ = child.communicate(timeout=30)
        latency = time.monotonic() - sent
        require(
            out.strip() == "interrupted" and child.returncode == 0 and latency < 2.0,
            f"Ctrl-C interrupts a native run ({latency:.2f} s, {out!r})",
        )
        rows.append({"check": "keyboard_interrupt", "passed": True})
    circuit = crn.load_example("repressilator")
    # The first species of the network model, as the 0.6.0 study used.
    species = json.loads((ROOT / "models/legacy/repressilator.model.json").read_text())[
        "species"
    ][0]["name"]
    study = []
    for split, seeds in [("development", [3101, 3102]), ("held_out", [3201, 3202])]:
        for seed in seeds:
            result = crn.ensemble(
                circuit,
                end=100.0,
                report_times=np.linspace(0, 100, 501).tolist(),
                replicates=4,
                seed=seed,
                keep_paths=4,
                method="static_direct",
            ).simulate()
            require(result.metadata["termination"] == "completed", "circuit incomplete")
            paths = result.paths.sel(variable=species)
            diagnostics = [
                oscillation(paths["time"].values, paths.sel(replicate=r).values)
                for r in paths["replicate"].values
            ]
            study.append(
                {
                    "split": split,
                    "seed": seed,
                    "events": int(result.records["events"].sum()),
                    "classification": [d["classification"] for d in diagnostics],
                    "scope": "synthetic fixed topology illustration; circuit-study source not supplied",
                }
            )
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.python-validation/0.2",
                "passed": True,
                "rows": rows,
                "circuit_study": study,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "checks": len(rows), "circuit_paths": 16}))


if __name__ == "__main__":
    main()
