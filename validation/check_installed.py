"""Acceptance from an installed wheel, independent of repository imports.

Run it with the environment's Python. The console scripts are taken from that
environment's scripts directory, never from PATH, so a native `crn-craft`
binary elsewhere on PATH cannot stand in for the installed entry point.
"""

import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import shutil
import subprocess
import sysconfig
import tempfile
import numpy as np
import pandas as pd
import crn_craft as c
from crn_craft import rules
from crn_craft.importers import import_bng_net
from crn_craft.observables import oscillation, sampled_first_passage

SBML = (
    b'<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" '
    b'version="2"><model id="Decay" substanceUnits="item" timeUnits="second">'
    b'<listOfCompartments><compartment id="cell" size="1" constant="true"/>'
    b'</listOfCompartments><listOfSpecies><species id="A" compartment="cell" '
    b'initialAmount="10" hasOnlySubstanceUnits="true" boundaryCondition="false" '
    b'constant="false"/></listOfSpecies><listOfParameters><parameter id="k" '
    b'value="0.5" constant="true"/></listOfParameters><listOfReactions>'
    b'<reaction id="decay" reversible="false"><listOfReactants><speciesReference '
    b'species="A" stoichiometry="1" constant="true"/></listOfReactants>'
    b'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/>'
    b"<ci>k</ci><ci>A</ci></apply></math></kineticLaw></reaction>"
    b"</listOfReactions></model></sbml>"
)


def script(name):
    """The console script installed with this interpreter's package."""
    found = shutil.which(name, path=sysconfig.get_path("scripts"))
    if found is None:
        raise RuntimeError(f"{name} is not installed next to this interpreter")
    return found


def require(v, m):
    if not v:
        raise RuntimeError(m)


def python_spelling(cargo):
    """Cargo's `0.3.0-dev.1` is Python's `0.3.0.dev1`; releases are equal."""
    release, _, pre = cargo.partition("-")
    if not pre:
        return release
    kind, _, number = pre.partition(".")
    return (
        release + {"dev": ".dev", "alpha": "a", "beta": "b", "rc": "rc"}[kind] + number
    )


def three_class_model():
    """models/rules/mtdna-three-class.rules.json, without its description.

    Written out here so that the check needs no file from the checkout.
    """

    def cls(clone, s, wt):
        return {"clone": clone, "s": s, "wt": wt}

    def pattern(clone):
        return [{"label": "x", "type": "mtDNA", "where": f"x.clone == {clone}"}]

    copy = {"type": "mtDNA", "copy": "x"}
    vacancy = "max(0, K - N) / D * "

    def rule(name, clone, rate, products):
        return {
            "name": name,
            "reactants": pattern(clone),
            "rate": rate,
            "products": products,
        }

    return {
        "schema": "crn-craft.rules-model/0.1",
        "name": "mtdna_three_class",
        "parameters": {
            "K": 1000,
            "delta": 0.0028881132523331052,
            "sw": 1.0,
            "sm": 1.5,
            "mu": 3e-05,
            "fracAdv": 0.2,
            "w0": 1000,
        },
        "types": [{"name": "mtDNA", "attributes": ["clone", "s", "wt"]}],
        "initial": [{"type": "mtDNA", "attributes": cls(0, "sw", 1), "count": "w0"}],
        "observables": {
            "N": {"type": "mtDNA"},
            "wild_type": {"type": "mtDNA", "where": "x.wt == 1"},
            "mutants": {"type": "mtDNA", "where": "x.wt == 0"},
            "advantaged_mutants": {"type": "mtDNA", "where": "x.clone == 2"},
            "D": {"type": "mtDNA", "weight": "x.wt + sm * (1 - x.wt)"},
            "neutral_mutants": {"type": "mtDNA", "where": "x.clone == 1"},
        },
        "rules": [
            rule("degradation_wt", 0, "delta", []),
            rule("replication_wt", 0, vacancy + "(sw * (1 - mu))", [copy, copy]),
            rule(
                "mutation_mt",
                0,
                vacancy + "sw * mu * (1 - fracAdv)",
                [copy, {"type": "mtDNA", "attributes": cls(1, "sw", 0)}],
            ),
            rule(
                "mutation_mtAdv",
                0,
                vacancy + "sw * mu * fracAdv",
                [copy, {"type": "mtDNA", "attributes": cls(2, "sm", 0)}],
            ),
            rule("degradation_mt", 1, "delta", []),
            rule("replication_mt", 1, vacancy + "(sw)", [copy, copy]),
            rule("degradation_mtAdv", 2, "delta", []),
            rule("replication_mtAdv", 2, vacancy + "(sm)", [copy, copy]),
        ],
    }


def check_event_log_default_budget():
    """The event-log request from the 0.2.0 feedback, at the default budget.

    Three-class model, 400 cells for three years, up to 50 logged events of
    the two mutation rules per cell, default max_output_bytes (64 MiB).
    Release 0.2.0 rejected it before running with CRN-R020 (observation
    reservation exceeds budget).
    """
    names = ["mutation_mt", "mutation_mtAdv"]
    try:
        run = rules.ensemble(
            three_class_model(),
            end=25_920.0,
            report_times=[0.0, 25_920.0],
            replicates=400,
            seed=4,
            threads=2,
            method="static_direct",
            event_log={"max_entries": 50, "rules": names},
        )
    except c.CrnError as error:
        raise RuntimeError(
            f"the reported event-log request is rejected: {error}"
        ) from None
    result = run.simulate()
    log = result.event_columns()
    omitted = sum(o["omitted_log_entries"] for o in result.observations)
    mutations = result.rule_events[names]
    counts = mutations.sum(axis=1)
    logged = log.groupby(level="replicate").size().reindex(counts.index, fill_value=0)
    require(
        result.metadata["termination"] == "completed"
        and bool((result.records["termination"] == "completed").all())
        and len(result.observations) == 400
        and set(log["rule_name"].unique()) <= set(names)
        and omitted == 0
        and len(log) == int(mutations.to_numpy().sum())
        and len(log) > 0
        and np.array_equal(logged.to_numpy(), counts.to_numpy()),
        "the reported event log records every mutation of every cell",
    )


def check_builder():
    """rules.Model: a string is a whole propensity, a number mass action."""
    model = rules.Model("birth_death")
    model.add_parameter("k", 10.0).add_species("A")
    model.add_reaction("k", [], ["A"], name="birth")
    model.add_reaction(1.0, ["A"], [], name="death")
    hand = {
        "schema": "crn-craft.rules-model/0.3",
        "name": "birth_death",
        "parameters": {"k": 10.0},
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
                "rate": "1.0",
                "products": [],
            },
        ],
    }
    require(
        model.to_dict() == hand and json.loads(model.to_json()) == hand,
        "the builder writes the hand-written rules-model/0.3 document",
    )
    times = [0.0, 1.0, 5.0]
    # Birth-death from zero: Poisson with mean k (1 - exp(-t)).
    exact = 10.0 * (1.0 - np.exp(-np.array(times)))
    results = {}
    for method in ("static_direct", None):
        results[method] = rules.ensemble(
            model,
            end=5.0,
            report_times=times,
            replicates=1000,
            seed=11,
            threads=2,
            method=method,
        ).simulate()
        mean = results[method].statistics["mean"]["A"].to_numpy()
        error = results[method].statistics["mean_standard_error"]["A"].to_numpy()
        require(
            mean[0] == 0.0
            and bool(np.all(np.abs(mean[1:] - exact[1:]) <= 5 * error[1:])),
            f"builder birth-death matches the Poisson law ({method or 'tree'})",
        )
    same = rules.ensemble(
        hand,
        end=5.0,
        report_times=times,
        replicates=1000,
        seed=11,
        threads=2,
        method="static_direct",
    ).simulate()
    require(
        same.metadata["model_sha256"]
        == results["static_direct"].metadata["model_sha256"]
        and same.records.equals(results["static_direct"].records)
        and results["static_direct"].rule_names == ("birth", "death"),
        "a builder model runs exactly as its JSON",
    )


def check_pandas_results(result):
    """Statistics by name, NumPy access and categorical terminations."""
    meta = result.metadata
    observables = list(meta["observables"])
    for key in ("mean", "sample_variance", "mean_standard_error"):
        frame = result.statistics[key]
        require(
            isinstance(frame, pd.DataFrame)
            and frame.index.name == "time"
            and list(frame.index) == meta["report_times"]
            and list(frame.columns) == observables
            and np.array_equal(
                frame.to_numpy(),
                np.array(meta["statistics"][key], dtype=np.float64).reshape(
                    frame.shape
                ),
            ),
            f"statistics[{key!r}] is the raw matrix by time and observable",
        )
    termination = result.records["termination"]
    codes = termination.cat.codes.to_numpy()
    require(
        list(termination.cat.categories) == list(rules.TERMINATIONS)
        and bool(np.all(codes == rules.TERMINATIONS.index("completed")))
        and isinstance(result.records["stop_time"].to_numpy(), np.ndarray),
        "termination codes and NumPy columns stay available",
    )


def check_hybrid():
    """A tiny hybrid run: x' = 2 from 0 drives immigration at rate x."""
    ramp = rules.Model("ramp")
    ramp.add_species("A")
    ramp.add_continuous("x", 0.0, "2")
    ramp.add_reaction("x", [], ["A"], name="immigration")
    try:
        rules.ensemble(ramp, end=1.0, report_times=[0.0, 1.0], method="static_direct")
    except c.CrnError as error:
        require(error.diagnostics[0]["code"] == "CRN-R020", "hybrid-only refusal")
    else:
        raise RuntimeError("static_direct accepted a continuous variable")
    times = [0.0, 1.5, 3.0]
    result = rules.ensemble(
        ramp,
        end=3.0,
        report_times=times,
        replicates=400,
        seed=13,
        threads=2,
        method="hybrid",
    ).simulate()
    meta = result.metadata
    mean = result.statistics["mean"]
    error = result.statistics["mean_standard_error"]["A"].to_numpy()
    # x = 2 t, and the count at t is Poisson with mean t^2.
    exact = np.array(times) ** 2
    require(
        meta["schema"] == "crn-craft.rules-result/0.3"
        and meta["method"] == {"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}
        and meta["continuous"] == ["x"]
        and meta["sampler"] == "chacha8-u53/rules-approximate-hybrid-dopri54/1"
        and list(mean.columns) == ["A", "x"]
        and bool((result.records["termination"] == "completed").all())
        and np.allclose(mean["x"], [0.0, 3.0, 6.0], rtol=1e-6, atol=1e-9)
        and mean["A"].iloc[0] == 0.0
        and bool(np.all(np.abs(mean["A"].to_numpy()[1:] - exact[1:]) <= 5 * error[1:])),
        "installed hybrid integrates x and matches the Poisson law",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    checks = []
    cli, importer = script("crn-craft"), script("crn-craft-import")
    version = importlib.metadata.version("crn-craft")
    require(c.__version__ == version, "package __version__ matches metadata")
    require(
        python_spelling(c._native.__version__) == version,
        "native extension and distribution versions agree",
    )
    reported = subprocess.run([cli, "--version"], capture_output=True, check=False)
    require(
        reported.returncode == 0
        and json.loads(reported.stdout) == {"version": version, "interface": "python"},
        "installed entry point is the Python CLI",
    )
    checks.append("version_and_entry_points")
    m = c.load_example("immigration-death")
    r = c.run_spec(m, end=2.0, replicates=32, seed=77, keep=32)
    a = c.simulate(m, r)
    b = c.simulate(m, r, cached=True)
    require(a == b, "installed cached identity")
    checks.append("native_API_exact_identity")
    require(
        c.describe()["schema"] == "crn-craft.describe/0.10", "installed capabilities"
    )
    checks.append("discovery")
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        (folder / "model.json").write_bytes(m)
        (folder / "run.json").write_bytes(c.model_bytes(r))
        proc = subprocess.run(
            [
                cli,
                "simulate-cached",
                str(folder / "model.json"),
                str(folder / "run.json"),
                str(folder / "out"),
            ],
            capture_output=True,
            check=False,
        )
        require(proc.returncode == 0, proc.stdout.decode() + proc.stderr.decode())
        require(
            json.loads((folder / "out/result.json").read_bytes()) == a,
            "installed CLI/API identity",
        )
        require(
            subprocess.run(
                [cli, "verify", str(folder / "out")],
                capture_output=True,
                check=False,
            ).returncode
            == 0,
            "installed verify",
        )
        for command in [
            ["describe", "tau-result"],
            ["analyze-v1", str(folder / "model.json")],
            ["export", str(folder / "model.json"), "gillespy2"],
        ]:
            extra = subprocess.run([cli, *command], capture_output=True, check=False)
            require(extra.returncode == 0, "installed discovery/legacy/export parity")
        tau_request = {
            "schema": "crn-craft.tau-request/0.1",
            "epsilon": 0.05,
            "max_leaps_per_replicate": 100000,
            "reference_run": r,
        }
        (folder / "tau.json").write_bytes(c.model_bytes(tau_request))
        extra = subprocess.run(
            [cli, "tau", str(folder / "model.json"), str(folder / "tau.json")],
            capture_output=True,
            check=False,
        )
        require(
            extra.returncode == 0
            and json.loads(extra.stdout)["assessment"] == "approximate_ensemble",
            "installed tau CLI parity",
        )
        (folder / "decay.xml").write_bytes(SBML)
        imported = subprocess.run(
            [
                importer,
                "sbml",
                str(folder / "decay.xml"),
                str(folder / "decay.json"),
                "--convention",
                "stochastic_count_propensity",
            ],
            capture_output=True,
            check=False,
        )
        require(imported.returncode == 0, imported.stdout.decode())
        decay = json.loads((folder / "decay.json").read_bytes())["model"]
        require(
            c.check(decay)["valid"]
            and decay["reactions"][0]["rate"]["coefficient"] == 0.5,
            "installed SBML import script",
        )
        checks.append("installed_CLI_archive_and_import_scripts")
    fsp = {k: v for k, v in r.items() if k != "ensemble"}
    fsp["method"] = {
        "kind": "fsp",
        "bounds": {"A": 100},
        "max_states": 101,
        "tolerance": 1e-12,
    }
    fsp["outputs"] = {"statistics": ["mean"], "keep_trajectories": 0}
    fsp["budget"] = {
        k: v for k, v in r["budget"].items() if k != "max_events_per_replicate"
    }
    result = c.fsp(m, fsp)
    require(not result["numerical_error_certified"], "installed FSP semantics")
    # Immigration-death from zero is Poisson with mean 20 (1 - exp(-0.2 t)).
    final = result["snapshots"][-1]
    mean = 20.0 * (1.0 - math.exp(-0.4))
    require(
        abs(final["projected_first_moment"][0] - mean) < 1e-9
        and abs(final["projected_second_raw_moment"][0] - mean - mean**2) < 1e-8,
        "installed FSP matches the Poisson moments",
    )
    checks.append("FSP")
    approx = c.tau(
        m,
        {
            "schema": "crn-craft.tau-request/0.1",
            "epsilon": 0.05,
            "max_leaps_per_replicate": 100000,
            "reference_run": r,
        },
    )
    require(
        approx["assessment"] == "approximate_ensemble",
        "installed approximate assessment",
    )
    checks.append("tau")
    net = b"begin species\n1 A() 2\nend species\nbegin reactions\n1 1 0 1.0\nend reactions\n"
    require(
        c.check(import_bng_net(net, rate_convention="binomial")["model"])["valid"],
        "installed importer",
    )
    from crn_craft.importers import ImportFailure

    try:
        import_bng_net(net.replace(b"1.0", b"1e-400"), rate_convention="binomial")
    except ImportFailure:
        pass
    else:
        raise RuntimeError(
            "installed importer silently underflowed a positive coefficient"
        )
    checks.append("importer")
    require(
        sampled_first_passage([0, 1], [[0, 0]], 1)["not_detected"] == 1,
        "installed observables",
    )
    t = np.linspace(0, 100, 501)
    require(
        oscillation(t, np.sin(2 * np.pi * t / 10))["classification"] == "oscillatory"
        and oscillation(t, np.full(len(t), 0.3))["classification"] == "nonoscillatory",
        "installed oscillation detector",
    )
    try:
        c.fsp(m, dict(fsp, schema="crn-craft.run/9"))
    except c.CrnError as error:
        require(error.diagnostics[0]["code"] == "CRN-E003", "installed CrnError")
    else:
        raise RuntimeError("installed FSP accepted an invalid run")
    checks.append("observables")
    # The rules engine (`crn_craft.rules`) ships in the wheel; the mtDNA manual
    # depends on it. Immigration-death has the Poisson law with mean
    # (k_in/k_out)(1 - exp(-k_out t)), so the means at t = 5 and 10 must lie
    # within 5 standard errors, and threads must not change any result.
    death = {
        "schema": "crn-craft.rules-model/0.1",
        "name": "immigration_death",
        "parameters": {"k_in": 10.0, "k_out": 1.0},
        "types": [{"name": "X", "attributes": []}],
        "initial": [],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {
                "name": "immigration",
                "reactants": [],
                "rate": "k_in",
                "products": [{"type": "X"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "X"}],
                "rate": "k_out",
                "products": [],
            },
        ],
    }

    def rules_run(threads):
        return rules.ensemble(
            death,
            end=10.0,
            report_times=[0.0, 5.0, 10.0],
            replicates=400,
            seed=7,
            threads=threads,
            method="static_direct",
            trackers=[{"observable": "X"}],
        ).simulate()

    serial, parallel = rules_run(1), rules_run(2)
    mean = serial.statistics["mean"]["X"].to_numpy()
    error = serial.statistics["mean_standard_error"]["X"].to_numpy()
    exact = 10.0 * (1.0 - np.exp(-np.array([0.0, 5.0, 10.0])))
    require(
        mean[0] == 0.0 and bool(np.all(np.abs(mean[1:] - exact[1:]) <= 5 * error[1:])),
        "installed rules ensemble matches the immigration-death law",
    )
    require(
        serial.records.equals(parallel.records)
        and serial.rule_events.equals(parallel.rule_events)
        and serial.statistics["mean"].equals(parallel.statistics["mean"])
        and set(serial.tracker_columns(0).columns)
        >= {"appearances", "first_passage", "integral"},
        "installed rules results are independent of threads and expose trackers",
    )
    # Results are pandas objects: records and rule counts by replicate,
    # statistics by report time, and every path with keep_paths="all".
    records = serial.records
    require(
        isinstance(records, pd.DataFrame)
        and records.index.name == "replicate"
        and list(records.index) == list(range(400))
        and isinstance(records["termination"].dtype, pd.CategoricalDtype)
        and bool((records["termination"] == "completed").all())
        and records["events"].dtype == np.uint64
        and serial.rule_names == ("immigration", "death")
        and list(serial.rule_events.columns) == ["immigration", "death"]
        and int(serial.rule_events.to_numpy().sum()) == int(records["events"].sum())
        and list(serial.statistics["mean"].index) == [0.0, 5.0, 10.0]
        and serial.statistics["samples"] == 400,
        "installed rules results are DataFrames",
    )
    kept = rules.ensemble(
        death,
        end=10.0,
        report_times=[0.0, 5.0, 10.0],
        replicates=50,
        seed=7,
        keep_paths="all",
    ).simulate()
    by_time = kept.paths["X"].groupby(level="time").mean().to_numpy()
    require(
        kept.paths.shape == (150, 1)
        and list(kept.paths.index.names) == ["replicate", "time"]
        and np.allclose(by_time, kept.statistics["mean"]["X"], rtol=1e-12, atol=0),
        "installed rules keep every path as a DataFrame",
    )
    try:
        import xarray  # noqa: F401  (optional: crn-craft[xarray])
    except ImportError:
        try:
            kept.to_xarray()
        except ImportError as error:
            require("crn-craft[xarray]" in str(error), "to_xarray names the extra")
        else:
            raise RuntimeError("to_xarray without xarray")
    else:
        data = kept.to_xarray()
        require(
            data["X"].dims == ("replicate", "time") and data["X"].shape == (50, 3),
            "installed to_xarray dimensions",
        )
    require(
        kept.paths.equals(
            rules.ensemble(
                death,
                end=10.0,
                report_times=[0.0, 5.0, 10.0],
                replicates=50,
                seed=7,
                keep_paths=True,
            )
            .simulate()
            .paths
        ),
        "keep_paths=True keeps every path, as 'all' does",
    )
    check_pandas_results(serial)
    require(
        {"mtdna_direct", "hybrid"} <= set(rules.describe()["methods"]),
        "installed rules describe lists mtdna_direct and hybrid",
    )
    checks.append("rules_engine")
    check_event_log_default_budget()
    checks.append("rules_event_log_default_budget")
    check_builder()
    checks.append("rules_builder")
    check_hybrid()
    checks.append("rules_hybrid")
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.installed-validation/0.1",
                "passed": True,
                "checks": checks,
                "package_path": str(Path(c.__file__).resolve()),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "checks": len(checks)}))


if __name__ == "__main__":
    main()
