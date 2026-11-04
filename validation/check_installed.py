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
    mean = serial.statistics["mean"][:, 0]
    error = serial.statistics["mean_standard_error"][:, 0]
    exact = 10.0 * (1.0 - np.exp(-np.array([0.0, 5.0, 10.0])))
    require(
        mean[0] == 0.0 and bool(np.all(np.abs(mean[1:] - exact[1:]) <= 5 * error[1:])),
        "installed rules ensemble matches the immigration-death law",
    )
    require(
        all(
            np.array_equal(serial.records[k], parallel.records[k])
            for k in serial.records
        )
        and np.array_equal(serial.statistics["mean"], parallel.statistics["mean"])
        and set(serial.tracker_columns(0))
        >= {"appearances", "first_passage", "integral"},
        "installed rules results are independent of threads and expose trackers",
    )
    require(
        "mtdna_direct" in str(rules.describe()["methods"]),
        "installed rules describe lists mtdna_direct",
    )
    checks.append("rules_engine")
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
