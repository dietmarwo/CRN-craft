"""The working directory the tutorial's commands run in.

The models are the packaged `crn-craft.model/1.0` examples; the requests are
`crn-craft.run/1.0` documents written here, each naming its model by file name
and SHA-256.
"""

import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
BUDGET = {
    "max_events_per_replicate": 1_000_000,
    "max_sampling_steps_per_replicate": 8_000_000,
    "max_instances": 100_000,
    "wall_time_s": 60.0,
    "threads": 1,
    "max_output_bytes": 1_048_576,
}
SAVE_TIMES = [0.0, 1.0, 2.0, 5.0, 10.0]


# An SBML Level 3 Version 2 model of the dimerization 2 A -> B.
SBML_DIMER = (
    '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" '
    'version="2">\n'
    '  <model id="Dimer" substanceUnits="item" timeUnits="second">\n'
    "    <listOfCompartments>\n"
    '      <compartment id="cell" size="1" constant="true"/>\n'
    "    </listOfCompartments>\n"
    "    <listOfSpecies>\n"
    '      <species id="A" compartment="cell" initialAmount="20"'
    ' hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>\n'
    '      <species id="B" compartment="cell" initialAmount="0"'
    ' hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>\n'
    "    </listOfSpecies>\n"
    "    <listOfParameters>\n"
    '      <parameter id="k" value="0.1" constant="true"/>\n'
    "    </listOfParameters>\n"
    "    <listOfReactions>\n"
    '      <reaction id="dimerize" reversible="false">\n'
    "        <listOfReactants>\n"
    '          <speciesReference species="A" stoichiometry="2" constant="true"/>\n'
    "        </listOfReactants>\n"
    "        <listOfProducts>\n"
    '          <speciesReference species="B" stoichiometry="1" constant="true"/>\n'
    "        </listOfProducts>\n"
    "        <kineticLaw>\n"
    '          <math xmlns="http://www.w3.org/1998/Math/MathML">\n'
    "            <apply><divide/>\n"
    "              <apply><times/><ci>k</ci><ci>A</ci>\n"
    "                <apply><minus/><ci>A</ci><cn>1</cn></apply></apply>\n"
    "              <cn>2</cn></apply>\n"
    "          </math>\n"
    "        </kineticLaw>\n"
    "      </reaction>\n"
    "    </listOfReactions>\n"
    "  </model>\n"
    "</sbml>\n"
)


def request(workdir, name, **changes):
    """A `run/1.0` Static Direct request of the model file `name` in `workdir`."""
    digest = hashlib.sha256((workdir / f"{name}.model.json").read_bytes()).hexdigest()
    run = {
        "schema": "crn-craft.run/1.0",
        "model": {"path": f"{name}.model.json", "sha256": digest},
        "method": {"kind": "static_direct"},
        "time": {"end": 10.0, "report_times": SAVE_TIMES},
        "ensemble": {"replicates": 64, "root_seed": 42},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 2,
            "keep_instances": 0,
        },
        "budget": dict(BUDGET),
    }
    run.update(changes)
    return run


def write(path, document):
    path.write_text(json.dumps(document, indent=2) + "\n")


def prepare_workdir(workdir):
    for path in (ROOT / "python/crn_craft/models").glob("*.model.json"):
        shutil.copy2(path, workdir / path.name)
    models = sorted(p.name.removesuffix(".model.json") for p in workdir.glob("*.json"))
    for name in models:
        write(workdir / f"{name}.run.json", request(workdir, name))
    # A model whose first rule names a parameter that does not exist.
    broken = json.loads((workdir / "immigration-death.model.json").read_text())
    broken["rules"][0]["rate"] = "missing"
    write(workdir / "missing-parameter.model.json", broken)
    run = request(workdir, "immigration-death")
    stale = dict(run, model=dict(run["model"], sha256="0" * 64))
    write(workdir / "stale-hash.run.json", stale)
    long = request(
        workdir,
        "immigration-death",
        time={"end": 1e6, "report_times": [0.0, 1e6]},
    )
    long["budget"] = dict(BUDGET, max_events_per_replicate=10**12)
    long["budget"].pop("wall_time_s")
    write(workdir / "long.run.json", long)
    fsp = request(
        workdir,
        "immigration-death",
        method={
            "kind": "fsp",
            "bounds": {"A": 100},
            "max_states": 102,
            "tolerance": 1e-10,
        },
    )
    fsp.pop("ensemble")
    fsp.pop("outputs")
    fsp["budget"] = {"wall_time_s": 60.0, "max_output_bytes": 1_048_576}
    write(workdir / "immigration-death.fsp.run.json", fsp)
    tau = request(workdir, "immigration-death", method={"kind": "tau", "epsilon": 0.03})
    write(workdir / "tau.run.json", tau)
    (workdir / "dimer.net").write_text(
        "begin parameters\n1 k 0.1\nend parameters\n"
        "begin species\n1 A() 20\n2 B() 0\nend species\n"
        "begin reactions\n1 1,1 2 k\nend reactions\n"
    )
    (workdir / "dimer.xml").write_text(SBML_DIMER)
    # The same model in the network format that release 0.6.0 wrote.
    shutil.copy2(
        ROOT / "models/legacy/immigration-death.model.json",
        workdir / "immigration-death.network.json",
    )
    (workdir / "runs").mkdir()
