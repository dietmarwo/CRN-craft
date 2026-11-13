"""The working directories that the documentation's commands run in.

`prepare_guide_workdir` writes the files of the guide's command-line chapters
(`docs/guide/10-command-line.md` and `11-analysis-import-export.md`). The
models are built with the same builder lines as the chapters and written
with `Model.to_json()`, so the hashes in the screenshots are those a reader
gets. The requests are `crn-craft.run/1.0` documents, each naming its model
by file name and SHA-256. `examples/tutorial/make_assets.py` runs its
terminal sessions there. `validation/check_guide.py` runs the chapters
themselves, each in an empty directory that its own blocks fill.
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


# An SBML Level 3 Version 2 model of the dimerization 2 A -> B, and the same
# reaction as an expanded BioNetGen network.
SBML_DIMER = ROOT / "examples/guide/dimer.xml"
BNG_DIMER = (
    "begin parameters\n1 k 0.1\nend parameters\n"
    "begin species\n1 A() 20\n2 B() 0\nend species\n"
    "begin reactions\n1 1,1 2 k\nend reactions\n"
)


def immigration_death():
    """The guide's running example, with the builder lines of its chapters."""
    import crn_craft as crn

    model = crn.Model("immigration_death")
    model.add_parameter({"birth": 4.0, "death": 0.2})
    model.add_species("A", 0)
    model.add_reaction("birth", [], ["A"], name="birth")
    model.add_reaction("death * A", ["A"], [], name="death")
    return model


def enzyme_mechanism():
    """S + E <-> ES -> E + P, as chapter 11 builds it."""
    import crn_craft as crn

    enzyme = crn.Model("enzyme_mechanism")
    enzyme.add_species({"S": 50, "E": 10, "ES": 0, "P": 0})
    enzyme.add_reaction(0.01, ["S", "E"], ["ES"], name="bind")
    enzyme.add_reaction(0.1, ["ES"], ["S", "E"], name="unbind")
    enzyme.add_reaction(0.5, ["ES"], ["E", "P"], name="catalysis")
    return enzyme


def repressilator():
    """Three genes repressing each other in a ring, as chapter 11 builds it."""
    import crn_craft as crn

    circuit = crn.Model("repressilator")
    circuit.add_parameter({"alpha0": 0.216, "alpha": 216.0, "K": 40.0, "n": 2.0})
    circuit.add_parameter({"beta": 5.0, "gamma_m": 1.0, "gamma_p": 0.2})
    circuit.add_species("p1", 5)
    for gene, repressor in ((1, 3), (2, 1), (3, 2)):
        mrna, protein = f"m{gene}", f"p{gene}"
        hill = f"alpha0 + alpha / (1 + (p{repressor} / K)^n)"
        circuit.add_reaction(hill, [], [mrna], name=f"tx{gene}")
        circuit.add_reaction(
            f"beta * {mrna}", [mrna], [mrna, protein], name=f"tl{gene}"
        )
        circuit.add_reaction(f"gamma_m * {mrna}", [mrna], [], name=f"dm{gene}")
        circuit.add_reaction(f"gamma_p * {protein}", [protein], [], name=f"dp{gene}")
    return circuit


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


def prepare_guide_workdir(workdir):
    """Write the files that chapters 10 and 11 of the guide write or copy."""
    import crn_craft as crn

    name = "immigration-death"
    (workdir / f"{name}.model.json").write_text(immigration_death().to_json())
    run = request(workdir, name)
    write(workdir / f"{name}.run.json", run)
    # A model whose first rate names a parameter that does not exist.
    broken = crn.Model("immigration_death")
    broken.add_species("A", 0)
    broken.add_reaction("missing", [], ["A"], name="birth")
    broken.add_reaction(0.2, ["A"], [], name="death")
    (workdir / "missing-parameter.model.json").write_text(broken.to_json())
    stale = dict(run, model=dict(run["model"], sha256="0" * 64))
    write(workdir / "stale-hash.run.json", stale)
    write(workdir / "tau.run.json", dict(run, method={"kind": "tau", "epsilon": 0.03}))
    fsp = {
        "schema": run["schema"],
        "model": run["model"],
        "method": {
            "kind": "fsp",
            "bounds": {"A": 100},
            "max_states": 102,
            "tolerance": 1e-10,
        },
        "time": run["time"],
        "budget": {"wall_time_s": 60.0, "max_output_bytes": 1_048_576},
    }
    write(workdir / f"{name}.fsp.run.json", fsp)
    long = request(workdir, name, time={"end": 1e6, "report_times": [0.0, 1e6]})
    long["budget"]["max_events_per_replicate"] = 10**12
    del long["budget"]["wall_time_s"]
    write(workdir / "long.run.json", long)
    (workdir / "runs").mkdir()
    # Chapter 11: one model file written from the builder, one written by a
    # block, and two repository files under their repository paths.
    (workdir / "enzyme-mechanism.model.json").write_text(enzyme_mechanism().to_json())
    (workdir / "dimer.net").write_text(BNG_DIMER)
    for fixture in (
        "models/legacy/immigration-death.model.json",
        "examples/guide/dimer.xml",
    ):
        (workdir / fixture).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / fixture, workdir / fixture)
