"""The working directory the tutorial's commands run in."""

import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def prepare_workdir(workdir):
    for path in (ROOT / "models/examples").iterdir():
        shutil.copy2(path, workdir / path.name)
    shutil.copy2(ROOT / "models/invalid/missing-parameter.model.json", workdir)
    run = json.loads((workdir / "immigration-death.run.json").read_text())
    stale = dict(run, model=dict(run["model"], sha256="0" * 64))
    (workdir / "stale-hash.run.json").write_text(json.dumps(stale, indent=2) + "\n")
    long = dict(run, time={"end": 1e6, "save_times": [0.0, 1e6]})
    long["budget"] = dict(run["budget"], max_events_per_replicate=10**12)
    long["budget"].pop("wall_time_s")
    (workdir / "long.run.json").write_text(json.dumps(long, indent=2) + "\n")
    tau = {
        "schema": "crn-craft.tau-request/0.1",
        "epsilon": 0.03,
        "max_leaps_per_replicate": 100000,
        "reference_run": run,
    }
    (workdir / "tau.json").write_text(json.dumps(tau, indent=2) + "\n")
    (workdir / "dimer.net").write_text(
        "begin parameters\n1 k 0.1\nend parameters\n"
        "begin species\n1 A() 20\n2 B() 0\nend species\n"
        "begin reactions\n1 1,1 2 k\nend reactions\n"
    )
    (workdir / "runs").mkdir()
