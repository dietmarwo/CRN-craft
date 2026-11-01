"""Acceptance from an installed wheel, independent of repository imports."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import crn_craft as c
from crn_craft.importers import import_bng_net
from crn_craft.observables import sampled_first_passage


def require(v, m):
    if not v:
        raise RuntimeError(m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    checks = []
    m = c.load_example("immigration-death")
    r = c.run_spec(m, end=2.0, replicates=32, seed=77, keep=32)
    a = c.simulate(m, r)
    b = c.simulate(m, r, cached=True)
    require(a == b, "installed cached identity")
    checks.append("native_API_exact_identity")
    require(
        c.describe()["schema"] == "crn-craft.describe/0.5", "installed capabilities"
    )
    checks.append("discovery")
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        (folder / "model.json").write_bytes(m)
        (folder / "run.json").write_bytes(c.model_bytes(r))
        proc = subprocess.run(
            [
                "crn-craft",
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
                ["crn-craft", "verify", str(folder / "out")],
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
            extra = subprocess.run(
                ["crn-craft", *command], capture_output=True, check=False
            )
            require(extra.returncode == 0, "installed discovery/legacy/export parity")
        tau_request = {
            "schema": "crn-craft.tau-request/0.1",
            "epsilon": 0.05,
            "max_leaps_per_replicate": 100000,
            "reference_run": r,
        }
        (folder / "tau.json").write_bytes(c.model_bytes(tau_request))
        extra = subprocess.run(
            ["crn-craft", "tau", str(folder / "model.json"), str(folder / "tau.json")],
            capture_output=True,
            check=False,
        )
        require(
            extra.returncode == 0
            and json.loads(extra.stdout)["assessment"] == "approximate_ensemble",
            "installed tau CLI parity",
        )
        checks.append("installed_CLI_and_archive")
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
        sampled_first_passage([0, 1], [[0, 0]], 1)["censored"] == 1,
        "installed observables",
    )
    checks.append("observables")
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
