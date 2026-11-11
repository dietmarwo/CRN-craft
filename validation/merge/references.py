"""Record the 0.6.0 network engine's laws for the 0.7.0 engine merge (task M0c).

For every model of `validation/merge/models.json` this writes, into the output
directory (default `validation/merge/reference/`):

- `<name>.moments.json`: mean, variance and their standard errors per species
  and report time from the declared replicates of the `direct` law (computed
  with `cached_direct`, whose paths are identical; checked per model), with the
  exact integer power sums they come from;
- `<name>.analysis.json.gz`: the output of `crn-craft analyze MODEL`;
- `<name>.kinetics.json`: the `kinetics` right-hand side and Jacobians at the
  initial and a probe state, or the refusal;
- for models with declared FSP bounds: `<name>.fsp.json.gz` (the full
  `crn-craft fsp` result), `<name>.fsp-moments.json` (the law's moments) and
  `<name>.self-check.json` (the SSA reference against the FSP law);
- `index.json`: SHA-256 of every file in the directory and the tool identities.

The declarations are in `validation/merge/RESULTS.md`. This process needs
only the standard library; the simulations run in a worker started with
`--python` (an environment with the 0.6.0 `crn_craft` wheel and numpy).
Run it under the CPU lock, for example:

    flock /home/claude/craft5/cpu.lock python validation/merge/references.py \
        --python /home/claude/m07/venv060/bin/python \
        --binary /home/claude/m07/bin060/crn-craft \
        --kinetics-binary /home/claude/m07/bin060/kinetics-reference

`--check-against DIR` compares the files written for the chosen models with
those in DIR byte for byte and exits 1 on a difference.
"""

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import (  # noqa: E402
    MANIFEST,
    ROOT,
    binary_identity,
    compare_to_law,
    encode,
    gzip_bytes,
    load_manifest,
    run_worker,
    sha256,
    ssa_request,
    statistics_from_power_sums,
    summarize,
)

DEFAULT_OUTPUT = ROOT / "validation/merge/reference"
REFERENCE_METHOD = {"aggregator": "cached_direct", "waiting_time": "inverse"}
IDENTITY_REPLICATES = 1000
FSP_ERROR_LIMIT = 1e-9


def dump(value):
    return (json.dumps(value, separators=(",", ":"), allow_nan=False) + "\n").encode()


# ---------------------------------------------------------------- worker side


def power_sums(x, np):
    """Exact sums of x, x^2, x^3, x^4 over replicates as Python-int object arrays."""
    top = int(x.max()) if x.size else 0
    sums = [np.zeros(x.shape[1:], dtype=object) for _ in range(4)]
    if top**4 >= 2**63:
        big = x.astype(object)
        for p in range(4):
            sums[p] += (big ** (p + 1)).sum(axis=0)
        return sums
    block = (2**63 - 1) // max(1, top**4)
    for start in range(0, x.shape[0], block):
        part = x[start : start + block]
        for p in range(4):
            sums[p] += (part ** (p + 1)).sum(axis=0, dtype=np.int64).astype(object)
    return sums


def worker_moments(job):
    import crn_craft as crn
    import numpy as np

    model = (ROOT / job["path"]).read_bytes()
    chunk = job["chunk"]
    times = job["report_times"]
    totals = None
    chunks = []
    sampler = None
    species = None
    identity = None

    def request(aggregator, replicates, seed):
        return ssa_request(
            model,
            aggregator=aggregator,
            waiting_time="inverse",
            fast=False,
            end=job["end"],
            save_times=times,
            replicates=replicates,
            seed=seed,
            threads=job["threads"],
            keep=replicates,
        )

    def completed(result, replicates):
        if (
            result["termination"] != "completed"
            or result["completed_replicates"] != replicates
            or len(result["trajectories"]) != replicates
            or [t["id"] for t in result["trajectories"]] != list(range(replicates))
        ):
            raise RuntimeError(
                f"{job['name']}: incomplete chunk {result['diagnostics']}"
            )
        return np.asarray([t["counts"] for t in result["trajectories"]], dtype=np.int64)

    for k in range(job["replicates"] // chunk):
        seed = job["reference_seed"] + k
        result = crn.simulate_ssa(
            model, request(REFERENCE_METHOD["aggregator"], chunk, seed)
        )
        x = completed(result, chunk)
        if k == 0:
            sampler = result["sampler"]
            species = result["species"]
            m = min(chunk, IDENTITY_REPLICATES)
            direct = crn.simulate_ssa(model, request("direct", m, seed))
            same = np.array_equal(completed(direct, m), x[:m]) and [
                r["events"] for r in direct["replicates"]
            ] == [r["events"] for r in result["replicates"][:m]]
            if not same:
                raise RuntimeError(
                    f"{job['name']}: cached_direct paths differ from direct"
                )
            identity = {
                "aggregator": "direct",
                "root_seed": seed,
                "replicates": m,
                "sampler": direct["sampler"],
                "identical": True,
            }
        if result["sampler"] != sampler or result["species"] != species:
            raise RuntimeError(f"{job['name']}: sampler or species changed")
        # Sanity check of the decoding against the engine's own statistics.
        for mine, theirs in [
            (x.mean(axis=0), result["statistics"]["mean"]),
            (x.var(axis=0, ddof=1), result["statistics"]["sample_variance"]),
        ]:
            theirs = np.asarray(theirs)
            if not np.all(
                np.abs(mine - theirs) <= 1e-8 * np.maximum(1.0, np.abs(theirs))
            ):
                raise RuntimeError(f"{job['name']}: statistics disagree with the paths")
        sums = power_sums(x, np)
        totals = (
            sums
            if totals is None
            else [a + b for a, b in zip(totals, sums, strict=True)]
        )
        chunks.append(
            {
                "root_seed": seed,
                "events": result["events"],
                "rng_words": result["work"]["rng_words"],
            }
        )
        print(
            f"  {job['name']} chunk {k} seed {seed} events {result['events']}",
            file=sys.stderr,
            flush=True,
        )
    return {
        "species": species,
        "sampler": sampler,
        "direct_identity": identity,
        "chunks": chunks,
        "power_sums": {
            "n": job["replicates"],
            **{f"s{p + 1}": totals[p].tolist() for p in range(4)},
        },
    }


def worker_identity():
    import crn_craft as crn
    import numpy as np

    return {
        "python": sys.version,
        "crn_craft": crn.__version__,
        "numpy": np.__version__,
    }


def worker(job_path, answer_path):
    job = json.loads(Path(job_path).read_text())
    if job["kind"] == "identity":
        answer = worker_identity()
    else:
        answer = worker_moments(job)
    Path(answer_path).write_text(json.dumps(answer))


# ----------------------------------------------------------------- outer side


def fsp_run(row, model_bytes):
    return {
        "schema": "crn-craft.run/0.1",
        "model": {"path": "model.json", "sha256": sha256(model_bytes)},
        "method": {"kind": "fsp", **row["fsp"]},
        "time": {"end": row["end"], "save_times": row["report_times"]},
        "outputs": {"statistics": ["mean", "variance"], "keep_trajectories": 0},
        "budget": {"threads": 1, "max_output_bytes": 200_000_000},
    }


def fsp_moments(result, row):
    """Law moments of each FSP snapshot, normalized by the retained mass."""
    states = result["states"]
    species = result["species"]
    out = {
        "mean": [],
        "variance": [],
        "fourth_central_moment": [],
        "retained_mass": [],
        "physical_leakage_estimate": [],
        "mass_defect": [],
        "series_tail_upper_bound": [],
    }
    for snapshot in result["snapshots"]:
        p = snapshot["probabilities"]
        mass = math.fsum(p)
        error = snapshot["physical_leakage_estimate"] + abs(snapshot["mass_defect"])
        if error > FSP_ERROR_LIMIT:
            raise RuntimeError(
                f"{row['name']}: FSP error {error} at {snapshot['time']}"
            )
        means, variances, fourths = [], [], []
        for j in range(len(species)):
            mean = math.fsum(q * x[j] for q, x in zip(p, states, strict=True)) / mass
            means.append(mean)
            variances.append(
                math.fsum(
                    q * (x[j] - mean) ** 2 for q, x in zip(p, states, strict=True)
                )
                / mass
            )
            fourths.append(
                math.fsum(
                    q * (x[j] - mean) ** 4 for q, x in zip(p, states, strict=True)
                )
                / mass
            )
        out["mean"].append(means)
        out["variance"].append(variances)
        out["fourth_central_moment"].append(fourths)
        out["retained_mass"].append(mass)
        for key in [
            "physical_leakage_estimate",
            "mass_defect",
            "series_tail_upper_bound",
        ]:
            out[key].append(snapshot[key])
    return out


def run_cli(binary, arguments, cwd):
    done = subprocess.run([str(binary), *arguments], capture_output=True, cwd=cwd)
    return done.returncode, done.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--python", type=Path, help="Python with the 0.6.0 wheel")
    parser.add_argument("--binary", type=Path, help="0.6.0 crn-craft CLI")
    parser.add_argument(
        "--kinetics-binary", type=Path, help="kinetics-reference printer"
    )
    parser.add_argument("--models", help="comma-separated model names (default: all)")
    parser.add_argument(
        "--manifest", type=Path, default=MANIFEST, help=argparse.SUPPRESS
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--check-against", type=Path)
    parser.add_argument("--worker", nargs=2, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
        return 0
    if not (args.python and args.binary and args.kinetics_binary):
        parser.error("--python, --binary and --kinetics-binary are required")
    binary = args.binary.resolve()
    rows = load_manifest(args.models.split(",") if args.models else None, args.manifest)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    written = []

    def write(name, data):
        (out / name).write_bytes(data)
        written.append(name)

    kinetics = subprocess.run(
        [str(args.kinetics_binary.resolve()), *(row["path"] for row in rows)],
        capture_output=True,
        check=True,
        cwd=ROOT,
    ).stdout
    kinetics = json.loads(kinetics)
    by_path = {entry["model"]: entry for entry in kinetics["models"]}
    with tempfile.TemporaryDirectory(prefix="crn-merge-ref-") as scratch:
        scratch = Path(scratch)
        identity = run_worker(args.python, __file__, {"kind": "identity"}, scratch)
        for row in rows:
            name = row["name"]
            print(f"{name}: analysis, kinetics, FSP", file=sys.stderr, flush=True)
            model_bytes = (ROOT / row["path"]).read_bytes()
            code, analysis = run_cli(binary, ["analyze", row["path"]], ROOT)
            if code != 0:
                raise RuntimeError(f"{name}: analyze exited {code}")
            write(f"{name}.analysis.json.gz", gzip_bytes(analysis))
            write(
                f"{name}.kinetics.json",
                dump(
                    {
                        "schema": kinetics["schema"],
                        "crn_craft": kinetics["crn_craft"],
                        "model": by_path[row["path"]],
                    }
                ),
            )
            law = None
            if row["fsp"]:
                folder = scratch / name
                folder.mkdir()
                (folder / "model.json").write_bytes(model_bytes)
                run = fsp_run(row, model_bytes)
                (folder / "run.json").write_bytes(encode(run))
                code, result = run_cli(
                    binary, ["fsp", "model.json", "run.json"], folder
                )
                if code != 0:
                    raise RuntimeError(f"{name}: fsp exited {code}: {result[:400]!r}")
                write(f"{name}.fsp.json.gz", gzip_bytes(result))
                law = fsp_moments(json.loads(result), row)
                write(
                    f"{name}.fsp-moments.json",
                    dump(
                        {
                            "schema": "crn-craft.merge-fsp-moments/0.1",
                            "model": {
                                "name": name,
                                "path": row["path"],
                                "sha256": row["sha256"],
                            },
                            "run": run,
                            "report_times": row["report_times"],
                            "species": json.loads(result)["species"],
                            **law,
                        }
                    ),
                )
            print(f"{name}: SSA reference", file=sys.stderr, flush=True)
            job = {"kind": "moments", **row, "threads": args.threads}
            answer = run_worker(args.python, __file__, job, scratch)
            stats = statistics_from_power_sums(answer["power_sums"])
            moments = {
                "schema": "crn-craft.merge-reference-moments/0.1",
                "model": {"name": name, "path": row["path"], "sha256": row["sha256"]},
                "law_of": "direct",
                "method": {"kind": "ssa", **REFERENCE_METHOD, "fast": False},
                "sampler": answer["sampler"],
                "direct_identity": answer["direct_identity"],
                "end": row["end"],
                "report_times": row["report_times"],
                "species": answer["species"],
                "replicates": row["replicates"],
                "chunk": row["chunk"],
                "chunks": answer["chunks"],
                "events": sum(c["events"] for c in answer["chunks"]),
                **stats,
                "power_sums": answer["power_sums"],
            }
            write(f"{name}.moments.json", dump(moments))
            if law is not None:
                rows_ = compare_to_law(
                    stats,
                    law,
                    row["replicates"],
                    row["report_times"],
                    answer["species"],
                )
                check = {
                    "schema": "crn-craft.merge-reference-self-check/0.1",
                    "model": name,
                    "declared": "RESULTS.md: FSP comparisons of the law gate, n = replicates",
                    **summarize(rows_),
                    "rows": rows_,
                }
                write(f"{name}.self-check.json", dump(check))
                print(
                    f"{name}: self-check {summarize(rows_)}",
                    file=sys.stderr,
                    flush=True,
                )
    files = {
        path.name: {"sha256": sha256(path.read_bytes()), "bytes": path.stat().st_size}
        for path in sorted(out.iterdir())
        if path.name != "index.json"
    }
    index = {
        "schema": "crn-craft.merge-reference-index/0.1",
        "models_json_sha256": sha256(args.manifest.read_bytes()),
        "tools": {
            "crn_craft_cli": binary_identity(binary),
            "kinetics_reference_sha256": sha256(args.kinetics_binary.read_bytes()),
            "worker": identity,
        },
        "files": files,
    }
    (out / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    if args.check_against:
        differ = [
            name
            for name in written
            if (args.check_against / name).read_bytes() != (out / name).read_bytes()
        ]
        print(
            json.dumps({"compared": len(written), "different": differ}, indent=2),
            flush=True,
        )
        return 1 if differ else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
