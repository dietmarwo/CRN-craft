"""References of the E1 law check: the builder form of the three-class model.

Declared in `validation/merge/RESULTS.md` ("E1 (`mtdna_direct` on the
builder form): law check declaration"). For every row of
`validation/merge/builder/models.json` this writes, into
`validation/merge/builder/reference/`:

- `<name>.moments.json`: mean, variance and their standard errors per
  species and report time from the row's replicates under `static_direct`
  (`crn-craft.run/1.0`) with the reference seeds, from the release CLI given
  with `--binary` (for E1 the CLI of `3d80c6d`, the 0.7.0 engine), with the
  exact integer power sums they come from;
- `<name>.fsp-moments.json`: the moments of the finite state projection of
  the same CLI (`fsp`, the row's bounds and tolerance), kept only if the
  leakage and the mass defect are below 1e-9 at every report time.

Standard library only; run under the CPU lock:

    flock /home/claude/craft5/cpu.lock python validation/merge/builder_references.py \\
        --binary /path/to/0.7.0/crn-craft
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_laws import chunk_plan, cli_model  # noqa: E402
from merge_common import (  # noqa: E402
    HERE,
    ROOT,
    binary_identity,
    encode,
    load_manifest,
    sha256,
    statistics_from_power_sums,
)
from references import fsp_moments  # noqa: E402

MANIFEST = HERE / "builder" / "models.json"
OUTPUT = HERE / "builder" / "reference"
METHOD = {"kind": "static_direct"}


def dump(value):
    return (json.dumps(value, indent=1, allow_nan=False) + "\n").encode()


def fsp_law(binary, row, folder):
    """The FSP result of the row's model through `crn-craft rules`."""
    model = (ROOT / row["path"]).read_bytes()
    (folder / "model.json").write_bytes(model)
    request = {
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": sha256(model)},
        "time": {"end": row["end"], "report_times": row["report_times"]},
        "method": {"kind": "fsp", **row["fsp"]},
        "budget": {"max_output_bytes": 1 << 30},
    }
    (folder / "request.json").write_bytes(encode(request))
    done = subprocess.run(
        [str(binary), "rules", "model.json", "request.json"],
        capture_output=True,
        cwd=folder,
        check=False,
    )
    result = json.loads(done.stdout)
    if done.returncode != 0 or result.get("schema") != "crn-craft.fsp-result/1.0":
        raise RuntimeError(f"{row['name']}: {result}")
    return request, result


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    binary = args.binary.resolve()
    identity = binary_identity(binary)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for row in load_manifest(None, MANIFEST):
        name = row["name"]
        species = row["species_names"]
        with tempfile.TemporaryDirectory(prefix="crn-builder-references-") as scratch:
            folder = Path(scratch)
            request, result = fsp_law(binary, row, folder)
            order = [result["species"].index(s) for s in species]
            reordered = dict(result)
            reordered["species"] = species
            reordered["states"] = [[x[i] for i in order] for x in result["states"]]
            law = fsp_moments(reordered, row)
            (OUTPUT / f"{name}.fsp-moments.json").write_bytes(
                dump(
                    {
                        "schema": "crn-craft.merge-fsp-moments/0.1",
                        "model": {
                            "name": name,
                            "path": row["path"],
                            "sha256": row["sha256"],
                        },
                        "engine": identity,
                        "run": request,
                        "report_times": row["report_times"],
                        "species": species,
                        "states": len(result["states"]),
                        **law,
                    }
                )
            )
            print(f"{name}: FSP {len(result['states'])} states", file=sys.stderr)
            chunk, chunks = chunk_plan(row, row["replicates"])
            job = {
                **row,
                "species": species,
                "chunk": chunk,
                "chunks": chunks,
                "threads": args.threads,
                "method": METHOD,
                # The reference seeds: chunk k uses reference_seed + k.
                "gate_seed": row["reference_seed"],
            }
            (folder / "sample").mkdir()
            answer = cli_model(binary, job, folder / "sample")
        stats = statistics_from_power_sums(answer["power_sums"])
        (OUTPUT / f"{name}.moments.json").write_bytes(
            dump(
                {
                    "schema": "crn-craft.merge-reference-moments/0.1",
                    "model": {
                        "name": name,
                        "path": row["path"],
                        "sha256": row["sha256"],
                    },
                    "engine": identity,
                    "method": METHOD,
                    "request_schema": answer["request_schema"],
                    "sampler": answer["sampler"],
                    "end": row["end"],
                    "report_times": row["report_times"],
                    "species": species,
                    "replicates": row["replicates"],
                    "chunk": chunk,
                    "chunks": answer["chunks"],
                    "events": sum(c["events"] for c in answer["chunks"]),
                    **{
                        k: stats[k]
                        for k in [
                            "mean",
                            "variance",
                            "mean_se",
                            "variance_se",
                            "fourth_central_moment",
                        ]
                    },
                    "power_sums": answer["power_sums"],
                }
            )
        )
        print(f"{name}: {answer['sampler']}", file=sys.stderr)


if __name__ == "__main__":
    main()
