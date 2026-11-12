"""FSP references of the higher-order exact-law tests (task M4b).

For every row of `validation/merge/higher-order/models.json` this runs the
0.6.0 network engine's `crn-craft fsp` (`run/0.1`, method `fsp`, the row's
bounds, `max_states` and tolerance, at the row's report times) and writes,
into `validation/merge/reference/higher-order/`:

- `<name>.fsp.json.gz`: the full FSP result (states and probabilities);
- `<name>.fsp-moments.json`: the law's mean, variance and fourth central
  moment per species and report time, normalized by the retained mass, with
  the same formulas and error limit (leakage + mass defect <= 1e-9) as the
  M0c references;
- `index.json`: SHA-256 of every file and the identity of the CLI.

The declarations are in `validation/merge/RESULTS.md` (M4b). This script
needs only the standard library and the kept 0.6.0 CLI:

    flock /home/claude/craft5/cpu.lock python \\
        validation/merge/higher_order_references.py \\
        --binary /home/claude/m07/bin060/crn-craft
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import (  # noqa: E402
    HERE,
    ROOT,
    binary_identity,
    encode,
    gzip_bytes,
    load_manifest,
    sha256,
)
from references import dump, fsp_moments, fsp_run, run_cli  # noqa: E402

MANIFEST = HERE / "higher-order" / "models.json"
OUTPUT = HERE / "reference" / "higher-order"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--binary", type=Path, required=True, help="0.6.0 CLI")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    binary = args.binary.resolve()
    rows = load_manifest(path=MANIFEST)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="crn-merge-ho-") as scratch:
        for row in rows:
            name = row["name"]
            model_bytes = (ROOT / row["path"]).read_bytes()
            folder = Path(scratch) / name
            folder.mkdir()
            (folder / "model.json").write_bytes(model_bytes)
            run = fsp_run(row, model_bytes)
            (folder / "run.json").write_bytes(encode(run))
            code, result = run_cli(binary, ["fsp", "model.json", "run.json"], folder)
            if code != 0:
                raise RuntimeError(f"{name}: fsp exited {code}: {result[:400]!r}")
            (out / f"{name}.fsp.json.gz").write_bytes(gzip_bytes(result))
            parsed = json.loads(result)
            law = fsp_moments(parsed, row)
            (out / f"{name}.fsp-moments.json").write_bytes(
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
                        "species": parsed["species"],
                        **law,
                    }
                )
            )
            print(
                f"{name}: {len(parsed['states'])} states, largest error "
                f"{max(a + abs(b) for a, b in zip(law['physical_leakage_estimate'], law['mass_defect'], strict=True)):.2e}",
                file=sys.stderr,
            )
    files = {
        path.name: {"sha256": sha256(path.read_bytes()), "bytes": path.stat().st_size}
        for path in sorted(out.iterdir())
        if path.name != "index.json"
    }
    index = {
        "schema": "crn-craft.merge-reference-index/0.1",
        "manifest_sha256": sha256(MANIFEST.read_bytes()),
        "tools": {"crn_craft_cli": binary_identity(binary)},
        "files": files,
    }
    (out / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
