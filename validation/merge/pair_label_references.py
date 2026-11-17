"""References of the 0.8 pair-label law check.

Declared in `validation/merge/RESULTS.md` ("0.8 pair labels: law check
declaration"). For every row of `validation/merge/pair-labels/models.json`
this writes, into `validation/merge/pair-labels/reference/`, the files of
`builder_references.py` computed from the row's attribute-free counts
network (`reference_path`) instead of the row's model:

- `<name>.moments.json`: sample moments under `static_direct`
  (`crn-craft.run/1.0`) with the reference seeds;
- `<name>.fsp-moments.json`: the moments of the finite state projection,
  kept only if the leakage and the mass defect are below 1e-9.

Standard library only; run under the CPU lock with the CLI under test:

    flock /home/claude/craft5/cpu.lock python validation/merge/pair_label_references.py \\
        --binary target/release/crn-craft
"""

import argparse
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import builder_references  # noqa: E402
from merge_common import HERE, ROOT, load_manifest, sha256  # noqa: E402

builder_references.MANIFEST = HERE / "pair-labels" / "models.json"
builder_references.OUTPUT = HERE / "pair-labels" / "reference"


def counts_rows(path):
    """The manifest rows with the counts network as their model."""
    rows = []
    for row in load_manifest(None, path):
        data = (ROOT / row["reference_path"]).read_bytes()
        if sha256(data) != row["reference_sha256"]:
            raise SystemExit(f"{row['reference_path']}: bytes differ from models.json")
        rows.append(
            {**row, "path": row["reference_path"], "sha256": row["reference_sha256"]}
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    rows = counts_rows(builder_references.MANIFEST)
    builder_references.load_manifest = lambda names, path: rows
    sys.argv = [
        sys.argv[0],
        "--binary",
        str(args.binary),
        "--threads",
        str(args.threads),
    ]
    with tempfile.TemporaryDirectory():
        builder_references.main()


if __name__ == "__main__":
    main()
