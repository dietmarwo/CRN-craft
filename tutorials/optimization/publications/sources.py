"""Fetch the pinned research fixtures. Tutorial runs read local files only."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
from urllib.request import urlopen

MANIFEST = Path(__file__).parents[3] / "examples/publication_tutorials/sources.json"
DEFAULT = Path(os.environ.get("CRN_PUBLICATION_DATA", ".publication-data"))


def fixtures(directory, chapter, *, fetch=False):
    spec = json.loads(MANIFEST.read_text())[chapter]
    paths = {}
    for name, (url, expected) in spec.items():
        path = directory / chapter / name
        if fetch and not path.exists():
            with urlopen(url, timeout=60) as response:
                data = response.read(64 << 20 | 1)
            if len(data) > 64 << 20 or sha256(data).hexdigest() != expected:
                raise ValueError(f"invalid source: {name}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if not path.is_file():
            raise ValueError(
                f"missing {path}. Run publications/sources.py with --data-dir."
            )
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"checksum mismatch: {path}")
        paths[name] = path
    return paths


def data_option(parser):
    parser.add_argument("--data-dir", type=Path, default=DEFAULT)


if __name__ == "__main__":
    manifest = json.loads(MANIFEST.read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    data_option(parser)
    parser.add_argument("chapters", nargs="*", default=list(manifest))
    args = parser.parse_args()
    for identifier in args.chapters:
        for path in fixtures(args.data_dir, identifier, fetch=True).values():
            print(path)
