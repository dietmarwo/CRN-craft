"""Fetch and verify the pinned research fixtures.

Tutorial runs read local files only. A manifest entry is either a pair
``[url, sha256]``, which pins the bytes of one file, or an object with
``"kind": "tar-content"``, which pins the CONTENT of selected members of a
gzip tar archive. Generated repository archives are not byte stable, so
their container hash is not pinned. The content digest is the SHA-256 of the
lines ``<sha256 of member>  <member path>\\n``, sorted by path (code point
order), for every regular member whose path, without the archive's first
component, starts with one of ``prefixes`` or equals one of ``files``. From a
clone of the pinned commit the same digest is

    git ls-files <prefixes and files> | LC_ALL=C sort | xargs sha256sum \\
        | sha256sum

Archives are read as a stream: nothing is extracted, executed or unpickled.
Links, devices, absolute or parent-relative member names and duplicate
members are refused, and every read is bounded.
"""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import tarfile
from tempfile import NamedTemporaryFile
from urllib.parse import urlsplit
import urllib.request

MANIFEST = Path(__file__).parents[3] / "examples/publication_tutorials/sources.json"
DEFAULT = Path(os.environ.get("CRN_PUBLICATION_DATA", ".publication-data"))
HOSTS = frozenset({"raw.githubusercontent.com", "codeload.github.com", "github.com"})
MAX_FILE = 64 << 20
MAX_ARCHIVE = 160 << 20
MAX_MEMBER = 4 << 20
MAX_UNPACKED = 1 << 30
MAX_MEMBERS = 100_000


def _check_url(url):
    parts = urlsplit(url)
    if parts.scheme != "https" or parts.hostname not in HOSTS or parts.port:
        raise ValueError(f"refused source URL: {url}")


class _Redirects(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _download(url, limit):
    _check_url(url)
    opener = urllib.request.build_opener(_Redirects)
    with opener.open(url, timeout=60) as response:
        data = response.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"source exceeds {limit} bytes: {url}")
    return data


def _member_path(name):
    """The member path below the archive's top folder, or a refusal."""
    if "\0" in name or "\\" in name:
        raise ValueError(f"unsafe archive member: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe archive member: {name!r}")
    return "/".join(path.parts[1:])


def _selected(spec, path):
    return path in spec.get("files", ()) or any(
        path.startswith(prefix) for prefix in spec.get("prefixes", ())
    )


def read_tar(path, spec, handle=None):
    """Stream the selected members of a pinned archive and verify them.

    ``handle(member_path, data)`` sees each selected member once, in archive
    order. Results a caller builds from it are valid only if this function
    returns: it raises ``ValueError`` after the stream when the content
    digest or member count differs from the manifest, and at once for an
    unsafe or oversized member.
    """
    if path.stat().st_size > MAX_ARCHIVE:
        raise ValueError(f"archive too large: {path}")
    digests = {}
    unpacked = 0
    with tarfile.open(path, mode="r|gz") as source:
        for count, entry in enumerate(source):
            if count >= MAX_MEMBERS:
                raise ValueError(f"too many archive members: {path}")
            member = _member_path(entry.name)
            if not (entry.isfile() or entry.isdir()):
                raise ValueError(f"refused archive member type: {entry.name!r}")
            unpacked += entry.size
            if unpacked > MAX_UNPACKED:
                raise ValueError(f"archive content too large: {path}")
            if not entry.isfile() or not _selected(spec, member):
                continue
            if member in digests:
                raise ValueError(f"duplicate archive member: {member}")
            if entry.size > MAX_MEMBER:
                raise ValueError(f"archive member too large: {member}")
            data = source.extractfile(entry).read(MAX_MEMBER + 1)
            if len(data) != entry.size:
                raise ValueError(f"truncated archive member: {member}")
            digests[member] = sha256(data).hexdigest()
            if handle is not None:
                handle(member, data)
    lines = "".join(f"{digests[m]}  {m}\n" for m in sorted(digests))
    if len(digests) != spec["members"] or (
        sha256(lines.encode()).hexdigest() != spec["sha256"]
    ):
        raise ValueError(f"archive content differs from the manifest: {path}")


def _verify(path, spec):
    if isinstance(spec, dict):
        if spec.get("kind") != "tar-content":
            raise ValueError(f"unknown fixture kind: {path}")
        read_tar(path, spec)
        return
    _, expected = spec
    if path.stat().st_size > MAX_FILE:
        raise ValueError(f"fixture too large: {path}")
    if sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"checksum mismatch: {path}")


def fixtures(directory, chapter, *, fetch=False):
    """Paths of a chapter's fixtures.

    Single files are verified by their SHA-256 here. Archives are verified
    by :func:`read_tar` while the tutorial reads them; with ``fetch=True``
    they are also verified before they are stored.
    """
    spec = json.loads(MANIFEST.read_text())[chapter]
    paths = {}
    for name, entry in spec.items():
        path = directory / chapter / name
        if fetch and not path.exists():
            archive = isinstance(entry, dict)
            url = entry["url"] if archive else entry[0]
            data = _download(url, MAX_ARCHIVE if archive else MAX_FILE)
            path.parent.mkdir(parents=True, exist_ok=True)
            with NamedTemporaryFile(dir=path.parent, delete=False) as handle:
                handle.write(data)
            staged = Path(handle.name)
            try:
                _verify(staged, entry)
                staged.replace(path)
            finally:
                staged.unlink(missing_ok=True)
        if not path.is_file():
            raise ValueError(
                f"missing {path}. Run publications/sources.py with --data-dir."
            )
        if not isinstance(entry, dict):
            _verify(path, entry)
        paths[name] = path
    return paths


def data_option(parser):
    parser.add_argument("--data-dir", type=Path, default=DEFAULT)


if __name__ == "__main__":
    manifest = json.loads(MANIFEST.read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    data_option(parser)
    parser.add_argument(
        "--verify", action="store_true", help="Verify present fixtures only."
    )
    parser.add_argument("chapters", nargs="*", default=list(manifest))
    args = parser.parse_args()
    for identifier in args.chapters:
        found = fixtures(args.data_dir, identifier, fetch=not args.verify)
        for name, path in found.items():
            entry = manifest[identifier][name]
            if isinstance(entry, dict):
                read_tar(path, entry)
            print(path)
