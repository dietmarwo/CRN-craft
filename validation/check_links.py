"""Check that every relative link of the documentation resolves.

Scans the Markdown files under `docs/` (without `docs/campaign/`, the
dated task logs), the READMEs of the repository and the tutorials. A link
`[text](target)` or `[text](target#fragment)` passes when the target file or
directory exists and, for a Markdown target, the fragment names a heading
(GitHub's slug rules). Links inside fenced code blocks and absolute URLs are
ignored. Exits 1 and lists every broken link otherwise.

It also checks the stubs that old paths keep after a move (a heading and one
sentence "This page moved to [title](path)."): every stub must point at an
existing file that is not itself a stub, and no tracked Markdown file except
the frozen `validation/compat/manual-0.6.0.md` may link to a stub. Links must
go to the real page.

    python validation/check_links.py [--root PATH] [--include-campaign]
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

LINK = re.compile(r"(?<!\!)\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
FENCE = re.compile(r"^(`{3,}|~{3,})")
STUB = re.compile(r"\A# [^\n]+\n\nThis page moved to \[[^\]]*\]\(([^)\s]+)\)\.\n?\Z")
FROZEN = "validation/compat/manual-0.6.0.md"


def tracked_markdown(root):
    out = subprocess.run(
        ["git", "-C", str(root), "ls-files", "*.md"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    return [root / name for name in out]


def wanted(root, path, include_campaign):
    relative = path.relative_to(root).as_posix()
    if relative.startswith("docs/campaign/") and not include_campaign:
        return False
    if relative.startswith(("docs/", "tutorials/", "examples/")):
        return True
    return "/" not in relative or relative.endswith("README.md")


def strip_fences(text):
    kept, fence = [], None
    for line in text.splitlines():
        match = FENCE.match(line)
        if fence is None and match:
            fence = match.group(1)
            kept.append("")
        elif fence is not None:
            if line.startswith(fence) and not line.strip(fence[0]):
                fence = None
            kept.append("")
        else:
            kept.append(line)
    return kept


def slug(heading):
    heading = re.sub(r"`", "", heading)
    heading = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", heading)
    heading = heading.strip().lower()
    heading = re.sub(r"[^\w\- ]", "", heading)
    return heading.replace(" ", "-")


_ANCHORS = {}


def anchors(path):
    cache = _ANCHORS
    if path not in cache:
        found, seen = set(), {}
        for line in strip_fences(path.read_text(encoding="utf-8")):
            match = re.match(r"^#{1,6}\s+(.*?)\s*#*\s*$", line)
            if match:
                base = slug(match.group(1))
                count = seen.get(base, 0)
                seen[base] = count + 1
                found.add(base if count == 0 else f"{base}-{count}")
        cache[path] = found
    return cache[path]


def check(root, include_campaign):
    broken, count = [], 0
    for path in tracked_markdown(root):
        if not wanted(root, path, include_campaign):
            continue
        lines = strip_fences(path.read_text(encoding="utf-8"))
        for number, line in enumerate(lines, 1):
            targets = LINK.findall(line) + IMAGE.findall(line)
            for target in targets:
                if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I):
                    continue
                count += 1
                name, _, fragment = target.partition("#")
                resolved = path if not name else (path.parent / unquote(name)).resolve()
                where = f"{path.relative_to(root)}:{number}: {target}"
                if not resolved.exists():
                    broken.append(f"{where}: no such file")
                elif fragment and resolved.suffix == ".md":
                    if unquote(fragment).lower() not in anchors(resolved):
                        broken.append(f"{where}: no such heading")
    return broken, count


def stub_target(path):
    """Return the link target of a stub page, or None if the page is no stub."""
    match = STUB.match(path.read_text(encoding="utf-8"))
    return match.group(1) if match else None


def check_stubs(root):
    """Check every stub's target and that no page links to a stub."""
    paths = tracked_markdown(root)
    stubs, broken = {}, []
    for path in paths:
        target = stub_target(path)
        if target is not None:
            stubs[path.resolve()] = target
    for path, target in stubs.items():
        name = target.partition("#")[0]
        resolved = (path.parent / unquote(name)).resolve()
        where = path.relative_to(root)
        if not resolved.exists():
            broken.append(f"{where}: stub target {target} does not exist")
        elif resolved in stubs:
            broken.append(f"{where}: stub target {target} is itself a stub")
    for path in paths:
        relative = path.relative_to(root).as_posix()
        if relative == FROZEN or path.resolve() in stubs:
            continue
        lines = strip_fences(path.read_text(encoding="utf-8"))
        for number, line in enumerate(lines, 1):
            for target in LINK.findall(line) + IMAGE.findall(line):
                name = target.partition("#")[0]
                if not name or re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I):
                    continue
                if (path.parent / unquote(name)).resolve() in stubs:
                    broken.append(
                        f"{relative}:{number}: {target}: links to a stub, "
                        "link the real page"
                    )
    return broken, len(stubs)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--include-campaign", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    broken, count = check(root, args.include_campaign)
    stub_broken, stubs = check_stubs(root)
    broken += stub_broken
    for line in broken:
        print(line, file=sys.stderr)
    print(
        f"{count} relative links checked, {stubs} stubs checked, {len(broken)} broken"
    )
    sys.exit(1 if broken else 0)


if __name__ == "__main__":
    main()
