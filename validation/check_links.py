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

It checks `docs/SUMMARY.md`, the sidebar of the rendered site: every entry
must exist, no page may be listed twice, and every Markdown file under
`docs/guide`, `docs/reference`, `docs/developer` and `docs/evidence` must be
listed.

With `--site PATH` it also checks the built site (`mdbook build`): every
relative `href` and `src` in the HTML pages must resolve to a file of the
site, and a fragment must name an `id` of the target page.

    python validation/check_links.py [--root PATH] [--include-campaign]
                                     [--site PATH]
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


SUMMARY_ENTRY = re.compile(r"^\s*(?:-\s+)?\[[^\]]*\]\(([^)\s]*)\)\s*$")
SUMMARY_DIRS = ("guide", "reference", "developer", "evidence")


def check_summary(root):
    """Check docs/SUMMARY.md against the tracked pages. Returns (broken, entries)."""
    summary = root / "docs" / "SUMMARY.md"
    if not summary.exists():
        return ["docs/SUMMARY.md: missing"], 0
    entries, broken = [], []
    for number, line in enumerate(strip_fences(summary.read_text("utf-8")), 1):
        match = SUMMARY_ENTRY.match(line)
        if match:
            entries.append((number, unquote(match.group(1))))
    seen = {}
    for number, name in entries:
        where = f"docs/SUMMARY.md:{number}: {name}"
        if not name:
            broken.append(f"{where}: draft entry without a file")
        elif not (root / "docs" / name).is_file():
            broken.append(f"{where}: no such file")
        if name in seen:
            broken.append(f"{where}: listed twice (first at line {seen[name]})")
        seen.setdefault(name, number)
    for path in tracked_markdown(root):
        relative = (
            path.relative_to(root / "docs").as_posix()
            if path.is_relative_to(root / "docs")
            else None
        )
        if relative and relative.startswith(tuple(d + "/" for d in SUMMARY_DIRS)):
            if relative not in seen:
                broken.append(f"docs/{relative}: not listed in docs/SUMMARY.md")
    return broken, len(entries)


HTML_REF = re.compile(r"""\b(?:href|src)=["']([^"'#][^"']*|#[^"']+)["']""")
HTML_ID = re.compile(r"""\bid=["']([^"']+)["']""")


def check_site(site):
    """Check the relative references of the built HTML pages."""
    site = site.resolve()
    broken, count, ids = [], 0, {}
    for page in sorted(site.rglob("*.html")):
        text = page.read_text(encoding="utf-8")
        for target in HTML_REF.findall(text):
            if re.match(r"^([a-z][a-z0-9+.-]*:|//|/)", target, re.I):
                continue
            count += 1
            name, _, fragment = target.partition("#")
            resolved = page if not name else (page.parent / unquote(name)).resolve()
            if resolved.is_dir():
                resolved = resolved / "index.html"
            where = f"{page.relative_to(site)}: {target}"
            if not resolved.exists():
                broken.append(f"{where}: no such file")
            elif fragment and resolved.suffix == ".html":
                if resolved not in ids:
                    ids[resolved] = set(
                        HTML_ID.findall(resolved.read_text(encoding="utf-8"))
                    )
                if unquote(fragment) not in ids[resolved]:
                    broken.append(f"{where}: no such id")
    return broken, count


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--include-campaign", action="store_true")
    parser.add_argument(
        "--site",
        type=Path,
        help="only check the built site in this directory, not the sources",
    )
    args = parser.parse_args()
    if args.site:
        broken, count = check_site(args.site)
        for line in broken:
            print(line, file=sys.stderr)
        print(f"{count} site references checked, {len(broken)} broken")
        sys.exit(1 if broken else 0)
    root = args.root.resolve()
    broken, count = check(root, args.include_campaign)
    stub_broken, stubs = check_stubs(root)
    broken += stub_broken
    summary_broken, entries = check_summary(root)
    broken += summary_broken
    for line in broken:
        print(line, file=sys.stderr)
    print(
        f"{count} relative links checked, {stubs} stubs checked, "
        f"{entries} SUMMARY.md entries checked, {len(broken)} broken"
    )
    sys.exit(1 if broken else 0)


if __name__ == "__main__":
    main()
