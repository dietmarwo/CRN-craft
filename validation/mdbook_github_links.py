"""mdBook preprocessor: links to pages outside the book point to GitHub.

mdBook renders only the pages named in `docs/SUMMARY.md`. The pages stay
plain Markdown with relative links, so a link to a file that is not in the
book (campaign logs, decisions, stubs, CHANGELOG.md, source files, a
directory) would be dead on the site. This preprocessor rewrites such links
to `https://github.com/dietmarwo/CRN-craft/blob/main/<path>` (`tree` for a
directory). Links to book pages and to assets under `docs/` are untouched,
as are absolute URLs and links in fenced code blocks.

Usage by mdBook: `python3 validation/mdbook_github_links.py` reads
`[context, book]` as JSON on stdin and writes the book to stdout.
"""

import json
from pathlib import Path
import posixpath
import re
import sys

REPO = "https://github.com/dietmarwo/CRN-craft"
BRANCH = "main"
LINK = re.compile(r"(\]\()([^)\s]+)((?:\s+\"[^\"]*\")?\))")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})")


def chapters(items):
    for item in items:
        chapter = item.get("Chapter") if isinstance(item, dict) else None
        if chapter:
            yield chapter
            yield from chapters(chapter["sub_items"])


def rewrite(content, chapter_path, book_pages, docs_dir):
    out, fence = [], None
    for line in content.split("\n"):
        match = FENCE.match(line)
        if fence is None and match:
            fence = match.group(1)
        elif fence is not None and line.strip().startswith(fence):
            fence = None
        elif fence is None:
            line = LINK.sub(
                lambda m: m.group(1)
                + convert(m.group(2), chapter_path, book_pages, docs_dir)
                + m.group(3),
                line,
            )
        out.append(line)
    return "\n".join(out)


def convert(target, chapter_path, book_pages, docs_dir):
    if re.match(r"^[a-z][a-z0-9+.-]*:|^#|^/", target, re.I):
        return target
    name, hash_, fragment = target.partition("#")
    relative = posixpath.normpath(posixpath.join(posixpath.dirname(chapter_path), name))
    if relative in book_pages:
        return target
    path = (docs_dir / relative).resolve()
    if not path.exists():
        return target
    if path.is_file() and path.suffix != ".md" and path.is_relative_to(docs_dir):
        return target  # an asset: mdBook copies it
    repo_path = path.relative_to(docs_dir.parent).as_posix()
    kind = "tree" if path.is_dir() else "blob"
    return f"{REPO}/{kind}/{BRANCH}/{repo_path}{hash_}{fragment}"


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "supports":
        sys.exit(0)
    context, book = json.load(sys.stdin)
    docs_dir = (Path.cwd() / context["config"]["book"]["src"]).resolve()
    items = book["items"] if "items" in book else book["sections"]
    pages = {c["path"] for c in chapters(items) if c.get("path")}
    for chapter in chapters(items):
        if chapter.get("path"):
            chapter["content"] = rewrite(
                chapter["content"], chapter["path"], pages, docs_dir
            )
    json.dump(book, sys.stdout)


if __name__ == "__main__":
    main()
