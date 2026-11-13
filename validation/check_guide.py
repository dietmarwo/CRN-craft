#!/usr/bin/env python3
"""Run the code blocks of the user guide and the READMEs, and compare outputs.

Checks every `docs/guide/*.md`, `README.md`, `python/README.md` and
`docs/index.md` (or the files named on the command line), one file per
process:

    python validation/check_guide.py --binary target/release/crn-craft \\
        --require-cargo --output guide.json
    python -O validation/check_guide.py ... --output guide-optimized.json
    diff guide.json guide-optimized.json

The rules, per file:

1. One file is one session. The file runs in a fresh Python process and a
   fresh, empty scratch directory. Nothing passes from one file to the next.
2. The first `python` block of a file is its setup block.
3. Python names are local to a level-2 section. At every level-2 heading
   after the setup block, the namespace is replaced by a fresh one in which
   the setup block has run again (its output is discarded). So "setup block
   plus one section" is what every section is checked as. The scratch
   directory is kept for the whole file: a file that one section writes is
   there for a later one (chapters 10, 11 and 14 rely on this).
4. Every file a page reads is written by one of its blocks, or is a
   repository file named in a comment `<!-- fixtures: PATH ... -->`. Those
   are copied into the scratch directory under the same relative paths.
5. Fence languages: `python` and `bash` run. `text` is the expected output
   of the runnable block before it and must follow one. `json` is strict
   JSON (no duplicate keys, no NaN or Infinity), and a document with a
   `crn-craft.model/` or `crn-craft.run/` schema must pass `crn-craft
   check`. `sh` passes `bash -n` and never runs. `r` is allowed only in a
   page marked `<!-- check: none -->`. Any other language, and a fence
   without one, is an error.
6. Output is compared exactly after stripping surrounding whitespace and
   trailing spaces. On lines that hold a timing (`<number> s`) the number is
   replaced on both sides and runs of whitespace are collapsed.
7. `<!-- check: standalone -->`: every `python` block runs alone in a fresh
   interpreter from the repository root. `<!-- check: none -->`: nothing
   runs; the blocks are listed as not run. `<!-- check: long -->` (with
   standalone): without `--long` the blocks run in smoke mode, with
   `CRN_CRAFT_SMOKE=1` (few replicates, three evaluations), and their
   output is not compared; with `--long` they run in full and their output
   is compared or written. A comment `<!-- file: PATH -->` before a block
   binds it to a repository file: the block must equal the file byte for
   byte, and a standalone page runs the file itself (`python PATH` from the
   repository root), so that a script can import its neighbors.
8. Blocks that call `compile_loop` need cargo, and are the last blocks of
   their page. Without cargo the session stops at the first such block and
   lists the rest of the file as skipped; with `--require-cargo` (CI) that
   is a failure. Compiled loops build against this checkout
   (`CRN_CRAFT_SOURCE` is set to the repository root unless already set).
9. `validation/guide-blocks.json` (`--counts`) maps every checked file to
   the least number of runnable blocks it has. A file with fewer, a checked
   file without an entry, and (when the default set of files is checked) an
   entry without a file all fail.
10. The example of `README.md` and its output equal those of
    `docs/index.md`, byte for byte.

`--max-replicates N` caps every `crn.ensemble` and `crn.rules.ensemble`
call at N replicates (and `keep_paths` at that count) without changing the
pages: a quick smoke run whose outputs are not compared. `bash` blocks run
unchanged.

`--binary PATH` puts that `crn-craft` first on PATH for `bash` blocks and
`json` checks; without it the console script of the installed wheel runs.
The checker also runs under `python -O` (the sessions then run with `-O`
too, which skips the `assert` statements of the blocks); the report has no
timings, so the two reports are equal.

With `--write`, every `text` block that follows a runnable block is replaced
by the real output of that block. Nothing else in the file changes. Expected
outputs are never written by hand.

Requires the installed wheel (`crn_craft`, `crn-craft-import`), numpy,
pandas, xarray and matplotlib; fcmaes, scipy and loguru for the fitting
page.
"""

import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "validation" / "guide-blocks.json"
FENCE = re.compile(r"^(`{3,})(\w*)\s*$")
SECONDS = re.compile(r"\b\d+(?:\.\d+)? s\b")
COMMENT = re.compile(r"^<!--\s*(check|fixtures|file)\s*:\s*(.*?)\s*-->\s*$")
RUNNABLE = ("python", "bash")
STATIC = ("text", "json", "sh")
MODES = {"session", "standalone", "none"}
COUNTED = ("ran", "skipped: needs cargo")
# The README's example must equal the landing page's (section 7.2 of the
# table of contents of the docs-0.7.1 campaign).
SAME_EXAMPLE = ("README.md", "docs/index.md")

# A launcher for capped runs: patch the ensemble functions, then run the
# script as __main__ with its own folder first on sys.path. `CAP` is filled in.
CAPPED = """\
import os
import runpy
import sys

import crn_craft as _crn
import crn_craft.rules as _rules


def _cap(ensemble):
    def capped(model, *, replicates=1, keep_paths=0, **kwargs):
        replicates = min(int(replicates), CAP)
        if not isinstance(keep_paths, (bool, str)):
            keep_paths = min(int(keep_paths), replicates)
        return ensemble(model, replicates=replicates, keep_paths=keep_paths, **kwargs)

    return capped


_rules.ensemble = _cap(_rules.ensemble)
_crn.ensemble = _cap(_crn.ensemble)
sys.argv = sys.argv[1:]
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
runpy.run_path(sys.argv[0], run_name="__main__")
"""


class CheckError(Exception):
    """A block failed, or the file breaks a rule."""


def parse(text):
    """Return (blocks, options) of a markdown file.

    A block is a dictionary with `language`, `line` (of the opening fence),
    `end` (line of the closing fence), `ticks`, `code`, `section` (the number
    of level-2 headings before it), `heading` (the text of the last one) and
    `file` (the path of a `file` comment right before it, or None). A fence
    is three or more backticks and closes with a line of at least as many,
    so a block can show markdown that holds fences itself. Options are the
    `check` words and the `fixtures` paths of the file's comments.
    """
    blocks, checks, fixtures = [], set(), []
    current, section, heading, bound = None, 0, "", None
    for number, line in enumerate(text.splitlines(), 1):
        if current is None:
            match = FENCE.match(line)
            if match:
                current = {
                    "language": match.group(2),
                    "line": number,
                    "ticks": len(match.group(1)),
                    "section": section,
                    "heading": heading,
                    "lines": [],
                    "file": bound,
                }
                bound = None
            elif line.startswith("## "):
                section += 1
                heading = line[3:].strip()
            else:
                comment = COMMENT.match(line)
                if comment and comment.group(1) == "file":
                    bound = comment.group(2)
                elif comment and comment.group(1) == "check":
                    checks.update(comment.group(2).replace(",", " ").split())
                elif comment:
                    fixtures.extend(comment.group(2).split())
        else:
            stripped = line.strip()
            if (
                stripped
                and stripped == "`" * len(stripped)
                and len(stripped) >= current["ticks"]
            ):
                current["end"] = number
                current["code"] = "\n".join(current.pop("lines")) + "\n"
                blocks.append(current)
                current = None
            else:
                current["lines"].append(line)
    if current is not None:
        raise CheckError(f"the fence at line {current['line']} is never closed")
    return blocks, {"checks": checks, "fixtures": fixtures}


def mode_of(options):
    unknown = options["checks"] - MODES - {"long"}
    if unknown:
        raise CheckError(f"unknown check option(s): {', '.join(sorted(unknown))}")
    modes = options["checks"] & MODES
    if len(modes) > 1:
        raise CheckError(f"more than one check mode: {', '.join(sorted(modes))}")
    mode = modes.pop() if modes else "session"
    if "long" in options["checks"] and mode != "standalone":
        raise CheckError("check: long needs check: standalone")
    return mode


def normalize(text):
    """The lines of an output, with timings and their spacing made equal."""
    lines = []
    for line in text.strip().splitlines():
        line = line.rstrip()
        if SECONDS.search(line):
            line = " ".join(SECONDS.sub("<t> s", line).split())
        lines.append(line)
    return lines


def strict_json(text):
    """Parse JSON, rejecting duplicate keys and non-finite constants."""

    def pairs(items):
        document = {}
        for key, value in items:
            if key in document:
                raise ValueError(f"duplicate key {key!r}")
            document[key] = value
        return document

    def constant(name):
        raise ValueError(f"non-standard constant {name}")

    return json.loads(text, object_pairs_hook=pairs, parse_constant=constant)


def check_json(block, folder, env):
    """A json block is strict JSON; models and requests also pass `check`."""
    line = block["line"]
    try:
        document = strict_json(block["code"])
    except ValueError as error:
        raise CheckError(
            f"json block at line {line} is not strict JSON: {error}"
        ) from None
    if folder is None:
        return
    schema = document.get("schema", "") if isinstance(document, dict) else ""
    if schema.startswith("crn-craft.model/"):
        target = folder / f"block-{line}.model.json"
        arguments = [target.name]
    elif schema.startswith("crn-craft.run/"):
        target = folder / f"block-{line}.run.json"
        arguments = [document["model"]["path"], target.name]
    else:
        return
    target.write_text(block["code"], encoding="utf-8")
    done = subprocess.run(
        ["crn-craft", "check", *arguments],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )
    target.unlink()
    try:
        valid = json.loads(done.stdout).get("valid")
    except ValueError:
        valid = False
    if done.returncode or not valid:
        raise CheckError(
            f"json block at line {line} fails crn-craft check:\n{done.stdout}{done.stderr}"
        )


def check_sh(block):
    done = subprocess.run(
        ["bash", "-n", "-c", block["code"]], capture_output=True, text=True, check=False
    )
    if done.returncode:
        raise CheckError(
            f"sh block at line {block['line']} has a syntax error:\n{done.stderr}"
        )


def copy_fixtures(paths, root, folder):
    for relative in paths:
        source = root / relative
        target = folder / relative
        if Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise CheckError(
                f"fixture {relative!r} must be a path inside the repository"
            )
        if source.is_dir():
            shutil.copytree(source, target, dirs_exist_ok=True)
        elif source.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        else:
            raise CheckError(f"fixture {relative!r} is not in the repository")


def environment(base, binary, root):
    """The environment of the blocks: the checked tools first on PATH."""
    path = [str(Path(sys.executable).parent), os.environ.get("PATH", "")]
    if binary is not None:
        # Beside the scratch directory, so the page starts in an empty one.
        bindir = base / "bin"
        bindir.mkdir()
        shutil.copy2(binary, bindir / "crn-craft")
        path.insert(0, str(bindir))
    env = dict(os.environ, PATH=os.pathsep.join(path), PYTHONDONTWRITEBYTECODE="1")
    # Compiled loops build against this checkout, which may be newer than the
    # release on crates.io.
    env.setdefault("CRN_CRAFT_SOURCE", str(root))
    return env


def cap_ensembles(limit):
    """In this process, cap the replicates of every ensemble call."""
    import crn_craft as crn
    import crn_craft.rules as rules

    def cap(ensemble):
        def capped(model, *, replicates=1, keep_paths=0, **kwargs):
            replicates = min(int(replicates), limit)
            if not isinstance(keep_paths, bool | str):
                keep_paths = min(int(keep_paths), replicates)
            return ensemble(
                model, replicates=replicates, keep_paths=keep_paths, **kwargs
            )

        return capped

    rules.ensemble = cap(rules.ensemble)
    crn.ensemble = cap(crn.ensemble)


def run_python(code, name, namespace):
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        exec(compile(code, name, "exec"), namespace)  # noqa: S102
    return buffer.getvalue()


def run_bash(code, folder, env):
    done = subprocess.run(
        ["bash", "-e", "-o", "pipefail", "-c", code],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=1800,
    )
    if done.returncode:
        raise CheckError(f"exit status {done.returncode}:\n{done.stdout}{done.stderr}")
    return done.stdout


def run_standalone(block, root, base, env, smoke, cap):
    """Run one block of a standalone page in a fresh interpreter."""
    flags = ["-O"] * sys.flags.optimize
    if cap is None:
        source = [block["file"]] if block["file"] else ["-c", block["code"]]
        command = [sys.executable, *flags, *source]
    else:
        script = block["file"]
        if script is None:
            script = base / f"block-{block['line']}.py"
            script.write_text(block["code"], encoding="utf-8")
        launcher = CAPPED.replace("CAP", str(cap))
        command = [sys.executable, *flags, "-c", launcher, str(script)]
    done = subprocess.run(
        command,
        cwd=root,
        env=env | {"CRN_CRAFT_SMOKE": "1"} if smoke else env,
        capture_output=True,
        text=True,
        check=False,
        timeout=7200,
    )
    if done.returncode:
        raise CheckError(done.stderr)
    return done.stdout


def compare(block, output):
    """Compare a `text` block with the output of the block before it."""
    shown, actual = normalize(block["code"]), normalize(output)
    if shown != actual:
        raise CheckError(
            f"text block at line {block['line']} differs from the output of the "
            f"block before it:\n--- printed\n{output}--- shown\n{block['code']}"
        )


def not_run(blocks):
    """The rows of a `check: none` page; its static blocks are still checked."""
    rows = []
    for block in blocks:
        language = block["language"]
        if language not in (*RUNNABLE, *STATIC, "r"):
            raise CheckError(
                f"block at line {block['line']} has the unchecked language "
                f"{language!r}; use python, bash, text, json, sh or r"
            )
        row = {"line": block["line"], "language": language, "status": "not run"}
        if language == "sh":
            check_sh(block)
            row["status"] = "syntax checked"
        elif language == "json":
            check_json(block, None, None)
            row["status"] = "valid"
        rows.append(row | {"section": block["heading"]})
    return rows


def session(path, root, args):
    """Run one file in this process. Return (rows, outputs by text-block line)."""
    name = path.name
    blocks, options = parse(path.read_text(encoding="utf-8"))
    mode = mode_of(options)
    rows, outputs = [], {}
    if mode == "none":
        return not_run(blocks), outputs
    base = Path(tempfile.mkdtemp(prefix="crn-guide-"))
    folder = base / "work"
    folder.mkdir()
    copy_fixtures(options["fixtures"], root, folder)
    env = environment(base, args.binary, root)
    os.environ.update({k: env[k] for k in ("PATH", "CRN_CRAFT_SOURCE")})
    if args.max_replicates is not None and mode == "session":
        cap_ensembles(args.max_replicates)
    cargo = shutil.which("cargo") is not None
    long = "long" in options["checks"]
    smoke = long and not args.long
    capped = args.max_replicates is not None
    setup = next((b for b in blocks if b["language"] == "python"), None)
    namespace = {"__name__": "__main__"}
    section = setup["section"] if setup else 0
    # What the text block after a runnable block is compared with: the output
    # of that block, a reason why it is not compared, or None (no such block).
    last = None
    os.chdir(folder)
    try:
        for index, block in enumerate(blocks):
            language, line = block["language"], block["line"]
            row = {"line": line, "language": language, "section": block["heading"]}
            if language in RUNNABLE:
                if "compile_loop" in block["code"] and not cargo:
                    if args.require_cargo:
                        raise CheckError(
                            f"block at line {line} needs cargo (https://rustup.rs)"
                        )
                    for rest in blocks[index:]:
                        rows.append(
                            {
                                "line": rest["line"],
                                "language": rest["language"],
                                "section": rest["heading"],
                                "status": "skipped: needs cargo",
                            }
                        )
                    break
                if block["file"] is not None:
                    bound = root / block["file"]
                    if (
                        not bound.is_file()
                        or bound.read_bytes().decode("utf-8") != block["code"]
                    ):
                        raise CheckError(
                            f"block at line {line} differs from its file {block['file']}"
                        )
                try:
                    if mode == "standalone":
                        if language != "python":
                            raise CheckError(
                                "a standalone page runs python blocks only"
                            )
                        output = run_standalone(
                            block, root, base, env, smoke, args.max_replicates
                        )
                    elif language == "bash":
                        output = run_bash(block["code"], folder, env)
                    else:
                        if block is not setup and block["section"] != section:
                            # A new level-2 section: only the setup block's names.
                            namespace = {"__name__": "__main__"}
                            run_python(
                                setup["code"], f"{name}:{setup['line']}", namespace
                            )
                        output = run_python(block["code"], f"{name}:{line}", namespace)
                    if language == "python":
                        section = block["section"]
                except CheckError as error:
                    raise CheckError(
                        f"{language} block at line {line} failed: {error}"
                    ) from None
                except Exception as error:
                    import traceback

                    raise CheckError(
                        f"{language} block at line {line} failed:\n"
                        f"{traceback.format_exc()}"
                    ) from error
                row["status"] = "ran"
                if smoke:
                    row["smoke"] = True
                if smoke:
                    last = ("not compared: smoke",)
                elif capped:
                    last = ("not compared: capped",)
                else:
                    last = output
            elif language == "text":
                if last is None:
                    raise CheckError(
                        f"text block at line {line} follows no runnable block"
                    )
                if isinstance(last, tuple):
                    row["status"] = last[0]
                elif args.write:
                    outputs[line] = last
                    row["status"] = "written"
                else:
                    compare(block, last)
                    row["status"] = "matches"
                last = None
            elif language == "json":
                check_json(block, folder, env)
                row["status"] = "valid"
            elif language == "sh":
                check_sh(block)
                row["status"] = "syntax checked"
            else:
                raise CheckError(
                    f"block at line {line} has the unchecked language {language!r}; "
                    "use python, bash, text, json or sh"
                )
            rows.append(row)
    finally:
        os.chdir(root)
    if args.keep:
        rows.append({"line": 0, "language": "", "status": f"scratch kept: {folder}"})
    else:
        shutil.rmtree(base, ignore_errors=True)
    return rows, outputs


def rewrite(path, outputs):
    """Replace the text blocks named by line with the given outputs."""
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    blocks, _ = parse(text)
    for block in reversed(blocks):
        if block["line"] not in outputs:
            continue
        body = outputs[block["line"]].strip("\n").splitlines()
        longest = max(
            (len(m.group(0)) for row in body for m in re.finditer(r"`{3,}", row)),
            default=0,
        )
        fence = "`" * max(3, longest + 1)
        lines[block["line"] - 1 : block["end"]] = [f"{fence}text", *body, fence]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def child(args):
    """Run one file and print its report as JSON on the last line of stdout."""
    path = args.files[0].resolve()
    root = args.root.resolve()
    report = {"file": path.relative_to(root).as_posix(), "passed": True}
    try:
        rows, outputs = session(path, root, args)
        if args.write and outputs:
            rewrite(path, outputs)
        report["blocks"] = rows
    except CheckError as error:
        report |= {"passed": False, "error": str(error)}
    # Blocks may leave threads or redirected streams; report on the real stdout.
    sys.__stdout__.write(json.dumps(report) + "\n")
    sys.__stdout__.flush()
    return 0 if report["passed"] else 1


def example_of(path):
    """The first python block of a page and the text block after it."""
    blocks, _ = parse(path.read_text(encoding="utf-8"))
    for index, block in enumerate(blocks):
        if block["language"] == "python":
            after = blocks[index + 1] if index + 1 < len(blocks) else None
            output = after["code"] if after and after["language"] == "text" else None
            return block["code"], output
    return None, None


def same_example(root):
    """Check that the README's example and output equal the landing page's."""
    first, second = (root / name for name in SAME_EXAMPLE)
    entry = {"files": list(SAME_EXAMPLE), "passed": True}
    if not first.is_file() or not second.is_file():
        return entry | {"passed": False, "error": "a page is missing"}
    code_a, output_a = example_of(first)
    code_b, output_b = example_of(second)
    if code_a is None or output_a is None:
        return entry | {"passed": False, "error": "README.md has no example"}
    if code_a != code_b:
        return entry | {"passed": False, "error": "the example blocks differ"}
    if output_a != output_b:
        return entry | {"passed": False, "error": "the shown outputs differ"}
    return entry


def default_files(root):
    pages = sorted((root / "docs" / "guide").glob("*.md"))
    pages += [root / "README.md", root / "python" / "README.md"]
    pages.append(root / "docs" / "index.md")
    return sorted(pages, key=lambda p: p.relative_to(root).as_posix())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "files", type=Path, nargs="*", help="markdown files (default: all)"
    )
    parser.add_argument("--binary", type=Path, help="the native crn-craft tool")
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root")
    parser.add_argument("--require-cargo", action="store_true")
    parser.add_argument("--max-replicates", type=int, help="quick smoke run")
    parser.add_argument(
        "--long", action="store_true", help="run long pages in full, not in smoke mode"
    )
    parser.add_argument("--write", action="store_true", help="rewrite the text blocks")
    parser.add_argument(
        "--keep", action="store_true", help="keep the scratch directory"
    )
    parser.add_argument(
        "--counts", type=Path, default=COUNTS, help="JSON: file to least block count"
    )
    parser.add_argument("--output", type=Path, help="write the report here")
    parser.add_argument("--session", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.max_replicates is not None and args.max_replicates < 1:
        parser.error("--max-replicates must be positive")
    if args.write and args.max_replicates is not None:
        parser.error("--write needs uncapped runs")
    if args.binary is not None:
        args.binary = args.binary.resolve()
    args.root = args.root.resolve()
    if args.session:
        return child(args)
    every = not args.files
    files = default_files(args.root) if every else [p.resolve() for p in args.files]
    counts = json.loads(args.counts.read_text(encoding="utf-8"))
    reports, failed = [], False
    for path in files:
        relative = path.relative_to(args.root).as_posix()
        command = [sys.executable, *["-O"] * sys.flags.optimize, __file__, "--session"]
        command += ["--root", str(args.root), str(path)]
        for flag in ("write", "require_cargo", "long", "keep"):
            if getattr(args, flag):
                command.append("--" + flag.replace("_", "-"))
        if args.binary is not None:
            command += ["--binary", str(args.binary)]
        if args.max_replicates is not None:
            command += ["--max-replicates", str(args.max_replicates)]
        start = time.perf_counter()
        done = subprocess.run(
            command, capture_output=True, text=True, check=False, cwd=args.root
        )
        try:
            report = json.loads(done.stdout.strip().splitlines()[-1])
        except (IndexError, ValueError):
            report = {
                "file": relative,
                "passed": False,
                "error": f"the session ended without a report:\n{done.stdout}{done.stderr}",
            }
        rows = report.get("blocks", [])
        runnable = sum(
            row["language"] in RUNNABLE and row["status"] in COUNTED for row in rows
        )
        least = counts.get(relative)
        if report["passed"] and least is None:
            report |= {
                "passed": False,
                "error": f"{relative} has no entry in {args.counts.name}",
            }
        elif report["passed"] and runnable < least:
            report |= {
                "passed": False,
                "error": f"{runnable} runnable blocks; at least {least} are expected",
            }
        report["runnable"] = runnable
        reports.append(report)
        failed |= not report["passed"]
        compared = sum(row["status"] in ("matches", "written") for row in rows)
        skipped = sum(row["status"].startswith(("skipped", "not ")) for row in rows)
        print(
            f"{relative}: {'ok' if report['passed'] else 'FAILED'}; "
            f"{runnable} blocks ran, {compared} outputs "
            f"{'written' if args.write else 'compared'}, {skipped} not run or "
            f"not compared, {time.perf_counter() - start:.1f} s",
            file=sys.stderr,
        )
        for row in rows:
            if row["status"].startswith("scratch kept"):
                print(f"  {row['status']}", file=sys.stderr)
        if not report["passed"]:
            print(report["error"], file=sys.stderr)
            print(done.stderr[-4000:], file=sys.stderr)
    problems = []
    if every:
        checked = {r["file"] for r in reports}
        problems += [
            f"{name} is in {args.counts.name} but is not a checked page"
            for name in sorted(set(counts) - checked)
        ]
    example = same_example(args.root)
    if not example["passed"]:
        problems.append(f"README.md and docs/index.md: {example['error']}")
    for problem in problems:
        print(problem, file=sys.stderr)
    failed |= bool(problems)
    document = {
        "schema": "crn-craft.guide-check/1.0",
        "passed": not failed,
        "capped_replicates": args.max_replicates,
        "long": args.long,
        "files": reports,
        "same_example": example,
        "problems": problems,
        "runnable_blocks": sum(r["runnable"] for r in reports),
    }
    if args.output:
        args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "passed": not failed,
                "files": len(reports),
                "runnable_blocks": document["runnable_blocks"],
            }
        )
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
