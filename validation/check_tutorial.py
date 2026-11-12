#!/usr/bin/env python3
"""Run every `bash` and `python` block of docs/tutorial.md, in order.

Blocks run in one scratch directory prepared like the tutorial's (copies of
models/examples plus a few inputs), with the native CLI first on PATH. Python
blocks share one namespace. Every fenced block is checked:

- `bash` and `python` blocks run, and a failure stops the check;
- a `text` block must equal the standard output of the runnable block before
  it, character for character after stripping surrounding whitespace, and
  must follow such a block;
- a `json` block must be strict JSON without duplicate keys; a document with a
  `crn-craft.model/` or `crn-craft.run/` schema must also pass
  `crn-craft check` (a request against the model file it names);
- an `sh` block (installation commands) must pass `bash -n` and is not run.

Any other fence language is an error. Requires the installed wheel (for
`crn_craft` and `crn-craft-import`), numpy, pandas and xarray.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples/tutorial"))
from workdir import prepare_workdir  # noqa: E402

FENCE = re.compile(r"^(`{3,})(\w*)\s*$")


def blocks(text):
    """Yield (language, line, code) of every fenced block.

    A fence is three or more backticks and closes with a fence of at least the
    same length, so a block can show markdown that holds fences itself.
    """
    language, lines, start, ticks = None, [], 0, 0
    for number, line in enumerate(text.splitlines(), 1):
        match = FENCE.match(line)
        if match and language is None:
            ticks, language = len(match.group(1)), match.group(2)
            lines, start = [], number
        elif (
            language is not None
            and line.strip() == "`" * len(line.strip())
            and len(line.strip()) >= ticks
        ):
            yield language, start, "\n".join(lines) + "\n"
            language = None
        elif language is not None:
            lines.append(line)


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


def check_json_block(code, line, folder, env):
    """A json block is strict JSON; models and requests also pass `check`."""
    try:
        document = strict_json(code)
    except ValueError as error:
        raise RuntimeError(
            f"json block at line {line} is not strict JSON: {error}"
        ) from None
    schema = document.get("schema", "") if isinstance(document, dict) else ""
    arguments = None
    if schema.startswith("crn-craft.model/"):
        target = folder / f"block-{line}.model.json"
        target.write_text(code)
        arguments = [target.name]
    elif schema.startswith("crn-craft.run/"):
        target = folder / f"block-{line}.run.json"
        target.write_text(code)
        arguments = [document["model"]["path"], target.name]
    if arguments is not None:
        proc = subprocess.run(
            ["crn-craft", "check", *arguments],
            cwd=folder,
            env=env,
            capture_output=True,
            text=True,
            check=False,
            timeout=600,
        )
        report = json.loads(proc.stdout or "{}")
        if proc.returncode or not report.get("valid"):
            raise RuntimeError(
                f"json block at line {line} fails crn-craft check:\n{proc.stdout}"
            )
        target.unlink()
    return "json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if shutil.which("crn-craft-import") is None:
        raise RuntimeError("install the crn-craft wheel first (crn-craft-import)")
    text = (ROOT / "docs/tutorial.md").read_text()
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-tutorial-check-") as temp:
        folder = Path(temp)
        prepare_workdir(folder)
        bindir = folder / "bin"
        bindir.mkdir()
        shutil.copy2(args.binary.resolve(), bindir / "crn-craft")
        env = dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ['PATH']}")
        namespace = {"__name__": "__tutorial__"}
        previous = Path.cwd()
        os.chdir(folder)
        # The output of the last runnable block, for the text block after it.
        last_output = None
        try:
            for language, line, code in blocks(text):
                if language == "bash":
                    proc = subprocess.run(
                        ["bash", "-e", "-o", "pipefail", "-c", code],
                        cwd=folder,
                        env=env,
                        capture_output=True,
                        text=True,
                        check=False,
                        timeout=600,
                    )
                    if proc.returncode:
                        raise RuntimeError(
                            f"bash block at line {line} failed:\n"
                            f"{proc.stdout}{proc.stderr}"
                        )
                    output = proc.stdout
                elif language == "python":
                    buffer = io.StringIO()
                    with contextlib.redirect_stdout(buffer):
                        exec(compile(code, f"tutorial.md:{line}", "exec"), namespace)
                    output = buffer.getvalue()
                elif language == "text":
                    if last_output is None:
                        raise RuntimeError(
                            f"text block at line {line} follows no runnable block"
                        )
                    if code.strip() != last_output.strip():
                        raise RuntimeError(
                            f"text block at line {line} differs from the output of "
                            f"the block before it:\n--- expected\n{last_output}"
                            f"--- found\n{code}"
                        )
                    last_output = None
                    rows.append(
                        {"line": line, "language": "text", "output": code.strip()}
                    )
                    continue
                elif language == "json":
                    check_json_block(code, line, folder, env)
                    rows.append({"line": line, "language": "json", "output": ""})
                    continue
                elif language == "sh":
                    proc = subprocess.run(
                        ["bash", "-n", "-c", code],
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    if proc.returncode:
                        raise RuntimeError(
                            f"sh block at line {line} has a syntax error:\n"
                            f"{proc.stderr}"
                        )
                    rows.append({"line": line, "language": "sh", "output": ""})
                    continue
                else:
                    raise RuntimeError(
                        f"block at line {line} has the unchecked language "
                        f"{language!r}; use bash, python, text, json or sh"
                    )
                last_output = output
                rows.append(
                    {"line": line, "language": language, "output": output.strip()}
                )
        finally:
            os.chdir(previous)
    runnable = sum(row["language"] in ("bash", "python") for row in rows)
    if runnable < 15:
        raise RuntimeError(f"only {runnable} runnable blocks found")
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.tutorial-check/0.1",
                "passed": True,
                "blocks": rows,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "blocks": len(rows)}))


if __name__ == "__main__":
    main()
