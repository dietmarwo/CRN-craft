#!/usr/bin/env python3
"""Run every `bash` and `python` block of docs/tutorial.md, in order.

Blocks run in one scratch directory prepared like the tutorial's (copies of
models/examples plus a few inputs), with the native CLI first on PATH. Python
blocks share one namespace. Blocks fenced as `sh`, `json` or `text` are
illustrations and are not run. Requires the installed wheel (for
`crn_craft` and `crn-craft-import`) and numpy.
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

FENCE = re.compile(r"^```(\w*)\s*$")


def blocks(text):
    language, lines, start = None, [], 0
    for number, line in enumerate(text.splitlines(), 1):
        match = FENCE.match(line)
        if match and language is None:
            language, lines, start = match.group(1), [], number
        elif line.strip() == "```" and language is not None:
            yield language, start, "\n".join(lines) + "\n"
            language = None
        elif language is not None:
            lines.append(line)


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
                else:
                    continue
                rows.append(
                    {"line": line, "language": language, "output": output.strip()}
                )
        finally:
            os.chdir(previous)
    if len(rows) < 15:
        raise RuntimeError(f"only {len(rows)} runnable blocks found")
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
