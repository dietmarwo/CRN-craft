"""Run every Python block of the mtDNA user manual on its own.

Each ```python block runs in a fresh interpreter, in manual order, from the
repository root, with a timeout. A block passes when it exits with status 0
and, if the manual shows its output (the first ```text block after it and
before the next Python block), when it prints that output. Timings
(`<number> s`) and runs of whitespace are normalized before comparing.

A fence may use more than three backticks; it ends at a line of at least as
many. Blocks inside a longer fence, such as an example of a markdown report,
are text and never run.

Blocks that read `sys.argv` get a fresh directory and, with `--binary`, the
native CLI as arguments. Without `--binary` their output is not compared.

`--max-replicates N` caps every `crn.ensemble` and `crn.rules.ensemble` call at N replicates (and
`keep_paths` at that count) without changing the manual text. It is a quick
smoke run: the numbers then differ from the manual, so output is not
compared.

Blocks that call `rules.compile_loop` need cargo. Without it they are
skipped with a message, or fail with `--require-cargo` (as in CI). They
build against this checkout, which may be newer than the release on
crates.io: the checker sets `CRN_CRAFT_SOURCE` to the repository root
unless it is already set.
"""

import argparse
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
MANUAL = ROOT / "docs" / "mtdna-user-manual.md"
FENCE = re.compile(r"^(`{3,})(\w*)\s*$")
SECONDS = re.compile(r"\b\d+(?:\.\d+)? s\b")

# Prepended in --max-replicates mode; `CAP` is filled in.
CAPPED = """\
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


# Both spellings: `crn.ensemble` and `crn.rules.ensemble`.
_rules.ensemble = _cap(_rules.ensemble)
_crn.ensemble = _cap(_crn.ensemble)
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
"""


def blocks(text):
    """Yield (line, heading, code, shown output or None) per Python block."""
    heading, found, current = "", [], None
    lines = text.splitlines(keepends=True)
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("#"):
            heading = line.strip("# \n")
        match = FENCE.match(line)
        if not match:
            index += 1
            continue
        ticks, language = match.group(1), match.group(2)
        start, body = index + 1, []
        index += 1
        closing = re.compile(r"^" + ticks + r"`*\s*$")
        while index < len(lines) and not closing.match(lines[index]):
            body.append(lines[index])
            index += 1
        index += 1
        if language == "python":
            current = {"line": start, "heading": heading, "code": "".join(body)}
            current["shown"] = None
            found.append(current)
        elif language == "text" and current is not None and current["shown"] is None:
            current["shown"] = "".join(body)
    return found


def normalize(text):
    out = [" ".join(SECONDS.sub("<t> s", line).split()) for line in text.splitlines()]
    while out and not out[-1]:
        out.pop()
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manual", type=Path, default=MANUAL)
    parser.add_argument("--binary", type=Path, help="native crn-craft CLI")
    parser.add_argument("--timeout", type=float, default=600.0, help="per block, s")
    parser.add_argument("--max-replicates", type=int, help="quick smoke run")
    parser.add_argument("--only", type=int, nargs="*", help="block numbers (1-based)")
    parser.add_argument(
        "--require-cargo",
        action="store_true",
        help="fail blocks that need cargo when it is missing, instead of skipping",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.max_replicates is not None and args.max_replicates < 1:
        parser.error("--max-replicates must be positive")
    found = blocks(args.manual.read_text(encoding="utf-8"))
    flags = ["-O"] * sys.flags.optimize
    environment = dict(os.environ)
    environment.setdefault("CRN_CRAFT_SOURCE", str(ROOT))
    cargo = shutil.which("cargo") is not None
    report = []
    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        for number, block in enumerate(found, 1):
            if args.only and number not in args.only:
                continue
            path = scratch / f"block_{number:02d}.py"
            path.write_text(block["code"], encoding="utf-8")
            argv = []
            compare = block["shown"] is not None and args.max_replicates is None
            reason = "no output shown" if block["shown"] is None else "capped run"
            if "sys.argv" in block["code"]:
                work = scratch / f"work_{number:02d}"
                binary = [str(args.binary.resolve())] if args.binary else []
                argv = [str(work), *binary]
                if args.binary is None and compare:
                    compare, reason = False, "needs --binary"
            if args.max_replicates is None:
                command = [sys.executable, *flags, str(path), *argv]
            else:
                launcher = CAPPED.replace("CAP", str(args.max_replicates))
                command = [sys.executable, *flags, "-c", launcher, str(path), *argv]
            entry = {
                "block": number,
                "line": block["line"],
                "heading": block["heading"],
            }
            start = time.perf_counter()
            if "compile_loop" in block["code"] and not cargo:
                status = "failed" if args.require_cargo else "skipped"
                entry |= {"status": status, "detail": "needs cargo (https://rustup.rs)"}
                print(f"block {number:2d}: {status}, it needs cargo", file=sys.stderr)
                report.append(entry)
                continue
            try:
                done = subprocess.run(
                    command,
                    cwd=ROOT,
                    env=environment,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                )
            except subprocess.TimeoutExpired:
                entry |= {"status": "timeout", "detail": f"over {args.timeout:g} s"}
            else:
                if done.returncode != 0:
                    tail = done.stderr.strip().splitlines()[-1:] or [""]
                    entry |= {"status": "failed", "detail": tail[0]}
                elif not compare:
                    entry |= {"status": "ran", "output": reason}
                elif normalize(done.stdout) == normalize(block["shown"]):
                    entry |= {"status": "passed", "output": "matches the manual"}
                else:
                    shown, actual = normalize(block["shown"]), normalize(done.stdout)
                    differ = [
                        {"shown": a, "printed": b}
                        for a, b in zip(shown, actual, strict=False)
                        if a != b
                    ]
                    if len(shown) != len(actual):
                        differ.append({"shown_lines": len(shown), "lines": len(actual)})
                    entry |= {"status": "mismatch", "differences": differ[:5]}
            seconds = time.perf_counter() - start
            print(
                f"block {number:2d} (line {block['line']:4d}, {block['heading']}): "
                f"{entry['status']} in {seconds:.1f} s",
                file=sys.stderr,
            )
            if entry["status"] in ("failed", "mismatch", "timeout"):
                print(json.dumps(entry, indent=2), file=sys.stderr)
            report.append(entry)
    passed = all(e["status"] in ("passed", "ran", "skipped") for e in report)
    document = {
        "schema": "crn-craft.manual-check/0.1",
        "manual": args.manual.resolve().relative_to(ROOT).as_posix()
        if args.manual.resolve().is_relative_to(ROOT)
        else str(args.manual),
        "capped_replicates": args.max_replicates,
        "passed": passed,
        "blocks": report,
    }
    text = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(json.dumps({"passed": passed, "blocks": len(report)}))
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
