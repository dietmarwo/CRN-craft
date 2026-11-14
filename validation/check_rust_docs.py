"""Compile and run the `rust` blocks of the Rust pages; check their `toml` blocks.

The pages are docs/developer/rust-api.md and the READMEs of the three
crates (`--page` names others). In them:

- A `rust` block is a complete program. The script writes every block as a
  binary of a scratch crate that depends on the workspace crate `crn-craft`
  by path, builds the crate with cargo and runs each binary, so a block that
  does not compile, or whose assertions fail, is reported. A block that is
  a fragment (it uses names it does not define) is marked with the comment
  `<!-- check: fragment -->` on the line before its fence; it is listed as
  a fragment and not compiled.
- A `toml` block is valid TOML, and a dependency on `crn-craft` or
  `crn-craft-rules` names a version that the workspace version satisfies
  (`"0.7"` for 0.7.0).
- An `sh` block passes `bash -n`.

Any other fence language is an error, and at least one complete `rust`
block must be compiled. The scratch crate lives in `--workdir` (a new
temporary directory by default), never in the repository.

    flock /home/claude/craft5/cpu.lock python validation/check_rust_docs.py

Set CARGO_TARGET_DIR to reuse a build directory. The script needs network
access to crates.io the first time, like any cargo build. It copies the
workspace Cargo.lock so that dependency versions match the workspace.
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
PAGES = [
    ROOT / "docs" / "developer" / "rust-api.md",
    ROOT / "crates" / "crn-craft" / "README.md",
    ROOT / "crates" / "crn-craft-cli" / "README.md",
    ROOT / "crates" / "crn-craft-rules" / "README.md",
]
FENCE = re.compile(r"^(`{3,})(\w*)\s*$")
FRAGMENT = re.compile(r"^<!--\s*check:\s*fragment\s*-->\s*$")
CRATES = ("crn-craft", "crn-craft-rules")


def blocks(text):
    """Yield (language, line, code, fragment) of every fenced block."""
    current, previous = None, ""
    for number, line in enumerate(text.splitlines(), 1):
        if current is None:
            match = FENCE.match(line)
            if match:
                fragment = bool(FRAGMENT.match(previous))
                current = (match.group(2), number, len(match.group(1)), fragment, [])
            elif line.strip():
                previous = line
        else:
            stripped = line.strip()
            if stripped and stripped == "`" * len(stripped):
                if len(stripped) >= current[2]:
                    language, start, _, fragment, lines = current
                    yield language, start, "\n".join(lines) + "\n", fragment
                    current, previous = None, ""
                    continue
            current[4].append(line)
    if current is not None:
        raise ValueError(f"the fence at line {current[1]} is never closed")


def satisfies(requirement, version):
    """A plain or caret requirement such as "0.7" or "^0.7.0"."""
    wanted = requirement.lstrip("^").split(".")
    return version.split(".")[: len(wanted)] == wanted


def check_toml(code, where, version):
    try:
        document = tomllib.loads(code)
    except tomllib.TOMLDecodeError as error:
        return [f"{where}: not valid TOML: {error}"]
    problems = []
    for name in CRATES:
        dependency = document.get("dependencies", {}).get(name)
        if dependency is None:
            continue
        requirement = (
            dependency if isinstance(dependency, str) else dependency.get("version")
        )
        if requirement is None or not satisfies(requirement, version):
            problems.append(
                f"{where}: {name} = {requirement!r} does not match the workspace "
                f"version {version}"
            )
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--page", type=Path, action="append", help="default: all")
    parser.add_argument(
        "--workdir", type=Path, help="scratch directory (must not exist)"
    )
    parser.add_argument("--cargo", default="cargo")
    args = parser.parse_args()

    workspace = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    version = workspace["workspace"]["package"]["version"]
    programs, problems = [], []
    for page in args.page or PAGES:
        relative = page.resolve().relative_to(ROOT).as_posix()
        counts = {"rust": 0, "fragment": 0, "toml": 0, "sh": 0}
        for language, line, code, fragment in blocks(page.read_text(encoding="utf-8")):
            where = f"{relative}:{line}"
            if language == "rust" and fragment:
                counts["fragment"] += 1
            elif language == "rust":
                programs.append((where, code))
                counts["rust"] += 1
            elif language == "toml":
                problems += check_toml(code, where, version)
                counts["toml"] += 1
            elif language == "sh":
                done = subprocess.run(
                    ["bash", "-n", "-c", code], capture_output=True, text=True
                )
                if done.returncode:
                    problems.append(f"{where}: sh syntax error: {done.stderr.strip()}")
                counts["sh"] += 1
            else:
                problems.append(
                    f"{where}: unchecked fence language {language!r}; "
                    "use rust, toml or sh"
                )
        print(
            f"{relative}: {counts['rust']} rust programs, {counts['fragment']} "
            f"fragments, {counts['toml']} toml, {counts['sh']} sh"
        )
    if not programs:
        problems.append("no complete rust block to compile")
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return 1

    work = args.workdir or Path(tempfile.mkdtemp(prefix="crn-rust-docs-"))
    if args.workdir:
        work.mkdir(parents=True)
    (work / "src" / "bin").mkdir(parents=True)
    crate = ROOT / "crates" / "crn-craft"
    (work / "Cargo.toml").write_text(
        '[package]\nname = "rust-api-blocks"\nversion = "0.0.0"\nedition = "2024"\n'
        f'publish = false\n\n[dependencies]\ncrn-craft = {{ path = "{crate}" }}\n',
        encoding="utf-8",
    )
    shutil.copy(ROOT / "Cargo.lock", work / "Cargo.lock")
    names = []
    for number, (where, code) in enumerate(programs, 1):
        name = f"block{number}"
        (work / "src" / "bin" / f"{name}.rs").write_text(code, encoding="utf-8")
        names.append((name, where))
    build = subprocess.run([args.cargo, "build", "--bins"], cwd=work)
    if build.returncode != 0:
        print("a rust block does not compile", file=sys.stderr)
        return 1
    failed = 0
    for name, where in names:
        run = subprocess.run([args.cargo, "run", "--quiet", "--bin", name], cwd=work)
        print(f"{name} ({where}): {'ok' if run.returncode == 0 else 'FAILED'}")
        failed += run.returncode != 0
    if not args.workdir:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(names)} rust blocks, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
