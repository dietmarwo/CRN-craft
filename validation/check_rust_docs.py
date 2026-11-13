"""Compile and run the `rust` blocks of docs/developer/rust-api.md.

Each block is a complete program. The script writes every block as a binary
of a scratch crate that depends on the workspace crate `crn-craft` by path,
builds the crate with cargo and runs each binary, so a block that does not
compile, or whose assertions fail, is reported. The scratch crate lives in
`--workdir` (a new temporary directory by default), never in the repository.

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

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "developer" / "rust-api.md"
BLOCK = re.compile(r"^```rust\n(.*?)^```$", re.S | re.M)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--page", type=Path, default=PAGE)
    parser.add_argument(
        "--workdir", type=Path, help="scratch directory (must not exist)"
    )
    parser.add_argument("--cargo", default="cargo")
    args = parser.parse_args()

    blocks = BLOCK.findall(args.page.read_text(encoding="utf-8"))
    if not blocks:
        print(f"{args.page}: no rust blocks", file=sys.stderr)
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
    for number, code in enumerate(blocks, 1):
        name = f"block{number}"
        (work / "src" / "bin" / f"{name}.rs").write_text(code, encoding="utf-8")
        names.append(name)
    build = subprocess.run([args.cargo, "build", "--bins"], cwd=work)
    if build.returncode != 0:
        print("a rust block does not compile", file=sys.stderr)
        return 1
    failed = 0
    for name in names:
        run = subprocess.run([args.cargo, "run", "--quiet", "--bin", name], cwd=work)
        print(f"{name}: {'ok' if run.returncode == 0 else 'FAILED'}")
        failed += run.returncode != 0
    if not args.workdir:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(names)} rust blocks, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
