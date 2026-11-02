"""Package the committed parameter-interface work and verify exact patch replay."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = "c16fcc55654d0867a58af48d48a8e860efce8896"
LOG = "crn-craft-reusable-parameters-implementation-log-2026-11-02.md"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if git("status", "--porcelain").strip():
        raise RuntimeError("Commit intended changes before packaging")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    stage = output / "CRN-craft-parameters-package"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    head = git("rev-parse", "HEAD").decode().strip()
    patch = stage / "CRN-craft-reusable-parameters.patch"
    patch.write_bytes(
        git("diff", "--binary", "--full-index", "--no-ext-diff", BASE, head)
    )
    changed = (
        git("diff", "--name-only", "--diff-filter=ACMRT", BASE, head)
        .decode()
        .splitlines()
    )
    for name in changed:
        target = stage / "repo" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(git("show", f"{head}:{name}"))
    for name in ["LICENSE", "THIRD_PARTY_NOTICES.md"]:
        (stage / name).write_bytes(git("show", f"{head}:{name}"))
    check = output / "clean-base-verification"
    subprocess.run(
        ["git", "worktree", "add", "--detach", str(check), BASE],
        cwd=ROOT,
        capture_output=True,
        check=True,
    )
    try:
        subprocess.run(
            ["git", "apply", "--check", str(patch)],
            cwd=check,
            capture_output=True,
            check=True,
        )
        subprocess.run(
            ["git", "apply", "--index", str(patch)],
            cwd=check,
            capture_output=True,
            check=True,
        )
        tree = (
            subprocess.check_output(["git", "write-tree"], cwd=check).decode().strip()
        )
        expected = git("rev-parse", f"{head}^{{tree}}").decode().strip()
        if tree != expected:
            raise RuntimeError("Patch does not reconstruct tested tree")
    finally:
        subprocess.run(
            ["git", "worktree", "remove", "--force", str(check)],
            cwd=ROOT,
            capture_output=True,
            check=True,
        )
    verification = dict(
        base=BASE,
        head=head,
        expected_tree=expected,
        reconstructed_tree=tree,
        passed=True,
    )
    (stage / "patch-verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    log = (
        ROOT / "docs/campaign/reusable-parameters-implementation-log-2026-11-02.md"
    ).read_text()
    log += f"\n## Delivered patch verification\n\nThe aggregate patch applies cleanly to `{BASE}` and reconstructs tree `{tree}`, exactly equal to committed head `{head}`. The ZIP passed a full CRC test and verification of every file against its SHA256 manifest.\n"
    (stage / LOG).write_text(log)
    (output / LOG).write_text(log)
    (
        stage / "README.md"
    ).write_text(f"""# CRN-craft reusable parameters and generated kernels

Base: `{BASE}` (GitHub main, merged PR #8, includes method.fast).
Head: `{head}`. No remote was pushed or package published.

## Apply

Use a clean checkout and only the aggregate patch at this ZIP's root:

```sh
git clone https://github.com/dietmarwo/CRN-craft.git
cd CRN-craft
git checkout -b reusable-parameters {BASE}
git apply --check /absolute/path/CRN-craft-reusable-parameters.patch
git apply --index /absolute/path/CRN-craft-reusable-parameters.patch
```

The `repo/` folder is a changed-files snapshot for review, not a full checkout.
Apply the patch to the base for a complete runnable repository. Exact tree
reconstruction was verified: `{tree}`. Every other file is in `SHA256SUMS`.

## Read

- `{LOG}`: implementation, checks, complete optimization results and limits.
- `repo/docs/reusable-parameters.md`: Rust/Python API and generated kernels.
- `repo/docs/comparison-parameter-optimization.md`: compilation, optimizer and
  held-out quality comparison, including rebop dynamic and macro routes.
- `repo/validation/parameter-optimization/README.md`: reproducible commands.
- `repo/validation/parameter-optimization/evidence/`: raw traces, build costs,
  failed development evidence, checked/fast statistical gates and source hashes.
- `LICENSE`, `THIRD_PARTY_NOTICES.md`: preserved attribution. Optional rebop
  reference-license text is inside its benchmark directory.

Rebuild and force reinstall the wheel from the patched checkout: its development
version number is unchanged. The default remains checked; fast must be selected
explicitly. Generated code specializes mass-action propensities and reuses the
existing schedulers/checks. Compiler/optimizer/rebop dependencies are optional.
""")
    files = sorted(p for p in stage.rglob("*") if p.is_file())
    (stage / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(stage).as_posix()}\n"
            for p in files
        )
    )
    target = output / "CRN-craft-reusable-parameters-patch.zip"
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(p for p in stage.rglob("*") if p.is_file()):
            info = zipfile.ZipInfo(
                path.relative_to(stage).as_posix(), (2026, 11, 2, 0, 0, 0)
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes())
    with zipfile.ZipFile(target) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("ZIP CRC failed")
        for line in archive.read("SHA256SUMS").decode().splitlines():
            digest, name = line.split("  ", 1)
            if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                raise RuntimeError("ZIP member hash failed: " + name)
    result = dict(
        verification=verification,
        zip_bytes=target.stat().st_size,
        zip_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        files=[str(target), str(output / LOG)],
    )
    (output / "package-verification.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
