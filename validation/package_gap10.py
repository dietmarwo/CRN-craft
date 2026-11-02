"""Package the committed ten-iteration campaign and verify its aggregate patch."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = "add6cca0fb146e8a40418432e963219514086024"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if git("status", "--porcelain").strip():
        raise RuntimeError("commit all intended changes before packaging")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    stage = out / "CRN-craft-gap10-package"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    head = git("rev-parse", "HEAD").decode().strip()
    patch = stage / "CRN-craft-julia-gap10.patch"
    patch.write_bytes(
        git("diff", "--binary", "--full-index", "--no-ext-diff", BASE, head)
    )
    for name in (
        git("diff", "--name-only", "--diff-filter=ACMRT", BASE, head)
        .decode()
        .splitlines()
    ):
        target = stage / "repo" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(git("show", f"{head}:{name}"))
    for name in ["LICENSE", "THIRD_PARTY_NOTICES.md"]:
        (stage / name).write_bytes(git("show", f"{head}:{name}"))
        (stage / "repo" / name).write_bytes(git("show", f"{head}:{name}"))
    check = out / "patch-base-check"
    subprocess.run(
        ["git", "worktree", "add", "--detach", str(check), BASE],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    try:
        subprocess.run(
            ["git", "apply", "--check", str(patch)],
            cwd=check,
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["git", "apply", "--index", str(patch)],
            cwd=check,
            check=True,
            capture_output=True,
        )
        actual = (
            subprocess.check_output(["git", "write-tree"], cwd=check).decode().strip()
        )
        expected = git("rev-parse", f"{head}^{{tree}}").decode().strip()
        if actual != expected:
            raise RuntimeError("patch tree differs from tested head")
    finally:
        subprocess.run(
            ["git", "worktree", "remove", "--force", str(check)],
            cwd=ROOT,
            check=True,
            capture_output=True,
        )
    verification = dict(
        base=BASE, head=head, expected_tree=expected, actual_tree=actual, passed=True
    )
    (stage / "patch-verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    log = (
        ROOT / "docs/campaign/julia-gap10-implementation-log-2026-11-02.md"
    ).read_text()
    log += f"\n## Delivered patch verification\n\nThe aggregate patch applies cleanly to GitHub base `{BASE}` and reconstructs tree `{actual}`, exactly equal to committed head `{head}`. The ZIP has a per-file SHA256 manifest and passed a complete CRC check.\n"
    logname = "crn-craft-julia-gap10-implementation-log-2026-11-02.md"
    (stage / logname).write_text(log)
    (out / logname).write_text(log)
    (stage / "README.md").write_text(f"""# Ten further CRN-craft Julia-gap iterations

Base: `{BASE}` (GitHub main, merged PR #7).
Head: `{head}`. The previous A01–A10 campaign is already in the base.
This ZIP contains the new C01–C10 campaign. No remote was pushed.

## Apply the final patch

Use a clean checkout and apply only the aggregate patch at the ZIP root:

```sh
git clone https://github.com/dietmarwo/CRN-craft.git
cd CRN-craft
git checkout -b julia-gap10 {BASE}
git apply --check /absolute/path/CRN-craft-julia-gap10.patch
git apply --index /absolute/path/CRN-craft-julia-gap10.patch
```

Only C05, C06 and C09 passed the fixed performance rule. Rejected CXX.patch
files nested under `repo/validation/julia-gap10` are evidence, not additional
patches to apply. The `repo/` snapshot contains only changed/added files for
inspection; it is not a full repository. Apply the aggregate patch to the base
for a complete source checkout.

The clean-base application was verified by comparing Git tree objects:
`{actual}`. `patch-verification.json` records the result. `SHA256SUMS` hashes
every other packaged file and the ZIP passed its CRC check.

## Read and reproduce

- `{logname}`: complete iteration decisions, tests, measured results and limits.
- `repo/docs/comparison-julia-gap10.md`: direct same-host baseline/final/Julia comparison.
- `repo/docs/campaign/julia-gap10-source-investigation.md`: inspected Julia mechanisms.
- `repo/validation/julia-gap10/README.md`: build and reproduction instructions.
- `repo/validation/julia-gap10/`: every observation, final gates, compressed samples and rejected experiment diffs.
- `LICENSE` and `THIRD_PARTY_NOTICES.md`: original MIT terms and retained attribution.

Build the Python wheel from the patched source and force reinstall it; an older
wheel with the same development version will not contain these optimizations.
All final statistical/contract checks are retained. Performance is a measured
single-host result with explicit API boundaries, not a universal speed claim.
""")
    members = sorted(f for f in stage.rglob("*") if f.is_file())
    (stage / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(stage).as_posix()}\n"
            for f in members
        )
    )
    target = out / "CRN-craft-julia-gap10-patch.zip"
    with zipfile.ZipFile(
        target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as z:
        for path in sorted(f for f in stage.rglob("*") if f.is_file()):
            info = zipfile.ZipInfo(
                path.relative_to(stage).as_posix(), date_time=(2026, 11, 2, 0, 0, 0)
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            z.writestr(info, path.read_bytes())
    with zipfile.ZipFile(target) as z:
        if z.testzip() is not None:
            raise RuntimeError("ZIP CRC failed")
        for line in z.read("SHA256SUMS").decode().splitlines():
            digest, name = line.split("  ", 1)
            if hashlib.sha256(z.read(name)).hexdigest() != digest:
                raise RuntimeError("ZIP hash failed: " + name)
    print(
        json.dumps(
            dict(
                verification=verification,
                zip_bytes=target.stat().st_size,
                zip_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                files=[str(target), str(out / logname)],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
