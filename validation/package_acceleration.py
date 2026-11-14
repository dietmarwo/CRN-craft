"""Package committed SSA campaign changes; no untracked files or build outputs."""

import argparse
import hashlib
import json
from pathlib import Path
import posixpath
import re
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = "baed0e4590f93e4c1dcc77cd7f039a2604a1594c"
PREVIOUS = "dbb680944301b620d3ad08070442a3723215d46a"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def standalone(text):
    text = re.sub(
        r"!\[[^\]]*\]\(([^)]+)\)",
        lambda m: (
            f"*Chart: `repo/{posixpath.normpath(posixpath.join('docs/evidence', m[1]))}` in the patch ZIP; all values also appear in the tables.*"
        ),
        text,
    )

    def link(match):
        name, target = match.groups()
        if "://" in target or target.startswith("#"):
            return match[0]
        return f"{name} (`repo/{posixpath.normpath(posixpath.join('docs/evidence', target))}` in the patch ZIP)"

    return re.sub(r"(?<!!)\[([^\]]+)\]\(([^)]+)\)", link, text)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if git("status", "--porcelain").strip():
        raise RuntimeError("commit the complete campaign before packaging")
    head = git("rev-parse", "HEAD").decode().strip()
    stage = out / "CRN-craft-ssa-acceleration-package"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    patches = [
        ("from-github-baed0e4.patch", BASE),
        ("from-previous-comparison-dbb6809.patch", PREVIOUS),
    ]
    for filename, base in patches:
        (stage / filename).write_bytes(
            git("diff", "--binary", "--full-index", "--no-ext-diff", base, head)
        )
    changed = (
        git("diff", "--name-only", "--diff-filter=ACMRT", BASE, head)
        .decode()
        .splitlines()
    )
    for name in changed:
        dest = stage / "repo" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(git("show", f"{head}:{name}"))
    readme = f"""# CRN-craft five-step SSA acceleration patch

Head commit: `{head}`. GitHub base: `{BASE}`.
Previous standalone Julia comparison base: `{PREVIOUS}`.
No remote was pushed and no package was published.

## Apply one patch

Use a clean checkout. If starting from GitHub main as fetched for this campaign:

```sh
git clone https://github.com/dietmarwo/CRN-craft.git
cd CRN-craft
git checkout -b ssa-acceleration {BASE}
git apply --check /absolute/path/from-github-baed0e4.patch
git apply --index /absolute/path/from-github-baed0e4.patch
```

If the previous comparison patch at `{PREVIOUS}` is already present, use
`from-previous-comparison-dbb6809.patch` instead. **Do not apply both patches.**
The first patch includes the earlier comparison prerequisites. The `repo/`
folder contains changed/added files for inspection; it is not a full checkout.
Both patches are verified by applying to clean base trees and comparing the
resulting Git tree with the head tree. `patch-verification.json` records this.

## Contents and validation

- `repo/docs/evidence/comparison-catalyst-jumpprocesses.md`: the new seven-workload comparison.
- `repo/docs/campaign/ssa-acceleration-implementation-log-2026-11-02.md`: implementation loops, tests, failures and limitations.
- `repo/docs/alternative-ssa.md`: new request/result contracts, Python and CLI usage, numerical policy.
- `repo/validation/ssa-acceleration/`: raw measurements, compressed samples, final/failed check logs and pinned environments.
- `repo/validation/ssa-acceleration/README.md`: full reproduction commands.
- `repo/THIRD_PARTY_NOTICES.md`: MIT attribution for the adapted exponential sampler.

Run `cargo test --release --workspace --locked` and the acceptance commands in
the evidence README after applying. Build the Python wheel from this source;
installing an older package from PyPI will not supply these methods.

`SHA256SUMS` hashes every other packaged file. ZIP integrity was checked before
delivery. The old Direct/cached schemas and numerical results remain unchanged;
new methods use standalone `ssa-result/0.1` JSON and separate sampler identities.
"""
    (stage / "README.md").write_text(readme)
    (stage / "LICENSE").write_bytes(git("show", f"{head}:LICENSE"))
    (stage / "repo" / "LICENSE").write_bytes(git("show", f"{head}:LICENSE"))
    # Verify patches independently; compare tree objects, not just line counts.
    verification = []
    for filename, base in patches:
        tree = out / ("verify-" + base[:7])
        subprocess.run(
            ["git", "worktree", "add", "--detach", str(tree), base],
            cwd=ROOT,
            check=True,
            capture_output=True,
        )
        try:
            subprocess.run(
                ["git", "apply", "--check", str(stage / filename)],
                cwd=tree,
                check=True,
                capture_output=True,
            )
            subprocess.run(
                ["git", "apply", "--index", str(stage / filename)],
                cwd=tree,
                check=True,
                capture_output=True,
            )
            actual = (
                subprocess.check_output(["git", "write-tree"], cwd=tree)
                .decode()
                .strip()
            )
            expected = git("rev-parse", f"{head}^{{tree}}").decode().strip()
            if actual != expected:
                raise RuntimeError("patch tree differs")
            verification.append(
                dict(
                    patch=filename,
                    base=base,
                    head=head,
                    expected_tree=expected,
                    actual_tree=actual,
                    passed=True,
                )
            )
        finally:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(tree)],
                cwd=ROOT,
                check=True,
                capture_output=True,
            )
    (stage / "patch-verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    members = sorted(p for p in stage.rglob("*") if p.is_file())
    (stage / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(stage).as_posix()}\n"
            for p in members
        )
    )
    target = out / "CRN-craft-ssa-acceleration-patch.zip"
    with zipfile.ZipFile(
        target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as z:
        for path in sorted(p for p in stage.rglob("*") if p.is_file()):
            info = zipfile.ZipInfo(
                path.relative_to(stage).as_posix(), date_time=(2026, 11, 2, 0, 0, 0)
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            z.writestr(info, path.read_bytes())
    with zipfile.ZipFile(target) as z:
        if z.testzip() is not None:
            raise RuntimeError("ZIP CRC failure")
    log = out / "CRN-craft-implementation-log.md"
    log.write_bytes(
        (
            ROOT / "docs/campaign/ssa-acceleration-implementation-log-2026-11-02.md"
        ).read_bytes()
    )
    log.write_text(
        log.read_text()
        + f"\n## Delivered patch verification\n\nBoth base patches apply cleanly and reconstruct Git tree `{verification[0]['actual_tree']}` for head `{head}`. The ZIP passed its CRC check; per-file SHA256 hashes and the full verification report are included.\n"
    )
    comparison = out / "comparison-catalyst-jumpprocesses.md"
    comparison.write_text(
        standalone(
            (ROOT / "docs/evidence/comparison-catalyst-jumpprocesses.md").read_text()
        )
    )
    print(
        json.dumps(
            dict(
                head=head,
                files=[str(target), str(log), str(comparison)],
                zip_bytes=target.stat().st_size,
                zip_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                patches=verification,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
