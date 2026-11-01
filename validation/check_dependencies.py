"""Fail closed on new dependencies, known native build tools and links metadata."""

import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def audit(packages, allowed, allowed_local):
    errors = []
    forbidden = {"cc", "cmake", "bindgen", "llvm-sys", "openblas-src", "lapack-sys"}
    for package in packages:
        name = package["name"]
        if name in forbidden or package.get("links"):
            errors.append(f"native build/link dependency: {name}")
        source = package.get("source")
        if source is None:
            if name not in allowed_local:
                errors.append(f"unreviewed local dependency: {name}")
        elif not source.startswith("registry+"):
            errors.append(f"unreviewed dependency source: {name}")
        elif name not in allowed:
            errors.append(f"unreviewed registry dependency: {name}")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cargo", default="cargo")
    args = parser.parse_args()
    policy = json.loads((ROOT / "validation/dependency-policy.json").read_text())
    allowed = set(policy["allowed_registry_packages"])
    allowed_local = set(policy["allowed_workspace_packages"])
    metadata = json.loads(
        subprocess.check_output(
            [args.cargo, "metadata", "--locked", "--format-version", "1"],
            cwd=ROOT,
            timeout=120,
        )
    )
    packages = metadata["packages"]
    errors = audit(packages, allowed, allowed_local)
    # Checker mutation probes, independent of cargo's actual clean tree.
    for bad in [
        {"name": "cmake", "source": "registry+fixture", "links": None},
        {"name": "serde", "source": "registry+fixture", "links": "hidden_c_library"},
        {"name": "surprise", "source": "registry+fixture", "links": None},
        {"name": "serde", "source": "git+fixture", "links": None},
        {"name": "local_surprise", "source": None, "links": None},
    ]:
        if not audit([bad], allowed, allowed_local):
            errors.append(f"dependency gate failed its mutation probe: {bad['name']}")
    print(
        json.dumps(
            {
                "schema": "crn-craft.dependency-audit/0.1",
                "passed": not errors,
                "errors": errors,
                "mutation_probes": 5,
                "packages": sorted(
                    [
                        {
                            "name": p["name"],
                            "version": p["version"],
                            "links": p["links"],
                        }
                        for p in packages
                    ],
                    key=lambda p: (p["name"], p["version"]),
                ),
            },
            indent=2,
        )
    )
    if errors:
        raise RuntimeError("dependency audit failed")


if __name__ == "__main__":
    main()
