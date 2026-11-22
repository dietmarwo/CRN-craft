"""Fail closed on new dependencies, known native build tools and links metadata.

The core workspace (``Cargo.toml``) is audited against
``dependency-policy.json``; the Python extension (``bindings/python``, its
own lockfile) against ``dependency-policy-python.json``, which also records
the declared license of every package and the two ``links`` keys of PyO3.
"""

import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


FORBIDDEN = {
    "cc",
    "cmake",
    "bindgen",
    "llvm-sys",
    "openblas-src",
    "lapack-sys",
    "blas-src",
    "lapack-src",
    "intel-mkl-src",
    "rayon",
}


def audit(packages, allowed, allowed_local, licenses=None, links=None):
    """Errors for packages outside the policy. With ``licenses``, each
    registry package must declare the recorded license; ``links`` maps the
    packages allowed to carry ``links`` metadata to its value."""
    errors = []
    links = links or {}
    for package in packages:
        name = package["name"]
        link = package.get("links")
        if name in FORBIDDEN or (link and links.get(name) != link):
            errors.append(f"native build/link dependency: {name}")
        source = package.get("source")
        if source is None:
            if name not in allowed_local:
                errors.append(f"unreviewed local dependency: {name}")
        elif not source.startswith("registry+"):
            errors.append(f"unreviewed dependency source: {name}")
        elif name not in allowed:
            errors.append(f"unreviewed registry dependency: {name}")
        elif licenses is not None and package.get("license") != licenses[name]:
            errors.append(f"license differs from the policy: {name}")
    return errors


def metadata(cargo, manifest):
    return json.loads(
        subprocess.check_output(
            [
                cargo,
                "metadata",
                "--locked",
                "--format-version",
                "1",
                "--manifest-path",
                str(ROOT / manifest),
            ],
            cwd=ROOT,
            timeout=120,
        )
    )["packages"]


def listing(packages):
    return sorted(
        [
            {
                "name": p["name"],
                "version": p["version"],
                "links": p["links"],
            }
            for p in packages
        ],
        key=lambda p: (p["name"], p["version"]),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cargo", default="cargo")
    args = parser.parse_args()
    policy = json.loads((ROOT / "validation/dependency-policy.json").read_text())
    allowed = set(policy["allowed_registry_packages"])
    allowed_local = set(policy["allowed_workspace_packages"])
    packages = metadata(args.cargo, "Cargo.toml")
    errors = audit(packages, allowed, allowed_local)
    python = json.loads((ROOT / "validation/dependency-policy-python.json").read_text())
    python_licenses = python["allowed_registry_packages"]
    python_args = (
        set(python_licenses),
        set(python["allowed_workspace_packages"]),
        python_licenses,
        python["allowed_links"],
    )
    python_packages = metadata(args.cargo, python["manifest"])
    errors += [
        f"bindings/python: {e}" for e in audit(python_packages, *python_args)
    ]  # Checker mutation probes, independent of cargo's actual clean tree.
    for bad in [
        {"name": "cmake", "source": "registry+fixture", "links": None},
        {"name": "serde", "source": "registry+fixture", "links": "hidden_c_library"},
        {"name": "surprise", "source": "registry+fixture", "links": None},
        {"name": "serde", "source": "git+fixture", "links": None},
        {"name": "local_surprise", "source": None, "links": None},
    ]:
        if not audit([bad], allowed, allowed_local):
            errors.append(f"dependency gate failed its mutation probe: {bad['name']}")
    for bad in [
        {"name": "rayon", "source": "registry+fixture", "links": None},
        {"name": "faer", "source": "registry+fixture", "license": "GPL-3.0"},
        {"name": "pyo3", "source": "registry+fixture", "links": "other"},
    ]:
        if not audit([{"links": None, **bad}], *python_args):
            errors.append(f"binding gate failed its mutation probe: {bad['name']}")
    print(
        json.dumps(
            {
                "schema": "crn-craft.dependency-audit/0.1",
                "passed": not errors,
                "errors": errors,
                "mutation_probes": 8,
                "packages": listing(packages),
                "python_extension_packages": listing(python_packages),
            },
            indent=2,
        )
    )
    if errors:
        raise RuntimeError("dependency audit failed")


if __name__ == "__main__":
    main()
