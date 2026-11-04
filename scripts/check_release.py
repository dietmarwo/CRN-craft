"""Fail before packaging on inconsistent versions, crate metadata or lockfiles.

Checks what `scripts/package_version.py` checks (the three version sources
agree), plus: every workspace member is publishable, inherits the workspace
version, ships a README, a LICENSE identical to the root one and a description;
path dependencies between the crates name the same version; and every lockfile
that contains a CRN-craft crate records that version. Prints a JSON summary.
"""

import importlib.util
import json
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
LOCKFILES = [
    "Cargo.lock",
    "bindings/python/Cargo.lock",
    "integrations/diffsol/Cargo.lock",
    "integrations/jobs/Cargo.lock",
    "validation/parameter-optimization/worker/Cargo.lock",
]


def read(path):
    return tomllib.loads((ROOT / path).read_text(encoding="utf-8"))


def require(value, message):
    if not value:
        sys.exit("check_release: " + message)


def main():
    spec = importlib.util.spec_from_file_location(
        "package_version", ROOT / "scripts/package_version.py"
    )
    versions = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(versions)
    python_version = versions.versions()  # exits when the sources disagree
    root = read("Cargo.toml")
    version = root["workspace"]["package"]["version"]
    names = []
    for member in root["workspace"]["members"]:
        manifest = read(member + "/Cargo.toml")
        package = manifest["package"]
        name = package["name"]
        require(package.get("publish") is not False, f"{name} has publish = false")
        require(
            package.get("version") == {"workspace": True},
            f"{name} must inherit the workspace version",
        )
        require(package.get("description"), f"{name} has no description")
        require((ROOT / member / "README.md").is_file(), f"{name} has no README.md")
        require(
            (ROOT / member / "LICENSE").read_bytes() == (ROOT / "LICENSE").read_bytes(),
            f"{name}/LICENSE differs from the root LICENSE",
        )
        for dep_name, dep in manifest.get("dependencies", {}).items():
            if isinstance(dep, dict) and dep.get("workspace"):
                dep = root["workspace"]["dependencies"][dep_name]
            if isinstance(dep, dict) and "path" in dep:
                require(
                    dep.get("version") == version,
                    f"{name}: path dependency {dep_name} must name version {version}",
                )
        names.append(name)
    first_party = set(names) | {read("bindings/python/Cargo.toml")["package"]["name"]}
    for lock in LOCKFILES:
        for package in read(lock)["package"]:
            if package["name"] in first_party and "source" not in package:
                require(
                    package["version"] == version,
                    f"{lock} records {package['name']} {package['version']}",
                )
    print(json.dumps({"version": python_version, "crates": names}))


if __name__ == "__main__":
    main()
