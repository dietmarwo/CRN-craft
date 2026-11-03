# Releasing

One tag publishes the Python package with the workspace version:

| Registry | Package | Workflow |
|---|---|---|
| PyPI | `crn-craft` (abi3 wheels for Linux x86_64, Windows x64, macOS x86_64 and arm64, plus the source distribution) | `python-release.yml` |

The Rust crates are not published (`publish = false`); the native CLI, the
rules engine, the diffsol integration and the job worker build from this
repository. The version is set in `[workspace.package]` of `Cargo.toml`, in
`[package]` of `bindings/python/Cargo.toml` and in `pyproject.toml`;
`scripts/package_version.py` fails when they disagree (Cargo's `0.3.0-dev.1`
equals Python's `0.3.0.dev1`). `crn_craft.__version__` reports the installed
metadata version and `crn_craft._native.__version__` the extension's Cargo
version. Publication is irreversible: never reuse a version. MIT project code
does not relicense dependencies or reference datasets; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## One-time setup

1. Create a GitHub environment named `release` (protect it if releases
   should need approval).
2. PyPI: create a pending Trusted Publisher for project `crn-craft`: owner
   `dietmarwo`, repository `CRN-craft`, workflow `python-release.yml`,
   environment `release`. The first tag then creates the project.

## Pre-release checks

```bash
python3 scripts/package_version.py
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
python validation/check_contracts.py --output contracts.json
```

Build the wheel with `maturin==1.9.6 build --release --locked`, install it in
a fresh environment and run `validation/check_installed.py --output
installed.json` with that environment's Python from outside the checkout. The
checker runs the `crn-craft` and `crn-craft-import` console scripts from that
environment's scripts directory (`sysconfig.get_path("scripts")`), never from
PATH, and requires `crn-craft --version` to report the Python interface and the
installed version. For the source distribution, build it with `maturin sdist`,
install it in a fresh environment (this compiles the extension) and run the
same checker.

`python-release.yml` does all of this on every platform: it builds the four
wheels and the source distribution and runs the checker on CPython 3.10 and
3.14, normally and under `-O`. It runs on pull requests that change the
packaging and manually (Actions, "Run workflow"); only a tag publishes.

The release commit must have a clean working tree, the intended version in
all three places and in every `Cargo.lock` that contains a CRN-craft crate
(root, `bindings/python`, `integrations/diffsol`, `integrations/jobs`,
`validation/parameter-optimization/worker`), a dated `CHANGELOG.md` entry with
its compare link, and passing CI.

## Tag and publish

Merge the release pull request, then GitHub → Releases → Draft a new release
→ new tag `vX.Y.Z` on `main` → publish. The tag push starts
`python-release.yml`, which refuses a tag that does not equal `v` plus the
version. From a computer:

```bash
version="$(python3 scripts/package_version.py)"
git tag -a "v$version" -m "Release CRN-craft $version"
git push origin "v$version"
```

## Post-release verification

```bash
python -m venv /tmp/crn-craft-check
/tmp/crn-craft-check/bin/python -m pip install crn-craft
cd /tmp && /tmp/crn-craft-check/bin/python -I /path/to/CRN-craft/validation/check_installed.py --output installed.json
```

Check the PyPI page (its description is `python/README.md`; run its example),
then the shared catalog of
[optimization-jobs](https://github.com/dietmarwo/optimization-jobs/blob/main/docs/discovery/README.md)
and the tutorial requirements that build CRN-craft from source.
