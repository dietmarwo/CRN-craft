# Releasing

One tag publishes every registry artifact with the workspace version:

| Registry | Package | Workflow |
|---|---|---|
| PyPI | `crn-craft` (abi3 wheels for Linux x86_64, Windows x64, macOS x86_64 and arm64, plus the source distribution) | `python-release.yml` |
| crates.io / docs.rs | `crn-craft`, `crn-craft-rules`, `crn-craft-cli` (binary `crn-craft`) | `publish-crates.yml` |

The Python bindings (`bindings/python`), the diffsol integration and the job
worker are separate workspaces with `publish = false`. The version is set in
`[workspace.package]` of `Cargo.toml`, in `[package]` of
`bindings/python/Cargo.toml` and in `pyproject.toml`;
`scripts/package_version.py` fails when they disagree (Cargo's `0.3.0-dev.1`
equals Python's `0.3.0.dev1`). `scripts/check_release.py` adds the crate
checks: publishable members inheriting the version, a description, a README
and a LICENSE identical to the root one, versioned path dependencies, and the
version recorded in every lockfile that contains a CRN-craft crate.
`crn_craft.__version__` reports the installed metadata version and
`crn_craft._native.__version__` the extension's Cargo version. Publication is
irreversible: never reuse a version. MIT project code does not relicense
dependencies or reference datasets; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## One-time setup

1. Create a GitHub environment named `release` (protect it if releases
   should need approval).
2. PyPI: create a pending Trusted Publisher for project `crn-craft`: owner
   `dietmarwo`, repository `CRN-craft`, workflow `python-release.yml`,
   environment `release`. The first tag then creates the project.
3. crates.io: the first version must be published with an API token, because
   crates.io accepts a Trusted Publisher only for an existing crate. On
   crates.io, Account Settings → API Tokens, create a token with the scopes
   `publish-new` and `publish-update`, restricted to `crn-craft*`, with a short
   expiry, and add it as the secret `CARGO_REGISTRY_TOKEN` of the `release`
   environment. The tag then publishes the crates in dependency order
   (`crn-craft`, `crn-craft-rules`, `crn-craft-cli`), skipping versions already
   on crates.io, so rerunning the job completes a partial upload. A tag run
   without a token and without Trusted Publishing fails instead of skipping.
   After the first upload, delete the secret, revoke the token, configure a
   Trusted Publisher for each crate (owner `dietmarwo`, repository
   `CRN-craft`, workflow `publish-crates.yml`, environment `release`) and set
   the repository variable `CRATES_IO_TRUSTED_PUBLISHING` to `true`.

## Pre-release checks

```bash
python3 scripts/check_release.py
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
python validation/check_contracts.py --output contracts.json
cargo publish -p crn-craft -p crn-craft-rules -p crn-craft-cli --dry-run --locked
```

The dry run packages the three crates and builds each from its own archive
(current stable Cargo stages the unpublished dependencies). A packaged crate
carries its own normalized `Cargo.lock`, which is what its archives'
`cargo_lock_sha256` records.

Build the wheel with `maturin==1.9.6 build --release --locked`, install it in
a fresh environment and run `validation/check_installed.py --output
installed.json` with that environment's Python from outside the checkout. The
checker runs the `crn-craft` and `crn-craft-import` console scripts from that
environment's scripts directory (`sysconfig.get_path("scripts")`), never from
PATH, and requires `crn-craft --version` to report the Python interface and the
installed version. For the source distribution, build it with `maturin sdist`,
install it in a fresh environment (this compiles the extension) and run the
same checker.

`publish-crates.yml` runs the Rust checks and the dry run on pull requests
that touch the crates and on every tag. `python-release.yml` does the Python
checks on every platform: it builds the four
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
`python-release.yml` and `publish-crates.yml`; both refuse a tag that does not
equal `v` plus the version. From a computer:

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

For the crates: `cargo install crn-craft-cli --version "$version" --locked`
and `crn-craft describe`, and a new project with `cargo add crn-craft@"$version"`
that builds. Check the PyPI page (its description is `python/README.md`; run its
example), the crate pages and docs.rs,
then the shared catalog of
[optimization-jobs](https://github.com/dietmarwo/optimization-jobs/blob/main/docs/discovery/README.md)
and the tutorial requirements that build CRN-craft from source.
