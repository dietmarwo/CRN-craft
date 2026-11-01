# Releasing

This is the `0.2.0-dev.1` campaign delivery (milestones M3–M10 and performance iterations P01–P20), with `publish = false` for every crate. The Python distribution is `crn-craft` `0.2.0.dev1`; `crn_craft.__version__` reports the installed metadata version and the native core's Cargo version is `crn_craft._native.__version__`. This remains a reviewable development source delivery. Public package release still requires owner approval and installed-package/platform acceptance.

For a future release: complete the declared milestone, regenerate/version schemas, run the documented checks on supported toolchains/platforms, review dependency licenses and reference-data provenance, update changelog/capabilities, and run installed-package acceptance. Preserve lockfiles and evidence. Package availability must be checked at release time; earlier planning searches do not reserve names.

The owner decides when to push, tag and publish. No release credentials or workflow that automatically publishes packages are included. MIT project code does not relicense dependencies or reference datasets.

## Campaign package checks

Build the isolated Python package with `maturin==1.9.6 build --release --locked`.
Use `maturin sdist`, extract it into a new directory, and rebuild from that
project directory to exercise its pyproject metadata and package resources.
Install the wheel in a fresh Python environment and run
`validation/check_installed.py --output installed.json` with that environment's
Python. The checker runs the `crn-craft` and `crn-craft-import` console scripts
from that environment's scripts directory (`sysconfig.get_path("scripts")`),
never from PATH, and requires `crn-craft --version` to report the Python
interface and the installed version. To check the sdist on its own, build it
with `maturin sdist`, extract it into a new empty directory, build a wheel there
in a fresh environment and run the same checker against that install. No
package index publication is part of these checks.

The `python` CI job builds the wheel on Linux and runs `check_python.py`,
`check_importers.py` and `check_installed.py`, each normally and under `-O`.

The optional wheel workflow is manual and covers Linux, macOS and Windows.
Linux x86_64/Python 3.12 is locally executed; an ABI3 tag alone is not proof of
runtime acceptance on every Python/OS version. Build the native CLI independently
with the root Cargo.lock. Optional jobs need Rust 1.88; optional diffsol needs
Rust 1.89. Preserve all three workspace lockfiles and the Python pins.
