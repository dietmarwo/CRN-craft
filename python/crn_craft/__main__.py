"""The crn-craft command of the Python package: the command set of the native CLI.

Every command prints JSON. ``run`` prints exactly what the native
``crn-craft run`` prints; the other reports are the native reports,
pretty-printed with sorted keys.
"""

import json
from pathlib import Path
import sys

from . import _native
from ._common import CrnError, __version__

MILESTONE = "unified-engine"
COMMANDS = ("describe", "check", "run", "verify", "analyze", "convert", "export")

HELP = """\
crn-craft (Python package) — reaction networks: exact and approximate
simulation, distributions, analysis

Usage:
  crn-craft describe [KIND]
  crn-craft check MODEL [REQUEST]
  crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]
  crn-craft verify DIRECTORY
  crn-craft analyze MODEL [STATE]
  crn-craft convert MODEL
  crn-craft export MODEL gillespy2
  crn-craft --version

Aliases: rules = run, rules-check = check.

The commands are those of the native CLI (`cargo install crn-craft-cli`), with
the same reports; `crn-craft --help` of the native CLI describes them.
Models are crn-craft.model/1.0 or crn-craft.rules-model/0.1-0.4; check,
analyze, convert and export also read a network crn-craft.model/0.1.
Requests are crn-craft.run/1.0 (method fsp included) or
crn-craft.rules-run/0.1-0.6. run exits 0 only for a completed ensemble;
verify exits 0 for a valid archive; other errors exit 1.
"""

USAGE_HINT = (
    "Use --help. Commands: describe, check, run, verify, analyze, convert and "
    "export (rules and rules-check are aliases of run and check)."
)


def _usage_error():
    """The native CLI's usage report (CRN-E031, crn-craft.check/1.0)."""
    return {
        "schema": "crn-craft.check/1.0",
        "valid": False,
        "model_sha256": None,
        "diagnostics": [
            {
                "code": "CRN-E031",
                "path": "",
                "message": "unknown command or wrong argument count",
                "hint": USAGE_HINT,
                "severity": "error",
            }
        ],
    }, False


def _print(report):
    print(json.dumps(report, indent=2, allow_nan=False, sort_keys=True))


def _run_options(options):
    """``--archive DIR`` and ``--cancel-file PATH``, each once; None otherwise."""
    found = {}
    if len(options) % 2:
        return None
    for name, value in zip(options[::2], options[1::2], strict=True):
        if name not in ("--archive", "--cancel-file") or name in found:
            return None
        found[name] = value
    return found.get("--archive"), found.get("--cancel-file")


def _run(args):
    """run / rules: the native CLI's output text and exit status."""
    archive, cancel_file = _run_options(args[3:])
    text, success = _native.run_files_json(args[1], args[2], archive, cancel_file)
    sys.stdout.write(text + "\n")
    return success


def _convert(path):
    try:
        data = Path(path).read_bytes()
    except OSError as error:
        _print({"success": False, "message": str(error)})
        return False
    try:
        converted = _native.convert_network(data)
    except _native.NativeError as error:
        _print({"success": False, "diagnostics": json.loads(str(error))})
        return False
    sys.stdout.flush()
    sys.stdout.buffer.write(converted)
    sys.stdout.buffer.flush()
    return True


def _command(args):
    """Runs one visible command: (report, success), or None for wrong arguments."""
    name, count = args[0], len(args)
    if name == "describe" and count == 1:
        return json.loads(_native.describe_json()), True
    if name == "describe" and count == 2:
        try:
            return json.loads(_native.schema_json(args[1])), True
        except ValueError:
            return _usage_error()
    if name in ("check", "rules-check") and count in (2, 3):
        report = json.loads(
            _native.check_any_files_json(args[1], args[2] if count == 3 else None)
        )
        return report, report["valid"]
    if name == "analyze" and count in (2, 3):
        report = json.loads(
            _native.analyze_state_json(
                Path(args[1]).read_bytes(),
                Path(args[2]).read_bytes() if count == 3 else None,
            )
        )
        return report, report["success"]
    if name == "export" and count == 3 and args[2] == "gillespy2":
        report = json.loads(_native.export_json(Path(args[1]).read_bytes()))
        return report, report["success"]
    if name == "verify" and count == 2:
        report = json.loads(_native.verify_json(args[1]))
        return report, report["valid"]
    return None


def main():
    args = sys.argv[1:]
    if args in (["--help"], ["-h"]):
        sys.stdout.write(HELP)
        return
    if args == ["--version"]:
        _print({"version": __version__, "interface": "python", "milestone": MILESTONE})
        return
    try:
        if len(args) == 2 and args[0] == "convert":
            success = _convert(args[1])
            if not success:
                raise SystemExit(1)
            return
        if (
            args
            and args[0] in ("run", "rules")
            and len(args) >= 3
            and _run_options(args[3:]) is not None
        ):
            if not _run(args):
                raise SystemExit(1)
            return
        outcome = _command(args) if args else None
        report, success = _usage_error() if outcome is None else outcome
    except CrnError as error:
        report, success = {"success": False, "diagnostics": error.diagnostics}, False
    except (OSError, ValueError) as error:
        report, success = {"success": False, "message": str(error)}, False
    except KeyboardInterrupt:
        print(json.dumps({"success": False, "message": "interrupted"}))
        raise SystemExit(130) from None
    _print(report)
    if not success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
