"""Installed native-backed CLI. Filesystem preflight and archives stay native."""

import argparse
import json
from pathlib import Path
import sys
from . import CrnError, __version__, _call, _native, describe, fsp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version",
        action="version",
        version=json.dumps({"version": __version__, "interface": "python"}),
    )
    sub = parser.add_subparsers(dest="command", required=True)
    describe_parser = sub.add_parser("describe")
    describe_parser.add_argument("schema", nargs="?")
    for command in ["check", "analyze", "analyze-v1", "fsp", "tau", "ssa", "export"]:
        p = sub.add_parser(command)
        p.add_argument("model", type=Path)
        if command == "check":
            p.add_argument("run", type=Path, nargs="?")
        if command in ["analyze", "analyze-v1"]:
            p.add_argument("state", type=Path, nargs="?")
        if command in ["fsp", "tau", "ssa"]:
            p.add_argument("run", type=Path)
            p.add_argument("--cancel-file")
        if command == "export":
            p.add_argument("target", choices=["gillespy2"])
    for command in ["simulate", "simulate-cached"]:
        p = sub.add_parser(command)
        p.add_argument("model", type=Path)
        p.add_argument("run", type=Path)
        p.add_argument("output", type=Path)
        p.add_argument("--cancel-file")
    p = sub.add_parser("verify")
    p.add_argument("directory", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "describe":
            report = (
                json.loads(_native.schema_json(args.schema))
                if args.schema
                else describe()
            )
            success = True
        elif args.command == "check":
            report = json.loads(
                _native.check_files_json(
                    str(args.model), str(args.run) if args.run else None
                )
            )
            success = report["valid"]
        elif args.command in ["analyze", "analyze-v1"]:
            report = json.loads(
                _native.analyze_state_json(
                    args.model.read_bytes(),
                    args.state.read_bytes() if args.state else None,
                    args.command == "analyze-v1",
                )
            )
            success = report["success"]
        elif args.command == "export":
            report = json.loads(_native.export_json(args.model.read_bytes()))
            success = report["success"]
        elif args.command == "ssa":
            report = _call(
                _native.ssa_files_json, str(args.model), str(args.run), args.cancel_file
            )
            success = report.get("termination") == "completed"
        elif args.command == "tau":
            report = json.loads(
                _native.tau_files_json(str(args.model), str(args.run), args.cancel_file)
            )
            success = report.get("termination") == "completed"
        elif args.command == "fsp":
            check = json.loads(_native.check_files_json(str(args.model), str(args.run)))
            if not check["valid"]:
                report = check
                success = False
            else:
                report = fsp(
                    args.model.read_bytes(),
                    args.run.read_bytes(),
                    cancel_file=args.cancel_file,
                )
                success = True
        elif args.command == "verify":
            report = json.loads(_native.verify_json(str(args.directory)))
            success = report["valid"]
        else:
            report = json.loads(
                _native.archive_json(
                    str(args.model),
                    str(args.run),
                    str(args.output),
                    args.command == "simulate-cached",
                    args.cancel_file,
                )
            )
            success = report["valid"] and report["simulation_complete"]
    except CrnError as error:
        report = {"success": False, "diagnostics": error.diagnostics}
        success = False
    except (OSError, ValueError) as error:
        report = {"success": False, "message": str(error)}
        success = False
    except KeyboardInterrupt:
        print(json.dumps({"success": False, "message": "interrupted"}))
        raise SystemExit(130) from None
    if args.command == "ssa" or (args.command in ["fsp", "tau"] and success):
        # Checked against max_output_bytes in compact form, like the native CLI.
        sys.stdout.write(
            json.dumps(report, allow_nan=False, separators=(",", ":"), sort_keys=True)
        )
        sys.stdout.write("\n")
        if not success:
            raise SystemExit(1)
        return
    print(json.dumps(report, indent=2, allow_nan=False, sort_keys=True))
    if not success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
