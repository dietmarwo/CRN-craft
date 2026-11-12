"""Write a model and a run request for the command line, and run them.

Usage: python examples/mtdna/11_cli_request.py [DIRECTORY [CRN_CRAFT_BINARY]]

DIRECTORY defaults to `takeover` in the temporary directory, /tmp/takeover
on Linux. With the native CLI binary, the script also runs the request.
"""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

if len(sys.argv) > 1:
    folder = Path(sys.argv[1])
else:
    folder = Path(tempfile.gettempdir()) / "takeover"
folder.mkdir(parents=True, exist_ok=True)
model = json.loads(Path("models/rules/mtdna-three-class.rules.json").read_text())
model["stop"] = "N > 0 && mutants >= 0.9 * N"
data = (json.dumps(model, indent=2) + "\n").encode()
(folder / "takeover.rules.json").write_bytes(data)
request = {
    "schema": "crn-craft.run/1.0",
    # The request binds the model by relative path and exact SHA-256.
    "model": {
        "path": "takeover.rules.json",
        "sha256": hashlib.sha256(data).hexdigest(),
    },
    "time": {"end": 25920.0, "report_times": [0.0, 8640.0, 17280.0, 25920.0]},
    "ensemble": {"replicates": 1000, "root_seed": 11},
    "outputs": {
        "statistics": ["mean", "variance"],
        "quantiles": [0.5],
        "keep_paths": 0,
        "keep_instances": 0,
    },
    "budget": {
        "threads": 2,
        "max_events_per_replicate": 10_000_000,
        "max_sampling_steps_per_replicate": 100_000_000,
        "max_instances": 10,
        "max_output_bytes": 64 << 20,
    },
    "method": {"kind": "mtdna_direct"},
    "trackers": [{"observable": "advantaged_mutants"}],
}
(folder / "takeover.run.json").write_text(json.dumps(request, indent=2) + "\n")
print("model sha256:", request["model"]["sha256"])
if len(sys.argv) > 2:
    done = subprocess.run(
        [sys.argv[2], "run", "takeover.rules.json", "takeover.run.json"],
        cwd=folder,
        capture_output=True,
        check=True,
    )
    result = json.loads(done.stdout)
    print(result["termination"], result["sampler"])
    print(
        "cells:",
        result["completed_replicates"],
        "stopped:",
        result["stopped_replicates"],
    )
    print("events:", result["events"])
    print("result fields:", sorted(result))
