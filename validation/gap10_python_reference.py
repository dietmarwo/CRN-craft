"""Isolated Python wheel reference process; IPC is outside its prepared timer."""

import hashlib
import json
from pathlib import Path
import sys
import time
import crn_craft as crn
from compare_julia_accelerated import RUST, rust_job

jobs = {}
for line in sys.stdin:
    request = json.loads(line)
    if request["command"] == "prepare":
        model = Path(request["model_path"]).read_bytes()
        jobs = {
            k: rust_job(
                k, model, request["end"], request["paths"], 4242, False, "prepared"
            )
            for k in RUST
        }
        result = dict(
            version=crn.__version__,
            module=crn.__file__,
            native_sha256=hashlib.sha256(
                Path(crn._native.__file__).read_bytes()
            ).hexdigest(),
        )
    elif request["command"] == "run":
        start = time.perf_counter()
        value = jobs[request["engine"]]()
        elapsed = time.perf_counter() - start
        result = dict(seconds=elapsed, value=value)
    else:
        raise RuntimeError("unknown command")
    print(json.dumps(result), flush=True)
