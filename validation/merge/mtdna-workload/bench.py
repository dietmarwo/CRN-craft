"""Time crn_craft.rules on the mtDNA manual workload in the running interpreter.

usage: python bench.py MODEL METHOD THREADS REPLICATES OUT_JSON
Prints one JSON line: version, method, threads, seconds (one run), digest of
the result (statistics, records and termination), sampler.
"""

import hashlib
import json
import sys
import time

import crn_craft
from crn_craft import rules

GRID = [float(t) for t in range(0, 25921, 720)]

model_path, method, threads, replicates, out = sys.argv[1:6]
model = rules.load(model_path)
method = None if method == "tree" else method
start = time.perf_counter()
run = rules.ensemble(
    model,
    end=GRID[-1],
    report_times=GRID,
    replicates=int(replicates),
    seed=1,
    threads=int(threads),
    method=method,
)
result = run.simulate()
seconds = time.perf_counter() - start
digest = hashlib.sha256()
for name in ("mean", "sample_variance"):
    digest.update(result.statistics[name].to_csv().encode())
digest.update(result.records.to_csv().encode())
row = {
    "version": crn_craft.__version__,
    "python": sys.executable,
    "method": method or "tree",
    "threads": int(threads),
    "replicates": int(replicates),
    "seconds": seconds,
    "digest": digest.hexdigest(),
    "sampler": result.metadata["sampler"],
    "run_sha256": result.metadata.get("run_sha256"),
}
with open(out, "a") as handle:
    handle.write(json.dumps(row) + "\n")
print(json.dumps(row))
