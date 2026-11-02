"""Development timing/replay smoke; not the final interleaved comparison."""

import argparse
import json
import time
from pathlib import Path
import crn_craft as crn
from benchmark_campaign import models


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--rssacr", action="store_true")
    args = p.parse_args()
    rows = []
    algorithms = ["direct", "cached_direct", "direct_cr", "sorting_direct"]
    if args.rssacr:
        algorithms.append("rssa_cr")
    for name, model in models():
        for method in algorithms:
            for wait in ["inverse", "ziggurat"]:
                request = crn.ssa_spec(
                    model,
                    aggregator=method,
                    waiting_time=wait,
                    end=100.0,
                    save_times=[0.0, 100.0],
                    replicates=16,
                    threads=1,
                    max_output_bytes=200_000_000,
                )
                prepared = crn.prepare_ssa(model, request)
                start = time.perf_counter()
                result = prepared.simulate()
                seconds = time.perf_counter() - start
                if result["termination"] != "completed":
                    raise RuntimeError(result)
                if result != crn.simulate_ssa(model, request):
                    raise RuntimeError("prepared equality")
                row = dict(
                    workload=name,
                    algorithm=method,
                    waiting_time=wait,
                    seconds=seconds,
                    events=result["events"],
                    work=result["work"],
                    prepared_equal=True,
                )
                rows.append(row)
                args.output.write_text(json.dumps(rows, indent=2) + "\n")
                print(name, method, wait, seconds, flush=True)


if __name__ == "__main__":
    main()
