"""Development timing/replay smoke of the Static Direct selectors; not the final comparison."""

import argparse
import json
import time
from pathlib import Path

import crn_craft as crn
from benchmark_campaign import models

SELECTORS = ["linear", "composition_rejection", "sorting", "rssa_cr"]


def prepared(binary_model, selector, wait):
    return crn.ensemble(
        binary_model,
        end=100.0,
        report_times=[0.0, 100.0],
        replicates=16,
        seed=4242,
        threads=1,
        max_events=10**9,
        max_sampling_steps=10**10,
        max_output_bytes=200_000_000,
        method={"kind": "static_direct", "waiting_time": wait},
        selector=selector,
        quantiles=(),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = []
    for name, model in models():
        converted = crn.convert((json.dumps(model) + "\n").encode())
        for selector in SELECTORS:
            for wait in ["inverse", "ziggurat"]:
                job = prepared(converted, selector, wait)
                start = time.perf_counter()
                result = job.simulate()
                seconds = time.perf_counter() - start
                meta = result.metadata
                if meta["termination"] != "completed":
                    raise RuntimeError(meta)
                again = prepared(converted, selector, wait).simulate().metadata
                if (
                    again["events"] != meta["events"]
                    or again["run_sha256"] != meta["run_sha256"]
                ):
                    raise RuntimeError("prepared equality")
                rows.append(
                    dict(
                        workload=name,
                        selector=selector,
                        waiting_time=wait,
                        seconds=seconds,
                        events=int(meta["events"]),
                        prepared_equal=True,
                    )
                )
                args.output.write_text(json.dumps(rows, indent=2) + "\n")
                print(name, selector, wait, seconds, flush=True)


if __name__ == "__main__":
    main()
