"""Use one prepared native concentration model per evaluation worker."""

import json
from pathlib import Path
import subprocess
from queue import Empty, LifoQueue
from threading import Lock
from contextlib import contextmanager

import numpy as np

from common import write


def document(
    states, parameters, initial, derivatives, *, assignments=(), observations=()
):
    return dict(
        schema="crn-craft.concentration-model/0.1",
        states=states,
        parameters=parameters,
        initial=initial,
        derivatives=derivatives,
        assignments=assignments,
        observations=observations,
    )


class Solver:
    def __init__(self, binary, model_path, model):
        self.binary = Path(binary).resolve(strict=True)
        self.model_path = Path(model_path).resolve()
        write(self.model_path, model)
        self.available = LifoQueue()
        self.lock = Lock()
        self.workers = []

    def __enter__(self):
        return self

    def __exit__(self, *_):
        for process in self.workers:
            if process.stdin:
                process.stdin.close()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            process.stdout.close()

    @contextmanager
    def worker(self):
        try:
            process = self.available.get_nowait()
        except Empty:
            process = subprocess.Popen(
                [self.binary, "concentrations", self.model_path],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                text=True,
                encoding="utf-8",
            )
            with self.lock:
                self.workers.append(process)
        try:
            yield process
        finally:
            if process.poll() is None:
                self.available.put(process)

    def batch(self, parameters, times, *, rtol=2e-8, atol=1e-10):
        request = dict(
            schema="crn-craft.concentration-request/0.1",
            cases=[
                dict(parameters=list(p), times=list(times), rtol=rtol, atol=atol)
                for p in parameters
            ],
        )
        try:
            with self.worker() as process:
                process.stdin.write(json.dumps(request, allow_nan=False) + "\n")
                process.stdin.flush()
                line = process.stdout.readline()
        except BrokenPipeError as error:
            raise RuntimeError("concentration worker closed its input") from error
        if not line:
            raise RuntimeError("concentration worker returned no result")
        result = json.loads(line)
        if result.get("schema") != "crn-craft.concentration-result/0.1":
            raise ValueError(result.get("error", "unexpected concentration result"))
        if len(result["cases"]) != len(parameters):
            raise ValueError("incomplete concentration batch")
        state, readout = [], []
        for i, item in enumerate(result["cases"]):
            if "error" in item:
                raise ValueError(f"condition {i}: {item['error']}")
            state.append(item["values"]["states"])
            readout.append(item["values"]["observations"])
        return np.asarray(state), np.asarray(readout)
