"""Client of the persistent `crn-diffsol concentrations` process.

One process per evaluation thread prepares the model once and answers one
request line with one result line (`integrations/diffsol/CONCENTRATIONS.md`).
A request holds 1 to 32 cases. A failed case is returned as data with a
diagnostic code; `Solver.batch` raises `CaseError` for it, so that the
search records the candidate as failed instead of scoring it.
"""

from contextlib import contextmanager
import json
from pathlib import Path
from queue import Empty, LifoQueue
import selectors
import subprocess
import tempfile
from threading import Lock

import numpy as np

from common import write

RESULT = "crn-craft.concentration-result/0.1"
MAX_CASES = 32


class CaseError(ValueError):
    """A case of a batch failed; `failures` holds (index, code, message)."""

    def __init__(self, failures):
        self.failures = failures
        index, code, message = failures[0]
        super().__init__(f"condition {index}: {code}: {message}")


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


def _failure(error):
    # A worker built before the coded failures reports a bare message.
    if isinstance(error, str):
        return "CRN-E080", error
    return error["code"], error["message"]


class Solver:
    """Pool of worker processes for one model; use it as a context manager.

    `timeout` bounds the wait for one result line in seconds. A worker that
    exceeds it, closes its output or answers out of protocol is killed and
    replaced on the next call.
    """

    def __init__(self, binary, model_path, model, *, timeout=600.0):
        self.binary = Path(binary).resolve(strict=True)
        self.model_path = Path(model_path).resolve()
        write(self.model_path, model)
        self.timeout = timeout
        self.available = LifoQueue()
        self.lock = Lock()
        self.workers = []

    def __enter__(self):
        return self

    def __exit__(self, *_):
        for process, errors in self.workers:
            if process.stdin and not process.stdin.closed:
                try:
                    process.stdin.close()
                except BrokenPipeError:
                    pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            process.stdout.close()
            errors.close()

    def _start(self):
        # The diagnostic of a refused model goes to stderr; a file instead of
        # a pipe keeps a talkative worker from blocking on a full buffer.
        errors = tempfile.TemporaryFile()
        process = subprocess.Popen(
            [self.binary, "concentrations", self.model_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=errors,
            text=True,
            encoding="utf-8",
        )
        with self.lock:
            self.workers.append((process, errors))
        return process, errors

    @contextmanager
    def worker(self):
        try:
            process, errors = self.available.get_nowait()
        except Empty:
            process, errors = self._start()
        try:
            yield process, errors
        finally:
            if process.poll() is None:
                self.available.put((process, errors))

    def _exchange(self, request):
        with self.worker() as (process, errors):
            try:
                process.stdin.write(request + "\n")
                process.stdin.flush()
            except BrokenPipeError:
                line = ""
            else:
                with selectors.DefaultSelector() as selector:
                    selector.register(process.stdout, selectors.EVENT_READ)
                    ready = selector.select(self.timeout)
                if not ready:
                    process.kill()
                    process.wait()
                    raise RuntimeError(
                        f"concentration worker gave no result in {self.timeout} s"
                    )
                line = process.stdout.readline()
            if not line:
                process.wait()
                errors.seek(0)
                detail = errors.read().decode("utf-8", "replace").strip()
                raise RuntimeError(
                    f"concentration worker exited with {process.returncode}: {detail}"
                )
            return line

    def batch(self, parameters, times, *, rtol=2e-8, atol=1e-10):
        """States and readouts, shape (cases, times, n), for each parameter vector."""
        parameters = [list(map(float, p)) for p in parameters]
        if not 1 <= len(parameters) <= MAX_CASES:
            raise ValueError(f"a batch holds 1 to {MAX_CASES} cases")
        request = dict(
            schema="crn-craft.concentration-request/0.1",
            cases=[
                dict(parameters=p, times=list(map(float, times)), rtol=rtol, atol=atol)
                for p in parameters
            ],
        )
        # allow_nan=False refuses NaN and infinity before they reach the worker.
        result = json.loads(self._exchange(json.dumps(request, allow_nan=False)))
        if result.get("schema") != RESULT:
            raise RuntimeError(f"unexpected concentration result: {result}")
        if "error" in result:
            raise ValueError("{}: {}".format(*_failure(result["error"])))
        if len(result["cases"]) != len(parameters):
            raise RuntimeError("incomplete concentration batch")
        failures = [
            (i, *_failure(item["error"]))
            for i, item in enumerate(result["cases"])
            if "error" in item
        ]
        if failures:
            raise CaseError(failures)
        state = [item["values"]["states"] for item in result["cases"]]
        readout = [item["values"]["observations"] for item in result["cases"]]
        return np.asarray(state), np.asarray(readout)
