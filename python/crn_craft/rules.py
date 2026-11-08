"""Prepared rules ensembles and their results as pandas and xarray objects.

``rules.ensemble(model, ...)`` builds and prepares a run request;
``simulate()`` returns a :class:`Result` whose records, statistics, rule
counts and trackers are pandas DataFrames, and whose kept paths and event
logs are xarray DataArrays (the 0.4.0 DataFrames are ``paths_frame()`` and
``event_frame()``). Every frame and array has ``.to_numpy()``, and
``Result.metadata`` keeps the full result document as the engine returned
it.

Opt-in, with a Rust toolchain: :func:`compile_loop` compiles a model's
Static Direct event loop, and ``compiled=`` runs it with byte-identical
results; :func:`loop_source` shows whether a model is supported.
"""

import hashlib
import json
import operator
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
import xarray

from . import CrnError, _crn_error, _native, _path, model_bytes
from ._rules_loop import CompiledLoop, fingerprint_of

# The compiled-loop API is part of this module.
from ._rules_loop import CompiledLoopError as CompiledLoopError
from ._rules_loop import compile_loop as compile_loop
from ._rules_loop import loop_source as loop_source
from .builder import Model

# Tolerances `ensemble` writes for method="hybrid" unless the method
# dictionary gives its own.
HYBRID_TOLERANCES = {"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}

TERMINATIONS = (
    "completed",
    "stopped",
    "event_limit",
    "sampling_limit",
    "instance_limit",
    "wall_time_limit",
    "cancelled",
    "failed",
    "output_limit",
)

# Record columns in DataFrame order, with the native buffers' dtypes.
_RECORD_DTYPES = {
    "started": "?",
    "termination": "u1",
    "stop_time": "<f8",
    "events": "<u8",
    "rng_words": "<u8",
    "peak_instances": "<u8",
    "created_instances": "<u8",
    "saved_points": "<u8",
}
_TRACKER_COUNTS = ("appearances", "extinctions", "omitted_episodes")
_TRACKER_VALUES = (
    "last_appearance",
    "last_extinction",
    "first_passage",
    "minimum",
    "maximum",
    "integral",
    "omitted_duration",
    "omitted_peak",
)
# Result metadata that Result.paths and Result.event_columns() carry as attrs.
_ARRAY_ATTRS = ("sampler", "method", "root_seed", "model_sha256", "run_sha256")
# The `rule` of an event-array entry past a replicate's logged length.
_PADDING_RULE = -2


def describe():
    """Return the native rules capabilities and limits."""
    return json.loads(_native.describe_rules_json())


def load(path):
    """Load and validate exact model bytes for reuse in multiple requests."""
    data = Path(path).read_bytes()
    check = json.loads(_native.check_rules_json(data))
    if not check["valid"]:
        from . import CrnError

        raise CrnError(check["diagnostics"])
    return data


def _replicate_index(ids):
    return pd.Index(np.asarray(ids, dtype=np.int64), name="replicate")


class Result:
    """The result of a rules ensemble.

    - ``records``: one row per replicate, indexed by replicate id. Columns
      ``started``, ``termination`` (categorical over :data:`TERMINATIONS`),
      ``stop_time``, ``events``, ``rng_words``, ``peak_instances``,
      ``created_instances`` and ``saved_points``.
    - ``rule_events``: events per replicate (rows) and rule (columns, in
      declaration order).
    - ``statistics``: None if withheld, else a dict with ``samples`` (the
      number of replicates included) and one DataFrame per statistic
      (``mean``, ``sample_variance``, ``mean_standard_error``), indexed by
      report time with one column per observable. A statistic that was not
      computed is None.
    - ``paths``: the kept series as an xarray DataArray with dimensions
      (``time``, ``variable``, ``replicate``); ``paths_frame()`` gives the
      DataFrame of release 0.4.0, ``to_xarray()`` a Dataset.
    - ``event_columns()``: event logs as a DataArray with dimensions
      (``replicate``, ``entry``, ``variable``); ``event_frame()`` gives the
      DataFrame of release 0.4.0.
    - ``tracker_columns(i)``: a tracker's fields as a DataFrame.
    - ``metadata``: the full result document; ``observations`` its trackers
      and logs; ``record_diagnostics`` replicate errors by id.

    Observables are in result order: sorted by name, then any continuous
    variable such as ``ATP``. Use ``.to_numpy()`` on any frame, column or
    array for NumPy arrays.
    """

    def __init__(self, metadata, columns):
        self.metadata = metadata["result"]
        self.record_diagnostics = dict(metadata["record_diagnostics"])
        self.observations = self.metadata.get("observations", [])
        self.rule_names = tuple(self.metadata["rules"])
        observables = list(self.metadata["observables"])
        index = _replicate_index(np.frombuffer(columns["id"], dtype="<u8"))
        data = {}
        for key, dtype in _RECORD_DTYPES.items():
            values = np.frombuffer(columns[key], dtype=dtype)
            if key == "termination":
                values = pd.Categorical.from_codes(values, categories=TERMINATIONS)
            data[key] = values
        self.records = pd.DataFrame(data, index=index)
        counts = np.frombuffer(columns["rule_events"], dtype="<u8")
        self.rule_events = pd.DataFrame(
            counts.reshape(len(index), len(self.rule_names)),
            index=index,
            columns=pd.Index(self.rule_names, name="rule"),
        )
        stats = self.metadata["statistics"]
        if stats is None:
            self.statistics = None
        else:
            times = pd.Index(self.metadata["report_times"], name="time")
            self.statistics = {}
            for key, value in stats.items():
                if isinstance(value, list):
                    matrix = np.array(value, dtype=np.float64)
                    value = pd.DataFrame(
                        matrix.reshape(len(times), len(observables)),
                        index=times,
                        columns=observables,
                    )
                self.statistics[key] = value

    @cached_property
    def _path_arrays(self):
        """Kept replicate ids, row counts, row positions and values."""
        kept = self.metadata["paths"]
        width = len(self.metadata["observables"])
        ids = np.fromiter((p["id"] for p in kept), dtype=np.int64, count=len(kept))
        lengths = np.fromiter(
            (len(p["values"]) for p in kept), dtype=np.intp, count=len(kept)
        )
        total = int(lengths.sum())
        values = np.array(
            [row for p in kept for row in p["values"]], dtype=np.float64
        ).reshape(total, width)
        # Position of each row in the report grid: rows are a prefix of it.
        positions = np.arange(total) - np.repeat(np.cumsum(lengths) - lengths, lengths)
        return ids, lengths, positions, values

    def _array_attrs(self):
        """Provenance for a DataArray: sampler, method, seed and hashes."""
        meta = self.metadata
        attrs = {key: meta.get(key) for key in _ARRAY_ATTRS}
        attrs["method"] = json.loads(json.dumps(attrs["method"]))
        return attrs

    @cached_property
    def paths(self):
        """Kept series as a DataArray: report times x variables x replicates.

        Dimensions are ``("time", "variable", "replicate")`` with the
        coordinates ``time`` (the report times), ``variable`` (the
        observables in result order: sorted names, then continuous
        variables) and ``replicate`` (the kept replicate ids, int64). Values
        are float64. ``attrs`` holds the result's ``sampler``, ``method``
        (the method object of rules-result/0.3, else None), ``root_seed``,
        ``model_sha256`` and ``run_sha256``.

        ``keep_paths`` decides which replicates are kept (``"all"`` keeps
        every one). A replicate that ended on a budget has values only up
        to its last saved report time and NaN after it. For example
        ``paths.sel(variable="mutants")`` is a time x replicate array and
        ``paths.sel(replicate=3)`` one replicate's series.
        ``paths_frame()`` returns the DataFrame of release 0.4.0.
        """
        ids, lengths, positions, values = self._path_arrays
        times = np.asarray(self.metadata["report_times"], dtype=np.float64)
        observables = list(self.metadata["observables"])
        cube = np.full((len(times), len(observables), len(ids)), np.nan)
        cube[positions, :, np.repeat(np.arange(len(ids)), lengths)] = values
        return xarray.DataArray(
            cube,
            dims=("time", "variable", "replicate"),
            coords={"time": times, "variable": observables, "replicate": ids},
            name="paths",
            attrs=self._array_attrs(),
        )

    def paths_frame(self):
        """Kept series as a DataFrame: (replicate, time) rows.

        One column per observable, in result order; this is ``paths`` of
        release 0.4.0. A replicate that ended on a budget has rows only up
        to its last saved report time.
        """
        ids, lengths, positions, values = self._path_arrays
        times = np.asarray(self.metadata["report_times"], dtype=np.float64)
        index = pd.MultiIndex.from_arrays(
            [np.repeat(ids, lengths), times[positions]], names=["replicate", "time"]
        )
        return pd.DataFrame(
            values, index=index, columns=list(self.metadata["observables"])
        )

    def to_xarray(self):
        """Return the kept paths as an xarray Dataset.

        Dimensions are ``replicate`` and ``time`` (the report times), with
        one variable per observable; rows a replicate did not save are NaN.
        It holds the values of ``paths``:
        ``paths.to_dataset("variable").transpose("replicate", "time")``.
        ``attrs`` holds ``sampler``, ``model_sha256``, ``run_sha256`` and
        ``root_seed``, all of which netCDF can store.
        """
        ids, lengths, positions, values = self._path_arrays
        times = np.asarray(self.metadata["report_times"], dtype=np.float64)
        observables = self.metadata["observables"]
        cube = np.full((len(ids), len(times), len(observables)), np.nan)
        cube[np.repeat(np.arange(len(ids)), lengths), positions] = values
        meta = self.metadata
        return xarray.Dataset(
            {
                name: (("replicate", "time"), cube[:, :, column])
                for column, name in enumerate(observables)
            },
            coords={"replicate": ids, "time": times},
            attrs={
                key: meta[key]
                for key in ("sampler", "model_sha256", "run_sha256", "root_seed")
            },
        )

    def _observed(self, replicate):
        """The observations of every replicate, or of one replicate id."""
        if replicate is None:
            return self.observations
        chosen = [o for o in self.observations if o["id"] == replicate][:1]
        if not chosen:
            raise KeyError(replicate)
        return chosen

    def event_columns(self, replicate=None):
        """Return retained event logs as a DataArray.

        Without a replicate id the dimensions are ``("replicate", "entry",
        "variable")``: one row per observed replicate (int64 ids), entries
        0, 1, ... in event order, and the observables after each event
        (result order), NaN past a replicate's last entry. Coordinates over
        (``replicate``, ``entry``): ``time`` (float64, NaN as padding),
        ``rule`` (the rule index, -1 for a leap, -2 as padding) and
        ``rule_name`` (the name; NaN for a leap and as padding); over
        ``replicate``: ``length``, the number of logged entries. With a
        replicate id the dimensions are ``("entry", "variable")``, exactly
        ``length`` entries, with scalar ``replicate`` and ``length``
        coordinates. ``attrs`` is that of ``paths``.

        A replicate without logged events has ``length`` 0. Participant
        instances and leap firings are in the observation entries.
        ``event_frame()`` returns the DataFrame of release 0.4.0.
        """
        chosen = self._observed(replicate)
        observables = list(self.metadata["observables"])
        lengths = np.fromiter(
            (len(o["log"]) for o in chosen), dtype=np.int64, count=len(chosen)
        )
        width = int(lengths.max(initial=0))
        rows = np.repeat(np.arange(len(chosen)), lengths)
        entries = np.arange(int(lengths.sum())) - np.repeat(
            np.cumsum(lengths) - lengths, lengths
        )
        log = [e for o in chosen for e in o["log"]]
        shape = (len(chosen), width)
        values = np.full((*shape, len(observables)), np.nan)
        values[rows, entries] = np.array(
            [e["observables"] for e in log], dtype=np.float64
        ).reshape(len(log), len(observables))
        time = np.full(shape, np.nan)
        time[rows, entries] = np.fromiter(
            (e["time"] for e in log), dtype=np.float64, count=len(log)
        )
        rule = np.full(shape, _PADDING_RULE, dtype=np.int64)
        rule[rows, entries] = np.fromiter(
            (-1 if e["rule"] is None else e["rule"] for e in log),
            dtype=np.int64,
            count=len(log),
        )
        # Index 0 is the missing name of leaps and padding; rule r is at r + 1.
        names = np.array([np.nan, *self.rule_names], dtype=object)
        rule_name = names[np.where(rule < 0, 0, rule + 1)]
        ids = np.fromiter((o["id"] for o in chosen), dtype=np.int64, count=len(chosen))
        array = xarray.DataArray(
            values,
            dims=("replicate", "entry", "variable"),
            coords={
                "replicate": ids,
                "entry": np.arange(width, dtype=np.int64),
                "variable": observables,
                "length": ("replicate", lengths),
                "time": (("replicate", "entry"), time),
                "rule": (("replicate", "entry"), rule),
                "rule_name": (("replicate", "entry"), rule_name),
            },
            name="event_log",
            attrs=self._array_attrs(),
        )
        return array if replicate is None else array.isel(replicate=0)

    def event_frame(self, replicate=None):
        """Return retained event logs as a DataFrame (release 0.4.0's form).

        With a replicate id, that replicate's log, indexed by ``entry``;
        without one, every logged event, indexed by (``replicate``,
        ``entry``). Rows are in event order: ``time``, ``rule`` (the rule
        index, -1 for a leap), ``rule_name`` (categorical, missing for a
        leap), then one column per observable after the event. Participant
        instances and leap firings are in the observation entries.
        """
        chosen = self._observed(replicate)
        entries = [e for o in chosen for e in o["log"]]
        size = len(entries)
        observables = list(self.metadata["observables"])
        rule = np.fromiter(
            (-1 if e["rule"] is None else e["rule"] for e in entries),
            dtype=np.int64,
            count=size,
        )
        if replicate is None:
            lengths = [len(o["log"]) for o in chosen]
            index = pd.MultiIndex.from_arrays(
                [
                    np.repeat([o["id"] for o in chosen], lengths).astype(np.int64),
                    np.concatenate([np.arange(n) for n in lengths] or [[]]).astype(
                        np.int64
                    ),
                ],
                names=["replicate", "entry"],
            )
        else:
            index = pd.RangeIndex(size, name="entry")
        frame = pd.DataFrame(
            {
                "time": np.fromiter(
                    (e["time"] for e in entries), dtype=np.float64, count=size
                ),
                "rule": rule,
                "rule_name": pd.Categorical.from_codes(
                    rule, categories=list(self.rule_names)
                ),
            },
            index=index,
        )
        values = np.array([e["observables"] for e in entries], dtype=np.float64)
        values = values.reshape(size, len(observables))
        return frame.join(pd.DataFrame(values, index=index, columns=observables))

    def tracker_columns(self, index=0):
        """Return a tracker's scalar fields, one row per observed replicate.

        Indexed by replicate id. Counts (``appearances``, ``extinctions``,
        ``omitted_episodes``) are unsigned integers; times and values are
        floats, NaN where absent (for example ``first_passage`` of a
        replicate that never met the predicate). Replicates that failed
        before tracking began carry no observation. Episodes are in
        ``observations``.
        """
        if index < 0:
            raise ValueError("tracker index must be nonnegative")
        count = max((len(ob["trackers"]) for ob in self.observations), default=0)
        if index >= count:
            raise IndexError(f"tracker index {index} with {count} trackers")
        rows = [ob["trackers"][index] for ob in self.observations]
        data = {}
        for key in _TRACKER_COUNTS:
            data[key] = np.array([t[key] for t in rows], dtype=np.uint64)
        for key in _TRACKER_VALUES:
            data[key] = np.array([t[key] for t in rows], dtype=np.float64)
        ids = [ob["id"] for ob in self.observations]
        return pd.DataFrame(data, index=_replicate_index(ids))


class PreparedEnsemble:
    def __init__(self, native, rule_names, documents=None, compiled=None):
        self._native = native
        self._rule_names = tuple(rule_names)
        # The exact (model, request) bytes, which a compiled loop runs.
        self._documents = documents
        self._fingerprint = None
        self._compiled = None
        if compiled is not None:
            self._compiled_request(compiled)
            self._compiled = compiled

    def _compiled_request(self, compiled):
        """Check that ``compiled`` runs this request; return the request."""
        if not isinstance(compiled, CompiledLoop):
            raise TypeError("compiled must be a CompiledLoop from rules.compile_loop")
        if self._fingerprint is None:
            try:
                source = self._native.loop_source()
            except _native.NativeError as error:
                raise _crn_error(error) from None
            self._fingerprint = fingerprint_of(source)
        if self._fingerprint != compiled.fingerprint or self._documents is None:
            raise CrnError(
                [
                    {
                        "code": "CRN-R020",
                        "path": "/model",
                        "message": "generated rules loop fingerprint mismatch",
                        "hint": "A compiled loop runs only the model bytes it was "
                        "compiled for; compile one for this model and parameter "
                        "set with rules.compile_loop.",
                    }
                ]
            )
        return self._documents[1]

    @property
    def rule_names(self):
        """Rule names in declaration order.

        Callback events and log entries carry a rule index into this tuple,
        so ``run.rule_names[event["rule"]]`` names the rule that fired.
        """
        return self._rule_names

    def bind_parameters(self, overrides):
        """Prepare an independent candidate with the same ensemble settings.

        The candidate's model and run hashes cover the re-serialized bound
        model and request, not the original file. It runs interpreted: a
        compiled loop serves only the model bytes it was compiled for.
        """
        try:
            native = self._native.bind_parameters(dict(overrides))
        except _native.NativeError as error:
            raise _crn_error(error) from None
        return PreparedEnsemble(native, self._rule_names)

    def simulate(
        self,
        *,
        cancel_file=None,
        callback=None,
        batch_callback=None,
        batch_size=256,
        compiled=None,
    ):
        """Run native workers.

        A live callback(id, event) returning False stops that replicate; any
        other return value continues. ``event`` is a log entry: ``time``,
        ``rule`` (an index into ``rule_names``, None for a leap),
        ``reactants``, ``products``, ``observables`` (in the result's column
        order) and ``firings``. Workers preserve per-replicate order, and
        calls from different replicates interleave. An exception in the
        callback cancels the run and is raised here; partial results are
        discarded. A callback that blocks also blocks Ctrl-C. A
        batch_callback(id, events) replays the retained log after execution
        and cannot affect the run.

        ``compiled``, a :class:`CompiledLoop` for this model, runs the
        request in the compiled program instead of the interpreter, with
        byte-identical results; it defaults to the loop given to
        ``ensemble`` or ``prepare``. The request must select
        ``static_direct`` (``CrnError`` ``CRN-R020`` at ``/method``
        otherwise) and use this loop's model bytes (``CRN-R020`` at
        ``/model``). Callbacks need the interpreter and are refused with
        ``ValueError``.
        """
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        compiled = self._compiled if compiled is None else compiled
        if compiled is not None:
            if callback is not None or batch_callback is not None:
                raise ValueError(
                    "callback and batch_callback need the interpreter; a compiled "
                    "loop runs in its own process. Use trackers or an event log, "
                    "or simulate without compiled."
                )
            request = self._compiled_request(compiled)
            metadata, columns = compiled.execute(request, cancel_file=cancel_file)
            return Result(json.loads(metadata), columns)

        def dispatch(replicate, event):
            return callback(replicate, json.loads(event))

        try:
            metadata, columns = self._native.simulate(
                _path(cancel_file), None if callback is None else dispatch
            )
        except _native.NativeError as error:
            raise _crn_error(error) from None
        result = Result(json.loads(metadata), columns)
        if batch_callback is not None:
            for observation in result.observations:
                log = observation["log"]
                for start in range(0, len(log), batch_size):
                    batch_callback(observation["id"], log[start : start + batch_size])
        return result


def prepare(model, request, *, compiled=None):
    """Prepare a rules-run contract, preserving the supplied model bytes.

    ``model`` is bytes, text, a dictionary or a ``Model``. ``compiled`` is
    an optional :class:`CompiledLoop` for this model that ``simulate`` then
    uses (see ``PreparedEnsemble.simulate``).
    """
    if isinstance(model, Model):
        model = model.to_dict()
    data = model_bytes(model)
    request = model_bytes(request)
    try:
        native = _native.prepare_rules(data, request)
    except _native.NativeError as error:
        raise _crn_error(error) from None
    # The engine has validated the model; results list rules in this order.
    names = [r["name"] for r in json.loads(data)["rules"]]
    return PreparedEnsemble(native, names, (data, request), compiled)


def ensemble(
    model,
    *,
    end,
    report_times,
    replicates=1,
    seed=0,
    threads=1,
    keep_paths=0,
    keep_instances=0,
    max_events=1_000_000,
    max_sampling_steps=8_000_000,
    max_instances=100_000,
    max_output_bytes=64 << 20,
    wall_time_s=None,
    trackers=(),
    event_log=None,
    method=None,
    compiled=None,
):
    """Build and prepare a rules request with explicit work/output budgets.

    ``model`` is model bytes (from :func:`load`), a JSON string, a
    dictionary or a :class:`Model`. ``report_times`` are increasing times in [0, ``end``].
    ``keep_paths`` is the number of replicates, from replicate 0, whose
    observable series are kept as ``Result.paths``; ``"all"`` or ``True``
    keeps every replicate. ``keep_instances`` keeps final instances the same
    way (a count). The budgets ``max_events``, ``max_sampling_steps`` and
    ``max_instances`` are per replicate; ``wall_time_s`` covers the run and
    ``max_output_bytes`` the result. ``trackers`` and ``event_log`` select
    the observation. The lowest request version that fits the arguments is
    used; ``"static_direct"`` selects rules-run/0.6 for its defaults.

    ``method`` is a name or a dictionary such as ``{"kind": "tau",
    "epsilon": 0.03}``. The generic methods run every model in their scope:

    - ``"static_direct"``: exact, for closed models (a fixed set of
      attribute vectors). It selects rules-run/0.6, whose defaults are
      Ziggurat waiting times and sparse integration of every eligible
      tracker (an unweighted single-class count without a first-passage
      predicate). The dictionary options ``"waiting_time"`` (``"inverse"``
      or ``"ziggurat"``) and ``"sparse_trackers"`` (``False`` or ``True``)
      override them; ``{"kind": "static_direct", "waiting_time":
      "inverse", "sparse_trackers": False}`` reproduces the
      ``"static_direct"`` numbers of 0.3.0;
    - ``"next_reaction"``: exact, for closed models;
    - ``None``: the exact tree method, for every model without an ATP
      block, including models whose instances change identity;
    - ``"tau"`` and ``"binomial"``: approximate, for closed models;
    - ``"hybrid"``: for models with continuous variables, integrating the
      hazards to a declared tolerance; ``"hybrid"`` or a hybrid dictionary
      without tolerances gets ``rtol=1e-6`` and ``atol=1e-9`` and selects
      rules-run/0.5.

    Two methods serve one structure each. ``"thinning"`` is required by,
    and only for, models with an ATP block. ``"mtdna_direct"`` is an opt-in
    exact fast path for the three-class fixed-population structure: three
    classes in eight unit channels, one shared constant degradation rate,
    replication and mutation rates ``max(0, K - N) / D`` times constants,
    observables that are weighted sums, and no stop or
    ``N > 0 && mutants >= f * N``. It is certified per model and parameter
    set at preparation, samples the same process with a different use of
    the random stream (sampler ``chacha8-u53-ziggurat/rules-mtdna-direct/1``)
    and is about 3.3x faster than the Static Direct default on that
    structure (3.9x faster than Static Direct with inverse waits). A model
    it cannot certify runs Static Direct with inverse waits and per-event
    trackers, with the numbers of ``{"kind": "static_direct",
    "waiting_time": "inverse", "sparse_trackers": False}`` and a
    ``CRN-R051`` warning naming the reason.

    ``compiled`` is an optional :class:`CompiledLoop` from
    :func:`compile_loop` for this model; ``simulate`` then runs the request
    in the compiled program, with results byte-identical to the
    interpreter's. It needs ``method="static_direct"`` (with any options).
    """
    if isinstance(model, Model):
        model = model.to_dict()
    data = model_bytes(model)
    # NumPy scalars and arrays are converted; counts must be integral.
    count = operator.index
    if isinstance(keep_paths, bool | np.bool_):
        keep_paths = count(replicates) if keep_paths else 0
    elif isinstance(keep_paths, str):
        if keep_paths != "all":
            raise ValueError(
                f"keep_paths must be a count, 'all' or True: {keep_paths!r}"
            )
        keep_paths = count(replicates)
    request = {
        "schema": "crn-craft.rules-run/0.1",
        "model": {
            "path": "model.rules.json",
            "sha256": hashlib.sha256(data).hexdigest(),
        },
        "time": {"end": float(end), "report_times": [float(t) for t in report_times]},
        "ensemble": {"replicates": count(replicates), "root_seed": count(seed)},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": count(keep_paths),
            "keep_instances": count(keep_instances),
        },
        "budget": {
            "threads": count(threads),
            "max_events_per_replicate": count(max_events),
            "max_sampling_steps_per_replicate": count(max_sampling_steps),
            "max_instances": count(max_instances),
            "max_output_bytes": count(max_output_bytes),
            "wall_time_s": None if wall_time_s is None else float(wall_time_s),
        },
    }
    if trackers or event_log is not None or method is not None:
        # Present fields must hold values; omit the ones that were not given.
        request["schema"] = "crn-craft.rules-run/0.2"
        if trackers:
            request["trackers"] = list(trackers)
        if event_log is not None:
            request["event_log"] = event_log
        if method is not None:
            request["method"] = {"kind": method} if isinstance(method, str) else method
        if (
            isinstance(request.get("method"), dict)
            and request["method"].get("kind") == "hybrid"
        ):
            # The request declares the tolerances; these are the defaults.
            request["method"] = {**HYBRID_TOLERANCES, **request["method"]}
            request["schema"] = "crn-craft.rules-run/0.5"
        elif (
            isinstance(request.get("method"), dict)
            and request["method"].get("kind") == "mtdna_direct"
        ):
            request["schema"] = "crn-craft.rules-run/0.3"
        elif (
            isinstance(request.get("method"), dict)
            and request["method"].get("kind") == "static_direct"
        ):
            # Options left out take the rules-run/0.6 defaults.
            request["schema"] = "crn-craft.rules-run/0.6"
    return prepare(data, request, compiled=compiled)
