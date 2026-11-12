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

from collections.abc import Iterable, Mapping
import copy
import hashlib
import json
import numbers
import operator
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
import xarray

from . import _native
from . import report as _report
from ._common import CrnError, _crn_error, _path, convert, model_bytes
from ._rules_loop import CompiledLoop, fingerprint_of

# The compiled-loop API is part of this module.
from ._rules_loop import CompiledLoopError as CompiledLoopError
from ._rules_loop import compile_loop as compile_loop
from ._rules_loop import loop_source as loop_source
from .builder import Model
from .report import load_report as load_report

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
# Quantile statistics and selectors need the 1.0 contracts (see ensemble).
MODEL_1_0 = "crn-craft.model/1.0"
RUN_1_0 = "crn-craft.run/1.0"
# What `fsp` requests when max_states or tolerance is not given: the engine's
# state cap and a total series truncation budget.
FSP_MAX_STATES = 100_000
FSP_TOLERANCE = 1e-12
_RULES_MODEL = "crn-craft.rules-model/"
_RULES_RUN = "crn-craft.rules-run/"
_NETWORK_MODEL = "crn-craft.model/0.1"
# The quantiles the top-level API (crn_craft.ensemble, crn_craft.prepare)
# requests unless told otherwise; crn_craft.rules requests none.
DEFAULT_QUANTILES = (0.025, 0.5, 0.975)
# Static Direct selectors of crn-craft.run/1.0 (method.selector).
SELECTORS = ("linear", "composition_rejection", "sorting", "rssa_cr", "auto")
# `quantiles=` of the top-level prepare: the default quantiles unless the
# request states its own.
_DEFAULT = object()
# Request versions whose Static Direct defaults differ from run/1.0's: inverse
# waits and per-event trackers unless the method says otherwise.
_BEFORE_ZIGGURAT_DEFAULTS = (
    "crn-craft.rules-run/0.2",
    "crn-craft.rules-run/0.3",
    "crn-craft.rules-run/0.4",
    "crn-craft.rules-run/0.5",
)


def _probabilities(quantiles):
    """The probabilities of ``quantiles=`` as floats; an empty list for none."""
    if quantiles is None:
        return []
    if isinstance(quantiles, str | bytes | Mapping) or not isinstance(
        quantiles, Iterable
    ):
        raise TypeError(
            "quantiles must be a sequence of probabilities, for example "
            f"(0.025, 0.5, 0.975): {quantiles!r}"
        )
    values = list(quantiles)
    for p in values:
        if isinstance(p, bool | np.bool_) or not isinstance(p, numbers.Real):
            raise TypeError(f"quantile probabilities must be real numbers: {p!r}")
    return [float(p) for p in values]


def _as_model_1_0(model, *, network=False):
    """The model's bytes as ``crn-craft.model/1.0``, only the schema changed.

    A rules-model document gets the schema string ``crn-craft.model/1.0``,
    which means what rules-model/0.4 means; every other byte is kept. A
    dictionary or ``Model`` is serialized as without quantiles, with the
    ``schema`` value replaced in place. With ``network=True`` (the top-level
    API) a network ``crn-craft.model/0.1`` document is converted
    (:func:`crn_craft.convert`). A model/1.0 document, and anything else
    (the engine refuses it), is returned unchanged.
    """
    if isinstance(model, Model):
        model = model.to_dict()
    if isinstance(model, Mapping):
        schema = model.get("schema")
        if isinstance(schema, str) and schema.startswith(_RULES_MODEL):
            model = {**model, "schema": MODEL_1_0}
        elif network and schema == _NETWORK_MODEL:
            return convert(model)
        return model_bytes(model)
    data = model_bytes(model)
    try:
        schema = json.loads(data).get("schema")
    except (ValueError, AttributeError):
        return data
    if network and schema == _NETWORK_MODEL:
        return convert(data)
    if not (isinstance(schema, str) and schema.startswith(_RULES_MODEL)):
        return data
    old = json.dumps(schema).encode()
    if data.count(old) != 1:
        raise ValueError(
            f"cannot rewrite the model's schema {schema!r} to {MODEL_1_0!r} in "
            f"place: the string occurs {data.count(old)} times in the bytes; "
            "pass the model as a dictionary"
        )
    return data.replace(old, json.dumps(MODEL_1_0).encode())


def _unique_keys(pairs):
    keys = [key for key, _ in pairs]
    if len(set(keys)) != len(keys):
        duplicate = next(k for k in keys if keys.count(k) > 1)
        raise CrnError(
            [
                {
                    "code": "CRN-E002",
                    "path": "",
                    "message": f"duplicate JSON key {duplicate!r} in the request",
                    "hint": "Remove the duplicate key.",
                }
            ]
        )
    return dict(pairs)


def _selector(selector):
    """``selector=`` checked: None or one of :data:`SELECTORS`."""
    if selector is not None and selector not in SELECTORS:
        raise ValueError(
            f"selector must be one of {', '.join(map(repr, SELECTORS))}: {selector!r}"
        )
    return selector


def _with_selector(method, selector):
    """A Static Direct method dictionary with ``selector`` added."""
    if not (isinstance(method, Mapping) and method.get("kind") == "static_direct"):
        raise ValueError(
            "selector= chooses the event selector of method='static_direct'; "
            f"the request's method is {method!r}"
        )
    if "selector" in method:
        raise ValueError(
            "the method already has a selector; give selector= or the field, not both"
        )
    return {**method, "selector": selector}


def _request_1_0(request, probabilities, selector, old_sha256, new_sha256):
    """A request as ``crn-craft.run/1.0``, with quantiles and a selector.

    A rules-run request becomes run/1.0 and keeps its meaning: a Static
    Direct method gets that version's defaults written out (rules-run/0.2 to 0.5:
    inverse waits and per-event trackers; every rules-run version: the
    linear selector, unless ``selector`` is given), and a model hash of the
    original model bytes is replaced by the hash of the model/1.0 bytes.
    Another version is left as it is, so the engine judges it.
    ``probabilities`` are added as ``outputs.quantiles`` (``ValueError`` if
    the request has the field already).
    """
    if isinstance(request, bytes | str):
        request = json.loads(request, object_pairs_hook=_unique_keys)
    else:
        request = copy.deepcopy(request)
    if not isinstance(request, dict):
        raise TypeError("the request must be a JSON object")
    outputs = request.get("outputs")
    if probabilities and isinstance(outputs, dict):
        if "quantiles" in outputs:
            raise ValueError(
                "the request already has outputs.quantiles; give quantiles= or "
                "the field, not both"
            )
        outputs["quantiles"] = probabilities
    method = request.get("method")
    if selector is not None:
        method = request["method"] = _with_selector(method, selector)
    schema = request.get("schema")
    if (
        isinstance(schema, str)
        and schema.startswith(_RULES_RUN)
        and isinstance(method, dict)
        and method.get("kind") == "static_direct"
    ):
        if schema in _BEFORE_ZIGGURAT_DEFAULTS:
            method.setdefault("waiting_time", "inverse")
            method.setdefault("sparse_trackers", False)
        # rules-run requests always scan linearly; run/1.0 defaults to auto.
        method.setdefault("selector", "linear")
    if isinstance(schema, str) and schema.startswith(_RULES_RUN):
        request["schema"] = RUN_1_0
    # Any other version (run/1.0 itself, or one the engine refuses) is kept.
    bound = request.get("model")
    if isinstance(bound, dict) and bound.get("sha256") == old_sha256:
        bound["sha256"] = new_sha256
    return request


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
      computed is None. A run with ``quantiles=`` also has ``quantiles``, a
      dict from each requested probability to a DataFrame of that shape,
      and, when 0.5 was requested, ``median``, the DataFrame of
      ``quantiles[0.5]`` (see :func:`ensemble`).
    - ``paths``: the kept series as an xarray DataArray with dimensions
      (``replicate``, ``time``, ``variable``); ``paths_frame()`` gives the
      DataFrame of release 0.4.0, ``to_xarray()`` a Dataset.
    - ``event_columns()``: event logs as a DataArray with dimensions
      (``replicate``, ``entry``, ``variable``); ``event_frame()`` gives the
      DataFrame of release 0.4.0.
    - ``tracker_columns(i)``: a tracker's fields as a DataFrame.
    - ``metadata``: the full result document; ``observations`` its trackers
      and logs; ``record_diagnostics`` replicate errors by id.
    - ``model`` and ``request``: the model and request documents that ran,
      and ``save(path)``, which writes them with the results to one
      markdown or JSON report.

    Observables are in result order: sorted by name, then any continuous
    variable such as ``ATP``. Use ``.to_numpy()`` on any frame, column or
    array for NumPy arrays.
    """

    # The documents that ran: (model bytes, request bytes, bound), where
    # `bound` marks a bind_parameters candidate (see PreparedEnsemble). None
    # for a Result built from a metadata document and columns.
    _ran = None

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

            def table(rows):
                matrix = np.array(rows, dtype=np.float64)
                return pd.DataFrame(
                    matrix.reshape(len(times), len(observables)),
                    index=times,
                    columns=observables,
                )

            self.statistics = {}
            for key, value in stats.items():
                if key == "quantiles":
                    # crn-craft.result/1.0: [{"p": ..., "values": ...}] in
                    # request order (increasing p).
                    value = {q["p"]: table(q["values"]) for q in value}
                elif isinstance(value, list):
                    value = table(value)
                self.statistics[key] = value
            if 0.5 in self.statistics.get("quantiles", {}):
                self.statistics["median"] = self.statistics["quantiles"][0.5]

    @property
    def model(self):
        """The model document that ran, as a new dictionary.

        For a run from ``ensemble`` or ``prepare`` (also with a compiled
        loop) it is the parsed model bytes whose SHA-256 is
        ``metadata["model_sha256"]``. For a ``bind_parameters`` candidate
        it is the template's model with the bound parameter values, all
        parameters as floats in name order: the document the engine
        re-serialized and ran, apart from fields that serialization adds
        with their default values. None for a Result built directly from
        a metadata document and columns.
        """
        return None if self._ran is None else json.loads(self._ran[0])

    @property
    def request(self):
        """The request document that ran, as a new dictionary.

        For a ``bind_parameters`` candidate it is the template's request
        bound to the candidate's model: ``model.sha256`` is
        ``metadata["model_sha256"]``. None for a Result built directly from
        a metadata document and columns.
        """
        if self._ran is None:
            return None
        request = json.loads(self._ran[1])
        if self._ran[2]:
            request["model"]["sha256"] = self.metadata["model_sha256"]
        return request

    def save(self, path, *, paths=True):
        """Write the model, the request and the results to one report file.

        The extension of ``path`` decides the format: ``.md`` or
        ``.markdown`` for markdown, ``.json`` for JSON (schema
        ``crn-craft.report/1.0``, read back with :func:`load_report`), in
        any case; another extension raises ``ValueError``. ``paths=False``
        leaves the single trajectories out. An existing file is replaced;
        the report is written to a temporary file in the same directory and
        renamed, so a failed write (``OSError``) leaves no partial file.
        Returns ``path`` as a ``pathlib.Path``.

        The same result always gives the same bytes: there is no timestamp,
        host name or timing. The markdown report lists the model (reactions,
        start values, parameters, stop condition, continuous variables),
        the simulation settings with the sampler, the package version and
        both hashes, and the results: terminations, mean events per rule,
        the kept trajectories and the summary statistics (quantiles, when
        the result has them, in ascending order). Observables follow the
        model's species declaration order when every observable is named
        after a species type, then continuous variables; otherwise result
        order. Each section's lines are in a fenced ``text`` block, so a
        markdown renderer keeps one line per entry. Numbers are written in
        shortest round-trip form, integral values without ``.0``.

        A report is a readable record of a run; an archive is the sealed,
        verifiable record of a run for replay. ``save`` changes nothing
        about the result.
        """
        return _report.save(self, path, paths=paths)

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
        """Provenance for a DataArray: sampler, method, seed and hashes.

        Every value is a string or a number, so netCDF can store it: the
        method is JSON text (``"null"`` when the result reports none).
        """
        meta = self.metadata
        attrs = {key: meta.get(key) for key in _ARRAY_ATTRS}
        attrs["method"] = json.dumps(attrs["method"], sort_keys=True)
        return attrs

    @cached_property
    def paths(self):
        """Kept series as a DataArray: replicates x report times x variables.

        Dimensions are ``("replicate", "time", "variable")``, the order of
        ``event_columns()``, with the coordinates ``replicate`` (the kept
        replicate ids, int64), ``time`` (the report times) and ``variable``
        (the observables in result order: sorted names, then continuous
        variables). Values
        are float64. ``attrs`` holds the result's ``sampler``, ``method``
        (the method object of rules-result/0.3 as JSON text, ``"null"``
        for older results), ``root_seed``, ``model_sha256`` and
        ``run_sha256``; netCDF can store all of them.

        ``keep_paths`` decides which replicates are kept (``"all"`` keeps
        every one). A replicate that ended on a budget has values only up
        to its last saved report time and NaN after it. For example
        ``paths.sel(variable="mutants")`` is a replicate x time array and
        ``paths.sel(replicate=3)`` one replicate's series. Release 0.5.0
        ordered the dimensions ``("time", "variable", "replicate")``;
        ``paths.transpose("time", "variable", "replicate")`` gives that
        order. ``paths_frame()`` returns the DataFrame of release 0.4.0.
        """
        ids, lengths, positions, values = self._path_arrays
        times = np.asarray(self.metadata["report_times"], dtype=np.float64)
        observables = list(self.metadata["observables"])
        cube = np.full((len(ids), len(times), len(observables)), np.nan)
        cube[np.repeat(np.arange(len(ids)), lengths), positions, :] = values
        return xarray.DataArray(
            cube,
            dims=("replicate", "time", "variable"),
            coords={"replicate": ids, "time": times, "variable": observables},
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
        It holds the values of ``paths``: ``paths.to_dataset("variable")``.
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
    def __init__(
        self,
        native,
        rule_names,
        documents=None,
        compiled=None,
        *,
        ran=None,
        unified_from=None,
        unified=False,
    ):
        self._native = native
        self._rule_names = tuple(rule_names)
        # The one flag of the unified API: prepared through the top-level
        # names (crn_craft.ensemble, crn_craft.prepare), which write the 1.0
        # documents. Candidates of bind_parameters keep it.
        self._unified = unified
        # SHA-256 of the rules-model bytes that quantiles= rewrote to
        # crn-craft.model/1.0, so a loop compiled for them is explained.
        self._unified_from = unified_from
        # The exact (model, request) bytes, which a compiled loop runs.
        self._documents = documents
        # What results report as Result.model and Result.request:
        # (model bytes, request bytes, bound); see Result._ran.
        if ran is None and documents is not None:
            ran = (*documents, False)
        self._ran = ran
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
            hint = (
                "A compiled loop runs only the model bytes it was compiled for; "
                "compile one for this model and parameter set with "
                "rules.compile_loop."
            )
            if self._unified_from is not None and (
                compiled.model_sha256 == self._unified_from
            ):
                hint = (
                    "The loop was compiled for the rules-model bytes, but "
                    "quantiles run the model as crn-craft.model/1.0. Compile the "
                    "loop for the model with its schema set to "
                    "crn-craft.model/1.0, or run without quantiles."
                )
            raise CrnError(
                [
                    {
                        "code": "CRN-R020",
                        "path": "/model",
                        "message": "generated rules loop fingerprint mismatch",
                        "hint": hint,
                    }
                ]
            )
        return self._documents[1]

    def _result(self, metadata, columns):
        """A Result of this run that knows the documents that ran."""
        result = Result(json.loads(metadata), columns)
        result._ran = self._ran
        return result

    def archive(self, directory, *, cancel_file=None):
        """Run this request once more and seal the run in ``directory``.

        As :func:`archive` with the exact model and request bytes this run
        was prepared from: the archived result equals the result of
        ``simulate()`` (completed runs replay byte for byte; a compiled loop
        gives the interpreter's result, which the archive holds). Returns
        the ``crn-craft.archive-check/1.0`` report. A ``bind_parameters``
        candidate has no original bytes and raises ``ValueError``: archive
        its model with the bound values and a request bound to it with
        :func:`archive`.
        """
        if self._documents is None:
            raise ValueError(
                "a bind_parameters candidate has no original model and request "
                "bytes to archive; pass the model with the bound parameter values "
                "and a request bound to it to rules.archive"
            )
        return archive(*self._documents, directory, cancel_file=cancel_file)

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

        Its results' ``model`` is this run's model with the bound values,
        and their ``request`` this run's request bound to it (see
        ``Result.model``).
        """
        overrides = dict(overrides)
        if self._unified:
            # The top-level API refuses what float() would turn into a number.
            for key, value in overrides.items():
                if isinstance(value, bool | np.bool_) or not isinstance(
                    value, numbers.Real
                ):
                    raise CrnError(
                        [
                            {
                                "code": "CRN-R020",
                                "path": f"/parameters/{key}",
                                "message": "a bound parameter value must be a real "
                                f"number, not {type(value).__name__}",
                                "hint": "Pass a float or an int for every parameter.",
                            }
                        ]
                    )
        try:
            native = self._native.bind_parameters(overrides)
        except _native.NativeError as error:
            raise _crn_error(error) from None
        ran = None
        if self._ran is not None:
            # The engine re-serializes the parsed model with the bound values
            # (parameters are a sorted map of floats) and binds the request to
            # the new bytes; results carry that model hash.
            model = json.loads(self._ran[0])
            parameters = {**(model.get("parameters") or {}), **overrides}
            model["parameters"] = {k: float(v) for k, v in sorted(parameters.items())}
            ran = (model_bytes(model), self._ran[1], True)
        return PreparedEnsemble(
            native, self._rule_names, ran=ran, unified=self._unified
        )

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

        ``compiled``, a :class:`CompiledLoop` for this model, runs this
        call in the compiled program instead of the interpreter, with
        byte-identical results. The request must select ``static_direct``
        (``CrnError`` ``CRN-R020`` at ``/method`` otherwise) and use this
        loop's model bytes (``CRN-R020`` at ``/model``); the first call
        checks this. ``None``, the default, means the loop attached with
        ``compiled=`` on ``ensemble`` or ``prepare``, and the interpreter
        if none was attached. So a run prepared without a loop can be
        simulated either way, call by call, while a run with an attached
        loop always runs compiled. Callbacks need the interpreter and are
        refused with ``ValueError`` whenever a loop runs the call.
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
            return self._result(metadata, columns)

        def dispatch(replicate, event):
            return callback(replicate, json.loads(event))

        try:
            metadata, columns = self._native.simulate(
                _path(cancel_file), None if callback is None else dispatch
            )
        except _native.NativeError as error:
            raise _crn_error(error) from None
        result = self._result(metadata, columns)
        if batch_callback is not None:
            for observation in result.observations:
                log = observation["log"]
                for start in range(0, len(log), batch_size):
                    batch_callback(observation["id"], log[start : start + batch_size])
        return result


def prepare(model, request, *, compiled=None, quantiles=None, selector=None):
    """Prepare a rules-run contract, preserving the supplied model bytes.

    ``model`` is bytes, text, a dictionary or a ``Model``. ``compiled`` is
    an optional :class:`CompiledLoop` for this model, attached as in
    :func:`ensemble`: the request is checked against it here, once, and
    every ``simulate()`` of the prepared run uses it.

    ``quantiles``, a sequence of increasing probabilities in (0, 1), adds
    ``outputs.quantiles`` to the request, as in :func:`ensemble`, and
    ``selector`` (one of :data:`SELECTORS`) adds the event selector to a
    ``static_direct`` method (``ValueError`` for another method). Both exist
    only in the 1.0 contracts: the request is written as
    ``crn-craft.run/1.0`` and a rules-model as ``crn-craft.model/1.0`` (only
    its schema string changes; the request's model hash follows when it
    was the hash of the supplied bytes). A Static Direct request of a
    rules-run version keeps its meaning: the defaults of its version are
    written out (``"waiting_time": "inverse"`` and ``"sparse_trackers":
    False`` for rules-run/0.2 to 0.5, and ``"selector": "linear"`` unless
    ``selector`` is given), so it runs the same paths. A request that
    already has ``outputs.quantiles`` (or a selector) raises ``ValueError``
    when the argument is given too. ``None`` or an empty sequence, the
    default, leaves both documents as supplied; a ``crn-craft.run/1.0``
    request may also carry ``outputs.quantiles`` and ``method.selector``
    itself.
    """
    return _prepare(model, request, compiled, quantiles, _selector(selector), False)


def _prepare(model, request, compiled, quantiles, selector, unified):
    """``prepare``; with ``unified`` (the top-level name) every request is
    written as run/1.0 and every model as model/1.0, ``quantiles`` may be
    :data:`_DEFAULT` (the default quantiles unless the request states its
    own) and network models are converted."""
    if quantiles is _DEFAULT:
        quantiles = () if _has_quantiles(request) else DEFAULT_QUANTILES
    probabilities = _probabilities(quantiles)
    if isinstance(model, Model):
        model = model.to_dict()
    unified_from = None
    original = model_bytes(model)
    if unified or probabilities or selector is not None:
        data = _as_model_1_0(model, network=unified)
        old, new = (hashlib.sha256(d).hexdigest() for d in (original, data))
        schema = _request_schema(request)
        if unified and (
            schema is _UNREADABLE
            or (
                schema == RUN_1_0
                and not probabilities
                and selector is None
                and old == new
            )
        ):
            # Already the 1.0 documents, or a request that does not parse:
            # its exact bytes go to the engine.
            pass
        else:
            request = _request_1_0(request, probabilities, selector, old, new)
        unified_from = old if old != new else None
    else:
        data = original
    return _prepared(data, request, compiled, unified_from, unified)


# `_request_schema` of a request that is no JSON object.
_UNREADABLE = object()


def _request_schema(request):
    """The request's ``schema`` value, read leniently; :data:`_UNREADABLE`
    for a request that is no JSON object."""
    if isinstance(request, bytes | str):
        try:
            request = json.loads(request)
        except ValueError:
            return _UNREADABLE
    return request.get("schema") if isinstance(request, Mapping) else _UNREADABLE


def _has_quantiles(request):
    """Whether a request (bytes, text or a dictionary) states outputs.quantiles."""
    if isinstance(request, bytes | str):
        try:
            request = json.loads(request)
        except ValueError:
            return False
    outputs = request.get("outputs") if isinstance(request, Mapping) else None
    return isinstance(outputs, Mapping) and "quantiles" in outputs


def _prepared(data, request, compiled, unified_from, unified):
    """Prepare exact model bytes with a request (bytes or a dictionary)."""
    request = model_bytes(request)
    try:
        native = _native.prepare_rules(data, request)
    except _native.NativeError as error:
        raise _crn_error(error) from None
    # The engine has validated the model; results list rules in this order.
    names = [r["name"] for r in json.loads(data)["rules"]]
    return PreparedEnsemble(
        native,
        names,
        (data, request),
        compiled,
        unified_from=unified_from,
        unified=unified,
    )


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
    quantiles=None,
    selector=None,
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
    used; ``"static_direct"`` selects rules-run/0.6 for its defaults, and a
    ``crn-craft.model/1.0`` model (from :func:`crn_craft.convert`) always
    selects ``crn-craft.run/1.0``, which has the fields and defaults of
    rules-run/0.6.

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
    :func:`compile_loop` for this model. It attaches the loop to the
    prepared run: preparation checks once that the loop was compiled for
    these model bytes (``CrnError`` ``CRN-R020`` at ``/model`` otherwise)
    and that the request selects ``static_direct`` with any options
    (``CRN-R020`` at ``/method``), and every ``simulate()`` then runs in the
    compiled program, with results byte-identical to the interpreter's.
    Attach it when every run of this request should use the loop, for
    example when ``simulate()`` is called repeatedly or by code that knows
    nothing of compiled loops; a mistake then fails before anything runs.
    ``simulate(compiled=loop)`` uses a loop for one call instead and leaves
    the run interpreted otherwise, which suits comparing the two. A run
    with an attached loop cannot run interpreted or take callbacks; prepare
    it without ``compiled`` for those. Candidates from
    ``bind_parameters`` run interpreted, since a loop serves only the
    parameter set it was compiled for.

    ``quantiles`` is a sequence of at most 16 strictly increasing
    probabilities in (0, 1), for example ``(0.025, 0.5, 0.975)``. The
    engine then also reports, for every report time and observable, the
    quantile of the replicates' values: ``result.statistics["quantiles"]``
    is a dict from each probability to a DataFrame shaped like
    ``statistics["mean"]``, and ``statistics["median"]`` is
    ``quantiles[0.5]`` when 0.5 was requested. The quantile is NumPy's
    default (``numpy.quantile(x, p)``, linear interpolation, type 7) and
    equals it bit for bit. Quantiles are withheld with the other statistics
    (``statistics`` is None when a replicate did not complete). The engine
    keeps replicates x report times x observables values for them during
    the run; above its 2 GiB workspace limit the request is refused with
    ``CRN-R020`` at ``/outputs/quantiles``.

    The default, ``None``, and an empty sequence request no quantiles: the
    documents are then exactly those written without the argument, as in
    release 0.6.0. Quantiles exist only in the 1.0 contracts: with them the request is written as
    ``crn-craft.run/1.0`` (the fields and defaults of rules-run/0.6, so the
    paths are those of the request without quantiles) and a rules-model,
    such as a :class:`Model`, as ``crn-craft.model/1.0``, which means what
    rules-model/0.4 means: only its schema string changes. The result is
    then ``crn-craft.result/1.0`` with new hashes, and ``result.model`` is
    the model/1.0 document. A compiled loop serves exact model bytes, so
    compile it for the model/1.0 document to use it with quantiles (for a
    :class:`Model`: ``d = m.to_dict(); d["schema"] = "crn-craft.model/1.0";
    loop = compile_loop(d)``); a loop compiled for the rules-model bytes is
    refused with ``CRN-R020`` at ``/model``.

    ``selector`` chooses the event selector of ``method="static_direct"``
    (``ValueError`` for another method): ``"linear"``, the scan of every
    rules-run version; ``"composition_rejection"``, ``"sorting"`` and
    ``"rssa_cr"``, ported from the network engine; or ``"auto"``, which
    picks one of them from the number of lowered channels and the initial
    populations (``crn_craft.describe()["selectors"]["auto"]`` lists the
    thresholds). The selector exists only in ``crn-craft.run/1.0``, so
    passing it writes the 1.0 documents, as ``quantiles`` does. In this
    module a Static Direct request written as ``crn-craft.run/1.0`` without
    a selector gets ``"selector": "linear"`` written out, so adding
    ``quantiles`` never changes the paths; the top-level
    :func:`crn_craft.ensemble` leaves the selector out and runs ``"auto"``.
    The selector is part of the sampler identity.
    """
    return _ensemble(
        model,
        end=end,
        report_times=report_times,
        replicates=replicates,
        seed=seed,
        threads=threads,
        keep_paths=keep_paths,
        keep_instances=keep_instances,
        max_events=max_events,
        max_sampling_steps=max_sampling_steps,
        max_instances=max_instances,
        max_output_bytes=max_output_bytes,
        wall_time_s=wall_time_s,
        trackers=trackers,
        event_log=event_log,
        method=method,
        compiled=compiled,
        quantiles=quantiles,
        selector=_selector(selector),
        unified=False,
    )


def _ensemble(
    model,
    *,
    end,
    report_times,
    replicates,
    seed,
    threads,
    keep_paths,
    keep_instances,
    max_events,
    max_sampling_steps,
    max_instances,
    max_output_bytes,
    wall_time_s,
    trackers,
    event_log,
    method,
    compiled,
    quantiles,
    selector,
    unified,
):
    """``ensemble``; with ``unified`` (the top-level name) the documents are
    always model/1.0 and run/1.0 and network models are converted."""
    probabilities = _probabilities(quantiles)
    if isinstance(model, Model):
        model = model.to_dict()
    unified_from = None
    if unified or probabilities or selector is not None:
        data = _as_model_1_0(model, network=unified)
        original = hashlib.sha256(model_bytes(model)).hexdigest()
        if original != hashlib.sha256(data).hexdigest():
            unified_from = original
    else:
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
    if selector is not None:
        request["method"] = _with_selector(request.get("method"), selector)
    if (
        json.loads(data).get("schema") == MODEL_1_0
        or probabilities
        or selector is not None
        or unified
    ):
        # A crn-craft.model/1.0 model runs only with crn-craft.run/1.0,
        # which has the fields and defaults of rules-run/0.6; quantiles and
        # selectors exist only there. Static Direct is already rules-run/0.6
        # above.
        request["schema"] = RUN_1_0
    if (
        not unified
        and request["schema"] == RUN_1_0
        and isinstance(request.get("method"), dict)
        and request["method"].get("kind") == "static_direct"
    ):
        # crn_craft.rules keeps the linear scan of every rules-run version;
        # a run/1.0 request without a selector would run "auto".
        request["method"].setdefault("selector", "linear")
    if probabilities:
        request["outputs"]["quantiles"] = probabilities
    return _prepared(data, request, compiled, unified_from, unified)


def fsp(model, end, report_times, bounds, max_states=None, tolerance=None):
    """Finite state projection: the distribution of a closed model's counts.

    ``model`` is model bytes, a JSON string, a dictionary or a
    :class:`Model`, of any rules-model version or ``crn-craft.model/1.0``
    (from :func:`crn_craft.convert`). The engine lowers it to its flat
    network, as the static simulation methods do, enumerates the states
    reachable from the initial counts inside the box ``bounds`` and
    propagates their probabilities by uniformization to each of
    ``report_times`` (increasing times in [0, ``end``]). Probability that
    leaves the box moves to an absorbing sink and is reported, not
    renormalized; a state where the model's ``stop`` condition holds is
    absorbing.

    ``bounds`` maps every species label to its inclusive upper count, at
    least the initial count: the type name for a type without attributes
    (``{"S": 999, "I": 1000, "R": 1000}``), ``"T{a=1,b=0.5}"`` for an
    attributed type (every attribute in declaration order). ``max_states``
    caps the enumerated states (default 100,000, which the engine never
    exceeds) and
    ``tolerance`` is the total truncation budget of the Poisson series over
    the run (default ``1e-12``); neither bounds floating-point error.

    Returns the parsed ``crn-craft.fsp-result/1.0`` document: ``species``
    (the labels), ``states`` (count vectors in discovery order, the initial
    state first), and one entry of ``snapshots`` per report time with
    ``probabilities`` (indexed like ``states``), ``physical_leakage_estimate``
    (the sink), the unnormalized ``projected_first_moment`` and
    ``projected_second_raw_moment`` per species, and the error terms
    ``series_tail_upper_bound``, ``discarded_weight_estimate``,
    ``mass_defect`` and ``roundoff_mass_residual``. A model that is not
    closed raises ``CrnError`` ``CRN-R017``; a request the engine refuses
    raises ``CRN-R020``; a projection above the state, edge, work or output
    limits raises ``CRN-E070``. Ctrl-C cancels the solve.
    """
    if isinstance(model, Model):
        model = model.to_dict()
    data = model_bytes(model)
    if not isinstance(bounds, Mapping):
        raise TypeError(
            f"bounds must map species labels to counts, for example {{'A': 100}}: {bounds!r}"
        )
    count = operator.index
    request = {
        "schema": RUN_1_0,
        "model": {
            "path": "model.rules.json",
            "sha256": hashlib.sha256(data).hexdigest(),
        },
        "time": {"end": float(end), "report_times": [float(t) for t in report_times]},
        "method": {
            "kind": "fsp",
            "bounds": {str(label): count(bound) for label, bound in bounds.items()},
            "max_states": count(FSP_MAX_STATES if max_states is None else max_states),
            "tolerance": float(FSP_TOLERANCE if tolerance is None else tolerance),
        },
        "budget": {"max_output_bytes": 1 << 30},
    }
    try:
        text = _native.rules_fsp_json(data, model_bytes(request))
    except _native.NativeError as error:
        raise _crn_error(error) from None
    return json.loads(text)


# The error terms of an FSP snapshot, the columns of FspResult.truncation.
_FSP_ERRORS = (
    "physical_leakage_estimate",
    "series_tail_upper_bound",
    "discarded_weight_estimate",
    "mass_defect",
    "roundoff_mass_residual",
)


class FspResult:
    """A finite state projection as pandas and xarray objects.

    :meth:`Model.fsp` returns it; :func:`fsp` returns the document it wraps.

    - ``states``: the projected states, a DataFrame indexed by ``state``
      (0 is the initial state) with one int64 column per species label.
    - ``probabilities``: a DataArray with dimensions ``("time", "state")``,
      the probability of each state at each report time (the mass that left
      the bounds is not in it; see ``truncation``).
    - ``mean``: the projected first moment of each species at each report
      time, a DataFrame indexed by ``time``: Σ p(x)·x over the retained
      states, not renormalized, so it is the mean of the full process only
      up to the leaked mass.
    - ``truncation``: the error terms per report time, a DataFrame with the
      columns ``physical_leakage_estimate`` (mass that left the bounds),
      ``series_tail_upper_bound``, ``discarded_weight_estimate``,
      ``mass_defect`` and ``roundoff_mass_residual``.
    - ``species``: the species labels; ``metadata``: the
      ``crn-craft.fsp-result/1.0`` document.
    """

    def __init__(self, document):
        self.metadata = document
        self.species = tuple(document["species"])
        snapshots = document["snapshots"]
        times = pd.Index([s["time"] for s in snapshots], dtype=np.float64, name="time")
        columns = pd.Index(self.species, name="species")
        count = len(document["states"])
        self.states = pd.DataFrame(
            np.array(document["states"], dtype=np.int64).reshape(
                count, len(self.species)
            ),
            index=pd.RangeIndex(count, name="state"),
            columns=columns,
        )
        self.probabilities = xarray.DataArray(
            np.array([s["probabilities"] for s in snapshots], dtype=np.float64).reshape(
                len(snapshots), count
            ),
            dims=("time", "state"),
            coords={"time": times.to_numpy(), "state": np.arange(count)},
            name="probability",
        )
        self.mean = pd.DataFrame(
            np.array(
                [s["projected_first_moment"] for s in snapshots], dtype=np.float64
            ).reshape(len(snapshots), len(self.species)),
            index=times,
            columns=columns,
        )
        self.truncation = pd.DataFrame(
            {key: [s[key] for s in snapshots] for key in _FSP_ERRORS},
            index=times,
            dtype=np.float64,
        )

    def __repr__(self):
        return (
            f"FspResult({len(self.states)} states of {', '.join(self.species)}, "
            f"{len(self.truncation)} report times)"
        )


def archive(model, request, directory, *, cancel_file=None):
    """Run a request and seal the run in a new archive directory.

    ``model`` is model bytes, a JSON string, a dictionary or a
    :class:`Model`; ``request`` is request bytes, a JSON string or a
    dictionary: any rules-run version, or ``crn-craft.run/1.0`` with any
    method, ``fsp`` included. The request must bind the exact model bytes
    (``model.sha256``); its ``model.path`` is kept but not resolved.
    ``directory`` must not exist; its parent must.

    The directory gets ``model.json`` and ``request.json`` (the exact
    bytes), ``result.json`` (the result as ``crn-craft rules`` prints it),
    ``manifest.json`` (``crn-craft.archive/1.0``: length and SHA-256 of the
    three) and, written last, ``COMPLETE``. This is
    ``crn-craft rules MODEL REQUEST --archive NEW_DIRECTORY`` without files
    for the inputs; :meth:`PreparedEnsemble.archive` archives a prepared
    run.

    Returns the ``crn-craft.archive-check/1.0`` report as a dictionary,
    also when the run or the write failed: read ``valid`` (the archive was
    written), ``simulation_complete`` (the run completed) and
    ``termination`` together, and ``diagnostics`` for the reason of a
    failure. Input errors create no directory. Check an archive later with
    :func:`verify_archive`. Ctrl-C cancels the run, which is archived as
    cancelled before ``KeyboardInterrupt`` is raised.

    An archive is the sealed, verifiable record of a run for replay; a
    report written by :meth:`Result.save` is a readable record.
    """
    if isinstance(model, Model):
        model = model.to_dict()
    return json.loads(
        _native.rules_archive_json(
            model_bytes(model),
            model_bytes(request),
            str(directory),
            _path(cancel_file),
        )
    )


def verify_archive(directory):
    """Check an archive directory without running anything.

    Reads ``crn-craft.archive/1.0`` directories (written by :func:`archive`
    and ``crn-craft run ... --archive``) and returns the report of
    ``crn-craft verify DIRECTORY`` as a dictionary
    (``crn-craft.archive-check/1.0``): ``valid`` is true for an intact,
    consistent archive, also of an incomplete run (``simulation_complete``
    false); ``diagnostics`` names the first inconsistency otherwise. The
    network engine's ``crn-craft.archive/0.1`` directories are no longer
    read (release 0.7.0).
    """
    return json.loads(_native.verify_json(str(directory)))
