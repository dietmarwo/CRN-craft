"""Readable reports of a rules result: ``Result.save`` and ``load_report``.

A report holds the model and request that ran and the results in one file,
as markdown (``.md``, ``.markdown``) or as JSON (``.json``, schema
``crn-craft.report/1.0``). The same result always gives the same bytes: a
report has no timestamp, host name or timing. It is a readable record; an
archive is the sealed, verifiable record of a run for replay.

Numbers in the results are written in their shortest round-trip form
(``repr`` of a float), except that a float with an integral value below 2**53
is written without ``.0``. Non-finite values are ``nan``, ``inf`` and
``-inf``; in JSON they are strings, so the file is strict JSON. The model and
request in a JSON report are the documents as they ran.
"""

import json
import math
import numbers
import os
from contextlib import suppress
from decimal import Decimal
from pathlib import Path

SCHEMA = "crn-craft.report/1.0"
MARKDOWN = (".md", ".markdown")
JSON = (".json",)
# Statistics in report order; quantiles follow in ascending order.
_STATISTICS = ("mean", "mean_standard_error", "sample_variance")
# Rule fields a reaction line can show; anything else is written as JSON.
_RULE_FIELDS = frozenset(("name", "reactants", "rate", "propensity", "products"))
_PATTERN_FIELDS = frozenset(("label", "type", "count"))
_PRODUCT_FIELDS = frozenset(("type", "attributes"))
# Fields of a continuous variable, in the order a report lists them.
_CONTINUOUS_FIELDS = ("initial", "derivative", "lower_bound", "on_bound", "independent")
_LIMIT = 2**53


def _non_finite(value):
    """``nan``, ``inf`` or ``-inf`` for a non-finite float, else None."""
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return None


def number(value):
    """The report text of a number (see the module documentation)."""
    if isinstance(value, numbers.Integral):
        return str(int(value))
    value = float(value)
    special = _non_finite(value)
    if special is not None:
        return special
    text = repr(value)
    if value.is_integer() and abs(value) < _LIMIT and text.endswith(".0"):
        return text[:-2]
    return text


def _document_number(value):
    """A number of the model or request: the standard JSON form."""
    if isinstance(value, numbers.Integral):
        return str(int(value))
    value = float(value)
    special = _non_finite(value)
    return json.dumps(special) if special is not None else repr(value)


def _result_number(value):
    """A number of the results in JSON: the report form, non-finite quoted."""
    text = number(value)
    return json.dumps(text) if text in ("nan", "inf", "-inf") else text


def _value(value):
    """A model value in markdown: an expression as written, or a number."""
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, numbers.Real):
        return number(value)
    return _compact(value)


def _compact(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


# ---------------------------------------------------------------------------
# JSON


def _holds_object(value):
    if isinstance(value, dict):
        return True
    return isinstance(value, list) and any(_holds_object(v) for v in value)


def _json(value, write_number, level):
    """JSON text: one object member per line, lists of numbers inline."""
    if isinstance(value, dict):
        if not value:
            return "{}"
        pad = "  " * (level + 1)
        members = [
            f"{pad}{json.dumps(str(key))}: {_json(item, write_number, level + 1)}"
            for key, item in value.items()
        ]
        return "{\n" + ",\n".join(members) + "\n" + "  " * level + "}"
    if isinstance(value, list | tuple):
        if all(type(v) is float or type(v) is int for v in value):
            return "[" + ", ".join(map(write_number, value)) + "]"
        if not any(_holds_object(v) for v in value):
            return "[" + ", ".join(_json(v, write_number, level) for v in value) + "]"
        pad = "  " * (level + 1)
        items = [pad + _json(v, write_number, level + 1) for v in value]
        return "[\n" + ",\n".join(items) + "\n" + "  " * level + "]"
    if isinstance(value, str):
        return json.dumps(value)
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, numbers.Real):
        return write_number(value)
    raise TypeError(f"a report cannot hold {type(value).__name__} values")


def _table(frame, observables):
    """A statistic as {observable: values over the report times}."""
    if frame is None:
        return None
    return {name: frame[name].to_numpy().tolist() for name in observables}


def _quantiles(stats):
    """(p, DataFrame) of the result's quantiles, in ascending order."""
    return sorted(
        ((float(p), frame) for p, frame in (stats.get("quantiles") or {}).items()),
        key=lambda item: item[0],
    )


def _kept(result):
    """Kept replicate ids (ascending) and their values, replicate x time x
    variable as nested lists, or None if no path was kept."""
    if not result.metadata["paths"]:
        return None
    cube = result.paths
    ids = cube["replicate"].values.tolist()
    values = cube.values.tolist()
    order = sorted(range(len(ids)), key=ids.__getitem__)
    return [ids[i] for i in order], [values[i] for i in order]


def _json_report(result, model, request, version, paths):
    meta = result.metadata
    observables = list(meta["observables"])
    stats = result.statistics
    if stats is None:
        statistics = None
    else:
        statistics = {"samples": stats["samples"]}
        for key in _STATISTICS:
            statistics[key] = _table(stats.get(key), observables)
        statistics["quantiles"] = [
            {"p": p, "values": _table(frame, observables)}
            for p, frame in _quantiles(stats)
        ]
    kept = _kept(result) if paths else None
    if kept is not None:
        ids, values = kept
        kept = {
            "replicates": ids,
            "values": {
                name: [[row[column] for row in rows] for rows in values]
                for column, name in enumerate(observables)
            },
        }
    body = {
        "model_sha256": meta["model_sha256"],
        "run_sha256": meta["run_sha256"],
        "sampler": meta["sampler"],
        "report_times": list(meta["report_times"]),
        "observables": observables,
        "terminations": _terminations(result),
        "rule_events_mean": _mean_events(result),
        "statistics": statistics,
        "paths": kept,
    }
    members = [
        ("schema", SCHEMA, _document_number),
        ("crn_craft_version", version, _document_number),
        ("model", model, _document_number),
        ("request", request, _document_number),
        ("result", body, _result_number),
    ]
    text = ",\n".join(
        f"  {json.dumps(key)}: {_json(value, write, 1)}"
        for key, value, write in members
    )
    return "{\n" + text + "\n}\n"


# ---------------------------------------------------------------------------
# Markdown


def _terminations(result):
    """Counts of the terminations that occurred, in category order."""
    column = result.records["termination"]
    counts = column.value_counts()
    return {
        str(name): int(counts[name])
        for name in column.cat.categories
        if int(counts[name]) > 0
    }


def _mean_events(result):
    """The mean over replicates of each rule's event count."""
    events = result.rule_events
    size = len(events)
    return {
        name: int(events[name].to_numpy().sum()) / size if size else math.nan
        for name in result.rule_names
    }


def _order(result, model):
    """Observables in report order (see Result.save)."""
    observables = list(result.metadata["observables"])
    continuous = set(result.metadata.get("continuous") or ())
    counted = [name for name in observables if name not in continuous]
    types = [t.get("name") for t in model.get("types", [])]
    if not all(name in types for name in counted):
        return observables
    declared = [name for name in types if name in counted]
    return declared + [name for name in observables if name in continuous]


def _attributes(model):
    return {
        t.get("name"): list(t.get("attributes") or ()) for t in model.get("types", [])
    }


def _reaction(rule, attributes):
    """One reaction line, or the rule's name and compact JSON."""
    reactants = rule.get("reactants", [])
    products = rule.get("products", [])
    plain = (
        set(rule) <= _RULE_FIELDS
        and ("rate" in rule) != ("propensity" in rule)
        and all(
            set(p) <= _PATTERN_FIELDS and not attributes.get(p.get("type"))
            for p in reactants
        )
        and all(
            set(p) <= _PRODUCT_FIELDS
            and not p.get("attributes")
            and not attributes.get(p.get("type"))
            for p in products
        )
    )
    if not plain:
        return f"{rule.get('name')}: {_compact(rule)}"
    left = [p["type"] for p in reactants for _ in range(int(p.get("count", 1)))]
    right = [p["type"] for p in products]
    if "propensity" in rule:
        law = _value(rule["propensity"])
    else:
        law = _value(rule["rate"]) + " (mass action)"
    return f"{' + '.join(left) or '0'} -> {' + '.join(right) or '0'}, {law}"


def _start_values(model, attributes):
    lines = []
    initial = model.get("initial", [])
    for name, declared in attributes.items():
        entries = [e for e in initial if e.get("type") == name]
        if declared and entries:
            for entry in entries:
                values = entry.get("attributes") or {}
                keys = [a for a in declared if a in values]
                keys += [k for k in values if k not in declared]
                inner = ", ".join(f"{k}={_value(values[k])}" for k in keys)
                lines.append(f"{name} {{{inner}}} = {_value(entry.get('count'))}")
            continue
        counts = [e.get("count") for e in entries]
        if not counts:
            text = "0"
        elif len(counts) == 1:
            text = _value(counts[0])
        elif all(isinstance(c, numbers.Real) for c in counts):
            text = number(sum(counts))
        else:
            text = " + ".join(_value(c) for c in counts)
        lines.append(f"{name} = {text}")
    return lines


def _plain_observable(name, spec):
    return (
        isinstance(spec, dict)
        and spec.get("type") == name
        and set(spec) <= {"type", "reduce"}
        and spec.get("reduce", "sum") == "sum"
    )


def _method(method):
    if method is None:
        return "tree"
    if isinstance(method, str):
        return method
    options = [
        f"{key} = {_value(value)}" for key, value in method.items() if key != "kind"
    ]
    return ", ".join([str(method.get("kind")), *options])


def _withheld(result, request):
    """Why a result has no statistics."""
    termination = result.metadata.get("termination")
    if termination not in (None, "completed"):
        return f"ensemble termination: {termination}"
    if not (request.get("outputs") or {}).get("statistics"):
        return "no statistics requested"
    return "not in the result"


def _percentile(p):
    """``median`` for 0.5, else ``<100 p>% percentile`` with the shortest
    exact decimal."""
    if p == 0.5:
        return "median"
    percent = (Decimal(repr(p)) * 100).normalize()
    return f"{percent:f}% percentile"


def _markdown_report(result, model, request, version, paths):
    meta = result.metadata
    attributes = _attributes(model)
    out = ["# Model", "", f"Name: {model.get('name')}", ""]

    def section(title, lines):
        # The data lines go in a fenced block, so a renderer keeps one line
        # per entry. Names are identifiers and the engine has validated every
        # expression, so no data line contains a backtick.
        out.extend([title, "", "```text", *lines, "```", ""])

    if any(attributes.values()):
        section(
            "## Species types",
            [
                f"{name}: {', '.join(declared)}" if declared else str(name)
                for name, declared in attributes.items()
            ],
        )
    section("## Reactions", [_reaction(r, attributes) for r in model.get("rules", [])])
    section("## Start values", _start_values(model, attributes))
    parameters = model.get("parameters") or {}
    if parameters:
        section(
            "## Parameter values",
            [f"{name} = {_value(parameters[name])}" for name in sorted(parameters)],
        )
    observables = model.get("observables") or {}
    if not all(_plain_observable(name, spec) for name, spec in observables.items()):
        section(
            "## Observables",
            [f"{name}: {_compact(spec)}" for name, spec in observables.items()],
        )
    if model.get("stop") is not None:
        section("## Stop condition", [_value(model["stop"])])
    if model.get("continuous"):
        lines = []
        for variable in model["continuous"]:
            keys = [k for k in _CONTINUOUS_FIELDS if k in variable]
            keys += [k for k in variable if k != "name" and k not in keys]
            fields = ", ".join(f"{k} = {_value(variable[k])}" for k in keys)
            lines.append(f"{variable.get('name')}: {fields}")
        section("## Continuous variables", lines)
    if model.get("atp") is not None:
        section("## ATP coupling", [_compact(model["atp"])])
    ensemble = request.get("ensemble") or {}
    section(
        "# Simulation",
        [
            f"seed = {_value(ensemble.get('root_seed'))}",
            f"end = {_value((request.get('time') or {}).get('end'))}",
            f"replicates = {_value(ensemble.get('replicates'))}",
            "report_times = " + ", ".join(map(number, meta["report_times"])),
            "method = " + _method(request.get("method") or meta.get("method")),
            f"sampler = {meta['sampler']}",
            f"crn-craft = {version}",
            f"model_sha256 = {meta['model_sha256']}",
            f"run_sha256 = {meta['run_sha256']}",
        ],
    )
    out += ["# Results", ""]
    section(
        "## Terminations",
        [f"{name} = {count}" for name, count in _terminations(result).items()],
    )
    section(
        "## rule_events (mean per replicate)",
        [f"{name} = {number(mean)}" for name, mean in _mean_events(result).items()],
    )
    order = _order(result, model)
    columns = [list(meta["observables"]).index(name) for name in order]
    kept = _kept(result) if paths else None
    if kept is not None:
        ids, values = kept
        width = max(3, len(str(max(ids))))
        lines = []
        for replicate, rows in zip(ids, values, strict=True):
            for name, column in zip(order, columns, strict=True):
                series = ", ".join(number(row[column]) for row in rows)
                lines.append(f"{name}{replicate:0{width}d}: {series}")
        section("## Single trajectories for " + ", ".join(order), lines)
    stats = result.statistics
    if stats is None:
        section(
            "## Summary statistics",
            [f"statistics withheld ({_withheld(result, request)})"],
        )
    else:
        section("## Summary statistics", [f"samples = {number(stats['samples'])}"])
        shown = [(f"### {key}", stats.get(key)) for key in _STATISTICS]
        shown += [(f"### {_percentile(p)}", frame) for p, frame in _quantiles(stats)]
        for title, frame in shown:
            if frame is None:
                continue
            section(
                title,
                [
                    f"{name}: "
                    + ", ".join(map(number, frame[name].to_numpy().tolist()))
                    for name in order
                ],
            )
    return "\n".join(out[:-1]) + "\n"


# ---------------------------------------------------------------------------
# Files


def _write(path, data):
    """Write ``data`` to ``path`` through a temporary file and a rename."""
    directory = path.parent
    for _ in range(64):
        temporary = directory / f".{path.name}.{os.urandom(6).hex()}.tmp"
        try:
            # Created like open(path, "wb") would, so the umask applies.
            handle = os.open(
                temporary,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0),
                0o666,
            )
        except FileExistsError:
            continue
        break
    else:
        raise FileExistsError(f"no free temporary file name next to {path}")
    try:
        with os.fdopen(handle, "wb") as file:
            file.write(data)
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary, path)
    except BaseException:
        with suppress(OSError):
            os.unlink(temporary)
        raise


def save(result, path, *, paths=True):
    """Write a report of ``result`` to ``path``; see ``Result.save``."""
    target = Path(path)
    suffix = target.suffix.lower()
    if suffix not in MARKDOWN + JSON:
        raise ValueError(
            f"unsupported report extension {suffix or '(none)'!r} in {str(path)!r}: "
            "use .md, .markdown or .json"
        )
    if not isinstance(paths, bool):
        raise TypeError(f"paths must be True or False, not {type(paths).__name__}")
    model, request = result.model, result.request
    if model is None or request is None:
        raise ValueError(
            "this Result does not know the model and request that ran; save the "
            "result of a prepared run's simulate()"
        )
    from . import __version__

    write = _markdown_report if suffix in MARKDOWN else _json_report
    _write(target, write(result, model, request, __version__, paths).encode("utf-8"))
    return target


def _unique(pairs):
    document = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate key {key!r} in a report")
        document[key] = value
    return document


def _constant(name):
    raise ValueError(f"non-standard JSON constant {name} in a report")


def load_report(path):
    """Read a JSON report written by ``Result.save``.

    Returns the parsed ``crn-craft.report/1.0`` document as a dictionary:
    ``schema``, ``crn_craft_version``, ``model`` and ``request`` (the
    documents that ran) and ``result``. Non-finite numbers in it are the
    strings ``"nan"``, ``"inf"`` and ``"-inf"``. Raises ``ValueError`` for
    a file that is not strict JSON (including duplicate keys) or whose
    ``schema`` is not ``crn-craft.report/1.0``.
    """
    document = json.loads(
        Path(path).read_bytes(), object_pairs_hook=_unique, parse_constant=_constant
    )
    schema = document.get("schema") if isinstance(document, dict) else None
    if schema != SCHEMA:
        raise ValueError(f"{str(path)!r} is not a {SCHEMA} report (schema {schema!r})")
    return document
