"""Helpers shared by the public modules: the error type, exact model bytes,
native calls and the converter.

``crn_craft`` and ``crn_craft.rules`` both use these objects, so
``crn_craft.CrnError is crn_craft.rules.CrnError``.
"""

import importlib.metadata
import json

from . import _native

try:
    __version__ = importlib.metadata.version("crn-craft")
except importlib.metadata.PackageNotFoundError:  # a source tree, not installed
    __version__ = _native.__version__


class CrnError(ValueError):
    """A request the native core rejected.

    ``diagnostics`` is the list of CRN-craft diagnostics (dictionaries with
    ``code``, ``path``, ``message`` and ``hint``). It is a ``ValueError``, so
    existing ``except ValueError`` handlers keep working.
    """

    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        first = diagnostics[0] if diagnostics else {}
        message = f"{first.get('code', 'CRN')}: {first.get('message', 'rejected')}"
        if first.get("path"):
            message += f" at {first['path']}"
        if len(diagnostics) > 1:
            message += f" (and {len(diagnostics) - 1} more)"
        super().__init__(message)


def _crn_error(error):
    text = str(error)
    try:
        value = json.loads(text)
    except ValueError:
        value = [{"code": "CRN-E000", "path": "", "message": text}]
    return CrnError(value if isinstance(value, list) else [value])


def model_bytes(model):
    """The exact bytes of a document: bytes as given, text as UTF-8, and a
    dictionary as compact JSON with a trailing newline."""
    if isinstance(model, bytes):
        return model
    if isinstance(model, str):
        return model.encode("utf-8")
    return (json.dumps(model, allow_nan=False, separators=(",", ":")) + "\n").encode()


def _path(value):
    return None if value is None else str(value)


def convert(model):
    """Convert a network model (``crn-craft.model/0.1``) into ``crn-craft.model/1.0``.

    ``model`` is bytes, text or a dictionary. Returns the converted model's
    exact bytes (pretty-printed JSON with a trailing newline), the same
    bytes ``crn-craft convert MODEL`` prints. Every species becomes a type
    and an observable of the same name, every reaction a rule with the same
    propensity; reactions of any order use the pattern field ``count``. A
    model the converter cannot hold (a population limit, a network without
    reactions, ...) raises ``CrnError`` with ``CRN-R016``; an invalid
    network model raises ``CrnError`` with its validation diagnostics.
    ``Model.from_dict(json.loads(convert(m)))`` gives a :class:`Model`.
    """
    try:
        return _native.convert_network(model_bytes(model))
    except _native.NativeError as error:
        raise _crn_error(error) from None
