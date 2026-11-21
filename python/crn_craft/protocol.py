"""Shared intervention descriptions for trajectory and distribution solvers."""

from copy import deepcopy
from operator import index


class Protocol:
    """An ordered calendar in the time unit of the model's rates.

    Use :meth:`at` for external interventions and :meth:`delay` to attach
    completions to the firings of a reaction. Equal times preserve insertion
    order. Reports see the state after actions and due completions.
    """

    def __init__(self, *, time_unit=None, max_pending=4096):
        self._document = {
            "actions": [],
            "delays": [],
            "max_pending": index(max_pending),
        }
        if time_unit is not None:
            if not isinstance(time_unit, str):
                raise TypeError("the time unit must be a string")
            self._document["time_unit"] = time_unit

    def at(
        self, time, *, parameters=None, change=None, set_counts=None, partition=None
    ):
        """Add an atomic group of rate edits and count operations.

        Count increments may be negative. Counts assigned or partitioned in
        the same group use their pre-group values. Give at most one count
        operation per species. Return the protocol for chaining.
        """
        item = {"time": float(time)}
        for key, value, convert in (
            ("parameters", parameters, float),
            ("change", change, index),
            ("set", set_counts, index),
            ("partition", partition, float),
        ):
            if value is not None:
                item[key] = {k: convert(v) for k, v in value.items()}
        self._document["actions"].append(item)
        return self

    def delay(self, reaction, duration, *, products=None, owner=None):
        """Add fixed delayed products to a named initiation reaction.

        When supplied, the owner is a busy species created once per
        initiation. Its loss cancels the completion; a completion consumes
        one surviving owner.
        """
        item = {"reaction": reaction, "duration": float(duration)}
        if products is not None:
            item["products"] = {k: index(v) for k, v in products.items()}
        if owner is not None:
            item["owner"] = owner
        self._document["delays"].append(item)
        return self

    def to_dict(self):
        """Return a detached copy of the protocol payload for run/1.1."""
        return deepcopy(self._document)
