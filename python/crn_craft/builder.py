"""Build rules models from Python, in the style of rebop's ``Gillespie``.

``Model`` writes a ``crn-craft.rules-model/0.3`` document with one
attribute-free species type per species, one observable per species (its
count, under the species name) and one rule per reaction::

    from crn_craft import rules

    sir = rules.Model("sir")
    sir.add_species("S", 990).add_species("I", 10).add_species("R")
    sir.add_reaction(1e-4, ["S", "I"], ["I", "I"], name="infection")
    sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    sir.add_parameter("gamma", 0.05)
    run = rules.ensemble(sir, end=100.0, report_times=range(101), replicates=1000)

A numeric rate is a mass-action constant k, counted like the rest of
CRN-craft: k times the number of distinct combinations of reactant molecules,
so ``A`` fires at k·A, ``A + B`` at k·A·B and ``A + A`` at k·A(A−1)/2. rebop
leaves out the 1/2 for identical reactants: its numeric k on ``["A", "A"]``
equals ``2 * k`` here. A string rate is the reaction's whole propensity, an
expression over parameters and species counts in the rules expression
language (``+ - * / ^``, comparisons, ``min max exp log sqrt abs floor ceil
if``), as in rebop. It applies while each reactant has enough molecules and
is 0 otherwise.

A reaction takes at most two reactant molecules, the limit of the rules
language. Species with attributes, observables with filters and the other
model fields remain available by editing ``to_dict()`` or writing JSON.
"""

import json
import math
import numbers
import operator
import re

SCHEMA = "crn-craft.rules-model/0.3"

_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_FUNCTIONS = frozenset(
    "min max exp log sqrt abs floor ceil if bernoulli choice uniform normal new_id".split()
)
_LABELS = ("x", "y")


def _check_name(kind, name):
    if not isinstance(name, str):
        raise TypeError(f"{kind} name must be a string, not {type(name).__name__}")
    if not _NAME.match(name) or len(name) > 64 or name in _FUNCTIONS:
        raise ValueError(
            f"invalid {kind} name {name!r}: use an ASCII identifier of at most "
            "64 bytes that is not a function name"
        )
    return name


def _number(value, what):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{what} must be a number, not {type(value).__name__}")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{what} must be finite, not {value}")
    return value


def _count(value, what):
    try:
        value = operator.index(value)
    except TypeError:
        raise TypeError(
            f"{what} must be an integer, not {type(value).__name__}"
        ) from None
    if not 0 <= value <= 1 << 53:
        raise ValueError(f"{what} must be in [0, 2**53], not {value}")
    return value


class Model:
    """A reaction network over attribute-free species.

    ``add_*`` methods return the model, so calls can be chained. ``to_dict``
    returns a fresh rules-model/0.3 document, which ``rules.ensemble`` and
    ``rules.prepare`` also accept as the ``Model`` itself.
    """

    def __init__(self, name="model"):
        self._name = _check_name("model", name)
        self._parameters = {}
        self._species = {}
        # Declared by a reaction before add_species; add_species may set them.
        self._implicit = set()
        self._reactions = []
        self._continuous = []

    @property
    def name(self):
        """The model's name."""
        return self._name

    @property
    def species(self):
        """Species names in declaration order."""
        return tuple(self._species)

    @property
    def parameters(self):
        """A copy of the parameter values by name."""
        return dict(self._parameters)

    @property
    def reaction_names(self):
        """Reaction (rule) names in declaration order."""
        return tuple(r["name"] for r in self._reactions)

    def _taken(self, name):
        return (
            name in self._parameters
            or name in self._species
            or any(c["name"] == name for c in self._continuous)
        )

    def add_parameter(self, name, value):
        """Declare a parameter that string rates may read."""
        _check_name("parameter", name)
        if self._taken(name):
            raise ValueError(f"name {name!r} is already declared")
        self._parameters[name] = _number(value, f"parameter {name}")
        return self

    def add_species(self, name, initial=0):
        """Declare a species with an initial count (default 0).

        A species that a reaction already introduced keeps its position and
        takes the count.
        """
        _check_name("species", name)
        if self._taken(name) and name not in self._implicit:
            raise ValueError(f"name {name!r} is already declared")
        self._species[name] = _count(initial, f"initial count of {name}")
        self._implicit.discard(name)
        return self

    def set_initial(self, counts=None, /, **kwargs):
        """Set initial counts of declared species, from a mapping or keywords.

        ``model.set_initial({"S": 990, "I": 10})`` and
        ``model.set_initial(S=990, I=10)`` are equivalent.
        """
        values = dict(counts or {})
        values.update(kwargs)
        for name in values:
            if name not in self._species:
                raise KeyError(f"unknown species {name!r}; declare it first")
        for name, value in values.items():
            self._species[name] = _count(value, f"initial count of {name}")
        return self

    def add_reaction(self, rate, reactants, products, name=None):
        """Add a reaction ``reactants -> products``.

        ``rate`` is a number (a mass-action constant, see the module
        documentation for the counting of identical reactants) or a string
        (the whole propensity). Species are lists of names, repeated for
        multiplicity, for example ``["A", "A"]``. Species not yet declared are
        declared with count 0, as in rebop. ``name`` defaults to ``r<index>``.
        """
        if isinstance(reactants, str) or isinstance(products, str):
            raise TypeError("reactants and products are lists of species names")
        reactants = [_check_name("species", s) for s in reactants]
        products = [_check_name("species", s) for s in products]
        if len(reactants) > 2:
            raise ValueError(
                "a reaction takes at most two reactant molecules in the rules "
                f"language, not {len(reactants)}; split it into steps"
            )
        if isinstance(rate, str):
            if not rate.strip():
                raise ValueError("a string rate must be a propensity expression")
            law = {"propensity": rate}
        else:
            k = _number(rate, "a numeric rate")
            if k < 0:
                raise ValueError(f"a numeric rate must be nonnegative, not {k}")
            law = {"rate": repr(k)}
        name = f"r{len(self._reactions)}" if name is None else name
        _check_name("reaction", name)
        if name in self.reaction_names:
            raise ValueError(f"reaction name {name!r} is already used")
        for species in reactants + products:
            if species not in self._species:
                if self._taken(species):
                    raise ValueError(f"{species!r} is a parameter or variable")
                self._species[species] = 0
                self._implicit.add(species)
        rule = {
            "name": name,
            "reactants": [
                {"label": label, "type": species}
                for label, species in zip(
                    _LABELS[: len(reactants)], reactants, strict=True
                )
            ],
        }
        rule.update(law)
        rule["products"] = [{"type": species} for species in products]
        self._reactions.append(rule)
        return self

    def to_dict(self):
        """Return the rules-model/0.3 document as a new dictionary."""
        document = {
            "schema": SCHEMA,
            "name": self._name,
            "parameters": dict(self._parameters),
            "types": [{"name": s} for s in self._species],
            "initial": [
                {"type": s, "count": n} for s, n in self._species.items() if n > 0
            ],
            "observables": {s: {"type": s} for s in self._species},
            "rules": json.loads(json.dumps(self._reactions)),
        }
        if self._continuous:
            document["continuous"] = json.loads(json.dumps(self._continuous))
        return document

    def to_json(self, *, indent=2):
        """Return the document as JSON text with a trailing newline."""
        return json.dumps(self.to_dict(), indent=indent, allow_nan=False) + "\n"

    def __repr__(self):
        return (
            f"Model({self._name!r}: {len(self._species)} species, "
            f"{len(self._reactions)} reactions)"
        )
