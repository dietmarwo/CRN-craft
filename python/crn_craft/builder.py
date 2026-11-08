"""Build rules models from Python, in the style of rebop's ``Gillespie``.

``Model`` writes a ``crn-craft.rules-model/0.3`` document (0.4 with
independent ODEs, see ``add_ode``) with one attribute-free species type per
species, one observable per species (its count, under the species name) and
one rule per reaction::

    from crn_craft import rules

    sir = rules.Model("sir")
    sir.add_species({"S": 990, "I": 10, "R": 0})
    sir.add_reaction(1e-4, ["S", "I"], ["I", "I"], name="infection")
    sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    sir.add_parameter("gamma", 0.05)
    run = rules.ensemble(sir, end=100.0, report_times=range(101), replicates=1000)

``add_species`` and ``add_parameter`` (also named ``add_parameters``) take
one name and value, or a mapping of several.

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
language. ``add_ode`` (also ``add_ODE``; ``add_continuous`` is its
deprecated name) declares a continuous variable, an ODE between events that
string rates may read; such models run with ``method="hybrid"``.
Species with attributes, observables with filters and the other model fields
remain available by editing ``to_dict()`` or writing JSON.
"""

from collections.abc import Mapping
import json
import math
import numbers
import operator
import re
import warnings

SCHEMA = "crn-craft.rules-model/0.3"
# Written instead of SCHEMA when a continuous variable is independent.
SCHEMA_INDEPENDENT = "crn-craft.rules-model/0.4"

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


class _Default:
    """A default argument value that tells "not given" from any given value."""

    def __init__(self, text):
        self._text = text

    def __repr__(self):
        return self._text


# `add_parameter`'s value has no default; `add_species`'s count defaults to 0.
_REQUIRED = _Default("<required>")
_ZERO = _Default("0")


def _entries(name, value, method, argument):
    """The (name, value) pairs of a one-entry call or of a mapping."""
    if isinstance(name, Mapping):
        if not isinstance(value, _Default):
            raise TypeError(
                f"{method}() takes a mapping or a name and its {argument}, not both"
            )
        return list(name.items())
    if value is _REQUIRED:
        raise TypeError(f"{method}() is missing the {argument} of {name!r}")
    return [(name, 0 if value is _ZERO else value)]


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
    returns a fresh rules-model/0.3 document (0.4 if an ODE is
    independent), which ``rules.ensemble`` and ``rules.prepare`` also
    accept as the ``Model`` itself.
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

    def add_parameter(self, name, value=_REQUIRED):
        """Declare parameters that string rates may read.

        ``add_parameter("gamma", 0.01)`` declares one;
        ``add_parameter({"gamma": 0.01, "beta": 2e-4})`` declares every
        entry of a mapping, in its order, as one call each would. A mapping
        takes no other argument (``TypeError``). If an entry is invalid,
        nothing is declared. ``add_parameters`` is the same method.
        """
        entries = _entries(name, value, "add_parameter", "value")
        checked = {}
        for key, number in entries:
            _check_name("parameter", key)
            if self._taken(key):
                raise ValueError(f"name {key!r} is already declared")
            checked[key] = _number(number, f"parameter {key}")
        self._parameters.update(checked)
        return self

    add_parameters = add_parameter

    def add_species(self, name, initial=_ZERO):
        """Declare species with initial counts (default 0).

        ``add_species("S", 999)`` declares one; ``add_species({"S": 999,
        "I": 1})`` declares every entry of a mapping, in its order, as one
        call each would. A mapping takes no other argument (``TypeError``).
        If an entry is invalid, nothing is declared. A species that a
        reaction already introduced keeps its position and takes the count.
        """
        entries = _entries(name, initial, "add_species", "initial")
        checked = {}
        for key, count in entries:
            _check_name("species", key)
            if self._taken(key) and key not in self._implicit:
                raise ValueError(f"name {key!r} is already declared")
            checked[key] = _count(count, f"initial count of {key}")
        for key, count in checked.items():
            self._species[key] = count
            self._implicit.discard(key)
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

    def add_ode(
        self,
        name,
        initial,
        derivative,
        *,
        independent=False,
        lower_bound=None,
        on_bound="stop",
    ):
        """Declare a continuous variable dx/dt = ``derivative`` between events.

        ``initial`` and ``lower_bound`` are numbers or parameter expressions;
        ``derivative`` is an expression over parameters, species counts and
        continuous variables. String rates may read the variable. With a
        ``lower_bound``, ``on_bound`` is ``"stop"`` (end the replicate when
        the variable reaches the bound) or ``"clamp"`` (hold it there). Run
        such a model with ``method="hybrid"``.

        ``independent=True`` declares that ``initial`` and ``derivative``
        read only parameters and other independent variables, never species
        counts, so the ODE can be solved once for the whole run instead of
        per replicate. It writes ``"independent": true`` on the variable and
        makes ``to_dict`` write ``crn-craft.rules-model/0.4``; the engine
        refuses an independent variable that reads anything else. Without
        independent variables the document is rules-model/0.3, as before.
        ``add_ODE`` is the same method.
        """
        if not isinstance(independent, bool):
            raise TypeError(
                f"independent must be True or False, not {type(independent).__name__}"
            )
        _check_name("variable", name)
        if self._taken(name):
            raise ValueError(f"name {name!r} is already declared")
        if not isinstance(derivative, str) or not derivative.strip():
            raise TypeError("derivative must be an expression string")
        variable = {
            "name": name,
            "initial": initial if isinstance(initial, str) else _number(initial, name),
            "derivative": derivative,
        }
        if independent:
            variable["independent"] = True
        if lower_bound is not None:
            if on_bound not in ("stop", "clamp"):
                raise ValueError(
                    f"on_bound must be 'stop' or 'clamp', not {on_bound!r}"
                )
            variable["lower_bound"] = (
                lower_bound
                if isinstance(lower_bound, str)
                else _number(lower_bound, f"lower bound of {name}")
            )
            variable["on_bound"] = on_bound
        self._continuous.append(variable)
        return self

    add_ODE = add_ode

    def add_continuous(
        self, name, initial, derivative, *, lower_bound=None, on_bound="stop"
    ):
        """Deprecated since 0.5.0: use :meth:`add_ode`, which takes the same
        arguments and ``independent=``. Writes the same document as before.
        """
        warnings.warn(
            "Model.add_continuous is deprecated; use Model.add_ode (or add_ODE), "
            "which takes the same arguments",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.add_ode(
            name, initial, derivative, lower_bound=lower_bound, on_bound=on_bound
        )

    def to_dict(self):
        """Return the document as a new dictionary.

        It is rules-model/0.3, or rules-model/0.4 if a continuous variable
        is independent (see :meth:`add_ode`).
        """
        independent = any(c.get("independent") for c in self._continuous)
        document = {
            "schema": SCHEMA_INDEPENDENT if independent else SCHEMA,
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
