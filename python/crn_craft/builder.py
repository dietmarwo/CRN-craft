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

A reaction with up to two reactant molecules is written as in release
0.6.0. Three or more molecules, such as ``["A", "A", "A"]``, give one
pattern per species with the pattern field ``count`` (``3A -> ...`` fires at
k·C(A, 3)); ``count`` exists only in ``crn-craft.model/1.0``, so such a
model writes that version. ``add_ode`` (also ``add_ODE``; ``add_continuous``
is its deprecated name) declares a continuous variable, an ODE between
events that string rates may read; such models run with ``method="hybrid"``.
``set_stop`` sets the model's stop condition, for example
``sir.set_stop("I == 0")``. Species with attributes, observables with filters
and the other model fields remain available by editing ``to_dict()`` or
writing JSON; :meth:`Model.from_dict` wraps any complete model document.

``crn_craft.Model`` (the top-level name) is this builder writing
``crn-craft.model/1.0`` always; ``crn_craft.rules.Model`` writes the 0.6.0
versions unless the model needs 1.0.
"""

from collections.abc import Mapping
import copy
import json
import math
import numbers
import operator
import re
import warnings

SCHEMA = "crn-craft.rules-model/0.3"
# Written instead of SCHEMA when a continuous variable is independent.
SCHEMA_INDEPENDENT = "crn-craft.rules-model/0.4"
# Written when the model uses a 1.0 feature (a pattern `count`, more than two
# patterns), and always by the top-level `crn_craft.Model`.
SCHEMA_1_0 = "crn-craft.model/1.0"
# Versions whose documents `Model.from_dict` can hand back to the builder.
_EDITABLE = (SCHEMA, SCHEMA_INDEPENDENT, SCHEMA_1_0)
_RULES_MODEL = "crn-craft.rules-model/"
_DOCUMENT_KEYS = frozenset(
    "schema name description parameters types initial observables rules stop "
    "continuous".split()
)

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


def _patterns(reactants):
    """The reactant patterns of a list of species names.

    Up to two molecules: one pattern per molecule, labeled ``x`` and ``y``
    (the encoding of release 0.6.0, also for ``["A", "A"]``). Three or
    more: one pattern per species in first-appearance order with its
    multiplicity as ``count`` (left out when 1), labeled ``x`` and ``y``
    for up to two species and ``x1``, ``x2``, ... for more, as
    ``crn_craft.convert`` labels network reactions.
    """
    if len(reactants) <= 2:
        return [
            {"label": label, "type": species}
            for label, species in zip(_LABELS[: len(reactants)], reactants, strict=True)
        ]
    counts = {}
    for species in reactants:
        counts[species] = counts.get(species, 0) + 1
    many = len(counts) > 2
    patterns = []
    for i, (species, count) in enumerate(counts.items()):
        pattern = {"label": f"x{i + 1}" if many else _LABELS[i], "type": species}
        if count > 1:
            pattern["count"] = count
        patterns.append(pattern)
    return patterns


def _needs_1_0(rule):
    """Whether a rule uses a feature of crn-craft.model/1.0 only."""
    reactants = rule.get("reactants")
    if not isinstance(reactants, list):
        return False
    return len(reactants) > 2 or any(
        isinstance(p, Mapping) and "count" in p for p in reactants
    )


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
    returns a fresh document: rules-model/0.3 (0.4 if an ODE is
    independent, ``crn-craft.model/1.0`` if a reaction takes three or more
    molecules), which ``rules.ensemble`` and ``rules.prepare`` also accept
    as the ``Model`` itself. ``crn_craft.Model``, the top-level name, writes
    ``crn-craft.model/1.0`` always.

    :meth:`from_dict` wraps a complete model document, for example a
    converted network model; :meth:`fsp`, :meth:`analyze` and
    :meth:`to_gillespy2` work on any model.
    """

    # The one flag of the unified API: True writes crn-craft.model/1.0. The
    # top-level `crn_craft.Model` sets it for the class; `from_dict` sets it
    # on an instance whose document is model/1.0.
    _unified = False

    def __init__(self, name="model"):
        self._name = _check_name("model", name)
        self._parameters = {}
        self._species = {}
        # Declared by a reaction before add_species; add_species may set them.
        self._implicit = set()
        self._reactions = []
        self._continuous = []
        self._stop = None
        # A document's description, kept by `from_dict`.
        self._description = None
        # From `from_dict`: the document, returned as it is until an edit.
        self._document = None
        # True when the builder cannot represent the document, so the model
        # cannot be edited.
        self._fixed = False

    @classmethod
    def from_dict(cls, document):
        """A model of a complete model document.

        ``document`` is a dictionary, JSON text or bytes: any rules-model
        version, ``crn-craft.model/1.0`` or a network ``crn-craft.model/0.1``,
        which is converted first (:func:`crn_craft.convert`). The engine
        checks it (``CrnError`` with its diagnostics if invalid).
        :meth:`to_dict` returns the document as given, key order included,
        so its hash is that of the document; ``crn_craft.Model.from_dict``
        writes it as ``crn-craft.model/1.0``, changing only the schema
        string of a rules-model.

        A document the builder could have written (attribute-free types,
        one count observable per type named after it, one numeric initial
        count per type; rules-model/0.3, 0.4 or model/1.0), which includes
        every converted network model and the packaged examples, stays
        editable: the ``add_*`` and ``set_*`` methods work, and after the
        first edit ``to_dict`` writes the builder's document (explicit
        default values such as ``"attributes": []`` and ``"reduce": "sum"``
        are then left out, which changes the bytes but not the meaning).
        Any other document (attributes, filtered or weighted observables, an
        ATP block, ...) is kept whole and cannot be edited: those methods
        raise ``ValueError``; edit the dictionary and call ``from_dict``
        again.
        """
        from ._common import _native, convert, model_bytes

        if isinstance(document, bytes | str):
            document = json.loads(document)
        if not isinstance(document, Mapping):
            raise TypeError(
                f"a model document is a dictionary, JSON text or bytes, not "
                f"{type(document).__name__}"
            )
        if document.get("schema") == "crn-craft.model/0.1":
            document = json.loads(convert(document))
        document = copy.deepcopy(dict(document))
        check = json.loads(_native.check_rules_json(model_bytes(document)))
        if not check["valid"]:
            from ._common import CrnError

            raise CrnError(check["diagnostics"])
        model = cls(document["name"])
        if document["schema"] == SCHEMA_1_0:
            model._unified = True
        if not model._load(document):
            model = cls(document["name"])
            if document["schema"] == SCHEMA_1_0:
                model._unified = True
            model._fixed = True
        model._document = document
        return model

    def _load(self, document):
        """Take over a document's content; True if the builder can write it.

        The builder writes attribute-free types, one numeric initial entry
        per type, one count observable per type named after it, and keeps
        rules and continuous variables as given. Explicit default values
        (``"attributes": []`` on a type, ``"attributes": {}`` on an initial
        entry, ``"reduce": "sum"``, an integral count written as a float),
        which the converter writes, are accepted.
        """
        if document.get("schema") not in _EDITABLE or not set(document) <= (
            _DOCUMENT_KEYS
        ):
            return False
        try:
            types = []
            for entry in document["types"]:
                if set(entry) - {"name", "attributes"} or entry.get("attributes", []):
                    return False
                types.append(entry["name"])
            counts = {}
            for entry in document["initial"]:
                count = entry["count"]
                if (
                    set(entry) - {"type", "attributes", "count"}
                    or entry.get("attributes", {})
                    or entry["type"] in counts
                    or entry["type"] not in types
                    or isinstance(count, str | bool)
                ):
                    return False
                if isinstance(count, float) and count.is_integer():
                    count = int(count)
                counts[entry["type"]] = _count(count, "an initial count")
            observables = document["observables"]
            if set(observables) != set(types) or any(
                set(o) - {"type", "reduce"}
                or o["type"] != name
                or o.get("reduce", "sum") != "sum"
                for name, o in observables.items()
            ):
                return False
            self._species = {name: counts.get(name, 0) for name in types}
            self._parameters = {
                key: _number(value, key)
                for key, value in document["parameters"].items()
            }
            self._reactions = copy.deepcopy(list(document["rules"]))
            self._continuous = copy.deepcopy(list(document.get("continuous", [])))
            self._stop = document.get("stop")
            self._description = document.get("description")
        except (KeyError, TypeError, ValueError, AttributeError):
            return False
        return True

    def _editable(self):
        if self._fixed:
            raise ValueError(
                f"model {self._name!r} holds a document the builder cannot "
                "write (attributes, filtered observables or other fields); edit "
                "the dictionary from to_dict() and call Model.from_dict instead"
            )

    def _edited(self):
        self._document = None

    @property
    def name(self):
        """The model's name."""
        return self._name

    @property
    def species(self):
        """Species (type) names in declaration order."""
        if self._fixed:
            return tuple(t["name"] for t in self._document["types"])
        return tuple(self._species)

    @property
    def parameters(self):
        """A copy of the parameter values by name."""
        if self._fixed:
            return dict(self._document.get("parameters", {}))
        return dict(self._parameters)

    @property
    def reaction_names(self):
        """Reaction (rule) names in declaration order."""
        if self._fixed:
            return tuple(r["name"] for r in self._document["rules"])
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
        self._editable()
        entries = _entries(name, value, "add_parameter", "value")
        checked = {}
        for key, number in entries:
            _check_name("parameter", key)
            if self._taken(key):
                raise ValueError(f"name {key!r} is already declared")
            checked[key] = _number(number, f"parameter {key}")
        self._parameters.update(checked)
        self._edited()
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
        self._editable()
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
        self._edited()
        return self

    def set_initial(self, counts=None, /, **kwargs):
        """Set initial counts of declared species, from a mapping or keywords.

        ``model.set_initial({"S": 990, "I": 10})`` and
        ``model.set_initial(S=990, I=10)`` are equivalent.
        """
        self._editable()
        values = dict(counts or {})
        values.update(kwargs)
        for name in values:
            if name not in self._species:
                raise KeyError(f"unknown species {name!r}; declare it first")
        for name, value in values.items():
            self._species[name] = _count(value, f"initial count of {name}")
        self._edited()
        return self

    def add_reaction(self, rate, reactants, products, name=None):
        """Add a reaction ``reactants -> products``.

        ``rate`` is a number (a mass-action constant k: the reaction fires
        at k times the number of distinct combinations of reactant
        molecules, so ``["A", "A"]`` fires at k·A(A−1)/2) or a string (the
        whole propensity). Species are lists of names, repeated for
        multiplicity, for example ``["A", "A"]``. Three or more reactant
        molecules, such as ``["A", "A", "A"]`` or ``["A", "A", "B"]``, are
        written with the pattern field ``count`` of ``crn-craft.model/1.0``
        (a numeric k fires at k·C(A, 2)·B for ``["A", "A", "B"]``; a string
        propensity applies while each species has its count). Species not
        yet declared are declared with count 0, as in rebop. ``name``
        defaults to ``r<index>``.
        """
        self._editable()
        if isinstance(reactants, str) or isinstance(products, str):
            raise TypeError("reactants and products are lists of species names")
        reactants = [_check_name("species", s) for s in reactants]
        products = [_check_name("species", s) for s in products]
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
        rule = {"name": name, "reactants": _patterns(reactants)}
        rule.update(law)
        rule["products"] = [{"type": species} for species in products]
        self._reactions.append(rule)
        self._edited()
        return self

    def set_stop(self, condition):
        """Set the stop condition, or remove it with ``None``.

        ``condition`` is an expression over parameters and species counts,
        for example ``"I == 0"`` or ``"N > 0 && mutants >= 0.9 * N"``. Once
        it is true (it is checked at time 0 and after every event), a
        replicate keeps its state to the end of the report grid and is
        recorded as ``"stopped"``, with its ``stop_time``. It cannot read
        continuous variables; a bound of ``add_ode`` stops on those.
        """
        self._editable()
        if condition is not None:
            if not isinstance(condition, str):
                raise TypeError(
                    "a stop condition is an expression string or None, "
                    f"not {type(condition).__name__}"
                )
            if not condition.strip():
                raise ValueError("a stop condition must be an expression")
        self._stop = condition
        self._edited()
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
        self._editable()
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
        self._edited()
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

    def _written(self):
        """The document the builder writes from its state."""
        independent = any(c.get("independent") for c in self._continuous)
        if self._unified or any(_needs_1_0(r) for r in self._reactions):
            schema = SCHEMA_1_0
        else:
            schema = SCHEMA_INDEPENDENT if independent else SCHEMA
        document = {
            "schema": schema,
            "name": self._name,
        }
        if self._description is not None:
            document["description"] = self._description
        document |= {
            "parameters": dict(self._parameters),
            "types": [{"name": s} for s in self._species],
            "initial": [
                {"type": s, "count": n} for s, n in self._species.items() if n > 0
            ],
            "observables": {s: {"type": s} for s in self._species},
            "rules": json.loads(json.dumps(self._reactions)),
        }
        if self._stop is not None:
            document["stop"] = self._stop
        if self._continuous:
            document["continuous"] = json.loads(json.dumps(self._continuous))
        return document

    def to_dict(self):
        """Return the document as a new dictionary.

        It is rules-model/0.3, or rules-model/0.4 if a continuous variable
        is independent (see :meth:`add_ode`), or ``crn-craft.model/1.0`` if
        a reaction takes three or more molecules; ``crn_craft.Model`` writes
        ``crn-craft.model/1.0`` always. A model from :meth:`from_dict`
        returns its document until it is edited.
        """
        if self._document is not None:
            document = copy.deepcopy(self._document)
            schema = document.get("schema")
            if (
                self._unified
                and isinstance(schema, str)
                and schema.startswith(_RULES_MODEL)
            ):
                # In place, so the key order (and the bytes) stay the same.
                document["schema"] = SCHEMA_1_0
            return document
        return self._written()

    def to_json(self, *, indent=2):
        """Return the document as JSON text with a trailing newline."""
        return json.dumps(self.to_dict(), indent=indent, allow_nan=False) + "\n"

    def fsp(self, end, report_times, bounds, *, max_states=None, tolerance=None):
        """The exact distribution of the counts by finite state projection.

        ``bounds`` maps every species label to its largest count (a type
        name for a type without attributes, for example ``{"S": 999, "I":
        1000, "R": 1000}``); ``end`` and ``report_times`` are as in
        ``ensemble``. ``max_states`` (default 100,000) caps the enumerated
        states and ``tolerance`` (default ``1e-12``) is the truncation
        budget of the series. The model must be closed (``CrnError``
        ``CRN-R017`` otherwise). Returns an :class:`FspResult`; see
        ``crn_craft.rules.fsp`` for the method and its error terms.
        """
        from .rules import FspResult, fsp

        return FspResult(
            fsp(
                self.to_dict(),
                end,
                report_times,
                bounds,
                max_states=max_states,
                tolerance=tolerance,
            )
        )

    def analyze(self):
        """The network analysis of this model (``crn-craft.analysis/1.0``).

        Returns the report as a dictionary, as ``crn-craft analyze MODEL``
        prints it: species and reactions of the lowered network, their
        dependencies, conservation laws, the cone of reachable directions,
        the complex graph and warnings. Read ``success`` first: a model
        without a flat view (continuous variables, ATP, products that
        create new attribute tuples) gives ``success`` false with
        ``CRN-R017``.
        """
        from ._common import _native, model_bytes

        return json.loads(_native.analyze_json(model_bytes(self.to_dict())))

    def to_gillespy2(self):
        """This model as a ``gillespy2.Model``.

        GillesPy2 must be installed, in the version the export names
        (``crn-craft describe`` lists it; another version raises
        ``RuntimeError``). A rule with a constant ``rate`` and at most two molecules becomes a
        mass-action reaction, any other rule a custom propensity, as
        ``crn-craft export MODEL gillespy2`` writes them. A model the export
        refuses raises ``CrnError`` (``CRN-E050`` for what GillesPy2 cannot
        evaluate, such as a stop condition or ``exp``).
        """
        from ._common import CrnError, _native, model_bytes

        report = json.loads(_native.export_json(model_bytes(self.to_dict())))
        if not report["success"]:
            raise CrnError(report["diagnostics"])
        namespace = {}
        # The engine's generated module; it imports gillespy2 and checks
        # its version when build_model() runs.
        exec(
            compile(report["artifact"]["source"], "<crn-craft export>", "exec"),
            namespace,
        )
        return namespace["build_model"]()

    def __repr__(self):
        return (
            f"Model({self._name!r}: {len(self.species)} species, "
            f"{len(self.reaction_names)} reactions)"
        )
