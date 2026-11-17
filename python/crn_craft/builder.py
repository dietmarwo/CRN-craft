"""Build models from Python, in the style of rebop's ``Gillespie``.

``Model`` writes a ``crn-craft.rules-model/0.3`` document (0.4 with
independent ODEs, see ``add_ode``; ``crn-craft.model/1.0`` when the model
uses a feature of that version) with one species type per species, one
observable per species (its count, under the species name), the
observables of ``add_observable`` and one rule per reaction::

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
is 0 otherwise. ``mass_action(expression)`` is a rate per combination of
reactant molecules, like a number, given as an expression.

A reaction with up to two reactant molecules is written as in release
0.6.0. Three or more molecules, such as ``["A", "A", "A"]``, give one
pattern per species with the pattern field ``count`` (``3A -> ...`` fires at
k·C(A, 3)); ``count`` exists only in ``crn-craft.model/1.0``, so such a
model writes that version. ``add_ode`` (also ``add_ODE``; ``add_continuous``
is its deprecated name) declares a continuous variable, an ODE between
events that string rates may read; such models run with ``method="hybrid"``.
``set_stop`` sets the model's stop condition, for example
``sir.set_stop("I == 0")``.

A species may carry numeric attributes, for example a clone identity:
``add_species("mtDNA", attributes=["clone"])``. Each attribute vector is
then an instance with its own count, created and removed during the run.
``add_initial("mtDNA", 1000, clone=0)`` gives the molecules at time 0;
reactions match molecules with a species name or :func:`pattern` (with a
``where`` condition or a ``factor``) and create them with :func:`copy`
(the attributes of a reactant) and :func:`new` (every attribute given);
``add_observable`` declares filtered, weighted and reduced observables::

    clones = rules.Model("clones")
    clones.add_parameter({"delta": 0.1, "pmut": 0.01})
    clones.add_species("mtDNA", attributes=["clone"])
    clones.add_initial("mtDNA", 100, clone=0)
    clones.add_observable("mutants", "mtDNA", where="x.clone > 0")
    clones.add_observable("types", "mtDNA", reduce="distinct")
    clones.add_reaction(mass_action("delta"), ["mtDNA"], [], name="degradation")
    clones.add_reaction(mass_action("delta * pmut"), ["mtDNA"],
                        [copy("x"), new("mtDNA", clone="new_id()")], name="mutation")

:meth:`Model.from_dict` wraps any complete model document; documents the
builder can write stay editable.

``crn_craft.Model`` (the top-level name) is this builder writing
``crn-craft.model/1.0`` always; ``crn_craft.rules.Model`` writes the 0.6.0
versions unless the model needs 1.0.
"""

from collections.abc import Iterable, Mapping
from copy import deepcopy
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
# patterns, a `max_count` observable), and always by the top-level
# `crn_craft.Model`.
SCHEMA_1_0 = "crn-craft.model/1.0"
# Versions whose documents `Model.from_dict` can hand back to the builder.
_EDITABLE = (SCHEMA, SCHEMA_INDEPENDENT, SCHEMA_1_0)
# Also editable for `crn_craft.Model`, which writes them as model/1.0, a
# version that means what they mean.
_EDITABLE_UNIFIED = (
    "crn-craft.rules-model/0.1",
    "crn-craft.rules-model/0.2",
)
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
# The reductions of an observable; `max_count` exists only in model/1.0.
REDUCTIONS = ("sum", "distinct", "max", "min", "max_count")
_REDUCTIONS_1_0 = frozenset({"max_count"})
_OBSERVABLE_KEYS = frozenset({"type", "label", "where", "weight", "reduce"})
_PATTERN_KEYS = frozenset({"label", "type", "where", "factor", "count"})
_PRODUCT_KEYS = frozenset({"type", "copy", "attributes"})
# The engine's limit of attributes per type.
_MAX_ATTRIBUTES = 64


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


def _value(value, what):
    """A document value: an integer, a finite float or an expression string."""
    if isinstance(value, str):
        if not value.strip():
            raise ValueError(f"{what} must be a number or an expression, not ''")
        return value
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(
            f"{what} must be a number or an expression string, not "
            f"{type(value).__name__}"
        )
    if isinstance(value, numbers.Integral):
        return int(value)
    return _number(value, what)


def _expression(value, what):
    """An expression string; a number is written as its shortest repr."""
    if isinstance(value, str):
        if not value.strip():
            raise ValueError(f"{what} must be an expression, not ''")
        return value
    return repr(_number(value, what))


def _condition(value, what):
    """None, or a nonempty condition string."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(
            f"{what} must be an expression string, not {type(value).__name__}"
        )
    if not value.strip():
        raise ValueError(f"{what} must be an expression, not ''")
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


def _attribute_names(attributes):
    """The attribute names of ``add_species``: a tuple, empty for None."""
    if attributes is None:
        return ()
    if isinstance(attributes, str | bytes | Mapping) or not isinstance(
        attributes, Iterable
    ):
        raise TypeError(
            "attributes is a list of attribute names, for example ['clone'], "
            f"not {type(attributes).__name__}"
        )
    names = tuple(attributes)
    for name in names:
        _check_name("attribute", name)
    for i, name in enumerate(names):
        if name in names[:i]:
            raise ValueError(f"attribute {name!r} is listed twice")
    if len(names) > _MAX_ATTRIBUTES:
        raise ValueError(
            f"a species has at most {_MAX_ATTRIBUTES} attributes, not {len(names)}"
        )
    return names


class _Pattern:
    """A reactant of :func:`pattern`."""

    __slots__ = ("species", "where", "factor", "label")

    def __init__(self, species, where, factor, label):
        self.species = species
        self.where = where
        self.factor = factor
        self.label = label

    def __repr__(self):
        given = [repr(self.species)]
        given += [
            f"{key}={getattr(self, key)!r}"
            for key in ("where", "factor", "label")
            if getattr(self, key) is not None
        ]
        return f"pattern({', '.join(given)})"


class _Copy:
    """A product of :func:`copy`."""

    __slots__ = ("label", "attributes")

    def __init__(self, label, attributes):
        self.label = label
        self.attributes = attributes

    def __repr__(self):
        given = [repr(self.label)]
        given += [f"{key}={value!r}" for key, value in self.attributes.items()]
        return f"copy({', '.join(given)})"


class _New:
    """A product of :func:`new`."""

    __slots__ = ("species", "attributes")

    def __init__(self, species, attributes):
        self.species = species
        self.attributes = attributes

    def __repr__(self):
        given = [repr(self.species)]
        given += [f"{key}={value!r}" for key, value in self.attributes.items()]
        return f"new({', '.join(given)})"


class _MassAction:
    """A rate of :func:`mass_action`."""

    __slots__ = ("expression",)

    def __init__(self, expression):
        self.expression = expression

    def __repr__(self):
        return f"mass_action({self.expression!r})"


def pattern(species, where=None, factor=None, label=None):
    """A reactant of ``add_reaction``: the molecules of ``species`` that match.

    ``where`` is a condition over parameters and the molecule's attributes,
    for example ``"x.clone > 0"``; ``factor`` is a weight h ≥ 0 per
    molecule over the same names (a number or an expression), so the
    reaction fires at its rate times the sum of h over the matching
    molecules. ``label`` names the molecule in ``where``, ``factor``,
    product attributes and :func:`copy`; it defaults to ``x`` for the first
    reactant and ``y`` for the second, so a condition on the second
    reactant reads ``y.clone`` (or give ``label``). A plain species name is
    the pattern without ``where`` and ``factor``. A reaction with a
    ``where`` or ``factor`` takes a rate per combination of molecules, a
    number or :func:`mass_action`, not a whole-propensity string. The
    model document writes it as ``{"label", "type", "where", "factor"}``;
    such a dictionary is accepted as a reactant too.
    """
    _check_name("species", species)
    where = _condition(where, "a pattern's where")
    if factor is not None:
        factor = _expression(factor, "a pattern's factor")
    if label is not None:
        _check_name("label", label)
    return _Pattern(species, where, factor, label)


def copy(label, /, **attributes):
    """A product of ``add_reaction``: a copy of the reactant labeled ``label``.

    The product has the reactant's species and attribute values, except the
    attributes given as keywords, which it overrides: for example
    ``copy("x", size="x.size + 1")``. Values are numbers or expressions over
    parameters, the reactants' attributes (``x.clone``) and random draws
    (``new_id()``, ``bernoulli(p)``, ``choice(p, a, b)``, ``uniform(a, b)``,
    ``normal(m, s)``). Reactants are consumed, so ``["mtDNA"] ->
    [copy("x"), copy("x")]`` is a replication that keeps the clone (one
    molecule more), and ``[copy("x")]`` alone gives the molecule back. The
    model document writes ``{"type": ..., "copy": label, "attributes":
    {...}}``; such a dictionary is accepted as a product too.
    """
    _check_name("label", label)
    values = {
        key: _value(value, f"attribute {key} of copy({label!r})")
        for key, value in attributes.items()
    }
    return _Copy(label, values)


def new(species, /, **attributes):
    """A product of ``add_reaction``: a molecule of ``species`` with the given
    attributes.

    Every attribute of the species is given as a keyword; values are numbers
    or expressions over parameters, the reactants' attributes and random
    draws, for example ``new("mtDNA", clone="new_id()")``, which founds a
    new clone, or ``new("Cluster", size="x.size + y.size")``. A species
    without attributes takes none: ``new("A")`` is the same as ``"A"``. The
    model document writes ``{"type": species, "attributes": {...}}``; such a
    dictionary is accepted as a product too.
    """
    _check_name("species", species)
    values = {
        key: _value(value, f"attribute {key} of new({species!r})")
        for key, value in attributes.items()
    }
    return _New(species, values)


def mass_action(expression):
    """A rate for ``add_reaction``: ``expression`` per combination of molecules.

    The reaction fires at the expression's value g times the number of
    distinct combinations of reactant molecules, as a numeric rate does: one
    reactant ``x`` of a species with N molecules fires at g·N, and with a
    :func:`pattern` at g times the sum of ``factor`` over the matching
    molecules. The expression reads parameters and observables, for example
    ``mass_action("r * max(0, K - N) / N")``, and is written as the rule's
    ``rate``. A plain string rate is different: it is the reaction's whole
    propensity (the rule's ``propensity``), and a firing then takes
    reactant molecules uniformly. For one reactant without ``where`` or
    ``factor`` the two describe the same process: ``mass_action("delta")``
    and ``"delta * N"`` (with ``N`` the count) on ``["mtDNA"]`` both fire at
    δN and take a molecule uniformly. A whole propensity cannot be combined
    with ``where`` or ``factor``, and the static methods run it only over a
    species with a single instance, so a per-molecule rate on a species
    with attributes is :func:`mass_action`. A number is accepted too and
    means what it means as a rate.
    """
    return _MassAction(_expression(expression, "a mass_action expression"))


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


def _pattern_of(item):
    """A pattern dictionary as (species, label, where, factor, count)."""
    extra = sorted(set(item) - _PATTERN_KEYS)
    if extra:
        raise ValueError(
            f"unknown reactant pattern field {extra[0]!r}; a pattern has label, "
            "type, where, factor and count"
        )
    if "type" not in item:
        raise ValueError("a reactant pattern dictionary needs its 'type'")
    species = _check_name("species", item["type"])
    label = item.get("label")
    if label is not None:
        _check_name("label", label)
    where = _condition(item.get("where"), "a pattern's where")
    factor = item.get("factor")
    if factor is not None:
        factor = _expression(factor, "a pattern's factor")
    count = item.get("count")
    if count is not None:
        count = _count(count, "a pattern's count")
        if count < 1:
            raise ValueError("a pattern's count must be at least 1")
    return species, label, where, factor, count


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
    """A model built call by call: species, parameters, reactions, observables.

    ``add_*`` methods return the model, so calls can be chained. ``to_dict``
    returns a fresh document: rules-model/0.3 (0.4 if an ODE is
    independent, ``crn-craft.model/1.0`` if a reaction takes three or more
    molecules or an observable reduces with ``max_count``), which
    ``rules.ensemble`` and ``rules.prepare`` also accept as the ``Model``
    itself. ``crn_craft.Model``, the top-level name, writes
    ``crn-craft.model/1.0`` always.

    Species are plain counts by default. A species with attributes
    (``add_species(name, attributes=[...])``) is a set of instances, one per
    attribute vector, that the run creates and removes; ``add_initial``
    gives its molecules at time 0, :func:`pattern`, :func:`copy` and
    :func:`new` describe how reactions match and create them, and
    ``add_observable`` reduces them to numbers that rates, the stop
    condition and results read.

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
        # Species in declaration order: the initial count of a species
        # without attributes, 0 for a species with attributes.
        self._species = {}
        # The attribute names of each species that has attributes.
        self._attributes = {}
        # The initial entries [attribute values, count] of those species.
        self._initial = {}
        # Observables of add_observable (and of a document), by name.
        self._observables = {}
        # Species without the count observable under their name; only a
        # document of `from_dict` can leave it out.
        self._uncounted = set()
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

        A document the builder could have written stays editable: the
        ``add_*`` and ``set_*`` methods work, and after the first edit
        ``to_dict`` writes the builder's document (explicit default values
        such as ``"attributes": []`` are then left out, and the entries of
        ``initial`` and ``observables`` are grouped by species, which changes
        the bytes but not the meaning). That is a rules-model/0.3, 0.4 or
        model/1.0 document (for ``crn_craft.Model`` also rules-model/0.1 and
        0.2 without an ATP block) whose types, initial entries and
        observables have the fields of :meth:`add_species`,
        :meth:`add_initial` and :meth:`add_observable`, with one numeric
        initial count per species without attributes; every converted
        network model, the packaged examples and documents that
        :meth:`to_dict` wrote are such documents. A species whose count
        observable (an observable named after it, over all its molecules) is
        missing in the document keeps it missing. Any other document (an
        ATP block, two initial entries of one species without attributes,
        ...) is kept whole and cannot be edited: those methods raise
        ``ValueError``; edit the dictionary and call ``from_dict`` again.
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
        document = deepcopy(dict(document))
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

        The builder writes types with attribute names, one numeric initial
        entry per species without attributes and any number of entries (with
        every attribute, counts as numbers or parameter expressions) per
        species with attributes, the observables of ``add_observable`` (a
        species' count observable is the one named after it with no other
        field than ``type``), and keeps rules and continuous variables as
        given. Explicit default values (``"attributes": []`` on a type,
        ``"attributes": {}`` on an initial entry, ``"reduce": "sum"``, an
        integral count written as a float), which the converter writes, are
        accepted.
        """
        editable = _EDITABLE + (_EDITABLE_UNIFIED if self._unified else ())
        if document.get("schema") not in editable or not set(document) <= (
            _DOCUMENT_KEYS
        ):
            return False
        try:
            types = []
            attributes = {}
            for entry in document["types"]:
                names = entry.get("attributes", [])
                if set(entry) - {"name", "attributes"} or not isinstance(names, list):
                    return False
                if names:
                    attributes[entry["name"]] = _attribute_names(names)
                types.append(entry["name"])
            counts = {}
            initial = {name: [] for name in attributes}
            for entry in document["initial"]:
                count = entry["count"]
                kind = entry["type"]
                values = entry.get("attributes", {})
                if (
                    set(entry) - {"type", "attributes", "count"}
                    or kind not in types
                    or isinstance(count, bool)
                    or not isinstance(values, Mapping)
                ):
                    return False
                if isinstance(count, float) and count.is_integer():
                    count = int(count)
                if kind in attributes:
                    if set(values) != set(attributes[kind]):
                        return False
                    if not isinstance(count, str):
                        count = _count(count, "an initial count")
                    initial[kind].append(
                        [{a: values[a] for a in attributes[kind]}, count]
                    )
                    continue
                if values or kind in counts or isinstance(count, str):
                    return False
                counts[kind] = _count(count, "an initial count")
            counted = set()
            observables = {}
            for name, spec in document["observables"].items():
                if set(spec) - _OBSERVABLE_KEYS or spec["type"] not in types:
                    return False
                if (
                    name == spec["type"]
                    and set(spec) <= {"type", "reduce"}
                    and spec.get("reduce", "sum") == "sum"
                ):
                    counted.add(name)
                else:
                    observables[name] = deepcopy(dict(spec))
            self._species = {name: counts.get(name, 0) for name in types}
            self._attributes = attributes
            self._initial = initial
            self._observables = observables
            self._uncounted = set(types) - counted
            self._parameters = {
                key: _number(value, key)
                for key, value in document["parameters"].items()
            }
            self._reactions = deepcopy(list(document["rules"]))
            self._continuous = deepcopy(list(document.get("continuous", [])))
            self._stop = document.get("stop")
            self._description = document.get("description")
        except (KeyError, TypeError, ValueError, AttributeError):
            return False
        return True

    def _editable(self):
        if self._fixed:
            raise ValueError(
                f"model {self._name!r} holds a document the builder cannot "
                "write (an ATP block or other fields); edit the dictionary from "
                "to_dict() and call Model.from_dict instead"
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
            or name in self._observables
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

    def add_species(self, name, initial=_ZERO, attributes=None):
        """Declare species with initial counts (default 0), or with attributes.

        ``add_species("S", 999)`` declares one; ``add_species({"S": 999,
        "I": 1})`` declares every entry of a mapping, in its order, as one
        call each would. A mapping takes no other argument (``TypeError``).
        If an entry is invalid, nothing is declared. A species that a
        reaction already introduced keeps its position and takes the count.

        ``attributes``, a list of attribute names, declares a species whose
        molecules carry numeric attributes, for example
        ``add_species("mtDNA", attributes=["clone"])``: each attribute
        vector is an instance with its own count, and the run creates and
        removes instances. Its molecules at time 0 come from
        :meth:`add_initial`, so ``initial`` must stay 0. Declare it before
        the reactions that use it. Every species, with attributes or not,
        has an observable under its own name: its total count of molecules.
        """
        self._editable()
        entries = _entries(name, initial, "add_species", "initial")
        names = _attribute_names(attributes)
        checked = {}
        for key, count in entries:
            _check_name("species", key)
            if self._taken(key) and key not in self._implicit:
                raise ValueError(f"name {key!r} is already declared")
            if not names:
                checked[key] = _count(count, f"initial count of {key}")
                continue
            if key in self._implicit:
                raise ValueError(
                    f"species {key!r} is already used by a reaction as a species "
                    "without attributes; declare it with its attributes before "
                    "the reactions that use it"
                )
            if isinstance(count, bool | str) or count != 0:
                raise ValueError(
                    f"species {key!r} has attributes, so its molecules at time 0 "
                    f"come from add_initial({key!r}, count, {names[0]}=...), one "
                    "call per attribute vector; initial= is for species without "
                    "attributes"
                )
            checked[key] = 0
        for key, count in checked.items():
            self._species[key] = count
            self._implicit.discard(key)
            if names:
                self._attributes[key] = names
                self._initial[key] = []
        self._edited()
        return self

    def _attribute_values(self, species, values, *, every, what):
        """Checked attribute values of ``species`` in declaration order."""
        declared = self._attributes.get(species, ())
        unknown = [key for key in values if key not in declared]
        if unknown:
            known = ", ".join(declared) if declared else "none"
            raise ValueError(
                f"{what}: species {species!r} has no attribute {unknown[0]!r} "
                f"(its attributes: {known})"
            )
        if every:
            missing = [key for key in declared if key not in values]
            if missing:
                raise ValueError(
                    f"{what}: give every attribute of {species!r}; missing "
                    f"{', '.join(missing)}"
                )
        return {
            key: _value(values[key], f"attribute {key} of {what}")
            for key in declared
            if key in values
        }

    def add_initial(self, name, count, /, **attributes):
        """Add ``count`` molecules of species ``name`` at time 0.

        For a species with attributes every attribute is a keyword, for
        example ``add_initial("mtDNA", 1000, clone=0)``: the molecules of one
        instance. Values are numbers or parameter expressions, and ``count``
        is a nonnegative integer or a parameter expression such as ``"K"``.
        Call it once per attribute vector; a second call with the same
        values adds to the count. For a species without attributes it adds
        ``count``, an integer, to the initial count.
        """
        self._editable()
        if not isinstance(name, str) or name not in self._species:
            raise KeyError(f"unknown species {name!r}; declare it with add_species")
        if name not in self._attributes:
            if attributes:
                raise ValueError(
                    f"species {name!r} has no attributes; declare them with "
                    f"add_species({name!r}, attributes=[...])"
                )
            added = _count(count, f"initial count of {name}")
            total = _count(self._species[name] + added, f"initial count of {name}")
            self._species[name] = total
            self._implicit.discard(name)
            self._edited()
            return self
        values = self._attribute_values(
            name, attributes, every=True, what=f"add_initial({name!r}, ...)"
        )
        if isinstance(count, str):
            if not count.strip():
                raise ValueError("an initial count must be an integer or an expression")
        else:
            count = _count(count, f"initial count of {name}")
        for entry in self._initial[name]:
            if (
                entry[0] == values
                and isinstance(entry[1], int)
                and isinstance(count, int)
            ):
                entry[1] = _count(entry[1] + count, f"initial count of {name}")
                break
        else:
            self._initial[name].append([values, count])
        self._edited()
        return self

    def set_initial(self, counts=None, /, **kwargs):
        """Set initial counts of declared species, from a mapping or keywords.

        ``model.set_initial({"S": 990, "I": 10})`` and
        ``model.set_initial(S=990, I=10)`` are equivalent. A species with
        attributes takes its molecules from :meth:`add_initial` instead.
        """
        self._editable()
        values = dict(counts or {})
        values.update(kwargs)
        for name in values:
            if name not in self._species:
                raise KeyError(f"unknown species {name!r}; declare it first")
            if name in self._attributes:
                raise ValueError(
                    f"species {name!r} has attributes; give its molecules with "
                    f"add_initial({name!r}, count, ...)"
                )
        for name, value in values.items():
            self._species[name] = _count(value, f"initial count of {name}")
        self._edited()
        return self

    def add_observable(
        self, name, species, where=None, weight=None, reduce="sum", label=None
    ):
        """Declare an observable: a number reduced from the molecules of ``species``.

        Over the present instances of ``species`` that satisfy ``where`` (a
        condition over parameters and the instance's attributes, for example
        ``"x.clone > 0"``), with ``weight`` per instance (a number or an
        expression over the same names, default 1), ``reduce`` gives:

        - ``"sum"`` (default): the sum of weight × count, so without
          ``weight`` the number of matching molecules;
        - ``"distinct"``: the sum of the weights, each instance once, so
          without ``weight`` the number of distinct matching instances (for
          example clones);
        - ``"max"``, ``"min"``: the largest or smallest weight, 0 when no
          instance matches;
        - ``"max_count"``: the largest copy number among the matching
          instances, times the weight, 0 when none matches (for example the
          size of the largest clone); it exists only in
          ``crn-craft.model/1.0``, so ``crn_craft.rules.Model`` then writes
          that version.

        ``label`` names the instance in ``where`` and ``weight`` (default
        ``x``). Rates (strings and :func:`mass_action`), the stop condition
        and trackers read observables by name, and results have a column for
        each. The name shares one namespace with parameters, species and
        continuous variables; every species already has an observable under
        its own name, its count.
        """
        self._editable()
        _check_name("observable", name)
        if self._taken(name):
            hint = (
                " (every species has an observable of its count under its name)"
                if name in self._species
                else ""
            )
            raise ValueError(f"name {name!r} is already declared{hint}")
        _check_name("species", species)
        if species not in self._species:
            raise KeyError(f"unknown species {species!r}; declare it with add_species")
        if reduce not in REDUCTIONS:
            raise ValueError(
                f"reduce must be one of {', '.join(map(repr, REDUCTIONS))}, not "
                f"{reduce!r}"
            )
        spec = {"type": species}
        if label is not None:
            spec["label"] = _check_name("label", label)
        if where is not None:
            spec["where"] = _condition(where, f"where of observable {name}")
        if weight is not None:
            spec["weight"] = _expression(weight, f"weight of observable {name}")
        if reduce != "sum":
            spec["reduce"] = reduce
        self._observables[name] = spec
        self._edited()
        return self

    def add_reaction(self, rate, reactants, products, name=None):
        """Add a reaction ``reactants -> products``.

        ``rate`` is a number (a mass-action constant k: the reaction fires
        at k times the number of distinct combinations of reactant
        molecules, so ``["A", "A"]`` fires at k·A(A−1)/2), a string (the
        whole propensity) or :func:`mass_action` (an expression that a
        reaction multiplies by the number of combinations, like a number).
        Species are lists of names, repeated for multiplicity, for example
        ``["A", "A"]``. Three or more reactant molecules, such as
        ``["A", "A", "A"]`` or ``["A", "A", "B"]``, are written with the
        pattern field ``count`` of ``crn-craft.model/1.0`` (a numeric k
        fires at k·C(A, 2)·B for ``["A", "A", "B"]``; a string propensity
        applies while each species has its count). Species not yet declared
        are declared with count 0, as in rebop. ``name`` defaults to
        ``r<index>``.

        For species with attributes, a reactant is a species name (any
        molecule; the first is labeled ``x``, the second ``y``) or a
        :func:`pattern` with a ``where`` condition or ``factor``; a product
        is :func:`copy` (``copy("x")`` keeps the attributes of reactant
        ``x``) or :func:`new` (every attribute given), never a bare name.
        The pattern and product dictionaries of the model document are
        accepted too. A reaction with ``where`` or ``factor`` takes a number
        or :func:`mass_action` as its rate. A plain string stays the whole
        propensity: for a species with attributes the reaction fires at that
        rate and takes a molecule uniformly, whatever its attributes, and
        only the ``tree`` method (``method=None``) runs it unless the species
        has a single instance. Reactions over species with attributes take
        at most two reactant molecules.
        """
        self._editable()
        if isinstance(reactants, str | Mapping | _Pattern) or isinstance(
            products, str | Mapping | _Copy | _New
        ):
            raise TypeError(
                "reactants and products are lists of species names (or patterns "
                "and products)"
            )
        entries, explicit = self._reactants(reactants)
        products = list(products)
        if isinstance(rate, _MassAction):
            law = {"rate": rate.expression}
        elif isinstance(rate, str):
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
        # Species a reaction declares; added only when the reaction is valid.
        introduced = []

        def declare(species):
            if species in self._species or species in introduced:
                return
            if self._taken(species):
                raise ValueError(f"{species!r} is a parameter, observable or variable")
            introduced.append(species)

        for species, *_ in entries:
            declare(species)
        patterns = self._reactant_patterns(entries, explicit, law)
        labels = {p["label"]: p["type"] for p in patterns}
        written = [self._product(item, labels, declare) for item in products]
        for species in introduced:
            self._species[species] = 0
            self._implicit.add(species)
        rule = {"name": name, "reactants": patterns}
        rule.update(law)
        rule["products"] = written
        self._reactions.append(rule)
        self._edited()
        return self

    def _reactants(self, reactants):
        """Reactants as (species, label, where, factor, count) tuples, and
        whether any is more than a species name."""
        entries = []
        explicit = False
        for item in reactants:
            if isinstance(item, str):
                entries.append((_check_name("species", item), None, None, None, None))
            elif isinstance(item, _Pattern):
                entries.append(
                    (item.species, item.label, item.where, item.factor, None)
                )
                explicit |= not (
                    item.label is None and item.where is None and item.factor is None
                )
            elif isinstance(item, Mapping):
                entries.append(_pattern_of(item))
                explicit = True
            else:
                raise TypeError(
                    "a reactant is a species name, crn_craft.pattern(...) or a "
                    f"pattern dictionary, not {type(item).__name__}"
                )
        return entries, explicit

    def _reactant_patterns(self, entries, explicit, law):
        """The rule's reactant patterns."""
        if len(entries) > 2 and (
            any(e[0] in self._attributes for e in entries)
            or any(e[2] is not None or e[3] is not None for e in entries)
        ):
            raise ValueError(
                "a reaction of three or more reactant molecules takes species "
                "without attributes and patterns without where or factor"
            )
        if not explicit:
            return _patterns([e[0] for e in entries])
        many = len(entries) > 2
        patterns = []
        for i, (species, label, where, factor, count) in enumerate(entries):
            if label is None:
                label = f"x{i + 1}" if many else _LABELS[i]
            if any(p["label"] == label for p in patterns):
                raise ValueError(
                    f"two reactants have the label {label!r}; give each its own "
                    "label (the defaults are x and y by position)"
                )
            pattern = {"label": label, "type": species}
            if where is not None:
                pattern["where"] = where
            if factor is not None:
                pattern["factor"] = factor
            if count is not None:
                pattern["count"] = count
            patterns.append(pattern)
        if "propensity" in law and any("where" in p or "factor" in p for p in patterns):
            raise ValueError(
                "a string rate is the reaction's whole propensity, which does not "
                "combine with a reactant's where or factor; give the rate per "
                "combination of reactant molecules as a number or "
                "crn_craft.mass_action(...)"
            )
        return patterns

    def _product(self, item, labels, declare):
        """One product of a rule, checked against the reactant labels."""
        if isinstance(item, str):
            species = _check_name("species", item)
            if species in self._attributes:
                first = self._attributes[species][0]
                raise ValueError(
                    f"species {species!r} has attributes, so a product of it says "
                    "where they come from: crn_craft.copy(label) copies a reactant, "
                    f"crn_craft.new({species!r}, {first}=...) gives every attribute"
                )
            declare(species)
            return {"type": species}
        if isinstance(item, _Copy):
            return self._copied(item.label, item.attributes, labels)
        if isinstance(item, _New):
            return self._created(item.species, item.attributes, declare)
        if isinstance(item, Mapping):
            extra = sorted(set(item) - _PRODUCT_KEYS)
            if extra:
                raise ValueError(
                    f"unknown product field {extra[0]!r}; a product has type, "
                    "copy and attributes"
                )
            if "type" not in item:
                raise ValueError("a product dictionary needs its 'type'")
            species = _check_name("species", item["type"])
            values = item.get("attributes", {})
            if not isinstance(values, Mapping):
                raise TypeError("a product's attributes are a dictionary")
            if item.get("copy") is None:
                return self._created(species, values, declare)
            product = self._copied(item["copy"], values, labels)
            if product["type"] != species:
                raise ValueError(
                    f"a copy has the species of its reactant: {item['copy']!r} is "
                    f"{product['type']!r}, not {species!r}"
                )
            return product
        raise TypeError(
            "a product is a species name, crn_craft.copy(...), crn_craft.new(...) "
            f"or a product dictionary, not {type(item).__name__}"
        )

    def _copied(self, label, overrides, labels):
        """A copy product of the reactant ``label``."""
        if label not in labels:
            known = ", ".join(labels) if labels else "none"
            raise ValueError(
                f"copy({label!r}): no reactant has the label {label!r} (the "
                f"reactant labels: {known})"
            )
        species = labels[label]
        values = self._attribute_values(
            species, overrides, every=False, what=f"copy({label!r})"
        )
        product = {"type": species, "copy": label}
        if values:
            product["attributes"] = values
        return product

    def _created(self, species, attributes, declare):
        """A product with every attribute of ``species`` given."""
        if species not in self._species:
            if attributes:
                raise KeyError(
                    f"unknown species {species!r}; declare it with "
                    f"add_species({species!r}, attributes=[...]) first"
                )
            declare(species)
        values = self._attribute_values(
            species, attributes, every=True, what=f"new({species!r})"
        )
        product = {"type": species}
        if values:
            product["attributes"] = values
        return product

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
        observables = {
            s: {"type": s} for s in self._species if s not in self._uncounted
        }
        observables |= json.loads(json.dumps(self._observables))
        if (
            self._unified
            or any(_needs_1_0(r) for r in self._reactions)
            or any(o.get("reduce") in _REDUCTIONS_1_0 for o in observables.values())
        ):
            schema = SCHEMA_1_0
        else:
            schema = SCHEMA_INDEPENDENT if independent else SCHEMA
        types = []
        initial = []
        for species, count in self._species.items():
            names = self._attributes.get(species)
            if not names:
                types.append({"name": species})
                if count > 0:
                    initial.append({"type": species, "count": count})
                continue
            types.append({"name": species, "attributes": list(names)})
            for values, number in self._initial[species]:
                if number != 0:
                    initial.append(
                        {"type": species, "attributes": dict(values), "count": number}
                    )
        document = {
            "schema": schema,
            "name": self._name,
        }
        if self._description is not None:
            document["description"] = self._description
        document |= {
            "parameters": dict(self._parameters),
            "types": types,
            "initial": json.loads(json.dumps(initial)),
            "observables": observables,
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
        a reaction takes three or more molecules or an observable reduces
        with ``max_count``; ``crn_craft.Model`` writes
        ``crn-craft.model/1.0`` always. A model from :meth:`from_dict`
        returns its document until it is edited.
        """
        if self._document is not None:
            document = deepcopy(self._document)
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
