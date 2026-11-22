"""Execute a documented stochastic SBML Core and SED-ML/OMEX subset."""

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from pathlib import Path, PurePosixPath
import json
import math
import re
import zipfile
import zlib

import numpy as np

from ._xml import _read_xml, _XmlRefused
from .importers import ImportFailure, _identifier, _number, _poly_add, _poly_mul, fail
from .protocol import Protocol

_SIZE = 16 << 20
# Work bounds for documents from untrusted archives. Rules and laws are
# re-evaluated once per event, so the products of these counts stay small.
_MAX_PARAMETERS = 4096
_MAX_REACTIONS = 1024
_MAX_EVENTS = 256
# Reported values (paths x report times x species) of one run; the series
# are returned as Python lists, about 32 bytes per value.
_MAX_VALUES = 1 << 24
# Device names that Windows reserves in every directory.
_RESERVED = {"CON", "PRN", "AUX", "NUL"} | {
    f"{kind}{i}" for kind in ("COM", "LPT") for i in range(1, 10)
}
_MATH = "{http://www.w3.org/1998/Math/MathML}"
_OMEX = "{http://identifiers.org/combine.specifications/omex-manifest}"
_SED = re.compile(r"\{http://sed-ml.org/sed-ml/level1/version([34])\}sedML")
_CORE = re.compile(r"\{http://www.sbml.org/sbml/level3/version([12])/core\}sbml")
_XPATH = re.compile(
    r"/(?:\w+:)?sbml/(?:\w+:)?model/(?:\w+:)?listOf(Parameters|Species)/"
    r"(?:\w+:)?(parameter|species)\[@id=(['\"])([A-Za-z_][A-Za-z0-9_]*)\3\](?:/@(\w+))?"
)
_EVENTS = "http://www.sbml.org/sbml/symbols/time"


def _xml(raw):
    try:
        return _read_xml(raw, _SIZE)
    except _XmlRefused as refused:
        if refused.reason == "size":
            fail("XML document exceeds 16 MiB")
        if refused.reason in {"nul", "dtd", "encoding"}:
            fail("XML must be UTF-8 without a document type or entities")
        fail(f"invalid XML: {refused.message}")


def _read(path, limit=_SIZE):
    path = Path(path)
    try:
        # A regular file only: a directory, pipe or device is refused.
        if not path.is_file():
            fail(f"not a readable regular file: {path.name}")
        with path.open("rb") as file:
            value = file.read(limit + 1)
    except OSError as error:
        fail(f"cannot read {path.name}: {error.strerror or error}")
    if len(value) > limit:
        fail(f"file exceeds its size budget: {Path(path).name}")
    return value


def _tag(node):
    return node.tag.rsplit("}", 1)[-1]


def _child(node, name, required=False):
    if node is None:
        return None
    namespace = node.tag[: -len(_tag(node))]
    values = node.findall(namespace + name)
    if len(values) > 1 or (required and not values):
        fail(f"{_tag(node)} needs a unique {name}")
    return values[0] if values else None


def _collection(node, container, item):
    listing = _child(node, container)
    if listing is None:
        return {}
    output = {}
    for member in listing:
        if (
            _tag(member) != item
            or member.tag[: -len(item)] != listing.tag[: -len(container)]
        ):
            fail(f"unsupported element in {container}")
        name = _identifier(member.get("id"))
        if name in output:
            fail(f"duplicate ID {name}")
        output[name] = member
    return output


def _allowed(node, allowed):
    for element in node:
        tag = _tag(element)
        if tag not in allowed and tag not in {"notes", "annotation"}:
            name = f" {node.get('id')!r}" if node.get("id") else ""
            fail(f"{_tag(node)}{name} contains unsupported {tag}")


def _math(node):
    values = node.findall(_MATH + "math")
    if len(values) != 1 or len(values[0]) != 1:
        fail(f"{_tag(node)} needs one MathML expression")
    return values[0][0]


def _literal(node):
    if len(node) == 0 and node.get("type", "real") in ("integer", "real"):
        return _number(node.text)
    if len(node) == 1 and node[0].tag == _MATH + "sep":
        a, b = _number(node.text), _number(node[0].tail)
        if node.get("type") == "rational" and b:
            return a / b
        if node.get("type") == "e-notation" and b.denominator == 1 and abs(b) <= 308:
            return a * Fraction(10) ** int(b)
    fail("unsupported MathML number")


def _poly(node, symbols, zero, depth=0):
    if depth > 32:
        fail("MathML expression depth limit")
    tag = _tag(node)
    if node.tag != _MATH + tag:
        fail("non-MathML expression element")
    if tag == "cn":
        n = _literal(node)
        return {zero: n} if n else {}
    if tag == "ci" and len(node) == 0:
        name = (node.text or "").strip()
        if name not in symbols:
            fail(f"undefined scalar or species {name!r}")
        return symbols[name]
    if tag != "apply" or len(node) < 2 or len(node[0]):
        fail(f"unsupported MathML expression {tag}")
    op = _tag(node[0])
    if node[0].tag != _MATH + op:
        fail("non-MathML operator")
    args = [_poly(arg, symbols, zero, depth + 1) for arg in node[1:]]
    if op == "times":
        result = {zero: Fraction(1)}
        for arg in args:
            result = _poly_mul(result, arg)
        return result
    if op == "plus":
        result = {}
        for arg in args:
            result = _poly_add(result, arg)
        return result
    if op == "minus" and len(args) in (1, 2):
        return _poly_add(args[0] if len(args) == 2 else {}, args[-1], Fraction(-1))
    if op == "divide" and len(args) == 2 and set(args[1]) == {zero}:
        return _poly_add({}, args[0], 1 / args[1][zero])
    if op == "power" and len(args) == 2 and args[1].keys() <= {zero}:
        n = args[1].get(zero, Fraction(0))
        if n.denominator == 1 and 0 <= n <= 8:
            result = {zero: Fraction(1)}
            for _ in range(int(n)):
                result = _poly_mul(result, args[0])
            return result
    fail(f"unsupported MathML operator or arguments: {op}")


def _scalar(node, environment):
    used = _symbols_used(node) & environment.keys()
    symbols = {k: {(): environment[k]} if environment[k] else {} for k in used}
    result = _poly(node, symbols, ())
    return result.get((), Fraction(0))


def _symbols_used(node):
    return {(x.text or "").strip() for x in node.iter(_MATH + "ci")}


def _unit(unit, definitions):
    """Dimensional signature in SI bases, with one numerical scale."""
    if unit in {"item", "second", "dimensionless", "metre", "litre"}:
        bases = {"litre": "metre"}
        return (
            {bases.get(unit, unit): 3 if unit == "litre" else 1}
            if unit != "dimensionless"
            else {}
        ), (1e-3 if unit == "litre" else 1.0)
    node = definitions.get(unit)
    if node is None:
        fail(f"undefined or unsupported unit {unit!r}")
    dims, scale = {}, 1.0
    for member in _child(node, "listOfUnits", required=True):
        like = member.get("kind")
        if like not in {"item", "second", "dimensionless", "metre", "litre"}:
            fail(f"unsupported base unit {like!r} for {unit!r}")
        exponent = _number(member.get("exponent", "1"))
        decade = _number(member.get("scale", "0"))
        factor = float(_number(member.get("multiplier", "1")))
        if (
            exponent.denominator != 1
            or abs(exponent) > 8
            or decade.denominator != 1
            or abs(decade) > 80
            or factor <= 0
        ):
            fail(f"unsupported scale in {unit!r}")
        base, value = _unit(like, {})
        try:
            scale *= (factor * 10.0 ** int(decade) * value) ** int(exponent)
        except (OverflowError, ZeroDivisionError):
            fail(f"unit scale exceeds f64 in {unit!r}")
        for name, power in base.items():
            dims[name] = dims.get(name, 0) + power * int(exponent)
    if not math.isfinite(scale) or scale <= 0:
        fail(f"invalid scale in unit {unit!r}")
    return {k: v for k, v in dims.items() if v}, scale


def _count(value):
    if value.denominator != 1 or not 0 <= value <= 2**53:
        fail(
            "stochastic amount must be a nonnegative integer representable in the engine"
        )
    return int(value)


def _side(reaction, list_name, species):
    side = {}
    listing = _child(reaction, list_name)
    for reference in () if listing is None else listing:
        name = reference.get("species")
        if (
            _tag(reference) != "speciesReference"
            or name not in species
            or len(reference)
        ):
            fail(f"invalid species reference in {reaction.get('id')}")
        n = _count(_number(reference.get("stoichiometry", "1")))
        if n == 0 or n > 8 or name in side:
            fail(
                "stoichiometry must be positive, bounded, and present once per species"
            )
        side[name] = n
    return side


def _basis(reactants, names):
    zero = (0,) * len(names)
    basis = {zero: Fraction(1)}
    for species, n in reactants.items():
        unit = tuple(int(s == species) for s in names)
        for j in range(n):
            basis = _poly_mul(
                basis, {unit: Fraction(1, j + 1), zero: Fraction(-j, j + 1)}
            )
    return basis


def _rate(value):
    try:
        result = float(value)
    except OverflowError:
        fail("rate coefficient exceeds the numerical range")
    if not math.isfinite(result) or (result == 0 and value > 0):
        fail("rate coefficient exceeds the numerical range")
    return result


def _coefficients(law, reactants, products, names, split):
    forward = _basis(reactants, names)
    key = next(iter(forward))
    coefficient = law.get(key, 0) / forward[key]
    if coefficient >= 0 and law == _poly_add({}, forward, coefficient):
        return _rate(coefficient), 0.0
    if split:
        reverse = _basis(products, names)
        for key2 in sorted(forward.keys() | reverse.keys()):
            a, b = forward[key], reverse.get(key, 0)
            c, d = forward.get(key2, 0), reverse.get(key2, 0)
            det = a * d - b * c
            if not det:
                continue
            x, y = law.get(key, 0), law.get(key2, 0)
            u, v = (x * d - y * b) / det, (a * y - c * x) / det
            if (
                u >= 0
                and v <= 0
                and law == _poly_add(_poly_add({}, forward, u), reverse, v)
            ):
                return _rate(u), _rate(-v)
            break
    fail("kinetic law does not match the declared stochastic reactant factors")


def _rule_order(known, rule_nodes):
    """Order assignment rules so that each follows the rules it reads.

    Kahn's algorithm, linear in the rules; a rule that reads an undefined
    symbol, a species or itself (directly or through other rules) fails.
    """
    uses = {name: _symbols_used(node) for name, node in rule_nodes.items()}
    users = {name: [] for name in rule_nodes}
    waiting = {}
    for name, symbols in uses.items():
        unknown = symbols - known - users.keys()
        if unknown:
            fail(
                f"assignment rule {name!r} reads a non-scalar or undefined symbol "
                f"{sorted(unknown)[0]!r}"
            )
        waiting[name] = len(symbols & users.keys())
        for symbol in symbols & users.keys():
            users[symbol].append(name)
    ready = sorted(name for name, n in waiting.items() if n == 0)
    order = []
    while ready:
        name = ready.pop()
        order.append(name)
        for user in users[name]:
            waiting[user] -= 1
            if waiting[user] == 0:
                ready.append(user)
    if len(order) != len(rule_nodes):
        cycle = sorted(set(rule_nodes) - set(order))
        fail(f"cyclic SBML assignment rules: {', '.join(cycle[:4])}")
    return tuple((name, rule_nodes[name]) for name in order)


def _resolve_rules(condition, raw, ordered_rules):
    env = {**condition, **raw}
    for name, math_node in ordered_rules:
        env[name] = _scalar(math_node, env)
    return env


@dataclass
class _SbmlModel:
    model: object
    events: tuple
    volume_divisors: dict
    concentrations: frozenset
    time_scale: float
    effective_rates: dict
    # Events whose trigger `time >= 0` is already true at time 0 with
    # initialValue="true": they never fire when the run starts at time 0.
    inactive_at_zero: tuple = ()
    # SED-ML parameter changes applied before compilation.
    changes: tuple = ()


def _translate_sbml(raw, overrides, *, split=False):
    from . import Model, mass_action

    root = _xml(raw)
    match = _CORE.fullmatch(root.tag)
    if not match or root.get("level") != "3" or root.get("version") != match[1]:
        fail("the experiment adapter supports SBML L3V1 and L3V2 Core")
    if any(k.endswith("}required") and v == "true" for k, v in root.attrib.items()):
        fail("required SBML extensions are outside the experiment subset")
    node = _child(root, "model", required=True)
    _allowed(
        node,
        {
            "listOfUnitDefinitions",
            "listOfCompartments",
            "listOfSpecies",
            "listOfParameters",
            "listOfRules",
            "listOfReactions",
            "listOfEvents",
        },
    )
    units = _collection(node, "listOfUnitDefinitions", "unitDefinition")
    time_dim, seconds = _unit(node.get("timeUnits"), units)
    if time_dim != {"second": 1}:
        fail("the model must explicitly declare a time unit")
    if node.get("extentUnits") is not None:
        if _unit(node.get("extentUnits"), units) != ({"item": 1}, 1.0):
            fail("reaction extent must be in items")
    compartments = _collection(node, "listOfCompartments", "compartment")
    volumes = {}
    dimension = {}
    for name, element in compartments.items():
        _allowed(element, set())
        if element.get("constant") != "true":
            fail(f"compartment {name!r} must be fixed")
        powers = _number(element.get("spatialDimensions", "3"))
        if powers not in (0, 1, 2, 3):
            fail("compartment dimensions must be 0, 1, 2 or 3")
        size = element.get("size")
        if size is not None:
            value = _number(size)
            if value <= 0:
                fail(f"compartment {name!r} must have positive size")
            units_name = element.get("units") or node.get(
                {1: "lengthUnits", 2: "areaUnits", 3: "volumeUnits"}.get(
                    int(powers), "dimensionless"
                )
            )
            if units_name:
                dims, _ = _unit(units_name, units)
                expected = {"metre": int(powers)} if powers else {}
                if dims != expected and units_name != "dimensionless":
                    fail(f"compartment {name!r} has inconsistent size units")
            volumes[name] = value
        dimension[name] = powers
    species = _collection(node, "listOfSpecies", "species")
    if not species or len(species) > 64:
        fail("experiments support one to 64 count species")
    divisors, counts, concentrations = {}, {}, set()
    for name, element in species.items():
        _allowed(element, set())
        if (
            element.get("constant") != "false"
            or element.get("boundaryCondition") != "false"
        ):
            fail(f"species {name!r} requires ordinary dynamic counts")
        if (
            element.get("conversionFactor") is not None
            or node.get("conversionFactor") is not None
        ):
            fail("conversion factors are outside the stochastic import subset")
        declared = element.get("substanceUnits") or node.get("substanceUnits")
        if _unit(declared, units) != ({"item": 1}, 1.0):
            fail(f"species {name!r} must have item amount units")
        region = element.get("compartment")
        if region not in compartments:
            fail(f"unknown compartment {region!r}")
        use_amount = element.get("hasOnlySubstanceUnits") == "true"
        divisor = (
            volumes.get(region, Fraction(1)) if use_amount else volumes.get(region)
        )
        if not use_amount and divisor is None:
            fail(f"concentration species {name!r} needs an explicit compartment size")
        divisors[name] = Fraction(1) if use_amount else divisor
        if not use_amount:
            concentrations.add(name)
        amt, conc = element.get("initialAmount"), element.get("initialConcentration")
        if amt is not None and conc is None:
            counts[name] = _count(_number(amt))
        elif (
            conc is not None
            and amt is None
            and region in volumes
            and dimension[region] > 0
        ):
            counts[name] = _count(_number(conc) * volumes[region])
        else:
            fail(f"ambiguous initial amount for {name!r}")
    parameters = _collection(node, "listOfParameters", "parameter")
    if len(parameters) + len(compartments) > _MAX_PARAMETERS:
        fail(f"experiments support at most {_MAX_PARAMETERS} parameters")
    if len(set(species) | set(parameters) | set(compartments)) != len(species) + len(
        parameters
    ) + len(compartments):
        fail("duplicate SBML symbol")
    rules = _child(node, "listOfRules")
    assignments = {}
    for rule in () if rules is None else rules:
        name = rule.get("variable")
        if (
            _tag(rule) != "assignmentRule"
            or name not in parameters
            or name in assignments
            or parameters[name].get("constant") != "false"
        ):
            fail("only acyclic scalar parameter assignment rules are supported")
        assignments[name] = _math(rule)
    if overrides.keys() - parameters.keys() or overrides.keys() & assignments.keys():
        fail("model change must target a declared independent SBML parameter")
    ordered_rules = _rule_order(
        volumes.keys() | (parameters.keys() - assignments.keys()), assignments
    )
    raw_parameters = {
        name: _number(p.get("value"))
        for name, p in parameters.items()
        if name not in assignments
    }
    raw_parameters.update(overrides)
    reactions = _collection(node, "listOfReactions", "reaction")
    if not reactions or len(reactions) > _MAX_REACTIONS:
        fail(f"SBML experiment model needs 1 to {_MAX_REACTIONS} reactions")
    names = tuple(species)
    zero = (0,) * len(names)
    rx = []
    for name, reaction in reactions.items():
        _allowed(
            reaction,
            {"listOfReactants", "listOfProducts", "kineticLaw", "listOfModifiers"},
        )
        if reaction.get("reversible") != "false" or reaction.get("fast") == "true":
            fail(f"reaction {name!r} requires an irreversible static jump")
        reactants = _side(reaction, "listOfReactants", species)
        products = _side(reaction, "listOfProducts", species)
        law = _child(reaction, "kineticLaw", required=True)
        _allowed(law, {"math", "listOfLocalParameters"})
        local = {
            k: _number(v.get("value"))
            for k, v in _collection(
                law, "listOfLocalParameters", "localParameter"
            ).items()
        }
        rx.append((name, reactants, products, _math(law), local))

    unit_counts = {
        s: {tuple(int(j == k) for j in range(len(names))): 1 / divisors[s]}
        for k, s in enumerate(names)
    }

    def rates(current):
        resolved = _resolve_rules(volumes, current, ordered_rules)
        shared = unit_counts | {k: {zero: v} if v else {} for k, v in resolved.items()}
        result = {}
        for i, (name, reactants, products, math_node, local) in enumerate(rx):
            context = (
                shared | {k: {zero: v} if v else {} for k, v in local.items()}
                if local
                else shared
            )
            try:
                forward, reverse = _coefficients(
                    _poly(math_node, context, zero),
                    reactants,
                    products,
                    names,
                    split,
                )
            except ImportFailure as error:
                fail(f"reaction {name!r}: {error}")
            result[f"sbml_rate_{i}"] = forward
            if split:
                result[f"sbml_reverse_{i}"] = reverse
        return resolved, result

    environment, initial_rates = rates(raw_parameters)
    title = node.get("id", "ImportedSBML")
    model = Model(_identifier(title))
    model.add_species(counts)
    model.add_parameter(initial_rates)
    for i, (name, reactants, products, _, _) in enumerate(rx):
        source = [s for s, n in reactants.items() for _ in range(n)]
        target = [s for s, n in products.items() for _ in range(n)]
        model.add_reaction(mass_action(f"sbml_rate_{i}"), source, target, name=name)
        if split:
            model.add_reaction(
                mass_action(f"sbml_reverse_{i}"),
                target,
                source,
                name=f"sbml_reverse_{i}",
            )

    schedule = _events(
        node,
        environment,
        raw_parameters,
        assignments,
        parameters,
        counts,
        divisors,
        rates,
    )
    schedule, inactive = schedule
    return _SbmlModel(
        model,
        schedule,
        divisors,
        frozenset(concentrations),
        seconds,
        initial_rates,
        inactive,
        tuple(sorted((k, float(v)) for k, v in overrides.items())),
    )


def _events(node, initial, parameters, rules, declarations, counts, divisors, rates):
    events = _collection(node, "listOfEvents", "event")
    if len(events) > _MAX_EVENTS:
        fail(f"experiments support at most {_MAX_EVENTS} events")
    timing, inactive = [], []
    for name, event in events.items():
        _allowed(event, {"trigger", "priority", "listOfEventAssignments"})
        trigger = _child(event, "trigger", required=True)
        math_node = _math(trigger)
        if (
            math_node.tag != _MATH + "apply"
            or len(math_node) != 3
            or math_node[0].tag != _MATH + "geq"
            or math_node[1].tag != _MATH + "csymbol"
            or math_node[1].get("definitionURL") != _EVENTS
            or math_node[2].tag != _MATH + "cn"
        ):
            fail(f"event {name!r} needs a literal time >= threshold trigger")
        when = _literal(math_node[2])
        if when < 0:
            fail(f"event {name!r} has a negative threshold")
        if when == 0 and trigger.get("initialValue", "true") == "true":
            inactive.append(name)
            continue
        priority = _child(event, "priority")
        order = None if priority is None else _scalar(_math(priority), {})
        assignments = _child(event, "listOfEventAssignments", required=True)
        keys = [item.get("variable") for item in assignments]
        if len(keys) != len(set(keys)) or any(
            _tag(item) != "eventAssignment" for item in assignments
        ):
            fail(f"ambiguous assignments in event {name!r}")
        # Simultaneous events read the state before the instant unless they
        # ask for values at execution (useValuesFromTriggerTime="false").
        at_trigger = event.get("useValuesFromTriggerTime", "true") != "false"
        timing.append(
            (when, order, name, tuple(zip(keys, assignments, strict=True)), at_trigger)
        )
    timing.sort(key=lambda item: (item[0], -(item[1] or 0), item[2]))
    output = []
    before = initial
    for i, (when, order, name, assignments, at_trigger) in enumerate(timing):
        if not i or when != timing[i - 1][0]:
            before = initial
        if (
            i
            and when == timing[i - 1][0]
            and (order is None or timing[i - 1][1] is None or order == timing[i - 1][1])
        ):
            fail("simultaneous events need distinct literal priorities")
        updates, copies = {}, {}
        for key, element in assignments:
            value = _scalar(_math(element), before if at_trigger else initial)
            if key in counts:
                copies[key] = _count(value * divisors[key])
            elif (
                key in parameters
                and key not in rules
                and declarations[key].get("constant") == "false"
            ):
                updates[key] = value
            else:
                fail(f"event {name!r} changes an unsupported or constant symbol")
        parameters = {**parameters, **updates}
        initial, effective = rates(parameters)
        output.append((float(when), effective if updates else {}, copies))
    return tuple(output), tuple(inactive)


def _member(name):
    if not isinstance(name, str) or "\\" in name or ":" in name or "\x00" in name:
        fail("invalid archive member or model URI")
    path = PurePosixPath(name)
    if not path.parts or path.is_absolute() or ".." in path.parts:
        fail("archive member must be relative to the experiment")
    return path.as_posix()


def _omex(path):
    raw = _read(path, 64 << 20)
    with zipfile.ZipFile(BytesIO(raw)) as archive:
        info = archive.infolist()
        if len(info) > 256 or sum(member.file_size for member in info) > 64 << 20:
            fail("archive exceeds member or uncompressed size limits")
        files = {}
        for member in info:
            name = _member(member.filename)
            if name in files:
                fail("duplicate archive member")
            if not member.is_dir():
                if ((member.external_attr >> 16) & 0o170000) == 0o120000:
                    fail("archive links are unsupported")
                files[name] = member
        if "manifest.xml" not in files:
            fail("OMEX archive needs a manifest")
        manifest = _xml(archive.read("manifest.xml"))
        if manifest.tag != _OMEX + "omexManifest":
            fail("unknown OMEX manifest")
        masters = []
        experiment_names = []
        for row in manifest:
            if row.tag != _OMEX + "content":
                fail("unknown OMEX manifest element")
            if "sed-ml" in row.get("format", "") and row.get("location") != ".":
                name = _member(row.get("location"))
                experiment_names.append(name)
                if row.get("master") == "true":
                    masters.append(name)
        chosen = masters or experiment_names
        if len(chosen) != 1 or chosen[0] not in files:
            fail("select an OMEX archive with one resolved master SED-ML experiment")
        data = {
            name: archive.read(member)
            for name, member in files.items()
            if member.file_size <= _SIZE
        }
    return data, chosen[0], sha256(raw).hexdigest()


def _target(text, kind, attribute=None):
    match = _XPATH.fullmatch(text or "")
    if not match or match[2] != kind or match[5] != attribute:
        fail(f"unsupported {kind} target {text!r}")
    return match[4]


def _report_series(metadata, variable, divisor, start):
    if variable is None:
        return [t + start for t in metadata["report_times"]]
    col = metadata["observables"].index(variable)
    return [
        [row[col] / divisor for row in path["values"]] for path in metadata["paths"]
    ]


@dataclass(frozen=True)
class Experiment:
    """A resolved experiment, ready for preflight or independently seeded replay."""

    name: str
    source_sha256: str
    member_sha256: dict
    convention: str
    tasks: dict
    outputs: dict

    def preflight(self):
        """List effective rates, units, report grids and output mappings."""
        return {
            "source": self.name,
            "source_sha256": self.source_sha256,
            "member_sha256": self.member_sha256,
            "convention": self.convention,
            "tasks": {
                k: {
                    "method": "KISAO:0000029",
                    "times": data[1].tolist(),
                    "time_scale_seconds": data[0].time_scale,
                    "parameter_changes": dict(data[0].changes),
                    "initial_rates": data[0].effective_rates,
                    "events": [
                        {"time": t, "parameters": p, "counts": c}
                        for t, p, c in data[0].events
                    ],
                    "events_inactive_at_time_zero": list(data[0].inactive_at_zero),
                    "repetitions": data[3],
                    "species_units": {
                        s: "item"
                        if s not in data[0].concentrations
                        else "item per model size unit"
                        for s in data[0].volume_divisors
                    },
                }
                for k, data in self.tasks.items()
            },
            "outputs": self.outputs,
        }

    def run(self, *, seed=0, replicates=1, threads=1, archive_directory=None):
        """Execute tasks in declaration order and return the requested series."""
        from . import ensemble

        values = sum(
            replicates * data[3] * len(data[1]) * len(data[0].volume_divisors)
            for data in self.tasks.values()
        )
        if values > _MAX_VALUES:
            fail(
                f"the run would return {values} values, above {_MAX_VALUES}; "
                "use fewer replicates"
            )
        archive = None if archive_directory is None else Path(archive_directory)
        if archive is not None:
            # Task IDs name the task directories: refuse names that collide on
            # a case-insensitive file system or that Windows reserves.
            folded = [name.upper() for name in self.tasks]
            if len(set(folded)) != len(folded) or _RESERVED & set(folded):
                fail("task IDs cannot name distinct archive directories")
            archive.mkdir()
        executed = {}
        for i, (name, (spec, times, start, repetitions)) in enumerate(
            self.tasks.items()
        ):
            protocol = Protocol()
            for t, pars, copies in spec.events:
                if start <= t <= times[-1]:
                    protocol.at(t - start, parameters=pars, set_counts=copies)
                elif t < start:
                    fail("experiment starts after a scheduled event")
            n = replicates * repetitions
            run = ensemble(
                spec.model,
                end=float(times[-1] - start),
                report_times=times - start,
                seed=seed + i,
                replicates=n,
                threads=threads,
                method="static_direct",
                protocol=protocol if spec.events else None,
                keep_paths=n,
                quantiles=(),
                strict=True,
            )
            if archive is None:
                result = run.simulate().metadata
            else:
                report = run.archive(archive / name)
                if not (report["valid"] and report["simulation_complete"]):
                    raise RuntimeError(
                        f"experiment task {name} did not finish: {report}"
                    )
                result = json.loads((archive / name / "result.json").read_bytes())
            executed[name] = result
        datasets = {}
        for report, series in self.outputs.items():
            datasets[report] = {}
            for label, (task, variable) in series.items():
                spec, _, start, _ = self.tasks[task]
                divisor = float(spec.volume_divisors.get(variable, 1))
                datasets[report][label] = _report_series(
                    executed[task], variable, divisor, start
                )
        result = {
            **self.preflight(),
            "schema": "crn-craft.experiment-result/1.0",
            "reports": datasets,
            "seed": seed,
            "replicates": replicates,
            "task_sampler": {k: m["sampler"] for k, m in executed.items()},
        }
        if archive is not None:
            (archive / "experiment.json").write_text(
                json.dumps(result, indent=2) + "\n"
            )
        return result


def read_experiment(path, *, convention):
    """Read an OMEX or SED-ML path and resolve its local stochastic experiment.

    The caller explicitly confirms stochastic count propensity semantics.
    Exact polynomial conversion verifies the falling-factorial reactants.
    """
    if convention not in {"stochastic_count_propensity", "split_signed_mass_action"}:
        fail("an explicit stochastic count convention is required")
    split = convention == "split_signed_mass_action"
    path = Path(path)
    if path.suffix.lower() in {".omex", ".zip"}:
        try:
            files, master, source_hash = _omex(path)
        except (
            ValueError,
            zipfile.BadZipFile,
            RuntimeError,
            NotImplementedError,
            EOFError,
            OSError,
            zlib.error,
        ) as error:
            fail(f"invalid OMEX archive: {error}")

        def read(member):
            if member not in files:
                fail(f"absent or oversized archive member {member!r}")
            return files[member]
    else:
        master, directory = path.name, path.resolve().parent
        source_hash = sha256(_read(path)).hexdigest()

        def read(member):
            candidate = (directory / member).resolve()
            if not candidate.is_relative_to(directory):
                fail("source escapes the SED-ML directory")
            return _read(candidate)

    document = _xml(read(master))
    match = _SED.fullmatch(document.tag)
    if not match or document.get("level") != "1" or document.get("version") != match[1]:
        fail("supported SED-ML versions are L1V3 and L1V4")
    _allowed(
        document,
        {
            "listOfModels",
            "listOfSimulations",
            "listOfTasks",
            "listOfDataGenerators",
            "listOfOutputs",
        },
    )
    model_nodes = _collection(document, "listOfModels", "model")
    cached_models = {}
    hashes = {master: sha256(read(master)).hexdigest()}

    def load_model(key, seen):
        if key in seen or key not in model_nodes:
            fail(f"cyclic or missing model reference {key!r}")
        if key in cached_models:
            return cached_models[key]
        node = model_nodes[key]
        if not node.get("language", "").startswith("urn:sedml:language:sbml"):
            fail(f"model {key!r} requires SBML")
        source = node.get("source", "")
        if source.startswith("#"):
            member, _, prior = load_model(source[1:], (*seen, key))
        else:
            member = _member(str(PurePosixPath(master).parent / _member(source)))
            prior = {}
        changes = _collection_changes(node)
        model_bytes = read(member)
        hashes[member] = sha256(model_bytes).hexdigest()
        overrides = {**prior, **changes}
        compiled = _translate_sbml(model_bytes, overrides, split=split)
        cached_models[key] = (member, compiled, overrides)
        return cached_models[key]

    simulations = _collection(document, "listOfSimulations", "uniformTimeCourse")
    model_tasks, task_aliases = {}, {}
    listing = _child(document, "listOfTasks", required=True)
    repeated = []
    for task in listing:
        key = _identifier(task.get("id"))
        if key in task_aliases:
            fail("duplicate SED-ML task ID")
        task_aliases[key] = key
        if _tag(task) == "repeatedTask":
            repeated.append((key, task))
            continue
        if _tag(task) != "task":
            fail("unsupported SED-ML task")
        simulation = simulations.get(task.get("simulationReference"))
        if simulation is None:
            fail("task needs a uniform time course")
        algorithm = _child(simulation, "algorithm", required=True)
        if algorithm.get("kisaoID") != "KISAO:0000029" or len(algorithm):
            fail(
                f"simulation {simulation.get('id')!r}: the stochastic adapter supports "
                "Gillespie direct (KISAO:0000029) without algorithm parameters"
            )
        initial, first, end = [
            float(_number(simulation.get(n)))
            for n in ("initialTime", "outputStartTime", "outputEndTime")
        ]
        if (
            simulation.get("numberOfPoints") is not None
            and simulation.get("numberOfSteps") is not None
        ):
            fail("uniform time course declares two grid sizes")
        steps = _count(
            _number(simulation.get("numberOfPoints") or simulation.get("numberOfSteps"))
        )
        if not initial <= first < end or not 1 <= steps <= 100_000:
            fail("invalid uniform time course")
        model = load_model(task.get("modelReference"), ())[1]
        if initial < 0 and model.inactive_at_zero:
            fail(
                f"event {model.inactive_at_zero[0]!r} at time 0 needs a run that "
                "starts at time 0"
            )
        model_tasks[key] = (model, np.linspace(first, end, steps + 1), initial, 1)
    for key, task in repeated:
        _allowed(task, {"listOfRanges", "listOfSubTasks"})
        subs = _child(task, "listOfSubTasks", required=True)
        ranges = _child(task, "listOfRanges", required=True)
        if len(subs) != 1 or len(ranges) != 1 or task.get("resetModel") != "true":
            fail("repeated task requires one subtask and a reset model")
        src, iteration = subs[0].get("task"), ranges[0]
        if (
            src not in model_tasks
            or _tag(iteration) != "uniformRange"
            or task.get("range") != iteration.get("id")
        ):
            fail("unsupported repeated task range or subtask")
        endpoint = _count(
            _number(iteration.get("numberOfPoints") or iteration.get("numberOfSteps"))
        )
        if not 1 <= endpoint <= 1000:
            fail("repeated task exceeds its run cap")
        model_tasks[key] = (*model_tasks[src][:3], endpoint + 1)

    outputs = _outputs(document, model_tasks)
    needed = {task for report in outputs.values() for task, _ in report.values()}
    return Experiment(
        master,
        source_hash,
        hashes,
        convention,
        {k: x for k, x in model_tasks.items() if k in needed},
        outputs,
    )


def _collection_changes(node):
    edits = _child(node, "listOfChanges")
    result = {}
    for edit in () if edits is None else edits:
        if _tag(edit) != "changeAttribute":
            fail("only explicit parameter value changes are supported")
        name = _target(edit.get("target"), "parameter", "value")
        if name in result:
            fail("parameter changed twice by one model")
        result[name] = _number(edit.get("newValue"))
    return result


def _outputs(document, tasks):
    fixtures = _collection(document, "listOfDataGenerators", "dataGenerator")
    data = {}
    for name, expression in fixtures.items():
        _allowed(expression, {"listOfVariables", "math"})
        variables = _collection(expression, "listOfVariables", "variable")
        expression = _math(expression)
        variable_id = (expression.text or "").strip()
        if (
            expression.tag != _MATH + "ci"
            or len(variables) != 1
            or variable_id not in variables
        ):
            fail(f"data generator {name!r} must select a time or species series")
        v = variables[variable_id]
        task = v.get("taskReference")
        if task not in tasks:
            fail(f"data generator {name!r} needs a task")
        if v.get("symbol") == "urn:sedml:symbol:time":
            quantity = None
        else:
            quantity = _target(v.get("target"), "species")
            if quantity not in tasks[task][0].volume_divisors:
                fail(f"data generator {name!r} has an unknown species")
        data[name] = (task, quantity)
    outputs = _child(document, "listOfOutputs", required=True)
    result = {}
    for item in outputs:
        name = _identifier(item.get("id"))
        if name in result:
            fail("duplicate experiment output")
        if _tag(item) == "report":
            sets = _collection(item, "listOfDataSets", "dataSet")
            series = [(k, v.get("dataReference")) for k, v in sets.items()]
        elif _tag(item) == "plot2D":
            curves = _collection(item, "listOfCurves", "curve")
            series = [
                (f"{key}.{axis}", c.get(f"{axis}DataReference"))
                for key, c in curves.items()
                for axis in ("x", "y")
            ]
        else:
            fail(f"unsupported experiment output {_tag(item)!r}")
        if any(source not in data for _, source in series) or not series:
            fail(f"output {name!r} has a missing data generator")
        result[name] = {key: data[source] for key, source in series}
    if not result:
        fail("SED-ML experiment has no report or curve data")
    return result


def main():
    """Run or preflight an SBML/SED-ML experiment from the command line."""
    import argparse

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("experiment", type=Path)
    parser.add_argument(
        "--convention",
        required=True,
        choices=("stochastic_count_propensity", "split_signed_mass_action"),
    )
    parser.add_argument(
        "--check", action="store_true", help="resolve the experiment without simulating"
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    try:
        study = read_experiment(args.experiment, convention=args.convention)
        result = (
            study.preflight()
            if args.check
            else study.run(
                seed=args.seed,
                replicates=args.replicates,
                threads=args.threads,
                archive_directory=args.archive,
            )
        )
    except (ValueError, RuntimeError, OSError) as error:
        # Refused input carries CRN-E090; a file system failure CRN-E030.
        code = getattr(error, "code", "CRN-E030" if isinstance(error, OSError) else "")
        parser.exit(2, f"{code + ': ' if code else ''}{error}\n")
    print(json.dumps(result, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
