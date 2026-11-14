"""Explicit finite-network import subsets. No silent SBML/rule expansion.

The importers build a network model (`crn-craft.model/0.1`) exactly as in
0.6.0, validate it, and return it converted to `crn-craft.model/1.0` by
:func:`crn_craft.convert`, inside a `crn-craft.import/1.0` report that keeps
the source's provenance. ``report["model"]`` is the converted model as a
dictionary, ready for the unified API (a later ``read_sbml``/``read_bngl``
wraps it in a ``Model``).
"""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET


class ImportFailure(ValueError):
    """Raised by an importer for input it does not support.

    A ``ValueError`` with the attribute ``code`` ``"CRN-E090"``.
    """

    code = "CRN-E090"


def fail(message):
    raise ImportFailure(message)


# Plain decimal literals with at most a four-digit exponent: Fraction would
# otherwise spend unbounded time on forms such as "1e-16000000", and accepts
# "1/10" and "1_0", which are not SBML or BioNetGen numbers.
_LITERAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]{1,4})?")


def _number(text):
    if not isinstance(text, str) or not _LITERAL.fullmatch(text.strip()):
        fail("finite decimal numeric literal required")
    try:
        value = Fraction(text.strip())
    except (ValueError, ZeroDivisionError, OverflowError):
        fail("finite numeric literal required")
    try:
        finite = math.isfinite(float(value))
    except OverflowError:
        finite = False
    if not finite:
        fail("numeric literal exceeds f64")
    return value


_INDEX = re.compile(r"[0-9]{1,9}")


def _identifier(name):
    if not name or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,63}", name):
        fail("unsupported or missing identifier")
    return name


def _model(name):
    return {
        "schema": "crn-craft.model/0.1",
        "name": _identifier(name),
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {},
        "species": [],
        "reactions": [],
    }


def _report(raw, source, model, **extra):
    # The converter applies the native strict domain validation of network
    # models; import is not biological validation.
    from ._common import CrnError, convert

    try:
        # The converter's canonical bytes, parsed in their key order.
        converted = json.loads(convert(model))
    except CrnError as error:
        # The network model's validation errors, or a model the converter
        # cannot hold, for example a network without reactions (CRN-R016).
        fail(json.dumps(error.diagnostics))
    return {
        "schema": "crn-craft.import/1.0",
        "source_format": source,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "model": converted,
        **extra,
    }


def _bounded(value):
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > 4096:
        fail("exact arithmetic work limit")
    return value


def _poly_add(a, b, scale=Fraction(1)):
    out = dict(a)
    for key, value in b.items():
        out[key] = _bounded(out.get(key, Fraction(0)) + _bounded(scale * value))
    return {k: v for k, v in out.items() if v}


def _poly_mul(a, b):
    out = {}
    if len(a) * len(b) > 4096:
        fail("polynomial work limit")
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(x + y for x, y in zip(ka, kb, strict=True))
            out[key] = _bounded(out.get(key, Fraction(0)) + _bounded(va * vb))
    if any(sum(k) > 8 for k in out):
        fail("polynomial degree limit")
    return {k: v for k, v in out.items() if v}


def import_sbml(raw, *, convention):
    """L3V2 Core, one fixed unit compartment, integer item amounts, polynomial jumps.

    convention must be explicitly 'stochastic_count_propensity'. Parameters are
    folded into reaction coefficients; the report names that import-time choice.
    Returns a ``crn-craft.import/1.0`` report as a dictionary:
    ``report["model"]`` is the converted ``crn-craft.model/1.0`` model, beside
    the source format and its SHA-256. Unsupported input raises
    :class:`crn_craft.importers.ImportFailure`.
    """
    if convention != "stochastic_count_propensity":
        fail("explicit stochastic_count_propensity convention required")
    if (
        len(raw) > 16 * 1024 * 1024
        or b"<!DOCTYPE" in raw.upper()
        or b"<!ENTITY" in raw.upper()
    ):
        fail("oversized XML or DTD/entities unsupported")
    try:
        raw.decode("utf-8")
        root = ET.fromstring(raw)
    except (ET.ParseError, UnicodeDecodeError, LookupError, ValueError) as e:
        fail(f"unreadable XML: {e}")
    sbml = "http://www.sbml.org/sbml/level3/version2/core"
    mathns = "http://www.w3.org/1998/Math/MathML"
    if (
        root.tag != f"{{{sbml}}}sbml"
        or root.get("level") != "3"
        or root.get("version") != "2"
    ):
        fail("only SBML Level 3 Version 2 Core is supported")

    def tag(x):
        if not x.tag.startswith("{" + sbml + "}"):
            fail("SBML extension namespaces unsupported")
        return x.tag.split("}", 1)[1]

    def attrs(x, allowed):
        if set(x.attrib) - set(allowed) - {"metaid", "sboTerm", "name"}:
            fail("unsupported SBML attributes on " + tag(x))

    def children(x, allowed):
        for child in x:
            if tag(child) not in allowed:
                fail("unsupported SBML element " + tag(child))

    attrs(root, ["level", "version"])
    children(root, {"model"})
    if len(root) != 1:
        fail("exactly one model required")
    node = root[0]
    attrs(
        node,
        [
            "id",
            "substanceUnits",
            "timeUnits",
            "volumeUnits",
            "areaUnits",
            "lengthUnits",
            "extentUnits",
        ],
    )
    children(
        node,
        {"listOfCompartments", "listOfSpecies", "listOfParameters", "listOfReactions"},
    )
    # L3V2 has no unit defaults: missing units are undefined, not items/seconds.
    if node.get("extentUnits", "item") != "item":
        fail("only item extent units supported")
    if node.get("substanceUnits") not in [None, "item"]:
        fail("only item substance units supported")
    if node.get("timeUnits") != "second":
        fail("explicit second time units required")
    model = _model(node.get("id", "ImportedSBML"))
    groups = {}
    for child in node:
        key = tag(child)
        if key in groups:
            fail("duplicate SBML list")
        attrs(child, [])
        groups[key] = child
    compartments = groups.get("listOfCompartments", [])
    if len(compartments) != 1:
        fail("one fixed compartment required")
    compartment = compartments[0]
    attrs(compartment, ["id", "size", "constant", "spatialDimensions", "units"])
    children(compartment, set())
    if (
        tag(compartment) != "compartment"
        or compartment.get("constant") != "true"
        or _number(compartment.get("size", "undefined")) != 1
    ):
        fail("unit-size constant compartment required")
    if compartment.get("units", "litre") not in ["litre", "dimensionless"]:
        fail("unsupported compartment units")
    compid = _identifier(compartment.get("id"))
    symbols = {compid: Fraction(1)}
    for species in groups.get("listOfSpecies", []):
        if tag(species) != "species":
            fail("non-species in species list")
        attrs(
            species,
            [
                "id",
                "compartment",
                "initialAmount",
                "hasOnlySubstanceUnits",
                "boundaryCondition",
                "constant",
                "substanceUnits",
            ],
        )
        children(species, set())
        name = _identifier(species.get("id"))
        count = _number(species.get("initialAmount", "undefined"))
        if name in symbols or any(s["name"] == name for s in model["species"]):
            fail("duplicate SBML identifier")
        if count.denominator != 1 or count < 0 or count > 2**63 - 1:
            fail("integer non-negative initialAmount required")
        if (
            species.get("compartment") != compid
            or species.get("hasOnlySubstanceUnits") != "true"
            or species.get("boundaryCondition") != "false"
            or species.get("constant") != "false"
            or species.get("substanceUnits", node.get("substanceUnits")) != "item"
        ):
            fail("unsupported species semantics or undefined substance units")
        model["species"].append({"name": name, "initial_count": int(count)})
    names = [s["name"] for s in model["species"]]
    if len(names) > 64:
        fail("SBML polynomial subset supports at most 64 species")
    zero = (0,) * len(names)
    for parameter in groups.get("listOfParameters", []):
        if tag(parameter) != "parameter":
            fail("non-parameter in parameter list")
        attrs(parameter, ["id", "value", "constant", "units"])
        children(parameter, set())
        name = _identifier(parameter.get("id"))
        if name in symbols or name in names:
            fail("duplicate SBML identifier")
        if (
            parameter.get("constant") != "true"
            or parameter.get("units", "dimensionless") != "dimensionless"
        ):
            fail(
                "only constant dimensionless parameter values supported; provide jump-rate units by explicit convention"
            )
        symbols[name] = _number(parameter.get("value", "undefined"))

    def polynomial(element, depth=0):
        if depth > 32:
            fail("MathML depth limit")
        if not element.tag.startswith("{" + mathns + "}"):
            fail("MathML namespace required")
        kind = element.tag.split("}", 1)[1]
        if element.attrib and not (
            kind == "cn" and element.attrib == {"type": "integer"}
        ):
            fail("unsupported MathML attributes")
        if kind == "ci":
            name = (element.text or "").strip()
            if len(element):
                fail("invalid identifier MathML")
            if name in symbols:
                return {zero: symbols[name]}
            if name in names:
                return {
                    tuple(
                        int(i == names.index(name)) for i in range(len(names))
                    ): Fraction(1)
                }
            fail("unknown MathML symbol " + name)
        if kind == "cn":
            if len(element):
                fail("structured MathML number unsupported")
            return {zero: _number((element.text or "").strip())}
        if kind != "apply" or len(element) < 2:
            fail("unsupported MathML node " + kind)
        op = element[0]
        if not op.tag.startswith("{" + mathns + "}") or len(op) or op.attrib:
            fail("unsupported MathML operator")
        op = op.tag.split("}", 1)[1]
        args = [polynomial(e, depth + 1) for e in list(element)[1:]]
        if op == "times":
            out = {zero: Fraction(1)}
            for arg in args:
                out = _poly_mul(out, arg)
            return out
        if op == "plus":
            out = {}
            for arg in args:
                out = _poly_add(out, arg)
            return out
        if op == "minus" and len(args) in [1, 2]:
            return (
                _poly_add({}, args[0], Fraction(-1))
                if len(args) == 1
                else _poly_add(args[0], args[1], Fraction(-1))
            )
        if (
            op == "divide"
            and len(args) == 2
            and set(args[1]) == {zero}
            and args[1][zero] != 0
        ):
            return {k: v / args[1][zero] for k, v in args[0].items()}
        if op == "power" and len(args) == 2 and set(args[1]) == {zero}:
            n = args[1][zero]
            if n.denominator != 1 or not 0 <= n <= 8:
                fail("integer MathML powers 0..8 only")
            out = {zero: Fraction(1)}
            for _ in range(int(n)):
                out = _poly_mul(out, args[0])
            return out
        fail("unsupported MathML operation " + op)

    for rx in groups.get("listOfReactions", []):
        if tag(rx) != "reaction":
            fail("non-reaction in reaction list")
        # L3V2 has no "fast" attribute; it is rejected as unsupported.
        attrs(rx, ["id", "reversible"])
        children(rx, {"listOfReactants", "listOfProducts", "kineticLaw"})
        if rx.get("reversible") != "false":
            fail("reversible reactions must be expanded explicitly")
        reaction_id = _identifier(rx.get("id"))
        if (
            reaction_id in symbols
            or reaction_id in names
            or any(r["name"] == reaction_id for r in model["reactions"])
        ):
            fail("duplicate SBML identifier")
        reaction = {"name": reaction_id, "reactants": {}, "products": {}}
        law = None
        seen = set()
        for child in rx:
            key = tag(child)
            if key in seen:
                fail("duplicate reaction list/law")
            seen.add(key)
            if key == "kineticLaw":
                attrs(child, [])
                if (
                    len(child) != 1
                    or child[0].tag != f"{{{mathns}}}math"
                    or len(child[0]) != 1
                    or child[0].attrib
                ):
                    fail(
                        "one pure MathML kinetic law required; local parameters unsupported"
                    )
                law = polynomial(child[0][0])
            else:
                attrs(child, [])
                side = reaction["reactants" if key == "listOfReactants" else "products"]
                for ref in child:
                    if tag(ref) != "speciesReference":
                        fail("unsupported reaction reference")
                    attrs(ref, ["species", "stoichiometry", "constant"])
                    children(ref, set())
                    name = ref.get("species")
                    # L3V2 requires both attributes; neither has a default.
                    if ref.get("stoichiometry") is None or ref.get("constant") is None:
                        fail("speciesReference stoichiometry and constant required")
                    stoich = _number(ref.get("stoichiometry"))
                    if (
                        name not in names
                        or name in side
                        or stoich.denominator != 1
                        or not 1 <= stoich <= 8
                        or ref.get("constant") != "true"
                    ):
                        fail("invalid species reference")
                    side[name] = int(stoich)
        if law is None:
            fail("missing kinetic law")
        basis = {zero: Fraction(1)}
        for s, n in reaction["reactants"].items():
            monomial = tuple(int(i == names.index(s)) for i in range(len(names)))
            for j in range(n):
                basis = _poly_mul(
                    basis, {monomial: Fraction(1, j + 1), zero: Fraction(-j, j + 1)}
                )
        key = next(iter(basis))
        coefficient = law.get(key, Fraction(0)) / basis[key]
        expected = {k: v * coefficient for k, v in basis.items() if v * coefficient}
        if coefficient < 0 or law != expected:
            fail(
                "kinetic law is not a non-negative multiple of the stochastic binomial reactant factors"
            )
        try:
            numeric_coefficient = float(coefficient)
        except OverflowError:
            fail("SBML coefficient exceeds f64")
        if not math.isfinite(numeric_coefficient) or (
            coefficient > 0 and numeric_coefficient == 0
        ):
            fail("SBML coefficient overflows or underflows f64")
        reaction["rate"] = {"kind": "mass_action", "coefficient": numeric_coefficient}
        model["reactions"].append(reaction)
    return _report(
        raw,
        "sbml_l3v2_core_subset",
        model,
        convention=convention,
        parameter_policy="constant parameters folded into literal coefficients at import",
        source_parameters={s: float(v) for s, v in symbols.items() if s != compid},
    )


def _expression(text, parameters):
    # Bound parser work before ast.parse: deep or long inputs exhaust the stack.
    if len(text) > 1024:
        fail("expression length limit")
    expression_text = text.replace("^", "**")
    try:
        tree = ast.parse(expression_text, mode="eval")
    except (SyntaxError, ValueError, RecursionError, MemoryError):
        fail("invalid finite-network scalar expression")
    nodes = list(ast.walk(tree))
    if len(nodes) > 128:
        fail("expression work limit")

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and type(node.value) in [int, float]:
            return _number(ast.get_source_segment(expression_text, node))
        if isinstance(node, ast.Name) and node.id in parameters:
            return parameters[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub | ast.UAdd):
            return (
                -visit(node.operand)
                if isinstance(node.op, ast.USub)
                else visit(node.operand)
            )
        if isinstance(node, ast.BinOp):
            a = visit(node.left)
            b = visit(node.right)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Div) and b:
                return a / b
            if isinstance(node.op, ast.Pow) and b.denominator == 1 and -8 <= b <= 8:
                if (
                    max(a.numerator.bit_length(), a.denominator.bit_length())
                    * abs(int(b))
                    > 4096
                ):
                    fail("scalar power work limit")
                return a ** int(b)
        fail(
            "unsupported scalar expression; calls and species-dependent rates are unsupported"
        )

    try:
        value = visit(tree)
        if (
            value.numerator.bit_length() > 4096
            or value.denominator.bit_length() > 4096
            or not math.isfinite(float(value))
        ):
            fail("scalar expression range limit")
        return value
    except (OverflowError, ZeroDivisionError, RecursionError):
        fail("scalar arithmetic range failure")


def import_bng_net(raw, *, rate_convention):
    """Expanded numeric .net only; explicit binomial or falling_factorial convention.

    Returns a ``crn-craft.import/1.0`` report as a dictionary:
    ``report["model"]`` is the converted ``crn-craft.model/1.0`` model, beside
    the source format and its SHA-256. Unsupported input raises
    :class:`crn_craft.importers.ImportFailure`.
    """
    if rate_convention not in ["binomial", "falling_factorial"]:
        fail("explicit rate convention required")
    if len(raw) > 16 * 1024 * 1024:
        fail("oversized .net input")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        fail("UTF-8 network required")
    blocks = {}
    current = None
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("begin "):
            name = line[6:].strip()
            if (
                current
                or name in blocks
                or name not in ["parameters", "species", "reactions", "groups"]
            ):
                fail("unsupported/duplicate/nested .net block")
            blocks[name] = []
            current = name
        elif line.startswith("end "):
            if current != line[4:].strip():
                fail("mismatched block")
            current = None
        elif current:
            blocks[current].append(line)
        else:
            fail("content outside .net block")
    if current or not {"species", "reactions"} <= blocks.keys():
        fail("incomplete expanded network")
    parameters = {}
    for line in blocks.get("parameters", []):
        parts = line.split(maxsplit=2)
        if len(parts) != 3 or not _INDEX.fullmatch(parts[0]):
            fail("parameter index/name/value required")
        name = _identifier(parts[1])
        if name in parameters:
            fail("duplicate parameter")
        parameters[name] = _expression(parts[2], parameters)
    model = _model("ImportedBioNetGen")
    patterns = {}
    indices = set()
    for line in blocks["species"]:
        parts = line.split()
        if len(parts) != 3 or not _INDEX.fullmatch(parts[0]):
            fail("species index/pattern/initial count required")
        if parts[1].startswith("$"):
            # BioNetGen holds "$" species at a fixed amount; that is not a
            # dynamic species and has no supported mapping here.
            fail("fixed ($) species are unsupported")
        i = int(parts[0])
        count = _expression(parts[2], parameters)
        if (
            i <= 0
            or i in indices
            or count.denominator != 1
            or count < 0
            or count > 2**63 - 1
        ):
            fail("invalid species index/count")
        indices.add(i)
        patterns[f"S_{i}"] = parts[1]
        model["species"].append({"name": f"S_{i}", "initial_count": int(count)})
    reaction_ids = set()

    def side(text):
        if text == "0":
            return {}
        items = text.split(",")
        if not all(_INDEX.fullmatch(x) for x in items):
            fail("numeric expanded species indices required")
        ids = [int(x) for x in items]
        if any(i not in indices for i in ids):
            fail("unknown reaction species index")
        return dict(Counter(f"S_{i}" for i in ids))

    for line in blocks["reactions"]:
        parts = line.split(maxsplit=3)
        if len(parts) != 4 or not _INDEX.fullmatch(parts[0]):
            fail("reaction index/reactants/products/rate required")
        i = int(parts[0])
        if i <= 0 or i in reaction_ids:
            fail("duplicate/invalid reaction index")
        reaction_ids.add(i)
        a = side(parts[1])
        b = side(parts[2])
        k = _expression(parts[3], parameters)
        if k < 0:
            fail("negative reaction coefficient")
        if rate_convention == "falling_factorial" and k:
            if any(v > 1024 for v in a.values()):
                fail("converted coefficient exceeds the supported numeric range")
            k *= math.prod(math.factorial(v) for v in a.values())
        try:
            coefficient = float(k)
        except OverflowError:
            fail("converted coefficient exceeds f64")
        if not math.isfinite(coefficient) or (k > 0 and coefficient == 0):
            fail("converted coefficient overflows or underflows f64")
        model["reactions"].append(
            {
                "name": f"R_{i}",
                "reactants": a,
                "products": b,
                "rate": {"kind": "mass_action", "coefficient": coefficient},
            }
        )
    # Observable groups are metadata, not reaction-rate modifiers.
    return _report(
        raw,
        "bionetgen_expanded_net_subset",
        model,
        rate_convention=rate_convention,
        species_patterns=patterns,
        groups=blocks.get("groups", []),
        parameter_policy="constant expressions folded at import",
    )


def main():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=["sbml", "bng-net"])
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--convention", required=True)
    args = parser.parse_args()
    try:
        with args.source.open("rb") as source:
            raw = source.read(16 * 1024 * 1024 + 1)
        report = (
            import_sbml(raw, convention=args.convention)
            if args.format == "sbml"
            else import_bng_net(raw, rate_convention=args.convention)
        )
        with args.destination.open("x") as file:
            json.dump(report, file, indent=2)
            file.write("\n")
    except (OSError, ValueError, LookupError, RecursionError, MemoryError) as error:
        # ImportFailure is a ValueError; other parser failures get the same code.
        message = str(error) or type(error).__name__
        print(json.dumps({"success": False, "code": "CRN-E090", "message": message}))
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
