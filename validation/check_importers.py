"""Importer semantics, combinatorial conversion and explicit rejection probes.

Since 0.7.0 (M5f) the importers return the converted `crn-craft.model/1.0`
document in a `crn-craft.import/1.0` report; the rate constants are read
from the converted rules, and a network without reactions, which a rules
model cannot hold, is rejected.
"""

import argparse
import json
from pathlib import Path
import re
import tempfile
import time
import zipfile

import crn_craft

# Validate new Python source against the installed native core before packaging.
crn_craft.__path__.insert(
    0, str(Path(__file__).resolve().parents[1] / "python/crn_craft")
)
from crn_craft._xml import _read_xml, _XmlRefused  # noqa: E402
from crn_craft.importers import ImportFailure, import_bng_net, import_sbml  # noqa: E402

SBML = b"""<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2"><model id="Dimer" substanceUnits="item" timeUnits="second"><listOfCompartments><compartment id="cell" size="1" constant="true"/></listOfCompartments><listOfSpecies><species id="A" compartment="cell" initialAmount="10" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/><species id="B" compartment="cell" initialAmount="0" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/></listOfSpecies><listOfParameters><parameter id="k" value="0.1" constant="true"/></listOfParameters><listOfReactions><reaction id="dimer" reversible="false"><listOfReactants><speciesReference species="A" stoichiometry="2" constant="true"/></listOfReactants><listOfProducts><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts><kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><divide/><apply><times/><ci>k</ci><ci>A</ci><apply><minus/><ci>A</ci><cn>1</cn></apply></apply><cn>2</cn></apply></math></kineticLaw></reaction></listOfReactions></model></sbml>"""
NET = b"""begin parameters
1 k 0.1
2 twice 2*k
end parameters
begin species
1 A(a) 10
2 B(b) 0
end species
begin reactions
1 1,1 2 k
2 2 1,1 twice
end reactions
begin groups
1 A_total 1,2*2
end groups
"""

SED = """<sedML xmlns="http://sed-ml.org/sed-ml/level1/version3" level="1" version="3">
<listOfModels><model id="m" source="model.xml" language="urn:sedml:language:sbml"/></listOfModels>
<listOfSimulations><uniformTimeCourse id="s" initialTime="0" outputStartTime="0" outputEndTime="1" numberOfPoints="2"><algorithm kisaoID="KISAO:0000029"/></uniformTimeCourse></listOfSimulations>
<listOfTasks><task id="t" modelReference="m" simulationReference="s"/></listOfTasks>
<listOfDataGenerators><dataGenerator id="g"><listOfVariables><variable id="x" taskReference="t" target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='A']"/></listOfVariables><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>x</ci></math></dataGenerator></listOfDataGenerators>
<listOfOutputs><report id="out"><listOfDataSets><dataSet id="A" dataReference="g"/></listOfDataSets></report></listOfOutputs></sedML>"""
MANIFEST = """<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">
<content location="study.sedml" format="http://identifiers.org/combine.specifications/sed-ml" master="true"/>
</omexManifest>"""
CONVENTION = "stochastic_count_propensity"
LIMIT = 16 << 20
UTF8_ONLY = "XML must be UTF-8 without a document type or entities"
DTD_OR_SIZE = "oversized XML or DTD/entities unsupported"
ONLY_UTF8 = "unreadable XML: only UTF-8 is supported"
TOO_DEEP = "XML nesting exceeds 256 levels"


def anchor(text, extra):
    """Add attribute text to the first element below the root."""
    return re.sub(r"<(model|content)\b", lambda m: m[0] + extra, text, count=1)


def pad(text, size):
    return text.encode() + b" " * (size - len(text.encode()))


def encoded(name):
    return lambda t: t.encode(name)


def declared(name, quote='"'):
    return (
        lambda t: f"<?xml version={quote}1.0{quote} encoding={quote}{name}{quote}?>{t}".encode()
    )


def prefixed(head):
    return lambda t: (head + t).encode()


def nested(levels):
    """Nest that many elements inside the root, as its last child."""

    def build(text):
        head, close, tail = text.rpartition("</")
        return (head + "<a>" * levels + "</a>" * levels + close + tail).encode()

    return build


DTD = '<!DOCTYPE a [<!ENTITY e "x">]>'
# case name -> (kind, builder). The kinds, with the refusal each one gets:
#   ok            read
#   type          not bytes (import_sbml only)
#   size          over 16 MiB
#   nul           NUL bytes (UTF-16, UTF-32)
#   dtd           a document type or entity declaration, anywhere in the bytes
#   enc_only      a declared encoding other than UTF-8 that the parser accepts
#   enc_parser    a declared encoding that the parser rejects as well
#   utf8          bytes that are not UTF-8
#   depth         elements nested deeper than 256, the root counted
#   syntax        not well-formed XML
# Not file contents: only import_sbml takes a value rather than a path.
IMPORT_ONLY = {"str", "none", "bytearray", "memoryview"}
GUARD_CASES = {
    "utf8": ("ok", encoded("utf-8")),
    "utf8_bom": ("ok", lambda t: b"\xef\xbb\xbf" + t.encode()),
    "declared_utf8_upper": ("ok", declared("UTF-8")),
    "declared_utf8_nodash": ("ok", declared("utf8")),
    "declared_none": ("ok", prefixed('<?xml version="1.0"?>')),
    "processing_instruction": ("ok", prefixed("<?target data?>")),
    "stylesheet_instruction": ("ok", prefixed('<?xml-stylesheet href="x.xsl"?>')),
    "comment": ("ok", prefixed("<!-- hello -->")),
    "cdata": ("ok", lambda t: t.replace(">", "><![CDATA[x]]>", 1).encode()),
    "predefined_entity": ("ok", lambda t: anchor(t, " name='&amp;&lt;'").encode()),
    "large_attribute": (
        "ok",
        lambda t: anchor(t, " name='" + "a" * (8 << 20) + "'").encode(),
    ),
    "size_limit_minus_one": ("ok", lambda t: pad(t, LIMIT - 1)),
    "size_limit": ("ok", lambda t: pad(t, LIMIT)),
    "size_limit_plus_one": ("size", lambda t: pad(t, LIMIT + 1)),
    "size_utf16": ("size", lambda t: t.encode("utf-16") + b" \x00" * LIMIT),
    "str": ("type", lambda t: t),
    "none": ("type", lambda t: None),
    "bytearray": ("ok", lambda t: bytearray(t.encode())),
    "memoryview": ("ok", lambda t: memoryview(t.encode())),
    "utf16_bom": ("nul", encoded("utf-16")),
    "utf16_le": ("nul", encoded("utf-16-le")),
    "utf16_be": ("nul", encoded("utf-16-be")),
    "utf32_bom": ("nul", encoded("utf-32")),
    "utf32_le": ("nul", encoded("utf-32-le")),
    "utf32_be": ("nul", encoded("utf-32-be")),
    "utf16_declared_dtd": (
        "nul",
        lambda t: ('<?xml version="1.0" encoding="UTF-16"?>' + DTD + t).encode(
            "utf-16-le"
        ),
    ),
    "nul_in_comment": ("nul", prefixed("<!--\x00-->")),
    "doctype": ("dtd", prefixed("<!DOCTYPE a>")),
    "doctype_lower": ("dtd", prefixed("<!doctype a>")),
    "doctype_mixed": ("dtd", prefixed("<!DocType a>")),
    "doctype_internal_subset": ("dtd", prefixed("<!DOCTYPE a [<!ELEMENT a ANY>]>")),
    "doctype_entity": ("dtd", prefixed(DTD)),
    "parameter_entity": ("dtd", prefixed('<!DOCTYPE a [<!ENTITY % p "x">%p;]>')),
    "external_entity": (
        "dtd",
        prefixed('<!DOCTYPE a [<!ENTITY e SYSTEM "file:///nonexistent/crn-craft">]>'),
    ),
    "external_subset": (
        "dtd",
        prefixed('<!DOCTYPE a SYSTEM "file:///nonexistent/crn-craft">'),
    ),
    "billion_laughs": (
        "dtd",
        lambda t: (
            '<!DOCTYPE l [<!ENTITY a "aaaaaaaaaa"><!ENTITY b "&a;&a;&a;&a;&a;&a;&a;&a;&a;&a;">'
            '<!ENTITY c "&b;&b;&b;&b;&b;&b;&b;&b;&b;&b;">]>' + anchor(t, " note='&c;'")
        ).encode(),
    ),
    "entity_declaration": ("dtd", prefixed("<!ENTITY e 'x'>")),
    "entity_declaration_lower": ("dtd", prefixed("<!entity e 'x'>")),
    "doctype_in_comment": ("dtd", prefixed("<!-- <!DOCTYPE a> -->")),
    "entity_in_comment": ("dtd", prefixed("<!-- <!ENTITY e 'x'> -->")),
    "doctype_in_cdata": (
        "dtd",
        lambda t: t.replace(">", "><![CDATA[<!DOCTYPE a>]]>", 1).encode(),
    ),
    "declared_utf7_doctype": (
        "dtd",
        lambda t: ('<?xml version="1.0" encoding="utf-7"?><!DOCTYPE a>' + t).encode(),
    ),
    "declared_latin1": ("enc_only", declared("ISO-8859-1")),
    "declared_latin1_single_quote": ("enc_only", declared("ISO-8859-1", "'")),
    "declared_latin1_after_bom": (
        "enc_only",
        lambda t: b"\xef\xbb\xbf" + declared("ISO-8859-1")(t),
    ),
    "declared_ascii": ("enc_only", declared("US-ASCII")),
    "declared_ascii_lower": ("enc_only", declared("ascii")),
    "declared_utf16": ("enc_parser", declared("UTF-16")),
    "declared_utf32": ("enc_parser", declared("UTF-32")),
    "declared_utf7": ("enc_parser", declared("utf-7")),
    "declared_ebcdic": ("enc_parser", declared("cp037")),
    "declared_unknown": ("enc_parser", declared("x-unknown")),
    "declared_empty": ("enc_parser", declared("")),
    "declared_latin1_bytes": (
        "enc_parser",
        lambda t: (
            '<?xml version="1.0" encoding="ISO-8859-1"?><!-- \xe9 -->' + t
        ).encode("latin-1"),
    ),
    "invalid_utf8": ("utf8", lambda t: anchor(t, " note='\xff'").encode("latin-1")),
    "binary": ("utf8", lambda t: bytes(range(1, 255)) * 8),
    "nesting_over_limit": ("depth", nested(256)),
    "nesting_far_over_limit": ("depth", nested(100_000)),
    "nesting_unclosed": ("depth", lambda t: nested(256)(t)[:-600]),
    "empty": ("syntax", lambda t: b""),
    "whitespace": ("syntax", lambda t: b"  \n"),
    "text": ("syntax", lambda t: b"hello"),
    "truncated": ("syntax", lambda t: t.encode()[:-12]),
    "undefined_entity": ("syntax", lambda t: anchor(t, " note='&e;'").encode()),
    "control_character": ("syntax", lambda t: anchor(t, " note='\x01'").encode()),
    "uppercase_declaration": ("syntax", prefixed('<?XML version="1.0"?>')),
    "space_in_declaration_open": ("syntax", prefixed("<! DOCTYPE a>")),
}


def guard_expectation(kind, channel):
    """The refusal expected from one entry point: (text, exact) or None to read."""
    if kind == "ok":
        return None
    if channel == "import":
        if kind in {"size", "nul", "dtd"}:
            return DTD_OR_SIZE, True
        if kind == "enc_only":
            return ONLY_UTF8, True
        if kind == "type":
            return "unreadable XML: XML must be given as bytes", True
        if kind == "enc_parser":
            return "unreadable XML: ", False
        if kind == "utf8":
            return "unreadable XML: 'utf-8' codec can't decode", False
        if kind == "depth":
            return "unreadable XML: " + TOO_DEEP, True
        return "unreadable XML: ", False
    # read_experiment; a manifest failure is wrapped as an archive failure
    wrap = "invalid OMEX archive: " if channel == "manifest" else ""
    if kind == "size":
        if channel == "manifest":
            return wrap + "XML document exceeds 16 MiB", True
        return f"file exceeds its size budget: {channel}", True
    if kind in {"nul", "dtd", "enc_only", "enc_parser"}:
        return wrap + UTF8_ONLY, True
    if kind == "utf8":
        return wrap + "invalid XML: 'utf-8' codec can't decode", False
    if kind == "depth":
        return wrap + "invalid XML: " + TOO_DEEP, True
    return wrap + "invalid XML: ", False


REASONS = {
    "type": "type",
    "size": "size",
    "nul": "nul",
    "dtd": "dtd",
    "enc_only": "encoding",
    "enc_parser": "encoding",
    "utf8": "utf8",
    "depth": "depth",
    "syntax": "syntax",
}


def guard_reason(data):
    """The reason and parser detail with which the guard itself refuses data."""
    try:
        _read_xml(data, LIMIT)
    except _XmlRefused as refused:
        return refused.reason, refused.detail
    return None, None


def refusal(read, *args, **kwargs):
    try:
        read(*args, **kwargs)
    except ImportFailure as error:
        return f"{error.code}: {error}"
    return None


def guard_rows():
    """Every guard condition through import_sbml and read_experiment.

    One transformation is applied to an SBML document read by import_sbml, to
    a SED-ML file, to a model file and to the manifest of an archive, so all
    four reading sites of the package meet the same bytes.
    """
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        for name, (kind, build) in GUARD_CASES.items():
            sites = {
                "import": SBML.decode(),
                "study.sedml": SED,
                "model.xml": SBML.decode(),
                "manifest": MANIFEST,
            }
            for channel, text in sites.items():
                if name in IMPORT_ONLY and channel != "import":
                    continue
                data = build(text)
                if channel == "import":
                    reason, detail = guard_reason(data)
                    require(
                        reason == REASONS.get(kind),
                        f"{name}: the guard says {reason!r}",
                    )
                    require(
                        (detail is None) == (kind != "enc_parser"),
                        f"{name}: parser detail {detail!r}",
                    )
                    result = refusal(import_sbml, data, convention=CONVENTION)
                elif channel == "manifest":
                    archive = temp / "case.omex"
                    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as arc:
                        arc.writestr("manifest.xml", data)
                        arc.writestr("study.sedml", SED)
                        arc.writestr("model.xml", SBML)
                    result = refusal(
                        crn_craft.read_experiment, archive, convention=CONVENTION
                    )
                else:
                    folder = temp / "files"
                    folder.mkdir(exist_ok=True)
                    (folder / "study.sedml").write_bytes(
                        data if channel == "study.sedml" else SED.encode()
                    )
                    (folder / "model.xml").write_bytes(
                        data if channel == "model.xml" else SBML
                    )
                    result = refusal(
                        crn_craft.read_experiment,
                        folder / "study.sedml",
                        convention=CONVENTION,
                    )
                expected = guard_expectation(kind, channel)
                if expected is None:
                    require(
                        result is None, f"{name} through {channel} refused: {result}"
                    )
                else:
                    text_, exact = expected
                    want = "CRN-E090: " + text_
                    require(
                        result is not None
                        and (result == want if exact else result.startswith(want))
                        # the parser's own message, not the guard's
                        and (kind != "enc_parser" or ONLY_UTF8 not in result),
                        f"{name} through {channel}: {result!r} not {want!r}",
                    )
            rows.append({"check": "xml_guard_" + name, "passed": True})
    # The nesting limit itself: 256 levels are read, 257 are not, and a
    # declared encoding is reported before the depth.
    at_limit = _read_xml(b"<a>" * 256 + b"</a>" * 256, LIMIT)
    require(
        sum(1 for _ in at_limit.iter()) == 256
        and guard_reason(b"<a>" * 257 + b"</a>" * 257) == ("depth", None)
        and guard_reason(b"<r>" + b"<a><b/></a>" * 1000 + b"</r>") == (None, None)
        and guard_reason(declared("ISO-8859-1")("<a>" * 257 + "</a>" * 257))
        == ("encoding", None),
        "the XML nesting limit is not 256 levels",
    )
    rows.append({"check": "xml_guard_nesting_limit", "passed": True})
    return rows


FORBIDDEN_XML = (
    "xml.etree",
    "xml.dom",
    "xml.sax",
    "expat",
    "lxml",
    "fromstring",
    "iterparse",
    "XMLParser",
    "ET.parse",
    "minidom",
    "defusedxml",
)


def xml_entry_points():
    """Only the guard module may name an XML parser."""
    package = Path(__file__).resolve().parents[1] / "python/crn_craft"
    guard = package / "_xml.py"
    require(guard.is_file(), "the XML guard module is missing")
    found = []
    for path in sorted(package.rglob("*.py")):
        if path == guard:
            continue
        source = path.read_text(encoding="utf-8")
        found += [f"{path.name}: {word}" for word in FORBIDDEN_XML if word in source]
    require(not found, "XML parsing outside the guard: " + ", ".join(found))
    require(
        "XMLParser" in guard.read_text(encoding="utf-8"), "the guard parses nothing"
    )
    return [{"check": "xml_parsing_only_in_guard", "passed": True}]


def require(v, message):
    if not v:
        raise RuntimeError(message)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    report = import_sbml(SBML, convention="stochastic_count_propensity")
    require(
        report["schema"] == "crn-craft.import/1.0"
        and report["source_format"] == "sbml_l3v2_core_subset",
        "import report contract",
    )
    model = report["model"]
    require(model["schema"] == "crn-craft.model/1.0", "converted model contract")
    analysis = crn_craft.Model.from_dict(model).analyze()
    require(analysis["success"], "imported model analysis")

    # Inspect pure compiled propensity through the public analysis report.
    def find_rates(value):
        if isinstance(value, dict):
            for k, v in value.items():
                if k == "propensities":
                    return v
                r = find_rates(v)
                if r is not None:
                    return r
        return None

    rates = find_rates(analysis)
    if rates is not None:
        require(rates == [4.5], f"expected binomial dimer propensity: {rates}")
    dimer = model["rules"][0]
    require(
        float(dimer["rate"]) == 0.1
        and [p["type"] for p in dimer["reactants"]] == ["A", "A"],
        "SBML falling/binomial divisor",
    )
    rows.append({"check": "SBML_dimer_binomial_semantics", "passed": True})
    for convention, k in [("binomial", 0.1), ("falling_factorial", 0.2)]:
        result = import_bng_net(NET, rate_convention=convention)
        require(
            float(result["model"]["rules"][0]["rate"]) == k,
            "explicit BNG conversion",
        )
        require(result["species_patterns"]["S_1"] == "A(a)", "species provenance")
        rows.append({"check": "BNG_" + convention, "passed": True})
    invalid_sbml = [
        ("events", SBML.replace(b"</model>", b"<listOfEvents/></model>")),
        (
            "concentration",
            SBML.replace(b'initialAmount="10"', b'initialConcentration="10"'),
        ),
        ("noninteger", SBML.replace(b'initialAmount="10"', b'initialAmount="1.5"')),
        (
            "variable_compartment",
            SBML.replace(b'constant="true"', b'constant="false"', 1),
        ),
        ("reversible", SBML.replace(b'reversible="false"', b'reversible="true"')),
        (
            "deterministic_homodimer",
            SBML.replace(b"<apply><minus/><ci>A</ci><cn>1</cn></apply>", b"<ci>A</ci>"),
        ),
        ("unsupported_symbol", SBML.replace(b"<ci>k</ci>", b"<ci>missing</ci>")),
        (
            "foreign_attribute",
            SBML.replace(b"<model id=", b'<model conversionFactor="x" id='),
        ),
        (
            "species_concentration_semantics",
            SBML.replace(
                b'hasOnlySubstanceUnits="true"', b'hasOnlySubstanceUnits="false"'
            ),
        ),
        ("units", SBML.replace(b'substanceUnits="item"', b'substanceUnits="mole"')),
        ("DTD", b"<!DOCTYPE sbml>" + SBML),
        (
            "local_parameters",
            SBML.replace(b"<kineticLaw>", b"<kineticLaw><listOfLocalParameters/>"),
        ),
    ]
    invalid_sbml.extend(
        [
            ("numeric_overflow", SBML.replace(b'value="0.1"', b'value="1e400"')),
            ("undefined_initial_amount", SBML.replace(b' initialAmount="10"', b"")),
            ("undefined_parameter", SBML.replace(b' value="0.1"', b"")),
            ("undefined_volume", SBML.replace(b' size="1"', b"")),
            ("undefined_reversibility", SBML.replace(b' reversible="false"', b"")),
            ("UTF16_XML", SBML.decode().encode("utf-16")),
            # L3V2 has no defaults for these: an absent value is undefined.
            (
                "fast_attribute",
                SBML.replace(b'reversible="false"', b'reversible="false" fast="false"'),
            ),
            ("missing_stoichiometry", SBML.replace(b' stoichiometry="2"', b"")),
            (
                "missing_reference_constant",
                SBML.replace(
                    b'stoichiometry="2" constant="true"', b'stoichiometry="2"'
                ),
            ),
            ("missing_substance_units", SBML.replace(b' substanceUnits="item"', b"")),
            ("missing_time_units", SBML.replace(b' timeUnits="second"', b"")),
            (
                "minute_time_units",
                SBML.replace(b'timeUnits="second"', b'timeUnits="minute"'),
            ),
            ("fraction_literal", SBML.replace(b'value="0.1"', b'value="1/10"')),
            ("underscore_literal", SBML.replace(b'value="0.1"', b'value="0_1"')),
            (
                "huge_exponent_literal",
                SBML.replace(b'value="0.1"', b'value="1e-16000000"'),
            ),
            (
                "reaction_id_collision",
                SBML.replace(b'reaction id="dimer"', b'reaction id="A"'),
            ),
            (
                "boundary_condition",
                SBML.replace(
                    b'boundaryCondition="false"', b'boundaryCondition="true"', 1
                ),
            ),
            (
                "duplicate_reference",
                SBML.replace(
                    b'<speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>',
                    b'<speciesReference species="B" stoichiometry="1" constant="true"/><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>',
                ),
            ),
            (
                "noninteger_stoichiometry",
                SBML.replace(b'stoichiometry="2"', b'stoichiometry="1.5"'),
            ),
            ("unknown_encoding", b'<?xml version="1.0" encoding="x-unknown"?>' + SBML),
            # A rules model needs a rule (CRN-R016 from the converter).
            (
                "no_reactions",
                SBML[: SBML.index(b"<listOfReactions>")]
                + SBML[
                    SBML.index(b"</listOfReactions>") + len(b"</listOfReactions>") :
                ],
            ),
        ]
    )
    for name, raw in invalid_sbml:
        started = time.perf_counter()
        try:
            import_sbml(raw, convention="stochastic_count_propensity")
        except ImportFailure:
            # Rejection must be cheap: no unbounded literal or parser work.
            elapsed = time.perf_counter() - started
            if elapsed > 1.0:
                raise RuntimeError(f"rejecting {name} took {elapsed:.1f} s") from None
            rows.append({"check": "reject_" + name, "passed": True})
        else:
            raise RuntimeError("accepted " + name)
    invalid_net = [
        ("unknown_block", NET + b"begin functions\nend functions\n"),
        ("species_reference", NET.replace(b"1 1,1 2 k", b"1 1,99 2 k")),
        ("expression_call", NET.replace(b"2*k", b'__import__("os")')),
        ("negative_count", NET.replace(b"A(a) 10", b"A(a) -1")),
        ("negative_rate", NET.replace(b"1 1,1 2 k", b"1 1,1 2 -k")),
        ("duplicate_reaction", NET.replace(b"2 2 1,1 twice", b"1 2 1,1 twice")),
    ]
    invalid_net.extend(
        [
            ("coefficient_underflow", NET.replace(b"1 1,1 2 k", b"1 1,1 2 1e-400")),
            ("coefficient_overflow", NET.replace(b"1 1,1 2 k", b"1 1,1 2 1e400")),
            ("fixed_species", NET.replace(b"1 A(a) 10", b"1 $A(a) 10")),
            ("non_ascii_index", NET.replace(b"1 A(a) 10", "\u00b9 A(a) 10".encode())),
            (
                "long_expression",
                NET.replace(b"2 twice 2*k", b"2 twice " + b"+".join([b"k"] * 10000)),
            ),
            (
                "deep_expression",
                NET.replace(b"2 twice 2*k", b"2 twice " + b"-" * 200000 + b"k"),
            ),
        ]
    )
    # The falling-factorial conversion multiplies by nu!, not nu: order three.
    trimer = NET.replace(b"1 1,1 2 k", b"1 1,1,1 2 k")
    converted = import_bng_net(trimer, rate_convention="falling_factorial")
    require(
        abs(float(converted["model"]["rules"][0]["rate"]) - 0.6) < 1e-15
        and converted["model"]["rules"][0]["reactants"][0]["count"] == 3,
        "falling factorial order three",
    )
    rows.append({"check": "BNG_falling_factorial_order_three", "passed": True})
    for name, raw in invalid_net:
        try:
            import_bng_net(raw, rate_convention="binomial")
        except ImportFailure:
            rows.append({"check": "reject_" + name, "passed": True})
        else:
            raise RuntimeError("accepted " + name)
    rows.extend(guard_rows())
    rows.extend(xml_entry_points())
    args.output.write_text(
        json.dumps(
            {"schema": "crn-craft.import-validation/0.1", "passed": True, "rows": rows},
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "checks": len(rows)}))


if __name__ == "__main__":
    main()
