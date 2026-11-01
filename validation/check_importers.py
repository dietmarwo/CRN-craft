"""Importer semantics, combinatorial conversion and explicit rejection probes."""

import argparse
import json
from pathlib import Path

import crn_craft

# Validate new Python source against the installed native core before packaging.
crn_craft.__path__.insert(
    0, str(Path(__file__).resolve().parents[1] / "python/crn_craft")
)
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


def require(v, message):
    if not v:
        raise RuntimeError(message)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    model = import_sbml(SBML, convention="stochastic_count_propensity")["model"]
    analysis = crn_craft.analyze(model)
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
    require(
        model["reactions"][0]["rate"]["coefficient"] == 0.1,
        "SBML falling/binomial divisor",
    )
    rows.append({"check": "SBML_dimer_binomial_semantics", "passed": True})
    for convention, k in [("binomial", 0.1), ("falling_factorial", 0.2)]:
        result = import_bng_net(NET, rate_convention=convention)
        require(
            result["model"]["reactions"][0]["rate"]["coefficient"] == k,
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
        ]
    )
    for name, raw in invalid_sbml:
        try:
            import_sbml(raw, convention="stochastic_count_propensity")
        except ImportFailure:
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
        ]
    )
    for name, raw in invalid_net:
        try:
            import_bng_net(raw, rate_convention="binomial")
        except ImportFailure:
            rows.append({"check": "reject_" + name, "passed": True})
        else:
            raise RuntimeError("accepted " + name)
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
