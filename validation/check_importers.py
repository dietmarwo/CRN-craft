"""Importer semantics, combinatorial conversion and explicit rejection probes.

Since 0.7.0 (M5f) the importers return the converted `crn-craft.model/1.0`
document in a `crn-craft.import/1.0` report; the rate constants are read
from the converted rules, and a network without reactions, which a rules
model cannot hold, is rejected.
"""

import argparse
import json
from pathlib import Path
import time

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
    report = import_sbml(SBML, convention="stochastic_count_propensity")
    require(
        report["schema"] == "crn-craft.import/1.0"
        and report["source_format"] == "sbml_l3v2_core_subset",
        "import report contract",
    )
    model = report["model"]
    require(model["schema"] == "crn-craft.model/1.0", "converted model contract")
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
