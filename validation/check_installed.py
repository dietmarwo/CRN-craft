"""Acceptance from an installed wheel, independent of repository imports.

Run it with the environment's Python. The console scripts are taken from that
environment's scripts directory, never from PATH, so a native `crn-craft`
binary elsewhere on PATH cannot stand in for the installed entry point.
"""

import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import shutil
import subprocess
import sysconfig
import tempfile
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import crn_craft as c
from crn_craft import rules
from crn_craft.importers import import_bng_net
from crn_craft.observables import oscillation, sampled_first_passage

SBML = (
    b'<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" '
    b'version="2"><model id="Decay" substanceUnits="item" timeUnits="second">'
    b'<listOfCompartments><compartment id="cell" size="1" constant="true"/>'
    b'</listOfCompartments><listOfSpecies><species id="A" compartment="cell" '
    b'initialAmount="10" hasOnlySubstanceUnits="true" boundaryCondition="false" '
    b'constant="false"/></listOfSpecies><listOfParameters><parameter id="k" '
    b'value="0.5" constant="true"/></listOfParameters><listOfReactions>'
    b'<reaction id="decay" reversible="false"><listOfReactants><speciesReference '
    b'species="A" stoichiometry="1" constant="true"/></listOfReactants>'
    b'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/>'
    b"<ci>k</ci><ci>A</ci></apply></math></kineticLaw></reaction>"
    b"</listOfReactions></model></sbml>"
)


def script(name):
    """The console script installed with this interpreter's package."""
    found = shutil.which(name, path=sysconfig.get_path("scripts"))
    if found is None:
        raise RuntimeError(f"{name} is not installed next to this interpreter")
    return found


def require(v, m):
    if not v:
        raise RuntimeError(m)


def python_spelling(cargo):
    """Cargo's `0.3.0-dev.1` is Python's `0.3.0.dev1`; releases are equal."""
    release, _, pre = cargo.partition("-")
    if not pre:
        return release
    kind, _, number = pre.partition(".")
    return (
        release + {"dev": ".dev", "alpha": "a", "beta": "b", "rc": "rc"}[kind] + number
    )


def three_class_model():
    """models/rules/mtdna-three-class.rules.json, without its description.

    Written out here so that the check needs no file from the checkout.
    """

    def cls(clone, s, wt):
        return {"clone": clone, "s": s, "wt": wt}

    def pattern(clone):
        return [{"label": "x", "type": "mtDNA", "where": f"x.clone == {clone}"}]

    copy = {"type": "mtDNA", "copy": "x"}
    vacancy = "max(0, K - N) / D * "

    def rule(name, clone, rate, products):
        return {
            "name": name,
            "reactants": pattern(clone),
            "rate": rate,
            "products": products,
        }

    return {
        "schema": "crn-craft.rules-model/0.1",
        "name": "mtdna_three_class",
        "parameters": {
            "K": 1000,
            "delta": 0.0028881132523331052,
            "sw": 1.0,
            "sm": 1.5,
            "mu": 3e-05,
            "fracAdv": 0.2,
            "w0": 1000,
        },
        "types": [{"name": "mtDNA", "attributes": ["clone", "s", "wt"]}],
        "initial": [{"type": "mtDNA", "attributes": cls(0, "sw", 1), "count": "w0"}],
        "observables": {
            "N": {"type": "mtDNA"},
            "wild_type": {"type": "mtDNA", "where": "x.wt == 1"},
            "mutants": {"type": "mtDNA", "where": "x.wt == 0"},
            "advantaged_mutants": {"type": "mtDNA", "where": "x.clone == 2"},
            "D": {"type": "mtDNA", "weight": "x.wt + sm * (1 - x.wt)"},
            "neutral_mutants": {"type": "mtDNA", "where": "x.clone == 1"},
        },
        "rules": [
            rule("degradation_wt", 0, "delta", []),
            rule("replication_wt", 0, vacancy + "(sw * (1 - mu))", [copy, copy]),
            rule(
                "mutation_mt",
                0,
                vacancy + "sw * mu * (1 - fracAdv)",
                [copy, {"type": "mtDNA", "attributes": cls(1, "sw", 0)}],
            ),
            rule(
                "mutation_mtAdv",
                0,
                vacancy + "sw * mu * fracAdv",
                [copy, {"type": "mtDNA", "attributes": cls(2, "sm", 0)}],
            ),
            rule("degradation_mt", 1, "delta", []),
            rule("replication_mt", 1, vacancy + "(sw)", [copy, copy]),
            rule("degradation_mtAdv", 2, "delta", []),
            rule("replication_mtAdv", 2, vacancy + "(sm)", [copy, copy]),
        ],
    }


def check_event_log_default_budget():
    """The event-log request from the 0.2.0 feedback, at the default budget.

    Three-class model, 400 cells for three years, up to 50 logged events of
    the two mutation rules per cell, default max_output_bytes (64 MiB).
    Release 0.2.0 rejected it before running with CRN-R020 (observation
    reservation exceeds budget).
    """
    names = ["mutation_mt", "mutation_mtAdv"]
    try:
        run = rules.ensemble(
            three_class_model(),
            end=25_920.0,
            report_times=[0.0, 25_920.0],
            replicates=400,
            seed=4,
            threads=2,
            method="static_direct",
            event_log={"max_entries": 50, "rules": names},
        )
    except c.CrnError as error:
        raise RuntimeError(
            f"the reported event-log request is rejected: {error}"
        ) from None
    result = run.simulate()
    log = result.event_columns()  # replicate x entry x variable, NaN-padded
    frame = result.event_frame()
    omitted = sum(o["omitted_log_entries"] for o in result.observations)
    mutations = result.rule_events[names]
    counts = mutations.sum(axis=1)
    logged = log["length"].to_pandas().reindex(counts.index, fill_value=0)
    fired = log["rule_name"].values[log["rule"].values >= 0]
    require(
        result.metadata["termination"] == "completed"
        and bool((result.records["termination"] == "completed").all())
        and len(result.observations) == 400
        and log.dims == ("replicate", "entry", "variable")
        and log.sizes["replicate"] == 400
        and 0 < log.sizes["entry"] <= 50
        and set(fired) <= set(names)
        and set(frame["rule_name"].unique()) <= set(names)
        and omitted == 0
        and len(fired) == len(frame) == int(mutations.to_numpy().sum())
        and len(frame) > 0
        and np.array_equal(logged.to_numpy(), counts.to_numpy()),
        "the reported event log records every mutation of every cell",
    )


def check_builder():
    """rules.Model: a string is a whole propensity, a number mass action."""
    model = rules.Model("birth_death")
    model.add_parameter("k", 10.0).add_species("A")
    model.add_reaction("k", [], ["A"], name="birth")
    model.add_reaction(1.0, ["A"], [], name="death")
    hand = {
        "schema": "crn-craft.rules-model/0.3",
        "name": "birth_death",
        "parameters": {"k": 10.0},
        "types": [{"name": "A"}],
        "initial": [],
        "observables": {"A": {"type": "A"}},
        "rules": [
            {
                "name": "birth",
                "reactants": [],
                "propensity": "k",
                "products": [{"type": "A"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "A"}],
                "rate": "1.0",
                "products": [],
            },
        ],
    }
    require(
        model.to_dict() == hand and json.loads(model.to_json()) == hand,
        "the builder writes the hand-written rules-model/0.3 document",
    )
    times = [0.0, 1.0, 5.0]
    # Birth-death from zero: Poisson with mean k (1 - exp(-t)).
    exact = 10.0 * (1.0 - np.exp(-np.array(times)))
    results = {}
    for method in ("static_direct", None):
        results[method] = rules.ensemble(
            model,
            end=5.0,
            report_times=times,
            replicates=1000,
            seed=11,
            threads=2,
            method=method,
        ).simulate()
        mean = results[method].statistics["mean"]["A"].to_numpy()
        error = results[method].statistics["mean_standard_error"]["A"].to_numpy()
        require(
            mean[0] == 0.0
            and bool(np.all(np.abs(mean[1:] - exact[1:]) <= 5 * error[1:])),
            f"builder birth-death matches the Poisson law ({method or 'tree'})",
        )
    same = rules.ensemble(
        hand,
        end=5.0,
        report_times=times,
        replicates=1000,
        seed=11,
        threads=2,
        method="static_direct",
    ).simulate()
    require(
        same.metadata["model_sha256"]
        == results["static_direct"].metadata["model_sha256"]
        and same.records.equals(results["static_direct"].records)
        and results["static_direct"].rule_names == ("birth", "death"),
        "a builder model runs exactly as its JSON",
    )


def raises(error, call, *args, **kwargs):
    """The exception `call` raises, required to be an `error`."""
    try:
        call(*args, **kwargs)
    except error as raised:
        return raised
    raise RuntimeError(f"{getattr(call, '__name__', call)} accepted {args} {kwargs}")


def check_builder_mappings():
    """add_species, add_parameter and add_parameters take mappings (0.5.0)."""
    single = rules.Model("sir")
    single.add_species("S", 999).add_species("I", 1).add_species("R")
    single.add_parameter("beta", 2e-4).add_parameter("gamma", 0.01)
    many = rules.Model("sir")
    many.add_species({"S": 999, "I": 1, "R": 0})
    many.add_parameter({"beta": 2e-4}).add_parameters({"gamma": 0.01})
    for model in (single, many):
        model.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
        model.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    require(
        many.to_json() == single.to_json()
        and many.species == ("S", "I", "R")
        and list(many.parameters) == ["beta", "gamma"],
        "a mapping adds every entry, in order, as one call each does",
    )
    require(
        rules.Model.add_parameters is rules.Model.add_parameter
        and rules.Model("x").add_parameters("k", 1.0).parameters == {"k": 1.0},
        "add_parameters is add_parameter",
    )
    # A species a reaction introduced keeps its place and takes the count.
    late = rules.Model("late").add_reaction(1.0, ["A"], ["B"])
    late.add_species({"B": 3, "C": 2})
    require(
        late.species == ("A", "B", "C")
        and late.to_dict()["initial"]
        == [
            {"type": "B", "count": 3},
            {"type": "C", "count": 2},
        ]
        and rules.Model("e").add_species({}).add_parameter({}).to_json()
        == rules.Model("e").to_json(),
        "mapping entries follow the one-at-a-time rules; an empty mapping adds none",
    )
    for call in (
        lambda m: m.add_species({"S": 1}, 5),
        lambda m: m.add_species({"S": 1}, initial=5),
        lambda m: m.add_parameter({"k": 1.0}, 2.0),
        lambda m: m.add_parameter({"k": 1.0}, value=2.0),
        lambda m: m.add_parameters({"k": 1.0}, value=2.0),
        lambda m: m.add_parameter("k"),
    ):
        model = rules.Model("mixed")
        raises(TypeError, call, model)
        require(model.to_json() == rules.Model("mixed").to_json(), "mixed arguments")
    # Invalid entries fail as they do one at a time, and leave the model as
    # it was.
    for build, bad, single_call in (
        (
            lambda m: m.add_species({"S": 1, "I": -1}),
            ValueError,
            lambda m: m.add_species("I", -1),
        ),
        (
            lambda m: m.add_species({"S": 1, "min": 2}),
            ValueError,
            lambda m: m.add_species("min", 2),
        ),
        (
            lambda m: m.add_species({"S": 1, 3: 2}),
            TypeError,
            lambda m: m.add_species(3, 2),
        ),
        (
            lambda m: m.add_species({"S": 1, "I": 1.5}),
            TypeError,
            lambda m: m.add_species("I", 1.5),
        ),
        (
            lambda m: m.add_species({"S": 1, "k": 1}),
            ValueError,
            lambda m: m.add_species("k", 1),
        ),
        (
            lambda m: m.add_parameter({"a": 1.0, "b": float("nan")}),
            ValueError,
            lambda m: m.add_parameter("b", float("nan")),
        ),
        (
            lambda m: m.add_parameter({"a": 1.0, "b": True}),
            TypeError,
            lambda m: m.add_parameter("b", True),
        ),
        (
            lambda m: m.add_parameters({"a": 1.0, "k": 2.0}),
            ValueError,
            lambda m: m.add_parameter("k", 2.0),
        ),
    ):

        def fresh():
            return rules.Model("strict").add_parameter("k", 1.0)

        model = fresh()
        mapped = raises(bad, build, model)
        alone = raises(bad, single_call, fresh())
        require(
            str(mapped) == str(alone) and model.to_json() == fresh().to_json(),
            f"mapping validation equals one call: {alone}",
        )


# rules.Model("odes") with two continuous variables, as release 0.4.0's
# add_continuous wrote it (to_json(indent=None)).
ODES_040 = (
    '{"schema": "crn-craft.rules-model/0.3", "name": "odes", "parameters": '
    '{"g": 0.01, "m0": 0.001}, "types": [{"name": "A"}], "initial": [{"type": '
    '"A", "count": 5}], "observables": {"A": {"type": "A"}}, "rules": [{"name": '
    '"death", "reactants": [{"label": "x", "type": "A"}], "propensity": '
    '"mortality * A", "products": []}], "continuous": [{"name": "mortality", '
    '"initial": "m0", "derivative": "g * mortality"}, {"name": "atp", "initial": '
    '450.0, "derivative": "0.1 * A - 0.2 * atp", "lower_bound": 0.0, "on_bound": '
    '"clamp"}]}\n'
)


def odes_model(add, first=None, second=None):
    """The ODES_040 model, its variables declared with method `add`."""
    model = rules.Model("odes")
    model.add_parameter({"g": 0.01, "m0": 1e-3}).add_species("A", 5)
    getattr(model, add)("mortality", "m0", "g * mortality", **(first or {}))
    getattr(model, add)(
        "atp",
        450.0,
        "0.1 * A - 0.2 * atp",
        lower_bound=0,
        on_bound="clamp",
        **(second or {}),
    )
    model.add_reaction("mortality * A", ["A"], [], name="death")
    return model


def check_builder_stop():
    """set_stop writes the model's stop condition; None removes it."""
    sir = rules.Model("sir").add_species({"S": 50, "I": 1, "R": 0})
    sir.add_parameter({"beta": 0.01, "gamma": 0.1})
    sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
    sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    plain = sir.to_json()
    require(
        sir.set_stop("I == 0") is sir
        and list(sir.to_dict())[-2:] == ["rules", "stop"]
        and sir.to_dict()["stop"] == "I == 0",
        "set_stop writes stop after the rules and returns the model",
    )
    result = rules.ensemble(
        sir, end=1000.0, report_times=[0.0, 1000.0], replicates=40, seed=3
    ).simulate()
    require(
        (result.records["termination"] == "stopped").all()
        and result.metadata["stopped_replicates"] == 40
        and (result.records["stop_time"] < 1000.0).all(),
        "every SIR replicate stops when the infection dies out",
    )
    require(
        sir.set_stop(None).to_json() == plain and "stop" not in sir.to_dict(),
        "set_stop(None) removes the stop and restores the bytes",
    )
    raises(TypeError, sir.set_stop, 0)
    raises(ValueError, sir.set_stop, "  ")
    unknown = rules.Model("u").add_species({"A": 1}).add_reaction(0.1, ["A"], [])
    unknown.set_stop("B == 0")
    error = raises(c.CrnError, rules.ensemble, unknown, end=1.0, report_times=[1.0])
    require(
        error.diagnostics[0]["path"] == "/stop",
        "a stop that reads an undeclared name is refused at /stop",
    )


def check_builder_odes():
    """add_ode, its alias add_ODE, independent=True and add_continuous."""
    plain = odes_model("add_ode")
    explicit = odes_model("add_ode", {"independent": False}, {"independent": False})
    require(
        plain.to_json(indent=None) == ODES_040
        and explicit.to_json(indent=None) == ODES_040
        and rules.model_bytes(plain.to_dict())
        == rules.model_bytes(json.loads(ODES_040))
        and plain.to_json() == json.dumps(json.loads(ODES_040), indent=2) + "\n",
        "add_ode without independent variables writes the bytes of 0.4.0",
    )
    require(rules.Model.add_ODE is rules.Model.add_ode, "add_ODE is add_ode")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        old = odes_model("add_continuous")
    require(
        old.to_json(indent=None) == ODES_040
        and [w.category for w in caught] == [DeprecationWarning] * 2
        and all("add_ode" in str(w.message) for w in caught)
        and all(Path(w.filename).resolve() == Path(__file__).resolve() for w in caught),
        "add_continuous warns, names add_ode, points at the caller, writes 0.4.0",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        raises(
            TypeError,
            rules.Model("x").add_continuous,
            "y",
            0.0,
            "1",
            independent=True,
        )
    # independent=True marks that variable and selects rules-model/0.4.
    for first, second in (({"independent": True}, None), (None, {"independent": True})):
        model = odes_model("add_ode", first, second)
        expected = json.loads(ODES_040)
        expected["schema"] = "crn-craft.rules-model/0.4"
        position = 0 if first else 1
        variable = expected["continuous"][position]
        expected["continuous"][position] = {
            "name": variable.pop("name"),
            "initial": variable.pop("initial"),
            "derivative": variable.pop("derivative"),
            "independent": True,
            **variable,
        }
        require(
            model.to_json(indent=None) == json.dumps(expected) + "\n"
            and model.to_dict() == expected
            and "independent" not in model.to_dict()["continuous"][1 - position],
            "independent=True writes independent: true and rules-model/0.4",
        )
    for bad in (1, 0, "yes", None, 1.0):
        model = rules.Model("x")
        raises(TypeError, model.add_ode, "y", 0.0, "1", independent=bad)
        require(
            model.to_json() == rules.Model("x").to_json(),
            f"independent={bad!r} is refused and declares nothing",
        )


def check_pandas_results(result):
    """Statistics by name, NumPy access and categorical terminations."""
    meta = result.metadata
    observables = list(meta["observables"])
    for key in ("mean", "sample_variance", "mean_standard_error"):
        frame = result.statistics[key]
        require(
            isinstance(frame, pd.DataFrame)
            and frame.index.name == "time"
            and list(frame.index) == meta["report_times"]
            and list(frame.columns) == observables
            and np.array_equal(
                frame.to_numpy(),
                np.array(meta["statistics"][key], dtype=np.float64).reshape(
                    frame.shape
                ),
            ),
            f"statistics[{key!r}] is the raw matrix by time and observable",
        )
    termination = result.records["termination"]
    codes = termination.cat.codes.to_numpy()
    require(
        list(termination.cat.categories) == list(rules.TERMINATIONS)
        and bool(np.all(codes == rules.TERMINATIONS.index("completed")))
        and isinstance(result.records["stop_time"].to_numpy(), np.ndarray),
        "termination codes and NumPy columns stay available",
    )


def check_hybrid():
    """Tiny hybrid runs: x' = 2 drives immigration; an independent ODE."""
    ramp = rules.Model("ramp")
    ramp.add_species("A")
    ramp.add_ode("x", 0.0, "2")
    ramp.add_reaction("x", [], ["A"], name="immigration")
    try:
        rules.ensemble(ramp, end=1.0, report_times=[0.0, 1.0], method="static_direct")
    except c.CrnError as error:
        require(error.diagnostics[0]["code"] == "CRN-R020", "hybrid-only refusal")
    else:
        raise RuntimeError("static_direct accepted a continuous variable")
    times = [0.0, 1.5, 3.0]
    result = rules.ensemble(
        ramp,
        end=3.0,
        report_times=times,
        replicates=400,
        seed=13,
        threads=2,
        method="hybrid",
    ).simulate()
    meta = result.metadata
    mean = result.statistics["mean"]
    error = result.statistics["mean_standard_error"]["A"].to_numpy()
    # x = 2 t, and the count at t is Poisson with mean t^2.
    exact = np.array(times) ** 2
    require(
        meta["schema"] == "crn-craft.rules-result/0.3"
        and meta["method"] == {"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}
        and meta["continuous"] == ["x"]
        and meta["sampler"] == "chacha8-u53/rules-approximate-hybrid-dopri54/1"
        and list(mean.columns) == ["A", "x"]
        and bool((result.records["termination"] == "completed").all())
        and np.allclose(mean["x"], [0.0, 3.0, 6.0], rtol=1e-6, atol=1e-9)
        and mean["A"].iloc[0] == 0.0
        and bool(np.all(np.abs(mean["A"].to_numpy()[1:] - exact[1:]) <= 5 * error[1:])),
        "installed hybrid integrates x and matches the Poisson law",
    )
    # An independent ODE (rules-model/0.4): Gompertz mortality m' = g m, so
    # each of 20 molecules survives to t with probability exp(-integral of m).
    gompertz = rules.Model("gompertz")
    gompertz.add_parameter({"g": 0.5, "m0": 0.1}).add_species({"A": 20})
    gompertz.add_ode("m", "m0", "g * m", independent=True)
    gompertz.add_reaction("m * A", ["A"], [], name="death")
    result = rules.ensemble(
        gompertz,
        end=4.0,
        report_times=[0.0, 4.0],
        replicates=400,
        seed=13,
        threads=2,
        method="hybrid",
    ).simulate()
    survival = np.exp(-0.1 * (np.exp(2.0) - 1.0) / 0.5)
    mean = result.statistics["mean"]
    error = result.statistics["mean_standard_error"]["A"].iloc[1]
    require(
        gompertz.to_dict()["schema"] == "crn-craft.rules-model/0.4"
        and result.metadata["continuous"] == ["m"]
        and np.isclose(mean["m"].iloc[1], 0.1 * np.exp(2.0), rtol=1e-5)
        and abs(mean["A"].iloc[1] - 20 * survival) <= 5 * error,
        "installed hybrid solves an independent ODE and matches the binomial law",
    )
    dependent = gompertz.to_dict()
    dependent["continuous"][0]["derivative"] = "g * m * A"
    try:
        rules.ensemble(dependent, end=1.0, report_times=[1.0], method="hybrid")
    except c.CrnError as error:
        first = error.diagnostics[0]
        require(
            (first["code"], first["path"]) == ("CRN-R012", "/continuous/0/derivative"),
            "an independent ODE that reads a species is refused",
        )
    else:
        raise RuntimeError("an independent ODE that reads a species was accepted")


def check_compiled_loop_refusals():
    """The compiled-loop API without cargo: source, driver and refusals."""
    model = three_class_model()
    source = rules.loop_source(model)
    require(
        "impl crn_craft::rules::RulesLoop for GeneratedRulesLoop" in source
        and "crn_craft_rules" not in source,
        "installed loop_source returns the generated loop",
    )
    rust_version, files = c._native.rules_loop_driver()
    require(
        rust_version and {"src/main.rs", "src/rules_payload.rs"} <= dict(files).keys(),
        "the installed extension carries the compiled-loop driver",
    )
    # Since 0.7.0 the driver, like the generated loop, builds against the
    # crn-craft crate alone.
    require(
        "use crn_craft::rules::" in dict(files)["src/main.rs"]
        and not any("crn_craft_rules" in text for _, text in files),
        "the compiled-loop driver names the crn-craft crate",
    )
    try:
        rules.compile_loop(model, cargo="crn-craft-no-such-cargo", progress=False)
    except rules.CompiledLoopError as error:
        require(error.reason == "cargo_missing", "missing cargo reason")
    else:
        raise RuntimeError("compile_loop without cargo was accepted")
    with tempfile.TemporaryDirectory() as folder:
        stranger = rules.CompiledLoop(
            fingerprint="0" * 64,
            model_sha256="0" * 64,
            version=c._native.__version__,
            directory=folder,
            cached=True,
        )
        try:
            rules.ensemble(
                model,
                end=1.0,
                report_times=[1.0],
                method="static_direct",
                compiled=stranger,
            )
        except c.CrnError as error:
            first = error.diagnostics[0]
            require(
                (first["code"], first["path"]) == ("CRN-R020", "/model"),
                "a loop for other model bytes is refused at /model",
            )
        else:
            raise RuntimeError("a loop for other model bytes was accepted")


# Expected output of Result.save, next to this script (not in the package).
SIR_REPORT = Path(__file__).resolve().parent / "fixtures" / "sir-report.md"


def check_report():
    """Result.save: the SIR example's markdown report, byte for byte.

    The expected file is validation/fixtures/sir-report.md. Its line
    `crn-craft = <version>` names the version it was written with; that one
    line is compared with the installed version instead, so the file stays
    valid across version bumps. Every other byte must be equal. The run asks
    for the quantiles 0.025, 0.5 and 0.975, so the documents are
    crn-craft.run/1.0 and crn-craft.model/1.0.
    """
    sir = c.Model("sir")
    sir.add_species({"S": 999, "I": 1, "R": 0})
    sir.add_parameter({"beta": 2e-4, "gamma": 0.01})
    sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
    sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    sir.set_stop("I == 0")
    settings = {
        "end": 250.0,
        "report_times": [0, 50, 100, 150, 200, 250],
        "replicates": 10,
        "seed": 1,
        "keep_paths": "all",
        "method": "static_direct",
    }
    # The top-level API asks for these quantiles by default; through
    # crn_craft.rules they are asked for, and the documents are the same
    # except that crn_craft.rules writes out the linear selector of every
    # rules-run version, where the top level leaves "auto" to the engine.
    prepared = c.ensemble(sir, **settings)
    old = rules.Model("sir")
    old.add_species({"S": 999, "I": 1, "R": 0})
    old.add_parameter({"beta": 2e-4, "gamma": 0.01})
    old.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
    old.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    old.set_stop("I == 0")
    old_model, old_request = rules.ensemble(
        old, quantiles=(0.025, 0.5, 0.975), **settings
    )._documents
    old_request = json.loads(old_request)
    require(
        old_model == prepared._documents[0]
        and old_request["method"].pop("selector") == "linear"
        and old_request == json.loads(prepared._documents[1]),
        "the top-level SIR documents equal rules.ensemble(quantiles=...) "
        "but for the linear selector",
    )
    result = prepared.simulate()
    require(
        list(result.statistics["quantiles"]) == [0.025, 0.5, 0.975]
        and result.statistics["median"] is result.statistics["quantiles"][0.5],
        "quantile statistics and the median",
    )
    expected = SIR_REPORT.read_bytes().split(b"\n")
    version = [i for i, line in enumerate(expected) if line.startswith(b"crn-craft = ")]
    require(len(version) == 1, "the SIR fixture has one crn-craft version line")
    expected[version[0]] = b"crn-craft = " + c.__version__.encode()
    expected = b"\n".join(expected)
    with tempfile.TemporaryDirectory() as tmp:
        first = result.save(Path(tmp) / "sir.md")
        second = result.save(Path(tmp) / "again.md")
        written = first.read_bytes()
        require(
            isinstance(first, Path)
            and first == Path(tmp) / "sir.md"
            and written == second.read_bytes(),
            "saving a result twice gives identical bytes",
        )
        require(
            written == expected,
            f"the SIR report differs from {SIR_REPORT}:\n{written.decode()}",
        )
        report = c.load_report(result.save(Path(tmp) / "sir.json"))
        require(
            report["schema"] == "crn-craft.report/1.0"
            and report["crn_craft_version"] == c.__version__
            and report["model"]
            == sir.to_dict()
            == {**old.to_dict(), "schema": "crn-craft.model/1.0"}
            == result.model
            and report["request"]["schema"] == "crn-craft.run/1.0"
            and report["request"] == result.request
            and report["result"]["observables"] == ["I", "R", "S"]
            and report["result"]["terminations"] == {"completed": 10}
            and rules.load_report is c.load_report,
            "the JSON report holds the model and request that ran",
        )


def poisson(n, mean):
    """The Poisson probability of n, in floating point."""
    return math.exp(n * math.log(mean) - mean - math.lgamma(n + 1)) if mean else n == 0


def check_unified():
    """The unified API from the installed wheel (M6a): the top-level names,
    Model.fsp against the Poisson law of immigration-death, Model.analyze of
    the enzyme model and a third-order reaction (FSP against its exact law,
    and an ensemble against the FSP mean within 5 standard errors)."""
    require(
        set(c.__all__)
        == {
            "CompiledLoopError",
            "CrnError",
            "Model",
            "archive",
            "check",
            "compile_loop",
            "convert",
            "describe",
            "ensemble",
            "load",
            "load_example",
            "load_report",
            "loop_source",
            "prepare",
            "read_bngl",
            "read_sbml",
            "verify_archive",
        }
        and not hasattr(c, "run_spec")
        and not hasattr(c, "simulate_ssa")
        and c.describe()["schema"] == "crn-craft.describe/1.0",
        "installed top level",
    )
    death = c.load_example("immigration-death")
    times = [0.0, 2.0, 5.0, 10.0]
    fsp = death.fsp(end=10.0, report_times=times, bounds={"A": 150})
    states = fsp.states["A"].to_numpy()
    for i, t in enumerate(times):
        # Immigration-death from zero: Poisson with mean 20 (1 - exp(-0.2 t)).
        mean = 20.0 * (1.0 - math.exp(-0.2 * t))
        exact = np.array([poisson(int(n), mean) for n in states])
        probabilities = fsp.probabilities.isel(time=i).to_numpy()
        require(
            float(np.max(np.abs(probabilities - exact))) <= 1e-10
            and abs(fsp.mean["A"].iloc[i] - mean) < 1e-9
            and fsp.truncation["physical_leakage_estimate"].iloc[i] < 1e-12,
            f"installed Model.fsp matches the Poisson law at t = {t}",
        )
    report = c.load_example("enzyme-mechanism").analyze()
    conservation = report["analysis"]["conservation"]
    require(
        report["schema"] == "crn-craft.analysis/1.0"
        and report["success"]
        and report["analysis"]["species"] == ["S", "E", "ES", "P"]
        and conservation["status"] == "exact"
        and [b["coefficients"] for b in conservation["basis"]]
        == [["0", "1", "1", "0"], ["1", "0", "1", "1"]]
        and [b["initial_total"] for b in conservation["basis"]] == ["10", "50"],
        "installed Model.analyze: the enzyme and substrate totals are conserved",
    )
    # 3A -> 0 at k C(A, 3) from A = 10: a pure-death chain 10 -> 7 -> 4 -> 1
    # with rates 120 k, 35 k and 4 k; 1 molecule is absorbing.
    k = 0.01
    tri = c.Model("third_order").add_species("A", 10)
    tri.add_reaction(k, ["A", "A", "A"], [], name="triple")
    require(
        tri.to_dict()["rules"][0]["reactants"]
        == [{"label": "x", "type": "A", "count": 3}]
        and c.check(tri)["valid"],
        "a third-order reaction is one pattern with count 3",
    )
    a, b, d = 120 * k, 35 * k, 4 * k
    times = [0.0, 1.0, 5.0, 20.0]
    law = c.Model.fsp(tri, end=20.0, report_times=times, bounds={"A": 10})
    for i, t in enumerate(times):
        p10 = math.exp(-a * t)
        p7 = a / (b - a) * (math.exp(-a * t) - math.exp(-b * t))
        p4 = (
            a
            * b
            * (
                math.exp(-a * t) / ((b - a) * (d - a))
                + math.exp(-b * t) / ((a - b) * (d - b))
                + math.exp(-d * t) / ((a - d) * (b - d))
            )
        )
        exact = {10: p10, 7: p7, 4: p4, 1: 1.0 - p10 - p7 - p4}
        got = dict(
            zip(
                law.states["A"].tolist(),
                law.probabilities.isel(time=i).to_numpy().tolist(),
                strict=True,
            )
        )
        require(
            set(got) == {10, 7, 4, 1}
            and max(abs(got[n] - exact[n]) for n in exact) <= 1e-10,
            f"installed third-order FSP matches the exact law at t = {t}",
        )
    ensemble = c.ensemble(
        tri,
        end=20.0,
        report_times=times,
        replicates=2000,
        seed=11,
        threads=2,
        keep_paths=50,
        method="static_direct",
    ).simulate()
    values = ensemble.paths.sel(variable="A").to_numpy()
    mean = ensemble.statistics["mean"]["A"].to_numpy()
    error = ensemble.statistics["mean_standard_error"]["A"].to_numpy()
    require(
        set(np.unique(values).tolist()) <= {10.0, 7.0, 4.0, 1.0}
        and bool(np.all(np.diff(values, axis=1) <= 0))
        and mean[0] == 10.0
        and bool(
            np.all(np.abs(mean[1:] - law.mean["A"].to_numpy()[1:]) <= 5 * error[1:])
        )
        and ensemble.statistics["median"]["A"].iloc[0] == 10.0,
        "installed third-order ensemble matches its exact law",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    checks = []
    cli, importer = script("crn-craft"), script("crn-craft-import")
    version = importlib.metadata.version("crn-craft")
    require(c.__version__ == version, "package __version__ matches metadata")
    require(
        python_spelling(c._native.__version__) == version,
        "native extension and distribution versions agree",
    )
    reported = subprocess.run([cli, "--version"], capture_output=True, check=False)
    require(
        reported.returncode == 0
        and json.loads(reported.stdout)
        == {"version": version, "interface": "python", "milestone": "unified-engine"},
        "installed entry point is the Python CLI",
    )
    checks.append("version_and_entry_points")
    m = c.load_example("immigration-death")
    prepared = c.ensemble(
        m,
        end=2.0,
        report_times=[0.0, 1.0, 2.0],
        replicates=32,
        seed=77,
        keep_paths=32,
        method={"kind": "static_direct", "waiting_time": "inverse"},
        selector="linear",
    )
    a = prepared.simulate()
    b = prepared.simulate()
    require(
        a.metadata == b.metadata and a.paths.equals(b.paths),
        "installed replay identity",
    )
    checks.append("native_API_exact_identity")
    require(
        c.describe()["schema"] == "crn-craft.describe/1.0",
        "installed capabilities",
    )
    checks.append("discovery")
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        model_file = folder / a.request["model"]["path"]
        model_file.write_bytes(rules.model_bytes(a.model))
        (folder / "run.json").write_bytes(rules.model_bytes(a.request))
        proc = subprocess.run(
            [
                cli,
                "run",
                str(model_file),
                str(folder / "run.json"),
                "--archive",
                str(folder / "out"),
            ],
            capture_output=True,
            check=False,
        )
        require(proc.returncode == 0, proc.stdout.decode() + proc.stderr.decode())
        archived = json.loads((folder / "out/result.json").read_bytes())
        require(
            {k: v for k, v in archived.items() if k != "replicates"}
            == {k: v for k, v in a.metadata.items() if k != "replicates"}
            and [r["events"] for r in archived["replicates"]]
            == a.records["events"].tolist(),
            "installed CLI/API identity",
        )
        require(
            subprocess.run(
                [cli, "verify", str(folder / "out")],
                capture_output=True,
                check=False,
            ).returncode
            == 0,
            "installed verify",
        )
        for command in [
            ["describe", "run-1.0"],
            ["check", str(model_file), str(folder / "run.json")],
            ["analyze", str(model_file)],
            ["export", str(model_file), "gillespy2"],
        ]:
            extra = subprocess.run([cli, *command], capture_output=True, check=False)
            require(extra.returncode == 0, "installed discovery/check/export parity")
        tau = c.ensemble(
            m,
            end=2.0,
            report_times=[0.0, 2.0],
            replicates=32,
            seed=77,
            method={"kind": "tau", "epsilon": 0.05},
        )
        approx = tau.simulate()
        (folder / "tau.json").write_bytes(rules.model_bytes(approx.request))
        extra = subprocess.run(
            [cli, "run", str(model_file), str(folder / "tau.json")],
            capture_output=True,
            check=False,
        )
        require(
            extra.returncode == 0
            and json.loads(extra.stdout)["method"]["kind"] == "tau"
            and json.loads(extra.stdout)["sampler"].startswith(
                "chacha8-u53/rules-approximate"
            ),
            "installed tau CLI parity",
        )
        (folder / "decay.xml").write_bytes(SBML)
        imported = subprocess.run(
            [
                importer,
                "sbml",
                str(folder / "decay.xml"),
                str(folder / "decay.json"),
                "--convention",
                "stochastic_count_propensity",
            ],
            capture_output=True,
            check=False,
        )
        require(imported.returncode == 0, imported.stdout.decode())
        decay = json.loads((folder / "decay.json").read_bytes())["model"]
        # The importers return the converted crn-craft.model/1.0 (M5f).
        require(
            c.Model.from_dict(decay).analyze()["success"]
            and decay["schema"] == "crn-craft.model/1.0"
            and float(decay["rules"][0]["rate"]) == 0.5,
            "installed SBML import script",
        )
        checks.append("installed_CLI_archive_and_import_scripts")
    result = rules.fsp(m.to_dict(), 2.0, [0.0, 2.0], {"A": 100}, 101, 1e-12)
    require(not result["numerical_error_certified"], "installed FSP semantics")
    # Immigration-death from zero is Poisson with mean 20 (1 - exp(-0.2 t)).
    final = result["snapshots"][-1]
    mean = 20.0 * (1.0 - math.exp(-0.4))
    require(
        abs(final["projected_first_moment"][0] - mean) < 1e-9
        and abs(final["projected_second_raw_moment"][0] - mean - mean**2) < 1e-8,
        "installed FSP matches the Poisson moments",
    )
    checks.append("FSP")
    require(
        approx.metadata["termination"] == "completed"
        and approx.metadata["method"]["kind"] == "tau"
        and approx.metadata["sampler"] == "chacha8-u53/rules-approximate-cgp-poisson/2",
        "installed approximate method identity",
    )
    checks.append("tau")
    net = b"begin species\n1 A() 2\nend species\nbegin reactions\n1 1 0 1.0\nend reactions\n"
    require(
        c.Model.from_dict(
            import_bng_net(net, rate_convention="binomial")["model"]
        ).analyze()["success"],
        "installed importer",
    )
    from crn_craft.importers import ImportFailure

    try:
        import_bng_net(net.replace(b"1.0", b"1e-400"), rate_convention="binomial")
    except ImportFailure:
        pass
    else:
        raise RuntimeError(
            "installed importer silently underflowed a positive coefficient"
        )
    checks.append("importer")
    require(
        sampled_first_passage([0, 1], [[0, 0]], 1)["not_detected"] == 1,
        "installed observables",
    )
    t = np.linspace(0, 100, 501)
    require(
        oscillation(t, np.sin(2 * np.pi * t / 10))["classification"] == "oscillatory"
        and oscillation(t, np.full(len(t), 0.3))["classification"] == "nonoscillatory",
        "installed oscillation detector",
    )
    try:
        c.prepare(m, dict(a.request, schema="crn-craft.run/9"))
    except c.CrnError as error:
        require(error.diagnostics[0]["code"] == "CRN-E003", "installed CrnError")
    else:
        raise RuntimeError("installed prepare accepted an invalid run")
    checks.append("observables")
    # The rules engine (`crn_craft.rules`) ships in the wheel; the mtDNA manual
    # depends on it. Immigration-death has the Poisson law with mean
    # (k_in/k_out)(1 - exp(-k_out t)), so the means at t = 5 and 10 must lie
    # within 5 standard errors, and threads must not change any result.
    death = {
        "schema": "crn-craft.rules-model/0.1",
        "name": "immigration_death",
        "parameters": {"k_in": 10.0, "k_out": 1.0},
        "types": [{"name": "X", "attributes": []}],
        "initial": [],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {
                "name": "immigration",
                "reactants": [],
                "rate": "k_in",
                "products": [{"type": "X"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "X"}],
                "rate": "k_out",
                "products": [],
            },
        ],
    }

    def rules_run(threads):
        return rules.ensemble(
            death,
            end=10.0,
            report_times=[0.0, 5.0, 10.0],
            replicates=400,
            seed=7,
            threads=threads,
            method="static_direct",
            trackers=[{"observable": "X"}],
        ).simulate()

    serial, parallel = rules_run(1), rules_run(2)
    mean = serial.statistics["mean"]["X"].to_numpy()
    error = serial.statistics["mean_standard_error"]["X"].to_numpy()
    exact = 10.0 * (1.0 - np.exp(-np.array([0.0, 5.0, 10.0])))
    require(
        mean[0] == 0.0 and bool(np.all(np.abs(mean[1:] - exact[1:]) <= 5 * error[1:])),
        "installed rules ensemble matches the immigration-death law",
    )
    require(
        serial.records.equals(parallel.records)
        and serial.rule_events.equals(parallel.rule_events)
        and serial.statistics["mean"].equals(parallel.statistics["mean"])
        and set(serial.tracker_columns(0).columns)
        >= {"appearances", "first_passage", "integral"},
        "installed rules results are independent of threads and expose trackers",
    )
    # Results are pandas objects: records and rule counts by replicate,
    # statistics by report time, and every path with keep_paths="all".
    records = serial.records
    require(
        isinstance(records, pd.DataFrame)
        and records.index.name == "replicate"
        and list(records.index) == list(range(400))
        and isinstance(records["termination"].dtype, pd.CategoricalDtype)
        and bool((records["termination"] == "completed").all())
        and records["events"].dtype == np.uint64
        and serial.rule_names == ("immigration", "death")
        and list(serial.rule_events.columns) == ["immigration", "death"]
        and int(serial.rule_events.to_numpy().sum()) == int(records["events"].sum())
        and list(serial.statistics["mean"].index) == [0.0, 5.0, 10.0]
        and serial.statistics["samples"] == 400,
        "installed rules results are DataFrames",
    )
    kept = rules.ensemble(
        death,
        end=10.0,
        report_times=[0.0, 5.0, 10.0],
        replicates=50,
        seed=7,
        keep_paths="all",
    ).simulate()
    paths = kept.paths  # replicates x report times x variables
    by_time = paths.sel(variable="X").mean("replicate").values
    frame = kept.paths_frame()
    require(
        isinstance(paths, xr.DataArray)
        and paths.dims == ("replicate", "time", "variable")
        and paths.shape == (50, 3, 1)
        and paths["time"].values.tolist() == [0.0, 5.0, 10.0]
        and paths["variable"].values.tolist() == ["X"]
        and paths["replicate"].values.tolist() == list(range(50))
        and paths.attrs["sampler"] == kept.metadata["sampler"]
        and paths.attrs["root_seed"] == 7
        and np.allclose(by_time, kept.statistics["mean"]["X"], rtol=1e-12, atol=0),
        "installed rules keep every path as a DataArray",
    )
    require(
        frame.shape == (150, 1)
        and list(frame.index.names) == ["replicate", "time"]
        and np.array_equal(
            frame["X"].unstack("time").to_numpy(), paths.sel(variable="X").values
        ),
        "paths_frame() is the 0.4.0 DataFrame",
    )
    data = kept.to_xarray()
    require(
        data["X"].dims == ("replicate", "time") and data["X"].shape == (50, 3),
        "installed to_xarray dimensions",
    )
    require(
        kept.paths.identical(
            rules.ensemble(
                death,
                end=10.0,
                report_times=[0.0, 5.0, 10.0],
                replicates=50,
                seed=7,
                keep_paths=True,
            )
            .simulate()
            .paths
        ),
        "keep_paths=True keeps every path, as 'all' does",
    )
    check_pandas_results(serial)
    require(
        {"mtdna_direct", "hybrid"} <= set(rules.describe()["methods"]),
        "installed rules describe lists mtdna_direct and hybrid",
    )
    checks.append("rules_engine")
    check_event_log_default_budget()
    checks.append("rules_event_log_default_budget")
    check_builder()
    checks.append("rules_builder")
    check_builder_mappings()
    checks.append("rules_builder_mappings")
    check_builder_odes()
    checks.append("rules_builder_odes")
    check_builder_stop()
    checks.append("rules_builder_stop")
    check_hybrid()
    checks.append("rules_hybrid")
    check_compiled_loop_refusals()
    checks.append("rules_compiled_loop_refusals")
    check_report()
    checks.append("rules_result_report")
    check_unified()
    checks.append("unified_api")
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.installed-validation/0.1",
                "passed": True,
                "checks": checks,
                "package_path": str(Path(c.__file__).resolve()),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "checks": len(checks)}))


if __name__ == "__main__":
    main()
