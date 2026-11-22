"""Fixed-work law checks for protocols, delay queues, projections and experiments."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import zipfile

import crn_craft as crn
from crn_craft.likelihood import (
    HistogramGroup,
    binned_probability,
    count_probability,
    histogram_score,
    mean_score,
    measured_counts,
)
from crn_craft.importers import ImportFailure
import numpy as np
from scipy.linalg import expm
from scipy.stats import binom, norm, poisson


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cdf_gap(sample, cdf):
    counts = np.bincount(np.asarray(sample, dtype=int))
    return float(
        np.max(np.abs(np.cumsum(counts) / counts.sum() - cdf(np.arange(len(counts)))))
    )


def np_gap(p, cdf):
    return float(np.max(np.abs(np.cumsum(p) - cdf(np.arange(len(p))))))


def immigration():
    m = crn.Model("timed_immigration")
    m.add_species("X", 0)
    m.add_parameter("birth", 2)
    m.add_reaction(crn.mass_action("birth"), [], ["X"], name="birth_event")
    return m


def timed_laws():
    m = immigration()
    policy = crn.Protocol(time_unit="hour").at(
        3, parameters={"birth": 4}, partition={"X": 0.5}
    )
    n = 20000
    bound = float(np.sqrt(np.log(400) / (2 * n)))
    rows = {}
    for method in ("static_direct", "next_reaction"):
        args = {
            "end": 5,
            "replicates": n,
            "keep_paths": n,
            "protocol": policy,
            "seed": 123,
            "method": method,
            "strict": True,
            "quantiles": (),
        }
        tested = crn.ensemble(m, report_times=[3, 5], threads=3, **args).simulate()
        grid = crn.ensemble(m, report_times=[0, 1, 3, 5], threads=1, **args).simulate()
        a = tested.paths.sel(variable="X").values
        b = grid.paths.sel(variable="X").values[:, -2:]
        gaps = [cdf_gap(a[:, i], poisson(mu).cdf) for i, mu in enumerate((3, 11))]
        require(max(gaps) < bound, "scheduled SSA law")
        require(np.array_equal(a, b), "report grid or thread count changes SSA")
        rows[method] = {"cdf": gaps, "bound": bound}
    g = crn.prepare_fsp(m, {"X": 48}, protocol=policy)
    law = g.solve(5, [0, 3, 5], tolerance=1e-12)
    gaps = [
        np_gap(law.marginal("X", at=i), poisson(mu).cdf) for i, mu in ((1, 3), (2, 11))
    ]
    require(max(gaps) < 1e-9, "partitioned FSP")
    resumed = g.solve(5, [5], initial=law.snapshots[1], tolerance=1e-12)
    require(
        np.max(abs(resumed.marginal("X") - law.marginal("X"))) < 1e-11,
        "FSP continuation",
    )
    forced = crn.prepare_fsp(m, {"X": 48}, initial={(0,): 0.25, (1,): 0.75})
    result = forced.solve(
        2,
        [0, 1, 2],
        signals={"birth_event": crn.Signal(lambda t: 2 + t)},
        tolerance=1e-10,
    )

    def expected(x):
        return 0.25 * poisson(6).cdf(x) + 0.75 * poisson(6).cdf(x - 1)

    error = np_gap(result.marginal("X"), expected)
    require(error < 1e-8, "smooth law with mixed initial state")
    for table in (
        crn.Signal.table([0, 2], [2, 4]),
        crn.Signal.table([0, 1], [2, 4], interpolation="step"),
    ):
        tabled = forced.solve(2, [2], signals={"birth_event": table}, tolerance=1e-10)
        require(np_gap(tabled.marginal("X"), expected) < 1e-8, "signal table law")
    received = []
    reduced = forced.solve(
        2, [0, 1, 2], on_report=received.append, keep_snapshots=False
    )
    require(len(received) == 3 and not reduced.snapshots, "report reductions")
    midpoint = forced.solve(
        2,
        [2],
        signals={"birth_event": crn.Signal(lambda t: 2 + t)},
        solver="midpoint",
        tolerance=1e-9,
    )
    require(np_gap(midpoint.marginal("X"), expected) < 1e-8, "exponential midpoint law")
    # Birth 2 + t and death 1 do not commute: Poisson mean (1 - e^-t) + t.
    decaying = immigration()
    decaying.add_parameter("death", 1.0)
    decaying.add_reaction(crn.mass_action("death"), ["X"], [], name="death_event")
    dg = crn.prepare_fsp(decaying, {"X": 48})
    noncommuting = {}
    for solver, tolerance, bound in (("bdf", 1e-10, 1e-8), ("midpoint", 1e-6, 1e-6)):
        law = dg.solve(
            2,
            [2],
            signals={"birth_event": crn.Signal(lambda t: 2 + t)},
            solver=solver,
            tolerance=tolerance,
        )
        gap = np_gap(law.marginal("X"), poisson(3 - np.exp(-2)).cdf)
        require(gap < bound and law.solver == solver, f"forced decay law {solver}")
        noncommuting[solver] = gap
    try:
        dg.solve(
            2,
            [2],
            signals={"birth_event": crn.Signal(lambda t: 2 + t)},
            solver="midpoint",
            tolerance=1e-14,
        )
    except RuntimeError as e:
        require("rounding level" in str(e), str(e))
    else:
        raise RuntimeError("unreachable midpoint tolerance accepted")
    small = crn.prepare_fsp(m, {"X": 3}, protocol=policy)
    overflow = small.solve(5, [3, 5], tolerance=1e-12)
    resumed = small.solve(5, [5], initial=overflow.snapshots[0], tolerance=1e-12)
    require(
        0
        < overflow.snapshots[0].physical_leakage
        < resumed.snapshots[0].physical_leakage,
        "escaped probability must persist",
    )

    # Dense reference: states 0..3 and the sink, with the partition kernel.
    def generator(rate):
        q = np.zeros((5, 5))
        for i in range(4):
            q[i, i] -= rate
            q[i, i + 1] += rate
        return q

    kernel = np.zeros((5, 5))
    kernel[4, 4] = 1
    for n in range(4):
        kernel[n, : n + 1] = binom.pmf(np.arange(n + 1), n, 0.5)
    hour3 = np.eye(5)[0] @ expm(3 * generator(2)) @ kernel
    for snap, reference in zip(
        overflow.snapshots, (hour3, hour3 @ expm(2 * generator(4))), strict=True
    ):
        got = np.r_[snap.probabilities, snap.physical_leakage]
        require(np.max(abs(got - reference)) < 1e-12, "sink across a partition")
    pulse = crn.Protocol().at(0, parameters={"birth": 0}, set_counts={"X": 4})
    pulse.at(1, parameters={"birth": 1}, partition={"X": 1}).at(2, partition={"X": 0})
    run = crn.ensemble(
        m,
        end=2,
        report_times=[0, 1, 2],
        protocol=pulse,
        keep_paths=1,
        strict=True,
        quantiles=(),
        event_log={"max_entries": 20},
    )
    got = run.simulate()
    require(
        got.paths.sel(variable="X").values.tolist() == [[4, 4, 0]],
        "right-continuous time zero or horizon",
    )
    with TemporaryDirectory() as temp:
        archive = Path(temp) / "run"
        require(
            run.archive(archive)["valid"] and crn.verify_archive(archive)["valid"],
            "protocol archive",
        )
    switch = crn.Model("ordered_switch")
    switch.add_species({"A": 1, "B": 0}).add_parameter({"f": 2.0, "r": 0.0})
    switch.add_reaction(crn.mass_action("f"), ["A"], ["B"], name="forward")
    switch.add_reaction(crn.mass_action("r"), ["B"], ["A"], name="reverse")
    order = crn.Protocol().at(0.5, parameters={"f": 0, "r": 0.2})
    ordered = crn.prepare_fsp(switch, {"A": 1, "B": 1}, protocol=order).solve(
        2, [2], tolerance=1e-12
    )
    require(
        abs(ordered.marginal("B")[1] - (1 - np.exp(-1)) * np.exp(-0.3)) < 1e-11,
        "ordered generators were averaged",
    )
    # A channel with zero rate at preparation switched on by a binding.
    dormant = crn.Model("dormant_switch")
    dormant.add_species({"A": 0, "B": 1}).add_parameter({"f": 0.0, "r": 0.0})
    dormant.add_reaction(crn.mass_action("f"), ["A"], ["B"], name="forward")
    dormant.add_reaction(crn.mass_action("r"), ["B"], ["A"], name="reverse")
    warm = crn.prepare_fsp(dormant, {"A": 1, "B": 1})
    switched = warm.solve(1, [1], parameters={"r": 1.0}, tolerance=1e-12)
    require(
        abs(switched.marginal("A")[1] - (1 - np.exp(-1))) < 1e-11,
        "channel switched on by a parameter binding",
    )
    cold = crn.prepare_fsp(dormant, {"A": 1, "B": 1}).solve(
        1, [1], parameters={"r": 1.0}, tolerance=1e-12
    )
    require(
        np.array_equal(
            cold.snapshots[0].probabilities, switched.snapshots[0].probabilities
        ),
        "cold and reused projections differ",
    )
    return {"ssa": rows, "fsp": gaps, "smooth": error, "noncommuting": noncommuting}


def delayed_laws():
    m = crn.Model("busy_copy")
    m.add_species({"Busy": 0, "Product": 0})
    m.add_parameter("start_rate", 2.5)
    m.add_reaction(crn.mass_action("start_rate"), [], ["Busy"], name="start")
    m.add_reaction(0.7, ["Busy"], [], name="loss")
    p = crn.Protocol().delay("start", 1.0, owner="Busy", products={"Product": 1})
    n = 12000
    v = (
        crn.ensemble(
            m,
            end=3,
            report_times=[1, 3],
            replicates=n,
            keep_paths=n,
            protocol=p,
            seed=1616,
            strict=True,
            quantiles=(),
            threads=3,
        )
        .simulate()
        .paths
    )
    a = v.sel(variable="Product").values
    require(a[:, 0].max() == 0, "fixed completion fired too early")
    gap = cdf_gap(a[:, 1], poisson(5 * np.exp(-0.7)).cdf)
    require(gap < 0.016, "completion/cancellation distribution")
    cancel = crn.Protocol().delay("start", 3, owner="Busy", products={"Product": 1})
    cancel.at(2, set_counts={"Busy": 0}, parameters={"start_rate": 0})
    lost = crn.ensemble(
        m,
        end=6,
        report_times=[2, 6],
        replicates=64,
        keep_paths=64,
        protocol=cancel,
        seed=131,
        strict=True,
    ).simulate()
    require(
        float(lost.paths.sel(variable="Product").values.max()) == 0,
        "scheduled removal did not cancel its pending products",
    )
    single = crn.Model("completion_boundary")
    single.add_species({"Free": 1, "Busy": 0, "Product": 0})
    single.add_reaction(1, ["Free"], ["Busy"], name="start")
    policy = crn.Protocol().delay("start", 1, owner="Busy", products={"Product": 1})
    tick = (
        crn.ensemble(
            single,
            end=5,
            report_times=[5],
            protocol=policy,
            seed=81,
            event_log={"max_entries": 16},
            strict=True,
        )
        .simulate()
        .observations[0]["log"][0]["time"]
    )
    for lead in (False, True):
        boundary = tick + 1
        stamp = np.nextafter(boundary, 0) if lead else boundary
        p = crn.Protocol().delay("start", 1, owner="Busy", products={"Product": 1})
        p.at(stamp, set_counts={"Busy": 0})
        observed = crn.ensemble(
            single,
            end=5,
            report_times=[boundary, 5],
            protocol=p,
            seed=81,
            keep_paths=1,
            event_log={"max_entries": 16},
            strict=True,
        ).simulate()
        log = [e["protocol"] for e in observed.observations[0]["log"]]
        values = observed.paths.sel(variable="Product").values
        require(
            np.array_equal(values, [[int(not lead), int(not lead)]]),
            "completion priority",
        )
        require(
            log[0]["kind"] == "initiation" and log[0]["pending"] == 1,
            "queued completion log",
        )
        require(
            (log[-1]["kind"], log[-1]["cancelled"]) == ("intervention", int(lead)),
            "intervention cancellation accounting",
        )
    return gap


def spatial_laws():
    n = 6000
    graph = crn.SpatialModel("two_sites", {"a": 1, "b": 1})
    graph.add_species("X", {"a": 10}).diffusion("X", "a", "b", 0.35, 1)
    run = crn.ensemble(
        graph.model,
        end=2,
        report_times=[2],
        replicates=n,
        keep_paths=n,
        seed=1199,
        method="next_reaction",
        quantiles=(),
        strict=True,
    ).simulate()
    a, b = (
        run.paths.sel(variable=graph.label("X", s)).values[:, 0] for s in ("a", "b")
    )
    require(np.all(a + b == 10), "transport must conserve the counts")
    theoretical = binom(10, (1 - np.exp(-1.4)) / 2)
    gap = cdf_gap(b, theoretical.cdf)
    require(gap < 0.022, "two-site transition")
    unequal = crn.SpatialModel("unequal_sites", {"a": 1, "b": 4})
    unequal.add_species("X", {"a": 10}).transport("X", "a", "b", 0.4).transport(
        "X", "b", "a", 0.1
    )
    balance = crn.ensemble(
        unequal.model,
        end=24,
        report_times=[24],
        replicates=2000,
        seed=122,
        method="static_direct",
        quantiles=(),
        strict=True,
    ).simulate()
    equilibrium = float(balance.statistics["mean"]["X_b"].iloc[0])
    require(abs(equilibrium - 8) < 0.12, "fixed-volume detailed balance")
    pop = []
    for seed in range(4000):
        result = crn.spatial_turnover(
            [(0, 10), (0, 0)],
            end=2,
            target=10,
            feedback=0,
            mutant_weight=1,
            death=0,
            hop=0.35,
            seed=seed,
            labels=True,
        )
        require(sum(result["mutant"]) == 10, "local sampler loses particles")
        require(
            [len(x) for x in result["labels"]] == result["mutant"], "clone bookkeeping"
        )
        pop.append(result["mutant"][1])
    local = cdf_gap(pop, theoretical.cdf)
    require(local < 0.027, "local scheduler transition law")
    chain = local_against_generic()
    return {
        "graph_gap": gap,
        "local_gap": local,
        "volume_equilibrium": equilibrium,
        "local_vs_generic_gaps": chain,
    }


def two_sample_gap(a, b):
    a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
    last = int(max(a.max(), b.max())) + 1
    fa = np.bincount(a, minlength=last).cumsum() / len(a)
    fb = np.bincount(b, minlength=last).cumsum() / len(b)
    return float(np.max(np.abs(fa - fb)))


def local_against_generic():
    """The local turnover sampler against the generic engine on a small chain.

    Four regions, a fixed left reservoir and feedback strong enough that the
    birth clamp at zero is reached (about 5% of endpoint regions). Five count
    CDFs, 20,000 paths per sampler, and a two-sample bound from DKW at
    alpha = 0.01 / 10 per sample: 2 * sqrt(ln(2000) / 40000) = 0.02757.
    """
    n, initial = 20_000, [(6, 4), (8, 2), (10, 0), (10, 0)]
    rates = {"target": 10, "feedback": 0.2, "mutant_weight": 0.5, "death": 0.5}
    model = crn.Model("turnover_chain")
    model.add_parameter(rates)
    for i, (w, m) in enumerate(initial):
        model.add_species({f"W{i}": w, f"M{i}": m})
    for i in range(4):
        for kind in "WM":
            x = f"{kind}{i}"
            if i:
                birth = (
                    "max(0, death + feedback * (target - W{0} - mutant_weight * M{0}))"
                )
                model.add_reaction(crn.mass_action(birth.format(i)), [x], [x, x])
                model.add_reaction(0.5, [x], [])
            for j in (i - 1, i + 1):
                if 0 <= j < 4:
                    source = [] if i == 0 else [x]
                    target = [] if j == 0 else [f"{kind}{j}"]
                    count = initial[0][kind == "M"]
                    model.add_reaction(0.3 * count if i == 0 else 0.3, source, target)
    run = crn.ensemble(
        model,
        end=5,
        report_times=[5],
        replicates=n,
        seed=777,
        method="next_reaction",
        keep_paths="all",
        quantiles=(),
        strict=True,
    ).simulate()
    names = [f"{k}{i}" for i in range(4) for k in "WM"]
    generic = run.paths.sel(variable=names).values[:, 0].reshape(n, 4, 2)
    local = np.empty_like(generic)
    for k in range(n):
        result = crn.spatial_turnover(
            initial, end=5, hop=0.3, fixed_left=True, seed=100_000 + k, **rates
        )
        local[k, :, 0], local[k, :, 1] = result["wild"], result["mutant"]
    require(np.all(local[:, 0] == initial[0]), "the reservoir holds its counts")

    def columns(x):
        w, m = x[..., 0], x[..., 1]
        return [w[:, 1], m[:, 1], m[:, 3], m[:, 1:].sum(1), (w + m)[:, 1:].sum(1)]

    gaps = [
        two_sample_gap(a, b)
        for a, b in zip(columns(generic), columns(local), strict=True)
    ]
    require(max(gaps) < 0.02757, "local turnover sampler against the generic engine")
    still = crn.spatial_turnover([(5, 5), (0, 0)], end=5, hop=0, seed=1, **rates)
    require(still["wild"][1] == still["mutant"][1] == 0, "zero transport")
    return gaps


def stationary_laws():
    m = immigration()
    m.add_reaction(1, ["X"], [], name="decay")
    result = crn.stationary_fsp(
        m,
        {"X": 15},
        return_state=(0,),
        expansions=4,
        tolerance=1e-8,
        flux_tolerance=1e-8,
    )
    error = np_gap(result.marginal("X"), poisson(2).cdf)
    require(result.converged and error < 1e-8, "stationary Poisson law")
    finite = crn.Model("several_classes")
    finite.add_species({"A": 1, "B": 0, "C": 0})
    finite.add_reaction(1, ["A"], ["B"])
    finite.add_reaction(1, ["A"], ["C"])
    try:
        crn.stationary_fsp(finite, {"A": 1, "B": 1, "C": 1})
    except ValueError as e:
        require("2 closed" in str(e), str(e))
    else:
        raise RuntimeError("ambiguous stationary law accepted")
    resolved = crn.stationary_fsp(
        finite, {"A": 1, "B": 1, "C": 1}, class_state=(0, 1, 0)
    )
    require(resolved.marginal("B")[1] == 1, "class selection")
    switch = crn.Model("two_state_switch")
    switch.add_species({"A": 1, "B": 0})
    switch.add_reaction(2.0, ["A"], ["B"])
    switch.add_reaction(3.0, ["B"], ["A"])
    exact = crn.stationary_fsp(switch, {"A": 1, "B": 1})
    require(
        abs(exact.marginal("B")[1] - 0.4) < 1e-12
        and exact.method == "finite-closed-class"
        and result.method == "stationary-fsp-return-state",
        "two-state stationary law",
    )
    # Plan 5.3: the same small box with two return states (descriptive).
    returns = {}
    for state in (0, 8):
        small = crn.stationary_fsp(m, {"X": 8}, return_state=(state,))
        returns[str(state)] = np_gap(small.marginal("X"), poisson(2).cdf)
    return {
        "poisson_cdf": error,
        "outward_flux": result.outward_flux,
        "residual": result.residual_l1,
        "return_state_cdf_at_bound_8": returns,
    }


def observations():
    score = histogram_score([HistogramGroup("ctr", "a", 1, [7, 2, 1], [0.6, 0.3, 0.1])])
    require(
        abs(score.nats_per_cell + np.dot([7, 2, 1], np.log([0.6, 0.3, 0.1])) / 10)
        < 1e-12,
        "histogram nats",
    )
    exact = count_probability(
        [[1, 2], [2, 1]], [0.2, 0.8], [[0, 0], [1, 2]], detection=[0, 1]
    )
    require(np.all(exact == 0), "binomial observation extremes")
    exact = count_probability(
        [[1, 2], [2, 1]], [0.2, 0.8], [[0, 2], [0, 1]], detection=[0, 1]
    )
    require(np.array_equal(exact, [0.2, 0.8]), "zero and perfect detection")
    thinned = count_probability(
        [[0], [1], [2]], [0.2, 0.5, 0.3], [[0], [1], [2]], detection=0.5
    )
    hand = [0.2 + 0.25 + 0.075, 0.25 + 0.15, 0.075]
    require(np.max(abs(thinned - hand)) < 1e-14, "binomial detection")
    zero = histogram_score([HistogramGroup("c", "a", 1, [7, 0, 1], [0.9, 0, 0.1])])
    require(
        abs(zero.nats_per_cell + (7 * np.log(0.9) + np.log(0.1)) / 8) < 1e-12
        and zero.impossible_bins == 0,
        "empty bins contribute zero",
    )
    blocked = histogram_score(
        [HistogramGroup("c", "a", 1, [7, 1, 1], [0.9, 0, 0.09], missing_mass=0.01)]
    )
    optimistic = -(7 * np.log(0.91) + np.log(0.01) + np.log(0.1)) / 9
    require(
        blocked.nats_per_cell == np.inf
        and blocked.impossible_bins == 1
        and abs(blocked.optimistic_nats_per_cell - optimistic) < 1e-12,
        "occupied impossible bin",
    )
    try:
        mean_score([1, 2], [0, 0], [[1, 1], [1, 1]])
    except ValueError:
        pass
    else:
        raise RuntimeError("singular covariance accepted")
    missing = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; sys.modules['scipy'] = None; import crn_craft as crn\n"
            "try:\n crn.prepare_fsp(crn.Model('m'), {})\n"
            "except ModuleNotFoundError as e:\n print(e)",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    require("crn-craft[projection]" in missing.stdout, "missing SciPy message")
    axes = measured_counts([[1, 2, 3], [5, 2, 3]], ("N", "C", "P"), (["N", "C"], "P"))
    require(np.array_equal(axes, [[3, 3], [7, 3]]), "joint and summed readout")
    noise = binned_probability([1], [1], [-np.inf, 0, 1, np.inf], gain=1, sigma=1)
    require(
        np.max(abs(noise - np.diff(norm.cdf([-np.inf, -1, 0, np.inf])))) < 1e-12,
        "Gaussian readout",
    )
    require(
        abs(mean_score([4], [2], [[2]]) - 1) < 1e-12,
        "mean score",
    )
    return score.nats_per_cell


def experiment_fixture(temp):
    source = Path(temp) / "fixture.omex"
    core = "http://www.sbml.org/sbml/level3/version2/core"
    mathml = "http://www.w3.org/1998/Math/MathML"
    sbml = f"""<sbml xmlns="{core}" level="3" version="2">
    <model id="pulse" timeUnits="second" substanceUnits="item">
      <listOfCompartments>
        <compartment id="left" size="2" constant="true" units="litre"/>
        <compartment id="right" size="4" constant="true" units="litre"/>
      </listOfCompartments>
      <listOfSpecies>
        <species id="A" compartment="left" initialConcentration="0"
          hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
        <species id="B" compartment="right" initialAmount="0"
          hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      </listOfSpecies>
      <listOfParameters><parameter id="k" constant="false" value="0"/></listOfParameters>
      <listOfReactions><reaction id="transfer" reversible="false">
        <listOfReactants><speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw><math xmlns="{mathml}"><apply><times/><ci>k</ci><ci>A</ci></apply></math></kineticLaw>
      </reaction></listOfReactions>
      <listOfEvents><event id="pulse_start" useValuesFromTriggerTime="true">
        <trigger initialValue="false" persistent="true"><math xmlns="{mathml}"><apply>
          <geq/><csymbol definitionURL="http://www.sbml.org/sbml/symbols/time">time</csymbol><cn>1</cn>
        </apply></math></trigger>
        <listOfEventAssignments>
          <eventAssignment variable="k"><math xmlns="{mathml}"><cn>1</cn></math></eventAssignment>
          <eventAssignment variable="A"><math xmlns="{mathml}"><cn>5</cn></math></eventAssignment>
        </listOfEventAssignments>
      </event></listOfEvents>
    </model></sbml>"""
    sed = f"""<sedML xmlns="http://sed-ml.org/sed-ml/level1/version3" level="1" version="3">
    <listOfModels><model id="m" source="model.xml" language="urn:sedml:language:sbml"/></listOfModels>
    <listOfSimulations><uniformTimeCourse id="s" initialTime="0" outputStartTime="0"
      outputEndTime="2" numberOfPoints="2"><algorithm kisaoID="KISAO:0000029"/></uniformTimeCourse></listOfSimulations>
    <listOfTasks><task id="t" modelReference="m" simulationReference="s"/></listOfTasks>
    <listOfDataGenerators><dataGenerator id="g"><listOfVariables>
      <variable id="x" taskReference="t" target="/sbml/model/listOfSpecies/species[@id='A']"/>
    </listOfVariables><math xmlns="{mathml}"><ci>x</ci></math></dataGenerator></listOfDataGenerators>
    <listOfOutputs><report id="out"><listOfDataSets><dataSet id="A_conc" dataReference="g"/>
    </listOfDataSets></report></listOfOutputs></sedML>"""
    manifest = """<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">
    <content location="study.sedml" format="http://identifiers.org/combine.specifications/sed-ml" master="true"/>
    <content location="model.xml" format="http://identifiers.org/combine.specifications/sbml"/>
    </omexManifest>"""
    with zipfile.ZipFile(source, "w") as arc:
        arc.writestr("manifest.xml", manifest)
        arc.writestr("study.sedml", sed)
        arc.writestr("model.xml", sbml)
    return source


_CORE = "http://www.sbml.org/sbml/level3/version2/core"
_MML = "http://www.w3.org/1998/Math/MathML"
_TIME = "http://www.sbml.org/sbml/symbols/time"
_MANIFEST = """<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">
<content location="study.sedml" format="http://identifiers.org/combine.specifications/sed-ml" master="true"/>
</omexManifest>"""


def _sbml(
    law="<apply><times/><ci>k</ci><ci>A</ci></apply>",
    parameters='<parameter id="k" value="1" constant="false"/>',
    rules="",
    events="",
    reactants='<speciesReference species="A" stoichiometry="1" constant="true"/>',
    initial="5",
):
    return f"""<sbml xmlns="{_CORE}" level="3" version="2">
<model id="m" timeUnits="second" substanceUnits="item">
<listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" initialAmount="{initial}"
 hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters>{parameters}</listOfParameters>{rules}
<listOfReactions><reaction id="r" reversible="false"><listOfReactants>{reactants}
</listOfReactants><kineticLaw><math xmlns="{_MML}">{law}</math></kineticLaw></reaction>
</listOfReactions>{events}</model></sbml>"""


def _sed(
    source="model.xml",
    changes="",
    algorithm='<algorithm kisaoID="KISAO:0000029"/>',
    start="0",
    end="1",
    task="t",
):
    return f"""<sedML xmlns="http://sed-ml.org/sed-ml/level1/version3" level="1" version="3">
<listOfModels><model id="m" source="{source}" language="urn:sedml:language:sbml">
{changes}</model></listOfModels>
<listOfSimulations><uniformTimeCourse id="s" initialTime="{start}" outputStartTime="0"
 outputEndTime="{end}" numberOfPoints="2">{algorithm}</uniformTimeCourse></listOfSimulations>
<listOfTasks><task id="{task}" modelReference="m" simulationReference="s"/></listOfTasks>
<listOfDataGenerators><dataGenerator id="g"><listOfVariables><variable id="x"
 taskReference="{task}" target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='A']"/>
</listOfVariables><math xmlns="{_MML}"><ci>x</ci></math></dataGenerator></listOfDataGenerators>
<listOfOutputs><report id="out"><listOfDataSets><dataSet id="A" dataReference="g"/>
</listOfDataSets></report></listOfOutputs></sedML>"""


def _event(name, time, assignments, priority=None, initial="true", extra=""):
    rank = (
        ""
        if priority is None
        else f'<priority><math xmlns="{_MML}"><cn>{priority}</cn></math></priority>'
    )
    body = "".join(
        f'<eventAssignment variable="{k}"><math xmlns="{_MML}">{v}</math></eventAssignment>'
        for k, v in assignments.items()
    )
    return f"""<event id="{name}" useValuesFromTriggerTime="{extra or "true"}">
<trigger initialValue="{initial}" persistent="true"><math xmlns="{_MML}"><apply><geq/>
<csymbol definitionURL="{_TIME}">time</csymbol><cn>{time}</cn></apply></math></trigger>
{rank}<listOfEventAssignments>{body}</listOfEventAssignments></event>"""


def _archive(path, members, *, compression=zipfile.ZIP_DEFLATED):
    members = {"manifest.xml": _MANIFEST, "study.sedml": _sed(), **members}
    with zipfile.ZipFile(path, "w", compression) as arc:
        for name, data in members.items():
            arc.writestr(name, data)
    return path


def _refused(read):
    try:
        read()
    except ImportFailure as error:
        return f"{error.code}: {error}"
    return "accepted"


def experiment_refusals(temp):
    """Hostile or unsupported inputs: each must end in a CRN-E090 refusal."""
    temp = Path(temp)
    model = {"model.xml": _sbml()}
    seed_parameter = (
        '<algorithm kisaoID="KISAO:0000029"><listOfAlgorithmParameters>'
        '<algorithmParameter kisaoID="KISAO:0000488" value="7"/>'
        "</listOfAlgorithmParameters></algorithm>"
    )
    change = (
        "<listOfChanges><changeAttribute newValue='2' target=\"{}\"/></listOfChanges>"
    )
    time_ge = f'<apply><geq/><csymbol definitionURL="{_TIME}">time</csymbol><cn>1</cn></apply>'
    k_rule = (
        '<parameter id="k" constant="false"/><parameter id="p" constant="false"/>',
        f'<listOfRules><assignmentRule variable="k"><math xmlns="{_MML}"><ci>p</ci>'
        f'</math></assignmentRule><assignmentRule variable="p"><math xmlns="{_MML}">'
        "<ci>k</ci></math></assignmentRule></listOfRules>",
    )
    utf16 = (
        '<?xml version="1.0" encoding="UTF-16"?><!DOCTYPE s [<!ENTITY e "x">]>' + _sed()
    ).encode("utf-16-le")
    cases = {
        "source_parent": {"study.sedml": _sed("../model.xml"), **model},
        "source_absolute": {"study.sedml": _sed("/model.xml"), **model},
        "source_url": {"study.sedml": _sed("https://example.org/m.xml"), **model},
        "source_urn": {"study.sedml": _sed("urn:miriam:biomodels.db:BIOMD1"), **model},
        "member_parent": {"../model.xml": _sbml(), **model},
        "member_missing": {},
        "member_count": {**model, **{f"pad{i}": "" for i in range(256)}},
        "doctype": {"study.sedml": '<!DOCTYPE s [<!ENTITY e "x">]>' + _sed(), **model},
        "utf16_doctype": {"study.sedml": utf16, **model},
        "encoding_latin1": {
            "study.sedml": '<?xml version="1.0" encoding="ISO-8859-1"?>' + _sed(),
            **model,
        },
        "nested_archive": {"model.xml": b"PK\x03\x04" + bytes(64)},
        "kisao_next_reaction": {
            "study.sedml": _sed(algorithm='<algorithm kisaoID="KISAO:0000027"/>'),
            **model,
        },
        "kisao_seed_parameter": {
            "study.sedml": _sed(algorithm=seed_parameter),
            **model,
        },
        "xpath_descendant": {
            "study.sedml": _sed(changes=change.format("//parameter[@id='k']/@value")),
            **model,
        },
        "xpath_species_amount": {
            "study.sedml": _sed(
                changes=change.format(
                    "/sbml/model/listOfSpecies/species[@id='A']/@initialAmount"
                )
            ),
            **model,
        },
        "change_derived_parameter": {
            "study.sedml": _sed(
                changes=change.format(
                    "/sbml/model/listOfParameters/parameter[@id='k']/@value"
                )
            ),
            "model.xml": _sbml(
                parameters='<parameter id="k" constant="false"/>'
                '<parameter id="p" value="1" constant="true"/>',
                rules=f'<listOfRules><assignmentRule variable="k"><math xmlns="{_MML}">'
                "<ci>p</ci></math></assignmentRule></listOfRules>",
            ),
        },
        "cyclic_rules": {"model.xml": _sbml(parameters=k_rule[0], rules=k_rule[1])},
        "rule_reads_species": {
            "model.xml": _sbml(
                parameters='<parameter id="k" constant="false"/>',
                rules=f'<listOfRules><assignmentRule variable="k"><math xmlns="{_MML}">'
                "<ci>A</ci></math></assignmentRule></listOfRules>",
            )
        },
        "pair_law_squared": {
            "model.xml": _sbml(
                law="<apply><times/><ci>k</ci><ci>A</ci><ci>A</ci></apply>",
                reactants='<speciesReference species="A" stoichiometry="2" constant="true"/>',
            )
        },
        "signed_law_strict": {
            "model.xml": _sbml(
                law="<apply><minus/><apply><times/><ci>k</ci><ci>A</ci></apply><cn>1</cn></apply>"
            )
        },
        "state_trigger": {
            "model.xml": _sbml(
                events=f'<listOfEvents><event id="e" useValuesFromTriggerTime="true">'
                f'<trigger initialValue="true" persistent="true"><math xmlns="{_MML}">'
                "<apply><geq/><ci>A</ci><cn>3</cn></apply></math></trigger>"
                f'<listOfEventAssignments><eventAssignment variable="k"><math xmlns="{_MML}">'
                "<cn>2</cn></math></eventAssignment></listOfEventAssignments></event>"
                "</listOfEvents>"
            )
        },
        "event_delay": {
            "model.xml": _sbml(
                events=f'<listOfEvents><event id="e" useValuesFromTriggerTime="true">'
                f'<trigger initialValue="true" persistent="true"><math xmlns="{_MML}">'
                f'{time_ge}</math></trigger><delay><math xmlns="{_MML}"><cn>1</cn>'
                '</math></delay><listOfEventAssignments><eventAssignment variable="k">'
                f'<math xmlns="{_MML}"><cn>2</cn></math></eventAssignment>'
                "</listOfEventAssignments></event></listOfEvents>"
            )
        },
        "event_reads_species": {
            "model.xml": _sbml(
                events="<listOfEvents>"
                + _event("e", 1, {"k": "<ci>A</ci>"})
                + "</listOfEvents>"
            )
        },
        "simultaneous_without_priority": {
            "model.xml": _sbml(
                events="<listOfEvents>"
                + _event("e1", 1, {"k": "<cn>2</cn>"})
                + _event("e2", 1, {"k": "<cn>3</cn>"})
                + "</listOfEvents>"
            )
        },
        "sbml_root_unknown_element": {
            "model.xml": _sbml().replace("<model ", "<extra/><model ", 1)
        },
        "time_zero_event_negative_start": {
            "study.sedml": _sed(start="-1"),
            "model.xml": _sbml(
                events="<listOfEvents>"
                + _event("e", 0, {"k": "<cn>2</cn>"})
                + "</listOfEvents>"
            ),
        },
    }
    result = {}
    for name, members in cases.items():
        path = _archive(temp / f"{name}.omex", members)
        result[name] = _refused(
            lambda path=path: crn.read_experiment(
                path, convention="stochastic_count_propensity"
            )
        )
    # Archive-level cases that need hand-edited bytes.
    symlink = zipfile.ZipInfo("model.xml")
    symlink.create_system = 3
    symlink.external_attr = 0o120777 << 16
    with zipfile.ZipFile(temp / "symlink.omex", "w") as arc:
        arc.writestr("manifest.xml", _MANIFEST)
        arc.writestr("study.sedml", _sed())
        arc.writestr(symlink, "/etc/hosts")
    raw = bytearray(_archive(temp / "bomb.omex", model).read_bytes())
    entry = raw.rfind(b"PK\x01\x02")  # last central record: model.xml
    raw[entry + 24 : entry + 28] = (65 << 20).to_bytes(4, "little")
    (temp / "bomb.omex").write_bytes(raw)
    raw = bytearray(_archive(temp / "corrupt.omex", model).read_bytes())
    start = raw.find(b"PK\x03\x04", raw.find(b"study.sedml")) + 30 + len("model.xml")
    raw[start : start + 16] = bytes(b ^ 0x5A for b in raw[start : start + 16])
    (temp / "corrupt.omex").write_bytes(raw)
    (temp / "folder").mkdir()
    (temp / "folder" / "study.sedml").write_text(_sed("absent.xml"))
    for name, path in (
        ("member_symlink", temp / "symlink.omex"),
        ("declared_size_bomb", temp / "bomb.omex"),
        ("corrupt_deflate", temp / "corrupt.omex"),
        ("folder_missing_model", temp / "folder" / "study.sedml"),
        ("directory_path", temp / "folder"),
    ):
        result[name] = _refused(
            lambda path=path: crn.read_experiment(
                path, convention="stochastic_count_propensity"
            )
        )
    # The single-file SBML importer applies the same UTF-8 rule. The UTF-8
    # form of this model imports; its BOM-less UTF-16 form with an entity
    # was imported with the entity expanded by 0.8.1.
    single = (
        f'<sbml xmlns="{_CORE}" level="3" version="2"><model id="m" '
        'substanceUnits="item" timeUnits="second" extentUnits="item">'
        '<listOfCompartments><compartment id="c" size="1" spatialDimensions="3" '
        'constant="true"/></listOfCompartments><listOfSpecies><species id="A" '
        'compartment="c" initialAmount="5" hasOnlySubstanceUnits="true" '
        'boundaryCondition="false" constant="false"/></listOfSpecies>'
        '<listOfParameters><parameter id="k" value="VALUE" constant="true"/>'
        '</listOfParameters><listOfReactions><reaction id="r" reversible="false">'
        '<listOfReactants><speciesReference species="A" stoichiometry="1" '
        'constant="true"/></listOfReactants><kineticLaw><math '
        f'xmlns="{_MML}"><apply><times/><ci>k</ci><ci>A</ci></apply></math>'
        "</kineticLaw></reaction></listOfReactions></model></sbml>"
    )
    convention = "stochastic_count_propensity"
    crn.importers.import_sbml(
        single.replace("VALUE", "0.25").encode(), convention=convention
    )
    hidden = '<?xml version="1.0" encoding="UTF-16"?><!DOCTYPE s [<!ENTITY e "0.25">]>'
    hidden += single.replace("VALUE", "&e;")
    result["importer_utf16_doctype"] = _refused(
        lambda: crn.importers.import_sbml(
            hidden.encode("utf-16-le"), convention=convention
        )
    )
    # Run-time budget and archive directory rules.
    plain = crn.read_experiment(
        _archive(temp / "plain.omex", model), convention="stochastic_count_propensity"
    )
    result["value_budget"] = _refused(lambda: plain.run(replicates=1 << 23))
    reserved = crn.read_experiment(
        _archive(temp / "reserved.omex", {"study.sedml": _sed(task="con"), **model}),
        convention="stochastic_count_propensity",
    )
    result["reserved_task_directory"] = _refused(
        lambda: reserved.run(archive_directory=temp / "reserved")
    )
    (temp / "exists").mkdir()
    try:
        plain.run(archive_directory=temp / "exists")
    except FileExistsError:
        result["existing_archive_directory"] = "FileExistsError"
    require(
        all(
            code.startswith("CRN-E090: ")
            for k, code in result.items()
            if "existing" not in k
        )
        and result["existing_archive_directory"] == "FileExistsError",
        f"experiment refusals: {result}",
    )
    return result


def experiment_semantics(temp):
    """Pair combinatorics, derived parameters after edits, and event timing."""
    temp = Path(temp)
    pair = '<speciesReference species="A" stoichiometry="2" constant="true"/>'
    falling = (
        "<apply><times/><ci>k</ci><ci>A</ci>"
        "<apply><minus/><ci>A</ci><cn>1</cn></apply></apply>"
    )
    study = crn.read_experiment(
        _archive(
            temp / "pair.omex",
            {
                "model.xml": _sbml(
                    law=falling,
                    reactants=pair,
                    initial="2",
                    parameters='<parameter id="k" value="0.5" constant="true"/>',
                )
            },
        ),
        convention="stochastic_count_propensity",
    )
    rate = study.preflight()["tasks"]["t"]["initial_rates"]["sbml_rate_0"]
    require(rate == 1.0, "k A (A - 1) is the hazard rate C(A, 2) x 2k")
    final = np.asarray(study.run(seed=1234, replicates=4000)["reports"]["out"]["A"])
    survived = float(np.mean(final[:, -1] == 2))
    require(abs(survived - np.exp(-1)) < 0.038, "pair hazard k A (A - 1)")
    # A derived rate follows a SED-ML change of the parameter it reads.
    derived = _sbml(
        parameters='<parameter id="k" constant="false"/>'
        '<parameter id="q" constant="false"/>'
        '<parameter id="p" value="1" constant="true"/>',
        rules=f'<listOfRules><assignmentRule variable="k"><math xmlns="{_MML}">'
        "<apply><times/><cn>2</cn><ci>q</ci></apply></math></assignmentRule>"
        f'<assignmentRule variable="q"><math xmlns="{_MML}"><ci>p</ci></math>'
        "</assignmentRule></listOfRules>",
    )
    edit = (
        "<listOfChanges><changeAttribute newValue='3' "
        "target=\"/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='p']"
        '/@value"/></listOfChanges>'
    )
    edited = crn.read_experiment(
        _archive(
            temp / "derived.omex",
            {"model.xml": derived, "study.sedml": _sed(changes=edit)},
        ),
        convention="stochastic_count_propensity",
    ).preflight()["tasks"]["t"]["initial_rates"]["sbml_rate_0"]
    require(edited == 6.0, "rule chain recalculated after an edit")
    # Events: time zero with both initial values, then two simultaneous
    # events whose second reads k at trigger time or at execution.
    zero = "<listOfEvents>{}</listOfEvents>"
    q = '<parameter id="k" value="1" constant="false"/><parameter id="q" value="1" constant="false"/>'
    law_q = "<apply><times/><ci>q</ci><ci>A</ci></apply>"
    timing = {}
    for name, events in {
        "zero_initial_true": _event("e", 0, {"A": "<cn>7</cn>"}),
        "zero_initial_false": _event("e", 0, {"A": "<cn>7</cn>"}, initial="false"),
        "same_time_trigger_values": _event("a", 0.5, {"k": "<cn>5</cn>"}, priority=2)
        + _event("b", 0.5, {"q": "<ci>k</ci>"}, priority=1),
        "same_time_execution_values": _event("a", 0.5, {"k": "<cn>5</cn>"}, priority=2)
        + _event("b", 0.5, {"q": "<ci>k</ci>"}, priority=1, extra="false"),
    }.items():
        loaded = crn.read_experiment(
            _archive(
                temp / f"{name}.omex",
                {
                    "model.xml": _sbml(
                        law=law_q, parameters=q, events=zero.format(events)
                    )
                },
            ),
            convention="stochastic_count_propensity",
        )
        task = loaded.preflight()["tasks"]["t"]
        first = np.asarray(loaded.run(seed=5, replicates=4)["reports"]["out"]["A"])
        timing[name] = {
            "events": [
                [e["time"], e["parameters"].get("sbml_rate_0"), e["counts"]]
                for e in task["events"]
            ],
            "inactive": task["events_inactive_at_time_zero"],
            "first_report": sorted({int(x) for x in first[:, 0]}),
        }
    require(
        timing["zero_initial_true"]["inactive"] == ["e"]
        and timing["zero_initial_true"]["first_report"] == [5]
        and timing["zero_initial_false"]["events"] == [[0.0, None, {"A": 7}]]
        and timing["zero_initial_false"]["first_report"] == [7]
        and [e[1] for e in timing["same_time_trigger_values"]["events"]] == [1.0, 1.0]
        and [e[1] for e in timing["same_time_execution_values"]["events"]]
        == [1.0, 5.0],
        f"SBML event timing: {timing}",
    )
    return {
        "pair_survival": survived,
        "derived_rate_after_edit": edited,
        "events": timing,
    }


def experiments(vilar):
    with TemporaryDirectory() as temp:
        refusals = experiment_refusals(temp)
        semantics = experiment_semantics(temp)
    with TemporaryDirectory() as temp:
        path = experiment_fixture(temp)
        loaded = crn.read_experiment(path, convention="stochastic_count_propensity")
        payload = loaded.run(
            seed=427, replicates=32, archive_directory=Path(temp) / "archive"
        )
        values = payload["reports"]["out"]["A_conc"]
        require(
            all(row[0] == 0 and row[1] == 5 for row in values),
            "pulse or concentration conversion",
        )
        require(
            crn.verify_archive(Path(temp) / "archive" / "t")["valid"],
            "archived experiment task",
        )
        m = crn.Model("pulse")
        m.add_species({"A": 0, "B": 0}).add_parameter("k", 0)
        m.add_reaction(crn.mass_action("k"), ["A"], ["B"], name="transfer")
        p = crn.Protocol().at(1, parameters={"k": 0.5}, set_counts={"A": 10})
        native = crn.ensemble(
            m,
            method="static_direct",
            protocol=p,
            end=2,
            report_times=[0, 1, 2],
            seed=427,
            replicates=32,
            keep_paths=32,
            quantiles=(),
        ).simulate()
        require(
            np.array_equal(
                np.asarray(values) * 2, native.paths.sel(variable="A").values
            ),
            "native and imported pulse differ",
        )
    checks = {"refusals": refusals, "semantics": semantics}
    if vilar is None:
        return {"pulse_replicates": len(values), **checks}
    expected = "d97b92fed06cc8b1d5cf152cdece7c30dee8458698f3c784876f68f0d4095c80"
    require(sha256(vilar.read_bytes()).hexdigest() == expected, "archive changed")
    try:
        crn.read_experiment(vilar, convention="stochastic_count_propensity")
    except ImportFailure:
        pass
    else:
        raise RuntimeError("signed reaction convention required for source archive")
    study = crn.read_experiment(vilar, convention="split_signed_mass_action")
    replay = study.run(seed=20, replicates=8)
    for k, report in replay["reports"].items():
        if k.startswith("report_"):
            require(
                all(
                    len(row) == 401
                    for v in report.values()
                    for row in (v if isinstance(v[0], list) else [v])
                ),
                "external report shape",
            )
    rates = [
        task["initial_rates"]["sbml_rate_5"]
        for task in study.preflight()["tasks"].values()
    ]
    require(
        abs(rates[0] / rates[1] - 4) < 1e-12,
        "condition change did not update derived rate",
    )
    changes = [t["parameter_changes"] for t in study.preflight()["tasks"].values()]
    require(changes == [{}, {"delta_R_R_decay": 0.05}], "recorded SED-ML change")
    return {
        "pulse_replicates": len(values),
        "vilar_reports": len(replay["reports"]),
        "effective_decay_rates": rates,
        **checks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vilar", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "protocols": timed_laws(),
        "delay_cdf": delayed_laws(),
        "spatial": spatial_laws(),
        "stationary": stationary_laws(),
        "histogram": observations(),
        "experiments": experiments(args.vilar),
    }
    result = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
