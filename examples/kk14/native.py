"""Small native CRN-craft models to compare against the study references."""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

from reference import Config, atp_equilibrium, simulate, size_capacity


def drift_pair(mutation, n=1000, initial=0, half_life=240):
    """The ordered survivor-victim law from chapter 16c."""
    import crn_craft as crn

    model = crn.Model("mtdna_clones_pairs")
    model.add_parameter({"delta": math.log(2) / half_life, "pmut": mutation})
    model.add_species("mtDNA", attributes=["clone"])
    if n > initial:
        model.add_initial("mtDNA", n - initial, clone=0)
    if initial:
        model.add_initial("mtDNA", initial, clone=-1)
    model.add_observable("mutants", "mtDNA", where="x.clone != 0")
    model.add_observable("mut_types", "mtDNA", where="x.clone != 0", reduce="distinct")
    model.add_observable("largest", "mtDNA", where="x.clone != 0", reduce="max_count")
    pair_rate = "2 * delta / (mtDNA - 1)"
    model.add_reaction(
        crn.mass_action(f"{pair_rate} * (1 - pmut)"),
        ["mtDNA", "mtDNA"],
        [crn.copy("x"), crn.copy("x")],
        name="replacement",
    )
    model.add_reaction(
        crn.mass_action(f"{pair_rate} * pmut"),
        ["mtDNA", "mtDNA"],
        [crn.copy("x"), crn.new("mtDNA", clone="new_id()")],
        name="new_lineage",
    )
    return model


def staged_size(stages=4, n=100, initial=10, half_life=48.0):
    """An Erlang approximation that runs on native Static Direct."""
    import crn_craft as crn

    if stages < 1:
        raise ValueError("the chain needs at least one busy stage")
    model = crn.Model(f"mtdna_size_erlang_{stages}")
    wild = ["W", *(f"W{i}" for i in range(1, stages + 1))]
    mutant = ["M", *(f"M{i}" for i in range(1, stages + 1))]
    species = dict.fromkeys(wild + mutant, 0)
    species.update(W=n - initial, M=initial)
    model.add_species(species)
    model.add_parameter(
        {
            "C": size_capacity(n, half_life, 2, 1),
            "delta": math.log(2) / half_life,
        }
    )
    total = "+".join(species)
    for labels, delay in ((wild, 2.0), (mutant, 1.0)):
        free, first = labels[:2]
        model.add_reaction(
            f"C * {free} / (C + {total})", [free], [first], name=f"start_{free}"
        )
        for origin, dest in zip(labels[1:-1], labels[2:], strict=True):
            model.add_reaction(
                stages / delay, [origin], [dest], name=f"advance_{origin}"
            )
        model.add_reaction(
            stages / delay,
            [labels[-1]],
            [free, free],
            name=f"finish_{free}",
        )
        for label in labels:
            model.add_reaction(
                crn.mass_action("delta"), [label], [], name=f"death_{label}"
            )
    return model


def coupled_atp(
    mutation, n=1000, initial=0, half_life=240, *, legacy=False, wild_only=False
):
    """Two genotypes sharing the same ATP pool, in either native solver."""
    import crn_craft as crn

    delta, atp, c = atp_equilibrium(n, half_life)
    if not legacy:
        model = crn.Model("mtdna_two_class_atp")
        model.add_species({"W": n - initial, "M": initial})
        model.add_parameter(
            {
                "delta": delta,
                "c": c,
                "mu": mutation,
                "mumut": 0.0 if wild_only else mutation,
            }
        )
        model.add_ode("ATP", atp, "0.1 * W - 0.2 * ATP - 0.01 * (W + M)", lower_bound=0)
        rate = "c / (c + ATP)"
        for parent, probability, synthesis in (("W", "mu", 1.0), ("M", "mumut", 1.5)):
            model.add_reaction(
                crn.mass_action(f"{synthesis} * ({rate}) * (1 - {probability})"),
                [parent],
                [parent, parent],
                name=f"copy_{parent}",
            )
            model.add_reaction(
                crn.mass_action(f"{synthesis} * ({rate}) * {probability}"),
                [parent],
                [parent, "M"],
                name=f"mutate_{parent}",
            )
            model.add_reaction(
                crn.mass_action("delta"), [parent], [], name=f"death_{parent}"
            )
        return model
    # The legacy ATP block takes two observables. A single attributed
    # species supplies W, M and the total with the same reaction laws.
    legacy_model = crn.Model("mtdna_two_class_thinning")
    legacy_model.add_parameter(
        {"delta": delta, "mu": mutation, "mumut": 0 if wild_only else mutation}
    )
    legacy_model.add_species("mtDNA", attributes=["mutant"])
    legacy_model.add_initial("mtDNA", n - initial, mutant=0)
    if initial:
        legacy_model.add_initial("mtDNA", initial, mutant=1)
    for group, attr, prob, speed in (("W", 0, "mu", 1.0), ("M", 1, "mumut", 1.5)):
        legacy_model.add_observable(group, "mtDNA", where=f"x.mutant == {attr}")
        parent = crn.pattern("mtDNA", where=f"x.mutant == {attr}")
        legacy_model.add_reaction(
            crn.mass_action(f"{speed} * (1 - {prob})"),
            [parent],
            [crn.copy("x"), crn.copy("x")],
            name=f"copy_{group}",
        )
        legacy_model.add_reaction(
            crn.mass_action(f"{speed} * {prob}"),
            [parent],
            [crn.copy("x"), crn.new("mtDNA", mutant=1)],
            name=f"mutate_{group}",
        )
        legacy_model.add_reaction(
            crn.mass_action("delta"), [parent], [], name=f"death_{group}"
        )
    doc = legacy_model.to_dict()
    doc["schema"] = "crn-craft.rules-model/0.2"
    doc["atp"] = {
        "initial": atp,
        "production": 0.1,
        "decay": 0.2,
        "maintenance": 0.01,
        "scale": c,
        "wild_type": "W",
        "total": "mtDNA",
        "replication_rules": ["copy_W", "mutate_W", "copy_M", "mutate_M"],
    }
    return doc


def marked_hybrid(n=100, initial=1, mutation=0.01, hours=240, reps=8, threads=2):
    """Mark aggregate mutant events without altering the native random stream."""
    import crn_craft as crn

    model = coupled_atp(mutation, n=n, initial=initial)
    names = [r["name"] for r in model.to_dict()["rules"]]
    pools = [Counter({0: initial}) if initial else Counter() for _ in range(reps)]
    next_ids = [1] * reps
    rngs = [np.random.default_rng(63188 + i) for i in range(reps)]

    def mark(replicate, event):
        counter = pools[replicate]
        rule = names[event["rule"]]
        if rule.startswith("mutate_"):
            counter[next_ids[replicate]] = 1
            next_ids[replicate] += 1
        elif rule in ("copy_M", "death_M"):
            rank = rngs[replicate].integers(sum(counter.values()))
            for label, count in counter.items():
                if rank < count:
                    clone = label
                    break
                rank -= count
            else:
                raise RuntimeError("native mutant event has no labelled template")
            if rule == "copy_M":
                counter[clone] += 1
            else:
                counter[clone] -= 1
                if not counter[clone]:
                    del counter[clone]
        # The two mass-action rules with an M daughter have different
        # identities even when their net change is identical.
        if sum(counter.values()) != event["observables"][mutant_column]:
            raise RuntimeError("marked mutants disagree with the native event")
        return True

    run = crn.ensemble(
        model,
        end=hours,
        report_times=[0, hours],
        replicates=reps,
        seed=207,
        threads=threads,
        method="hybrid",
        keep_paths="all",
        strict=True,
    )
    plain = run.simulate()
    mutant_column = plain.paths.coords["variable"].to_numpy().tolist().index("M")
    marked = run.simulate(callback=mark)
    same = bool(np.array_equal(plain.paths, marked.paths))
    if not same:
        raise RuntimeError("marking changed native demographic trajectories")
    return {
        "replicates": reps,
        "matched_counts": same,
        "events": int(marked.records["events"].sum()),
        "mutant_counts": [sum(c.values()) for c in pools],
        "clone_counts": [len(c) for c in pools],
        "sampler": marked.metadata["sampler"],
        "method": marked.metadata["method"],
    }


def _end(model, method, hours, reps, threads, seed):
    import crn_craft as crn

    return crn.ensemble(
        model,
        end=hours,
        report_times=[0, hours],
        replicates=reps,
        seed=seed,
        threads=threads,
        method=method,
        keep_paths="all",
        strict=True,
    ).simulate()


def _mean_se(values):
    values = np.asarray(values, dtype=float)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(len(values)))


def staged_comparison(threads=2, reps=320, seed=None, stages=(1, 4, 16)):
    """Endpoint mutant fractions of the hourly queue and of Erlang chains.

    Without ``seed`` the seeds of the saved 320-cell table are used. With a
    seed, the hourly queue takes ``seed, seed + 1, ...`` and the chain with
    ``k`` stages takes ensemble seed ``seed + 1000 * k``.
    """
    n, m0, hours, half_life = 100, 10, 24 * 120, 48
    cfg = Config(
        "size", hours, 0, population=n, initial_mutants=m0, half_life=half_life
    )
    first = 1010 if seed is None else seed
    base = simulate(
        cfg, np.arange(first, first + reps), threads=threads, labelled=False
    )
    total = base["wild"] + base["mutant"]
    delta = math.log(2) / half_life
    mean, se = _mean_se(base["mutant"] / total)
    rows = [
        {
            "sampler": "fixed-delay hourly",
            "stages": 0,
            "cells": reps,
            "mean_mutant_fraction": mean,
            "mean_mutant_fraction_se": se,
            "mutant_delay_variance_h2": 0.0,
            "mutant_template_survival_continuous": math.exp(-delta),
            "mutant_template_survival_sampler": 1.0,
        }
    ]
    # With a one-hour delay a new template completes at the next tick
    # before degradation draws for that tick.
    for k in stages:
        result = _end(
            staged_size(k, n, m0, half_life),
            "static_direct",
            hours,
            reps,
            threads,
            420 + k if seed is None else seed + 1000 * k,
        )
        final = result.paths.isel(time=-1)
        w = final.sel(variable=["W", *(f"W{i}" for i in range(1, k + 1))]).sum(
            "variable"
        )
        m = final.sel(variable=["M", *(f"M{i}" for i in range(1, k + 1))]).sum(
            "variable"
        )
        p = (1 + delta / k) ** -k
        mean, se = _mean_se(m / (m + w))
        rows.append(
            {
                "sampler": result.metadata["sampler"],
                "stages": k,
                "cells": reps,
                "mean_mutant_fraction": mean,
                "mean_mutant_fraction_se": se,
                "mutant_delay_variance_h2": 1.0 / k,
                "mutant_template_survival_continuous": p,
                "mutant_template_survival_sampler": p,
            }
        )
    return rows


def atp_comparison(threads=2, reps=320, seed=None):
    """COX-negative cells and ATP stops of hybrid, thinning and the hourly scheme.

    Without ``seed`` the seeds of the saved 320-cell table are used. With a
    seed, hybrid takes ensemble seed ``seed + 1``, thinning ``seed + 2`` and
    the hourly scheme the cell seeds ``seed + 10000, seed + 10001, ...``.
    """
    # Reduced population and a single founder give a mixture of loss and takeover.
    # 8,760 h is a fixed test horizon for this small cell, not a year of the
    # tissue study (which counts 360 days per year).
    n, initial, hours, half_life = 100, 1, 8760, 120
    mutation = 0.0001
    rows = []
    for method in ("hybrid", "thinning"):
        model = coupled_atp(
            mutation, n, initial, half_life, legacy=method == "thinning"
        )
        result = _end(
            model,
            method,
            hours,
            reps,
            threads,
            (521 if method == "hybrid" else 644)
            if seed is None
            else seed + (1 if method == "hybrid" else 2),
        )
        final = result.paths.isel(time=-1)
        w = final.sel(variable="W").to_numpy()
        m = final.sel(variable="M").to_numpy()
        collapsed = int((result.records["termination"] == "stopped").sum())
        rows.append(
            {
                "method": method,
                "sampler": result.metadata["sampler"],
                "cells": reps,
                "cox_negative": int(((w + m > 0) & (5 * m >= 3 * (w + m))).sum()),
                "collapsed": collapsed,
            }
        )
    for step in (1.0, 0.5):
        cfg = Config(
            "transcription",
            hours,
            mutation,
            population=n,
            initial_mutants=initial,
            half_life=half_life,
            step=step,
        )
        first = 3300 if seed is None else seed + 10000
        batch = simulate(
            cfg, np.arange(first, first + reps), labelled=False, threads=threads
        )
        total = batch["wild"] + batch["mutant"]
        rows.append(
            {
                "method": f"hourly reconstruction, step {step:g} h",
                "sampler": "kk14-hourly",
                "cells": reps,
                "cox_negative": int(
                    ((total > 0) & (5 * batch["mutant"] >= 3 * total)).sum()
                ),
                "collapsed": int((batch["status"] == 1).sum()),
            }
        )
    return rows
