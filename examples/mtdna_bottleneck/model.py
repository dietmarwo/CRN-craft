"""Linear birth–death phases and deterministic times of binomial partition."""

import math

import crn_craft as crn
import numpy as np

try:
    from numba import njit
except ImportError:

    def njit(function):
        return function


def build(initial, birth, death):
    model = crn.Model("birth_death_partition")
    model.add_species({"W": int(initial[0]), "M": int(initial[1])})
    for s in ("W", "M"):
        model.add_reaction(birth, [s], [s, s], name=f"replicate_{s}")
        model.add_reaction(death, [s], [], name=f"degrade_{s}")
    return model


def native_lineage(initial, phases, seed):
    """Advance one lineage using the Python builder and seeded external division."""
    rng = np.random.default_rng(seed)
    counts = np.array(initial, dtype=int)
    saved = []
    for hours, birth, death, divide in phases:
        model = build(counts, birth, death)
        run = crn.ensemble(
            model,
            end=hours,
            report_times=[hours],
            keep_paths=1,
            method="static_direct",
            strict=True,
            seed=int(rng.integers(2**63)),
        )
        counts = run.simulate().paths.sel(variable=["W", "M"]).values[0, 0].astype(int)
        if divide:
            counts = rng.binomial(counts, 0.5)
        saved.append(counts.copy())
    return np.asarray(saved)


def birth_death_draw(n, birth, death, hours, rng):
    """Sample the independent founder families, each zero or geometric."""
    net = birth - death
    growth = math.exp(net * hours)
    integral = -math.expm1(-net * hours) / net if net else hours
    # 1 + birth * (exp(net*t)-1)/net is the conditional clone size.
    geometric = 1 / (1 + birth * growth * integral)
    survival = min(1.0, growth * geometric)
    alive = rng.binomial(n, survival)
    return alive + rng.negative_binomial(np.maximum(1, alive), geometric) * (alive > 0)


@njit
def predictions(theta, ages):
    """Exact log mean and leading normalized heteroplasmy variance.

    theta holds initial copies, six rate differences, six removal rates
    and three quiescent phase lengths. Phase six has zero growth and
    removal as in the steady terminal-oocyte convention.
    """
    copy_log = math.log(theta[0])
    growth, death = theta[1:7], theta[7:13]
    ends = np.array(
        [
            203.0,
            315.0,
            315.0 + 24 * theta[13],
            315.0 + 24 * (theta[13] + theta[14]),
            315.0 + 24 * (theta[13] + theta[14] + theta[15]),
            np.inf,
        ]
    )
    result = np.empty((len(ages), 2))
    current, phase, divisions, next_division, z = 0.0, 0, 0, 7.0, 0.0
    peak = copy_log
    for i in range(len(ages)):
        target = ages[i] * 24
        while current < target:
            end = min(target, ends[phase], next_division)
            dt = end - current
            g = growth[phase]
            d = death[phase]
            increment = -math.expm1(-g * dt) / g if g else dt
            z += (g + 2 * d) * math.exp(-copy_log) * increment
            copy_log += g * dt
            peak = max(peak, copy_log)
            if end == next_division:
                copy_log -= math.log(2)
                z += math.exp(-copy_log) / 2
                divisions += 1
                if divisions < 29:
                    next_division += 7
                elif divisions < 36:
                    next_division += 16
                else:
                    next_division = np.inf
            if end == ends[phase] and phase < 5:
                phase += 1
            current = end
        result[i, 0] = copy_log
        result[i, 1] = z
    return result, peak


@njit
def accept_distance(theta, ages, kinds, measured):
    for i in range(6):
        if not 0 <= theta[7 + i] <= 1 or not 0 <= theta[1 + i] + theta[7 + i] <= 1:
            return 1e20
    if (
        theta[0] <= 0
        or theta[0] > 1e6
        or min(theta[13:16]) < 0
        or max(theta[13:16]) > 50
    ):
        return 1e20
    theory, peak = predictions(theta, ages)
    if peak >= math.log(500000) or np.max(theory[:, 1]) > 1:
        return 1e20
    score = 0.0
    for i in range(len(kinds)):
        j = kinds[i]
        score += (theory[i, j] - measured[i]) ** 2 * (1000 if j else 1)
    return score


@njit
def chain(initial, steps, changes, indices, ages, kinds, measured, threshold, seed):
    np.random.seed(seed)
    thin = 10
    positions = np.empty((steps // thin, len(initial)))
    scores = np.empty(steps // thin)
    x = initial.copy()
    score = accept_distance(x, ages, kinds, measured)
    accepted = 0
    for i in range(steps):
        direction = indices[i % len(indices)]
        candidate = x.copy()
        candidate[direction] += changes[direction] * np.random.randn()
        candidate_score = accept_distance(candidate, ages, kinds, measured)
        if candidate_score <= threshold:
            x, score = candidate, candidate_score
            accepted += 1
        if i % thin == thin - 1:
            positions[i // thin] = x
            scores[i // thin] = score
    return positions, scores, accepted


def forecast_sample(initial_count, h, log_mean, relative_var, n, rng):
    """Exact one-age founder distribution from two count moments."""
    founders = max(1, int(round(initial_count)))
    mean = math.exp(log_mean) / founders
    var = relative_var * founders * mean * mean
    geometric = min(1.0, 2 * mean / (var + mean + mean * mean))
    survival = min(1.0, mean * geometric)
    initials = np.array([round(founders * (1 - h)), 0], dtype=int)
    initials[1] = founders - initials[0]
    families = rng.binomial(initials, survival, size=(n, 2))
    counts = families + rng.negative_binomial(np.maximum(1, families), geometric) * (
        families > 0
    )
    total = counts.sum(axis=1)
    return np.divide(counts[:, 1], total, out=np.full(n, np.nan), where=total > 0)


@njit
def ensemble_chain(
    initial, sweeps, changes, indices, ages, kinds, measured, threshold, seed
):
    """Symmetric differential evolution proposals on the hard ABC target."""
    np.random.seed(seed)
    width, dims = initial.shape
    every = 10
    states = np.empty((sweeps // every, width, dims))
    distances = np.empty((sweeps // every, width))
    population = initial.copy()
    scores = np.array([accept_distance(x, ages, kinds, measured) for x in population])
    accepted = 0
    scales = (0.06, 0.12, 0.23, 0.4, 0.08, 0.16, 0.27, 0.6)
    for sweep in range(sweeps):
        for i in range(width):
            a = np.random.randint(width - 1)
            if a >= i:
                a += 1
            b = np.random.randint(width)
            while b == i or b == a:
                b = np.random.randint(width)
            jump = scales[(sweep + i) % len(scales)]
            proposal = population[i] + jump * (population[a] - population[b])
            for k in indices:
                proposal[k] += 0.025 * changes[k] * np.random.randn()
            score = accept_distance(proposal, ages, kinds, measured)
            if score <= threshold:
                population[i], scores[i] = proposal, score
                accepted += 1
        if sweep % every == every - 1:
            states[sweep // every] = population
            distances[sweep // every] = scores
    return states, distances, accepted
