"""Compare the public experiment adapter with two external reference fixtures.

Inputs, pinned by SHA-256 (CI downloads them by commit; they are not in the
repository):

- the BioSimulators discrete-SSA oscillator archive (`--vilar`),
- `00001-sbml-l3v2.xml` and `00001-results.csv` of the SBML Test Suite
  stochastic case 00001 (`--dsmts`, a directory).
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import crn_craft as crn
from numba import njit
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from check_paper_capabilities import require

VILAR_SHA256 = "d97b92fed06cc8b1d5cf152cdece7c30dee8458698f3c784876f68f0d4095c80"
DSMTS_SHA256 = {
    "00001-sbml-l3v2.xml": "1ba61a28537918456018ee3fb162ffcc52775c6504ff4e97ff1847f16c9c5c3d",
    "00001-results.csv": "bc7164b674ac8788223f3a3f74401fc5d8f4920ad727845c7fd432b471511e52",
}


def pinned(path, expected):
    require(sha256(path.read_bytes()).hexdigest() == expected, f"{path.name} changed")


def stochastic_case(folder, threads):
    for name, digest in DSMTS_SHA256.items():
        pinned(folder / name, digest)
    raw = (folder / "00001-sbml-l3v2.xml").read_bytes()
    reference = pd.read_csv(folder / "00001-results.csv")
    with TemporaryDirectory() as temp:
        model = Path(temp) / "model.xml"
        model.write_bytes(raw)
        sedml = Path(temp) / "study.sedml"
        sedml.write_text("""<sedML xmlns="http://sed-ml.org/sed-ml/level1/version3"
          level="1" version="3">
        <listOfModels><model id="model" language="urn:sedml:language:sbml"
          source="model.xml"/></listOfModels>
        <listOfSimulations><uniformTimeCourse id="series" initialTime="0"
          outputStartTime="0" outputEndTime="50" numberOfPoints="50">
          <algorithm kisaoID="KISAO:0000029"/>
        </uniformTimeCourse></listOfSimulations>
        <listOfTasks><task id="task" modelReference="model"
          simulationReference="series"/></listOfTasks>
        <listOfDataGenerators><dataGenerator id="g"><listOfVariables>
          <variable id="x" taskReference="task"
            target="/sbml/model/listOfSpecies/species[@id='X']"/>
        </listOfVariables><math xmlns="http://www.w3.org/1998/Math/MathML">
          <ci>x</ci></math></dataGenerator></listOfDataGenerators>
        <listOfOutputs><report id="report"><listOfDataSets>
          <dataSet id="X" dataReference="g"/></listOfDataSets></report></listOfOutputs>
        </sedML>""")
        experiment = crn.read_experiment(
            sedml, convention="stochastic_count_propensity"
        )
        data = experiment.run(seed=843, replicates=5000, threads=threads)
    x = np.asarray(data["reports"]["report"]["X"], dtype=float)
    n = len(x)
    mu = reference["X-mean"].to_numpy()[1:]
    sigma = reference["X-sd"].to_numpy()[1:]
    # The suite's own statistics and ranges (00001-settings.txt): at each
    # report time after 0, Z = sqrt(n) (mean - mu) / sigma within (-3, 3)
    # and Y = sqrt(n / 2) (S^2 / sigma^2 - 1) within (-5, 5).
    z = np.sqrt(n) * (x.mean(axis=0)[1:] - mu) / sigma
    y = np.sqrt(n / 2) * (x.var(axis=0, ddof=1)[1:] / sigma**2 - 1)
    z_max, y_max = float(np.max(np.abs(z))), float(np.max(np.abs(y)))
    require(z_max < 3 and y_max < 5, "DSMTS case 00001 Z and Y ranges")
    return {
        "source_sha256": sha256(raw).hexdigest(),
        "replicates": n,
        "report_times": len(z),
        "max_abs_z": z_max,
        "max_abs_y": y_max,
    }


@njit(cache=False)
def oscillator(initial, needs, changes, rates, t, n, seed):
    np.random.seed(seed)
    out = np.zeros((n, len(t), len(initial)), dtype=np.int64)
    a = np.zeros(len(rates))
    for i in range(n):
        x = initial.copy()
        at = 0.0
        j = 0
        while j < len(t):
            for r in range(len(rates)):
                a[r] = rates[r]
                for s in needs[r]:
                    if s >= 0:
                        a[r] *= x[s]
            total = a.sum()
            next_time = at + np.random.exponential(1 / total) if total else np.inf
            while j < len(t) and t[j] < next_time:
                out[i, j] = x
                j += 1
            if j == len(t):
                break
            draw = np.random.random() * total
            r = 0
            while r < len(a) - 1 and draw >= a[r]:
                draw -= a[r]
                r += 1
            x += changes[r]
            at = next_time
    return out


def vilar(path, n, threads):
    pinned(path, VILAR_SHA256)
    t = np.array([0.5, 1, 2])
    experiment = crn.read_experiment(path, convention="split_signed_mass_action")
    # An independent literal reconstruction of the deposited archive's
    # stoichiometries and numeric rates. A and R are monomer counts.
    labels = (
        "A",
        "R",
        "C",
        "mRNA_A_",
        "mRNA_R",
        "PrmA",
        "PrmR",
        "PrmA_bound",
        "PrmR_bound",
    )
    statements = (
        (("C",), ("R",), 1),
        (("A",), (), 1),
        (("mRNA_A_",), (), 10),
        (("mRNA_R",), (), 0.5),
        (("A", "R"), ("C",), 2),
        (("R",), (), 0.2),
        (("PrmA",), ("mRNA_A_", "PrmA"), 50),
        (("PrmA_bound",), ("mRNA_A_", "PrmA_bound"), 500),
        (("PrmR",), ("mRNA_R", "PrmR"), 0.01),
        (("PrmR_bound",), ("mRNA_R", "PrmR_bound"), 50),
        (("PrmA", "A"), ("PrmA_bound",), 1),
        (("PrmR", "A"), ("PrmR_bound",), 1),
        (("mRNA_R",), ("mRNA_R", "R"), 5),
        (("mRNA_A_",), ("mRNA_A_", "A"), 50),
        (("PrmA_bound",), ("PrmA", "A"), 50),
        (("PrmR_bound",), ("PrmR", "A"), 100),
    )
    initial = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0], dtype=np.int64)
    needs = np.full((len(statements), 2), -1, dtype=np.int64)
    changes = np.zeros((len(statements), len(labels)), dtype=np.int64)
    rates = np.zeros(len(statements))
    for i, (reactants, products, rate) in enumerate(statements):
        req = [labels.index(s) for s in reactants]
        needs[i, : len(req)] = req
        changes[i, req] -= 1
        changes[i, [labels.index(s) for s in products]] += 1
        rates[i] = rate * 4.1887902047863905
    output = {}
    # Two-sample DKW bound for equal sizes, Bonferroni over 12 comparisons,
    # family-wise level 0.01: 0.1233 at 512 paths, 0.0299 at 8,704.
    bound = float(np.sqrt(np.log(2 * 12 / 0.01) / n))
    for i, task in enumerate(("tsk_1_0", "tsk_1_1")):
        expected = rates.copy()
        if i == 1:
            expected[5] *= 0.25
        oracle = oscillator(initial, needs, changes, expected, t, n, 5351 + i)
        # The imported task model, sampled at the three comparison times
        # only, so that large ensembles stay small in memory. The 401-point
        # reports are checked by check_paper_capabilities.py --vilar.
        actual = crn.ensemble(
            experiment.tasks[task][0].model,
            end=2.0,
            report_times=t,
            seed=4391 + i,
            replicates=n,
            threads=threads,
            method="static_direct",
            keep_paths=n,
            quantiles=(),
            strict=True,
        ).simulate()
        for label in ("A", "R"):
            trace = actual.paths.sel(variable=label).values
            statistics = [
                float(
                    ks_2samp(
                        trace[:, j], oracle[:, j, labels.index(label)], method="asymp"
                    ).statistic
                )
                for j in range(len(t))
            ]
            output[task + "/" + label] = statistics
    require(
        max(max(v) for v in output.values()) < bound,
        "published archive vs independent SSA",
    )
    return {
        "replicates_per_task": n,
        "time": t.tolist(),
        "cdf_gaps": output,
        "simultaneous_bound": bound,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vilar", type=Path, required=True)
    parser.add_argument("--dsmts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--vilar-paths",
        type=int,
        default=512,
        help="paths per task and sampler; 8704 resolves CDF gaps of 0.03",
    )
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    result = {
        "vilar": vilar(args.vilar, args.vilar_paths, args.threads),
        "dsmts_00001": stochastic_case(args.dsmts, args.threads),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
