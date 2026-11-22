"""Checks of C10 (the 36-parameter three-step pathway)."""

import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from publications import pathway
from publications.concentrations import CaseError, Solver
from publications.search import PENALTY

HERE = Path(__file__).parent
RESULTS = HERE / "results" / "C10"
BINARY = Path(
    os.environ.get(
        "CRN_DIFFSOL",
        HERE.parents[1] / "integrations/diffsol/target/release/crn-diffsol",
    )
)


def rows(name):
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def hand_rhs(x, theta, s, p, literal=False):
    """The eight equations written out term by term (Mendes 2001 model)."""
    g1, g2, g3, e1, e2, e3, m1, m2 = x
    t = dict(zip(pathway.NAMES, theta, strict=True))

    def gene(i, regulator):
        return t[f"V{i}"] / (
            1
            + (p / t[f"Ki{i}"]) ** t[f"ni{i}"]
            + (t[f"Ka{i}"] / regulator) ** t[f"na{i}"]
        )

    def flux(i, a, b):
        km_a, km_b = t[f"Km{2 * i - 1}"], t[f"Km{2 * i}"]
        return t[f"kcat{i}"] * e[i - 1] / km_a * (a - b) / (1 + a / km_a + b / km_b)

    e = (e1, e2, e3)
    return np.array(
        [
            gene(1, s) - t["k1"] * g1,
            gene(2, m1) - t["k2"] * g2,
            gene(3, m2) - t["k3"] * g3,
            t["V4"] * g1 / (t["K4"] + g1) - t["k4"] * e1,
            t["V5"] * g2 / (t["K5"] + g2) - t["k5"] * e2,
            t["V6"] * (g2 if literal else g3) / (t["K6"] + g3) - t["k6"] * e3,
            flux(1, s, m1) - flux(2, m1, m2),
            flux(2, m1, m2) - flux(3, m2, p),
        ]
    )


class Objective(unittest.TestCase):
    def test_weights(self):
        self.assertEqual(
            pathway.weighted_error([1.0, 2.0], [2.0, 4.0], [2.0, 4.0]), 0.5
        )
        self.assertEqual(pathway.weighted_error([1.0], [1.0], [3.0]), 0.0)

    def test_nominal_values_and_conditions(self):
        theta = dict(zip(pathway.NAMES, pathway.TRUTH, strict=True))
        self.assertEqual(len(pathway.NAMES), 36)
        self.assertEqual(
            (theta["ni2"], theta["V5"], theta["k6"], theta["Km6"]), (2.0, 0.1, 0.1, 1.0)
        )
        self.assertEqual(len(set(pathway.CONDITIONS)), 16)

    def test_best_so_far_skips_failures(self):
        scores = np.array([PENALTY, 5.0, PENALTY, 3.0, 4.0])
        best = [r[-1] for r in pathway.best_so_far(scores, ("x",))]
        self.assertEqual(best[0], np.inf)
        self.assertEqual(best[-1], 3.0)
        self.assertTrue(np.all(np.diff(best[1:]) <= 0))

    def test_reference_against_hand_equations(self):
        # The reference integrator must follow the written-out equations.
        rng = np.random.default_rng(5)
        theta = pathway.TRUTH * np.exp(rng.uniform(-0.3, 0.3, 36))
        times = np.array([0.0, 6.0])
        for literal in (False, True):
            solved = pathway.reference(theta, times, literal)
            s, p = pathway.CONDITIONS[5]
            hand = solve_ivp(
                lambda _, x, s=s, p=p, literal=literal: hand_rhs(
                    x, theta, s, p, literal
                ),
                [0, 6.0],
                pathway.INITIAL,
                method="LSODA",
                rtol=1e-11,
                atol=1e-12,
            ).y[:, -1]
            np.testing.assert_allclose(solved[5, -1], hand, atol=1e-8)


class SavedResults(unittest.TestCase):
    def test_runs_match_reported_numbers(self):
        table = {
            (r["stage"], r["region"], r["method"], r["retry"]): r
            for r in rows("runs.csv")
        }
        de = table[("20000", "broad", "de", "0")]
        cma = table[("20000", "broad", "cmaes", "0")]
        self.assertAlmostEqual(float(de["best_proposed"]), 601.918, places=3)
        self.assertAlmostEqual(float(de["verified"]), 959.6845, places=4)
        self.assertEqual((de["failures"], de["rejected_finalists"]), ("2612", "24"))
        self.assertAlmostEqual(float(cma["verified"]), 667.4722, places=4)
        self.assertEqual(cma["failures"], "539")
        best = min(float(r["verified"]) for r in table.values() if r["verified"])
        self.assertAlmostEqual(best, 62.148705, places=5)

    def test_local_fits_reach_the_floor(self):
        for row in rows("local.csv"):
            self.assertLess(float(row["final"]), 1e-10)
            self.assertLess(float(row["max_log10_error"]), 1e-5)

    def test_nominal_table(self):
        table = rows("nominal.csv")
        self.assertEqual(len(table), 16 * 21)
        first = [float(table[0][s]) for s in pathway.STATES]
        np.testing.assert_allclose(first, pathway.INITIAL)

    def test_figures_draw(self):
        with TemporaryDirectory() as folder:
            pathway.draw(RESULTS, Path(folder))
            for name in ("C10-best-so-far.png", "C10-nominal.png"):
                self.assertLess((Path(folder) / name).stat().st_size, 250_000)


@unittest.skipUnless(BINARY.is_file(), "crn-diffsol not built; set CRN_DIFFSOL")
class NativeModel(unittest.TestCase):
    def test_native_matches_reference(self):
        times = np.arange(0.0, 31.0, 6.0)
        theta = pathway.TRUTH * np.exp(np.linspace(-0.2, 0.2, 36))
        with TemporaryDirectory() as folder:
            with Solver(
                BINARY, Path(folder) / "m.json", pathway.model(False)
            ) as engine:
                native, _ = engine.batch(
                    [(*theta, s, p) for s, p in pathway.CONDITIONS],
                    times,
                    rtol=1e-10,
                    atol=1e-12,
                )
                with self.assertRaises(CaseError):
                    engine.batch([(*np.full(36, 1e12), 0.1, 0.05)], times)
        np.testing.assert_allclose(native, pathway.reference(theta, times), atol=1e-7)


if __name__ == "__main__":
    unittest.main()
