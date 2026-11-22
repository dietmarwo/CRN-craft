"""Checks of C7 (Summit SnAr reactor) and of the concentration solver client."""

import csv
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.stats import mannwhitneyu

from publications import reactor
from publications.concentrations import CaseError, Solver, document

HERE = Path(__file__).parent
RESULTS = HERE / "results" / "C7"
ROOT = HERE.parents[1]
BINARY = Path(
    os.environ.get(
        "CRN_DIFFSOL",
        ROOT / "integrations/diffsol/target/release/crn-diffsol",
    )
)
DATA = Path(os.environ.get("CRN_PUBLICATION_DATA", ROOT / ".publication-data"))


def brute_hypervolume(points, ref):
    """Union of dominated rectangles by slabs between distinct abscissae."""
    inside = [p for p in points if p[0] < ref[0] and p[1] < ref[1]]
    xs = sorted({p[0] for p in inside} | {ref[0]})
    area = 0.0
    for x0, x1 in zip(xs[:-1], xs[1:], strict=True):
        low = min((p[1] for p in inside if p[0] <= x0), default=ref[1])
        area += (x1 - x0) * (ref[1] - low)
    return area


class Hypervolume(unittest.TestCase):
    def test_hand_front(self):
        # (3, 4) is dominated by (2, 2); (6, 0) lies outside the rectangle.
        points = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 4.0], [6.0, 0.0]])
        self.assertEqual(reactor.hypervolume(points, np.array([5.0, 5.0])), 11.0)
        self.assertEqual(reactor.hypervolume(points[3:], np.array([5.0, 5.0])), 0.0)

    def test_against_slabs(self):
        rng = np.random.default_rng(3)
        for _ in range(50):
            points = rng.uniform(-1, 1, (rng.integers(1, 30), 2))
            self.assertAlmostEqual(
                reactor.hypervolume(points, np.array([0.5, 0.5])),
                brute_hypervolume(points, (0.5, 0.5)),
                places=12,
            )

    def test_nondominated(self):
        points = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 4.0], [1.0, 5.0]])
        self.assertEqual(sorted(reactor.nondominated(points)), [0, 1])


class Emulator(unittest.TestCase):
    def test_rates_and_kinetics_variants(self):
        a0, b0, ka, kb, kc, kd, kq = reactor.inputs((1.0, 2.0, 0.2, 90.0))
        self.assertEqual((a0, b0), (0.2, 0.4))
        np.testing.assert_allclose(
            (ka, kb, kc, kd), 0.6 * np.array([57.9, 2.70, 0.865, 1.63])
        )
        self.assertEqual(kq, ka)
        self.assertEqual(
            reactor.inputs((1.0, 2.0, 0.2, 60.0), "pinned")[-1],
            reactor.inputs((1.0, 2.0, 0.2, 60.0), "pinned")[3],
        )

    def test_objectives(self):
        outlet = np.array([0.01, 0.2, 0.1, 0.02, 0.005])
        sty, waste = reactor.metrics(2.0, outlet)
        product = 210.21 * 0.1
        self.assertAlmostEqual(-sty, 60 * product / 2.0)
        other = 159.09 * 0.01 + 71.12 * 0.2 + 210.21 * 0.02 + 261.33 * 0.005
        self.assertAlmostEqual(waste, (789.0 + other) / product)
        self.assertEqual(list(reactor.metrics(1.0, np.zeros(5))), [-1e-6, 1000.0])

    def test_archived_summit_outputs(self):
        path = DATA / "C7" / "pareto_front_snar.csv"
        if not path.is_file():
            self.skipTest(f"fixture {path} absent; run publications/sources.py C7")
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))[1:]
        for row in rows[:: len(rows) // 25][:25]:
            design = [
                float(row[k]) for k in ("tau", "equiv_pldn", "conc_dfnb", "temperature")
            ]
            ours = reactor.metrics(
                design[0], reactor.outlets_reference(design, loose=True)
            )
            stored = (-float(row["sty"]), float(row["e_factor"]))
            np.testing.assert_allclose(ours, stored, rtol=1e-6)


class SavedResults(unittest.TestCase):
    def test_terminal_table(self):
        with (RESULTS / "terminal.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        stats = json.loads((RESULTS / "stats.json").read_text())
        arms = {}
        for row in rows:
            self.assertEqual(row["evaluations"], row["budget"])
            self.assertEqual(row["failures"], "0")
            arms.setdefault(row["arm"], []).append(float(row["hypervolume"]))
        self.assertEqual(
            {k: len(v) for k, v in arms.items()},
            {"mode": 20, "random": 20, "mode-500": 20},
        )
        for arm, values in arms.items():
            self.assertAlmostEqual(
                np.median(values), stats["archive"]["terminal"][arm]["median"]
            )
        p = mannwhitneyu(arms["mode"], arms["random"], alternative="two-sided").pvalue
        self.assertAlmostEqual(p, stats["archive"]["mode_vs_random"]["p_two_sided"])

    def test_figures_draw(self):
        with TemporaryDirectory() as folder:
            reactor.draw(RESULTS, Path(folder))
            for name in ("C7-pareto.png", "C7-hypervolume.png"):
                self.assertLess((Path(folder) / name).stat().st_size, 250_000)


@unittest.skipUnless(BINARY.is_file(), "crn-diffsol not built; set CRN_DIFFSOL")
class ConcentrationSolver(unittest.TestCase):
    def setUp(self):
        self.folder = TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)

    def solver(self, model, **options):
        return Solver(BINARY, Path(self.folder.name) / "model.json", model, **options)

    def test_decay_and_transport(self):
        model = document(
            ("C1", "C2"),
            ("k", "V1", "V2"),
            ("2", "0.5"),
            ("-J/V1", "J/V2"),
            assignments=(("J", "k*(C1-C2)"),),
            observations=("V1*C1+V2*C2",),
        )
        times = np.linspace(0, 10, 21)
        with self.solver(model) as engine:
            states, amount = engine.batch(
                [(0.8, 1.0, 4.0)], times, rtol=1e-10, atol=1e-13
            )
        mean = (2.0 + 0.5 * 4.0) / 5.0
        exact = mean + (2.0 - mean) * np.exp(-0.8 * 1.25 * times)
        np.testing.assert_allclose(states[0, :, 0], exact, atol=1e-8)
        np.testing.assert_allclose(amount[0, :, 0], 4.0, rtol=1e-12)

    def test_failed_case_is_data(self):
        model = document(("A",), ("k",), ("1",), ("-k*A",), observations=("A/k",))
        with self.solver(model) as engine:
            with self.assertRaises(CaseError) as caught:
                engine.batch([(1.0,), (0.0,)], [0.0, 1.0])
            self.assertEqual(caught.exception.failures[0][:2], (1, "CRN-E080"))
            with self.assertRaises(ValueError):
                engine.batch([(1.0,)], [0.0, 1.0, 1.0])
            states, _ = engine.batch([(1.0,)], [0.0, 1.0])  # the worker survives
            self.assertAlmostEqual(states[0, 1, 0], np.exp(-1.0), places=7)

    def test_refused_model_reports_diagnostic(self):
        with self.solver(document(("A",), ("k",), ("1",), ("-q*A",))) as engine:
            with self.assertRaisesRegex(RuntimeError, "CRN-E080"):
                engine.batch([(1.0,)], [0.0, 1.0])


if __name__ == "__main__":
    unittest.main()
