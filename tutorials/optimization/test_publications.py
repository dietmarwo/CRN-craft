"""Regression checks for the publication tutorial scoring and replay."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

try:
    import crn_craft as crn
    from crn_craft.sensitivities import forward_sensitivities
    from publications.experiment_design import analytic_check
    from publications.pathway import weighted_error
    from publications.petab import gaussian_nll
    from publications.reactor import hypervolume
    from publications.rna_counts import CdfScore
    from publications.search import search
    from publications.stimulus_fit import information
except ImportError as error:  # pragma: no cover - depends on the environment
    raise unittest.SkipTest(
        f"publication tutorials need the packages of tutorials/optimization/README.md: {error}"
    ) from None


class PublicationTests(unittest.TestCase):
    def test_cdf_repeated_thresholds(self):
        sample = np.array([[0, 0], [0, 1], [0, 1], [0, 2]])
        score = CdfScore(sample, midrank=True)
        self.assertAlmostEqual(score([[0.5, 1.0, 1.0], [0.125, 0.5, 0.875]]), 0.0)
        self.assertGreater(score([[0.6, 1.0, 1.0], [0.125, 0.5, 0.875]]), 0.0)

    def test_units_and_scoring(self):
        self.assertEqual(
            hypervolume(
                np.array([[1.0, 3.0], [2.0, 2.0], [6.0, 0.0]]), np.array([5.0, 5.0])
            ),
            11.0,
        )
        self.assertAlmostEqual(
            float(gaussian_nll(3.0, 3.0, 2.0) - gaussian_nll(3.0, 3.0, 1.0)),
            np.log(2.0),
        )
        self.assertEqual(weighted_error([1.0], [2.0], [2.0]), 0.25)
        self.assertEqual(information(6.0, 2, 1), (16.0, 12.0))

    def test_forward_sensitivity_with_and_without_cutoff(self):
        report = {"checks": []}
        analytic_check(report)
        self.assertTrue(report["checks"][0]["passed"], report)
        model = crn.Model("retained_probability")
        model.add_species({"A": 0, "B": 0})
        model.add_parameter("a", 2.0)
        model.add_reaction(crn.mass_action("a"), [], ["A"])
        model.add_reaction(1.0, [], ["B"])
        initial = {"A": 1, "B": 1}
        q = crn.prepare_fsp(model, initial, retained=[(0, 0), (1, 0), (0, 1)])
        solved = forward_sensitivities(
            q, [(1.0, {"a": 2.0})], [1.0], ["a"], tolerance=2e-8
        )
        self.assertEqual(len(solved.states), 3)
        self.assertGreater(float(solved.probabilities[0, -1]), 0.5)
        self.assertAlmostEqual(float(solved.derivatives[0, 0].sum()), 0.0, places=9)
        with self.assertRaises(ValueError):
            forward_sensitivities(q, [(1.0, {"a": 2.0})], [2.0], ["a"])

    def test_checkpoint_replay_accounts_for_final_partial_population(self):
        with TemporaryDirectory() as root:
            dirs = [Path(root) / name for name in ("first", "second", "direct")]
            for path in dirs:
                path.mkdir()
            visited = []

            def loss(x):
                visited.append(x.copy())
                return np.sum(x * x)

            short = SimpleNamespace(threads=1, output=dirs[0])
            continued = SimpleNamespace(threads=1, output=dirs[1], resume_from=dirs[0])
            direct = SimpleNamespace(threads=1, output=dirs[2])
            search(
                short,
                "fixture",
                "de",
                loss,
                [(-1.0, 1.0)] * 2,
                12,
                11,
                pop=8,
                context="objective",
            )
            self.assertEqual(len(visited), 11)
            extended = search(
                continued,
                "fixture",
                "de",
                loss,
                [(-1.0, 1.0)] * 2,
                12,
                21,
                pop=8,
                context="objective",
            )
            self.assertEqual(len(visited), 21)
            fresh = search(
                direct,
                "fixture",
                "de",
                loss,
                [(-1.0, 1.0)] * 2,
                12,
                21,
                pop=8,
                context="objective",
            )
            np.testing.assert_array_equal(extended.candidates, fresh.candidates)
            np.testing.assert_array_equal(extended.scores, fresh.scores)
