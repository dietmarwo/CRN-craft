"""Independent small references and negative acceptance controls (also under -O)."""

import copy
import unittest
import numpy as np
from calibration import binding_model, moments
from common import checked, require
from oscillators import generator_check, grid_niche, nondominated, score
from rules import exact_law


class ScienceChecks(unittest.TestCase):
    def test_require_survives_optimization(self):
        with self.assertRaises(ValueError):
            require(False, "negative control")

    def test_partial_ensemble_is_not_complete(self):
        r = dict(
            termination="completed",
            requested_replicates=4,
            completed_replicates=3,
            statistics={},
        )
        with self.assertRaises(ValueError):
            checked(r)
        r.update(completed_replicates=4, termination="event_limit")
        with self.assertRaises(ValueError):
            checked(r)

    def test_two_monomer_survival_probability(self):
        times = [0, 0.3, 1.0]
        law = exact_law(2, 1.7, times)
        # The sole pair waits Exp(1.7); at most one cluster iff it has merged.
        np.testing.assert_allclose(
            law[:, 0, 1], 1 - np.exp(-1.7 * np.asarray(times)), atol=1e-12
        )
        # Largest cluster remains <= 1 iff it has not merged.
        np.testing.assert_allclose(
            law[:, 1, 1], np.exp(-1.7 * np.asarray(times)), atol=1e-12
        )

    def test_wrong_kernel_has_a_different_law(self):
        a = exact_law(6, 1.2, [0.25, 0.7])
        b = exact_law(6, 1.2, [0.25, 0.7], True)
        self.assertGreater(float(np.max(np.abs(a - b))), 0.1)

    def test_binding_initial_and_conserved_moments(self):
        self.assertEqual(len(binding_model()["reactions"]), 2)
        mean, variance = moments("binding", [0.4, 0.2], [0, 1, 100])
        np.testing.assert_allclose(mean.sum(axis=1), 30)
        np.testing.assert_array_equal(mean[0], [30, 0])
        np.testing.assert_array_equal(variance[0], [0, 0])
        np.testing.assert_allclose(mean[-1], [10, 20])

    def test_hill_generator_and_archive_boundaries(self):
        self.assertLess(generator_check(), 1e-10)
        self.assertEqual(grid_niche([0, 0]), 0)
        self.assertEqual(grid_niche([80, 80]), 15)
        self.assertEqual(grid_niche([20, 0]), 1)
        self.assertEqual(grid_niche([0, 20]), 4)

    def test_constant_paths_are_not_oscillations(self):
        _, _, metrics = score(np.linspace(0, 80, 129), np.ones((3, 129, 3)))
        self.assertTrue(all(m["failed"] for m in metrics))

    def test_numerical_failure_cannot_enter_pareto_set(self):
        good = dict(valid=True, objectives=[1, 2, 3])
        bad = dict(valid=False, objectives=[0, 0, 0])
        dominated = copy.deepcopy(good)
        dominated["objectives"] = [2, 3, 4]
        self.assertEqual(nondominated([bad, dominated, good]), [good])


if __name__ == "__main__":
    unittest.main()
