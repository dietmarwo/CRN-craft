"""Checks of the C5 sparse-data CDF score and its saved recovery tables."""

import csv
from pathlib import Path
import unittest

import numpy as np
from scipy.special import ndtr

from publications import rna_counts as r

RESULTS = Path(__file__).parent / "results" / "C5"

# Two times, four cells. Ties at count 0 (time 1) and count 2 (time 2).
DATA = np.array([[0, 2], [0, 1], [1, 2], [3, 2]])
# A predicted CDF on the support 0..3 at each time.
PREDICTED = np.array([[0.3, 0.6, 0.7, 0.9], [0.1, 0.2, 0.5, 0.8]])


def term(rank, value, cells=4):
    """One Equation 7 term with the Equation 8 weight cells / (F (1 - F))."""
    return cells * (value - rank) ** 2 / (rank * (1 - rank))


def rows(name):
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


class Score(unittest.TestCase):
    def test_published_ranks_by_hand(self):
        # Equation 9: each sorted cell l has rank (l - 0.5) / m; the last
        # cell is left out (L - 1 terms); ties keep separate ranks.
        expected = (
            term(1 / 8, 0.3)
            + term(3 / 8, 0.3)
            + term(5 / 8, 0.6)  # time 1: cells 0, 0, 1
            + term(1 / 8, 0.2)
            + term(3 / 8, 0.5)
            + term(5 / 8, 0.5)  # time 2: cells 1, 2, 2
        )
        self.assertAlmostEqual(r.CdfScore(DATA)(PREDICTED), expected, places=12)
        self.assertAlmostEqual(expected, 1.9657142857142857, places=12)

    def test_midranks_by_hand(self):
        # One term per distinct count at its midrank (below + half of ties) / m.
        expected = (
            term(1 / 4, 0.3)
            + term(5 / 8, 0.6)
            + term(7 / 8, 0.9)  # time 1: 0, 1, 3
            + term(1 / 8, 0.2)
            + term(5 / 8, 0.5)  # time 2: 1, 2
        )
        score = r.CdfScore(DATA, midrank=True)
        self.assertAlmostEqual(score(PREDICTED), expected, places=12)

    def test_population_ranks_by_hand(self):
        expected = term(1 / 2, 0.3) + term(3 / 4, 0.6) + term(1 / 4, 0.2)
        score = r.CdfScore(DATA, population=True)
        self.assertAlmostEqual(score(PREDICTED), expected, places=12)

    def test_counts_route_matches_the_cdf_route(self):
        # The empirical midrank CDF of the data itself scores zero.
        score = r.CdfScore(DATA, midrank=True)
        self.assertAlmostEqual(score.from_counts(DATA), 0.0, places=12)

    def test_midrank_score_is_consistent_with_tied_counts(self):
        # A large sample from a known integer law: the true law must score
        # lower than a law with less mass at zero. The individually ranked
        # reading prefers the wrong law, the tie bias found in review.
        rng = np.random.default_rng(11)
        law = np.array([0.5, 0.3, 0.15, 0.05])
        sample = rng.choice(4, size=(20_000, 1), p=law)
        shifted = np.array([0.42, 0.36, 0.16, 0.06])

        def mid(p):
            return (np.cumsum(p) - 0.5 * p)[None, :]

        def upper(p):
            return np.cumsum(p)[None, :]

        midrank = r.CdfScore(sample, midrank=True)
        self.assertLess(midrank(mid(law)), midrank(mid(shifted)))
        published = r.CdfScore(sample)
        self.assertGreater(published(upper(law)), published(upper(shifted)))


class Readout(unittest.TestCase):
    def test_rounded_noise_on_a_silent_promoter(self):
        # With k1 = 0 the RNA count stays zero, so an observation is the
        # rounded normal error: P(obs <= x) = Phi((x + 0.5) / 0.5).
        import crn_craft as crn

        times = np.array([1.0, 2.0])
        graph = crn.prepare_fsp(r.model(), {"silent": 1, "active": 1, "rna": 20})
        data = np.array([[-1, 0], [0, 0], [0, 1], [1, 2]])
        for midrank in (False, True):
            score = r.CdfScore(data, midrank=midrank)
            x = np.arange(-1, 3)
            cdf = ndtr((x + 0.5) / 0.5)  # N(0, 0.25) read as a variance
            if midrank:
                cdf = 0.5 * (cdf + ndtr((x - 0.5) / 0.5))
            expected = score(np.array([cdf, cdf]))
            got = r.objective(graph, score, times)((0.0, 1.0, 1.0))
            self.assertAlmostEqual(got, expected, delta=1e-7 * expected)


@unittest.skipUnless(RESULTS.is_dir(), "saved C5 tables are absent")
class SavedTables(unittest.TestCase):
    def test_delivered_medians_and_search_shortfall(self):
        table = rows("recovery.csv")
        for variant, median in (
            ("midrank", (0.030791, 0.56681, 0.85696)),
            ("published", (0.13053, 5.0, 2.2138)),
        ):
            fits = [
                row
                for row in table
                if row["run"] == "delivered-300"
                and row["score"] != ""
                and row["sparse_score"] == variant
                and row["cells"] == "100"
            ]
            self.assertEqual(len(fits), 16)
            best = {}
            for row in fits:
                key = row["data_index"]
                if key not in best or float(row["score"]) < float(best[key]["score"]):
                    best[key] = row
            rates = np.median(
                [[float(row[k]) for k in r.NAMES] for row in best.values()], axis=0
            )
            np.testing.assert_allclose(rates, median, rtol=5e-5)
        short = [
            row
            for row in table
            if row["run"] == "delivered-300"
            and row["sparse_score"] == "midrank"
            and row["cells"] == "100"
            and float(row["score"]) > float(row["truth_score"])
        ]
        self.assertEqual(len(short), 11)

    def test_reached_flags_match_scores(self):
        for row in rows("recovery.csv"):
            reached = float(row["score"]) <= float(row["truth_score"])
            self.assertEqual(row["reached_truth_score"] == "true", reached, row)

    def test_review_fits_reach_the_truth_score(self):
        review = [row for row in rows("recovery.csv") if row["run"] == "review-3000"]
        self.assertGreaterEqual(len(review), 4)
        self.assertTrue(all(row["reached_truth_score"] == "true" for row in review))

    def test_contour_grid(self):
        table = rows("contour-m10.csv")
        self.assertEqual(len(table), 31 * 31)
        k2 = sorted({float(row["k2"]) for row in table})
        np.testing.assert_allclose(k2, np.linspace(0.1, 1.0, 31))


if __name__ == "__main__":
    unittest.main()
