"""Checks of the C6 BayFish reader, generator and information criteria."""

import csv
import os
from pathlib import Path
import unittest

import numpy as np

from publications import stimulus_fit as s

RESULTS = Path(__file__).parent / "results" / "C6"
DATA = Path(os.environ.get("CRN_PUBLICATION_DATA", ".publication-data"))
FIXTURES = DATA / "C6"


def rows(name):
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


class Reader(unittest.TestCase):
    def test_two_state_site_rule_and_encodings(self):
        text = "0.00\t0.00\t3\n3.56\t0.00\t25\n0.00\t2.10\t4\n1.5\t9.25\t120\n\n"
        expected = [[0, 3], [1, 25], [1, 4], [2, 120]]
        for raw in (
            text.encode(),
            text.encode("utf-8-sig"),
            text.replace("\n", "\r\n").encode("utf-16"),  # with a BOM
        ):
            np.testing.assert_array_equal(s.parse_counts(raw, "t"), expected)

    def test_rejects_malformed_rows(self):
        for bad in (
            "0 0\n",  # two fields
            "0 0 1 2\n",  # four fields
            "a 0 1\n",
            "0 0 1.5\n",  # RNA count must be an integer
            "nan 0 1\n",
            "0 nan 1\n",
            "0 inf 1\n",
            "-1 0 1\n",
            "0 0 -1\n",
            "0 0 20000\n",
            "\n\n",  # no cells
        ):
            with self.assertRaises(ValueError, msg=bad):
                s.parse_counts(bad.encode(), "t")

    def test_rejects_unmarked_utf16_and_large_files(self):
        with self.assertRaises(ValueError):
            s.parse_counts("0 0 1\n".encode("utf-16-le"), "t")
        with self.assertRaises(ValueError):
            s.parse_counts(b"0 0 1\n" * (s.MAX_BYTES // 6 + 1), "t")

    def test_counts_above_the_cutoff_are_refused(self):
        groups = [np.array([[0, 5]]), np.array([[1, 11]])]
        with self.assertRaises(ValueError):
            s.StimulusFit(groups, 10)


class Generator(unittest.TestCase):
    def test_columns_conserve_probability(self):
        rates = (0.0093, 0.1108, 0.0078, 8.14)
        for reflect in (True, False):
            q = s.cme(rates, 40, reflect=reflect).toarray()
            np.testing.assert_allclose(q.sum(axis=0), 0.0, atol=1e-12)
        # Without reflection the top RNA row leaks into the sink state.
        q = s.cme(rates, 40, reflect=False).toarray()
        self.assertAlmostEqual(q[-1, 40], 2 * 0.0078)

    def test_stationary_law(self):
        rates = (0.0093, 0.1108, 0.0078, 8.14)
        p = s.normalized_law(rates, 60)
        q = s.cme(rates, 60, reflect=True)[:-1, :-1]
        self.assertLess(np.abs(q @ p).sum(), 1e-12)
        sites = p.reshape(3, 61).sum(axis=1)
        f = rates[0] / (rates[0] + rates[1])
        np.testing.assert_allclose(sites, [(1 - f) ** 2, 2 * f * (1 - f), f**2])


class Criteria(unittest.TestCase):
    def test_information_criteria(self):
        aic, bic = s.information(2098.0, 6, 376)
        self.assertAlmostEqual(aic, 4208.0)
        self.assertAlmostEqual(bic, 4196.0 + 6 * np.log(376))
        self.assertAlmostEqual(s.corrected_aic(2098.0, 6, 376), 4208.0 + 84 / 369)

    def test_nested_pairs(self):
        models = [
            dict(changed="kon", changed_indices=(0,), nll=10.0),
            dict(changed="koff", changed_indices=(1,), nll=12.0),
            dict(changed="kon-koff", changed_indices=(0, 1), nll=11.0),
        ]
        pairs = {(p["superset"], p["subset"]): p for p in s.nested_pairs(models)}
        self.assertEqual(set(pairs), {("kon-koff", "kon"), ("kon-koff", "koff")})
        self.assertFalse(pairs["kon-koff", "kon"]["passed"])
        self.assertTrue(pairs["kon-koff", "koff"]["passed"])


@unittest.skipUnless(RESULTS.is_dir(), "saved C6 tables are absent")
class SavedTables(unittest.TestCase):
    def test_table_follows_from_the_saved_nlls(self):
        table = [row for row in rows("models.csv") if row["run"] == "delivered"]
        self.assertEqual(len(table), 7)
        criteria = {
            row["changed"]: s.information(float(row["nll"]), int(row["k"]), 376)
            for row in table
        }
        best_aic = min(a for a, _ in criteria.values())
        best_bic = min(b for _, b in criteria.values())
        for row in table:
            aic, bic = criteria[row["changed"]]
            self.assertEqual(int(row["k"]), 4 + len(row["changed"].split("-")))
            self.assertAlmostEqual(float(row["delta_aic"]), aic - best_aic, places=3)
            self.assertAlmostEqual(float(row["delta_bic"]), bic - best_bic, places=3)
        self.assertEqual(min(criteria, key=lambda k: criteria[k][1]), "kon-koff")

    def test_saved_nested_pairs_hold(self):
        table = [
            dict(
                changed=row["changed"],
                changed_indices=tuple(
                    s.NAMES.index(n) for n in row["changed"].split("-")
                ),
                nll=float(row["nll"]),
            )
            for row in rows("models.csv")
            if row["run"] == "delivered"
        ]
        pairs = s.nested_pairs(table)
        self.assertEqual(len(pairs), 12)
        self.assertTrue(all(p["passed"] for p in pairs))


@unittest.skipUnless(
    FIXTURES.is_dir(), "BayFish fixtures absent (CRN_PUBLICATION_DATA)"
)
class Fixtures(unittest.TestCase):
    def test_cells_and_fixed_vector(self):
        from publications.sources import fixtures

        groups = s.read_counts(fixtures(DATA, "C6"))
        self.assertEqual(tuple(map(len, groups)), s.TOTALS)
        self.assertEqual(
            [np.bincount(g[:, 0], minlength=3).tolist() for g in groups],
            [[105, 23, 9], [53, 44, 38], [6, 11, 21], [6, 17, 43]],
        )
        case = s.StimulusFit(groups, 300)
        native, *_ = case.probabilities(s.TABLE_BASAL, s.TABLE_STIMULATED)
        reference = s.reference(s.TABLE_BASAL, s.TABLE_STIMULATED, 300)
        total = case.score(native)
        nll = total.cells * total.nats_per_cell
        self.assertAlmostEqual(nll, 2098.742654619, places=5)
        check = case.score(reference)
        self.assertAlmostEqual(check.cells * check.nats_per_cell, nll, delta=1e-5)


if __name__ == "__main__":
    unittest.main()
