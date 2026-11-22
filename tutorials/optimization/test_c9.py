"""Tests of the C9 PEtab reader, the archive score and the STAT5 model mapping."""

import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import yaml

from publications import petab

DATA = Path(os.environ.get("CRN_PUBLICATION_DATA", ".publication-data")) / "C9"
ODE = Path(
    os.environ.get(
        "CRN_DIFFSOL",
        Path(__file__).parents[2] / "integrations/diffsol/target/release/crn-diffsol",
    )
)

YAML = """format_version: 1
parameter_file: parameters.tsv
problems:
- condition_files: [conditions.tsv]
  measurement_files: [measurements.tsv]
  observable_files: [observables.tsv]
  sbml_files: [model.xml]
"""
TABLES = {
    "parameters.tsv": [
        "parameterId\tparameterScale\tlowerBound\tupperBound\tnominalValue\testimate",
        "k\tlog10\t1E-3\t1000\t0.5\t1",
        "sd_x\tlog10\t1E-3\t1000\t2\t1",
        "c\tlin\t0\t5\t1\t0",
    ],
    "conditions.tsv": ["conditionId\tconditionName", "c1\tone"],
    "observables.tsv": [
        "observableId\tobservableFormula\tnoiseFormula"
        "\tobservableTransformation\tnoiseDistribution",
        "x\t2 * A + c\tnoiseParameter1_x\tlin\tnormal",
    ],
    "measurements.tsv": [
        "observableId\tsimulationConditionId\tmeasurement\ttime\tnoiseParameters"
        "\tdatasetId",
        "x\tc1\t1.0\t0\tsd_x\td",
        "x\tc1\t3.0\t2.5\tsd_x\td",
        "x\tc1\t4.0\t2.5\tsd_x\td",
    ],
    "simulations.tsv": [
        "observableId\tsimulationConditionId\tsimulation\ttime\tnoiseParameters"
        "\tdatasetId",
        "x\tc1\t2.0\t2.50\t2\td",
        "x\tc1\t1.5\t0.0\t2\td",
        "x\tc1\t2.5\t2.5\t2\td",
    ],
}


def gauss(y, mean, sigma):
    return (
        0.5 * ((y - mean) / sigma) ** 2 + math.log(sigma) + 0.5 * math.log(2 * math.pi)
    )


class Fixture(unittest.TestCase):
    def setUp(self):
        self.folder = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.folder)
        (self.folder / "problem.yaml").write_text(YAML)
        (self.folder / "model.xml").write_text("<sbml/>")
        for name, lines in TABLES.items():
            (self.folder / name).write_text("\n".join(lines) + "\n")

    def edit(self, name, old, new):
        path = self.folder / name
        content = path.read_text()
        self.assertIn(old, content)
        path.write_text(content.replace(old, new, 1))

    def refuses(self, pattern):
        with self.assertRaisesRegex(ValueError, pattern):
            petab.read(self.folder / "problem.yaml")


class ReaderAccepts(Fixture):
    def test_reads_and_binds(self):
        problem = petab.read(self.folder / "problem.yaml")
        self.assertEqual(problem.fits, ("k", "sd_x"))
        np.testing.assert_allclose(problem.bounds, [[-3, 3], [-3, 3]])
        np.testing.assert_allclose(problem.times, [0.0, 2.5])
        bound = problem.bind([1.0, 0.0])
        self.assertEqual((bound["k"], bound["sd_x"], bound["c"]), (10.0, 1.0, 1.0))
        np.testing.assert_allclose(problem.scaled(bound), [1.0, 0.0])

    def test_archive_joins_repeats_in_order(self):
        problem = petab.read(self.folder / "problem.yaml")
        costs = problem.score_archive(self.folder / "simulations.tsv")
        expected = [gauss(1.0, 1.5, 2), gauss(3.0, 2.0, 2), gauss(4.0, 2.5, 2)]
        np.testing.assert_allclose(costs, expected, rtol=0, atol=1e-14)

    def test_archive_missing_and_unused_rows(self):
        problem = petab.read(self.folder / "problem.yaml")
        self.edit("simulations.tsv", "x\tc1\t2.5\t2.5\t2\td\n", "")
        with self.assertRaisesRegex(ValueError, "no archived simulation row"):
            problem.score_archive(self.folder / "simulations.tsv")
        self.edit(
            "simulations.tsv", "x\tc1\t1.5\t0.0", "x\tc1\t1.5\t0.0\t2\td\nx\tc1\t1.5\t9"
        )
        with self.assertRaisesRegex(ValueError, "no archived|unused"):
            problem.score_archive(self.folder / "simulations.tsv")

    def test_score_uses_row_noise_and_normalization(self):
        problem = petab.read(self.folder / "problem.yaml")
        readout = np.array([[1.5], [2.0]])
        for sigma in (0.5, 2.0, 7.0):
            p = problem.bind([0.0, math.log10(sigma)])
            expected = sum(
                gauss(y, m, sigma) for y, m in [(1.0, 1.5), (3.0, 2.0), (4.0, 2.0)]
            )
            self.assertAlmostEqual(float(problem.score(readout, p).sum()), expected, 12)
        for bad in (0.0, -1.0, math.nan, math.inf):
            with self.assertRaises(ValueError):
                petab.gaussian_nll([1.0], [1.0], [bad])


class ReaderRefuses(Fixture):
    def test_unsafe_yaml(self):
        # An unsafe loader would build the integer 1 and accept the problem.
        (self.folder / "problem.yaml").write_text(
            YAML.replace(
                "format_version: 1",
                "format_version: !!python/object/apply:builtins.int ['1']",
            )
        )
        with self.assertRaises(yaml.constructor.ConstructorError):
            petab.read(self.folder / "problem.yaml")
        (self.folder / "problem.yaml").write_text(YAML + "format_version: 1\n")
        self.refuses("duplicate PEtab YAML key")
        (self.folder / "problem.yaml").write_text(
            YAML.replace("[conditions.tsv]", "&a [conditions.tsv]").replace(
                "[observables.tsv]", "*a"
            )
        )
        self.refuses("anchors and aliases")

    def test_yaml_fields(self):
        self.edit("problem.yaml", "format_version: 1", "format_version: 2")
        self.refuses("one PEtab v1 problem")
        (self.folder / "problem.yaml").write_text(YAML + "extra: 1\n")
        self.refuses("unsupported PEtab problem fields")
        (self.folder / "problem.yaml").write_text(
            YAML.replace("[model.xml]", "[model.xml, model.xml]")
        )
        self.refuses("one file per PEtab table")

    def test_paths_stay_inside(self):
        outside = self.folder.parent / f"{self.folder.name}-outside.tsv"
        shutil.copy(self.folder / "parameters.tsv", outside)
        self.addCleanup(outside.unlink)
        for name in (
            str(outside),
            f"../{outside.name}",
            "https://example.org/p.tsv",
            "sub\\p.tsv",
        ):
            (self.folder / "problem.yaml").write_text(
                YAML.replace("parameters.tsv", name)
            )
            self.refuses("relative to the problem folder")
        (self.folder / "link.tsv").symlink_to(outside)
        (self.folder / "problem.yaml").write_text(
            YAML.replace("parameters.tsv", "link.tsv")
        )
        self.refuses("symbolic link")

    def test_size_and_encoding(self):
        with (self.folder / "measurements.tsv").open("a") as stream:
            stream.write("#" * (petab.LIMIT + 1))
        self.refuses("exceeds")
        (self.folder / "measurements.tsv").write_bytes(
            "\n".join(TABLES["measurements.tsv"]).encode("utf-16")
        )
        self.refuses("not UTF-8")

    def test_columns_and_fields(self):
        self.edit("parameters.tsv", "\testimate\n", "\n")
        self.refuses("lacks columns")

    def test_extra_field(self):
        self.edit("measurements.tsv", "sd_x\td\n", "sd_x\td\textra\n")
        self.refuses("fields, expected")

    def test_missing_field(self):
        self.edit("measurements.tsv", "sd_x\td\n", "sd_x\n")
        self.refuses("fields, expected")

    def test_unknown_column(self):
        self.edit("measurements.tsv", "\tdatasetId", "\tdatasetIdx")
        self.refuses("unsupported columns")

    def test_condition_override(self):
        self.edit("conditions.tsv", "conditionName", "conditionName\tk")
        self.edit("conditions.tsv", "one", "one\t2")
        self.refuses("parameter overrides")

    def test_parameter_values(self):
        cases = [
            ("0.5\t1", "abc\t1", "not a number"),
            ("0.5\t1", "nan\t1", "not finite"),
            ("1000\t0.5", "1000\t5000", "outside its bounds"),
            ("1E-3\t1000\t0.5", "1000\t1E-3\t0.5", "not increasing"),
            ("k\tlog10", "k\tln", "parameterScale 'ln'"),
            ("0.5\t1", "0.5\tyes", "estimate of k"),
            ("c\tlin\t0", "c\tlog10\t0", "needs bounds above 0"),
        ]
        for old, new, message in cases:
            with self.subTest(message):
                self.setUp()
                self.edit("parameters.tsv", old, new)
                self.refuses(message)

    def test_priors_and_duplicates(self):
        self.edit("parameters.tsv", "estimate\n", "estimate\tobjectivePriorType\n")
        for row in ("0.5\t1\n", "2\t1\n", "1\t0\n"):
            self.edit("parameters.tsv", row, row[:-1] + "\tnormal\n")
        self.refuses("priors are not supported")
        self.setUp()
        self.edit("parameters.tsv", "c\tlin", "k\tlin")
        self.refuses("duplicate parameterId")

    def test_observables(self):
        cases = [
            ("\tlin\tnormal", "\tlog\tnormal", "observableTransformation"),
            ("\tlin\tnormal", "\tlin\tlaplace", "noiseDistribution"),
            ("noiseParameter1_x", "noiseParameter1_x * 2", "noiseFormula"),
            ("2 * A + c", "__import__('os')", "must be arithmetic"),
        ]
        for old, new, message in cases:
            with self.subTest(message):
                self.setUp()
                self.edit("observables.tsv", old, new)
                self.refuses(message)

    def test_measurements(self):
        cases = [
            ("\t1.0\t0\t", "\tnan\t0\t", "not finite"),
            ("\t1.0\t0\t", "\t1.0\t-1\t", "negative time"),
            ("\t1.0\t0\t", "\t1.0\tinf\t", "not finite"),
            ("1.0\t0\tsd_x", "1.0\t0\t0.5", "noiseParameters must name"),
            ("1.0\t0\tsd_x", "1.0\t0\tc", "positive lower bound"),
            ("x\tc1\t1.0", "y\tc1\t1.0", "unknown observable"),
            ("x\tc1\t1.0", "x\tc2\t1.0", "unknown simulation condition"),
        ]
        for old, new, message in cases:
            with self.subTest(message):
                self.setUp()
                self.edit("measurements.tsv", old, new)
                self.refuses(message)
        for column in ("observableParameters", "preequilibrationConditionId"):
            with self.subTest(column):
                self.setUp()
                lines = TABLES["measurements.tsv"]
                (self.folder / "measurements.tsv").write_text(
                    "\n".join(
                        [lines[0] + f"\t{column}", lines[1] + "\tc1"]
                        + [line + "\t" for line in lines[2:]]
                    )
                    + "\n"
                )
                self.refuses("not supported")


@unittest.skipUnless(DATA.joinpath("problem.yaml").is_file(), "C9 fixtures absent")
class Archive(unittest.TestCase):
    def test_archived_nll(self):
        problem = petab.read(DATA / "problem.yaml")
        costs = problem.score_archive(DATA / "simulations.tsv")
        self.assertEqual(len(costs), 48)
        self.assertEqual(len(problem.fits), 9)
        self.assertAlmostEqual(float(costs.sum()), 138.2219997062, delta=1e-9)
        per_signal = [float(costs[problem.o_index == i].sum()) for i in range(3)]
        np.testing.assert_allclose(
            per_signal, [44.2806969, 52.8702405, 41.0710623], rtol=0, atol=5e-8
        )


@unittest.skipUnless(
    DATA.joinpath("problem.yaml").is_file() and ODE.is_file(),
    "C9 fixtures or crn-diffsol binary absent",
)
class Model(unittest.TestCase):
    def test_hand_mapping_matches_sbml_at_three_vectors(self):
        from publications import stat5
        from publications.concentrations import Solver

        problem = petab.read(DATA / "problem.yaml")
        with tempfile.TemporaryDirectory() as folder:
            with Solver(ODE, Path(folder) / "m.json", stat5.model(problem)) as solver:
                for name, parameters in stat5.vectors(problem).items():
                    with self.subTest(name):
                        _, reference = stat5.reference(problem, parameters)
                        _, native = solver.batch(
                            [[parameters[k] for k in stat5.KINETIC]],
                            problem.times,
                            rtol=2e-10,
                            atol=1e-12,
                        )
                        self.assertLess(np.max(np.abs(native[0] - reference)), 1e-5)
                        gap = (
                            problem.score(native[0], parameters).sum()
                            - problem.score(reference, parameters).sum()
                        )
                        self.assertLess(abs(gap), 1e-4)


RESULTS = Path(__file__).parent / "results" / "C9"


class SavedResults(unittest.TestCase):
    def test_tables_agree_and_draw(self):
        from publications import stat5

        rows = {r["source"]: r for r in stat5.read_csv(RESULTS / "nll.csv")}
        self.assertAlmostEqual(float(rows["archive"]["total"]), 138.2219997062, 9)
        for row in rows.values():
            parts = sum(
                float(row[k]) for k in ("pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel")
            )
            self.assertAlmostEqual(parts, float(row["total"]), 8)
        for row in stat5.read_csv(RESULTS / "verification.csv"):
            self.assertLess(float(row["signal_max_abs"]), 0.01)
            self.assertLess(
                abs(float(row["nll_native"]) - float(row["nll_reference"])), 0.05
            )
        with tempfile.TemporaryDirectory() as folder:
            stat5.draw(RESULTS, Path(folder))
            for name in ("C9-signals.png", "C9-fits.png"):
                self.assertGreater((Path(folder) / name).stat().st_size, 10_000)


if __name__ == "__main__":
    unittest.main()
