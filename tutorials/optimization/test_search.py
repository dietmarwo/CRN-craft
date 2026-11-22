"""The shared search driver: accounting, failures, checkpoint replay, refusals.

Every check runs on a cheap analytic objective and needs only NumPy, SciPy,
loguru and fcmaes; the module is skipped with the missing package named.
"""

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

try:
    from publications.search import PENALTY, search
except ImportError as error:  # pragma: no cover - depends on the environment
    raise unittest.SkipTest(
        f"search driver needs fcmaes, scipy and loguru: {error}"
    ) from None

BOUNDS = [(-1.0, 1.0)] * 3
SEED = 12
POP = 8
# Keys of run.json that depend on the host, the budget or the interruption.
VOLATILE = {"elapsed_seconds", "requested_evaluations", "replayed", "workers"}


def sphere(x):
    return float(np.sum(np.square(x - 0.25)))


def pareto(x):
    return np.array([np.sum(np.square(x)), np.sum(np.square(x - 1.0))])


def make_args(root, name, threads=1, resume=None):
    output = Path(root) / name
    output.mkdir()
    return SimpleNamespace(threads=threads, output=output, resume_from=resume)


def run(root, name, algorithm, budget, function=sphere, resume=None, threads=1, **kw):
    """Run one arm; returns the result, the arguments and the evaluation count."""
    calls = []

    def counted(x):
        calls.append(1)
        return function(x)

    args = make_args(root, name, threads, resume)
    kw.setdefault("pop", POP)
    kw.setdefault("observables", 2 if algorithm == "mode" else 1)
    if algorithm == "mode":
        counted = lambda x: pareto(x)  # noqa: E731
    found = search(args, "arm", algorithm, counted, BOUNDS, SEED, budget, **kw)
    return found, args, len(calls)


def manifest(args):
    meta = json.loads((Path(args.output) / "arm.run.json").read_text())
    return {k: v for k, v in meta.items() if k not in VOLATILE}


def trace(args):
    return (Path(args.output) / "arm.csv").read_bytes()


class SearchReplay(unittest.TestCase):
    def check_same(self, left, right):
        np.testing.assert_array_equal(left[0].candidates, right[0].candidates)
        np.testing.assert_array_equal(left[0].scores, right[0].scores)
        np.testing.assert_array_equal(left[0].best, right[0].best)
        self.assertEqual(left[0].value, right[0].value)
        self.assertEqual(left[0].failures, right[0].failures)
        self.assertEqual(trace(left[1]), trace(right[1]))
        self.assertEqual(manifest(left[1]), manifest(right[1]))

    def test_interrupted_search_equals_uninterrupted_for_every_optimizer(self):
        # 19 = two full populations and a partial one of three; 16 ends on a
        # population boundary; 11 then 19 then 37 chains two continuations.
        for algorithm in ("de", "cmaes", "mode", "random"):
            for first in (19, 16):
                with self.subTest(algorithm=algorithm, first=first):
                    with TemporaryDirectory() as root:
                        direct = run(root, "direct", algorithm, 37)
                        one = run(root, "one", algorithm, first)
                        resumed = run(
                            root, "resumed", algorithm, 37, resume=one[1].output
                        )
                        self.check_same(resumed, direct)
                        self.assertEqual(direct[2], 37 if algorithm != "mode" else 0)
                        if algorithm != "mode":
                            self.assertEqual(resumed[2], 37 - first)
                        meta = json.loads(
                            (resumed[1].output / "arm.run.json").read_text()
                        )
                        self.assertEqual(meta["replayed"], first)
            with self.subTest(algorithm=algorithm, chain=True):
                with TemporaryDirectory() as root:
                    direct = run(root, "direct", algorithm, 37)
                    a = run(root, "a", algorithm, 11)
                    b = run(root, "b", algorithm, 19, resume=a[1].output)
                    c = run(root, "c", algorithm, 37, resume=b[1].output)
                    self.check_same(c, direct)

    def test_threads_do_not_change_the_proposals(self):
        for algorithm in ("de", "cmaes"):
            with self.subTest(algorithm=algorithm), TemporaryDirectory() as root:
                one = run(root, "one", algorithm, 29)
                three = run(root, "three", algorithm, 29, threads=3)
                np.testing.assert_array_equal(one[0].scores, three[0].scores)
                np.testing.assert_array_equal(one[0].candidates, three[0].candidates)

    def test_budget_is_exact_and_partial_batch_is_only_scored(self):
        for algorithm in ("de", "cmaes", "random"):
            for budget in (8, 9, 11, 16, 17):
                with self.subTest(algorithm=algorithm, budget=budget):
                    with TemporaryDirectory() as root:
                        found, args, calls = run(root, "x", algorithm, budget)
                        self.assertEqual(calls, budget)
                        self.assertEqual(len(found.scores), budget)
                        rows = list(
                            csv.reader(
                                (args.output / "arm.csv").read_text().splitlines()
                            )
                        )
                        self.assertEqual(len(rows) - 1, budget)
                        meta = json.loads((args.output / "arm.run.json").read_text())
                        self.assertEqual(meta["actual_evaluations"], budget)
                        self.assertEqual(meta["requested_evaluations"], budget)
        with TemporaryDirectory() as root:
            # A budget below the population shrinks the population to it.
            found, args, calls = run(root, "x", "de", 5, pop=20)
            self.assertEqual((calls, len(found.scores)), (5, 5))
            with self.assertRaises(ValueError):
                run(root, "y", "de", 3)

    def test_failures_count_as_evaluations_and_replay_identically(self):
        def hostile(x):
            kind = int(abs(x[0]) * 1000) % 5
            if kind == 0:
                raise ValueError("domain")
            if kind == 1:
                return float("nan")
            if kind == 2:
                return float("inf")
            if kind == 3:
                raise ArithmeticError("overflow")
            return sphere(x)

        for algorithm in ("de", "cmaes", "random"):
            with self.subTest(algorithm=algorithm), TemporaryDirectory() as root:
                direct = run(root, "direct", algorithm, 37, function=hostile)
                part = run(root, "part", algorithm, 19, function=hostile)
                resumed = run(
                    root,
                    "resumed",
                    algorithm,
                    37,
                    function=hostile,
                    resume=part[1].output,
                )
                self.check_same(resumed, direct)
                found = direct[0]
                self.assertGreater(found.failures, 0)
                self.assertEqual(len(found.scores), 37)
                self.assertEqual(
                    int(np.sum(found.scores[:, 0] == PENALTY)), found.failures
                )
                self.assertTrue(np.isfinite(found.scores).all())
                self.assertLess(found.value, PENALTY)

    def test_other_exceptions_are_not_hidden(self):
        def broken(x):
            raise KeyError("bug")

        with TemporaryDirectory() as root, self.assertRaises(KeyError):
            run(root, "x", "de", 16, function=broken)

    def test_replay_trusts_saved_scores_only_for_the_same_context(self):
        with TemporaryDirectory() as root:
            first = run(root, "first", "de", 16, context="objective A")
            run(root, "same", "de", 24, resume=first[1].output, context="objective A")
            with self.assertRaises(ValueError):
                run(
                    root,
                    "other",
                    "de",
                    24,
                    resume=first[1].output,
                    context="objective B",
                )


class ResumeRefusals(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.root = self.temporary.name
        self.saved = run(self.root, "saved", "de", 19, context="c", sigma=0.3)[1]
        self.count = 0

    def tearDown(self):
        self.temporary.cleanup()

    def attempt(self, resume=None, **change):
        self.count += 1
        settings = dict(algorithm="de", budget=27, context="c", sigma=0.3)
        settings.update(change)
        return run(
            self.root,
            f"try{self.count}",
            settings.pop("algorithm"),
            settings.pop("budget"),
            resume=resume or self.saved.output,
            **settings,
        )

    def copy_saved(self):
        self.count += 1
        target = Path(self.root) / f"copy{self.count}"
        target.mkdir()
        for path in self.saved.output.iterdir():
            (target / path.name).write_bytes(path.read_bytes())
        return target

    def test_the_unchanged_configuration_continues(self):
        self.assertEqual(self.attempt()[2], 27 - 19)

    def test_another_configuration_is_refused(self):
        for change in (
            dict(algorithm="random"),
            dict(algorithm="cmaes"),
            dict(pop=10),
            dict(context="d"),
            dict(sigma=0.2),
            dict(observables=2),
        ):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.attempt(**change)

    def test_another_seed_or_bounds_is_refused(self):
        args = make_args(self.root, "seed", resume=self.saved.output)
        with self.assertRaises(ValueError):
            search(
                args, "arm", "de", sphere, BOUNDS, SEED + 1, 27, pop=POP, context="c"
            )
        args = make_args(self.root, "box", resume=self.saved.output)
        with self.assertRaises(ValueError):
            search(
                args, "arm", "de", sphere, [(-2, 1)] * 3, SEED, 27, pop=POP, context="c"
            )

    def test_missing_saved_output_is_refused(self):
        with self.assertRaises(ValueError):
            self.attempt(resume=Path(self.root) / "absent")

    def test_incomplete_saved_run_is_refused(self):
        copy = self.copy_saved()
        (copy / "arm.run.json").unlink()
        with self.assertRaises(ValueError):
            self.attempt(resume=copy)

    def test_truncated_trace_is_refused(self):
        lines = (self.saved.output / "arm.csv").read_bytes().splitlines(keepends=True)
        for name, content in {
            "fewer rows": b"".join(lines[:-3]),
            "cut row": b"".join(lines[:-1]) + lines[-1][:25],
            "header only": lines[0],
        }.items():
            copy = self.copy_saved()
            (copy / "arm.csv").write_bytes(content)
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.attempt(resume=copy)

    def test_edited_or_foreign_trace_is_refused(self):
        copy = self.copy_saved()
        rows = list(csv.reader((copy / "arm.csv").read_text().splitlines()))
        rows[5][1] = repr(float(rows[5][1]) + 1e-9)
        with (copy / "arm.csv").open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)
        with self.assertRaises(ValueError):
            self.attempt(resume=copy)
        other = self.copy_saved()
        meta = json.loads((other / "arm.run.json").read_text())
        meta["status"] = "running"
        (other / "arm.run.json").write_text(json.dumps(meta))
        with self.assertRaises(ValueError):
            self.attempt(resume=other)

    def test_columns_of_another_search_are_refused(self):
        copy = self.copy_saved()
        text = (copy / "arm.csv").read_text().replace("x2", "y2", 1)
        (copy / "arm.csv").write_text(text)
        with self.assertRaises(ValueError):
            self.attempt(resume=copy)


class RunRecord(unittest.TestCase):
    def test_run_json_records_backend_options_and_no_hidden_state(self):
        for algorithm, backend in (("de", "fcmaes"), ("random", "numpy")):
            with self.subTest(algorithm=algorithm), TemporaryDirectory() as root:
                _, args, _ = run(root, "x", algorithm, 17, context="ctx", sigma=0.4)
                meta = json.loads((args.output / "arm.run.json").read_text())
                self.assertEqual(meta["backend"], backend)
                self.assertRegex(meta["backend_version"], r"^\d+\.\d+")
                for key, value in dict(
                    algorithm=algorithm,
                    population=POP,
                    seed=SEED,
                    sigma=0.4,
                    context="ctx",
                    failures=0,
                    status="completed",
                    observables=1,
                    bounds=[[-1.0, 1.0]] * 3,
                ).items():
                    self.assertEqual(meta[key], value, key)

    def test_trace_holds_no_timing(self):
        with TemporaryDirectory() as root:
            _, args, _ = run(root, "x", "de", 17)
            header = (args.output / "arm.csv").read_text().splitlines()[0]
            self.assertEqual(header, "call,x0,x1,x2,f0,error")

    def test_unknown_optimizer_is_refused(self):
        with TemporaryDirectory() as root, self.assertRaises(ValueError):
            run(root, "x", "annealing", 16)


if __name__ == "__main__":
    unittest.main()
