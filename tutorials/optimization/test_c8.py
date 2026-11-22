"""Checks of C8: fixture verification, Fisher information, sensitivities."""

from hashlib import sha256
import io
import json
import os
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

import crn_craft as crn
from crn_craft.sensitivities import fisher_information, forward_sensitivities
import numpy as np

from publications import sources
from publications.experiment_design import A, B, analytic_check

FILES = {
    "keep/a.txt": b"1 2 3\n",
    "keep/b.txt": b"4.5\n",
    "ref.txt": b"7\n",
    "other/script.py": b"raise SystemExit\n",
}


def digest(files):
    lines = "".join(
        f"{sha256(files[name]).hexdigest()}  {name}\n" for name in sorted(files)
    )
    return sha256(lines.encode()).hexdigest()


SPEC = {
    "kind": "tar-content",
    "prefixes": ["keep/"],
    "files": ["ref.txt"],
    "members": 3,
    "sha256": digest({k: v for k, v in FILES.items() if k != "other/script.py"}),
}


def archive(path, files=FILES, prefix="repo-abc", extra=(), order=None, mtime=0):
    """A gzip tar with a top folder, like GitHub tarballs and git archive."""
    with tarfile.open(path, "w:gz", format=tarfile.PAX_FORMAT) as tar:
        top = tarfile.TarInfo(prefix)
        top.type = tarfile.DIRTYPE
        tar.addfile(top)
        for name in order or files:
            info = tarfile.TarInfo(f"{prefix}/{name}")
            info.size = len(files[name])
            info.mtime = mtime
            tar.addfile(info, io.BytesIO(files[name]))
        for info, data in extra:
            tar.addfile(info, None if data is None else io.BytesIO(data))
    return path


class FixtureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def read(self, path, spec=SPEC):
        seen = {}
        sources.read_tar(path, spec, lambda name, data: seen.update({name: data}))
        return seen

    def test_content_digest_ignores_the_container(self):
        one = archive(self.dir / "github.tar.gz")
        # Other top folder, member order and timestamps: a git archive.
        two = archive(
            self.dir / "git.tar.gz",
            prefix="Repo-0123",
            order=list(reversed(FILES)),
            mtime=1_700_000_000,
        )
        self.assertNotEqual(one.read_bytes(), two.read_bytes())
        self.assertEqual(self.read(one), self.read(two))
        self.assertEqual(set(self.read(one)), {"keep/a.txt", "keep/b.txt", "ref.txt"})

    def test_modified_missing_and_extra_members_are_refused(self):
        changed = dict(FILES, **{"keep/a.txt": b"1 2 4\n"})
        with self.assertRaisesRegex(ValueError, "differs"):
            self.read(archive(self.dir / "changed.tar.gz", changed))
        missing = {k: v for k, v in FILES.items() if k != "ref.txt"}
        with self.assertRaisesRegex(ValueError, "differs"):
            self.read(archive(self.dir / "missing.tar.gz", missing))
        extra = dict(FILES, **{"keep/c.txt": b"0\n"})
        with self.assertRaisesRegex(ValueError, "differs"):
            self.read(archive(self.dir / "extra.tar.gz", extra))

    def test_unsafe_members_are_refused(self):
        def entry(name, kind=tarfile.REGTYPE, link=""):
            info = tarfile.TarInfo(name)
            info.type = kind
            info.linkname = link
            return info, (b"" if kind == tarfile.REGTYPE else None)

        cases = {
            "member type": entry("repo-abc/keep/l", tarfile.SYMTYPE, "/etc/passwd"),
            "member type ": entry("repo-abc/x", tarfile.LNKTYPE, "repo-abc/ref.txt"),
            "member type  ": entry("repo-abc/dev", tarfile.CHRTYPE),
            "unsafe": entry("repo-abc/../../evil.txt"),
            "unsafe ": entry("/abs/evil.txt"),
            "duplicate": entry("repo-abc/ref.txt"),
        }
        for message, item in cases.items():
            path = archive(self.dir / "bad.tar.gz", extra=[item])
            with (
                self.subTest(message),
                self.assertRaisesRegex(ValueError, message.strip()),
            ):
                self.read(path)

    def test_sizes_are_bounded(self):
        path = archive(self.dir / "big.tar.gz")
        with mock.patch.object(sources, "MAX_MEMBER", 3):
            with self.assertRaisesRegex(ValueError, "too large"):
                self.read(path)
        with mock.patch.object(sources, "MAX_UNPACKED", 8):
            with self.assertRaisesRegex(ValueError, "too large"):
                self.read(path)
        with mock.patch.object(sources, "MAX_ARCHIVE", 10):
            with self.assertRaisesRegex(ValueError, "too large"):
                self.read(path)

    def test_download_urls_are_restricted(self):
        for url in (
            "http://raw.githubusercontent.com/a/b/c",
            "file:///etc/passwd",
            "https://example.org/a.tar.gz",
            "https://raw.githubusercontent.com:8443/a",
        ):
            with self.subTest(url), self.assertRaisesRegex(ValueError, "refused"):
                sources._download(url, 10)
        handler = sources._Redirects()
        with self.assertRaisesRegex(ValueError, "refused"):
            handler.redirect_request(None, None, 302, "", {}, "https://evil.test/x")

    def test_manifest_pins_every_fixture(self):
        manifest = json.loads(sources.MANIFEST.read_text())
        for chapter, entries in manifest.items():
            for name, entry in entries.items():
                with self.subTest(chapter=chapter, name=name):
                    if isinstance(entry, dict):
                        self.assertEqual(entry["kind"], "tar-content")
                        self.assertRegex(entry["sha256"], "^[0-9a-f]{64}$")
                        url = entry["url"]
                    else:
                        url, pin = entry
                        self.assertRegex(pin, "^[0-9a-f]{64}$")
                    sources._check_url(url)
                    self.assertIn("2e2f76f" if chapter == "C8" else "", url)

    def test_a_byte_pin_still_refuses_another_file(self):
        path = self.dir / "C9" / "x.tsv"
        path.parent.mkdir()
        path.write_bytes(b"x")
        with self.assertRaisesRegex(ValueError, "checksum"):
            sources._verify(path, ["https://raw.githubusercontent.com/x", "0" * 64])


class FisherTests(unittest.TestCase):
    def test_formula_symmetry_and_semidefiniteness(self):
        rng = np.random.default_rng(1)
        p = rng.dirichlet(np.ones(6), size=3)
        p[1, 2] = 0.0  # unreachable bin: no information, no mass
        g = rng.normal(size=(3, 4, 6))
        g -= g.mean(axis=2, keepdims=True)
        g[1, :, 2] = 0.0
        cells = np.array([100.0, 50.0, 10.0])
        info, excluded = fisher_information(p, g, cells)
        expected = sum(
            cells[t] * g[t][:, b, None] * g[t][None, :, b] / p[t, b]
            for t in range(3)
            for b in range(6)
            if p[t, b] > 0
        )
        np.testing.assert_allclose(info, expected, rtol=1e-12)
        np.testing.assert_array_equal(info, info.T)
        self.assertGreaterEqual(np.linalg.eigvalsh(info)[0], -1e-9 * np.abs(info).max())
        np.testing.assert_array_equal(excluded, 0.0)
        floored, excluded = fisher_information(p, g, cells, probability_floor=0.1)
        self.assertGreater(excluded.max(), 0.0)
        self.assertTrue(np.all(np.linalg.eigvalsh(info - floored) > -1e-9))

    def test_invalid_input_is_refused(self):
        p = np.full((1, 2), 0.5)
        g = np.zeros((1, 1, 2))
        for args in (
            (p, np.zeros((1, 1, 3)), 1.0),
            (np.array([[0.5, np.nan]]), g, 1.0),
            (np.array([[1.1, -0.1]]), g, 1.0),
            (p, g, -1.0),
            (p, g, [1.0, 2.0]),
        ):
            with self.assertRaises(ValueError):
                fisher_information(*args)


class SensitivityTests(unittest.TestCase):
    def test_birth_death_against_poisson_and_differences(self):
        report = {"checks": []}
        analytic_check(report)
        row = report["checks"][0]
        self.assertTrue(row["passed"], row)
        self.assertLess(row["analytic_gap"], 1e-7)
        self.assertLess(row["mass_error"], 1e-12)
        self.assertLess(row["derivative_mass_error"], 1e-12)

    def test_default_initial_law_is_the_prepared_one(self):
        model = crn.Model("mixture")
        model.add_species("x", 0)
        model.add_parameter("d", 1.0)
        model.add_reaction(crn.mass_action("d"), ["x"], [])
        law = {(2,): 0.25, (3,): 0.75}
        graph = crn.prepare_fsp(model, {"x": 3}, initial=law)
        result = forward_sensitivities(graph, [(1.0, {"d": 1.0})], [0.0, 1.0], ["d"])
        start = dict(
            zip(map(tuple, graph.states), result.probabilities[0, :-1], strict=True)
        )
        self.assertEqual(start[(2,)], 0.25)
        self.assertEqual(start[(3,)], 0.75)
        mean = sum(
            x * p
            for (x,), p in zip(graph.states, result.probabilities[1, :-1], strict=True)
        )
        self.assertAlmostEqual(mean, 2.75 * np.exp(-1.0), places=6)

    def test_hill_rate_through_a_switch_matches_differences(self):
        model = crn.Model("hill_switch")
        model.add_species({"x": 0, "y": 0})
        model.add_observable("cy", "y")
        base = dict(b=0.5, k=6.0, a=0.3, c=2.0, d=1.0)
        for name, value in base.items():
            model.add_parameter(name, value)
        model.add_reaction("b+k/(1+a*cy^2)", [], ["x"])
        model.add_reaction(crn.mass_action("c"), [], ["y"])
        model.add_reaction(crn.mass_action("d"), ["x"], [])
        model.add_reaction(crn.mass_action("d"), ["y"], [])
        graph = crn.prepare_fsp(model, {"x": 30, "y": 25})
        segments = [(1.0, base), (2.0, base | {"c": 6.0})]
        result = forward_sensitivities(
            graph, segments, [2.0], ["a", "c"], tolerance=1e-10
        )
        self.assertLess(abs(result.derivatives.sum(axis=2)).max(), 1e-9)

        def law(name, delta):
            initial = None
            for (end, binding), start in zip(segments, (0.0, 1.0), strict=True):
                binding = binding | {name: binding[name] + delta}
                solved = graph.solve(
                    end - start,
                    [end - start],
                    initial=initial,
                    parameters=binding,
                    tolerance=1e-12,
                )
                p = np.maximum(solved.snapshots[-1].probabilities, 0.0)
                initial = {
                    tuple(map(int, s)): q
                    for s, q in zip(graph.states, p / p.sum(), strict=True)
                    if q > 0
                }
            return p

        for j, name in enumerate(("a", "c")):
            h = 1e-3 * base[name]
            difference = (law(name, h) - law(name, -h)) / (2 * h)
            gap = np.abs(result.derivatives[0, j, :-1] - difference).max()
            self.assertLess(gap, 1e-6, name)

    def test_retained_sets_are_validated(self):
        model = crn.Model("pair")
        model.add_species({"A": 0, "B": 0})
        model.add_reaction(1.0, [], ["A"])
        for retained in ([], [(0, -1)], [(0, 0, 0)], [(3, 0)], [(1, 0)]):
            with self.subTest(retained=retained), self.assertRaises(ValueError):
                crn.prepare_fsp(model, {"A": 2, "B": 2}, retained=retained)
        cut = crn.prepare_fsp(model, {"A": 2, "B": 2}, retained=[(0, 0), (1, 0)])
        self.assertEqual(sorted(map(tuple, cut.states)), [(0, 0), (1, 0)])


DATA = Path(os.environ.get("CRN_PUBLICATION_DATA", ".publication-data"))


@unittest.skipUnless(
    (DATA / "C8" / "fsp-fim.tar.gz").is_file(),
    "C8 archive absent: run publications/sources.py C8 or set CRN_PUBLICATION_DATA",
)
class ArchiveTests(unittest.TestCase):
    """Needs the 65 MB author archive; reads it once (about 10 s)."""

    @classmethod
    def setUpClass(cls):
        from publications.experiment_design import archive

        spec = json.loads(sources.MANIFEST.read_text())["C8"]["fsp-fim.tar.gz"]
        cls.designs, cls.matrices, cls.references = archive(
            DATA / "C8" / "fsp-fim.tar.gz", spec
        )

    def test_archive_ranking(self):
        self.assertEqual(len(self.designs), 672)
        self.assertEqual(max(self.designs, key=lambda k: self.designs[k].mean()), A)
        self.assertEqual(int(np.sum(self.designs[A] > self.designs[B])), 97)

    def test_native_matrix_reproduces_the_archive(self):
        from publications.experiment_design import ToggleDesign

        toggle = ToggleDesign(self.references)
        native, _ = toggle.matrix(A, 1)
        archived = self.matrices[A][1]
        gap = np.linalg.norm(native - archived) / np.linalg.norm(archived)
        self.assertLess(gap, 1e-4)


if __name__ == "__main__":
    unittest.main()
