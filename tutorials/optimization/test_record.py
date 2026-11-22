"""The record of the publication smoke runs refuses incomplete evidence."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import record


def summary(**change):
    value = dict(
        state="complete",
        passed=True,
        preset="smoke",
        versions={"crn-craft": "0.9.0", "numpy": "2.3.5"},
        checks=[dict(name="gate", passed=True, detail=1.5)],
    )
    return value | change


def write(root, **changes):
    for name in record.PUBLICATIONS:
        folder = Path(root) / name
        folder.mkdir(exist_ok=True)
        (folder / "summary.json").write_text(
            json.dumps(summary(**changes.get(name, {})))
        )


class PublicationRecord(unittest.TestCase):
    def test_complete_smoke_runs_keep_check_names_and_verdicts_only(self):
        with TemporaryDirectory() as root:
            write(root, C9=dict(versions={"crn-craft": "0.9.0", "PyYAML": "6.0.4"}))
            made = record.publications(Path(root))
        self.assertEqual(list(made["chapters"]), record.PUBLICATIONS)
        self.assertEqual(made["chapters"]["C5"], [dict(name="gate", passed=True)])
        self.assertEqual(
            made["versions"],
            {"crn-craft": "0.9.0", "numpy": "2.3.5", "PyYAML": "6.0.4"},
        )

    def test_incomplete_failed_or_mixed_evidence_is_refused(self):
        cases = {
            "failed": dict(C7=dict(passed=False)),
            "running": dict(C6=dict(state="running")),
            "study": dict(C8=dict(preset="study")),
            "no checks": dict(C5=dict(checks=[])),
            "versions": dict(C10=dict(versions={"crn-craft": "0.8.0"})),
        }
        for label, changes in cases.items():
            with self.subTest(label), TemporaryDirectory() as root:
                write(root, **changes)
                with self.assertRaises(ValueError):
                    record.publications(Path(root))
        with TemporaryDirectory() as root:
            write(root)
            (Path(root) / "C8" / "summary.json").unlink()
            with self.assertRaises(ValueError):
                record.publications(Path(root))


if __name__ == "__main__":
    unittest.main()
