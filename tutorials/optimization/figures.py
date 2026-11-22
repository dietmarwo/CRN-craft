#!/usr/bin/env python3
"""Draw the figures of the publication lessons C5 to C10 from saved tables.

Each lesson page shows figures from ``img/``. They are drawn from the
tables in ``results/C#/`` and from nothing else, so they can be redrawn
without running a study:

    python tutorials/optimization/figures.py            # all six lessons
    python tutorials/optimization/figures.py C7 C9      # two of them

``--tables DIR`` reads the tables of one lesson from another folder, for
example from a study run of your own, and ``--folder DIR`` writes the
figures somewhere else. Matplotlib is required.
"""

import argparse
from importlib import import_module
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
LESSONS = {
    "C5": "publications.rna_counts",
    "C6": "publications.stimulus_fit",
    "C7": "publications.reactor",
    "C8": "publications.experiment_design",
    "C9": "publications.stat5",
    "C10": "publications.pathway",
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("lessons", nargs="*", metavar="LESSON", help="C5 to C10")
    parser.add_argument("--tables", type=Path, help="tables of the one lesson named")
    parser.add_argument("--folder", type=Path, default=HERE / "img")
    args = parser.parse_args(argv)
    lessons = args.lessons or list(LESSONS)
    for lesson in lessons:
        if lesson not in LESSONS:
            parser.error(f"unknown lesson {lesson!r}; choose from {', '.join(LESSONS)}")
    if args.tables is not None and len(lessons) != 1:
        parser.error("--tables needs exactly one lesson")
    sys.path.insert(0, str(HERE))
    for lesson in lessons:
        tables = args.tables or HERE / "results" / lesson
        if not tables.is_dir():
            parser.error(f"no tables in {tables}")
        import_module(LESSONS[lesson]).draw(tables, args.folder)
        print(f"{lesson}: drew the figures in {args.folder} from {tables}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
