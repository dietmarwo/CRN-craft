#!/usr/bin/env python3
"""Dispatch one simulator-owned optimization tutorial."""

from importlib import import_module
import sys
from common import main

CHAPTERS = {
    "C1": "calibration",
    "C2": "sensitivity",
    "C3": "oscillators",
    "C4": "rules",
    "C5": "publications.rna_counts",
    "C6": "publications.stimulus_fit",
    "C7": "publications.reactor",
    "C8": "publications.experiment_design",
    "C9": "publications.stat5",
    "C10": "publications.pathway",
}

USAGE = (
    "usage: run.py {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10}"
    " --output NEW_DIRECTORY [--preset smoke|study]\n"
    "       run.py CHAPTER --help lists the options of one chapter"
)

if len(sys.argv) == 2 and sys.argv[1] in {"-h", "--help"}:
    print(USAGE)
    raise SystemExit(0)
if len(sys.argv) < 2 or sys.argv[1] not in CHAPTERS:
    unknown = f"unknown chapter {sys.argv[1]!r}\n" if len(sys.argv) > 1 else ""
    raise SystemExit(unknown + USAGE)
chapter = sys.argv.pop(1)
module = import_module(CHAPTERS[chapter])
raise SystemExit(
    main(
        chapter,
        module.run,
        getattr(module, "options", None),
        packages=getattr(module, "PACKAGES", None),
    )
)
