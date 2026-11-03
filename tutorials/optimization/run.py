#!/usr/bin/env python3
"""Dispatch one simulator-owned optimization tutorial."""

import sys
from common import main

if len(sys.argv) < 2 or sys.argv[1] not in {"C1", "C2", "C3", "C4"}:
    raise SystemExit(
        "usage: run.py {C1,C2,C3,C4} --output NEW_DIRECTORY [--preset smoke|study]"
    )
chapter = sys.argv.pop(1)
module = __import__(
    {"C1": "calibration", "C2": "sensitivity", "C3": "oscillators", "C4": "rules"}[
        chapter
    ]
)
raise SystemExit(main(chapter, module.run))
