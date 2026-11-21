"""Convert reviewed pixel anchors from the 1,000-pixel eLife Figure 3 render."""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# From elife-07464-fig3-v3 at 1000 px wide. Each tuple is
# (day, center y, outer top y, inner top y, inner bottom y, outer bottom y).
PIXELS = {
    "copies": [
        (10, 151, 141, 147, 154, 158),
        (20, 116, 103, 112, 123, 134),
        (35, 68, 41, 64, 88, 105),
        (50, 49, 24, 43, 69, 90),
        (70, 40, 15, 36, 47, 78),
        (100, 81, 24, 62, 136, 175),
    ],
    "normalized_variance": [
        (10, 329, 311, 326, 334, 342),
        (20, 304.5, 286, 297, 311, 323),
        (35, 289.5, 263, 281, 292, 313),
        (50, 274, 236, 261, 286, 303),
        (70, 256, 212, 252, 264, 287),
        (100, 327.5, 308, 322, 334, 343),
    ],
}


def scale(observable, day, y):
    if observable == "copies":
        return 10 ** (5 + (45 - y) / 61) if day < 100 else 1e5 * 2 ** ((105 - y) / 40)
    return (366 - y) * 0.02 / 47.5 if day < 100 else (380 - y) * 0.1 / 94


def digitize():
    with (ROOT / "reference-figure3.csv").open("w", newline="") as f:
        out = csv.writer(f, lineterminator="\n")
        out.writerow(
            [
                "day",
                "observable",
                "center",
                "outer_low",
                "inner_low",
                "inner_high",
                "outer_high",
                "center_pixel_low",
                "center_pixel_high",
            ]
        )
        for key, items in PIXELS.items():
            for day, center, top, top_i, bottom_i, bottom in items:
                values = [
                    scale(key, day, y)
                    for y in (
                        center,
                        bottom,
                        bottom_i,
                        top_i,
                        top,
                        center + 3,
                        center - 3,
                    )
                ]
                out.writerow([day, key, *[f"{v:.9g}" for v in values]])


if __name__ == "__main__":
    digitize()
