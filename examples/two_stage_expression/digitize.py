"""Convert marker centers from page 6 of the authors' arXiv PDF."""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# 1400-pixel-high PDF rendering, cropped to the upper 860 pixels.
# The plot spans y=499 to 318 for divergence 0 to 1.5.
EXPONENTS = (-1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2)
PIXELS = {
    (20.0, 2.5): (459.5, 466.5, 474.5, 483.5, 492.0, 497.5, 499.0, 499.0, 499.0),
    (0.5, 100.0): (328.0, 381.0, 432.5, 470.0, 488.0, 494.5, 497.0, 498.5, 498.5),
}


def digitize():
    with (ROOT / "reference-figure1d.csv").open("w", newline="") as f:
        out = csv.writer(f, lineterminator="\n")
        out.writerow(["a", "b", "gamma", "published_kl_display", "readout_error"])
        for (a, b), centers in PIXELS.items():
            for x, y in zip(EXPONENTS, centers, strict=True):
                out.writerow(
                    (a, b, 10**x, round((499 - y) * 1.5 / 181, 5), round(3 / 181, 5))
                )


if __name__ == "__main__":
    digitize()
