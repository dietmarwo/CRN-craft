"""Convert the two small numerical exhibits in the authors' archive."""

import csv
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
REVISION = "45462390b4c27cf2ffe4860c4f526e67d49ad6c8"
URL = f"https://raw.githubusercontent.com/ferdinando17/Survival-of-the-Densest/{REVISION}/Simulated_data/Main_text"

PANELS = {
    "deterministic": "Upper_left",
    "isolated": "Upper_right",
    "equal": "Bottom_left",
    "dense": "Bottom_right",
}


def read_panel(url, local):
    if local is not None:
        return local.read_text().splitlines()
    with urlopen(url, timeout=30) as response:
        return response.read().decode("utf-8-sig").splitlines()


def convert(local=None):
    with (ROOT / "reference-1c.csv").open("w", newline="") as file:
        out = csv.writer(file, lineterminator="\n")
        out.writerow(["condition", "degradation", "day", "mean", "sem_upper"])
        for condition, panel in PANELS.items():
            for degradation in ("neutral", "extra"):
                kind = (
                    "preferential_elimination"
                    if degradation == "extra"
                    else "deterministic"
                    if condition == "deterministic"
                    else "neutral"
                )
                name = f"ssd-{panel}-{kind}.dat"
                lines = read_panel(
                    f"{URL}/Figure1/1C/{panel}/{kind}.dat",
                    local / name if local else None,
                )
                for line in lines:
                    if not line.strip() or line.startswith("#"):
                        continue
                    data = [float(x) for x in line.replace(",", " ").split()]
                    # The average of marginal SEMs upper-bounds the SEM
                    # of the paired-region mean for any covariance.
                    error = (
                        sum(data[3::2]) / 2
                        if len(data) == 6
                        else data[2]
                        if len(data) == 3
                        else 0
                    )
                    out.writerow([condition, degradation, data[0], data[1], error])
    with (ROOT / "reference-2e.csv").open("w", newline="") as file:
        out = csv.writer(file, lineterminator="\n")
        out.writerow(["day", "position_um", "mean"])
        for day in (0, 1480):
            lines = read_panel(
                f"{URL}/Figure2/2E/Day{day}.csv",
                local / f"ssd-day{day}.csv" if local else None,
            )
            x, h = list(csv.reader(lines))
            out.writerows(
                (day, position, fraction)
                for position, fraction in zip(x[1:], h[1:], strict=True)
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local", type=Path)
    convert(parser.parse_args().local)
