"""Prepare parameter inputs and read published online variance summaries."""

import argparse
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
REV = json.loads((ROOT / "manifest.json").read_text())["upstream"]
BASE = f"https://raw.githubusercontent.com/ImperialCollegeLondon/MitoNetworksGenetics/{REV}/Linear_Feedback_Control/Stochastic_analysis/"
TABLES = {
    "E": ("mu_sweep/param_sweep_vals.csv", "mu_sweep@param_sweep_vals.csv", "mu_sweep"),
    "F": (
        "kappa_sweep/param_sweep_vals.csv",
        "kappa_sweep@param_sweep_vals.csv",
        "kappa_sweep",
    ),
    "G": (
        "delta_sweep/Make_sweep_vals/delta_sw_params.csv",
        "delta_sweep@params.csv",
        "delta_sweep/Post_processing",
    ),
    "H": (
        "network_and_b_sweep/param_sweep_vals.csv",
        "network_and_b_sweep@param_sweep_vals.csv",
        "network_and_b_sweep/Post_processing",
    ),
}
PARAMETERS = (
    "ws_init",
    "ms_init",
    "wf_init",
    "mf_init",
    "mu",
    "beta",
    "gamma",
    "kappa",
    "b",
    "delta",
    "net_mag",
    "net_rat",
)


def get(path):
    with urlopen(BASE + path, timeout=30) as response:
        return response.read().decode()


def prepare(local=None):
    rows = []
    for panel, (path, filename, _) in TABLES.items():
        source = (local / filename).read_text() if local else get(path)
        for index, data in enumerate(csv.DictReader(io.StringIO(source))):
            if not data["mu"]:
                continue
            if panel == "H" and (
                index >= 121
                or (float(data["net_mag"]) != 1 and float(data["net_rat"]) != 1)
            ):
                continue
            rows.append(
                {
                    "panel": panel,
                    "index": index,
                    **{k: float(data.get(k, 1)) for k in PARAMETERS},
                }
            )
    with (ROOT / "conditions.csv").open("w", newline="") as f:
        out = csv.DictWriter(f, ["panel", "index", *PARAMETERS], lineterminator="\n")
        out.writeheader()
        out.writerows(rows)
    return rows


def summaries(row):
    panel, index = row["panel"], row["index"]
    try:
        source = get(f"{TABLES[panel][2]}/Data/online_stats_{index}.csv")
    except HTTPError as e:
        if e.code != 404:
            raise
        return []
    lines = list(csv.DictReader(io.StringIO(source)))
    # Keep the author's column names in a log so an accidental schema
    # change cannot silently change the variance observation.
    keys = {k.strip().lower(): k for k in lines[0]}
    time = next((keys[k] for k in ("t", "time", "t_store") if k in keys), None)
    variance = next(
        (
            keys[k]
            for k in ("var_h", "v_h", "h_var", "varh", "h_v", "v(h)")
            if k in keys
        ),
        None,
    )
    if time is None or variance is None:
        raise ValueError(f"Unknown online summary columns: {list(keys)}")
    return [(panel, index, float(a[time]), float(a[variance])) for a in lines]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local", type=Path)
    parser.add_argument("--online-summaries", action="store_true")
    args = parser.parse_args()
    rows = prepare(args.local)
    if args.online_summaries:
        with (
            ThreadPoolExecutor(max_workers=8) as pool,
            (ROOT / "reference-variance.csv").open("w", newline="") as f,
        ):
            writer = csv.writer(f, lineterminator="\n")
            writer.writerow(["panel", "index", "day", "variance"])
            for result in pool.map(summaries, rows):
                writer.writerows(result)


if __name__ == "__main__":
    main()
