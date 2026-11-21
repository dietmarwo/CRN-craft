"""Compile compact count histograms from the pinned public research archive."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent
REV = json.loads((ROOT / "manifest.json").read_text())["upstream"]
UPSTREAM = (
    f"https://raw.githubusercontent.com/MunskyGroup/Munsky_et_al_2018/{REV}/Codes"
)
KEYS = (
    "k12",
    "k21",
    "k21_signal",
    "k23",
    "k32",
    "k34",
    "k43",
    "i1",
    "i2",
    "i3",
    "i4",
    "decay",
    "delay",
    "export",
    "nuclear_decay",
)
FILES = {}
SITE_DOSES = (0.4, 0.4, 0.4, 0.2, 0.2)


def readmat(path, local):
    if local:
        raw = (local / path).read_bytes()
    else:
        with urlopen(f"{UPSTREAM}/{path}", timeout=45) as response:
            raw = response.read()
    FILES[path] = hashlib.sha256(raw).hexdigest()
    return loadmat(io.BytesIO(raw), simplify_cells=True)


def write(name, header, rows):
    with (ROOT / name).open("w", newline="") as file:
        out = csv.writer(file, lineterminator="\n")
        out.writerow(header.split(","))
        out.writerows(rows)


def mature_counts(experiments):
    histogram, audit = [], []
    for dose, experiment in zip((0.2, 0.4), experiments, strict=True):
        for rep, values in zip(
            experiment["rep"], experiment["gene"][0]["rep"], strict=True
        ):
            identifier = rep["Name"].split("date")[1]
            data = values["DAPIth"][1]["type"]
            total, nuclear, cytoplasm = (
                values["RNArawtot"],
                data[1]["RNAraw"],
                data[2]["RNAraw"],
            )
            for i, minute in enumerate(rep["tt"]):
                observed = np.column_stack(
                    [x[:, i] for x in (total, nuclear, cytoplasm)]
                )
                observed = observed[np.isfinite(observed[:, 0])]
                if (
                    len(observed) != rep["NoCells"][i]
                    or not np.all(observed >= 0)
                    or np.any(observed != observed.astype(int))
                ):
                    raise ValueError(f"Inconsistent counts for {identifier}, {minute}")
                rows, repeats = np.unique(
                    observed.astype(int), axis=0, return_counts=True
                )
                histogram.extend(
                    (dose, identifier, int(minute), *counts, n)
                    for counts, n in zip(rows, repeats, strict=True)
                )
                audit.append(
                    (
                        dose,
                        identifier,
                        int(minute),
                        len(observed),
                        int(
                            np.count_nonzero(
                                observed[:, 0] != observed[:, 1] + observed[:, 2]
                            )
                        ),
                    )
                )
    return histogram, audit


def figure_traces(node):
    if node["type"] in ("line", "graph2d.lineseries"):
        p = node["properties"]
        if "XData" in p and "YData" in p and np.size(p["XData"]) > 3:
            yield (
                np.asarray(p["XData"]),
                np.asarray(p["YData"]),
                p.get("LineStyle", "-"),
            )
    children = node["children"]
    if isinstance(children, dict):
        children = [children]
    for item in children:
        yield from figure_traces(item)


def published_sites(local):
    activity, distribution = [], []
    for dose, tag in ((0.2, "0p2"), (0.4, "0p4")):
        prefix = f"Figures/1_TS_Figures/STL1_{tag}NaCl_Spatial_TS"
        a = readmat(prefix + "_Fraction_vs_t.fig", local)["hgS_070000"]
        e = readmat(prefix + "_Ave_Distributions.fig", local)["hgS_070000"]
        for i, (x, y, style) in enumerate(figure_traces(a), 1):
            kind = f"measurement-{i}" if style == "none" else "published-model"
            activity.extend(
                (dose, kind, float(t), float(p)) for t, p in zip(x, y, strict=True)
            )
        for i, (x, y, _) in enumerate(figure_traces(e), 1):
            staircase = np.any(np.diff(x) == 0)
            kind = f"measurement-{i}" if staircase else "published-model"
            for k, first in zip(*np.unique(x, return_index=True), strict=True):
                if staircase:
                    # At a staircase edge, both the previous and the next
                    # heights are serialized. Keep the right-hand height.
                    first = np.flatnonzero(x == k)[-1]
                distribution.append((dose, kind, int(k), float(y[first])))
    write("reference-site-activity.csv", "dose,series,minute,active_fraction", activity)
    write(
        "reference-site-distribution.csv",
        "dose,series,upper_intensity,probability",
        distribution,
    )


def setup(local):
    data = readmat("Hog_Data/Analyzed_FISH_Data_2016_06.mat", local)
    histogram, audit = mature_counts(data["RNAFISHexp"])
    write(
        "source-counts.csv",
        "dose,replicate,minute,total,nuclear,cytoplasm,cells",
        histogram,
    )
    write("source-audit.csv", "dose,replicate,minute,cells,total_disagrees", audit)
    published_sites(local)
    table = readmat("Hog_Data/RNA_TS_Brian.mat", local)["RNA_TS"]
    ts = []
    for irow, (dose, info) in enumerate(zip(SITE_DOSES, table, strict=True), 1):
        rep = f"site-{irow}"
        for i, minute in enumerate(info["tt"]):
            raw = np.array(info["TS"][:, i], dtype=float)
            # Archive is NaN-padded. The denominator includes the cells
            # whose brightest-spot intensity was unavailable.
            n = int(info["Nc"][i])
            detected = raw[np.isfinite(raw)]
            active = detected[detected > 2]
            values, counts = np.unique(np.floor(active).astype(int), return_counts=True)
            ts.append((dose, rep, int(minute), 0, n - len(active), n))
            ts.extend(
                (dose, rep, int(minute), int(k), int(c), n)
                for k, c in zip(values, counts, strict=True)
            )
    write(
        "source-sites.csv", "dose,replicate,minute,intensity_bin,cells,total_cells", ts
    )
    anchors = {}
    for space in ("NonSpatial", "Spatial"):
        for fit in ("Means", "Distributions"):
            path = f"Fit_Results/STL1_0p2and0p4/{space}/{fit}/4State_k21only_noloop/Best_Pars_AC.mat"
            v = readmat(path, local)["BEST_PARS"]
            if space == "Spatial":
                sites = [v[34], v[38], v[36], v[32]]
            else:
                sites = [v[32], v[38], 0.0, 0.0]
            k = [v[i] for i in (0, 3, 15, 4, 7, 8, 11, 24, 25, 26, 27)]
            k[2] = -k[2]
            anchors[f"{space}-{fit}"] = dict(zip(KEYS, [*k, *sites], strict=True))
    (ROOT / "source-parameters.json").write_text(
        json.dumps({"parameters": anchors, "sha256": FILES}, indent=2) + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--local", type=Path, help="Path to the upstream Codes directory")
    setup(p.parse_args().local)
