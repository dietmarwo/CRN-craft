"""The PEtab v1 fields needed for a single-model linear Gaussian fit."""

from collections import defaultdict, deque
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path

import numpy as np
import yaml

from common import require


class UniqueYaml(yaml.SafeLoader):
    def construct_mapping(self, node, deep=False):
        pairs = self.construct_pairs(node, deep)
        result = dict(pairs)
        require(len(result) == len(pairs), "duplicate PEtab YAML key")
        return result


def table(path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(
            reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
            f"duplicate or absent TSV header in {path.name}",
        )
        rows = list(reader)
        require(
            rows and all(None not in r and None not in r.values() for r in rows),
            f"invalid rows in {path.name}",
        )
        return rows


def unique(rows, field):
    result = {row[field]: row for row in rows}
    require(len(result) == len(rows), f"duplicate {field}")
    return result


def gaussian_nll(value, prediction, sigma):
    sigma = np.asarray(sigma, dtype=float)
    if np.any(sigma <= 0) or not np.isfinite(sigma).all():
        raise ValueError("Gaussian standard deviations must be positive and finite")
    return (
        0.5 * ((np.asarray(value) - prediction) / sigma) ** 2
        + np.log(sigma)
        + 0.5 * np.log(2 * np.pi)
    )


def key(row):
    return (
        row["observableId"],
        row["datasetId"],
        row["simulationConditionId"],
        float(row["time"]),
    )


@dataclass
class Problem:
    model: Path
    parameters: dict
    fits: tuple
    bounds: np.ndarray
    observations: dict
    measurements: list
    condition: str
    times: np.ndarray
    y: np.ndarray
    t_index: np.ndarray
    o_index: np.ndarray
    noise: tuple

    @property
    def nominal(self):
        return {
            name: float(row["nominalValue"]) for name, row in self.parameters.items()
        }

    def bind(self, x):
        result = self.nominal
        for name, scaled in zip(self.fits, x, strict=True):
            result[name] = (
                10**scaled
                if self.parameters[name]["parameterScale"] == "log10"
                else scaled
            )
        return result

    def score(self, readout, parameters):
        return gaussian_nll(
            self.y,
            readout[self.t_index, self.o_index],
            [parameters[name] for name in self.noise],
        )

    def score_archive(self, path):
        indexed = defaultdict(deque)
        for row in table(path):
            indexed[key(row)].append(
                (float(row["simulation"]), float(row["noiseParameters"]))
            )
        predictions = []
        for row in self.measurements:
            bucket = indexed[key(row)]
            require(bool(bucket), "missing reference simulation row")
            predictions.append(bucket.popleft())
        require(not any(indexed.values()), "unused reference simulation rows")
        predicted, sigma = np.asarray(predictions).T
        return gaussian_nll(self.y, predicted, sigma)


def read(problem_path):
    base = problem_path.resolve().parent
    config = yaml.load(problem_path.read_text(encoding="utf-8"), Loader=UniqueYaml)
    require(
        config.get("format_version") == 1 and len(config.get("problems", [])) == 1,
        "expected one PEtab v1 problem",
    )
    require(
        set(config) == {"format_version", "parameter_file", "problems"},
        "unsupported PEtab problem fields",
    )
    sub = config["problems"][0]
    keys = {"condition_files", "measurement_files", "observable_files", "sbml_files"}
    require(
        keys <= set(sub) <= keys | {"visualization_files"},
        "unsupported PEtab file groups",
    )
    require(
        all(isinstance(sub[k], list) and len(sub[k]) == 1 for k in keys),
        "expected one file per PEtab table and model",
    )

    def local(name):
        path = (base / name).resolve()
        require(
            path.is_relative_to(base) and not Path(name).is_absolute(),
            "PEtab paths must stay in the input folder",
        )
        require(
            path.is_file() and path.stat().st_size <= 4 << 20,
            f"invalid PEtab input {name}",
        )
        return path

    params = unique(table(local(config["parameter_file"])), "parameterId")
    bounds, free = [], []
    for name, row in params.items():
        a, b, x = map(
            float, [row["lowerBound"], row["upperBound"], row["nominalValue"]]
        )
        scale = row["parameterScale"]
        require(
            all(map(isfinite, (a, b, x)))
            and a <= x <= b
            and a < b
            and scale in ("lin", "log10")
            and (a > 0 or scale == "lin"),
            f"invalid parameter domain for {name}",
        )
        require(
            row["estimate"] in ("0", "1")
            and not any(
                value for field, value in row.items() if "prior" in field.lower()
            ),
            f"unsupported estimation flag or prior for {name}",
        )
        if row["estimate"] == "1":
            free.append(name)
            bounds.append(np.log10([a, b]) if scale == "log10" else [a, b])
    conditions = table(local(sub["condition_files"][0]))
    require(
        len(conditions) == 1 and set(conditions[0]) <= {"conditionId", "conditionName"},
        "this adapter supports one condition without parameter overrides",
    )
    condition = conditions[0]["conditionId"]
    observations = unique(table(local(sub["observable_files"][0])), "observableId")
    for name, row in observations.items():
        require(
            row["observableTransformation"] == "lin"
            and row["noiseDistribution"] == "normal"
            and row["noiseFormula"] == f"noiseParameter1_{name}",
            f"expected a linear Gaussian observation with a named scale for {name}",
        )
    rows = table(local(sub["measurement_files"][0]))
    for row in rows:
        require(
            not row.get("preequilibrationConditionId")
            and not row.get("observableParameters")
            and row["simulationConditionId"] == condition
            and row["observableId"] in observations
            and row["noiseParameters"] in params
            and isfinite(float(row["measurement"]))
            and isfinite(float(row["time"]))
            and float(row["time"]) >= 0.0,
            "unsupported measurement row",
        )
    times = np.unique([0.0, *[float(row["time"]) for row in rows]])
    ids = list(observations)
    return Problem(
        local(sub["sbml_files"][0]),
        params,
        tuple(free),
        np.asarray(bounds),
        observations,
        rows,
        condition,
        times,
        np.array([float(r["measurement"]) for r in rows]),
        np.searchsorted(times, [float(r["time"]) for r in rows]),
        np.array([ids.index(r["observableId"]) for r in rows]),
        tuple(r["noiseParameters"] for r in rows),
    )
