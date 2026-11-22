"""Read the subset of PEtab v1 needed for a single-model linear Gaussian fit.

The files are untrusted input. Every file is at most `LIMIT` bytes of UTF-8
text inside the folder of `problem.yaml` (no absolute paths, `..`, URLs or
symbolic links). The YAML file is read with a safe loader that also refuses
duplicate keys, anchors and aliases. The TSV tables must have exactly the
supported columns, one value per column, finite numbers and unique ids.
Formulas are never evaluated in Python: the observation formulas go to the
native expression parser after a character and length check, and the noise
formula must be the single placeholder `noiseParameter1_<observableId>`.

Supported: one condition without parameter overrides, one SBML file,
parameters on `lin` or `log10` scale without priors, `estimate` 0 or 1,
linear observations with normal noise, one named noise parameter per
measurement row, nonnegative finite times and repeated measurements.
Everything else is refused with a message that names the feature.
"""

from collections import defaultdict, deque
import csv
from dataclasses import dataclass
import io
from math import isfinite
from pathlib import Path
import re

import numpy as np
import yaml

from common import require

LIMIT = 4 << 20
ROWS = 100_000
FORMULA = 2_000
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
ARITHMETIC = re.compile(r"[A-Za-z0-9_ .+\-*/^()]*\Z")
PRIORS = (
    "initializationPriorType",
    "initializationPriorParameters",
    "objectivePriorType",
    "objectivePriorParameters",
)
COLUMNS = {
    "parameters": (
        {
            "parameterId",
            "parameterScale",
            "lowerBound",
            "upperBound",
            "nominalValue",
            "estimate",
        },
        {"parameterName", *PRIORS},
    ),
    "conditions": ({"conditionId"}, {"conditionName"}),
    "observables": (
        {"observableId", "observableFormula", "noiseFormula"},
        {"observableName", "observableTransformation", "noiseDistribution"},
    ),
    "measurements": (
        {
            "observableId",
            "simulationConditionId",
            "measurement",
            "time",
            "noiseParameters",
        },
        {
            "preequilibrationConditionId",
            "observableParameters",
            "datasetId",
            "replicateId",
        },
    ),
    "simulations": (
        {
            "observableId",
            "simulationConditionId",
            "simulation",
            "time",
            "noiseParameters",
        },
        {
            "preequilibrationConditionId",
            "observableParameters",
            "datasetId",
            "replicateId",
            "observableTransformation",
        },
    ),
}


class UniqueYaml(yaml.SafeLoader):
    """Safe YAML without duplicate keys, anchors or aliases."""

    def compose_node(self, parent, index):
        event = self.peek_event()
        require(
            not isinstance(event, yaml.AliasEvent)
            and getattr(event, "anchor", None) is None,
            "PEtab YAML anchors and aliases are not supported",
        )
        return super().compose_node(parent, index)

    def construct_mapping(self, node, deep=False):
        pairs = self.construct_pairs(node, deep)
        result = dict(pairs)
        require(len(result) == len(pairs), "duplicate PEtab YAML key")
        return result


def text(path):
    """Return the UTF-8 text of a regular file of at most `LIMIT` bytes."""
    require(
        path.is_file() and not path.is_symlink(),
        f"PEtab input {path.name} is not a regular file",
    )
    require(
        path.stat().st_size <= LIMIT, f"PEtab input {path.name} exceeds {LIMIT} bytes"
    )
    raw = path.read_bytes()
    require(len(raw) <= LIMIT, f"PEtab input {path.name} exceeds {LIMIT} bytes")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"PEtab input {path.name} is not UTF-8") from error


def number(value, what):
    try:
        result = float(value)
    except ValueError:
        raise ValueError(f"{what} is not a number: {value!r}") from None
    require(isfinite(result), f"{what} is not finite: {value!r}")
    return result


def table(path, kind):
    """Read one TSV table with exactly the supported columns of `kind`."""
    required, optional = COLUMNS[kind]
    content = text(path)
    if content.startswith("﻿"):
        content = content[1:]
    reader = csv.reader(io.StringIO(content, newline=""), delimiter="\t", strict=True)
    try:
        lines = [line for line in reader if line]
    except csv.Error as error:
        raise ValueError(f"malformed TSV in {path.name}: {error}") from None
    require(lines, f"empty TSV file {path.name}")
    header, body = lines[0], lines[1:]
    require(len(set(header)) == len(header), f"duplicate TSV column in {path.name}")
    missing = required - set(header)
    unknown = set(header) - required - optional
    require(not missing, f"{path.name} lacks columns {sorted(missing)}")
    require(
        not unknown,
        f"{path.name} has unsupported columns {sorted(unknown)}"
        + (" (condition parameter overrides)" if kind == "conditions" else ""),
    )
    require(body, f"{path.name} has no rows")
    require(len(body) <= ROWS, f"{path.name} has more than {ROWS} rows")
    for i, line in enumerate(body, start=2):
        require(
            len(line) == len(header),
            f"{path.name} line {i} has {len(line)} fields, expected {len(header)}",
        )
    return [dict(zip(header, line, strict=True)) for line in body]


def unique(rows, field, source):
    result = {}
    for row in rows:
        name = row[field]
        require(
            IDENTIFIER.match(name) is not None,
            f"{field} {name!r} in {source} is not an identifier",
        )
        require(name not in result, f"duplicate {field} {name!r} in {source}")
        result[name] = row
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
    """Join key of a measurement or simulation row."""
    return (
        row["observableId"],
        row.get("datasetId", ""),
        row["simulationConditionId"],
        number(row["time"], "time"),
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

    def scaled(self, parameters):
        """Inverse of `bind`: the estimated coordinates of a parameter mapping."""
        return np.array(
            [
                np.log10(parameters[name])
                if self.parameters[name]["parameterScale"] == "log10"
                else parameters[name]
                for name in self.fits
            ]
        )

    def score(self, readout, parameters):
        return gaussian_nll(
            self.y,
            readout[self.t_index, self.o_index],
            [parameters[name] for name in self.noise],
        )

    def score_archive(self, path):
        """Score archived predictions joined by observable, dataset, condition, time.

        Repeated measurements take the archived rows of their key in file order.
        """
        indexed = defaultdict(deque)
        for row in table(path, "simulations"):
            indexed[key(row)].append(
                (
                    number(row["simulation"], "archived simulation"),
                    number(row["noiseParameters"], "archived noise scale"),
                )
            )
        predictions = []
        for row in self.measurements:
            bucket = indexed[key(row)]
            require(bool(bucket), f"no archived simulation row for {key(row)}")
            predictions.append(bucket.popleft())
        require(not any(indexed.values()), "unused archived simulation rows")
        predicted, sigma = np.asarray(predictions).T
        return gaussian_nll(self.y, predicted, sigma)


def read(problem_path):
    problem_path = Path(problem_path)
    base = problem_path.resolve().parent
    config = yaml.load(text(problem_path), Loader=UniqueYaml)
    require(isinstance(config, dict), "the PEtab YAML file must be a mapping")
    require(
        set(config) == {"format_version", "parameter_file", "problems"},
        "unsupported PEtab problem fields: expected exactly format_version,"
        " parameter_file and problems",
    )
    require(
        config["format_version"] in (1, "1")
        and isinstance(config["problems"], list)
        and len(config["problems"]) == 1
        and isinstance(config["problems"][0], dict),
        "expected one PEtab v1 problem",
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
        require(isinstance(name, str) and name, "PEtab file names must be text")
        relative = Path(name)
        require(
            not relative.is_absolute()
            and ".." not in relative.parts
            and ":" not in name
            and "\\" not in name,
            f"PEtab path {name!r} must be relative to the problem folder",
        )
        path = base / relative
        for part in [path, *path.parents]:
            if part == base:
                break
            require(not part.is_symlink(), f"PEtab path {name!r} is a symbolic link")
        require(
            path.resolve().is_relative_to(base),
            f"PEtab path {name!r} leaves the problem folder",
        )
        return path

    params = unique(
        table(local(config["parameter_file"]), "parameters"),
        "parameterId",
        "parameters",
    )
    bounds, free = [], []
    for name, row in params.items():
        a = number(row["lowerBound"], f"lowerBound of {name}")
        b = number(row["upperBound"], f"upperBound of {name}")
        x = number(row["nominalValue"], f"nominalValue of {name}")
        scale = row["parameterScale"]
        require(
            scale in ("lin", "log10"),
            f"parameterScale {scale!r} of {name} is not supported (lin, log10)",
        )
        require(a < b, f"bounds of {name} are not increasing")
        require(a <= x <= b, f"nominalValue of {name} lies outside its bounds")
        require(a > 0 or scale == "lin", f"log10 parameter {name} needs bounds above 0")
        require(row["estimate"] in ("0", "1"), f"estimate of {name} must be 0 or 1")
        require(
            not any(row.get(field) for field in PRIORS),
            f"priors are not supported ({name})",
        )
        if row["estimate"] == "1":
            free.append(name)
            bounds.append(np.log10([a, b]) if scale == "log10" else [a, b])
    conditions = table(local(sub["condition_files"][0]), "conditions")
    require(len(conditions) == 1, "this adapter supports exactly one condition")
    condition = conditions[0]["conditionId"]
    observations = unique(
        table(local(sub["observable_files"][0]), "observables"),
        "observableId",
        "observables",
    )
    for name, row in observations.items():
        require(
            row.get("observableTransformation", "lin") in ("", "lin"),
            f"observableTransformation of {name} must be lin",
        )
        require(
            row.get("noiseDistribution", "normal") in ("", "normal"),
            f"noiseDistribution of {name} must be normal",
        )
        require(
            row["noiseFormula"] == f"noiseParameter1_{name}",
            f"noiseFormula of {name} must be noiseParameter1_{name}",
        )
        formula = row["observableFormula"]
        require(
            0 < len(formula) <= FORMULA and ARITHMETIC.match(formula) is not None,
            f"observableFormula of {name} must be arithmetic of at most"
            f" {FORMULA} characters",
        )
    rows = table(local(sub["measurement_files"][0]), "measurements")
    for i, row in enumerate(rows, start=2):
        where = f"measurement line {i}"
        require(
            not row.get("preequilibrationConditionId"),
            f"{where}: preequilibration is not supported",
        )
        require(
            not row.get("observableParameters"),
            f"{where}: observableParameters are not supported",
        )
        require(
            row["simulationConditionId"] == condition,
            f"{where}: unknown simulation condition",
        )
        require(row["observableId"] in observations, f"{where}: unknown observable")
        noise = row["noiseParameters"]
        require(
            noise in params,
            f"{where}: noiseParameters must name one parameter of the parameter table",
        )
        require(
            float(params[noise]["lowerBound"]) > 0,
            f"{where}: noise parameter {noise} needs a positive lower bound",
        )
        number(row["measurement"], f"{where}: measurement")
        require(number(row["time"], f"{where}: time") >= 0.0, f"{where}: negative time")
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
