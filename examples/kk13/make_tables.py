"""Tables of tutorial S2 from the result files of kowald_kirkwood_2013.py.

    python examples/kk13/make_tables.py                    # print the tables
    python examples/kk13/make_tables.py --write PAGE ...   # fill the pages
    python examples/kk13/make_tables.py --check PAGE ...   # exit 1 if a page differs

A table goes between `<!-- kk13-table: NAME -->` and `<!-- /kk13-table -->`
in a page. The tables show the run in examples/kk13/results/<size>/ (default
`draft`), with the author's 0.7.1 run and the paper for comparison.
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kowald_kirkwood_2013 as kk  # noqa: E402

# The author's run with CRN-craft 0.7.1 (his tutorial) and the paper (his
# tutorial's "Paper" columns).
AUTHOR_1A = {1e-4: "22.6", 1e-5: "2.4", 1e-6: "0.20"}
PAPER_1A = {1e-4: "22.3", 1e-5: "2.7", 1e-6: "≈ 0.25"}
AUTHOR_1B = {
    3: ("6.41e-3", "7.53e-3", "8.32e-3"),
    10: ("9.69e-4", "1.39e-3", "1.74e-3"),
    40: ("8.77e-5", "1.69e-4", "2.53e-4"),
    80: ("3.38e-5", "6.55e-5", "9.86e-5"),
    120: ("1.97e-5", "4.03e-5", "6.27e-5"),
}
PAPER_1B = {
    3: ("6.5e-3", "7.6e-3", "8.5e-3"),
    10: ("9.5e-4", "1.45e-3", "1.8e-3"),
    40: ("9e-5", "1.7e-4", "2.7e-4"),
    80: ("3e-5", "7e-5", "1.1e-4"),
    120: ("2e-5", "4e-5", "8e-5"),
}
AUTHOR_1C = {
    3: "31.3 / 35.9 / 38.4",
    10: "8.45 / 10.34 / 12.02",
    40: "2.25 / 2.95 / 3.29",
    80: "1.56 / 1.78 / 2.06",
    120: "1.31 / 1.51 / 1.66",
}
PAPER_1C = {
    3: "31.5 / 35.7 / 38.7",
    10: "7.8 / 10.4 / 12.2",
    40: "2.3 / 2.8 / 3.5",
    80: "1.44 / 1.75 / 2.08",
    120: "1.35 / 1.46 / 1.61",
}
AUTHOR_1D = {
    3: "17.3 / 15.6 / 14.3",
    10: "44.9 / 38.9 / 36.7",
    40: "80.2 / 76.4 / 74.4",
    80: "90.1 / 88.1 / 87.2",
    120: "95.5 / 92.3 / 91.6",
}
PAPER_1D = {
    3: "17 / 16 / 15",
    10: "44 / 40 / 37",
    40: "81 / 76 / 72",
    80: "92 / 90 / 86",
    120: "95 / 93 / 92",
}


def _sci(x):
    mantissa, exponent = f"{x:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def _mean_se(stat, digits):
    mean, se = stat
    if math.isnan(mean):
        return "–"
    return f"{mean:.{digits}f} ± {se:.{digits}f}"


def table_fig1a(out):
    rows = [
        "| P_mut | Cells | COX-negative at 120 y (this run) | Author's run (0.7.1) | Paper |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for pmut in kk.FIG1A_PMUT:
        path = out / f"fig1a_pmut{pmut:.0e}.csv"
        if not path.exists():
            continue
        info = json.loads(path.with_suffix(".json").read_text())
        last = kk._readCsv(path)[-1]
        n, k = info["cells"], int(last["cox_negative"])
        se = 100 * math.sqrt(max(k, 1) * (n - k) / n) / n
        exponent = round(math.log10(pmut))
        rows.append(
            f"| 10^{exponent} | {n:,} | {100 * k / n:.2f} % ± {se:.2f} | {AUTHOR_1A[pmut]} % | {PAPER_1A[pmut]} % |"
        )
    return "\n".join(rows)


def table_fig1b(out):
    rows = [
        "| Lifespan | 5 % (this run / author / paper) | 10 % | 15 % |",
        "| --- | --- | --- | --- |",
    ]
    for years in kk.LIFESPANS:
        path = out / f"fig1b_{years}y.json"
        if not path.exists():
            continue
        est = json.loads(path.read_text())["estimates"]
        cells = [
            f"{_sci(est[t]['pmut'])} / {a} / {p}"
            for t, a, p in zip(
                ["5", "10", "15"], AUTHOR_1B[years], PAPER_1B[years], strict=True
            )
        ]
        rows.append(f"| {years} y | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def _table_clones(out, column, digits, author, paper, extra=None):
    stats = kk.cloneStatistics(out)
    head = "| Lifespan | 5 % (this run) | 10 % | 15 % | COX-negative cells | Author's run | Paper |"
    if extra:
        head = head[:-1] + "| Mutant types only (this run) |"
    rows = [head, "|" + " --- |" * (head.count("|") - 1)]
    for years in kk.LIFESPANS:
        keys = [(years, p) for p in (5, 10, 15)]
        if not all(k in stats for k in keys):
            continue
        values = [_mean_se(stats[k][column], digits) for k in keys]
        counts = " / ".join(str(stats[k]["n"]) for k in keys)
        row = (
            f"| {years} y | "
            + " | ".join(values)
            + f" | {counts} | {author[years]} | {paper[years]} |"
        )
        if extra:
            row += (
                " "
                + " / ".join(f"{stats[k][extra][0]:.{digits}f}" for k in keys)
                + " |"
            )
        rows.append(row)
    return "\n".join(rows)


def table_fig1c(out):
    return _table_clones(out, "types", 2, AUTHOR_1C, PAPER_1C, extra="mut_types")


def table_fig1d(out):
    return _table_clones(out, "largest_pct", 1, AUTHOR_1D, PAPER_1D)


TABLES = {
    "fig1a": table_fig1a,
    "fig1b": table_fig1b,
    "fig1c": table_fig1c,
    "fig1d": table_fig1d,
}
MARKER = re.compile(r"(<!-- kk13-table: (\w+) -->\n)(.*?)(<!-- /kk13-table -->)", re.S)


def fill(text, out):
    return MARKER.sub(
        lambda m: m.group(1) + TABLES[m.group(2)](out) + "\n" + m.group(4), text
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--size", default="draft")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", nargs="+", type=Path, metavar="PAGE")
    mode.add_argument("--check", nargs="+", type=Path, metavar="PAGE")
    args = parser.parse_args()
    out = HERE / "results" / args.size
    if not (args.write or args.check):
        for name, table in TABLES.items():
            print(f"{name}:\n{table(out)}\n")
        return 0
    status = 0
    for page in args.write or args.check:
        text = page.read_text()
        new = fill(text, out)
        if args.write and new != text:
            page.write_text(new)
            print(f"{page}: tables written")
        elif args.check and new != text:
            print(
                f"{page}: a table differs from the result files; run --write",
                file=sys.stderr,
            )
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
