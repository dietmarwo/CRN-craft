"""Markdown tables of page O3b, rendered from the saved result files."""

import pandas as pd

from model import ROOT

RESULTS = ROOT / "results"


def score_table():
    scores = pd.read_csv(RESULTS / "scores.csv")
    archived = pd.read_csv(RESULTS / "review-data-choices.csv")
    archived = archived[archived.vector == "archived histogram fit"].set_index("cells")
    rows = [
        "| Evaluation set | Cells | Run size | Mean fit | Histogram fit | Archived histogram vector |",
        "| --- | ---: | --- | ---: | ---: | ---: |",
    ]
    full = scores[(scores.partition == "all") & (scores.preset == "full")]
    full = full.set_index("objective").histogram_nll_per_cell
    rows.append(
        f"| All usable cells, both doses (in sample) | 59,557 | full "
        f"| {full['mean']:.5f} | {full['histogram']:.5f} "
        f"| {archived.nll_per_cell['all replicates, recorded total']:.5f} |"
    )
    held = scores[scores.partition != "all"]
    for dose, part in held.groupby("held_dose"):
        part = part.set_index("objective")
        rep = str(int(part.held_rep.iloc[0]))
        n = int(part.held_n.iloc[0])
        reference = archived.nll_per_cell[f"replicate {rep}, recorded total"]
        rows.append(
            f"| {dose} M replicate {rep}, held out | {n:,} | draft "
            f"| {part.held_nll_per_cell['mean']:.5f} "
            f"| {part.held_nll_per_cell['histogram']:.5f} | {reference:.5f} |"
        )
    return "\n".join(rows) + "\n"


def start_table():
    runs = pd.read_csv(RESULTS / "review-optimizer-start.csv")
    rows = [
        "| Box | Seed | Evaluations | Invalid | Change vs archived (nats/cell) | Farthest parameter (decades) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in runs.itertuples():
        rows.append(
            f"| {r.box} | {r.seed} | {r.evaluations:,} | {r.invalid} "
            f"| {r.change_vs_archived:+.4f} | {r.max_abs_log10_from_archived:.2f} |"
        )
    return "\n".join(rows) + "\n"


TABLES = {"yeast-scores": score_table, "yeast-start": start_table}


def region(page, name):
    return page.partition(f"<!-- {name} -->\n")[2].partition(f"<!-- /{name} -->")[0]


if __name__ == "__main__":
    for name, make in TABLES.items():
        print(f"<!-- {name} -->\n{make()}<!-- /{name} -->\n")
