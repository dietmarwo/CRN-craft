#!/usr/bin/env python3
"""Regenerate the tutorial's terminal transcripts and figures.

Every terminal screenshot in docs/tutorial.md is a real command run against
the repository's examples with the native CLI, and every figure is computed
with the installed `crn_craft` package. Run from the repository root after
`cargo build --release --workspace --locked` and installing the wheel:

    python examples/tutorial/make_assets.py --binary target/release/crn-craft \
        --output docs/images/tutorial

Transcripts are written as HTML pages (`*.terminal.html`) next to the
figures; examples/tutorial/render_terminals.py turns them into PNG
screenshots with a headless Chromium.
"""

import argparse
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import crn_craft as crn  # noqa: E402
from crn_craft.observables import oscillation, sampled_first_passage  # noqa: E402
from workdir import prepare_workdir  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
# Categorical slots in fixed order (dataviz reference palette, light mode).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.titlecolor": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "legend.frameon": False,
            "font.size": 10,
            "lines.linewidth": 1.6,
            "savefig.dpi": 130,
        }
    )


# --------------------------------------------------------------------------
# Terminal transcripts


TERMINAL_CSS = """
body { margin: 0; background: #ffffff; font-family: 'DejaVu Sans Mono', monospace; }
.window { margin: 12px; border-radius: 10px; overflow: hidden; width: 940px;
  box-shadow: 0 4px 18px rgba(0,0,0,.25); background: #1e1f22; }
.bar { background: #2b2d31; padding: 9px 14px; display: flex; align-items: center; }
.dot { width: 12px; height: 12px; border-radius: 6px; margin-right: 7px; }
.title { color: #b8b9bf; font-family: 'DejaVu Sans', sans-serif; font-size: 13px;
  margin-left: 12px; }
pre { margin: 0; padding: 14px 18px 16px; color: #d4d4d8; font-size: 13.5px;
  line-height: 1.42; white-space: pre-wrap; word-break: break-all; }
.prompt { color: #6ad49b; } .cmd { color: #ffffff; font-weight: bold; }
.more { color: #8b8d94; font-style: italic; }
.key { color: #8fc1ff; } .str { color: #f2c27b; } .num { color: #c6a0f6; }
.lit { color: #f28fad; }
"""


def colorize(text):
    """Light JSON highlighting on already-escaped text."""
    import re

    out = []
    for line in text.split("\n"):
        line = re.sub(
            r"(&quot;[^&]*?&quot;)(:)", r'<span class="key">\1</span>\2', line
        )
        line = re.sub(
            r"(: )(&quot;.*?&quot;)(,?)$", r'\1<span class="str">\2</span>\3', line
        )
        line = re.sub(
            r"(: |^\s+)(-?\d[\d.eE+-]*)(,?)$", r'\1<span class="num">\2</span>\3', line
        )
        line = re.sub(
            r"(: |^\s+)(true|false|null)(,?)$", r'\1<span class="lit">\2</span>\3', line
        )
        out.append(line)
    return "\n".join(out)


class Terminal:
    def __init__(self, workdir, binary):
        self.workdir = workdir
        self.env = dict(os.environ)
        bindir = workdir / "bin"
        bindir.mkdir(exist_ok=True)
        shutil.copy2(binary, bindir / "crn-craft")
        self.env["PATH"] = f"{bindir}{os.pathsep}{self.env['PATH']}"
        self.env["PYTHONDONTWRITEBYTECODE"] = "1"

    def run(self, command, expect=0):
        proc = subprocess.run(
            ["bash", "-c", command],
            cwd=self.workdir,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
            timeout=600,
        )
        if expect is not None and proc.returncode != expect:
            raise RuntimeError(
                f"{command!r} exited {proc.returncode}:\n{proc.stdout}{proc.stderr}"
            )
        return proc.stdout + proc.stderr, proc.returncode

    def session(self, output, name, title, commands):
        """Run commands and write an HTML terminal page plus a text copy."""
        parts, text = [], []
        for item in commands:
            command, limit, expect = item if isinstance(item, tuple) else (item, 40, 0)
            out, code = self.run(command, expect)
            lines = out.rstrip("\n").split("\n") if out.strip() else []
            shown = lines[:limit]
            text.append(f"$ {command}\n" + "\n".join(lines))
            block = (
                f'<span class="prompt">$ </span>'
                f'<span class="cmd">{html.escape(command)}</span>\n'
            )
            block += colorize(html.escape("\n".join(shown)))
            if len(lines) > limit:
                block += (
                    f'\n<span class="more">… {len(lines) - limit} more lines</span>'
                )
            if code and expect:
                block += f'\n<span class="more">(exit status {code})</span>'
            parts.append(block)
        page = (
            f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<style>{TERMINAL_CSS}</style></head><body><div class='window'>"
            f"<div class='bar'><span class='dot' style='background:#ff5f57'></span>"
            f"<span class='dot' style='background:#febc2e'></span>"
            f"<span class='dot' style='background:#28c840'></span>"
            f"<span class='title'>{html.escape(title)}</span></div>"
            f"<pre>{chr(10).join(parts)}</pre></div></body></html>"
        )
        (output / f"{name}.terminal.html").write_text(page)
        (output / f"{name}.txt").write_text("\n\n".join(text) + "\n")


def sessions(terminal, output):
    t = terminal
    t.session(
        output,
        "01-help",
        "crn-craft --help",
        [("crn-craft --help", 60, 0), "crn-craft --version"],
    )
    t.session(
        output,
        "02-check",
        "Checking a model and its run",
        [
            "sha256sum immigration-death.model.json",
            "crn-craft check immigration-death.model.json immigration-death.run.json",
        ],
    )
    t.session(
        output,
        "03-diagnostics",
        "Diagnostics name the code, the JSON path and a hint",
        [
            ("crn-craft check missing-parameter.model.json", 40, 1),
            ("crn-craft check immigration-death.model.json stale-hash.run.json", 40, 1),
        ],
    )
    t.session(
        output,
        "04-simulate",
        "An exact ensemble written as a sealed archive",
        [
            "crn-craft run immigration-death.model.json"
            " immigration-death.run.json --archive runs/first",
            "ls runs/first && cat runs/first/COMPLETE",
            ("crn-craft verify runs/first", 12, 0),
        ],
    )
    t.session(
        output,
        "05-tamper",
        "verify detects an edited result",
        [
            'cp -r runs/first runs/edited && sed -i \'s/"events":[0-9]*/"events":1/\''
            " runs/edited/result.json",
            ("crn-craft verify runs/edited", 30, 1),
        ],
    )
    t.session(
        output,
        "06-replay",
        "The same request gives the same bytes",
        [
            "crn-craft run repressilator.model.json repressilator.run.json"
            " --archive runs/replay-a > /dev/null",
            "crn-craft run repressilator.model.json repressilator.run.json"
            " --archive runs/replay-b > /dev/null",
            "sha256sum runs/replay-a/result.json runs/replay-b/result.json",
        ],
    )
    t.session(
        output,
        "07-fsp",
        "Finite state projection of the immigration-death law",
        [
            "crn-craft run immigration-death.model.json immigration-death.fsp.run.json"
            " > fsp.json",
            "python3 -c \"import json; r = json.load(open('fsp.json'));"
            " s = r['snapshots'][-1]; s.pop('probabilities');"
            " print(len(r['states']), 'states,', r['sparse_edges'], 'edges');"
            ' print(json.dumps(s, indent=2))"',
        ],
    )
    t.session(
        output,
        "08-tau",
        "Approximate tau leaping is labelled as such",
        [
            "crn-craft run immigration-death.model.json tau.run.json > tau-result.json",
            "python3 -c \"import json; r = json.load(open('tau-result.json'));"
            " keys = ['method', 'sampler', 'termination', 'completed_replicates',"
            " 'events'];"
            ' print(json.dumps({k: r.get(k) for k in keys}, indent=2))"',
        ],
    )
    t.session(
        output,
        "09-cancel",
        "Cooperative cancellation with a cancel file",
        [
            "touch stop",
            (
                "crn-craft run immigration-death.model.json long.run.json"
                " --archive runs/stopped --cancel-file stop",
                30,
                1,
            ),
            ("crn-craft verify runs/stopped", 12, 0),
        ],
    )
    t.session(
        output,
        "10-import",
        "Importing a BioNetGen network",
        [
            "cat dimer.net",
            "crn-craft-import bng-net dimer.net dimer.import.json --convention binomial",
            "python3 -c \"import json; r = json.load(open('dimer.import.json'));"
            " print(json.dumps(r['model']['rules'], indent=2))\"",
        ],
    )
    t.session(
        output,
        "11-export",
        "Exporting to GillesPy2 source",
        [
            "crn-craft export enzyme-mechanism.model.json gillespy2"
            " | python3 -c \"import json, sys; print(json.load(sys.stdin)['artifact']['source'])\"",
        ],
    )


# --------------------------------------------------------------------------
# Figures


def ensemble(model, end, replicates, report_times=None, **options):
    """One exact Static Direct ensemble; the statistics have no quantiles."""
    if report_times is None:
        report_times = [0.0, end]
    return crn.ensemble(
        model,
        end=end,
        report_times=report_times,
        replicates=replicates,
        method=options.pop("method", "static_direct"),
        quantiles=(),
        **options,
    )


def final_counts(result, species):
    """The last saved count of every kept path."""
    return result.paths.sel(variable=species).to_numpy()[:, -1]


def immigration_death(output):
    model = crn.load_example("immigration-death")
    grid = np.linspace(0, 30, 301).round(10).tolist()
    result = ensemble(
        model, 30.0, 400, grid, keep_paths=8, seed=42, max_events=10_000_000
    ).simulate()
    t = np.array(grid)
    mean = result.statistics["mean"]["A"].to_numpy()
    variance = result.statistics["sample_variance"]["A"].to_numpy()
    exact = 20 * (1 - np.exp(-0.2 * t))
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 3.8))
    for path in result.paths.sel(variable="A").to_numpy():
        left.step(t, path, where="post", color=GRID, lw=1)
    left.plot(t, mean, color=BLUE, label="ensemble mean (400 paths)")
    left.plot(
        t, exact, color=ORANGE, ls="--", label=r"exact mean $20\,(1 - e^{-0.2t})$"
    )
    left.set(title="Eight paths and the ensemble mean", xlabel="time", ylabel="A")
    left.legend(loc="lower right")
    right.plot(t, variance, color=BLUE, label="sample variance")
    right.plot(t, exact, color=ORANGE, ls="--", label="exact variance (Poisson)")
    right.set(title="Variance equals the mean", xlabel="time", ylabel="variance of A")
    right.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(output / "immigration-death.png")
    plt.close(fig)
    return {"events": result.metadata["events"], "final_mean": float(mean[-1])}


def fsp_vs_ssa(output):
    model = crn.load_example("immigration-death")
    saves = [0.0, 2.0, 5.0, 20.0]
    projection = model.fsp(
        end=20.0,
        report_times=saves,
        bounds={"A": 80},
        max_states=81,
        tolerance=1e-10,
    )
    states = projection.states["A"].to_numpy()
    order = np.argsort(states)
    ssa = ensemble(
        model, 20.0, 2000, saves, keep_paths=2000, seed=42, max_events=10_000_000
    ).simulate()
    counts = ssa.paths.sel(variable="A").to_numpy().astype(int)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for axis, index in zip(axes, [1, 2, 3], strict=True):
        time = saves[index]
        p = projection.probabilities.isel(time=index).to_numpy()[order]
        mu = 20 * (1 - math.exp(-0.2 * time))
        k = states[order]
        poisson = np.exp(
            -mu + k * math.log(mu) - np.array([math.lgamma(x + 1) for x in k])
        )
        hist = np.bincount(counts[:, index], minlength=len(k))[: len(k)] / len(counts)
        axis.bar(k, hist, width=0.9, color=GRID, label="SSA, 2000 paths")
        axis.plot(k, p, color=BLUE, marker="o", ms=3, lw=1.2, label="FSP")
        axis.plot(k, poisson, color=ORANGE, ls="--", label="exact Poisson")
        axis.set(title=f"t = {time:g}", xlabel="A", xlim=(-0.5, 40))
    axes[0].set_ylabel("probability")
    axes[0].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(output / "fsp-vs-ssa.png")
    plt.close(fig)
    last = projection.probabilities.isel(time=-1).to_numpy()
    return {
        "max_abs_error_vs_poisson": float(
            max(
                abs(
                    pi
                    - math.exp(
                        -20 * (1 - math.exp(-4))
                        + s * math.log(20 * (1 - math.exp(-4)))
                        - math.lgamma(s + 1)
                    )
                )
                for s, pi in zip(states, last, strict=True)
            )
        ),
        "mass_defect": float(projection.truncation["mass_defect"].iloc[-1]),
    }


def repressilator(output):
    model = crn.load_example("repressilator")
    grid = np.linspace(0, 200, 1001).round(10).tolist()
    result = ensemble(
        model, 200.0, 1, grid, keep_paths=1, seed=42, max_events=10_000_000
    ).simulate()
    t = np.array(grid)
    path = result.paths.isel(replicate=0)

    def counts(name):
        return path.sel(variable=name).to_numpy()

    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(10, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1.2]}
    )
    for color, name in zip([BLUE, ORANGE, AQUA], ["p1", "p2", "p3"], strict=True):
        top.plot(t, counts(name), color=color, label=name, lw=1.2)
    verdict = oscillation(t, counts("p1"))
    top.set(
        title=(
            f"Repressilator, one exact path; p1 is {verdict['classification']}"
            f" (period ~{verdict['period']:.1f})"
        ),
        ylabel="protein count",
    )
    top.legend(loc="upper right", ncol=3)
    for color, name in zip([BLUE, ORANGE, AQUA], ["m1", "m2", "m3"], strict=True):
        bottom.plot(t, counts(name), color=color, lw=1)
    bottom.set(xlabel="time", ylabel="mRNA")
    fig.tight_layout()
    fig.savefig(output / "repressilator.png")
    plt.close(fig)
    return verdict


def tau_accuracy(output):
    """Tau leaping against exact SSA as epsilon varies (large-mean model)."""
    model = crn.Model("large_mean")
    model.add_species({"X": 1000})
    model.add_reaction(1000.0, [], ["X"], name="birth")
    model.add_reaction(1.0, ["X"], [], name="death")
    paths = 1000

    def sample(method):
        return ensemble(
            model,
            5.0,
            paths,
            keep_paths=paths,
            seed=7,
            max_events=10_000_000,
            method=method,
        ).simulate()

    exact_final = final_counts(sample("static_direct"), "X")
    rows = []
    approximate = {}
    for epsilon in [0.003, 0.01, 0.03, 0.1]:
        approx = sample({"kind": "tau", "epsilon": epsilon})
        final = final_counts(approx, "X")
        approximate[epsilon] = final
        rows.append(
            {
                "epsilon": epsilon,
                "variance_ratio": float(final.var(ddof=1) / exact_final.var(ddof=1)),
                "mean": float(final.mean()),
                "events": approx.metadata["events"],
            }
        )
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 3.6))
    eps = [r["epsilon"] for r in rows]
    left.plot(eps, [r["variance_ratio"] for r in rows], color=BLUE, marker="o")
    left.axhline(1.0, color=MUTED, lw=1, ls=":")
    left.axvline(1 / math.sqrt(1000), color=ORANGE, lw=1, ls="--")
    left.annotate(
        "relative s.d. of X\n= 1/sqrt(1000)",
        (1 / math.sqrt(1000), 1.0),
        xytext=(0.0335, 1.45),
        color=MUTED,
        fontsize=9,
    )
    left.set(
        xscale="log",
        title="Tau variance / exact variance at t = 5",
        xlabel="epsilon",
        ylabel="variance ratio",
    )
    bins = np.linspace(860, 1140, 41)
    right.hist(exact_final, bins=bins, color=GRID, label="exact SSA")
    for color, row, label in [(BLUE, 0, "tau, eps 0.003"), (ORANGE, 3, "tau, eps 0.1")]:
        right.hist(
            approximate[rows[row]["epsilon"]],
            bins=bins,
            histtype="step",
            color=color,
            lw=1.6,
            label=label,
        )
    right.set(title="Final counts, 1000 paths each", xlabel="X(5)", ylabel="paths")
    right.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(output / "tau-accuracy.png")
    plt.close(fig)
    return rows


def first_passage(output):
    model = crn.load_example("immigration-death")
    fine = np.linspace(0, 10, 201).round(10).tolist()
    result = ensemble(
        model, 10.0, 500, fine, keep_paths=500, seed=42, max_events=10_000_000
    ).simulate()
    paths = list(result.paths.sel(variable="A").to_numpy())
    rows = []
    for step in [1, 4, 20]:
        times = fine[::step]
        report = sampled_first_passage(times, [p[::step] for p in paths], 15)
        rows.append({"grid_spacing": 0.05 * step, **report})
    fig, axis = plt.subplots(figsize=(6.5, 3.4))
    spacing = [r["grid_spacing"] for r in rows]
    axis.plot(
        spacing, [r["restricted_mean_upper"] for r in rows], color=BLUE, marker="o"
    )
    axis.set(
        xscale="log",
        title="First passage to A >= 15: coarser grids give looser bounds",
        xlabel="saved-grid spacing",
        ylabel="upper bound on E[min(T, 10)]",
    )
    fig.tight_layout()
    fig.savefig(output / "first-passage.png")
    plt.close(fig)
    return rows


def threads_identity():
    model = crn.load_example("repressilator")
    digests = {}
    for threads in [1, 2, 8]:
        result = ensemble(
            model, 50.0, 64, keep_paths=4, threads=threads, seed=42
        ).simulate()
        # The run hash covers `threads`; everything numerical must agree.
        numerical = {k: v for k, v in result.metadata.items() if k != "run_sha256"}
        digests[threads] = hashlib.sha256(
            json.dumps(numerical, sort_keys=True).encode()
        ).hexdigest()[:16]
    return digests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    with tempfile.TemporaryDirectory(prefix="crn-tutorial-") as temp:
        workdir = Path(temp)
        prepare_workdir(workdir)
        terminal = Terminal(workdir, args.binary.resolve())
        importer = shutil.which("crn-craft-import")
        if importer is None:
            raise RuntimeError("install the crn-craft wheel first")
        sessions(terminal, args.output)
    numbers = {
        "immigration_death": immigration_death(args.output),
        "fsp": fsp_vs_ssa(args.output),
        "repressilator": repressilator(args.output),
        "tau": tau_accuracy(args.output),
        "first_passage": first_passage(args.output),
        "threads": threads_identity(),
    }
    (args.output / "numbers.json").write_text(
        json.dumps(numbers, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(numbers, indent=1)[:3000])


if __name__ == "__main__":
    main()
