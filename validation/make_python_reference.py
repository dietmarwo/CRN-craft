"""Generate docs/reference/python-api.md from the installed `crn_craft` package.

The page lists the public names of `crn_craft` (`__all__`), of `crn_craft.Model`,
of the result and ensemble classes (`PreparedEnsemble`, `Result`, `FspResult`,
`CompiledLoop`) and of `crn_craft.importers`, grouped by task. Each entry shows
the signature, the first paragraph of the docstring, a table of the arguments
and the rest of the docstring. The argument descriptions are the table
`ARGUMENTS` below: the docstrings give the arguments of a call in prose, not
per argument.

    python validation/make_python_reference.py            # print the page
    python validation/make_python_reference.py --write    # write the page
    python validation/make_python_reference.py --check    # exit 1 on a difference

`--check` never writes. A public name that no group lists, a group that lists
a name that does not exist, and an argument without a description all fail the
run, so the page cannot fall behind the package. Problems with docstrings
(missing, or a reference that does not resolve) are printed to stderr with
`--lint` and never change the page.

The page is generated from the installed package: install the wheel built from
this checkout first (CI does).
"""

import argparse
import difflib
import inspect
from pathlib import Path
import re
import sys

import crn_craft
from crn_craft import importers

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "reference" / "python-api.md"

# ---------------------------------------------------------------------------
# What the page lists, by task. A path is a dotted name below `crn_craft`;
# "importers.x" is a name of `crn_craft.importers`. Result classes live in
# `crn_craft.rules` and are reached through `crn_craft.rules`.
# ---------------------------------------------------------------------------

GROUPS = [
    (
        "Models",
        "Build a model in Python, or take a packaged one. `Model` writes "
        "`crn-craft.model/1.0`. The language of the document is in the "
        "[model reference](model-language.md); chapter 2 of the guide teaches "
        "the builder.",
        [
            "Model",
            "Model.add_species",
            "Model.add_parameter",
            "Model.add_parameters",
            "Model.add_reaction",
            "Model.set_initial",
            "Model.set_stop",
            "Model.add_ode",
            "Model.add_ODE",
            "Model.add_continuous",
            "Model.from_dict",
            "Model.to_dict",
            "Model.to_json",
            "Model.name",
            "Model.species",
            "Model.parameters",
            "Model.reaction_names",
            "load_example",
        ],
    ),
    (
        "Running",
        "Prepare a request and simulate it, check documents without running "
        "them, ask for the finite state projection, and read the discovery "
        "document. `ensemble` builds the request from keyword arguments; "
        "`prepare` takes a request you wrote. Both return a "
        "`PreparedEnsemble`. The guide teaches this in chapters 3 and 5.",
        [
            "ensemble",
            "prepare",
            "PreparedEnsemble",
            "PreparedEnsemble.simulate",
            "PreparedEnsemble.bind_parameters",
            "PreparedEnsemble.rule_names",
            "Model.fsp",
            "check",
            "describe",
            "CrnError",
        ],
    ),
    (
        "Results",
        "What `simulate()` and `Model.fsp` return. The fields of the result "
        "document and the terminations are in the [results reference]"
        "(results.md). `Result.save` is under "
        "[files and archives](#files-and-archives). Chapter 4 of the guide "
        "teaches the data frames.",
        [
            "Result",
            "Result.paths",
            "Result.paths_frame",
            "Result.to_xarray",
            "Result.tracker_columns",
            "Result.event_columns",
            "Result.event_frame",
            "Result.model",
            "Result.request",
            "FspResult",
        ],
    ),
    (
        "Files and archives",
        "Read model and report files, write reports, and seal a run in an "
        "archive directory that can be verified later. The formats are in "
        "the [archive reference](archives.md) and the "
        "[contracts](contracts.md).",
        [
            "load",
            "load_report",
            "Result.save",
            "archive",
            "PreparedEnsemble.archive",
            "verify_archive",
        ],
    ),
    (
        "Analysis and import/export",
        "Analyze the structure of a network, export it, and import SBML and "
        "BioNetGen files. The reports are described in the [analysis "
        "reference](analysis.md) and the [GillesPy2 export]"
        "(gillespy2-export.md). Chapter 11 of the guide teaches this.",
        [
            "Model.analyze",
            "Model.to_gillespy2",
            "convert",
            "read_sbml",
            "read_bngl",
            "importers.import_sbml",
            "importers.import_bng_net",
            "importers.ImportFailure",
        ],
    ),
    (
        "Compiled loops",
        "Build the Static Direct event loop of a model into a program and run "
        "requests with it. This needs a Rust toolchain; the requirements, "
        "the cache and the limits are in the [methods reference]"
        "(methods.md#compiled-event-loops). Chapter 5 of the guide teaches "
        "when to use it.",
        [
            "compile_loop",
            "loop_source",
            "CompiledLoop",
            "CompiledLoop.execute",
            "CompiledLoop.verify",
            "CompiledLoopError",
        ],
    ),
]

# Public names that are not part of the page, with the reason.
EXCLUDED = {
    "importers.main": "command-line entry point of `python -m crn_craft.importers`",
    "importers.fail": "internal helper",
}

# ---------------------------------------------------------------------------
# Text for what a docstring does not give.
# ---------------------------------------------------------------------------

# Descriptions of arguments. A key "qualified.name.argument" wins over "argument".
ARGUMENTS = {
    "model": "A model: bytes from `load`, a JSON string, a dictionary or a `Model`.",
    "request": "A request: bytes, a JSON string or a dictionary.",
    "end": "Horizon of the run. Report times lie in [0, `end`].",
    "report_times": "Increasing report times in [0, `end`].",
    "replicates": "Number of replicates.",
    "seed": "Root seed. Replicate `i` draws from its own stream, so results do not depend on `threads`.",
    "threads": "Native worker threads, 1 to 64.",
    "keep_paths": 'How many replicates, from replicate 0, keep their full observable series (`result.paths`). `"all"` or `True` keeps every replicate.',
    "keep_instances": 'How many replicates keep their final instances (`metadata["final_instances"]`).',
    "max_events": "Event budget per replicate.",
    "max_sampling_steps": "Budget per replicate for random words and leap iterations.",
    "max_instances": "Most instances alive at once, per replicate.",
    "max_output_bytes": "Budget for the size of the result. Event logs reserve space for every replicate before the run.",
    "wall_time_s": "Wall-clock limit of the whole run in seconds, or `None`.",
    "trackers": 'A list of `{"observable", "first_passage"?, "keep_episodes"?}`.',
    "event_log": '`{"max_entries": n, "rules": [...]}`, or `None`.',
    "method": 'A method name or a dictionary such as `{"kind": "tau", "epsilon": 0.03}`; `None` runs the tree method. See [methods](methods.md).',
    "compiled": "A loop from `compile_loop` for this model. Checked at preparation and used by every `simulate()`.",
    "quantiles": 'Probabilities for `statistics["quantiles"]` and `statistics["median"]`. `()` or `None` asks for none.',
    "selector": 'The event selector of `static_direct`: `"linear"`, `"composition_rejection"`, `"sorting"`, `"rssa_cr"` or `"auto"`. See [methods](methods.md#selectors).',
    "cancel_file": "Path of a file; the run is cancelled when the file exists.",
    "directory": "The archive directory.",
    "path": "A file path.",
    "name": "A name.",
    "kind": 'A contract name such as `"model-1.0"` or `"run-1.0"`. `None` returns the discovery document.',
    "convention": "The rate convention. The importer refuses an absent or other value.",
    "raw": "The file's bytes.",
    "rate_convention": '`"binomial"` or `"falling_factorial"`.',
    # Model
    "Model.name": "The model's name.",
    "Model.add_species.name": "A species name, or a mapping from names to initial counts.",
    "Model.add_species.initial": "Initial count. Not allowed with a mapping.",
    "Model.add_parameter.name": "A parameter name, or a mapping from names to values.",
    "Model.add_parameter.value": "The value. Not allowed with a mapping.",
    "Model.add_parameters.name": "A parameter name, or a mapping from names to values.",
    "Model.add_parameters.value": "The value. Not allowed with a mapping.",
    "Model.add_reaction.rate": "A number (a mass-action constant) or a string (the whole propensity).",
    "Model.add_reaction.reactants": 'Species names, repeated for multiplicity: `["A", "A"]`.',
    "Model.add_reaction.products": "Species names, repeated for multiplicity.",
    "Model.add_reaction.name": "The rule's name. Default `r<index>`.",
    "Model.set_initial.counts": "A mapping from species names to counts.",
    "Model.set_initial.kwargs": "The same, as keywords: `set_initial(S=990, I=10)`.",
    "Model.set_stop.condition": "An expression over parameters and species counts, or `None` to remove the stop.",
    "Model.to_json.indent": "JSON indentation.",
    "Model.from_dict.document": "A model document: `crn-craft.model/1.0`, a rules-model or a network `model/0.1`.",
    "Model.fsp.bounds": "A mapping from every species label to its largest count.",
    "Model.fsp.max_states": "Cap on the enumerated states. `None` means 100,000.",
    "Model.fsp.tolerance": "Truncation budget of the series. `None` means `1e-12`.",
    "Model.add_ode.name": "The variable's name.",
    "Model.add_ode.initial": "A number or a parameter expression.",
    "Model.add_ode.derivative": "An expression over parameters, species counts and continuous variables.",
    "Model.add_ode.independent": "`True` if `initial` and `derivative` read only parameters and other independent variables.",
    "Model.add_ode.lower_bound": "A number or a parameter expression, or `None`.",
    "Model.add_ode.on_bound": '`"stop"` ends the replicate at the bound, `"clamp"` holds the variable there.',
    "load_example.name": '`"enzyme-mechanism"`, `"immigration-death"`, `"michaelis-menten"`, `"repressilator"` or `"vilar-noise"`.',
    # running
    "check.model": "A model: bytes, text, a dictionary or a `Model`.",
    "check.request": "A request for the model, or `None` to check the model alone.",
    "prepare.request": "A request you wrote: bytes, a JSON string or a dictionary.",
    "prepare.quantiles": "Unset: `(0.025, 0.5, 0.975)` unless the request states `outputs.quantiles` itself. `()` asks for none.",
    "PreparedEnsemble.simulate.callback": "`callback(id, event)` on the native workers; returning `False` stops that replicate. Needs the interpreter.",
    "PreparedEnsemble.simulate.batch_callback": "`batch_callback(id, events)` receives the retained log after the run.",
    "PreparedEnsemble.simulate.batch_size": "Entries per call of `batch_callback`.",
    "PreparedEnsemble.simulate.compiled": "A `CompiledLoop` for this call only. `None` uses the loop attached at preparation, if any.",
    "PreparedEnsemble.bind_parameters.overrides": "A mapping from parameter names to new finite values.",
    "PreparedEnsemble.archive.directory": "The new archive directory. Its parent must exist.",
    "archive.directory": "The new archive directory. It must not exist; its parent must.",
    "verify_archive.directory": "The archive directory to check.",
    # results
    "Result.event_columns.replicate": "A replicate id, or `None` for every log.",
    "Result.event_frame.replicate": "A replicate id, or `None` for every logged event.",
    "Result.tracker_columns.index": "The index of the tracker in the request's `trackers`.",
    "Result.save.path": "The file. `.md` or `.markdown` writes markdown, `.json` writes `crn-craft.report/1.0`.",
    "Result.save.paths": "`False` leaves the single trajectories out.",
    "load.path": "A model file.",
    "load_report.path": "A JSON file written by `Result.save`.",
    "read_sbml.path": "An SBML file.",
    "read_sbml.convention": 'Must be `"stochastic_count_propensity"`.',
    "read_bngl.path": "An expanded `.net` file.",
    "read_bngl.convention": '`"binomial"` or `"falling_factorial"`.',
    "importers.import_sbml.raw": "The SBML document's bytes.",
    "importers.import_sbml.convention": 'Must be `"stochastic_count_propensity"`.',
    "importers.import_bng_net.raw": "The `.net` file's bytes.",
    # compiled loops
    "compile_loop.cache_dir": "Where to keep builds. Default: the user cache directory, or `CRN_CRAFT_CACHE_DIR`.",
    "compile_loop.cargo": "The cargo executable.",
    "compile_loop.source": "A CRN-craft checkout of the installed version, for an offline build or a version that is not on crates.io. Default: `CRN_CRAFT_SOURCE`.",
    "compile_loop.progress": "`False` silences cargo's progress on stderr.",
    "loop_source.model": "A model: bytes, text, a dictionary or a `Model`.",
    "CompiledLoop.execute.request": "A rules-run request for this loop's model.",
    "CompiledLoop.verify.run": "A `PreparedEnsemble` for this loop's model.",
    "CrnError.diagnostics": "The list of diagnostics.",
    "CompiledLoopError.message": "The message.",
    "CompiledLoopError.reason": '`"cargo_missing"`, `"build_failed"` or `"mismatch"`.',
    "CompiledLoopError.log": "The path of the build log, if a build started.",
}

# Methods that take the arguments of another.
ALIASES = {"Model.add_ODE": "Model.add_ode", "Model.add_continuous": "Model.add_ode"}

# First paragraphs for names whose class or function has no docstring of its own.
MISSING_DOC = {
    "PreparedEnsemble": (
        "A prepared request, returned by `ensemble` and `prepare`. `simulate()` "
        "runs it, `bind_parameters` makes a candidate with other parameter "
        "values, and `archive` seals a run. It can be simulated again with "
        "independent execution state."
    ),
    "importers.ImportFailure": (
        "Raised by an importer for input it does not support. A `ValueError` "
        'with the attribute `code` `"CRN-E090"`.'
    ),
}

KINDS = {
    "class": "class",
    "function": "function",
    "method": "method",
    "property": "property",
}

REFERENCE = re.compile(r":(?:func|class|meth|data|exc|attr):`~?([^`]+)`")
LITERAL = re.compile(r"``([^`]+)``")


def resolve(path):
    """The object and its kind for a path of GROUPS."""
    if path.startswith("importers."):
        parent, name = importers, path.split(".", 1)[1]
        obj = getattr(parent, name)
        return obj, ("class" if inspect.isclass(obj) else "function"), None
    if path in ("PreparedEnsemble", "Result", "FspResult", "CompiledLoop"):
        from crn_craft import rules

        return getattr(rules, path), "class", None
    head, _, tail = path.partition(".")
    if not tail:
        obj = getattr(crn_craft, head)
        return obj, ("class" if inspect.isclass(obj) else "function"), None
    classes = {
        "Model": crn_craft.Model,
        "PreparedEnsemble": _rules().PreparedEnsemble,
        "Result": _rules().Result,
        "CompiledLoop": _rules().CompiledLoop,
    }
    cls = classes[head]
    raw = inspect.getattr_static(cls, tail)
    if isinstance(raw, property):
        return raw, "property", cls
    return getattr(cls, tail), "method", cls


def _rules():
    from crn_craft import rules

    return rules


def public_names():
    """Every public name that the page has to list or exclude."""
    names = {n for n in crn_craft.__all__ if n != "Model"}
    names.add("Model")
    for cls_name, cls in (
        ("Model", crn_craft.Model),
        ("PreparedEnsemble", _rules().PreparedEnsemble),
        ("Result", _rules().Result),
        ("FspResult", _rules().FspResult),
        ("CompiledLoop", _rules().CompiledLoop),
    ):
        names.add(cls_name)
        for member, value in inspect.getmembers(cls):
            if member.startswith("_"):
                continue
            raw = inspect.getattr_static(cls, member)
            if inspect.isroutine(value) or isinstance(raw, property):
                names.add(f"{cls_name}.{member}")
    for member, value in vars(importers).items():
        if member.startswith("_"):
            continue
        if getattr(value, "__module__", None) == importers.__name__:
            names.add(f"importers.{member}")
    return names


def doc_of(obj):
    """The object's own docstring, cleaned, or None. Inherited ones do not count."""
    if isinstance(obj, property):
        obj = obj.fget
    doc = obj.__doc__
    if doc is None:
        return None
    if inspect.isclass(obj) and doc == getattr(obj.__mro__[1], "__doc__", None):
        return None
    return inspect.cleandoc(doc)


def clean_signature(obj, kind, drop_self=True):
    try:
        signature = inspect.signature(obj)
    except (TypeError, ValueError):
        return None
    params = list(signature.parameters.values())
    if kind in ("method",) and params and params[0].name in ("self", "cls"):
        params = params[1:]
    return signature.replace(parameters=params)


def default_text(parameter):
    if parameter.kind in (
        inspect.Parameter.VAR_KEYWORD,
        inspect.Parameter.VAR_POSITIONAL,
    ):
        return "none"
    if parameter.default is inspect.Parameter.empty:
        return "required"
    text = repr(parameter.default)
    if text == "<required>":
        return "required"
    if re.fullmatch(r"<object object at 0x[0-9a-f]+>", text):
        return "unset"
    return f"`{text}`"


def render_signature(display, signature):
    text = re.sub(r"<object object at 0x[0-9a-f]+>", "<unset>", f"{display}{signature}")
    if len(text) <= 78:
        return text
    lines = [f"{display}("]
    for parameter in signature.parameters.values():
        piece = str(parameter)
        piece = re.sub(r"<object object at 0x[0-9a-f]+>", "<unset>", piece)
        lines.append(f"    {piece},")
    lines.append(")")
    return "\n".join(lines)


def slug(heading):
    text = heading.lower().replace("`", "")
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def convert_text(text, anchors, context, unresolved):
    """Docstring reStructuredText to markdown: literals, references, literal blocks."""

    def reference(match):
        target = match.group(1).strip()
        key = target
        if key.startswith("crn_craft."):
            key = key[len("crn_craft.") :]
        candidates = [key]
        if key.startswith("rules."):
            candidates.append(key[len("rules.") :])
        if context:
            candidates.append(f"{context}.{key}")
        for candidate in candidates:
            if candidate in anchors:
                return f"[`{target.lstrip('~')}`](#{anchors[candidate]})"
        unresolved.append(target)
        return f"`{target}`"

    out, lines, i = [], text.split("\n"), 0
    previous_blank = True
    while i < len(lines):
        line = lines[i]
        if line.rstrip().endswith("::") and not line.strip() == "::":
            line = line.rstrip()[:-1]
            out.append(line)
            out.append("")
            i += 1
            while i < len(lines) and lines[i].strip() == "":
                i += 1
            block = []
            while i < len(lines) and (
                lines[i].startswith(" ") or lines[i].strip() == ""
            ):
                block.append(lines[i])
                i += 1
            while block and block[-1].strip() == "":
                block.pop()
            out.append("```text")
            out.extend(inspect.cleandoc("\n".join(block)).split("\n"))
            out.append("```")
            previous_blank = False
            continue
        if (
            line.startswith("- ")
            and not previous_blank
            and not out[-1].startswith(("- ", "  "))
        ):
            out.append("")
        out.append(line)
        previous_blank = line.strip() == ""
        i += 1
    joined = "\n".join(out)
    joined = REFERENCE.sub(reference, joined)
    joined = LITERAL.sub(lambda m: f"`{m.group(1)}`", joined)
    return joined


def first_paragraph(doc):
    parts = re.split(r"\n\s*\n", doc, maxsplit=1)
    return parts[0], (parts[1] if len(parts) > 1 else "")


def build(lint=None):
    """The page text. `lint` is a list that collects docstring problems."""
    problems = lint if lint is not None else []
    listed = [path for _, _, paths in GROUPS for path in paths]
    if len(listed) != len(set(listed)):
        duplicate = sorted({p for p in listed if listed.count(p) > 1})
        raise SystemExit(f"listed twice: {duplicate}")
    existing = public_names()
    missing = sorted(existing - set(listed) - set(EXCLUDED))
    if missing:
        raise SystemExit(
            "public names that no group lists (add them to GROUPS or EXCLUDED): "
            + ", ".join(missing)
        )
    unknown = sorted(set(listed) - existing)
    if unknown:
        raise SystemExit("listed names that do not exist: " + ", ".join(unknown))

    anchors = {}
    for path in listed:
        anchors[path] = slug(path)
    unresolved = []
    undocumented_args = []

    body = []
    for title, intro, paths in GROUPS:
        body.append(f"## {title}\n")
        body.append(intro + "\n")
        for path in paths:
            obj, kind, owner = resolve(path)
            context = (
                path.split(".")[0]
                if "." in path and not path.startswith("importers.")
                else None
            )
            display = path
            doc = doc_of(obj)
            if doc is None:
                if path in MISSING_DOC:
                    doc = MISSING_DOC[path]
                    problems.append(
                        f"{path}: no docstring (the page uses its own sentence)"
                    )
                else:
                    problems.append(f"{path}: no docstring")
                    doc = "(no description)"
            head, rest = first_paragraph(doc)
            body.append(f"### `{display}`\n")
            if kind == "property":
                body.append(f"```text\n{display}\n```\n")
                body.append("Read-only property.\n")
            elif kind == "class" and path in (
                "PreparedEnsemble",
                "Result",
                "FspResult",
                "CompiledLoop",
            ):
                body.append(f"```text\nclass {display}\n```\n")
                body.append(
                    "Created by the functions of this page, not by calling it.\n"
                )
                signature = None
            else:
                target = obj
                signature = clean_signature(
                    target, "method" if kind == "method" else kind
                )
                label = ("class " if kind == "class" else "") + display
                if signature is None:
                    body.append(f"```text\n{label}\n```\n")
                else:
                    body.append(
                        "```text\n" + render_signature(label, signature) + "\n```\n"
                    )
            body.append(convert_text(head, anchors, context, unresolved) + "\n")
            # Arguments.
            if kind in ("function", "method", "class") and path not in (
                "PreparedEnsemble",
                "Result",
                "FspResult",
                "CompiledLoop",
            ):
                signature = clean_signature(obj, "method" if kind == "method" else kind)
                if signature is not None and signature.parameters:
                    rows = []
                    for parameter in signature.parameters.values():
                        key = f"{path}.{parameter.name}"
                        alias = f"{ALIASES.get(path, path)}.{parameter.name}"
                        note = ARGUMENTS.get(
                            key, ARGUMENTS.get(alias, ARGUMENTS.get(parameter.name))
                        )
                        if note is None:
                            undocumented_args.append(key)
                            note = ""
                        name = parameter.name
                        if parameter.kind is inspect.Parameter.VAR_KEYWORD:
                            name = f"**{name}"
                        elif parameter.kind is inspect.Parameter.VAR_POSITIONAL:
                            name = f"*{name}"
                        only = (
                            " (keyword only)"
                            if parameter.kind is inspect.Parameter.KEYWORD_ONLY
                            else ""
                        )
                        rows.append(
                            f"| `{name}`{only} | {default_text(parameter)} | {note} |"
                        )
                    body.append(
                        "| Argument | Default | Meaning |\n| --- | --- | --- |\n"
                        + "\n".join(rows)
                        + "\n"
                    )
            if rest.strip():
                body.append(convert_text(rest, anchors, context, unresolved) + "\n")
    if undocumented_args:
        raise SystemExit(
            "arguments without a description in ARGUMENTS: "
            + ", ".join(sorted(set(undocumented_args)))
        )
    for target in sorted(set(unresolved)):
        problems.append(f"reference not on this page: {target}")

    module_doc = inspect.cleandoc(importers.__doc__ or "")
    version = crn_craft.__version__
    header = [
        "# Python API\n",
        f"<!-- Generated by validation/make_python_reference.py from crn-craft {version}. "
        "Do not edit: run the script with --write. -->\n",
        f"The public names of `crn_craft` {version}: everything in `crn_craft.__all__`, the methods of "
        "`Model`, `PreparedEnsemble`, `Result`, `FspResult` and `CompiledLoop`, and "
        "`crn_craft.importers`. The page is generated from the signatures and "
        "docstrings of the installed package, so it lists what the package does. "
        "Examples are in the [guide](../guide/01-install.md); this page has none.\n",
        "Import the package as `import crn_craft as crn`. `crn_craft.rules` is the "
        "API of release 0.6.0 under its old name. It has the same names with the "
        "defaults and document versions of 0.6.0; the [migration page]"
        "(../guide/migrating.md) lists the differences. Code written for it "
        "keeps running.\n",
        "`crn_craft.rules.ensemble` takes the arguments of `ensemble` below. Its "
        "`quantiles` defaults to `None` (no quantiles), and an absent `selector` "
        "means the linear scan.\n",
        "Contents: " + ", ".join(f"[{t}](#{slug(t)})" for t, _, _ in GROUPS) + ".\n",
    ]
    footer = [
        "## The importers module\n",
        convert_text(module_doc, anchors, None, []) + "\n",
    ]
    return "\n".join(header + body + footer).rstrip() + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Generate docs/reference/python-api.md from the installed crn_craft package."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="write the page")
    mode.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if the committed page differs from the generated one",
    )
    parser.add_argument(
        "--lint", action="store_true", help="print docstring problems to stderr"
    )
    parser.add_argument(
        "--output", type=Path, default=PAGE, help="the page (default: %(default)s)"
    )
    args = parser.parse_args()

    problems = []
    text = build(problems)
    if args.lint:
        for problem in problems:
            print(f"docstring: {problem}", file=sys.stderr)
    if args.check:
        current = (
            args.output.read_text(encoding="utf-8") if args.output.exists() else ""
        )
        if current != text:
            diff = difflib.unified_diff(
                current.splitlines(),
                text.splitlines(),
                fromfile=f"{args.output} (committed)",
                tofile=f"{args.output} (generated)",
                lineterm="",
            )
            for count, line in enumerate(diff):
                if count >= 60:
                    print("...", file=sys.stderr)
                    break
                print(line, file=sys.stderr)
            print(
                f"{args.output} differs from the generated page: run "
                "python validation/make_python_reference.py --write",
                file=sys.stderr,
            )
            return 1
        print(f"{args.output}: up to date")
        return 0
    if args.write:
        args.output.write_text(text, encoding="utf-8")
        print(f"wrote {args.output}")
        return 0
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
