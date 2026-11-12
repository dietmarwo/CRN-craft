"""CRN-craft: stochastic simulation, exact distributions and analysis of
reaction networks, with one engine and one API.

::

    import crn_craft as crn

    m = crn.Model("sir")
    m.add_species({"S": 999, "I": 1, "R": 0})
    m.add_parameter({"beta": 2e-4, "gamma": 0.01})
    m.add_reaction("beta * S * I", ["S", "I"], ["I", "I"])
    m.add_reaction("gamma * I", ["I"], ["R"])
    result = crn.ensemble(m, end=200.0, report_times=[0, 100, 200],
                          replicates=200).simulate()
    result.statistics["median"]      # also "quantiles": {0.025: ..., 0.975: ...}
    result.save("sir.md")            # or "sir.json"
    dist = m.fsp(end=10.0, report_times=[10.0],
                 bounds={"S": 999, "I": 1000, "R": 1000})
    report = m.analyze()

The names here write the 1.0 contracts: ``Model.to_dict()`` is
``crn-craft.model/1.0``, ``ensemble`` and ``prepare`` write
``crn-craft.run/1.0`` and ask for the quantiles 0.025, 0.5 and 0.975 by
default, and results are ``crn-craft.result/1.0``. ``crn_craft.rules`` is
the API of release 0.6.0 under the same names: it writes the 0.6.0 contract
versions unless a call needs a 1.0 feature, so 0.6.0 code keeps its
documents, hashes and results.
"""

import json as _json
from pathlib import Path as _Path

from . import _native
from . import rules
from ._common import CrnError, convert, model_bytes as _model_bytes
from ._common import __version__ as __version__
from . import _rules_loop
from ._rules_loop import CompiledLoopError
from .builder import Model as _RulesModel
from .report import load_report
from .rules import DEFAULT_QUANTILES as _DEFAULT_QUANTILES
from .rules import archive, verify_archive

__all__ = [
    "CompiledLoopError",
    "CrnError",
    "Model",
    "archive",
    "check",
    "compile_loop",
    "convert",
    "describe",
    "ensemble",
    "load",
    "load_example",
    "load_report",
    "loop_source",
    "prepare",
    "read_bngl",
    "read_sbml",
    "verify_archive",
]

# The examples packaged as crn-craft.model/1.0 (python/crn_craft/models/).
_EXAMPLES = (
    "enzyme-mechanism",
    "immigration-death",
    "michaelis-menten",
    "repressilator",
    "vilar-noise",
)


class Model(_RulesModel):
    """A reaction network over attribute-free species, written as
    ``crn-craft.model/1.0``.

    The builder of ``crn_craft.rules.Model`` (the same methods: ``add_species``,
    ``add_parameter``, ``add_reaction``, ``set_initial``, ``set_stop``,
    ``add_ode``, ``to_dict``, ``to_json``, ``from_dict``, ``fsp``,
    ``analyze``, ``to_gillespy2``), whose documents are always
    ``crn-craft.model/1.0``. ``add_reaction`` accepts any number of
    reactant molecules; ``["A", "A", "A"]`` is a third-order reaction.
    """

    _unified = True


def ensemble(
    model,
    *,
    end,
    report_times,
    replicates=1,
    seed=0,
    threads=1,
    keep_paths=0,
    keep_instances=0,
    max_events=1_000_000,
    max_sampling_steps=8_000_000,
    max_instances=100_000,
    max_output_bytes=64 << 20,
    wall_time_s=None,
    trackers=(),
    event_log=None,
    method=None,
    compiled=None,
    quantiles=_DEFAULT_QUANTILES,
    selector=None,
):
    """Build and prepare a ``crn-craft.run/1.0`` request; ``simulate()`` runs it.

    The arguments are those of :func:`crn_craft.rules.ensemble`, with two
    differences: the documents are always ``crn-craft.model/1.0`` (a
    rules-model gets that schema string, a network ``crn-craft.model/0.1``
    is converted) and ``crn-craft.run/1.0``, and ``quantiles`` defaults to
    ``(0.025, 0.5, 0.975)``, so ``result.statistics`` has ``"median"`` and
    ``"quantiles"``; ``quantiles=()`` or ``None`` asks for none.
    ``selector`` chooses the event selector of ``method="static_direct"``
    (absent means ``"auto"``). ``method=None`` runs the exact tree method;
    higher-order reactions (``count``) need ``"static_direct"``,
    ``"next_reaction"``, ``"tau"`` or ``"hybrid"``.
    """
    return rules._ensemble(
        model,
        end=end,
        report_times=report_times,
        replicates=replicates,
        seed=seed,
        threads=threads,
        keep_paths=keep_paths,
        keep_instances=keep_instances,
        max_events=max_events,
        max_sampling_steps=max_sampling_steps,
        max_instances=max_instances,
        max_output_bytes=max_output_bytes,
        wall_time_s=wall_time_s,
        trackers=trackers,
        event_log=event_log,
        method=method,
        compiled=compiled,
        quantiles=quantiles,
        selector=rules._selector(selector),
        unified=True,
    )


def loop_source(model):
    """Return the generated Static Direct event loop of a model as Rust.

    As :func:`crn_craft.rules.loop_source`, for the document that
    :func:`ensemble` runs: a rules-model gets the schema string
    ``crn-craft.model/1.0`` and a network ``crn-craft.model/0.1`` is
    converted, so the loop's fingerprint is that of the model/1.0 bytes.
    """
    return _rules_loop.loop_source(rules._as_model_1_0(model, network=True))


def compile_loop(model, *, cache_dir=None, cargo="cargo", source=None, progress=True):
    """Compile the Static Direct event loop of a model; opt-in, needs Rust.

    As :func:`crn_craft.rules.compile_loop`, for the document that
    :func:`ensemble` runs: a rules-model gets the schema string
    ``crn-craft.model/1.0`` and a network ``crn-craft.model/0.1`` is
    converted. So ``ensemble(model, ..., compiled=compile_loop(model))``
    works for the same ``model`` argument. Generated loops run the linear
    scan: pass ``selector="linear"`` to :func:`ensemble` unless ``"auto"``
    picks it for the model.
    """
    return _rules_loop.compile_loop(
        rules._as_model_1_0(model, network=True),
        cache_dir=cache_dir,
        cargo=cargo,
        source=source,
        progress=progress,
    )


def prepare(model, request, *, compiled=None, quantiles=rules._DEFAULT, selector=None):
    """Prepare a request you wrote, as ``crn-craft.run/1.0``.

    As :func:`crn_craft.rules.prepare`, but the documents are always the
    1.0 contracts: a rules-run request is rewritten as ``crn-craft.run/1.0``
    with the defaults of its version written out (so it runs the same
    paths), a rules-model gets the schema ``crn-craft.model/1.0`` and a
    network ``crn-craft.model/0.1`` is converted; the request's model hash
    follows when it was the hash of the supplied model. ``quantiles``
    defaults to ``(0.025, 0.5, 0.975)`` unless the request states
    ``outputs.quantiles`` itself; ``quantiles=()`` asks for none. A
    ``crn-craft.run/1.0`` request for ``crn-craft.model/1.0`` bytes that
    needs no change keeps its exact bytes.
    """
    return rules._prepare(
        model, request, compiled, quantiles, rules._selector(selector), True
    )


def load(path):
    """Read a model file and check it; return its exact bytes.

    Every model version is read: ``crn-craft.model/1.0``, the rules-model
    versions and the network ``crn-craft.model/0.1``. An invalid model
    raises ``CrnError`` with the diagnostics of :func:`check`. Pass the
    bytes to :func:`ensemble` or :func:`prepare`, or to
    :meth:`Model.from_dict`.
    """
    data = _Path(path).read_bytes()
    report = check(data)
    if not report["valid"]:
        raise CrnError(report["diagnostics"])
    return data


def describe(kind=None):
    """The discovery document, or the JSON Schema of one contract.

    Without ``kind``: ``crn-craft.describe/1.0``, as ``crn-craft describe``
    prints it, with the commands, the 1.0 contracts, methods, selectors,
    quantiles and limits first and the older contracts that are still read
    under ``legacy``. With ``kind`` (for example ``"model-1.0"`` or
    ``"run-1.0"``), the schema ``crn-craft describe KIND`` prints;
    ``ValueError`` for an unknown kind.
    """
    if kind is None:
        return _json.loads(_native.describe_json())
    return _json.loads(_native.schema_json(str(kind)))


def check(model, request=None):
    """Check a model, and optionally a request for it, without running anything.

    ``model`` and ``request`` are bytes, text, dictionaries or a
    :class:`Model`. Returns the report ``crn-craft check MODEL [REQUEST]``
    prints, as a dictionary; read ``valid`` and ``diagnostics``:

    - a ``crn-craft.model/1.0`` model, with or without a request, and a
      rules-model with a request: ``crn-craft.check/1.0``; with a request,
      the report says whether ``run`` would accept the two documents
      (simulation and ``fsp`` requests) and adds ``run_sha256``;
    - a rules-model alone: ``crn-craft.rules-check/0.1``, as release 0.6.0;
    - a network ``crn-craft.model/0.1`` alone: ``crn-craft.check/1.0`` of
      the conversion (``CRN-R016`` and validation diagnostics if it cannot
      be converted).
    """
    if isinstance(model, _RulesModel):
        model = model.to_dict()
    return _json.loads(
        _native.check_any_json(
            _model_bytes(model), None if request is None else _model_bytes(request)
        )
    )


def load_example(name):
    """One of the packaged example models as a :class:`Model`.

    ``name`` is one of ``"enzyme-mechanism"``,
    ``"immigration-death"``, ``"michaelis-menten"``, ``"repressilator"``
    and ``"vilar-noise"``. The packaged files are the converted
    ``crn-craft.model/1.0`` documents of ``models/examples/`` in the
    repository (``crn_craft.convert`` of each); ``to_dict()`` returns that
    document until the model is edited.
    """
    from importlib.resources import files

    if name not in _EXAMPLES:
        raise ValueError(
            f"unknown example {name!r}; the examples are {', '.join(_EXAMPLES)}"
        )
    data = files("crn_craft").joinpath("models", name + ".model.json").read_bytes()
    return Model.from_dict(_json.loads(data))


def _imported(report):
    model = Model.from_dict(report["model"])
    model.import_report = {k: v for k, v in report.items() if k != "model"}
    return model


def read_sbml(path, *, convention):
    """Import an SBML file as a :class:`Model`.

    The SBML subset and the required ``convention``
    (``"stochastic_count_propensity"``) are those of
    :func:`crn_craft.importers.import_sbml`. The network is converted to
    ``crn-craft.model/1.0``; ``model.import_report`` keeps the import
    report's provenance (``source_format``, ``source_sha256``, the
    convention and parameter policy) without the model. Unsupported input
    raises ``crn_craft.importers.ImportFailure`` (a ``ValueError``).
    """
    from .importers import import_sbml

    return _imported(import_sbml(_Path(path).read_bytes(), convention=convention))


def read_bngl(path, *, convention):
    """Import an expanded BioNetGen network (``.net``) as a :class:`Model`.

    ``convention`` is the rate convention of
    :func:`crn_craft.importers.import_bng_net` (``"binomial"`` or
    ``"falling_factorial"``). The network is converted to
    ``crn-craft.model/1.0``; ``model.import_report`` keeps the import
    report's provenance (species patterns, groups, the convention) without
    the model. Unsupported input raises
    ``crn_craft.importers.ImportFailure`` (a ``ValueError``).
    """
    from .importers import import_bng_net

    return _imported(
        import_bng_net(_Path(path).read_bytes(), rate_convention=convention)
    )
