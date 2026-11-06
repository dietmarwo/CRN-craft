"""Optional compiled Static Direct event loops for rules models.

``compile_loop`` asks the extension for the model's generated event loop,
writes a small Cargo project around it and builds it with the user's Rust
toolchain. The project depends on ``crn-craft-rules`` at the installed
version, from crates.io or from a local checkout, and builds a standalone
program; Python runs it as a subprocess, so no Rust object crosses a
library boundary. Builds are cached, and every handle is checked against
the interpreter before it is returned.
"""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import time

import numpy as np

from . import _crn_error, _native, model_bytes

SOURCE_ENV = "CRN_CRAFT_SOURCE"
CACHE_ENV = "CRN_CRAFT_CACHE_DIR"
HEADER = b"crn-craft-rules-loop/1\n"
# Exit status of the program for a rejected request (diagnostics on stdout).
_REJECTED = 2
_FINGERPRINT = re.compile(r'fn fingerprint\(&self\) -> &str \{ "([0-9a-f]{64})" \}')
# Settings that change the compiled program; they are part of the cache key.
_BUILD_ENV = (
    "RUSTFLAGS",
    "CARGO_ENCODED_RUSTFLAGS",
    "CARGO_BUILD_RUSTFLAGS",
    "CARGO_BUILD_TARGET",
)
_EXECUTABLE = "crn-craft-loop" + (".exe" if os.name == "nt" else "")


class CompiledLoopError(RuntimeError):
    """A compiled event loop could not be built, or it disagreed with the
    interpreter.

    ``reason`` is ``"cargo_missing"``, ``"build_failed"`` or ``"mismatch"``;
    ``log`` is the path of the build log, if a build started.
    """

    def __init__(self, message, reason, log=None):
        super().__init__(message)
        self.reason = reason
        self.log = log


def _document(model):
    from .builder import Model

    if isinstance(model, Model):
        model = model.to_dict()
    return model_bytes(model)


def _generator_request(data):
    """A minimal Static Direct request; the loop does not depend on it."""
    return model_bytes(
        {
            "schema": "crn-craft.rules-run/0.2",
            "model": {
                "path": "model.rules.json",
                "sha256": hashlib.sha256(data).hexdigest(),
            },
            "time": {"end": 1.0, "report_times": [0.0, 1.0]},
            "ensemble": {"replicates": 1, "root_seed": 0},
            "outputs": {"statistics": [], "keep_paths": 0, "keep_instances": 0},
            "budget": {
                "threads": 1,
                "max_events_per_replicate": 1,
                "max_sampling_steps_per_replicate": 1,
                "max_instances": 100_000,
                "max_output_bytes": 1 << 20,
            },
            "method": {"kind": "static_direct"},
        }
    )


def _prepare(data, request):
    try:
        return _native.prepare_rules(data, request)
    except _native.NativeError as error:
        raise _crn_error(error) from None


def _source(native):
    try:
        return native.loop_source()
    except _native.NativeError as error:
        raise _crn_error(error) from None


def fingerprint_of(source):
    """The fingerprint compiled into a generated loop's source."""
    match = _FINGERPRINT.search(source)
    if match is None:
        raise RuntimeError("generated rules loop without a fingerprint")
    return match.group(1)


def loop_source(model):
    """Return the generated Static Direct event loop of a model as Rust.

    Needs no Rust toolchain. ``model`` is bytes, text, a dictionary or a
    ``rules.Model``. The generator supports closed models with constant or
    vacancy rates and sum observables, up to 1 MiB of source; it refuses
    other models with ``CrnError`` ``CRN-R020`` at ``/method``: expression
    rates, ``propensity`` rules, max or distinct observables, densely
    coupled networks, and models that only another method runs (continuous
    variables, ATP blocks, instances that change identity).
    """
    data = _document(model)
    return _source(_prepare(data, _generator_request(data)))


def _cache_root(cache_dir):
    if cache_dir is not None:
        return Path(cache_dir).expanduser().resolve()
    if os.environ.get(CACHE_ENV):
        return Path(os.environ[CACHE_ENV]).expanduser().resolve()
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base).resolve() / "crn-craft" / "rules-loops"


def _checkout(source):
    """The local CRN-craft checkout to build against, or None for crates.io."""
    if source is None:
        source = os.environ.get(SOURCE_ENV) or None
    if source is None:
        return None
    root = Path(source).expanduser().resolve()
    if not (root / "crates" / "crn-craft-rules" / "Cargo.toml").is_file():
        raise ValueError(
            f"source={str(root)!r} is not a CRN-craft checkout: "
            "crates/crn-craft-rules/Cargo.toml is missing"
        )
    return root


def _tree_digest(root):
    """SHA-256 over the checkout files that make up crn-craft-rules."""
    digest = hashlib.sha256()
    files = [root / "Cargo.toml", root / "Cargo.lock"]
    for crate in ("crn-craft", "crn-craft-rules"):
        folder = root / "crates" / crate
        files += [folder / "Cargo.toml", folder / "build.rs"]
        files += sorted((folder / "src").rglob("*"))
    for path in files:
        if path.is_file():
            digest.update(path.relative_to(root).as_posix().encode() + b"\0")
            digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def _rustc(cargo, cwd):
    """``rustc -vV`` of the toolchain next to cargo, else of RUSTC or PATH."""
    sibling = Path(cargo).with_name("rustc" + Path(cargo).suffix)
    rustc = os.environ.get("RUSTC") or (
        str(sibling) if sibling.is_file() else shutil.which("rustc")
    )
    if rustc is None:
        raise CompiledLoopError(
            "rustc was not found next to cargo or on PATH; compile_loop needs a "
            "Rust toolchain (https://rustup.rs)",
            "cargo_missing",
        )
    done = subprocess.run(
        [rustc, "-vV"], cwd=cwd, capture_output=True, text=True, check=False
    )
    if done.returncode != 0:
        raise CompiledLoopError(
            f"{rustc} -vV failed: {done.stderr.strip()}", "cargo_missing"
        )
    return done.stdout


def _manifest(name, version, rust_version, checkout):
    rules = {"version": f"={version}"}
    if checkout is not None:
        rules["path"] = (checkout / "crates" / "crn-craft-rules").as_posix()
    # JSON strings are valid TOML basic strings.
    spec = ", ".join(
        f"{key} = {json.dumps(value, ensure_ascii=False)}"
        for key, value in rules.items()
    )
    return (
        "# Generated by crn_craft.rules.compile_loop.\n"
        "[package]\n"
        f'name = "{name}"\n'
        'version = "0.0.0"\n'
        'edition = "2024"\n'
        f'rust-version = "{rust_version}"\n'
        "publish = false\n\n"
        "[[bin]]\n"
        f'name = "{name}"\n'
        'path = "src/main.rs"\n\n'
        "[dependencies]\n"
        f"crn-craft-rules = {{ {spec} }}\n"
        'serde_json = "1"\n\n'
        "[workspace]\n"
    )


def _write(path, data):
    """Write bytes or text atomically, leaving equal files untouched."""
    data = data.encode() if isinstance(data, str) else data
    if path.is_file() and path.read_bytes() == data:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def _build(entry, root, cargo, name, progress, checkout):
    """Run cargo; return the built executable and the seconds it took."""
    log = entry / "build.log"
    messages = entry / "build-messages.json"
    env = dict(os.environ)
    env.setdefault("CARGO_TARGET_DIR", str(root / "target"))
    command = [
        cargo,
        "build",
        "--release",
        "--manifest-path",
        str(entry / "Cargo.toml"),
        "--message-format=json-render-diagnostics",
    ]
    if progress:
        print(
            "crn_craft.rules.compile_loop: building the event loop with cargo. "
            "The first build downloads and compiles crn-craft-rules and can take "
            f"several minutes; later models reuse {env['CARGO_TARGET_DIR']}. "
            f"Log: {log}",
            file=sys.stderr,
            flush=True,
        )
    start = time.perf_counter()
    with open(log, "wb") as text, open(messages, "wb") as out:
        process = subprocess.Popen(
            command, cwd=entry, env=env, stdout=out, stderr=subprocess.PIPE
        )
        try:
            for line in process.stderr:
                text.write(line)
                text.flush()
                if progress:
                    sys.stderr.write("  " + line.decode("utf-8", "replace"))
                    sys.stderr.flush()
            status = process.wait()
        except BaseException:
            process.kill()
            process.wait()
            raise
    seconds = time.perf_counter() - start
    if status != 0:
        lines = log.read_text("utf-8", "replace").strip().splitlines()
        errors = [line for line in lines if line.startswith("error")] or lines
        hint = (
            ""
            if checkout is not None
            else f" If crn-craft-rules {_native.__version__} is not on crates.io, "
            f"pass source= a CRN-craft checkout of that version or set {SOURCE_ENV}."
        )
        raise CompiledLoopError(
            f"cargo build failed with status {status}: "
            f"{errors[0] if errors else 'no output'}. Log: {log}.{hint}",
            "build_failed",
            log,
        )
    built = None
    for line in messages.read_text("utf-8").splitlines():
        try:
            message = json.loads(line)
        except ValueError:
            continue
        if (
            message.get("reason") == "compiler-artifact"
            and message.get("executable")
            and message.get("target", {}).get("name") == name
        ):
            built = Path(message["executable"])
    if built is None or not built.is_file():
        raise CompiledLoopError(
            f"cargo reported no executable for {name}. Log: {log}",
            "build_failed",
            log,
        )
    return built, seconds


class CompiledLoop:
    """A compiled Static Direct event loop for one model and parameter set.

    Returned by :func:`compile_loop`. Pass it as ``compiled=`` to
    ``rules.ensemble``, ``rules.prepare`` or ``PreparedEnsemble.simulate``
    to run Static Direct requests for exactly this model with it.

    - ``fingerprint``: the generated loop's fingerprint, which covers the
      emitted source and the model bytes;
    - ``model_sha256``: SHA-256 of the model bytes it was compiled for;
    - ``version``: the ``crn-craft-rules`` version it was built against;
    - ``path``: the program; ``directory`` its cache entry with the Cargo
      project and ``log`` the build log;
    - ``cached``: whether a cached build was reused; ``build_seconds`` the
      cargo time of a new build, else None.
    """

    def __init__(self, *, fingerprint, model_sha256, version, directory, cached):
        self.fingerprint = fingerprint
        self.model_sha256 = model_sha256
        self.version = version
        self.directory = Path(directory)
        self.path = self.directory / _EXECUTABLE
        self.log = self.directory / "build.log"
        self.cached = cached
        self.build_seconds = None

    def __repr__(self):
        return (
            f"CompiledLoop(fingerprint={self.fingerprint[:16]}..., "
            f"version={self.version!r}, path={str(self.path)!r})"
        )

    def execute(self, request, *, cancel_file=None):
        """Run a rules-run request for this loop's model in its program.

        Returns the metadata text and the record columns, as the extension
        does. A rejected request raises ``CrnError``.
        """
        command = [str(self.path)]
        if cancel_file is not None:
            command += ["--cancel-file", os.fspath(cancel_file)]
        done = subprocess.run(command, input=request, capture_output=True, check=False)
        if done.returncode == _REJECTED:
            raise _crn_error(done.stdout.decode("utf-8", "replace")) from None
        if done.returncode != 0 or not done.stdout.startswith(HEADER):
            stderr = done.stderr.decode("utf-8", "replace").strip()
            raise RuntimeError(
                f"compiled rules loop {self.path} failed with status "
                f"{done.returncode}: {stderr or 'no message'}"
            )
        return _frames(done.stdout)

    def verify(self, run):
        """Run a prepared ensemble interpreted and compiled, and compare.

        Raises ``CompiledLoopError`` (reason ``"mismatch"``) unless the
        metadata and every record column are byte-identical. Returns the
        result.
        """
        from .rules import Result

        request = run._compiled_request(self)
        interpreted = run._native.simulate(None, None)
        compiled = self.execute(request)
        _compare(interpreted, compiled, self.directory / "mismatch")
        return Result(json.loads(compiled[0]), compiled[1])


def _frames(data):
    view = memoryview(data)
    position, frames = len(HEADER), {}
    while position < len(view):
        (size,) = struct.unpack_from("<I", view, position)
        name = bytes(view[position + 4 : position + 4 + size]).decode()
        position += 4 + size
        (size,) = struct.unpack_from("<Q", view, position)
        frames[name] = bytes(view[position + 8 : position + 8 + size])
        position += 8 + size
    return frames.pop("metadata").decode(), frames


def _compare(interpreted, compiled, folder):
    if interpreted[0] == compiled[0] and interpreted[1] == compiled[1]:
        return
    folder.mkdir(parents=True, exist_ok=True)
    for label, (metadata, columns) in (
        ("interpreted", interpreted),
        ("compiled", compiled),
    ):
        (folder / f"{label}.json").write_text(metadata, "utf-8")
        for key, value in columns.items():
            (folder / f"{label}.{key}.bin").write_bytes(value)
    differ = [key for key in interpreted[1] if interpreted[1][key] != compiled[1][key]]
    if interpreted[0] != compiled[0]:
        differ.insert(0, "metadata")
    raise CompiledLoopError(
        f"the compiled loop's result differs from the interpreter's in "
        f"{', '.join(differ)}; both are saved in {folder}",
        "mismatch",
    )


def _self_check(handle, data):
    """Short requests through the program and the interpreter; identical bytes.

    The requests are built by ``rules.ensemble``, so they use the request
    version it writes; a program built against another engine version
    rejects them or differs.
    """
    from .rules import ensemble

    def run(**settings):
        return ensemble(data, replicates=8, seed=1, threads=2, **settings)

    # A horizon most replicates reach in about 2000 events.
    probe = run(end=1e300, report_times=[0.0, 1e300], max_events=2000)
    _, columns = probe._native.simulate(None, None)
    stop = np.frombuffer(columns["stop_time"], dtype="<f8")
    limited = stop[np.frombuffer(columns["termination"], dtype="u1") == 2]
    end = float(limited.min()) if len(limited) else 1.0
    times = [0.0, end / 4, end / 2, 3 * end / 4, end]
    if any(b <= a for a, b in zip(times, times[1:], strict=False)):
        times = [0.0, end]
    names = sorted(json.loads(data).get("observables", {}))
    trackers = []
    if names:
        trackers = [
            {"observable": names[0], "keep_episodes": 2},
            {"observable": names[-1], "first_passage": f"{names[-1]} > 0"},
        ]
    for method in (
        "static_direct",
        {"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False},
        {"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": True},
    ):
        check = run(
            end=end,
            report_times=times,
            max_events=4000,
            method=method,
            trackers=trackers,
            event_log={"max_entries": 16},
            keep_paths=2,
            keep_instances=1,
        )
        interpreted = check._native.simulate(None, None)
        try:
            compiled = handle.execute(check._documents[1])
        except (ValueError, RuntimeError) as error:
            raise CompiledLoopError(
                f"the compiled loop failed its self-check: {error}",
                "mismatch",
                handle.log,
            ) from error
        _compare(interpreted, compiled, handle.directory / "self-check")


def compile_loop(model, *, cache_dir=None, cargo="cargo", source=None, progress=True):
    """Compile the Static Direct event loop of a model; opt-in, needs Rust.

    ``model`` is bytes, text, a dictionary or a ``rules.Model``; the loop
    serves Static Direct requests for exactly these model bytes, with any
    run settings. Parameters are compiled in, so another parameter set
    needs another loop.

    The build writes a Cargo project that depends on ``crn-craft-rules`` at
    the installed version (``crn_craft._native.__version__``) from
    crates.io, or from the checkout ``source`` (default: the
    ``CRN_CRAFT_SOURCE`` environment variable), and runs ``cargo build
    --release``. It needs Rust at least at the crate's MSRV, and network
    access to crates.io unless cargo is offline with the crates cached.
    The first build can take several minutes; ``progress`` prints cargo's
    status lines and the log path to stderr.

    Builds are cached under ``cache_dir`` (default: ``CRN_CRAFT_CACHE_DIR``,
    else the user cache directory, for example
    ``~/.cache/crn-craft/rules-loops``), keyed by the loop fingerprint, the
    installed version, ``rustc -vV``, the target and Rust flags, and the
    checkout's contents. Builds share one Cargo target directory there
    (``CARGO_TARGET_DIR`` if set), so later models compile in seconds.

    Before returning, short Static Direct requests run through the program
    and the interpreter, with both waiting-time policies, trackers and an
    event log; the handle is returned only if their result bytes are
    identical.

    Raises ``CrnError`` (``CRN-R020``) for models the generator refuses
    (see :func:`loop_source`) and ``CompiledLoopError`` if cargo is missing
    (reason ``"cargo_missing"``), the build fails (``"build_failed"``, with
    ``log``) or the self-check differs (``"mismatch"``).
    """
    data = _document(model)
    source_code = _source(_prepare(data, _generator_request(data)))
    fingerprint = fingerprint_of(source_code)
    version = _native.__version__
    rust_version, driver = _native.rules_loop_driver()
    found = shutil.which(os.fspath(cargo))
    if found is None:
        raise CompiledLoopError(
            f"cargo was not found ({os.fspath(cargo)!r}); compile_loop needs a Rust "
            f"toolchain {rust_version} or newer (https://rustup.rs)",
            "cargo_missing",
        )
    checkout = _checkout(source)
    root = _cache_root(cache_dir)
    root.mkdir(parents=True, exist_ok=True)
    rustc = _rustc(found, root)
    host = re.search(r"^host: (\S+)$", rustc, re.MULTILINE)
    key_fields = {
        "fingerprint": fingerprint,
        "version": version,
        "rustc": rustc,
        "target": os.environ.get("CARGO_BUILD_TARGET")
        or (host.group(1) if host else ""),
        "environment": {name: os.environ.get(name) for name in _BUILD_ENV},
        "driver": hashlib.sha256(json.dumps(driver).encode()).hexdigest(),
        "source": "crates.io"
        if checkout is None
        else {"path": str(checkout), "sha256": _tree_digest(checkout)},
    }
    key = hashlib.sha256(json.dumps(key_fields, sort_keys=True).encode()).hexdigest()
    entry = root / key[:32]
    record = entry / "loop.json"
    handle = CompiledLoop(
        fingerprint=fingerprint,
        model_sha256=hashlib.sha256(data).hexdigest(),
        version=version,
        directory=entry,
        cached=False,
    )
    try:
        handle.cached = (
            json.loads(record.read_text("utf-8"))["key"] == key
            and handle.path.is_file()
        )
    except (OSError, ValueError, KeyError):
        handle.cached = False
    if not handle.cached:
        name = f"crn-craft-loop-{fingerprint[:16]}"
        entry.mkdir(parents=True, exist_ok=True)
        _write(entry / "Cargo.toml", _manifest(name, version, rust_version, checkout))
        for path, text in driver:
            _write(entry / path, text)
        _write(entry / "src" / "generated_loop.rs", source_code)
        _write(entry / "model.rules.json", data)
        if checkout is not None and (checkout / "Cargo.lock").is_file():
            # Start from the checkout's pinned versions; cargo prunes the rest.
            _write(entry / "Cargo.lock", (checkout / "Cargo.lock").read_bytes())
        built, handle.build_seconds = _build(
            entry, root, found, name, progress, checkout
        )
        temporary = handle.path.with_name(f"{_EXECUTABLE}.{os.getpid()}.tmp")
        shutil.copy2(built, temporary)
        os.replace(temporary, handle.path)
        _write(record, json.dumps({"key": key, **key_fields}, indent=2) + "\n")
    try:
        _self_check(handle, data)
    except CompiledLoopError:
        record.unlink(missing_ok=True)
        raise
    return handle
