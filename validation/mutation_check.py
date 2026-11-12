"""Exercise source mutations in an isolated copy; never edits the checkout.

The unmutated copy must pass first, so a mutant counts as caught only when it
turns a passing suite into named test failures.
"""

import argparse
import contextlib
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MUTATIONS = [
    (
        "missing_binomial_divisor",
        "crates/crn-craft/src/rules/combination.rs",
        "(n - order + i) as f64 / i as f64",
        "(n - order + i) as f64 / 1.0",
    ),
    (
        "wrong_product_sign",
        "crates/crn-craft/src/rules/flat.rs",
        ".chain(products.iter().map(|&s| (s, 1)))",
        ".chain(products.iter().map(|&s| (s, -1)))",
    ),
    (
        "missing_custom_rate_read",
        "crates/crn-craft/src/rules/flat.rs",
        "rule_reads[c.rule]\n                .iter()\n                .copied()\n                .chain(",
        "std::iter::empty::<usize>()\n                .chain(",
    ),
    (
        "missing_reactant_gate",
        "crates/crn-craft/src/rules/flat.rs",
        "if reactants.iter().any(|&(s, k)| counts[s] < k) {\n            return 0.0;\n        }\n        let global = self.global(net, c, env);",
        "if reactants.iter().any(|&(s, k)| counts[s] < k.min(0)) {\n            return 0.0;\n        }\n        let global = self.global(net, c, env);",
    ),
    (
        "hill_denominator_without_one",
        "crates/crn-craft/src/rules/expr.rs",
        "self.c + self.v / (self.one + z)",
        "self.c + self.v / z",
    ),
    (
        "bound_not_divided_by_weight",
        "crates/crn-craft/src/conservation.rs",
        "let bound = total / weight;",
        "let bound = total;",
    ),
    (
        "binomial_without_symmetry",
        "crates/crn-craft/src/rules/combination.rs",
        "        let order = k.min(n - k);\n        for i in 1..=order {",
        "        let order = k;\n        for i in 1..=order {",
    ),
    (
        "cone_missing_positive_contribution",
        "crates/crn-craft/src/cone.rs",
        "let a = minus_q / divisor;",
        "let a = 0;",
    ),
    (
        "cone_without_adjacency",
        "crates/crn-craft/src/cone.rs",
        "if !adjacent {",
        "if !adjacent && false {",
    ),
    (
        "cone_bound_without_weight",
        "crates/crn-craft/src/cone.rs",
        "let bound = total / weight;",
        "let bound = total;",
    ),
    (
        "linkage_reversed_outgoing",
        "crates/crn-craft/src/linkage.rs",
        "out[edge.source].push(edge.target);",
        "out[edge.target].push(edge.source);",
    ),
    (
        "complex_identity_ignores_weight",
        "crates/crn-craft/src/linkage.rs",
        ".map(|t| (t.species, t.coefficient))",
        ".map(|t| (t.species, 1))",
    ),
    (
        "export_without_homodimer_factor",
        "crates/crn-craft/src/export.rs",
        "coefficient / 2.0",
        "coefficient",
    ),
    (
        "analysis_drops_fallback_bounds",
        "crates/crn-craft/src/analysis.rs",
        "\n        conservation.initial_class_upper_bounds.clone()",
        "\n        vec![None; view.species().len()]",
    ),
    (
        "cone_rays_not_primitive",
        "crates/crn-craft/src/cone.rs",
        "*value /= common;",
        "*value /= 1;",
    ),
    (
        "cone_initial_ray_cap_ignored",
        "crates/crn-craft/src/cone.rs",
        "if n > limits.max_rays {\n        return Err(ConeStatus::RayLimit);\n    }",
        "",
    ),
    (
        "linkage_directed_edges_only",
        "crates/crn-craft/src/linkage.rs",
        "        undirected[edge.target].push(edge.source);\n",
        "",
    ),
    (
        "ssa_wait_wrong_scale",
        "crates/crn-craft/src/rules/flat_path.rs",
        "-libm::log(1.0 - words.open()?) / total\n            };\n            let next = clock.advance(*time, wait)?;\n            let channel = scan(",
        "-libm::log(1.0 - words.open()?) * total\n            };\n            let next = clock.advance(*time, wait)?;\n            let channel = scan(",
    ),
    (
        "ssa_scan_does_not_subtract",
        "crates/crn-craft/src/rules/flat_path.rs",
        "                break;\n            }\n            target -= a;",
        "                break;\n            }\n            target -= 0.0;",
    ),
    (
        "ssa_reuses_replicate_stream",
        "crates/crn-craft/src/rules/mod.rs",
        "rng.set_stream(replicate);",
        "rng.set_stream(0);",
    ),
    (
        "ssa_save_row_not_written",
        "crates/crn-craft/src/rules/flat_path.rs",
        "out.values[row..row + n].copy_from_slice(values);",
        "let _ = values;",
    ),
    (
        "ssa_event_cap_off_by_one",
        "crates/crn-craft/src/rules/flat_path.rs",
        "events <= self.limits.max_events.saturating_sub(out.events)",
        "events < self.limits.max_events.saturating_sub(out.events)",
    ),
    (
        "ssa_biased_variance",
        "crates/crn-craft/src/rules/mod.rs",
        "n.saturating_sub(1) as f64",
        "n as f64",
    ),
    (
        "ssa_conditional_statistics",
        "crates/crn-craft/src/rules/mod.rs",
        "let complete = replicates.iter().all(|r| r.termination.complete());",
        "let complete = replicates.iter().any(|r| r.termination.complete());",
    ),
    (
        "archive_ignores_payload_digest",
        "crates/crn-craft/src/archive/v1.rs",
        "entry.sha256 != sha256(&bytes)",
        "false",
    ),
    (
        "archive_ignores_completion_digest",
        "crates/crn-craft/src/archive/v1.rs",
        'marker != format!("{}\\n", sha256(&manifest_bytes)).as_bytes()',
        "false",
    ),
    (
        "ssa_generic_poll_removed",
        "crates/crn-craft/src/rules/engine.rs",
        "if out.events % POLL_EVENTS == 0 && out.events > 0 {\n                self.poll()?;\n            }\n            let total = state.propensities",
        "if false {\n                self.poll()?;\n            }\n            let total = state.propensities",
    ),
    (
        "ssa_cancel_file_condition_inverted",
        "crates/crn-craft/src/rules/mod.rs",
        "self.cancel_file.as_ref().is_some_and(|p| p.exists())",
        "self.cancel_file.as_ref().is_some_and(|p| !p.exists())",
    ),
    (
        "ssa_extra_retained_trajectory",
        "crates/crn-craft/src/rules/mod.rs",
        "if id < keep_paths {",
        "if id <= keep_paths {",
    ),
    (
        "direct_scan_accepts_unresolved_rate",
        "crates/crn-craft/src/rules/flat.rs",
        "rates.iter().any(|&a| a > 0.0 && a < floor)",
        "rates.iter().any(|&a| a > 0.0 && a < 0.0)",
    ),
    (
        "linear_accepts_unresolved_rate",
        "crates/crn-craft/src/rules/flat_path.rs",
        "if !total.is_finite() || least < floor {\n                return Err(unresolved(total));",
        "if !total.is_finite() || least < 0.0 {\n                return Err(unresolved(total));",
    ),
    (
        "tree_accepts_unresolved_rate",
        "crates/crn-craft/src/rules/select.rs",
        "root.minimum.is_finite() && root.minimum / root.sum < RESOLUTION",
        "false",
    ),
    (
        "ssa_accepts_zero_wait",
        "crates/crn-craft/src/rules/engine.rs",
        "if next <= time {",
        "if next < time {",
    ),
    (
        "archive_ignores_byte_budget",
        "crates/crn-craft/src/archive/v1.rs",
        "if bytes.len() as u64 > output_limit(max_output_bytes) {",
        "if false {",
    ),
    (
        "archive_ignores_unstarted_record_rules",
        "crates/crn-craft/src/archive/v1.rs",
        "|| (!r.started\n",
        "|| (false && !r.started\n",
    ),
    (
        "archive_ignores_initial_state",
        "crates/crn-craft/src/archive/v1.rs",
        ".any(|row| !same_bits(row, &initial))",
        ".any(|_| false)",
    ),
    (
        "archive_accepts_state_after_stop_time",
        "crates/crn-craft/src/archive/v1.rs",
        "|| r.stop_time > end",
        "|| false",
    ),
    (
        "cached_dependency_updates_removed",
        "crates/crn-craft/src/rules/flat.rs",
        ".map(|c| union(&readers, c.changes.iter().map(|&(s, _)| s)))",
        ".map(|_| Vec::new())",
    ),
    (
        "rule_resolution_check_removed",
        "crates/crn-craft/src/rules/engine.rs",
        "let unresolved = a[r] < total * RESOLUTION\n",
        "let unresolved = false\n",
    ),
    (
        "fsp_reflects_sink",
        "crates/crn-craft/src/fsp.rs",
        "let target = if escaped {\n                usize::MAX",
        "let target = if escaped {\n                i",
    ),
    (
        "fsp_wrong_sparse_orientation",
        "crates/crn-craft/src/fsp.rs",
        "next[if j == usize::MAX { n } else { j }] += term[i] * rate;",
        "next[if j == usize::MAX { n } else { j }] += term[n] * rate;",
    ),
    (
        "fsp_tail_bound_zeroed",
        "crates/crn-craft/src/fsp.rs",
        "tail_bound += next_weight / (1.0 - ratio);",
        "tail_bound += 0.0;",
    ),
    (
        "tau_claims_exact",
        "crates/crn-craft/src/rules/mod.rs",
        '"chacha8-u53/rules-approximate-cgp-poisson/2"',
        '"chacha8-u53/rules-exact-cgp-poisson/2"',
    ),
    (
        "tau_proposal_event_cap_removed",
        "crates/crn-craft/src/rules/leap.rs",
        "self.check_events(out, events)?;",
        "let _ = events;",
    ),
    (
        "kinetics_parameter_derivative_zeroed",
        "crates/crn-craft/src/kinetics.rs",
        "Input::Parameter(j) => dp[j] += g * scale,",
        "Input::Parameter(j) => dp[j] += 0.0,",
    ),
    (
        "kinetics_combination_derivative_dropped",
        "crates/crn-craft/src/kinetics.rs",
        "dy[s] += e.value * channel.factor * (others * factors[i].1);",
        "dy[s] += 0.0;",
    ),
    (
        "dual_power_base_factor_dropped",
        "crates/crn-craft/src/rules/dual.rs",
        "y * libm::pow(x, y - 1.0)",
        "libm::pow(x, y - 1.0)",
    ),
]

# Integration tests of the unified engine and of the kept analysis, structure,
# FSP, tau and kinetics code, plus the library tests. A name is used only if
# its file exists, so the list survives the removal of the old engine's tests.
TEST_FILES = [
    "boundaries_v1",
    "archives_v1",
    "clock",
    "contracts_v1",
    "convert",
    "flat_view",
    "higher_order",
    "limits_v1",
    "quantiles",
    "selectors",
    "analysis",
    "export",
    "fsp",
    "kinetics",
    "kinetics_references",
    "structure",
    "tau",
]


def test_command(root):
    command = ["test", "--release", "--locked", "--offline", "-p", "crn-craft", "--lib"]
    for name in TEST_FILES:
        if (root / "crates/crn-craft/tests" / f"{name}.rs").exists():
            command += ["--test", name]
    return command


@contextlib.contextmanager
def locked(path):
    """Holds an exclusive advisory lock on `path` (a no-op without one)."""
    if path is None:
        yield
        return
    with open(path, "a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cargo", default="cargo")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="keep the mutants already recorded in --output and run the rest",
    )
    parser.add_argument(
        "--lock", type=Path, help="advisory lock file taken around every cargo run"
    )
    parser.add_argument("--mutant-timeout", type=int, default=900)
    args = parser.parse_args()
    rows = []
    if args.resume and args.output and args.output.exists():
        rows = json.loads(args.output.read_text())["mutations"]
    done = {row["mutation"] for row in rows}
    for name, relative, old, _ in MUTATIONS:
        if (ROOT / relative).read_text().count(old) != 1:
            raise RuntimeError(f"{name}: mutation anchor not unique before compilation")

    def report(complete):
        return {
            "schema": "crn-craft.mutation-check/0.1",
            "baseline_passed": True,
            "complete": complete,
            "mutations": rows,
            "caught": sum(row["caught_by_test_failure"] for row in rows),
            "total": len(MUTATIONS),
        }

    with tempfile.TemporaryDirectory(prefix="crn-mutations-") as folder:
        clone = Path(folder) / "project"
        shutil.copytree(
            ROOT,
            clone,
            ignore=shutil.ignore_patterns(
                ".git",
                "target",
                "__pycache__",
                ".ruff_cache",
                "iteration-02",
                "iteration-03",
                "iteration-04",
            ),
        )
        command = [args.cargo, *test_command(clone)]
        # Each mutant changes one file of one crate: incremental compilation
        # makes the rebuilds a fraction of a clean release build.
        os.environ.setdefault("CARGO_INCREMENTAL", "1")
        with locked(args.lock):
            baseline = subprocess.run(
                command,
                cwd=clone,
                capture_output=True,
                text=True,
                timeout=3600,
                check=False,
            )
        if baseline.returncode != 0:
            raise RuntimeError(
                "unmutated tests fail; no mutant can count as caught\n"
                + baseline.stdout[-2000:]
                + baseline.stderr[-2000:]
            )
        print("baseline passed", file=sys.stderr, flush=True)
        for name, relative, old, new in MUTATIONS:
            if name in done:
                continue
            path = clone / relative
            original = path.read_text()
            if original.count(old) != 1:
                raise RuntimeError(f"{name}: mutation anchor not unique")
            path.write_text(original.replace(old, new))
            try:
                try:
                    with locked(args.lock):
                        result = subprocess.run(
                            command,
                            cwd=clone,
                            capture_output=True,
                            text=True,
                            timeout=args.mutant_timeout,
                            check=False,
                        )
                    stdout, stderr, code = (
                        result.stdout,
                        result.stderr,
                        result.returncode,
                    )
                except subprocess.TimeoutExpired:
                    stdout, stderr, code = "", "", None
                failed_tests = [
                    line
                    for line in stdout.splitlines()
                    if line.startswith("test ") and line.endswith(" ... FAILED")
                ]
                caught = (
                    code == 101
                    and bool(failed_tests)
                    and "could not compile" not in stdout + stderr
                )
                print(
                    f"{len(rows) + 1}/{len(MUTATIONS)} {name}: {caught}",
                    file=sys.stderr,
                    flush=True,
                )
                rows.append(
                    {
                        "mutation": name,
                        "file": relative,
                        "caught_by_test_failure": caught,
                        "exit_code": code,
                        "failed_tests": failed_tests,
                    }
                )
                if args.output:
                    args.output.write_text(json.dumps(report(False), indent=2) + "\n")
            finally:
                path.write_text(original)
    final = report(True)
    text = json.dumps(final, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    if final["caught"] != final["total"]:
        raise RuntimeError("a mutation survived or did not produce a test failure")


if __name__ == "__main__":
    main()
