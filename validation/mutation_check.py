"""Exercise source mutations in an isolated copy; never edits the checkout.

The unmutated copy must pass first, so a mutant counts as caught only when it
turns a passing suite into named test failures.
"""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MUTATIONS = [
    (
        "missing_binomial_divisor",
        "crates/crn-craft/src/rates.rs",
        "(count - order + i) as f64 / i as f64",
        "(count - order + i) as f64 / 1.0",
    ),
    (
        "wrong_product_sign",
        "crates/crn-craft/src/network.rs",
        "*changes.entry(term.species).or_insert(0) += term.coefficient;",
        "*changes.entry(term.species).or_insert(0) -= term.coefficient;",
    ),
    (
        "missing_custom_rate_read",
        "crates/crn-craft/src/network.rs",
        "reads.insert(species);",
        "let _ = species;",
    ),
    (
        "missing_reactant_gate",
        "crates/crn-craft/src/rates.rs",
        "state[term.species] < term.coefficient",
        "state[term.species] < 0",
    ),
    (
        "hill_without_log1p_near_threshold",
        "crates/crn-craft/src/rates.rs",
        "relative.abs() <= 0.5",
        "relative.abs() < 0.0",
    ),
    (
        "bound_not_divided_by_weight",
        "crates/crn-craft/src/conservation.rs",
        "let bound = total / weight;",
        "let bound = total;",
    ),
    (
        "binomial_without_symmetry",
        "crates/crn-craft/src/rates.rs",
        "let order = term.coefficient.min(count - term.coefficient);",
        "let order = term.coefficient;",
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
        "\n            conservation.initial_class_upper_bounds.clone()",
        "\n            vec![None; network.species().len()]",
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
        "crates/crn-craft/src/sampling.rs",
        "-libm::log(uniform) / total",
        "-libm::log(uniform) * total",
    ),
    (
        "ssa_selection_wrong_boundary",
        "crates/crn-craft/src/sampling.rs",
        "threshold < cumulative && rate > 0.0",
        "threshold <= cumulative && rate > 0.0",
    ),
    (
        "ssa_reuses_replicate_stream",
        "crates/crn-craft/src/sampling.rs",
        "rng.set_stream(replicate);",
        "rng.set_stream(0);",
    ),
    (
        "ssa_save_before_equal_event",
        "crates/crn-craft/src/simulation.rs",
        "run.time.save_times[save] < next",
        "run.time.save_times[save] <= next",
    ),
    (
        "ssa_event_cap_off_by_one",
        "crates/crn-craft/src/simulation.rs",
        "record.events >= run.budget.max_events_per_replicate.unwrap()",
        "record.events > run.budget.max_events_per_replicate.unwrap()",
    ),
    (
        "ssa_biased_variance",
        "crates/crn-craft/src/simulation.rs",
        "(completed - 1).max(1)",
        "completed.max(1)",
    ),
    (
        "ssa_conditional_statistics",
        "crates/crn-craft/src/simulation.rs",
        "completed == ensemble.replicates && !run.outputs.statistics.is_empty()",
        "completed > 0 && !run.outputs.statistics.is_empty()",
    ),
    (
        "archive_ignores_payload_digest",
        "crates/crn-craft/src/archive.rs",
        "entry.sha256 != sha256(&bytes)",
        "false",
    ),
    (
        "archive_ignores_completion_digest",
        "crates/crn-craft/src/archive.rs",
        'marker != format!("{}\\n", sha256(&manifest_bytes)).as_bytes()',
        "false",
    ),
    (
        "ssa_poller_skips_first_check",
        "crates/crn-craft/src/simulation.rs",
        "            countdown: 0,\n",
        "            countdown: CLOCK_POLL_STRIDE,\n",
    ),
    (
        "ssa_cancel_file_interval_inverted",
        "crates/crn-craft/src/simulation.rs",
        "now.duration_since(last) >= CANCEL_FILE_INTERVAL",
        "now.duration_since(last) < CANCEL_FILE_INTERVAL",
    ),
    (
        "ssa_extra_retained_trajectory",
        "crates/crn-craft/src/simulation.rs",
        "if id < run.outputs.keep_trajectories {",
        "if id <= run.outputs.keep_trajectories {",
    ),
    (
        "ssa_accepts_lost_positive_increment",
        "crates/crn-craft/src/sampling.rs",
        "|| (rate > 0.0 && next <= total)",
        "|| false",
    ),
    (
        "ssa_accepts_zero_wait",
        "crates/crn-craft/src/sampling.rs",
        "if !wait.is_finite() || wait <= 0.0 {",
        "if !wait.is_finite() {",
    ),
    (
        "archive_ignores_byte_budget",
        "crates/crn-craft/src/archive.rs",
        "if total_bytes > run.budget.max_output_bytes\n",
        "if false\n",
    ),
    (
        "archive_ignores_unstarted_record_rules",
        "crates/crn-craft/src/archive.rs",
        "            || (!record.started\n",
        "            || (false && !record.started\n",
    ),
    (
        "archive_ignores_initial_state",
        "crates/crn-craft/src/archive.rs",
        ".is_some_and(|counts| counts != prepared.network.initial_counts())",
        ".is_some_and(|_| false)",
    ),
    (
        "archive_accepts_state_after_stop_time",
        "crates/crn-craft/src/archive.rs",
        "&& run.time.save_times[record.saved_points - 1] > record.stop_time",
        "&& false",
    ),
]
MUTATIONS.extend(
    [
        (
            "cached_dependency_updates_removed",
            "crates/crn-craft/src/simulation.rs",
            "rates[reaction] = network.propensity_in_valid_state(reaction, &state)?;",
            "let _ = network.propensity_in_valid_state(reaction, &state)?;",
        ),
        (
            "minimum_resolution_check_removed",
            "crates/crn-craft/src/sampling.rs",
            "minimum.is_finite() && minimum / total < UNIFORM_STEP",
            "false",
        ),
        (
            "fsp_reflects_sink",
            "crates/crn-craft/src/fsp.rs",
            "let target = if escaped {\n                    usize::MAX",
            "let target = if escaped {\n                    i",
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
            "crates/crn-craft/src/tau.rs",
            '"approximate_ensemble"',
            '"exact_ensemble"',
        ),
        (
            "tau_proposal_event_cap_removed",
            "crates/crn-craft/src/tau.rs",
            "if events > run.budget.max_events_per_replicate.unwrap() - record.events {",
            "if false {",
        ),
        (
            "kinetics_parameter_derivative_zeroed",
            "crates/crn-craft/src/kinetics.rs",
            "dp[j] += product;",
            "dp[j] += 0.0;",
        ),
    ]
)

TEST_COMMAND = [
    "test",
    "--release",
    "--locked",
    "--offline",
    "-p",
    "crn-craft",
    "--lib",
    "--test",
    "ssa",
    "--test",
    "archives",
    "--test",
    "network",
    "--test",
    "structure",
    "--test",
    "fsp",
    "--test",
    "tau",
    "--test",
    "kinetics",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cargo", default="cargo")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = []
    for name, relative, old, _ in MUTATIONS:
        if (ROOT / relative).read_text().count(old) != 1:
            raise RuntimeError(f"{name}: mutation anchor not unique before compilation")
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
        baseline = subprocess.run(
            [args.cargo, *TEST_COMMAND],
            cwd=clone,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        if baseline.returncode != 0:
            raise RuntimeError(
                "unmutated tests fail; no mutant can count as caught\n"
                + baseline.stdout[-2000:]
                + baseline.stderr[-2000:]
            )
        for name, relative, old, new in MUTATIONS:
            path = clone / relative
            original = path.read_text()
            if original.count(old) != 1:
                raise RuntimeError(f"{name}: mutation anchor not unique")
            path.write_text(original.replace(old, new))
            try:
                result = subprocess.run(
                    [args.cargo, *TEST_COMMAND],
                    cwd=clone,
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=False,
                )
                combined = result.stdout + result.stderr
                failed_tests = [
                    line
                    for line in result.stdout.splitlines()
                    if line.startswith("test ") and line.endswith(" ... FAILED")
                ]
                caught = (
                    result.returncode == 101
                    and bool(failed_tests)
                    and "could not compile" not in combined
                )
                print(
                    f"{len(rows) + 1}/{len(MUTATIONS)} {name}: {caught}",
                    file=__import__("sys").stderr,
                    flush=True,
                )
                rows.append(
                    {
                        "mutation": name,
                        "file": relative,
                        "caught_by_test_failure": caught,
                        "exit_code": result.returncode,
                        "failed_tests": failed_tests,
                    }
                )
            finally:
                path.write_text(original)
    report = {
        "schema": "crn-craft.mutation-check/0.1",
        "baseline_passed": True,
        "mutations": rows,
        "caught": sum(row["caught_by_test_failure"] for row in rows),
        "total": len(MUTATIONS),
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    if report["caught"] != report["total"]:
        raise RuntimeError("a mutation survived or did not produce a test failure")


if __name__ == "__main__":
    main()
