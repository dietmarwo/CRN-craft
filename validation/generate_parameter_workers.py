"""Generate CRN propensity modules and rebop macro modules from frozen JSON."""

import json
import math
import time
import crn_craft as crn
from parameter_optimization_models import FOLDER, guard_fixture


def rebop_source(model):
    species = [s["name"] for s in model["species"]]
    names = sorted(model["parameters"])
    lines = [
        "use rebop::define_system;",
        "define_system! {",
        " ".join(names) + ";",
        "Model { " + ", ".join(species) + " }",
    ]
    for r in model["reactions"]:

        def side(d):
            return " + ".join(
                (str(n) + " " if n != 1 else "") + s for s, n in d.items()
            )

        coefficient = r["rate"]["coefficient"]
        coefficient = coefficient if isinstance(coefficient, str) else repr(coefficient)
        divisor = math.prod(math.factorial(n) for n in r["reactants"].values())
        rate = coefficient if divisor == 1 else f"{coefficient} / {divisor}.0"
        lines.append(
            f"{r['name']}: {side(r['reactants'])} => {side(r['products'])} @ {rate}"
        )
    lines.extend(
        [
            "}",
            "pub fn run(values: &[f64], times: &[f64], paths: u64, seeds: &[u64]) -> Result<serde_json::Value, String> {",
            "let mut prototype = Model::new();",
        ]
    )
    for i, name in enumerate(names):
        lines.append(f"prototype.{name} = values[{i}];")
    lines.append(f"let mut moments = super::Moments::new(times.len(), {len(species)});")
    lines.append(
        "for id in 0..paths { let mut p = prototype.clone(); p.seed(seeds[id as usize]);"
    )
    for s in model["species"]:
        lines.append(f"p.{s['name']} = {s['initial_count']};")
    lines.append("let mut states = Vec::new(); for &t in times { p.advance_until(t);")
    lines.append(
        "let state = vec![" + ", ".join("p." + s + " as f64" for s in species) + "];"
    )
    lines.append(
        'if state.iter().any(|x| *x < 0.) { return Err("negative rebop population".into()); } states.push(state); }'
    )
    lines.append("moments.add(&states); } Ok(moments.finish()) }")
    return "\n".join(lines) + "\n"


def main():
    report = {}
    for row in json.loads((FOLDER / "fixtures.json").read_text())["rows"]:
        case = row["case"]
        data = (FOLDER / "models" / f"{case}.model.json").read_bytes()
        start = time.perf_counter()
        source = crn.mass_action_kernel_source(data)
        (FOLDER / "worker/src" / f"crn_{case}.rs").write_text(source)
        report[case + "/crn"] = time.perf_counter() - start
        start = time.perf_counter()
        source = rebop_source(json.loads(data))
        (FOLDER / "worker/src" / f"rebop_{case}.rs").write_text(source)
        report[case + "/rebop"] = time.perf_counter() - start
    data = json.dumps(guard_fixture()).encode()
    (FOLDER / "worker/src/crn_guard.rs").write_text(crn.mass_action_kernel_source(data))
    (FOLDER / "evidence/generation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
