use crn_craft_rules::{generate_rules_loop, prepare};

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let [model, run] = args.as_slice() else {
        return Err("usage: generate_rules_loop MODEL REQUEST".into());
    };
    let read = |path| std::fs::read(path).map_err(|e| e.to_string());
    let prepared = prepare(&read(model)?, &read(run)?).map_err(|e| format!("{e:?}"))?;
    print!(
        "{}",
        generate_rules_loop(&prepared).map_err(|e| format!("{e:?}"))?
    );
    Ok(())
}
