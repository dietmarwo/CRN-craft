//! Emit an optional Rust propensity module: generate_kernel MODEL.json OUTPUT.rs
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: generate_kernel MODEL.json OUTPUT.rs".into());
    }
    let model = std::fs::read(&args[1])?;
    let source = crn_craft::codegen::mass_action_source(&model)
        .map_err(|e| serde_json::to_string(&e).unwrap())?;
    std::fs::write(&args[2], source)?;
    Ok(())
}
