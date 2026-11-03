//! Emit the normal rules result plus phase timings on stderr.
use crn_craft_rules::{Control, execute, prepare};
use std::io::{BufWriter, Write};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let [model, request] = &std::env::args_os().skip(1).collect::<Vec<_>>()[..] else {
        return Err("usage: benchmark_rules MODEL REQUEST".into());
    };
    let start = Instant::now();
    let model = std::fs::read(model)?;
    let request = std::fs::read(request)?;
    let prepared = prepare(&model, &request).map_err(|e| format!("{e:?}"))?;
    let preparation = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let result = execute(&prepared, &Control::default()).map_err(|e| format!("{e:?}"))?;
    let simulation = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let mut out = BufWriter::new(std::io::stdout().lock());
    serde_json::to_writer(&mut out, &result)?;
    out.flush()?;
    eprintln!(
        "{{\"preparation\":{preparation},\"simulation\":{simulation},\"output\":{}}}",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
