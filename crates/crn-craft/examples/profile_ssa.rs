//! Diagnostic microbenchmarks; public rate APIs include their validation cost.
use crn_craft::{CompiledNetwork, parse_model, sampling, simulation};
use rand_core::RngCore;
use serde_json::json;
use std::{hint::black_box, time::Instant};

fn measure(mut f: impl FnMut(), count: usize) -> f64 {
    let start = Instant::now();
    for _ in 0..count {
        f();
    }
    start.elapsed().as_secs_f64() / count as f64
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let model = std::fs::read(&args[1]).unwrap();
    let run = std::fs::read(&args[2]).unwrap();
    let prepared = simulation::prepare(&model, &run).unwrap();
    let network = CompiledNetwork::compile(&parse_model(&model).unwrap()).unwrap();
    let state: Vec<_> = network
        .initial_counts()
        .iter()
        .map(|&x| x.max(10))
        .collect();
    let mut rng = sampling::stream(4242, 0);
    let word = measure(
        || {
            black_box(rng.next_u64());
        },
        1_000_000,
    );
    let mut rng = sampling::stream(4242, 0);
    let inverse = measure(
        || {
            black_box(sampling::waiting_time(100.0, sampling::uniform_open(&mut rng)).unwrap());
        },
        1_000_000,
    );
    let rates = measure(
        || {
            black_box(network.propensities(&state).unwrap());
        },
        10_000,
    );
    let preparation = measure(
        || {
            black_box(simulation::prepare(&model, &run).unwrap());
        },
        100,
    );
    let control = simulation::Control::default();
    let mut direct = Vec::new();
    let mut cached = Vec::new();
    for _ in 0..5 {
        direct.push(measure(
            || {
                black_box(simulation::execute(&prepared, &control).unwrap());
            },
            1,
        ));
        cached.push(measure(
            || {
                black_box(simulation::execute_cached(&prepared, &control).unwrap());
            },
            1,
        ));
    }
    let result = simulation::execute_cached(&prepared, &control).unwrap();
    let encoding = measure(
        || {
            black_box(serde_json::to_vec(&result.result).unwrap());
        },
        100,
    );
    println!(
        "{}",
        json!({"rng_word_ns":word*1e9,"inverse_wait_and_rng_ns":inverse*1e9,
        "public_propensity_batch_ns":rates*1e9,"prepare_seconds":preparation,
        "result_encoding_seconds":encoding,"direct_seconds":direct,"cached_seconds":cached,
        "events":result.result.events,"species":state.len(),"reactions":network.reactions().len()})
    );
}
