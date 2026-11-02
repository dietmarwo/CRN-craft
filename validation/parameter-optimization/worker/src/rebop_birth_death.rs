use rebop::define_system;
define_system! {
birth death;
Model { X }
birth:  => X @ birth
death: X =>  @ death
}
pub fn run(
    values: &[f64],
    times: &[f64],
    paths: u64,
    seeds: &[u64],
) -> Result<serde_json::Value, String> {
    let mut prototype = Model::new();
    prototype.birth = values[0];
    prototype.death = values[1];
    let mut moments = super::Moments::new(times.len(), 1);
    for id in 0..paths {
        let mut p = prototype.clone();
        p.seed(seeds[id as usize]);
        p.X = 0;
        let mut states = Vec::new();
        for &t in times {
            p.advance_until(t);
            let state = vec![p.X as f64];
            if state.iter().any(|x| *x < 0.) {
                return Err("negative rebop population".into());
            }
            states.push(state);
        }
        moments.add(&states);
    }
    Ok(moments.finish())
}
