use rebop::define_system;
define_system! {
bind unbind;
Model { A, B, C }
bind: A + B => C @ bind
unbind: C => A + B @ unbind
}
pub fn run(
    values: &[f64],
    times: &[f64],
    paths: u64,
    seeds: &[u64],
) -> Result<serde_json::Value, String> {
    let mut prototype = Model::new();
    prototype.bind = values[0];
    prototype.unbind = values[1];
    let mut moments = super::Moments::new(times.len(), 3);
    for id in 0..paths {
        let mut p = prototype.clone();
        p.seed(seeds[id as usize]);
        p.A = 24;
        p.B = 12;
        p.C = 0;
        let mut states = Vec::new();
        for &t in times {
            p.advance_until(t);
            let state = vec![p.A as f64, p.B as f64, p.C as f64];
            if state.iter().any(|x| *x < 0.) {
                return Err("negative rebop population".into());
            }
            states.push(state);
        }
        moments.add(&states);
    }
    Ok(moments.finish())
}
