use rebop::define_system;
define_system! {
birth death;
Model { X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27, X28, X29, X30, X31, X32 }
birth_X0:  => X0 @ birth
death_X0: X0 =>  @ death
birth_X1:  => X1 @ birth
death_X1: X1 =>  @ death
birth_X2:  => X2 @ birth
death_X2: X2 =>  @ death
birth_X3:  => X3 @ birth
death_X3: X3 =>  @ death
birth_X4:  => X4 @ birth
death_X4: X4 =>  @ death
birth_X5:  => X5 @ birth
death_X5: X5 =>  @ death
birth_X6:  => X6 @ birth
death_X6: X6 =>  @ death
birth_X7:  => X7 @ birth
death_X7: X7 =>  @ death
birth_X8:  => X8 @ birth
death_X8: X8 =>  @ death
birth_X9:  => X9 @ birth
death_X9: X9 =>  @ death
birth_X10:  => X10 @ birth
death_X10: X10 =>  @ death
birth_X11:  => X11 @ birth
death_X11: X11 =>  @ death
birth_X12:  => X12 @ birth
death_X12: X12 =>  @ death
birth_X13:  => X13 @ birth
death_X13: X13 =>  @ death
birth_X14:  => X14 @ birth
death_X14: X14 =>  @ death
birth_X15:  => X15 @ birth
death_X15: X15 =>  @ death
birth_X16:  => X16 @ birth
death_X16: X16 =>  @ death
birth_X17:  => X17 @ birth
death_X17: X17 =>  @ death
birth_X18:  => X18 @ birth
death_X18: X18 =>  @ death
birth_X19:  => X19 @ birth
death_X19: X19 =>  @ death
birth_X20:  => X20 @ birth
death_X20: X20 =>  @ death
birth_X21:  => X21 @ birth
death_X21: X21 =>  @ death
birth_X22:  => X22 @ birth
death_X22: X22 =>  @ death
birth_X23:  => X23 @ birth
death_X23: X23 =>  @ death
birth_X24:  => X24 @ birth
death_X24: X24 =>  @ death
birth_X25:  => X25 @ birth
death_X25: X25 =>  @ death
birth_X26:  => X26 @ birth
death_X26: X26 =>  @ death
birth_X27:  => X27 @ birth
death_X27: X27 =>  @ death
birth_X28:  => X28 @ birth
death_X28: X28 =>  @ death
birth_X29:  => X29 @ birth
death_X29: X29 =>  @ death
birth_X30:  => X30 @ birth
death_X30: X30 =>  @ death
birth_X31:  => X31 @ birth
death_X31: X31 =>  @ death
birth_X32:  => X32 @ birth
death_X32: X32 =>  @ death
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
    let mut moments = super::Moments::new(times.len(), 33);
    for id in 0..paths {
        let mut p = prototype.clone();
        p.seed(seeds[id as usize]);
        p.X0 = 0;
        p.X1 = 0;
        p.X2 = 0;
        p.X3 = 0;
        p.X4 = 0;
        p.X5 = 0;
        p.X6 = 0;
        p.X7 = 0;
        p.X8 = 0;
        p.X9 = 0;
        p.X10 = 0;
        p.X11 = 0;
        p.X12 = 0;
        p.X13 = 0;
        p.X14 = 0;
        p.X15 = 0;
        p.X16 = 0;
        p.X17 = 0;
        p.X18 = 0;
        p.X19 = 0;
        p.X20 = 0;
        p.X21 = 0;
        p.X22 = 0;
        p.X23 = 0;
        p.X24 = 0;
        p.X25 = 0;
        p.X26 = 0;
        p.X27 = 0;
        p.X28 = 0;
        p.X29 = 0;
        p.X30 = 0;
        p.X31 = 0;
        p.X32 = 0;
        let mut states = Vec::new();
        for &t in times {
            p.advance_until(t);
            let state = vec![
                p.X0 as f64,
                p.X1 as f64,
                p.X2 as f64,
                p.X3 as f64,
                p.X4 as f64,
                p.X5 as f64,
                p.X6 as f64,
                p.X7 as f64,
                p.X8 as f64,
                p.X9 as f64,
                p.X10 as f64,
                p.X11 as f64,
                p.X12 as f64,
                p.X13 as f64,
                p.X14 as f64,
                p.X15 as f64,
                p.X16 as f64,
                p.X17 as f64,
                p.X18 as f64,
                p.X19 as f64,
                p.X20 as f64,
                p.X21 as f64,
                p.X22 as f64,
                p.X23 as f64,
                p.X24 as f64,
                p.X25 as f64,
                p.X26 as f64,
                p.X27 as f64,
                p.X28 as f64,
                p.X29 as f64,
                p.X30 as f64,
                p.X31 as f64,
                p.X32 as f64,
            ];
            if state.iter().any(|x| *x < 0.) {
                return Err("negative rebop population".into());
            }
            states.push(state);
        }
        moments.add(&states);
    }
    Ok(moments.finish())
}
