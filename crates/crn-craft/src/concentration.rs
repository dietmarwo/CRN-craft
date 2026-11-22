//! Prepared continuous rate expressions for deterministic integrations.
use crate::rules::{
    dual::{Input, Program},
    expr::{Env, Expr, Scope},
};
use std::collections::BTreeMap;

#[derive(Clone)]
struct Formula {
    expression: Expr,
    derivative: Program,
}

impl Formula {
    fn compile(source: &str, slots: &BTreeMap<String, usize>) -> Result<Self, String> {
        if source.len() > 4096 {
            return Err("concentration formula exceeds 4096 bytes".into());
        }
        let expression = Expr::parse(
            source,
            &Scope {
                params: slots,
                observables: None,
                reactants: &[],
                random: false,
            },
        )?;
        let derivative = Program::new(&expression);
        Ok(Self {
            expression,
            derivative,
        })
    }

    fn eval(&self, values: &[f64]) -> f64 {
        self.expression.eval(&Env {
            params: values,
            observables: &[],
            attributes: [&[], &[]],
        })
    }

    fn product(&self, values: &[f64], direction: &[f64]) -> f64 {
        let gradient = self.derivative.evaluate(values, &[]);
        self.derivative
            .inputs()
            .iter()
            .zip(gradient.grad)
            .map(|(slot, value)| match *slot {
                Input::Parameter(index) if direction[index] != 0.0 => value * direction[index],
                _ => 0.0,
            })
            .sum()
    }
}

/// Deterministic concentration laws in model order, with sequential named
/// intermediates. The name time is available in derivatives and assignments.
#[derive(Clone)]
pub struct Concentration {
    states: usize,
    parameters: usize,
    initial: Vec<Formula>,
    assignments: Vec<Formula>,
    derivatives: Vec<Formula>,
    observables: Vec<Formula>,
}

fn reserve(slots: &mut BTreeMap<String, usize>, name: &str) -> Result<(), String> {
    if name.is_empty()
        || name.len() > 64
        || !name
            .bytes()
            .enumerate()
            .all(|(i, c)| c.is_ascii_alphabetic() || c == b'_' || (i > 0 && c.is_ascii_digit()))
        || slots.contains_key(name)
    {
        return Err(format!("invalid or repeated concentration name {name:?}"));
    }
    slots.insert(name.to_owned(), slots.len());
    Ok(())
}

impl Concentration {
    /// Compile initial expressions, named assignments, state derivatives and
    /// readouts. Initial expressions can reference only the parameters.
    ///
    /// # Errors
    /// Invalid names, size or formulas, or inconsistent state dimensions.
    pub fn new(
        states: &[String],
        parameters: &[String],
        initial: &[String],
        assignments: &[(String, String)],
        derivatives: &[String],
        observables: &[String],
    ) -> Result<Self, String> {
        if !(1..=64).contains(&states.len())
            || initial.len() != states.len()
            || derivatives.len() != states.len()
            || parameters.len() > 128
            || assignments.len() > 128
            || observables.len() > 64
        {
            return Err("invalid concentration dimensions".into());
        }
        let mut slots = BTreeMap::new();
        for name in parameters {
            reserve(&mut slots, name)?;
        }
        let compile = |sources: &[String], slots: &_| {
            sources
                .iter()
                .map(|s| Formula::compile(s, slots))
                .collect::<Result<Vec<_>, _>>()
        };
        let initial = compile(initial, &slots)?;
        reserve(&mut slots, "time")?;
        for name in states {
            reserve(&mut slots, name)?;
        }
        let mut named = Vec::with_capacity(assignments.len());
        for (name, source) in assignments {
            named.push(Formula::compile(source, &slots)?);
            reserve(&mut slots, name)?;
        }
        Ok(Self {
            states: states.len(),
            parameters: parameters.len(),
            initial,
            assignments: named,
            derivatives: compile(derivatives, &slots)?,
            observables: compile(observables, &slots)?,
        })
    }

    /// Initial concentrations for a parameter vector.
    ///
    /// # Errors
    /// Invalid parameter or initial state values.
    pub fn initial(&self, parameters: &[f64]) -> Result<Vec<f64>, String> {
        if parameters.len() != self.parameters || parameters.iter().any(|v| !v.is_finite()) {
            return Err("invalid concentration parameters".into());
        }
        let initial: Vec<_> = self.initial.iter().map(|f| f.eval(parameters)).collect();
        if initial.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err("invalid initial concentrations".into());
        }
        Ok(initial)
    }

    fn inputs(&self, time: f64, state: &[f64], parameters: &[f64]) -> Vec<f64> {
        let mut values =
            Vec::with_capacity(self.parameters + 1 + self.states + self.assignments.len());
        values.extend_from_slice(parameters);
        values.push(time);
        values.extend(state.iter().map(|x| x.max(0.0)));
        values
    }

    /// Evaluate the ODE right-hand side. Trial states use max(state, 0).
    pub fn rhs(&self, t: f64, y: &[f64], p: &[f64], out: &mut [f64]) {
        if y.len() != self.states
            || p.len() != self.parameters
            || out.len() != self.states
            || y.iter().any(|x| !x.is_finite())
        {
            out.fill(f64::NAN);
            return;
        }
        let mut values = self.inputs(t, y, p);
        for f in &self.assignments {
            values.push(f.eval(&values));
        }
        for (v, f) in out.iter_mut().zip(&self.derivatives) {
            *v = f.eval(&values);
        }
    }

    /// State Jacobian times a direction, including the named assignments.
    pub fn jacobian_product(&self, t: f64, y: &[f64], p: &[f64], v: &[f64], out: &mut [f64]) {
        if y.len() != self.states
            || p.len() != self.parameters
            || out.len() != self.states
            || v.len() != self.states
            || y.iter().any(|x| !x.is_finite())
        {
            out.fill(f64::NAN);
            return;
        }
        let mut values = self.inputs(t, y, p);
        let mut direction = vec![0.; self.parameters + 1];
        direction.extend(
            y.iter()
                .zip(v)
                .map(|(y, v)| if *y < 0.0 { 0.0 } else { *v }),
        );
        for f in &self.assignments {
            let d = f.product(&values, &direction);
            values.push(f.eval(&values));
            direction.push(d);
        }
        for (result, f) in out.iter_mut().zip(&self.derivatives) {
            *result = f.product(&values, &direction);
        }
    }

    /// Readouts in declared order for a validated trajectory state.
    ///
    /// # Errors
    /// Invalid sizes or a non-finite readout.
    pub fn observe(&self, t: f64, y: &[f64], p: &[f64]) -> Result<Vec<f64>, String> {
        if y.len() != self.states || p.len() != self.parameters || y.iter().any(|x| !x.is_finite())
        {
            return Err("invalid concentration readout state".into());
        }
        let mut values = self.inputs(t, y, p);
        for f in &self.assignments {
            values.push(f.eval(&values));
        }
        let result: Vec<_> = self.observables.iter().map(|f| f.eval(&values)).collect();
        if result.iter().any(|x| !x.is_finite()) {
            return Err(format!("non-finite observation at {t}"));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::Concentration;

    fn strings(values: &[&str]) -> Vec<String> {
        values.iter().map(|s| (*s).to_owned()).collect()
    }

    fn transport() -> Concentration {
        Concentration::new(
            &strings(&["C1", "C2"]),
            &strings(&["k", "V1", "V2", "c0"]),
            &strings(&["c0", "c0/4"]),
            &[("J".into(), "k*(C1-C2)".into())],
            &strings(&["-J/V1", "J/V2 + time*0"]),
            &strings(&["V1*C1+V2*C2"]),
        )
        .unwrap()
    }

    #[test]
    fn right_hand_side_uses_assignments_and_parameters() {
        let model = transport();
        let p = [0.8, 1.0, 4.0, 2.0];
        assert_eq!(model.initial(&p).unwrap(), vec![2.0, 0.5]);
        let mut out = [0.0; 2];
        model.rhs(0.0, &[2.0, 0.5], &p, &mut out);
        assert!((out[0] + 1.2).abs() < 1e-15 && (out[1] - 0.3).abs() < 1e-15);
        // Amount is conserved: V1 dC1 + V2 dC2 = 0.
        assert!((out[0] + 4.0 * out[1]).abs() < 1e-15);
        assert_eq!(model.observe(1.0, &[2.0, 0.5], &p).unwrap(), vec![4.0]);
    }

    #[test]
    fn jacobian_product_matches_central_differences() {
        let model = Concentration::new(
            &strings(&["A", "B"]),
            &strings(&["k", "K"]),
            &strings(&["1", "0.5"]),
            &[
                ("v".into(), "k*A*B/(K+A)".into()),
                ("w".into(), "v^2".into()),
            ],
            &strings(&["-v+w", "v-2*w*B"]),
            &[],
        )
        .unwrap();
        let (p, y, v) = ([1.3, 0.4], [0.7, 0.9], [0.3, -1.1]);
        let mut product = [0.0; 2];
        model.jacobian_product(0.0, &y, &p, &v, &mut product);
        let h = 1e-6;
        let (mut plus, mut minus) = ([0.0; 2], [0.0; 2]);
        let shift = |s: f64| [y[0] + s * v[0], y[1] + s * v[1]];
        model.rhs(0.0, &shift(h), &p, &mut plus);
        model.rhs(0.0, &shift(-h), &p, &mut minus);
        for i in 0..2 {
            let difference = (plus[i] - minus[i]) / (2.0 * h);
            assert!((product[i] - difference).abs() < 1e-8, "{i}");
        }
    }

    #[test]
    fn negative_trial_states_are_read_as_zero() {
        let model = Concentration::new(
            &strings(&["A"]),
            &strings(&["k"]),
            &strings(&["1"]),
            &[],
            &strings(&["-k*A"]),
            &[],
        )
        .unwrap();
        let mut out = [1.0];
        model.rhs(0.0, &[-0.5], &[2.0], &mut out);
        assert_eq!(out, [0.0]);
        model.jacobian_product(0.0, &[-0.5], &[2.0], &[1.0], &mut out);
        assert_eq!(out, [0.0]);
    }

    #[test]
    fn invalid_models_and_inputs_are_refused() {
        let new = |states: &[&str], initial: &[&str], derivatives: &[&str]| {
            Concentration::new(
                &strings(states),
                &strings(&["k"]),
                &strings(initial),
                &[],
                &strings(derivatives),
                &[],
            )
        };
        assert!(new(&["A", "A"], &["1", "1"], &["0", "0"]).is_err());
        assert!(new(&["time"], &["1"], &["0"]).is_err());
        assert!(new(&["k"], &["1"], &["0"]).is_err());
        assert!(new(&["1A"], &["1"], &["0"]).is_err());
        assert!(new(&["A"], &["A"], &["0"]).is_err());
        assert!(new(&["A"], &["1"], &["-q*A"]).is_err());
        assert!(new(&["A"], &["1", "2"], &["0"]).is_err());
        assert!(new(&[], &[], &[]).is_err());
        let long = format!("{}k", "0+".repeat(2100));
        assert!(new(&["A"], &["1"], &[long.as_str()]).is_err());
        let model = new(&["A"], &["k-1"], &["-k*A"]).unwrap();
        assert!(model.initial(&[0.5]).is_err());
        assert!(model.initial(&[f64::NAN]).is_err());
        assert!(model.initial(&[1.0, 2.0]).is_err());
        assert_eq!(model.initial(&[1.5]).unwrap(), vec![0.5]);
        let mut out = [0.0];
        model.rhs(0.0, &[f64::INFINITY], &[1.0], &mut out);
        assert!(out[0].is_nan());
        assert!(model.observe(0.0, &[0.5, 0.5], &[1.0]).is_err());
        let ratio = Concentration::new(
            &strings(&["A"]),
            &strings(&["k"]),
            &strings(&["1"]),
            &[],
            &strings(&["0"]),
            &strings(&["A/k"]),
        )
        .unwrap();
        assert!(ratio.observe(0.0, &[1.0], &[0.0]).is_err());
    }
}
