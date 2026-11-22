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
