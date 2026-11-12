//! Dual-number evaluation of rule expressions, for ODE Jacobians.
//!
//! A [`Program`] evaluates a rate or propensity expression together with its
//! partial derivatives with respect to its inputs: the parameters and
//! observables it reads. The value is computed by the same operations as
//! [`Expr::eval`], so it is bit for bit the value the engine uses; the
//! derivatives follow the chain rule through every operator and function of
//! the expression language:
//!
//! - `+ - * /`, unary `-`, `^` (`∂(x^y) = y·x^(y−1)·∂x + x^y·ln x·∂y`),
//!   `exp`, `log`, `sqrt` and `abs` by their derivatives;
//! - points where an expression is not differentiable use the derivative of
//!   the branch the evaluation selects: `if` and `&&`/`||` take the branch
//!   their condition selects, `min(x, y)` takes `x` when x ≤ y and `max(x,
//!   y)` takes `x` when x ≥ y, `abs(x)` takes `x` when x ≥ 0 (as
//!   `if(x < 0, -x, x)`); comparisons, `!`, `&&`, `||`, `floor` and `ceil`
//!   have derivative 0;
//! - a partial derivative whose chain-rule factor is zero is zero, also
//!   where the other factor is infinite: a constant exponent contributes
//!   nothing at a base of 0, where `ln 0 = −∞`;
//! - a power, square root or logarithm whose derivative is infinite at 0
//!   (`x^y` with 0 < y < 1, `sqrt(x)`, `log(x)` at x = 0) gives an infinite
//!   or undefined derivative and marks the evaluation singular. At x = 0,
//!   `x^y` has derivative 0 for y > 1 and 1 for y = 1, and `∂(x^y)/∂y` is 0
//!   for y > 0 (the limits from above).
//!
//! An intermediate value can be infinite while the result is finite, as in
//! the converted Hill activation `vmax / (1 + (k / R)^n)` at R = 0, where
//! `k / R = ∞`. Such a value carries the derivatives of its reciprocal
//! (whose value is ±0), so the result's derivatives are the limits of the
//! quotient rule: for the Hill activation at R = 0 they are 0 for n > 1,
//! vmax / k with respect to R for n = 1, and infinite for n < 1, the
//! one-sided derivatives of the law.
use super::expr::{Expr, Func, Node, Op, Slot, deterministic_call, eval_binary, flag, truth};

/// An input of an expression: a parameter or an observable, by the index
/// the expression uses (sorted names).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Input {
    Parameter(usize),
    Observable(usize),
}

/// An expression prepared for dual-number evaluation with respect to the
/// parameters and observables it reads.
#[derive(Clone, Debug)]
pub(crate) struct Program {
    expr: Expr,
    inputs: Vec<Input>,
}

/// The value of a [`Program`] and its partial derivatives.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Evaluation {
    /// The value, as [`Expr::eval`] computes it.
    pub value: f64,
    /// ∂value/∂input in [`Program::inputs`] order; NaN when the value is
    /// not finite.
    pub grad: Vec<f64>,
    /// Whether a power, square root or logarithm had an infinite derivative
    /// at 0 on the way to a derivative.
    pub singular: bool,
}

/// A value and the derivatives of the value, or, when the value is ±∞, the
/// derivatives of its reciprocal 1/value = ±0.
struct Dual {
    v: f64,
    g: Vec<f64>,
}

/// `c · g`, but 0 where `g` is 0 whatever `c` is.
fn scaled(c: f64, g: f64) -> f64 {
    if g == 0.0 { 0.0 } else { c * g }
}

/// ∂(x^y)/∂x for finite x and y with a finite or zero power.
fn base_coefficient(x: f64, y: f64) -> f64 {
    if y == 0.0 {
        // x^0 = 1 for every x.
        0.0
    } else {
        y * libm::pow(x, y - 1.0)
    }
}

/// ∂(x^y)/∂y = x^y · ln x for finite x, with the limit 0 at x = 0 for y > 0.
fn exponent_coefficient(x: f64, y: f64, power: f64) -> f64 {
    if x == 0.0 && y > 0.0 {
        0.0
    } else {
        power * libm::log(x)
    }
}

struct Evaluator<'a> {
    params: &'a [f64],
    observables: &'a [f64],
    inputs: &'a [Input],
    singular: bool,
}

impl Evaluator<'_> {
    fn constant(&self, v: f64) -> Dual {
        Dual {
            v,
            g: vec![0.0; self.inputs.len()],
        }
    }

    fn undefined(&self) -> Dual {
        Dual {
            v: f64::NAN,
            g: vec![f64::NAN; self.inputs.len()],
        }
    }

    fn seed(&self, input: Input, v: f64) -> Dual {
        let mut g = vec![0.0; self.inputs.len()];
        if let Some(k) = self.inputs.iter().position(|&i| i == input) {
            g[k] = 1.0;
        }
        Dual { v, g }
    }

    /// Records an infinite derivative coefficient met by a nonzero input
    /// derivative at a zero argument.
    fn note(&mut self, coefficient: f64, g: &[f64]) {
        if !coefficient.is_finite() && g.iter().any(|&x| x != 0.0) {
            self.singular = true;
        }
    }

    fn node(&mut self, node: &Node) -> Dual {
        match node {
            Node::Num(x) => self.constant(*x),
            Node::Var(Slot::Param(i)) => self.seed(Input::Parameter(*i), self.params[*i]),
            Node::Var(Slot::Observable(i)) => {
                self.seed(Input::Observable(*i), self.observables[*i])
            }
            // Rates and propensities read no reactant attributes.
            Node::Var(Slot::Attribute { .. }) => self.undefined(),
            Node::Neg(a) => {
                let a = self.node(a);
                Dual {
                    v: -a.v,
                    g: a.g.iter().map(|x| -x).collect(),
                }
            }
            Node::Not(a) => {
                let x = self.node(a).v;
                self.constant(if x.is_nan() { x } else { flag(!truth(x)) })
            }
            Node::Bin(op, a, b) => {
                let x = self.node(a);
                // Short-circuit logic as `Expr::eval`; logic has derivative 0.
                match op {
                    Op::And | Op::Or if x.v.is_nan() => return self.constant(x.v),
                    Op::And if !truth(x.v) => return self.constant(0.0),
                    Op::Or if truth(x.v) => return self.constant(1.0),
                    _ => {}
                }
                let y = self.node(b);
                self.binary(*op, x, y)
            }
            Node::Call(Func::If, args) => {
                let c = self.node(&args[0]).v;
                if c.is_nan() {
                    return self.undefined();
                }
                if truth(c) {
                    self.node(&args[1])
                } else {
                    self.node(&args[2])
                }
            }
            Node::Call(f @ (Func::Min | Func::Max), args) => {
                let a = self.node(&args[0]);
                let b = self.node(&args[1]);
                let v = deterministic_call(*f, a.v, b.v);
                if v.is_nan() {
                    return self.undefined();
                }
                let first = if *f == Func::Min {
                    a.v <= b.v
                } else {
                    a.v >= b.v
                };
                Dual {
                    v,
                    g: if first { a.g } else { b.g },
                }
            }
            Node::Call(
                f @ (Func::Exp | Func::Log | Func::Sqrt | Func::Abs | Func::Floor | Func::Ceil),
                args,
            ) => {
                let a = self.node(&args[0]);
                self.unary(*f, a)
            }
            // Random draws are refused in rates at parse time.
            Node::Call(..) => self.undefined(),
        }
    }

    fn unary(&mut self, f: Func, a: Dual) -> Dual {
        let v = deterministic_call(f, a.v, 0.0);
        if v.is_nan() {
            return self.undefined();
        }
        let g = match f {
            Func::Floor | Func::Ceil => vec![0.0; a.g.len()],
            Func::Abs => {
                if a.v < 0.0 {
                    a.g.iter().map(|x| -x).collect()
                } else {
                    a.g
                }
            }
            Func::Exp => {
                if a.v.is_finite() && v.is_finite() {
                    a.g.iter().map(|x| v * x).collect()
                } else if a.v.is_finite() {
                    // Overflow: the reciprocal exp(-a) has derivative -exp(-a).
                    let w = libm::exp(-a.v);
                    a.g.iter().map(|x| -w * x).collect()
                } else {
                    // exp(±∞) and its reciprocal are flat in the limit.
                    a.g.iter().map(|&x| scaled(0.0, x)).collect()
                }
            }
            Func::Log => {
                if a.v.is_finite() && a.v > 0.0 {
                    a.g.iter().map(|x| x / a.v).collect()
                } else if a.v == 0.0 {
                    // log 0 = -∞; 1/log x has derivative -1/(x log² x) → -∞.
                    self.note(f64::NEG_INFINITY, &a.g);
                    a.g.iter().map(|&x| scaled(f64::NEG_INFINITY, x)).collect()
                } else {
                    // log(+∞) = +∞; its reciprocal has an infinite derivative.
                    a.g.iter().map(|&x| scaled(f64::INFINITY, x)).collect()
                }
            }
            Func::Sqrt => {
                if a.v.is_finite() && a.v > 0.0 {
                    a.g.iter().map(|x| x / (2.0 * v)).collect()
                } else {
                    // At 0, and for the reciprocal sqrt(1/a) at a = +∞.
                    if a.v == 0.0 {
                        self.note(f64::INFINITY, &a.g);
                    }
                    a.g.iter().map(|&x| scaled(f64::INFINITY, x)).collect()
                }
            }
            _ => unreachable!("not a unary function"),
        };
        Dual { v, g }
    }

    #[allow(clippy::too_many_lines)]
    fn binary(&mut self, op: Op, a: Dual, b: Dual) -> Dual {
        let v = eval_binary(op, a.v, b.v);
        let n = a.g.len();
        if matches!(
            op,
            Op::Eq | Op::Ne | Op::Lt | Op::Le | Op::Gt | Op::Ge | Op::And | Op::Or
        ) {
            return self.constant(v);
        }
        if v.is_nan() {
            return self.undefined();
        }
        let (fa, fb) = (a.v.is_finite(), b.v.is_finite());
        // The reciprocal derivatives of a result that overflowed from finite
        // operands: -g / v² with v = ±∞, so 0 for every finite g.
        let overflow = |g: Vec<f64>| -> Vec<f64> {
            g.into_iter()
                .map(|x| if x.is_finite() { 0.0 } else { f64::NAN })
                .collect()
        };
        // Both operands infinite and the result infinite: the reciprocal
        // derivative is defined only where neither operand moves.
        let both = |a: &[f64], b: &[f64]| -> Vec<f64> {
            a.iter()
                .zip(b)
                .map(|(&x, &y)| if x == 0.0 && y == 0.0 { 0.0 } else { f64::NAN })
                .collect()
        };
        let g = match op {
            Op::Add | Op::Sub => {
                let sign = if op == Op::Sub { -1.0 } else { 1.0 };
                match (fa, fb) {
                    (true, true) => {
                        let g = a.g.iter().zip(&b.g).map(|(x, y)| x + sign * y).collect();
                        if v.is_finite() { g } else { overflow(g) }
                    }
                    // ∞ + c: 1/(a + c) = r/(1 + c r) has derivative ∂r at r = 0.
                    (false, true) => a.g,
                    (true, false) => b.g.iter().map(|y| sign * y).collect(),
                    (false, false) => both(&a.g, &b.g),
                }
            }
            Op::Mul => match (fa, fb) {
                (true, true) => {
                    let g =
                        a.g.iter()
                            .zip(&b.g)
                            .map(|(x, y)| b.v * x + a.v * y)
                            .collect();
                    if v.is_finite() { g } else { overflow(g) }
                }
                // 1/(a c) = r / c.
                (false, true) => a.g.iter().map(|x| x / b.v).collect(),
                (true, false) => b.g.iter().map(|y| y / a.v).collect(),
                // 1/(a b) = r s has derivative r ∂s + s ∂r = 0 at r = s = 0.
                (false, false) => {
                    a.g.iter()
                        .zip(&b.g)
                        .map(|(&x, &y)| scaled(0.0, x) + scaled(0.0, y))
                        .collect()
                }
            },
            Op::Div => match (fa, fb) {
                (true, true) if v.is_finite() => {
                    a.g.iter()
                        .zip(&b.g)
                        .map(|(x, y)| (x - v * y) / b.v)
                        .collect()
                }
                // a / b = ±∞ (b = 0 or overflow): the reciprocal b / a has
                // derivative (∂b - (b / a) ∂a) / a.
                (true, true) => {
                    let r = b.v / a.v;
                    a.g.iter()
                        .zip(&b.g)
                        .map(|(x, y)| (y - r * x) / a.v)
                        .collect()
                }
                // a / ±∞ = a · s with s = ±0: derivative a ∂s + s ∂a.
                (true, false) => {
                    let s = 1.0 / b.v;
                    a.g.iter()
                        .zip(&b.g)
                        .map(|(&x, &y)| a.v * y + scaled(s, x))
                        .collect()
                }
                // ±∞ / c: the reciprocal c · r has derivative c ∂r + r ∂c.
                (false, true) => {
                    let r = 1.0 / a.v;
                    a.g.iter()
                        .zip(&b.g)
                        .map(|(&x, &y)| b.v * x + scaled(r, y))
                        .collect()
                }
                (false, false) => unreachable!("∞/∞ is NaN"),
            },
            Op::Pow => self.power(a, b, v),
            _ => unreachable!("comparisons and logic handled above"),
        };
        debug_assert_eq!(g.len(), n);
        Dual { v, g }
    }

    fn power(&mut self, a: Dual, b: Dual, v: f64) -> Vec<f64> {
        let (x, y) = (a.v, b.v);
        let combine = |this: &mut Self, ca: f64, cb: f64, ga: &[f64], gb: &[f64]| -> Vec<f64> {
            if x == 0.0 {
                this.note(ca, ga);
                this.note(cb, gb);
            }
            ga.iter()
                .zip(gb)
                .map(|(&p, &q)| scaled(ca, p) + scaled(cb, q))
                .collect()
        };
        match (x.is_finite(), y.is_finite()) {
            (true, true) if v.is_finite() => {
                let ca = base_coefficient(x, y);
                let cb = exponent_coefficient(x, y, v);
                combine(self, ca, cb, &a.g, &b.g)
            }
            (true, true) => {
                // x^y = ±∞ (x = 0 with y < 0, or overflow): the reciprocal
                // x^(-y) has derivative -y x^(-y-1) ∂x - x^(-y) ln x ∂y.
                let w = libm::pow(x, -y);
                let ca = base_coefficient(x, -y);
                let cb = -exponent_coefficient(x, -y, w);
                combine(self, ca, cb, &a.g, &b.g)
            }
            (false, true) if x > 0.0 => {
                // a = +∞ with reciprocal r = +0 and derivative ∂r: a^y is
                // (1/r)^y = r^(-y). For y > 0 the result is +∞ and its
                // reciprocal r^y; for y < 0 the result is r^(-y). The
                // derivative with respect to y vanishes in the limit,
                // except at y = 0 where it is ln ∞.
                let ca = if y > 0.0 {
                    base_coefficient(0.0, y)
                } else {
                    base_coefficient(0.0, -y)
                };
                let cb = if y == 0.0 { f64::INFINITY } else { 0.0 };
                if !ca.is_finite() && a.g.iter().any(|&p| p != 0.0) {
                    self.singular = true;
                }
                a.g.iter()
                    .zip(&b.g)
                    .map(|(&p, &q)| scaled(ca, p) + scaled(cb, q))
                    .collect()
            }
            // A base of -∞ or an infinite exponent: defined only where
            // neither operand moves.
            _ => {
                a.g.iter()
                    .zip(&b.g)
                    .map(|(&p, &q)| if p == 0.0 && q == 0.0 { 0.0 } else { f64::NAN })
                    .collect()
            }
        }
    }
}

impl Program {
    /// Prepares `expr`; its inputs are the parameters and observables it
    /// reads, parameters first, each in increasing index order.
    pub(crate) fn new(expr: &Expr) -> Self {
        fn visit(node: &Node, out: &mut Vec<Input>) {
            match node {
                Node::Var(Slot::Param(i)) => out.push(Input::Parameter(*i)),
                Node::Var(Slot::Observable(i)) => out.push(Input::Observable(*i)),
                Node::Num(_) | Node::Var(Slot::Attribute { .. }) => {}
                Node::Neg(a) | Node::Not(a) => visit(a, out),
                Node::Bin(_, a, b) => {
                    visit(a, out);
                    visit(b, out);
                }
                Node::Call(_, args) => args.iter().for_each(|a| visit(a, out)),
            }
        }
        let mut inputs = Vec::new();
        visit(&expr.root, &mut inputs);
        inputs.sort_by_key(|i| match *i {
            Input::Parameter(k) => (0, k),
            Input::Observable(k) => (1, k),
        });
        inputs.dedup();
        Self {
            expr: expr.clone(),
            inputs,
        }
    }

    /// The parameters and observables the derivatives are taken with
    /// respect to.
    pub(crate) fn inputs(&self) -> &[Input] {
        &self.inputs
    }

    /// The value and its derivatives at the given parameter and observable
    /// values (indexed as the expression's slots).
    pub(crate) fn evaluate(&self, params: &[f64], observables: &[f64]) -> Evaluation {
        let mut evaluator = Evaluator {
            params,
            observables,
            inputs: &self.inputs,
            singular: false,
        };
        let dual = evaluator.node(&self.expr.root);
        let grad = if dual.v.is_finite() {
            dual.g
        } else {
            vec![f64::NAN; self.inputs.len()]
        };
        Evaluation {
            value: dual.v,
            grad,
            singular: evaluator.singular,
        }
    }

    /// Whether the expression is 0 whenever every observable for which
    /// `zero` holds is 0, judged from its structure alone: a flagged
    /// observable or the literal 0; a product with such a factor; a
    /// quotient, negation, `sqrt`, `abs`, `floor`, `ceil` or a power with a
    /// positive literal exponent of such an expression; a sum, difference,
    /// `min`, `max` or both branches of `if` of such expressions. A value
    /// that becomes undefined (`0 · ∞`, `0 / 0`) is not finite and is caught
    /// by the evaluation.
    pub(crate) fn vanishes(&self, zero: &dyn Fn(usize) -> bool) -> bool {
        fn check(node: &Node, zero: &dyn Fn(usize) -> bool) -> bool {
            match node {
                Node::Num(x) => *x == 0.0,
                Node::Var(Slot::Observable(o)) => zero(*o),
                Node::Var(_) | Node::Not(_) => false,
                Node::Neg(a) => check(a, zero),
                Node::Bin(Op::Mul, a, b) => check(a, zero) || check(b, zero),
                Node::Bin(Op::Div, a, _) => check(a, zero),
                Node::Bin(Op::Add | Op::Sub, a, b) => check(a, zero) && check(b, zero),
                Node::Bin(Op::Pow, a, b) => {
                    matches!(b.as_ref(), Node::Num(y) if *y > 0.0) && check(a, zero)
                }
                Node::Bin(..) => false,
                Node::Call(Func::Sqrt | Func::Abs | Func::Floor | Func::Ceil, args) => {
                    check(&args[0], zero)
                }
                Node::Call(Func::Min | Func::Max, args) => {
                    check(&args[0], zero) && check(&args[1], zero)
                }
                Node::Call(Func::If, args) => check(&args[1], zero) && check(&args[2], zero),
                Node::Call(..) => false,
            }
        }
        check(&self.expr.root, zero)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::expr::{Env, Scope};
    use rand_core::{RngCore, SeedableRng};
    use std::collections::BTreeMap;

    /// Parses `source` over the parameters `k`, `n` and the observables
    /// `x`, `y`.
    fn program(source: &str) -> Program {
        let params = BTreeMap::from([("k".to_string(), 0), ("n".to_string(), 1)]);
        let names = BTreeMap::from([("x".to_string(), 0), ("y".to_string(), 1)]);
        let scope = Scope {
            params: &params,
            observables: Some(&names),
            reactants: &[],
            random: false,
        };
        Program::new(&Expr::parse(source, &scope).unwrap())
    }

    /// The derivatives in the fixed order k, n, x, y (0 for an input the
    /// expression does not read).
    fn full(p: &Program, k: f64, n: f64, x: f64, y: f64) -> (f64, [f64; 4], bool) {
        let e = p.evaluate(&[k, n], &[x, y]);
        let mut out = [0.0; 4];
        for (input, g) in p.inputs().iter().zip(&e.grad) {
            let at = match input {
                Input::Parameter(i) => *i,
                Input::Observable(i) => 2 + i,
            };
            out[at] = *g;
        }
        (e.value, out, e.singular)
    }

    fn eval(p: &Program, k: f64, n: f64, x: f64, y: f64) -> f64 {
        p.expr.eval(&Env {
            params: &[k, n],
            observables: &[x, y],
            attributes: [&[], &[]],
        })
    }

    #[test]
    fn values_equal_the_tree_evaluation_bit_for_bit() {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0x5e);
        fn expression(rng: &mut impl RngCore, depth: u8) -> String {
            let draw = rng.next_u32();
            if depth == 0 {
                return ["x", "y", "k", "n", "0", "-0", "0/0", "1/0", "2", "0.5"]
                    [(draw % 10) as usize]
                    .into();
            }
            let a = expression(rng, depth - 1);
            let b = expression(rng, depth - 1);
            match draw % 23 {
                0 => format!("max(({a}),({b}))"),
                1 => format!("min(({a}),({b}))"),
                2 => format!("if(({a}),({b}),-({b}))"),
                3 => format!("!({a})"),
                4 => format!("exp({a})"),
                5 => format!("log({a})"),
                6 => format!("sqrt({a})"),
                7 => format!("abs({a})"),
                8 => format!("floor({a})"),
                9 => format!("ceil({a})"),
                k => {
                    let op = [
                        "+", "-", "*", "/", "^", "<", "!=", "&&", "||", ">=", "==", "<=", ">",
                    ][(k - 10) as usize];
                    format!("({a}){op}({b})")
                }
            }
        }
        for _ in 0..2048 {
            let source = expression(&mut rng, 4);
            let p = program(&source);
            for [k, n, x, y] in [
                [0.0, -0.0, 0.0, -0.0],
                [2.0, 0.5, -2.0, 2.0],
                [3.0, 2.5, f64::INFINITY, 1.0],
                [1e300, 1e-300, 1.0, 0.4],
            ] {
                let (a, b) = (eval(&p, k, n, x, y), p.evaluate(&[k, n], &[x, y]).value);
                assert!(
                    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
                    "{source}: {a:?} != {b:?}"
                );
            }
        }
    }

    /// Central differences of smooth expressions at points away from their
    /// kinks agree with the dual derivatives.
    #[test]
    fn derivatives_match_central_differences() {
        let sources = [
            "k * x * y",
            "k * x^2 / (n + y)",
            "x^y",
            "k^n",
            "exp(-k * x) * y",
            "log(x + y) * sqrt(k * y)",
            "max(k * x, y) - min(x, n * y)",
            "abs(x - 3 * y) + floor(x) * y + ceil(y) * k",
            "if(x > y, k * x, n / y)",
            "(x < 5) * k + (x >= 5 && y > 0) * n * y",
            "k + n / (1 + (k / x)^n)",
            "k + n / (1 + (x / k)^n)",
            "-x^2 / -(y - 9)",
        ];
        let points = [[1.7, 2.3, 3.1, 0.7], [0.4, 3.9, 7.3, 4.6]];
        for source in sources {
            let p = program(source);
            for point in points {
                let (_, exact, singular) = full(&p, point[0], point[1], point[2], point[3]);
                assert!(!singular, "{source}");
                for j in 0..4 {
                    let h = 1e-6 * point[j].abs().max(1.0);
                    let (mut plus, mut minus) = (point, point);
                    plus[j] += h;
                    minus[j] -= h;
                    let f = |q: [f64; 4]| eval(&p, q[0], q[1], q[2], q[3]);
                    let numeric = (f(plus) - f(minus)) / (2.0 * h);
                    assert!(
                        (numeric - exact[j]).abs() <= 1e-6 * numeric.abs().max(1.0),
                        "{source} at {point:?}, input {j}: {} vs {numeric}",
                        exact[j]
                    );
                }
            }
        }
    }

    #[test]
    fn branch_rules_at_points_without_a_derivative() {
        // min and max take the first argument on ties; abs takes x at 0.
        assert_eq!(
            full(&program("min(x, y)"), 0., 0., 2., 2.).1,
            [0., 0., 1., 0.]
        );
        assert_eq!(
            full(&program("max(y, x)"), 0., 0., 2., 2.).1,
            [0., 0., 0., 1.]
        );
        assert_eq!(
            full(&program("min(x, y)"), 0., 0., 3., 2.).1,
            [0., 0., 0., 1.]
        );
        assert_eq!(full(&program("abs(x)"), 0., 0., 0., 0.).1, [0., 0., 1., 0.]);
        assert_eq!(
            full(&program("abs(x)"), 0., 0., -1., 0.).1,
            [0., 0., -1., 0.]
        );
        // floor, ceil, comparisons and logic are flat; if takes its branch.
        assert_eq!(full(&program("floor(k * x)"), 2., 0., 1., 0.).1, [0.; 4]);
        assert_eq!(
            full(&program("ceil(x) + (x > y)"), 0., 0., 1., 1.).1,
            [0.; 4]
        );
        assert_eq!(
            full(&program("x && k * y || n"), 1., 1., 1., 1.),
            (1.0, [0.0; 4], false)
        );
        let branch = program("if(x >= 1, k * x, n * y)");
        assert_eq!(full(&branch, 2., 3., 1., 5.).1, [1., 0., 2., 0.]);
        assert_eq!(full(&branch, 2., 3., 0.5, 5.).1, [0., 5., 0., 3.]);
        // An undefined condition gives an undefined value.
        let e = program("if(0/0, x, y)").evaluate(&[1., 1.], &[1., 1.]);
        assert!(e.value.is_nan() && e.grad.iter().all(|g| g.is_nan()));
    }

    #[test]
    fn powers_at_a_zero_base_take_the_limits_from_above() {
        // x^n at x = 0: derivative 0 for n > 1, 1 for n = 1, singular below.
        let p = program("x^n");
        assert_eq!(full(&p, 0., 2.5, 0., 0.), (0.0, [0.; 4], false));
        assert_eq!(full(&p, 0., 1., 0., 0.), (0.0, [0., 0., 1., 0.], false));
        let (v, g, singular) = full(&p, 0., 0.5, 0., 0.);
        assert!(v == 0.0 && g[2].is_infinite() && singular);
        // A literal exponent contributes nothing although ln 0 = -∞.
        assert_eq!(full(&program("x^2"), 0., 0., 0., 0.).1, [0.; 4]);
        // x^0 = 1 is flat in x; its exponent derivative at 0 is ln 0.
        let (v, g, singular) = full(&p, 0., 0., 0., 0.);
        assert!(v == 1.0 && g[2] == 0.0 && g[1] == f64::NEG_INFINITY && singular);
        // sqrt and log at 0 are singular; a branch that discards them is not.
        assert!(full(&program("sqrt(x)"), 0., 0., 0., 0.).2);
        assert!(full(&program("log(x) * 0 + y"), 0., 0., 0., 1.).1[2].is_nan());
        let e = program("max(sqrt(x), 1)").evaluate(&[0., 0.], &[0., 0.]);
        assert_eq!((e.value, e.grad.clone()), (1.0, vec![0.0]));
    }

    #[test]
    fn infinite_intermediates_give_the_one_sided_hill_derivatives() {
        // The converted Hill activation at R = 0: k / R = ∞.
        // vmax = 30, K = 8, basal = 1: the derivative is vmax / K = 3.75 for
        // n = 1 and 0 for n > 1.
        let activation = program("y + k / (1 + (n / x)^2)");
        assert_eq!(
            full(&activation, 30., 8., 0., 1.),
            (1.0, [0., 0., 0., 1.], false)
        );
        let linear = program("y + k / (1 + (n / x)^1)");
        let (v, g, singular) = full(&linear, 30., 8., 0., 1.);
        assert_eq!((v, g, singular), (1.0, [0., 0., 3.75, 1.], false));
        // Exponent as a parameter: d/dn vanishes in the limit for n > 0.
        let general = program("k / (1 + (4 / x)^n)");
        assert_eq!(full(&general, 30., 2.5, 0., 0.), (0.0, [0.; 4], false));
        assert_eq!(full(&general, 30., 1., 0., 0.).1, [0., 0., 7.5, 0.]);
        let (v, g, singular) = full(&general, 30., 0.5, 0., 0.);
        assert!(v == 0.0 && g[2].is_infinite() && singular);
        // The repression form at R = 0 needs no infinity.
        let repression = program("k / (1 + (x / 4)^n)");
        assert_eq!(
            full(&repression, 30., 2.5, 0., 0.),
            (30.0, [1., 0., 0., 0.], false)
        );
        assert_eq!(full(&repression, 30., 1., 0., 0.).1, [1., 0., -7.5, 0.]);
        // Overflow of an intermediate: (x / k)^n = ∞ for a huge ratio.
        let (v, g, singular) = full(&repression, 1e-300, 400., 1e300, 0.);
        assert!(
            v == 0.0 && g.iter().all(|x| x.is_finite()) && !singular,
            "{g:?}"
        );
        // a / 0 with both operands moving, then back to a finite value.
        let p = program("1 / (k / x)");
        assert_eq!(full(&p, 2., 0., 0., 0.), (0.0, [0., 0., 0.5, 0.], false));
        // exp(-∞) is flat; exp(log x) at x = 0 is undefined, not 0.
        assert_eq!(
            full(&program("exp(-k / x)"), 1., 0., 0., 0.),
            (0.0, [0.; 4], false)
        );
        assert!(full(&program("exp(log(x))"), 0., 0., 0., 0.).1[2].is_nan());
        // An infinite value has no derivatives.
        let e = program("k / x").evaluate(&[1., 0.], &[0., 0.]);
        assert!(e.value.is_infinite() && e.grad.iter().all(|g| g.is_nan()));
    }

    #[test]
    fn vanishing_is_judged_from_the_structure() {
        let zero_x = |o: usize| o == 0;
        for (source, expected) in [
            ("x", true),
            ("k * (x / (n + x))", true),
            ("k * x * y", true),
            ("x^2 * k", true),
            ("x^n", false),
            ("x + y", false),
            ("x - 2 * x", true),
            ("k + n / (1 + (k / x)^n)", false),
            ("if(y > 1, x, 0)", true),
            ("min(x, y)", false),
            ("max(x, -x) / k", true),
            ("sqrt(abs(floor(ceil(-x))))", true),
            ("(x > 0) * k", false),
            ("exp(x) - 1", false),
            ("0 * y", true),
        ] {
            assert_eq!(program(source).vanishes(&zero_x), expected, "{source}");
        }
    }
}
