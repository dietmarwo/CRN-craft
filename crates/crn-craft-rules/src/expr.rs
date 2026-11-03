//! Arithmetic expressions for rates, conditions and product attributes.
//!
//! Grammar (lowest precedence first): `||`, `&&`, `== !=`, `< <= > >=`,
//! `+ -`, `* /`, unary `- !`, `^` (right associative). Comparisons and logic
//! return 1 or 0; any nonzero value is true. Names are parameters,
//! observables or `label.attribute` of a reactant. Functions: `min max exp
//! log sqrt abs floor ceil if`, and in product attributes only the random
//! draws `bernoulli(p) choice(p, a, b) uniform(a, b) normal(m, s) new_id()`.

use rand_core::RngCore;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Op {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    And,
    Or,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Func {
    Min,
    Max,
    Exp,
    Log,
    Sqrt,
    Abs,
    Floor,
    Ceil,
    If,
    Bernoulli,
    Choice,
    Uniform,
    Normal,
    NewId,
}

impl Func {
    fn parse(name: &str) -> Option<(Self, usize, bool)> {
        // (function, arity, random)
        Some(match name {
            "min" => (Self::Min, 2, false),
            "max" => (Self::Max, 2, false),
            "exp" => (Self::Exp, 1, false),
            "log" => (Self::Log, 1, false),
            "sqrt" => (Self::Sqrt, 1, false),
            "abs" => (Self::Abs, 1, false),
            "floor" => (Self::Floor, 1, false),
            "ceil" => (Self::Ceil, 1, false),
            "if" => (Self::If, 3, false),
            "bernoulli" => (Self::Bernoulli, 1, true),
            "choice" => (Self::Choice, 3, true),
            "uniform" => (Self::Uniform, 2, true),
            "normal" => (Self::Normal, 2, true),
            "new_id" => (Self::NewId, 0, true),
            _ => return None,
        })
    }
}

/// A resolved variable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Slot {
    Param(usize),
    Observable(usize),
    Attribute { reactant: usize, index: usize },
}

#[derive(Debug, Clone, PartialEq)]
enum Node {
    Num(f64),
    Var(Slot),
    Neg(Box<Node>),
    Not(Box<Node>),
    Bin(Op, Box<Node>, Box<Node>),
    Call(Func, Vec<Node>),
}

/// Names an expression may use, decided by where it appears.
pub(crate) struct Scope<'a> {
    pub params: &'a BTreeMap<String, usize>,
    pub observables: Option<&'a BTreeMap<String, usize>>,
    /// (label, attribute names) per reactant.
    pub reactants: &'a [(&'a str, &'a [String])],
    pub random: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct Expr {
    root: Node,
}

/// Values an expression is evaluated against.
pub(crate) struct Env<'a> {
    pub params: &'a [f64],
    pub observables: &'a [f64],
    pub attributes: [&'a [f64]; 2],
}

/// Random draws and the replicate's identity counter, for product attributes.
pub(crate) struct Draw<'a, R: RngCore> {
    pub rng: &'a mut R,
    pub next_id: &'a mut f64,
    pub words: &'a mut u64,
    /// Set when `new_id()` would return 2^53 or more.
    pub id_overflow: &'a mut bool,
}

/// First identity `new_id()` cannot return exactly.
pub(crate) const ID_LIMIT: f64 = 9_007_199_254_740_992.0;

fn uniform(word: u64) -> f64 {
    (word >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
}

impl Expr {
    /// Recognize the guarded fixed-population takeover predicate.
    pub(crate) fn fraction_stop(&self, total: usize, mutants: usize) -> Option<f64> {
        let Node::Bin(Op::And, guard, predicate) = &self.root else {
            return None;
        };
        let total = Node::Var(Slot::Observable(total));
        let mutants = Node::Var(Slot::Observable(mutants));
        let Node::Bin(Op::Gt, left, right) = guard.as_ref() else {
            return None;
        };
        if **left != total || **right != Node::Num(0.0) {
            return None;
        }
        let Node::Bin(Op::Ge, left, right) = predicate.as_ref() else {
            return None;
        };
        let Node::Bin(Op::Mul, fraction, population) = right.as_ref() else {
            return None;
        };
        match fraction.as_ref() {
            Node::Num(f)
                if **left == mutants && **population == total && (0.0..=1.0).contains(f) =>
            {
                Some(*f)
            }
            _ => None,
        }
    }

    pub(crate) fn parse(source: &str, scope: &Scope) -> Result<Self, String> {
        let tokens = lex(source)?;
        let mut parser = Parser {
            tokens,
            pos: 0,
            scope,
        };
        let root = parser.expression(0)?;
        if parser.pos != parser.tokens.len() {
            return Err(format!("unexpected `{}`", parser.tokens[parser.pos].text()));
        }
        Ok(Self { root })
    }

    pub(crate) fn constant(value: f64) -> Self {
        Self {
            root: Node::Num(value),
        }
    }

    /// Structural equality of optional expressions (labels already resolved
    /// to reactant positions).
    pub(crate) fn same_option(a: Option<&Self>, b: Option<&Self>) -> bool {
        match (a, b) {
            (None, None) => true,
            (Some(a), Some(b)) => a.root == b.root,
            _ => false,
        }
    }

    /// Whether the expression calls `new_id()`.
    pub(crate) fn uses_new_id(&self) -> bool {
        fn walk(node: &Node) -> bool {
            match node {
                Node::Num(_) | Node::Var(_) => false,
                Node::Neg(a) | Node::Not(a) => walk(a),
                Node::Bin(_, a, b) => walk(a) || walk(b),
                Node::Call(f, args) => *f == Func::NewId || args.iter().any(walk),
            }
        }
        walk(&self.root)
    }

    /// Deterministic evaluation; random functions are rejected at parse time
    /// wherever this is used.
    pub(crate) fn eval(&self, env: &Env) -> f64 {
        eval_node::<NoRng>(&self.root, env, &mut None)
    }

    pub(crate) fn eval_random<R: RngCore>(&self, env: &Env, draw: &mut Draw<R>) -> f64 {
        let mut d = Some(Draw {
            rng: &mut *draw.rng,
            next_id: &mut *draw.next_id,
            words: &mut *draw.words,
            id_overflow: &mut *draw.id_overflow,
        });
        eval_node(&self.root, env, &mut d)
    }
}

struct NoRng;
impl RngCore for NoRng {
    fn next_u32(&mut self) -> u32 {
        unreachable!("deterministic expression drew a random number")
    }
    fn next_u64(&mut self) -> u64 {
        unreachable!("deterministic expression drew a random number")
    }
    fn fill_bytes(&mut self, _: &mut [u8]) {
        unreachable!()
    }
    fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), rand_core::Error> {
        unreachable!()
    }
}

fn truth(x: f64) -> bool {
    x != 0.0
}

fn flag(x: bool) -> f64 {
    f64::from(u8::from(x))
}

fn eval_node<R: RngCore>(node: &Node, env: &Env, draw: &mut Option<Draw<R>>) -> f64 {
    match node {
        Node::Num(x) => *x,
        Node::Var(Slot::Param(i)) => env.params[*i],
        Node::Var(Slot::Observable(i)) => env.observables[*i],
        Node::Var(Slot::Attribute { reactant, index }) => env.attributes[*reactant][*index],
        Node::Neg(a) => -eval_node(a, env, draw),
        // An undefined (NaN) condition stays NaN through logic, so `where`
        // and `stop` can reject it instead of treating it as true.
        Node::Not(a) => {
            let x = eval_node(a, env, draw);
            if x.is_nan() { x } else { flag(!truth(x)) }
        }
        Node::Bin(op, a, b) => {
            let x = eval_node(a, env, draw);
            // Short-circuit logic, so a guarded division is not evaluated.
            match op {
                Op::And | Op::Or if x.is_nan() => return x,
                Op::And if !truth(x) => return 0.0,
                Op::Or if truth(x) => return 1.0,
                _ => {}
            }
            let y = eval_node(b, env, draw);
            match op {
                Op::Add => x + y,
                Op::Sub => x - y,
                Op::Mul => x * y,
                Op::Div => x / y,
                Op::Pow => libm::pow(x, y),
                Op::Eq => f64::from(u8::from(x == y)),
                Op::Ne => f64::from(u8::from(x != y)),
                Op::Lt => f64::from(u8::from(x < y)),
                Op::Le => f64::from(u8::from(x <= y)),
                Op::Gt => f64::from(u8::from(x > y)),
                Op::Ge => f64::from(u8::from(x >= y)),
                Op::And | Op::Or if y.is_nan() => y,
                Op::And | Op::Or => flag(truth(y)),
            }
        }
        Node::Call(f, args) => {
            if *f == Func::If {
                let c = eval_node(&args[0], env, draw);
                if c.is_nan() {
                    return c;
                }
                return if truth(c) {
                    eval_node(&args[1], env, draw)
                } else {
                    eval_node(&args[2], env, draw)
                };
            }
            if *f == Func::NewId {
                let d = draw.as_mut().expect("random context");
                let id = *d.next_id;
                if id >= ID_LIMIT {
                    *d.id_overflow = true;
                    return f64::NAN;
                }
                *d.next_id += 1.0;
                return id;
            }
            if *f == Func::Choice {
                // Only the chosen branch is evaluated (and draws).
                let p = eval_node(&args[0], env, draw);
                let d = draw.as_mut().expect("random context");
                *d.words += 1;
                let u = uniform(d.rng.next_u64());
                return if u < p {
                    eval_node(&args[1], env, draw)
                } else {
                    eval_node(&args[2], env, draw)
                };
            }
            let v: Vec<f64> = args.iter().map(|a| eval_node(a, env, draw)).collect();
            let mut word = || {
                let d = draw.as_mut().expect("random context");
                *d.words += 1;
                d.rng.next_u64()
            };
            match f {
                Func::Min | Func::Max if v[0].is_nan() || v[1].is_nan() => f64::NAN,
                Func::Min => v[0].min(v[1]),
                Func::Max => v[0].max(v[1]),
                Func::Exp => libm::exp(v[0]),
                Func::Log => libm::log(v[0]),
                Func::Sqrt => libm::sqrt(v[0]),
                Func::Abs => v[0].abs(),
                Func::Floor => libm::floor(v[0]),
                Func::Ceil => libm::ceil(v[0]),
                Func::Bernoulli => f64::from(u8::from(uniform(word()) < v[0])),
                Func::Uniform => v[0] + (v[1] - v[0]) * uniform(word()),
                Func::Normal => {
                    // Box–Muller with an open first uniform.
                    let mut u1 = uniform(word());
                    while u1 == 0.0 {
                        u1 = uniform(word());
                    }
                    let u2 = uniform(word());
                    let r = libm::sqrt(-2.0 * libm::log(u1));
                    v[0] + v[1] * r * libm::cos(2.0 * std::f64::consts::PI * u2)
                }
                Func::If | Func::NewId | Func::Choice => unreachable!(),
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Num(f64),
    Name(String),
    Sym(&'static str),
}
impl Token {
    fn text(&self) -> String {
        match self {
            Self::Num(x) => format!("{x}"),
            Self::Name(n) => n.clone(),
            Self::Sym(s) => (*s).to_string(),
        }
    }
}

fn lex(source: &str) -> Result<Vec<Token>, String> {
    let chars: Vec<char> = source.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_ascii_whitespace() {
            i += 1;
        } else if c.is_ascii_digit()
            || (c == '.' && chars.get(i + 1).is_some_and(char::is_ascii_digit))
        {
            let start = i;
            while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                i += 1;
            }
            if i < chars.len() && (chars[i] == 'e' || chars[i] == 'E') {
                i += 1;
                if i < chars.len() && (chars[i] == '+' || chars[i] == '-') {
                    i += 1;
                }
                while i < chars.len() && chars[i].is_ascii_digit() {
                    i += 1;
                }
            }
            let text: String = chars[start..i].iter().collect();
            let value: f64 = text
                .parse()
                .map_err(|_| format!("invalid number `{text}`"))?;
            if !value.is_finite() {
                return Err(format!("number `{text}` is not finite"));
            }
            tokens.push(Token::Num(value));
        } else if c.is_ascii_alphabetic() || c == '_' {
            let start = i;
            while i < chars.len()
                && (chars[i].is_ascii_alphanumeric() || chars[i] == '_' || chars[i] == '.')
            {
                i += 1;
            }
            tokens.push(Token::Name(chars[start..i].iter().collect()));
        } else {
            let two: String = chars[i..chars.len().min(i + 2)].iter().collect();
            let sym = ["<=", ">=", "==", "!=", "&&", "||"]
                .into_iter()
                .find(|s| *s == two)
                .or_else(|| {
                    ["+", "-", "*", "/", "^", "<", ">", "!", "(", ")", ","]
                        .into_iter()
                        .find(|s| s.starts_with(c))
                })
                .ok_or_else(|| format!("unexpected character `{c}`"))?;
            i += sym.len();
            tokens.push(Token::Sym(sym));
        }
    }
    if tokens.is_empty() {
        return Err("empty expression".into());
    }
    Ok(tokens)
}

struct Parser<'a, 'b> {
    tokens: Vec<Token>,
    pos: usize,
    scope: &'a Scope<'b>,
}

fn binary(sym: &str) -> Option<(Op, u8, bool)> {
    // (operator, precedence, right associative)
    Some(match sym {
        "||" => (Op::Or, 1, false),
        "&&" => (Op::And, 2, false),
        "==" => (Op::Eq, 3, false),
        "!=" => (Op::Ne, 3, false),
        "<" => (Op::Lt, 4, false),
        "<=" => (Op::Le, 4, false),
        ">" => (Op::Gt, 4, false),
        ">=" => (Op::Ge, 4, false),
        "+" => (Op::Add, 5, false),
        "-" => (Op::Sub, 5, false),
        "*" => (Op::Mul, 6, false),
        "/" => (Op::Div, 6, false),
        "^" => (Op::Pow, 8, true),
        _ => return None,
    })
}

impl Parser<'_, '_> {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }
    fn expect(&mut self, sym: &str) -> Result<(), String> {
        match self.peek() {
            Some(Token::Sym(s)) if *s == sym => {
                self.pos += 1;
                Ok(())
            }
            other => Err(format!(
                "expected `{sym}`, found {}",
                other.map_or("end".into(), |t| format!("`{}`", t.text()))
            )),
        }
    }
    fn expression(&mut self, min: u8) -> Result<Node, String> {
        let mut left = self.unary()?;
        while let Some(Token::Sym(s)) = self.peek() {
            let Some((op, prec, right)) = binary(s) else {
                break;
            };
            if prec < min {
                break;
            }
            self.pos += 1;
            let rhs = self.expression(if right { prec } else { prec + 1 })?;
            left = Node::Bin(op, Box::new(left), Box::new(rhs));
        }
        Ok(left)
    }
    fn unary(&mut self) -> Result<Node, String> {
        match self.peek() {
            Some(Token::Sym("-")) => {
                self.pos += 1;
                // Unary minus binds looser than `^`: -x^2 = -(x^2).
                Ok(Node::Neg(Box::new(self.expression(7)?)))
            }
            Some(Token::Sym("!")) => {
                self.pos += 1;
                Ok(Node::Not(Box::new(self.expression(7)?)))
            }
            _ => self.primary(),
        }
    }
    fn primary(&mut self) -> Result<Node, String> {
        let token = self
            .tokens
            .get(self.pos)
            .cloned()
            .ok_or("unexpected end of expression")?;
        self.pos += 1;
        match token {
            Token::Num(x) => Ok(Node::Num(x)),
            Token::Sym("(") => {
                let inner = self.expression(0)?;
                self.expect(")")?;
                Ok(inner)
            }
            Token::Name(name) => {
                if matches!(self.peek(), Some(Token::Sym("("))) {
                    let (f, arity, random) =
                        Func::parse(&name).ok_or_else(|| format!("unknown function `{name}`"))?;
                    if random && !self.scope.random {
                        return Err(format!(
                            "random function `{name}` is only allowed in product attributes"
                        ));
                    }
                    self.pos += 1;
                    let mut args = Vec::new();
                    if !matches!(self.peek(), Some(Token::Sym(")"))) {
                        loop {
                            args.push(self.expression(0)?);
                            if matches!(self.peek(), Some(Token::Sym(","))) {
                                self.pos += 1;
                            } else {
                                break;
                            }
                        }
                    }
                    self.expect(")")?;
                    if args.len() != arity {
                        return Err(format!("`{name}` takes {arity} arguments"));
                    }
                    return Ok(Node::Call(f, args));
                }
                self.resolve(&name).map(Node::Var)
            }
            Token::Sym(s) => Err(format!("unexpected `{s}`")),
        }
    }
    fn resolve(&self, name: &str) -> Result<Slot, String> {
        if let Some((label, attribute)) = name.split_once('.') {
            for (reactant, (l, attributes)) in self.scope.reactants.iter().enumerate() {
                if *l == label {
                    return attributes
                        .iter()
                        .position(|a| a == attribute)
                        .map(|index| Slot::Attribute { reactant, index })
                        .ok_or_else(|| format!("`{label}` has no attribute `{attribute}`"));
                }
            }
            return Err(format!("unknown reactant label `{label}`"));
        }
        if let Some(&i) = self.scope.params.get(name) {
            return Ok(Slot::Param(i));
        }
        if let Some(&i) = self.scope.observables.and_then(|o| o.get(name)) {
            return Ok(Slot::Observable(i));
        }
        let allowed = if self.scope.observables.is_some() {
            "parameter or observable"
        } else {
            "parameter (or `label.attribute`)"
        };
        Err(format!("unknown name `{name}`; expected a {allowed}"))
    }
}

impl Expr {
    /// Whether all inputs are available without a reactant or a random draw.
    pub(crate) fn parameter_only(&self) -> bool {
        fn check(node: &Node) -> bool {
            match node {
                Node::Num(_) | Node::Var(Slot::Param(_)) => true,
                Node::Var(_) => false,
                Node::Neg(a) | Node::Not(a) => check(a),
                Node::Bin(_, a, b) => check(a) && check(b),
                Node::Call(f, args) => {
                    !matches!(
                        f,
                        Func::Bernoulli | Func::Choice | Func::Uniform | Func::Normal | Func::NewId
                    ) && args.iter().all(check)
                }
            }
        }
        check(&self.root)
    }

    pub(crate) fn deterministic(&self) -> bool {
        fn check(node: &Node) -> bool {
            match node {
                Node::Num(_) | Node::Var(_) => true,
                Node::Neg(a) | Node::Not(a) => check(a),
                Node::Bin(_, a, b) => check(a) && check(b),
                Node::Call(f, args) => {
                    !matches!(
                        f,
                        Func::Bernoulli | Func::Choice | Func::Uniform | Func::Normal | Func::NewId
                    ) && args.iter().all(check)
                }
            }
        }
        check(&self.root)
    }

    pub(crate) fn reads(&self, flags: &mut [bool]) {
        fn visit(node: &Node, flags: &mut [bool]) {
            match node {
                Node::Var(Slot::Observable(i)) => flags[*i] = true,
                Node::Neg(a) | Node::Not(a) => visit(a, flags),
                Node::Bin(_, a, b) => {
                    visit(a, flags);
                    visit(b, flags);
                }
                Node::Call(_, args) => {
                    for arg in args {
                        visit(arg, flags);
                    }
                }
                Node::Num(_) | Node::Var(_) => {}
            }
        }
        visit(&self.root, flags);
    }

    pub(crate) fn kernel(&self, params: &[f64]) -> RateKernel {
        if self.parameter_only() {
            return RateKernel::Constant(self.eval(&Env {
                params,
                observables: &[],
                attributes: [&[], &[]],
            }));
        }
        let mut root = &self.root;
        let mut factors = Vec::new();
        while let Node::Bin(Op::Mul, left, right) = root {
            let right = Self {
                root: (**right).clone(),
            };
            if !right.parameter_only() {
                break;
            }
            factors.push(right.eval(&Env {
                params,
                observables: &[],
                attributes: [&[], &[]],
            }));
            root = left;
        }
        let vacancy = || {
            let Node::Bin(Op::Div, num, den) = root else {
                return None;
            };
            let Node::Var(Slot::Observable(d)) = den.as_ref() else {
                return None;
            };
            let Node::Call(Func::Max, args) = num.as_ref() else {
                return None;
            };
            if args[0] != Node::Num(0.0) {
                return None;
            }
            let Node::Bin(Op::Sub, k, n) = &args[1] else {
                return None;
            };
            let Node::Var(Slot::Observable(n)) = n.as_ref() else {
                return None;
            };
            let limit = Self {
                root: (**k).clone(),
            };
            if !limit.parameter_only() {
                return None;
            }
            let limit = limit.eval(&Env {
                params,
                observables: &[],
                attributes: [&[], &[]],
            });
            if !limit.is_finite() {
                return None;
            }
            factors.reverse();
            Some(RateKernel::Vacancy {
                limit,
                count: *n,
                denominator: *d,
                factors,
            })
        };
        if let Some(kernel) = vacancy() {
            return kernel;
        }
        RateKernel::Expression(self.clone())
    }
}

pub(crate) enum RateKernel {
    Constant(f64),
    Vacancy {
        limit: f64,
        count: usize,
        denominator: usize,
        factors: Vec<f64>,
    },
    Expression(Expr),
}

impl RateKernel {
    pub(crate) fn eval(&self, env: &Env<'_>) -> f64 {
        match self {
            Self::Constant(value) => *value,
            Self::Vacancy {
                limit,
                count,
                denominator,
                factors,
            } => {
                let base =
                    (limit - env.observables[*count]).max(0.0) / env.observables[*denominator];
                factors.iter().fold(base, |x, f| x * f)
            }
            Self::Expression(expr) => expr.eval(env),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vacancy_kernel_keeps_undefined_limit_undefined() {
        let params = BTreeMap::new();
        let names = BTreeMap::from([("N".into(), 0), ("D".into(), 1)]);
        let scope = Scope {
            params: &params,
            observables: Some(&names),
            reactants: &[],
            random: false,
        };
        let expr = Expr::parse("max(0, log(-1) - N) / D", &scope).unwrap();
        let env = Env {
            params: &[],
            observables: &[1.0, 1.0],
            attributes: [&[], &[]],
        };
        assert!(expr.eval(&env).is_nan());
        assert!(expr.kernel(&[]).eval(&env).is_nan());
    }

    fn eval(source: &str) -> f64 {
        let params = BTreeMap::from([("k".to_string(), 0usize)]);
        let attrs = vec!["len".to_string(), "s".to_string()];
        let reactants = [("x", attrs.as_slice())];
        let scope = Scope {
            params: &params,
            observables: None,
            reactants: &reactants,
            random: false,
        };
        let e = Expr::parse(source, &scope).unwrap();
        e.eval(&Env {
            params: &[2.0],
            observables: &[],
            attributes: [&[3.0, 0.5], &[]],
        })
    }

    #[test]
    fn precedence_and_functions() {
        assert_eq!(eval("1 + 2 * 3"), 7.0);
        assert_eq!(eval("-2 ^ 2"), -4.0);
        assert_eq!(eval("2 ^ 3 ^ 2"), 512.0);
        assert_eq!(eval("k * x.len - x.s"), 5.5);
        assert_eq!(eval("x.len >= 3 && x.s < 1"), 1.0);
        assert_eq!(eval("!(x.len == 3) || 0"), 0.0);
        assert_eq!(eval("if(x.len > 2, max(1, k), 7)"), 2.0);
        assert_eq!(eval("0 && 1 / 0"), 0.0);
        assert_eq!(eval("1.5e1 / 3"), 5.0);
        // Undefined conditions stay undefined; min/max do not drop NaN.
        assert!(eval("(0/0) && 1").is_nan());
        assert!(eval("1 || (0/0)") == 1.0 && eval("0 || (0/0)").is_nan());
        assert!(eval("!(0/0)").is_nan() && eval("if(0/0, 1, 2)").is_nan());
        assert!(eval("max(0/0, 1)").is_nan() && eval("min(1, 0/0)").is_nan());
    }

    #[test]
    fn names_and_random_functions_are_checked() {
        let params = BTreeMap::from([("k".to_string(), 0usize)]);
        let observables = BTreeMap::from([("N".to_string(), 0usize)]);
        let attrs = vec!["len".to_string()];
        let reactants = [("x", attrs.as_slice())];
        let local = Scope {
            params: &params,
            observables: None,
            reactants: &reactants,
            random: false,
        };
        assert!(
            Expr::parse("N", &local)
                .unwrap_err()
                .contains("unknown name")
        );
        assert!(
            Expr::parse("x.size", &local)
                .unwrap_err()
                .contains("no attribute")
        );
        assert!(Expr::parse("y.len", &local).unwrap_err().contains("label"));
        assert!(
            Expr::parse("uniform(0, 1)", &local)
                .unwrap_err()
                .contains("random")
        );
        assert!(
            Expr::parse("max(1)", &local)
                .unwrap_err()
                .contains("arguments")
        );
        assert!(Expr::parse("1 +", &local).is_err());
        assert!(Expr::parse("(1", &local).is_err());
        assert!(Expr::parse("1e999", &local).is_err());
        assert!(Expr::parse("", &local).is_err());
        let global = Scope {
            params: &params,
            observables: Some(&observables),
            reactants: &[],
            random: false,
        };
        assert!(Expr::parse("k * N", &global).is_ok());
        assert!(Expr::parse("x.len", &global).is_err());
    }

    #[test]
    fn random_draws_and_new_ids() {
        let params = BTreeMap::new();
        let scope = Scope {
            params: &params,
            observables: None,
            reactants: &[],
            random: true,
        };
        let id = Expr::parse("new_id()", &scope).unwrap();
        let coin = Expr::parse("choice(0.25, 10, 20)", &scope).unwrap();
        let mut rng = crate::stream(1, 2);
        let (mut next, mut words, mut overflow) = (5.0, 0, false);
        let env = Env {
            params: &[],
            observables: &[],
            attributes: [&[], &[]],
        };
        let mut draw = Draw {
            rng: &mut rng,
            next_id: &mut next,
            words: &mut words,
            id_overflow: &mut overflow,
        };
        assert_eq!(id.eval_random(&env, &mut draw), 5.0);
        assert_eq!(id.eval_random(&env, &mut draw), 6.0);
        let tens = (0..4000)
            .filter(|_| coin.eval_random(&env, &mut draw) == 10.0)
            .count();
        assert!((800..1200).contains(&tens), "{tens}");
        assert_eq!(words, 4000);
    }
}
