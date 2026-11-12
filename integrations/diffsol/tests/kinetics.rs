//! Analytic Jacobians of the core deterministic kinetics at states whose components all differ,
//! so that a derivative taken with respect to the wrong species cannot coincide with the right one.
use crn_craft::kinetics::Kinetics;

const MODELS: [&str; 7] = [
    include_str!("../../../models/legacy/repressilator.model.json"),
    include_str!("../../../models/legacy/enzyme-mechanism.model.json"),
    include_str!("../../../models/legacy/michaelis-menten.model.json"),
    include_str!("../../../models/legacy/vilar-noise.model.json"),
    include_str!("../fixtures/dimerization.model.json"),
    include_str!("../fixtures/hill-constant.model.json"),
    include_str!("../fixtures/hill-dynamic.model.json"),
];

/// Three distinct, well-separated states per model; the third puts a regulator near K = 10.
fn states(n: usize) -> Vec<Vec<f64>> {
    vec![
        (0..n).map(|i| 3.7 + 2.9 * i as f64).collect(),
        (0..n).map(|i| 41.3 - 3.1 * i as f64).collect(),
        (0..n).map(|i| 9.6 + 0.7 * (i * i) as f64).collect(),
    ]
}

#[test]
fn state_and_parameter_jacobians_match_central_differences_at_distinct_states() {
    for (label, text) in MODELS.iter().enumerate() {
        // Any model version: the kinetics read the model through its flat view.
        let kinetics = Kinetics::from_bytes(text.as_bytes()).unwrap();
        let p = kinetics.parameter_values().to_vec();
        let n = kinetics.species().len();
        for y in states(n) {
            assert!(y.iter().all(|v| *v > 0.5));
            let exact = kinetics.evaluate(&y, &p).unwrap();
            for j in 0..n {
                let h = 1e-5 * y[j];
                let (mut plus, mut minus) = (y.clone(), y.clone());
                plus[j] += h;
                minus[j] -= h;
                let a = kinetics.evaluate(&plus, &p).unwrap().rhs;
                let b = kinetics.evaluate(&minus, &p).unwrap().rhs;
                for i in 0..n {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    let analytic = exact.state_jacobian[i][j];
                    assert!(
                        (numeric - analytic).abs() <= 1e-6 * numeric.abs().max(1.),
                        "model {label} d rhs[{i}] / d y[{j}] at {y:?}: {analytic} vs {numeric}"
                    );
                }
            }
            for j in 0..p.len() {
                let h = 1e-5 * p[j];
                let (mut plus, mut minus) = (p.clone(), p.clone());
                plus[j] += h;
                minus[j] -= h;
                let a = kinetics.evaluate(&y, &plus).unwrap().rhs;
                let b = kinetics.evaluate(&y, &minus).unwrap().rhs;
                for i in 0..n {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    let analytic = exact.parameter_jacobian[i][j];
                    assert!(
                        (numeric - analytic).abs() <= 1e-6 * numeric.abs().max(1.),
                        "model {label} d rhs[{i}] / d {} at {y:?}: {analytic} vs {numeric}",
                        kinetics.parameters[j]
                    );
                }
            }
        }
    }
}

#[test]
fn mass_action_right_hand_side_uses_powers_over_factorials() {
    // X + 2Y -> Z at rate c X Y^2 / 2!; dimerization 2A -> B at rate c A^2 / 2!.
    let kinetics = Kinetics::from_bytes(MODELS[6].as_bytes()).unwrap();
    let p = kinetics.parameter_values().to_vec();
    let y = [7.0, 3.0, 5.0, 11.0]; // R, X, Y, Z
    let rhs = kinetics.evaluate(&y, &p).unwrap().rhs;
    let c = p[kinetics
        .parameters
        .iter()
        .position(|name| name == "c")
        .unwrap()];
    assert!((rhs[3] - c * 3.0 * 25.0 / 2.0).abs() <= 1e-12, "{rhs:?}");
    let kinetics = Kinetics::from_bytes(MODELS[4].as_bytes()).unwrap();
    let rhs = kinetics.evaluate(&[12.0, 0.0], &[0.01]).unwrap().rhs;
    assert!((rhs[0] + 0.01 * 144.0).abs() <= 1e-12 && (rhs[1] - 0.01 * 72.0).abs() <= 1e-12);
}
