use crn_craft::{kinetics::Kinetics, parse_model};
#[test]
fn analytic_state_and_parameter_derivatives_match_independent_central_differences() {
    for bytes in [
        include_bytes!("../../../models/examples/repressilator.model.json").as_slice(),
        include_bytes!("../../../models/examples/enzyme-mechanism.model.json").as_slice(),
        include_bytes!("../../../models/examples/michaelis-menten.model.json").as_slice(),
    ] {
        let model = parse_model(bytes).unwrap();
        let k = Kinetics::compile(&model).unwrap();
        let p: Vec<_> = k
            .parameters
            .iter()
            .map(|name| model.parameters[name])
            .collect();
        // Distinct components: a derivative taken with respect to the wrong species must not
        // coincide with the right one.
        let n = model.species.len();
        let distinct = [
            (0..n).map(|i| 3.7 + 2.9 * i as f64).collect::<Vec<_>>(),
            (0..n).map(|i| 41.3 - 3.1 * i as f64).collect(),
        ];
        for y in distinct {
            let result = k.evaluate(&y, &p).unwrap();
            for j in 0..y.len() {
                let h = 1e-4;
                let mut plus = y.clone();
                let mut minus = y.clone();
                plus[j] += h;
                minus[j] -= h;
                let a = k.evaluate(&plus, &p).unwrap().rhs;
                let b = k.evaluate(&minus, &p).unwrap().rhs;
                for i in 0..y.len() {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    assert!(
                        (numeric - result.state_jacobian[i][j]).abs()
                            < 1e-7 * numeric.abs().max(1.)
                    );
                }
            }
            for j in 0..p.len() {
                let h = p[j] * 1e-5;
                let mut plus = p.clone();
                let mut minus = p.clone();
                plus[j] += h;
                minus[j] -= h;
                let a = k.evaluate(&y, &plus).unwrap().rhs;
                let b = k.evaluate(&y, &minus).unwrap().rhs;
                for i in 0..y.len() {
                    let numeric = (a[i] - b[i]) / (2. * h);
                    assert!(
                        (numeric - result.parameter_jacobian[i][j]).abs()
                            < 1e-6 * numeric.abs().max(1.)
                    );
                }
            }
        }
        assert!(k.evaluate(&vec![-1.; model.species.len()], &p).is_err());
    }
}
