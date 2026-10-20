use thermopack_rs::fmt::{FmtDensity, FmtModel};

const MODELS: [FmtModel; 3] = [
    FmtModel::Rosenfeld,
    FmtModel::WhiteBear,
    FmtModel::WhiteBearII,
];

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn values(s: FmtDensity) -> Vec<f64> {
    [s.value]
        .into_iter()
        .chain(s.gradient)
        .chain(s.hessian.into_iter().flatten())
        .collect()
}

#[test]
fn source_functionals_and_optional_output_branches_match() {
    let mut count = 0;
    for row in include_str!("fixtures/fmt.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let model = FmtModel::from_code(f[1].parse().unwrap()).unwrap();
        let n: Vec<f64> = f[2].split(',').map(|v| v.parse().unwrap()).collect();
        let s = model
            .baseline_energy_density(n.try_into().unwrap())
            .unwrap();
        let flags: u8 = f[3].parse().unwrap();
        for (i, a) in values(s).into_iter().enumerate() {
            let active = if i == 0 {
                flags & 1 != 0
            } else if i < 7 {
                flags & 2 != 0
            } else {
                flags & 4 != 0
            };
            if active {
                close(a, c[i + 1].parse().unwrap(), 1e-10, 2e-10);
            } else {
                assert_eq!(c[i + 1].parse::<f64>().unwrap(), 0.);
            }
        }
        count += 1;
    }
    assert_eq!(count, 1288);
}

#[test]
fn analytic_gradients_and_hessians_close_independent_differences() {
    for model in MODELS {
        for eta in [0.01, 0.2, 0.7] {
            let n = [0.1, 0.2, 0.3, eta, 0.02, -0.04];
            let s = model.energy_density(n).unwrap();
            for i in 0..6 {
                let delta = 1e-5;
                let mut hi = n;
                hi[i] += delta;
                let mut lo = n;
                lo[i] -= delta;
                let a = model.energy_density(hi).unwrap();
                let b = model.energy_density(lo).unwrap();
                close(
                    s.gradient[i],
                    (a.value - b.value) / (2. * delta),
                    1e-8,
                    2e-5,
                );
                for j in 0..6 {
                    assert_eq!(s.hessian[i][j], s.hessian[j][i]);
                    close(
                        s.hessian[j][i],
                        (a.gradient[j] - b.gradient[j]) / (2. * delta),
                        1e-7,
                        2e-4,
                    );
                }
            }
            if eta < 1. {
                let reference = model.baseline_energy_density(n).unwrap();
                for (a, b) in values(s).into_iter().zip(values(reference)) {
                    close(a, b, 1e-9, 1e-8);
                }
            }
        }
    }
}

#[test]
fn uniform_sphere_limit_has_the_expected_excess_energy() {
    for eta in [0.01, 0.1, 0.4, 0.8] {
        let rho = 3. * eta / (4. * std::f64::consts::PI);
        let n = [rho, rho, 3. * eta, eta, 0., 0.];
        let u = 1. - eta;
        for model in MODELS {
            let excess = match model {
                FmtModel::Rosenfeld => -u.ln() + 3. * eta / u + 1.5 * eta * eta / (u * u),
                _ => (4. * eta - 3. * eta * eta) / (u * u),
            };
            close(
                model.energy_density(n).unwrap().value,
                rho * excess,
                1e-12,
                1e-11,
            );
        }
    }
}

#[test]
fn grids_invalid_inputs_and_source_derivative_floor_are_explicit() {
    for model in MODELS {
        let n = [0.1, 0.2, 0.3, 0.2, 0.02, -0.04];
        let initial = model.energy_density(n).unwrap();
        let inputs = [n, [0.2, 0.3, 0.4, 0.3, 0., 0.]];
        let mut out = [initial; 2];
        model.energy_density_grid(&inputs, &mut out).unwrap();
        for (s, n) in out.into_iter().zip(inputs) {
            assert_eq!(values(s), values(model.energy_density(n).unwrap()));
        }
        assert!(
            model
                .energy_density_grid(&inputs, &mut [initial; 1])
                .is_err()
        );
        model.energy_density_grid(&[], &mut []).unwrap();
        model
            .baseline_energy_density_grid(&inputs, &mut out)
            .unwrap();
        for (s, n) in out.into_iter().zip(inputs) {
            assert_eq!(values(s), values(model.baseline_energy_density(n).unwrap()));
        }
        model.baseline_energy_density_grid(&[], &mut []).unwrap();
        assert!(
            model
                .baseline_energy_density_grid(&inputs, &mut [initial; 1])
                .is_err()
        );
        let mut invalid = inputs;
        invalid[1][0] = f64::NAN;
        for baseline in [false, true] {
            let mut partial = [initial; 2];
            let mut first = inputs[1];
            first[3] = 0.4;
            invalid[0] = first;
            let error = if baseline {
                model.baseline_energy_density_grid(&invalid, &mut partial)
            } else {
                model.energy_density_grid(&invalid, &mut partial)
            };
            assert!(error.is_err());
            let expected = if baseline {
                model.baseline_energy_density(first)
            } else {
                model.energy_density(first)
            }
            .unwrap();
            assert_eq!(values(partial[0]), values(expected));
            assert_eq!(values(partial[1]), values(initial));
        }
        for i in 0..6 {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut invalid = n;
                invalid[i] = bad;
                assert!(model.energy_density(invalid).is_err());
                assert!(model.baseline_energy_density(invalid).is_err());
            }
        }
    }
    for model in [FmtModel::WhiteBear, FmtModel::WhiteBearII] {
        assert!(model.energy_density([0.; 6]).is_err());
        assert!(model.baseline_energy_density([0.; 6]).is_err());
    }
    assert_eq!(
        FmtModel::Rosenfeld.energy_density([0.; 6]).unwrap().value,
        0.
    );
    assert!(FmtModel::from_code(0).is_err());
    assert!(FmtModel::from_code(4).is_err());
    for eta in [1., 1.1] {
        let s = FmtModel::Rosenfeld
            .energy_density([0.1, 0.2, 0.3, eta, 0.02, -0.04])
            .unwrap();
        assert_eq!(s.gradient[3], 0.);
        assert!(s.hessian[3].iter().all(|v| *v == 0.));
    }
    let n = [0.1, 0.2, 0.3, 1. - 1e-9, 0.02, -0.04];
    let analytic = FmtModel::Rosenfeld.energy_density(n).unwrap();
    let baseline = FmtModel::Rosenfeld.baseline_energy_density(n).unwrap();
    close(analytic.value, baseline.value, 1e-10, 1e-14);
    assert!((analytic.gradient[3] / baseline.gradient[3]).abs() > 10.);
}
