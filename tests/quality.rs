use thermopack_rs::envelope::EnvelopeConstraint;
use thermopack_rs::quality::QualityGuess;
use thermopack_rs::{Cubic, CubicKind};

#[test]
fn quality_states_close_fugacity_material_balance_and_branch_slopes() {
    for eos in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(eos, &["C1", "C2", "C3"]).unwrap();
        let z = [0.4, 0.6, 0.];
        for beta in [0., 0.01, 0.25, 0.5, 0.75, 0.99, 1.] {
            let a = model.quality_temperature(1e6, &z, beta).unwrap();
            assert!(a.fugacity_residual < 1e-9);
            assert!(a.material_balance_residual < 1e-9);
            assert_eq!(a.vapor_fraction, beta);
            assert_eq!(a.feed, z);
            assert_eq!(a.liquid[2], 0.);
            assert_eq!(a.vapor[2], 0.);
            assert!(a.liquid_volume < a.vapor_volume);
            let l = model.prepare(a.temperature, &a.liquid).unwrap();
            let v = model.prepare(a.temperature, &a.vapor).unwrap();
            let mut fl = vec![0.; 3];
            let mut fv = fl.clone();
            l.ln_fugacity_tp_into(a.pressure, a.liquid_volume, &mut fl)
                .unwrap();
            v.ln_fugacity_tp_into(a.pressure, a.vapor_volume, &mut fv)
                .unwrap();
            for i in 0..2 {
                assert!(((1. - beta) * a.liquid[i] + beta * a.vapor[i] - z[i]).abs() < 1e-9);
                assert!((a.liquid[i].ln() + fl[i] - a.vapor[i].ln() - fv[i]).abs() < 1e-9);
            }
            for (state, vol) in [(l, a.liquid_volume), (v, a.vapor_volume)] {
                assert!((state.pressure(vol).unwrap().value - a.pressure).abs() < 0.1);
            }
            let h = 1e-3;
            let plus = model
                .refine_quality_state(&a, EnvelopeConstraint::Temperature(a.temperature + h))
                .unwrap();
            let minus = model
                .refine_quality_state(&a, EnvelopeConstraint::Temperature(a.temperature - h))
                .unwrap();
            let fd = (plus.pressure - minus.pressure) / (2. * h);
            assert!(
                (fd - a.dpdt).abs() < 1e-3 + 1e-6 * a.dpdt.abs(),
                "{eos:?} beta={beta}: {fd} vs {}",
                a.dpdt
            );
            let b = model.quality_pressure(a.temperature, &z, beta).unwrap();
            assert!((a.pressure - b.pressure).abs() < 0.1);
            if beta == 0. {
                assert!(
                    (model.bubble_temperature(1e6, &z).unwrap().temperature - a.temperature).abs()
                        < 1e-6
                );
            } else if beta == 1. {
                assert!(
                    (model.dew_temperature(1e6, &z).unwrap().temperature - a.temperature).abs()
                        < 1e-6
                );
            }
        }
    }
}

#[test]
fn embedded_pure_quality_and_invalid_inputs_have_explicit_behavior() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let a = model.quality_temperature(1e6, &[0., 1.], 0.3).unwrap();
    assert_eq!(a.liquid, [0., 1.]);
    assert_eq!(a.vapor, [0., 1.]);
    let sat = model.pseudo_saturation_temperature(1e6, &[0., 1.]).unwrap();
    assert_eq!(a.temperature, sat.temperature);
    assert_eq!(a.dpdt, sat.dpdt);
    for beta in [-1., 1.01, f64::NAN, f64::INFINITY] {
        assert!(model.quality_temperature(1e6, &[0.5, 0.5], beta).is_err());
    }
    for z in [vec![], vec![0., 0.], vec![1., -0.1], vec![f64::NAN, 1.]] {
        assert!(model.quality_pressure(200., &z, 0.5).is_err());
    }
    let guess = QualityGuess {
        temperature: 500.,
        pressure: 1e6,
        ratios: vec![1., 1.],
    };
    assert!(
        model
            .quality_state_with_guess(&[0.5, 0.5], 0.5, EnvelopeConstraint::Pressure(1e6), &guess)
            .is_err()
    );
    assert!(model.quality_temperature(-1., &[0.5, 0.5], 0.5).is_err());
}

#[test]
fn source_quality_states_retain_strict_matches_and_declared_fidelity_gaps() {
    use thermopack_rs::{CubicOptions, MixingRule};
    let mut total = 0;
    let mut matching = 0;
    let mut cold_matches = 0;
    let mut cold_gaps = 0;
    for row in include_str!("fixtures/quality.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let names: Vec<_> = f[2].split(',').collect();
        let model = if f[3] == "quantum" {
            Cubic::quantum(&names).unwrap()
        } else {
            Cubic::with_options(
                match f[1] {
                    "SRK" => CubicKind::SoaveRedlichKwong,
                    "RK" => CubicKind::RedlichKwong,
                    "VDW" => CubicKind::VanDerWaals,
                    _ => CubicKind::PengRobinson,
                },
                &names,
                CubicOptions {
                    mixing_rule: match f[3] {
                        "unifac" => MixingRule::Unifac,
                        "hv1" => MixingRule::HuronVidal1,
                        "hv2" => MixingRule::HuronVidal2,
                        "nrtl" => MixingRule::Nrtl,
                        _ => MixingRule::VdW,
                    },
                    volume_translation: f[3] == "translated",
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let z: Vec<f64> = f[4].split(',').map(|v| v.parse().unwrap()).collect();
        let beta = f[5].parse::<f64>().unwrap();
        let o: Vec<_> = f[6].split(':').collect();
        let value = o[1].parse::<f64>().unwrap();
        let constraint = match o[0] {
            "T" => EnvelopeConstraint::Temperature(value),
            "P" => EnvelopeConstraint::Pressure(value),
            _ => EnvelopeConstraint::LogK {
                component: o[2].parse().unwrap(),
                value,
            },
        };
        let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
        let guess = QualityGuess {
            temperature: g[0],
            pressure: g[1],
            ratios: g[2..].to_vec(),
        };
        let a = model
            .baseline_quality_state_with_guess(&z, beta, constraint, &guess)
            .unwrap();
        let expected: Vec<f64> = c[1].split(',').map(|v| v.parse().unwrap()).collect();
        let encode = |a: &thermopack_rs::quality::QualityState| {
            let mut v = vec![
                a.temperature,
                a.pressure,
                a.vapor_fraction,
                a.liquid_volume,
                a.vapor_volume,
            ];
            v.extend(&a.liquid);
            v.extend(&a.vapor);
            v.push(a.dpdt);
            v
        };
        let compare = |actual: Vec<f64>, allow_pending: bool| {
            let mut misses = 0;
            for (j, (&a, b)) in expected.iter().zip(actual).enumerate() {
                let floor = if j == 0 {
                    2e-6
                } else if j == 1 {
                    0.1
                } else if j == 2 {
                    0.
                } else if j < 5 {
                    1e-12
                } else if j == expected.len() - 1 {
                    1e-3
                } else {
                    2e-8
                };
                let rtol = if j == expected.len() - 1 { 2e-6 } else { 2e-7 };
                if (a - b).abs() > floor + rtol * a.abs().max(b.abs()) {
                    assert!(
                        allow_pending && [6, 8].contains(&j) && (a - b).abs() < 1.2e-7,
                        "{} {j}: {a} vs {b}",
                        c[0]
                    );
                    misses += 1;
                }
                assert!(b.is_finite());
            }
            misses
        };
        let pending = c[7] == "0";
        if pending {
            assert_eq!(f[2], "C1,H2O");
            assert!(
                a.liquid
                    .iter()
                    .zip(&a.vapor)
                    .all(|(x, y)| (x - y).abs() < 1e-4)
            );
        }
        let misses = compare(encode(&a), pending);
        if misses == 0 {
            matching += 1;
        }
        assert!(a.fugacity_residual < 1e-8);
        assert!(a.material_balance_residual < 1e-8);
        if c[7] != "-1" {
            let cold = match constraint {
                EnvelopeConstraint::Temperature(t) => model.quality_pressure(t, &z, beta),
                EnvelopeConstraint::Pressure(p) => model.quality_temperature(p, &z, beta),
                _ => unreachable!(),
            };
            if c[7] == "1" {
                assert_eq!(compare(encode(&cold.unwrap()), false), 0);
                cold_matches += 1;
            } else {
                assert!(cold.is_err());
                cold_gaps += 1;
            }
        }
        total += 1;
    }
    assert_eq!((total, cold_matches, cold_gaps), (577, 491, 8));
    // Eight nearly coincident states remain a conditioning/fidelity audit.
    // All other 569 states must match strictly on every build. The release
    // qualification records its additional strict matches individually.
    assert!((569..=577).contains(&matching));
}
