use thermopack_rs::quality::{QualityGuess, QualitySpecification};
use thermopack_rs::{Cubic, CubicKind};

fn decode(case: &str) -> (Cubic, Vec<f64>, f64, QualitySpecification, QualityGuess) {
    use thermopack_rs::{CubicOptions, MixingRule};
    let f: Vec<_> = case.split('|').collect();
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
    let z = f[4].split(',').map(|v| v.parse().unwrap()).collect();
    let o: Vec<_> = f[6].split(':').collect();
    let target = o[1].parse().unwrap();
    let spec = if o[0] == "H" {
        QualitySpecification::Enthalpy(target)
    } else {
        QualitySpecification::Entropy(target)
    };
    let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
    (
        model,
        z,
        f[5].parse().unwrap(),
        spec,
        QualityGuess {
            temperature: g[0],
            pressure: g[1],
            ratios: g[2..].to_vec(),
        },
    )
}

fn check_target(
    model: &Cubic,
    state: &thermopack_rs::quality::QualityState,
    spec: QualitySpecification,
) {
    assert!(state.fugacity_residual < 1e-9);
    assert!(state.material_balance_residual < 1e-9);
    let p = model.quality_properties(state).unwrap();
    let (actual, target) = match spec {
        QualitySpecification::Enthalpy(h) => (p.enthalpy, h),
        QualitySpecification::Entropy(s) => (p.entropy, s),
        _ => unreachable!(),
    };
    assert!((actual - target).abs() < 1.1e-9 * target.abs().max(1.));
}

#[test]
fn source_caloric_quality_states_match_and_three_physical_differences_stay_explicit() {
    let mut total = 0;
    let mut physical_matches = 0;
    let mut gaps = Vec::new();
    for row in include_str!("fixtures/quality_caloric.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let (model, z, beta, spec, guess) = decode(c[0]);
        let baseline = model
            .baseline_quality_state_with_specification(&z, beta, spec, &guess)
            .unwrap();
        let physical = model
            .quality_state_with_specification(&z, beta, spec, &guess)
            .unwrap();
        check_target(&model, &physical, spec);
        let expected: Vec<f64> = c[1].split(',').map(|v| v.parse().unwrap()).collect();
        for (state, physical_mode) in [(&baseline, false), (&physical, true)] {
            assert_eq!(state.feed, z);
            let p = model.quality_properties(state).unwrap();
            let mut values = vec![
                state.temperature,
                state.pressure,
                beta,
                state.liquid_volume,
                state.vapor_volume,
            ];
            values.extend(&state.liquid);
            values.extend(&state.vapor);
            values.extend([state.dpdt, p.enthalpy, p.entropy]);
            let mut misses = Vec::new();
            for (j, (&a, b)) in expected.iter().zip(values).enumerate() {
                let floor = if j == 0 {
                    2e-5
                } else if j == 1 {
                    2.
                } else if j == 2 {
                    0.
                } else if j < 5 {
                    1e-10
                } else if j == 5 + 2 * z.len() {
                    1e-3
                } else if j == 6 + 2 * z.len() {
                    1e-5
                } else if j == 7 + 2 * z.len() {
                    1e-8
                } else {
                    2e-8
                };
                let rtol = if physical_mode || j == 5 + 2 * z.len() {
                    2e-6
                } else {
                    2e-7
                };
                assert!(b.is_finite());
                if (a - b).abs() > floor + rtol * a.abs().max(b.abs()) {
                    assert!(physical_mode, "{} {j}: {a} vs {b}", c[0]);
                    match (f[1], f[2], f[3], f[5], f[6].split(':').next().unwrap()) {
                        ("RK", "C1,C3", "classic", "0.75", "S") => {
                            assert!([4, 9].contains(&j));
                            assert!((a - b).abs() < 5e-6 * a.abs().max(b.abs()));
                        }
                        ("PR", "C1,C2,C3", "unifac", "0.75", "S") => {
                            assert_eq!(j, 4);
                            assert!((a - b).abs() < 3e-6 * a.abs().max(b.abs()));
                        }
                        ("SRK", "C1,H2O", "hv1", "0.25", "H") => {
                            assert!([0, 1, 3, 4, 9, 11].contains(&j));
                            assert!((state.temperature - 135.9567589314).abs() < 2e-5);
                            assert!((state.pressure - 1e6).abs() < 2.);
                        }
                        _ => panic!("unexpected physical difference {} {j}: {a} vs {b}", c[0]),
                    }
                    misses.push(j);
                }
            }
            if physical_mode {
                if misses.is_empty() {
                    physical_matches += 1;
                } else {
                    gaps.push((c[0], misses));
                }
            }
        }
        total += 1;
    }
    assert_eq!((total, physical_matches, gaps.len()), (203, 200, 3));
}

#[test]
fn physical_caloric_solver_recovers_all_seventy_seven_source_solver_failures() {
    let successful: std::collections::HashSet<_> = include_str!("fixtures/quality_caloric.csv")
        .lines()
        .skip(1)
        .map(|v| v.split(';').next().unwrap())
        .collect();
    let mut count = 0;
    for row in include_str!("fixtures/quality_caloric_jacobian.csv")
        .lines()
        .skip(1)
    {
        let case = row.split(';').next().unwrap();
        if successful.contains(case) {
            continue;
        }
        let (model, z, beta, spec, guess) = decode(case);
        assert!(
            model
                .baseline_quality_state_with_specification(&z, beta, spec, &guess)
                .is_err(),
            "{case}"
        );
        let state = model
            .quality_state_with_specification(&z, beta, spec, &guess)
            .unwrap();
        check_target(&model, &state, spec);
        assert_eq!(state.feed, z);
        count += 1;
    }
    assert_eq!(count, 77);
}

#[test]
fn caloric_reference_shifts_and_zero_targets_keep_the_same_quality_state() {
    use thermopack_rs::EntropyReference;
    for entropy in [false, true] {
        let mut model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
        let reference = model.quality_temperature(1e6, &[0.5, 0.5], 0.25).unwrap();
        let props = model.quality_properties(&reference).unwrap();
        for i in 0..2 {
            let ideal = model.ideal_gas_mut(i).unwrap();
            if entropy {
                ideal
                    .set_standard_entropy(
                        ideal.standard_entropy() - props.entropy,
                        EntropyReference::Bar,
                    )
                    .unwrap();
            } else {
                ideal
                    .set_formation_enthalpy(ideal.formation_enthalpy() - props.enthalpy)
                    .unwrap();
            }
        }
        let mut guess = reference.guess();
        guess.temperature *= 1.005;
        guess.pressure *= 0.99;
        let spec = if entropy {
            QualitySpecification::Entropy(0.)
        } else {
            QualitySpecification::Enthalpy(0.)
        };
        let state = model
            .quality_state_with_specification(&reference.feed, 0.25, spec, &guess)
            .unwrap();
        check_target(&model, &state, spec);
        assert!((state.temperature - reference.temperature).abs() < 2e-6);
        assert!((state.pressure - reference.pressure).abs() < 0.1);
        let refined = model.refine_quality_specification(&state, spec).unwrap();
        assert_eq!(refined.feed, reference.feed);
    }
}

#[test]
fn caloric_quality_corrections_close_h_s_and_recover_independent_tp_states() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.4, 0.6, 0.];
        for beta in [0., 0.1, 0.25, 0.5, 0.75, 0.9, 1.] {
            let reference = model.quality_temperature(1e6, &z, beta).unwrap();
            let properties = model.quality_properties(&reference).unwrap();
            assert_eq!(
                properties.enthalpy,
                properties.internal_energy + reference.pressure * properties.volume
            );
            let mut guess = reference.guess();
            guess.temperature *= 1.005;
            guess.pressure *= 0.99;
            for k in &mut guess.ratios {
                *k = k.powf(0.98);
            }
            for spec in [
                QualitySpecification::Enthalpy(properties.enthalpy),
                QualitySpecification::Entropy(properties.entropy),
            ] {
                let actual = model
                    .quality_state_with_specification(&z, beta, spec, &guess)
                    .unwrap();
                assert!(
                    (actual.temperature - reference.temperature).abs() < 2e-6,
                    "{kind:?} {beta} {spec:?}"
                );
                assert!((actual.pressure - reference.pressure).abs() < 0.1);
                assert_eq!(actual.feed, z);
                assert_eq!(actual.liquid[2], 0.);
                assert_eq!(actual.vapor[2], 0.);
                assert!(actual.fugacity_residual < 1e-9);
                assert!(actual.material_balance_residual < 1e-9);
                let p = model.quality_properties(&actual).unwrap();
                let residual = match spec {
                    QualitySpecification::Enthalpy(h) => (p.enthalpy - h).abs() / h.abs().max(1.),
                    QualitySpecification::Entropy(s) => (p.entropy - s).abs() / s.abs().max(1.),
                    _ => unreachable!(),
                };
                assert!(residual < 1e-9);
                assert!((actual.dpdt - reference.dpdt).abs() < 1e-3 + 1e-7 * reference.dpdt.abs());
            }
        }
    }
}

#[test]
fn pure_caloric_quality_uses_coexistence_and_baseline_retains_pure_rejection() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let z = [0., 1.];
    for beta in [0., 0.3, 1.] {
        let reference = model.quality_temperature(1e6, &z, beta).unwrap();
        let p = model.quality_properties(&reference).unwrap();
        let guess = QualityGuess {
            temperature: reference.temperature * 1.01,
            pressure: reference.pressure * 0.95,
            ratios: vec![1., 1.],
        };
        for spec in [
            QualitySpecification::Enthalpy(p.enthalpy),
            QualitySpecification::Entropy(p.entropy),
        ] {
            let a = model
                .quality_state_with_specification(&z, beta, spec, &guess)
                .unwrap();
            assert!((a.temperature - reference.temperature).abs() < 2e-6);
            assert!((a.pressure - reference.pressure).abs() < 0.1);
            assert!(
                model
                    .baseline_quality_state_with_specification(&z, beta, spec, &guess)
                    .is_err()
            );
        }
    }
    let guess = QualityGuess {
        temperature: 200.,
        pressure: 1e6,
        ratios: vec![2., 0.5],
    };
    assert!(
        model
            .quality_state_with_specification(
                &[0.5, 0.5],
                0.4,
                QualitySpecification::Entropy(f64::NAN),
                &guess
            )
            .is_err()
    );
}
