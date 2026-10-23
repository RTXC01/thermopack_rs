use thermopack_rs::pure_saturation_curve::{
    PureSaturationEnd, PureSaturationOptions, PureSaturationStart,
};
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule};

#[test]
fn complete_source_grids_match_with_explicit_critical_root_conditioning() {
    let mut cases = std::collections::BTreeMap::<&str, Vec<Vec<f64>>>::new();
    for row in include_str!("fixtures/pure_curve.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        cases
            .entry(columns[0])
            .or_default()
            .push(columns[2..].iter().map(|s| s.parse().unwrap()).collect());
    }
    let mut count = 0;
    let mut critical_bounds = 0;
    assert_eq!(cases.len(), 201);
    for (case, rows) in cases {
        let f: Vec<_> = case.split('|').collect();
        let o: Vec<_> = f[4].split(':').collect();
        let kind = match f[1] {
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => CubicKind::PengRobinson,
        };
        let names: Vec<_> = f[2].split(',').collect();
        let model = if f[3] == "quantum" {
            Cubic::quantum(&names).unwrap()
        } else {
            Cubic::with_options(
                kind,
                &names,
                CubicOptions {
                    volume_translation: f[3] == "translated",
                    parameter_reference: f[5],
                    mixing_rule: match f[3] {
                        "unifac" => MixingRule::Unifac,
                        "hv1" => MixingRule::HuronVidal1,
                        "hv2" => MixingRule::HuronVidal2,
                        "nrtl" => MixingRule::Nrtl,
                        _ => MixingRule::VdW,
                    },
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let component = o[0].parse().unwrap();
        let options = PureSaturationOptions {
            start: if o[1] == "T" {
                PureSaturationStart::Temperature(o[2].parse().unwrap())
            } else {
                PureSaturationStart::Pressure(o[2].parse().unwrap())
            },
            maximum_pressure_step: o[3].parse().unwrap(),
            maximum_points: o[4].parse().unwrap(),
            logarithmic: o[5] == "1",
        };
        let curve = model
            .baseline_pure_saturation_curve(component, options)
            .unwrap();
        assert_eq!(
            curve.end,
            PureSaturationEnd::SourceCriticalParameters,
            "{case}: {:?}",
            curve.failure
        );
        assert_eq!(curve.points.len(), rows.len(), "{case}");
        for (index, (source, native)) in rows.iter().zip(&curve.points).enumerate() {
            let values = [
                native.temperature,
                native.pressure,
                native.liquid_volume,
                native.vapor_volume,
            ];
            for j in 0..4 {
                let floor = [2e-6, 0.1, 1e-12, 1e-12][j];
                let tolerance = floor + 2e-7 * source[j].abs().max(values[j].abs());
                let mut bound = tolerance;
                if index == rows.len() - 1 && j >= 2 {
                    let (a, b, c, shift, eta, radius) = (
                        source[7], source[8], source[9], source[10], source[11], source[12],
                    );
                    let center = -a / 3.;
                    let p = b - a * a / 3.;
                    let q = 2. * a * a * a / 27. - a * b / 3. + c;
                    let scale = center.abs().powi(3)
                        + a.abs() * center.powi(2)
                        + (b * center).abs()
                        + c.abs();
                    assert!((eta - 32. * f64::EPSILON * scale).abs() < 1e-28);
                    let mut low = 0.;
                    let mut high = 1.;
                    for _ in 0..100 {
                        let r = (low + high) * 0.5;
                        if r * r * r - p.abs() * r - q.abs() - eta > 0. {
                            high = r;
                        } else {
                            low = r;
                        }
                    }
                    let expected = 2. * high * thermopack_rs::R * source[0] / source[1];
                    assert!((radius - expected).abs() < 1e-12 * radius);
                    for volume in [source[j], values[j]] {
                        let z = (volume + shift) * source[1] / (thermopack_rs::R * source[0]);
                        let backward = (z * z * z + a * z * z + b * z + c).abs();
                        assert!(
                            backward <= eta,
                            "{case}: polynomial residual {backward} > {eta}"
                        );
                    }
                    critical_bounds += 1;
                    bound += radius;
                }
                assert!(
                    (source[j] - values[j]).abs() <= bound,
                    "{case} point {index} quantity {j}: {} vs {}, bound {bound}",
                    values[j],
                    source[j]
                );
            }
            count += 1;
        }
        let physical = model.pure_saturation_curve(component, options).unwrap();
        assert_eq!(
            physical.end,
            PureSaturationEnd::CriticalPoint,
            "{case}: {:?}",
            physical.failure
        );
        let mut z = vec![0.; names.len()];
        z[component] = 1.;
        for (i, point) in physical.points.iter().enumerate() {
            let state = model.prepare(point.temperature, &z).unwrap();
            let mut l = vec![0.; z.len()];
            let mut v = l.clone();
            state
                .ln_fugacity_tp_into(point.pressure, point.liquid_volume, &mut l)
                .unwrap();
            state
                .ln_fugacity_tp_into(point.pressure, point.vapor_volume, &mut v)
                .unwrap();
            assert!((l[component] - v[component]).abs() < 3e-11, "{case} {i}");
            if i + 1 < physical.points.len() {
                assert!(point.liquid_volume < point.vapor_volume);
            } else {
                assert_eq!(point.liquid_volume, point.vapor_volume);
            }
        }
    }
    assert_eq!(count, 19352);
    assert_eq!(critical_bounds, 402);
}

#[test]
fn physical_grids_close_coexistence_and_finish_at_pressure_inflection() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for name in ["C1", "CO2", "H2O"] {
            let model = Cubic::new(kind, &[name]).unwrap();
            let c = &model.components()[0];
            for start in [
                PureSaturationStart::Pressure(0.05 * c.critical_pressure),
                PureSaturationStart::Temperature(0.65 * c.critical_temperature),
            ] {
                for logarithmic in [false, true] {
                    let curve = model
                        .pure_saturation_curve(
                            0,
                            PureSaturationOptions {
                                start,
                                logarithmic,
                                maximum_points: 17,
                                ..Default::default()
                            },
                        )
                        .unwrap();
                    assert_eq!(
                        curve.end,
                        PureSaturationEnd::CriticalPoint,
                        "{kind:?} {name} {start:?} {logarithmic}: {:?}",
                        curve.failure
                    );
                    assert!(curve.failure.is_none());
                    assert_eq!(curve.points.len(), 17);
                    for point in &curve.points[..curve.points.len() - 1] {
                        assert!(point.liquid_volume < point.vapor_volume);
                        let state = model.prepare(point.temperature, &[1.]).unwrap();
                        let mut l = [0.];
                        let mut v = [0.];
                        state
                            .ln_fugacity_tp_into(point.pressure, point.liquid_volume, &mut l)
                            .unwrap();
                        state
                            .ln_fugacity_tp_into(point.pressure, point.vapor_volume, &mut v)
                            .unwrap();
                        assert!((l[0] - v[0]).abs() < 3e-11);
                        for volume in [point.liquid_volume, point.vapor_volume] {
                            assert!(
                                (state.pressure(volume).unwrap().value - point.pressure).abs()
                                    < 1e-3 + 1e-7 * point.pressure
                            );
                        }
                        assert!(point.dpdt > 0. && point.dpdt.is_finite());
                    }
                    let last = curve.points.last().unwrap();
                    assert_eq!(last.liquid_volume, last.vapor_volume);
                    let state = model.prepare(last.temperature, &[1.]).unwrap();
                    let p = state.pressure(last.liquid_volume).unwrap();
                    assert!((p.value - last.pressure).abs() < 1e-5);
                    assert!((p.dv * last.liquid_volume / last.pressure).abs() < 1e-8);
                    assert!(
                        (state
                            .pressure_second_derivative(last.liquid_volume)
                            .unwrap()
                            * last.liquid_volume.powi(2)
                            / last.pressure)
                            .abs()
                            < 1e-8
                    );
                    assert_eq!(last.dpdt, p.dt);
                    let coordinates: Vec<_> = curve
                        .points
                        .iter()
                        .map(|p| match start {
                            PureSaturationStart::Pressure(_) => p.pressure,
                            PureSaturationStart::Temperature(_) => p.temperature,
                        })
                        .map(|x| if logarithmic { x.ln() } else { x })
                        .collect();
                    let step = (coordinates.last().unwrap() - coordinates[0]) / 16.;
                    for (i, x) in coordinates.iter().enumerate() {
                        assert!(
                            (x - coordinates[0] - step * i as f64).abs() < 1e-9 * x.abs().max(1.)
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn selected_component_retains_model_and_absent_species() {
    for mixing_rule in [
        MixingRule::VdW,
        MixingRule::Unifac,
        MixingRule::HuronVidal1,
        MixingRule::HuronVidal2,
        MixingRule::Nrtl,
    ] {
        let model = Cubic::with_options(
            CubicKind::PengRobinson,
            &["C1", "C2", "C3"],
            CubicOptions {
                mixing_rule,
                ..Default::default()
            },
        )
        .unwrap();
        for (component, name) in ["C1", "C2", "C3"].iter().enumerate() {
            let pure = Cubic::new(CubicKind::PengRobinson, &[name]).unwrap();
            let options = PureSaturationOptions {
                maximum_points: 17,
                ..Default::default()
            };
            let curve = model.pure_saturation_curve(component, options).unwrap();
            let expected = pure.pure_saturation_curve(0, options).unwrap();
            assert_eq!(
                curve.end,
                PureSaturationEnd::CriticalPoint,
                "{mixing_rule:?} {name}"
            );
            assert_eq!(curve.component, component);
            assert_eq!(curve.points.len(), expected.points.len());
            for (a, b) in curve.points.iter().zip(expected.points) {
                for (a, b) in [a.temperature, a.pressure, a.liquid_volume, a.vapor_volume]
                    .iter()
                    .zip([b.temperature, b.pressure, b.liquid_volume, b.vapor_volume])
                {
                    assert!(
                        (a - b).abs() < 1e-12 + 1e-8 * a.abs().max(b.abs()),
                        "{mixing_rule:?} {name}: {a} vs {b}"
                    );
                }
            }
        }
    }
}

#[test]
fn point_limit_resamples_the_entire_grid_and_invalid_input_is_rejected() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let expected = model
        .baseline_pure_saturation_curve(
            0,
            PureSaturationOptions {
                logarithmic: true,
                ..Default::default()
            },
        )
        .unwrap();
    for ignored_step in [0., -1., f64::NAN] {
        let actual = model
            .baseline_pure_saturation_curve(
                0,
                PureSaturationOptions {
                    logarithmic: true,
                    maximum_pressure_step: ignored_step,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(actual.points.len(), expected.points.len());
        for (a, b) in actual.points.iter().zip(&expected.points) {
            assert_eq!(
                [
                    a.temperature,
                    a.pressure,
                    a.liquid_volume,
                    a.vapor_volume,
                    a.dpdt
                ],
                [
                    b.temperature,
                    b.pressure,
                    b.liquid_volume,
                    b.vapor_volume,
                    b.dpdt
                ]
            );
        }
    }
    for maximum_points in [2, 3, 17] {
        let options = PureSaturationOptions {
            maximum_points,
            ..Default::default()
        };
        let curve = model.pure_saturation_curve(0, options).unwrap();
        assert_eq!(curve.points.len(), maximum_points);
        assert_eq!(curve.end, PureSaturationEnd::CriticalPoint);
        assert!((curve.points.first().unwrap().pressure - 1e5).abs() < 1e-4);
        assert!(
            (curve.points.last().unwrap().pressure - model.components()[0].critical_pressure).abs()
                < 1e-4
        );
    }
    for options in [
        PureSaturationOptions {
            maximum_points: 0,
            ..Default::default()
        },
        PureSaturationOptions {
            maximum_points: 1,
            ..Default::default()
        },
        PureSaturationOptions {
            maximum_pressure_step: 0.,
            ..Default::default()
        },
        PureSaturationOptions {
            maximum_pressure_step: f64::NAN,
            ..Default::default()
        },
        PureSaturationOptions {
            start: PureSaturationStart::Pressure(-1.),
            ..Default::default()
        },
        PureSaturationOptions {
            start: PureSaturationStart::Temperature(f64::NAN),
            ..Default::default()
        },
        PureSaturationOptions {
            start: PureSaturationStart::Pressure(1e8),
            ..Default::default()
        },
        PureSaturationOptions {
            start: PureSaturationStart::Temperature(200.),
            ..Default::default()
        },
    ] {
        assert!(model.pure_saturation_curve(0, options).is_err());
    }
    assert!(model.pure_saturation_curve(1, Default::default()).is_err());
}

#[test]
fn quantum_and_translated_curves_have_physical_critical_volumes() {
    for model in [
        Cubic::quantum(&["H2"]).unwrap(),
        Cubic::quantum(&["HE"]).unwrap(),
        Cubic::with_options(
            CubicKind::PengRobinson,
            &["CO2"],
            CubicOptions {
                volume_translation: true,
                ..Default::default()
            },
        )
        .unwrap(),
    ] {
        let c = &model.components()[0];
        let curve = model
            .pure_saturation_curve(
                0,
                PureSaturationOptions {
                    start: PureSaturationStart::Pressure(0.05 * c.critical_pressure),
                    maximum_points: 17,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(
            curve.end,
            PureSaturationEnd::CriticalPoint,
            "{:?}",
            curve.failure
        );
        let last = curve.points.last().unwrap();
        assert!(last.liquid_volume > 0.);
        let state = model.prepare(last.temperature, &[1.]).unwrap();
        assert!(
            (state.pressure(last.liquid_volume).unwrap().dv * last.liquid_volume / last.pressure)
                .abs()
                < 1e-8
        );
    }
}
