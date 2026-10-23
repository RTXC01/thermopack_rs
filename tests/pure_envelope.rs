use thermopack_rs::envelope::{EnvelopeConstraint, EnvelopeEnd, EnvelopeOptions};
use thermopack_rs::pure_saturation_curve::{
    PureSaturationEnd, PureSaturationOptions, PureSaturationStart,
};
use thermopack_rs::{Cubic, CubicKind, SaturationKind};

#[test]
fn source_pure_dispatch_and_critical_initializers_match_all_finite_grids() {
    use thermopack_rs::{CubicOptions, MixingRule};
    let mut cases = std::collections::BTreeMap::<&str, Vec<Vec<f64>>>::new();
    for row in include_str!("fixtures/pure_envelope.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        cases
            .entry(c[0])
            .or_default()
            .push(c[2..].iter().map(|v| v.parse().unwrap()).collect());
    }
    assert_eq!(cases.len(), 169);
    let mut points = 0;
    let mut critical_nodes = 0;
    for (case, rows) in cases {
        let f: Vec<_> = case.split('|').collect();
        let o: Vec<_> = f[4].split(':').collect();
        let names: Vec<_> = f[2].split(',').collect();
        let feed: Vec<f64> = f[6].split(',').map(|v| v.parse().unwrap()).collect();
        let model = if f[3] == "quantum" {
            Cubic::quantum(&names).unwrap()
        } else {
            let kind = match f[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            Cubic::with_options(
                kind,
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
                    parameter_reference: f[5],
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let start = if o[1] == "T" {
            EnvelopeConstraint::Temperature(o[2].parse().unwrap())
        } else {
            EnvelopeConstraint::Pressure(o[2].parse().unwrap())
        };
        let actual: Vec<_> = if f[0] == "PURE_ENV" {
            let e: Vec<f64> = f[7].split(':').map(|v| v.parse().unwrap()).collect();
            let curve = model
                .baseline_phase_envelope(
                    &feed,
                    EnvelopeOptions {
                        start,
                        maximum_points: o[4].parse().unwrap(),
                        maximum_pressure: e[0],
                        minimum_temperature: Some(e[1]),
                        maximum_step: e[2],
                        initial_kind: if e[3] == 0. {
                            SaturationKind::Bubble
                        } else {
                            SaturationKind::Dew
                        },
                    },
                )
                .unwrap();
            assert_eq!(
                curve.end,
                EnvelopeEnd::SourceCriticalParameters,
                "{case}: {:?}",
                curve.failure
            );
            assert!(curve.failure.is_none());
            curve
                .points
                .into_iter()
                .map(|p| {
                    assert_eq!(p.kind, SaturationKind::Bubble);
                    assert_eq!(p.liquid, feed);
                    assert_eq!(p.vapor, feed);
                    (
                        p.temperature,
                        p.pressure,
                        p.liquid_volume,
                        p.vapor_volume,
                        Some(p.fugacity_residual),
                    )
                })
                .collect()
        } else {
            let curve = model
                .baseline_pure_saturation_curve(
                    o[0].parse().unwrap(),
                    PureSaturationOptions {
                        start: match start {
                            EnvelopeConstraint::Temperature(t) => {
                                PureSaturationStart::Temperature(t)
                            }
                            EnvelopeConstraint::Pressure(p) => PureSaturationStart::Pressure(p),
                            _ => unreachable!(),
                        },
                        maximum_points: o[4].parse().unwrap(),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(
                curve.end,
                PureSaturationEnd::SourceCriticalParameters,
                "{case}: {:?}",
                curve.failure
            );
            curve
                .points
                .into_iter()
                .map(|p| {
                    (
                        p.temperature,
                        p.pressure,
                        p.liquid_volume,
                        p.vapor_volume,
                        None,
                    )
                })
                .collect()
        };
        assert_eq!(actual.len(), rows.len(), "{case}");
        for (index, (s, a)) in rows.iter().zip(actual).enumerate() {
            let values = [a.0, a.1, a.2, a.3];
            let critical = s[11] > 0.;
            for j in 0..4 {
                let mut tolerance =
                    [2e-6, 0.1, 1e-12, 1e-12][j] + 2e-7 * s[j].abs().max(values[j].abs());
                if critical && j >= 2 {
                    let (pa, pb, pc, shift, eta, radius) = (s[7], s[8], s[9], s[10], s[11], s[12]);
                    let center = -pa / 3.;
                    let scale = center.abs().powi(3)
                        + pa.abs() * center.powi(2)
                        + (pb * center).abs()
                        + pc.abs();
                    assert!((eta - 32. * f64::EPSILON * scale).abs() < 1e-28);
                    let p = pb - pa * pa / 3.;
                    let q = 2. * pa * pa * pa / 27. - pa * pb / 3. + pc;
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
                    assert!(
                        (radius - 2. * high * thermopack_rs::R * s[0] / s[1]).abs()
                            < 1e-12 * radius
                    );
                    for v in [s[j], values[j]] {
                        let z = (v + shift) * s[1] / (thermopack_rs::R * s[0]);
                        let error = (z * z * z + pa * z * z + pb * z + pc).abs();
                        assert!(error <= eta, "{case} {index} {j}: {error} > {eta}");
                    }
                    tolerance += radius;
                }
                assert!(
                    (s[j] - values[j]).abs() <= tolerance,
                    "{case} {index} {j}: {} vs {}, tolerance {tolerance}",
                    values[j],
                    s[j]
                );
            }
            if let Some(residual) = a.4 {
                let state = model.prepare(a.0, &feed).unwrap();
                let mut l = vec![0.; feed.len()];
                let mut v = l.clone();
                state.ln_fugacity_tp_into(a.1, a.2, &mut l).unwrap();
                state.ln_fugacity_tp_into(a.1, a.3, &mut v).unwrap();
                let expected = feed
                    .iter()
                    .enumerate()
                    .filter(|(_, x)| **x > 0.)
                    .map(|(i, _)| (l[i] - v[i]).abs())
                    .fold(0., f64::max);
                assert!((residual - expected).abs() < 1e-12);
                if !critical {
                    assert!(
                        (residual - s[13]).abs() < 2e-9 + 2e-7 * residual.abs().max(s[13].abs())
                    );
                }
            }
            points += 1;
            critical_nodes += usize::from(critical);
        }
    }
    assert_eq!(points, 13917);
    assert_eq!(critical_nodes, 175);
}

#[test]
fn pure_and_embedded_components_reach_critical_points_with_full_phase_arrays() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        for index in 0..3 {
            let mut z = vec![0.; 3];
            z[index] = 1.;
            let curve = model
                .phase_envelope(
                    &z,
                    EnvelopeOptions {
                        maximum_pressure: 1e8,
                        maximum_points: 17,
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(
                curve.end,
                EnvelopeEnd::CriticalPoint,
                "{kind:?} {index}: {:?}",
                curve.failure
            );
            assert_eq!(curve.points.len(), 17);
            for point in &curve.points {
                assert_eq!(point.kind, SaturationKind::Bubble);
                assert_eq!(point.liquid, z);
                assert_eq!(point.vapor, z);
                assert!(point.fugacity_residual < 3e-11);
                for volume in [point.liquid_volume, point.vapor_volume] {
                    let state = model.prepare(point.temperature, &z).unwrap();
                    assert!(
                        (state.pressure(volume).unwrap().value - point.pressure).abs()
                            < 1e-3 + 1e-7 * point.pressure
                    );
                }
            }
            let last = curve.points.last().unwrap();
            assert_eq!(last.liquid_volume, last.vapor_volume);
            let state = model.prepare(last.temperature, &z).unwrap();
            assert!(
                (state.pressure(last.liquid_volume).unwrap().dv * last.liquid_volume
                    / last.pressure)
                    .abs()
                    < 1e-8
            );
        }
    }
}

#[test]
fn physical_limits_choose_first_boundary_and_single_point_prefix() {
    let mut model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    for (pmax, tmax, end) in [
        (2e6, 999., EnvelopeEnd::MaximumPressure),
        (2e6, 150., EnvelopeEnd::MaximumTemperature),
        (5e5, 150., EnvelopeEnd::MaximumPressure),
    ] {
        model.set_temperature_bounds(80., tmax).unwrap();
        let curve = model
            .phase_envelope(
                &[1.],
                EnvelopeOptions {
                    maximum_pressure: pmax,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(curve.end, end, "{:?}", curve.failure);
        let last = curve.points.last().unwrap();
        if end == EnvelopeEnd::MaximumPressure {
            assert!((last.pressure - pmax).abs() < 1e-3);
        } else {
            assert!((last.temperature - tmax).abs() < 1e-8);
        }
        assert!(
            curve
                .points
                .iter()
                .all(|p| p.pressure <= pmax * (1. + 1e-10) && p.temperature <= tmax * (1. + 1e-10))
        );
    }
    let curve = model
        .phase_envelope(
            &[1.],
            EnvelopeOptions {
                maximum_points: 1,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(curve.end, EnvelopeEnd::PointLimit);
    assert_eq!(curve.points.len(), 1);
    assert!((curve.points[0].pressure - 1e5).abs() < 1e-4);
    assert!(curve.failure.is_none());
}

#[test]
fn baseline_controls_and_roundoff_trace_classification_follow_source() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let z = [1., 1e-18];
    let options = EnvelopeOptions {
        maximum_points: 17,
        ..Default::default()
    };
    let a = model.baseline_phase_envelope(&z, options).unwrap();
    let b = model
        .baseline_phase_envelope(
            &z,
            EnvelopeOptions {
                maximum_pressure: 0.,
                minimum_temperature: Some(-1.),
                maximum_step: -1.,
                initial_kind: SaturationKind::Bubble,
                ..options
            },
        )
        .unwrap();
    assert_eq!(a.end, EnvelopeEnd::SourceCriticalParameters);
    assert_eq!(b.end, a.end);
    assert_eq!(a.points.len(), 17);
    assert_eq!(b.points.len(), 17);
    assert!(a.points.iter().any(|p| p.fugacity_residual > 1e-3));
    for (a, b) in a.points.iter().zip(&b.points) {
        assert_eq!(a.liquid, z);
        assert_eq!(a.vapor, z);
        assert_eq!(a.kind, SaturationKind::Bubble);
        assert_eq!(
            [a.temperature, a.pressure, a.liquid_volume, a.vapor_volume],
            [b.temperature, b.pressure, b.liquid_volume, b.vapor_volume]
        );
    }
    if let Ok(physical) = model.phase_envelope(
        &z,
        EnvelopeOptions {
            maximum_points: 1,
            ..Default::default()
        },
    ) {
        assert_ne!(physical.end, EnvelopeEnd::SourceCriticalParameters);
        assert_ne!(physical.end, EnvelopeEnd::CriticalPoint);
    }
}

#[test]
fn one_point_baseline_and_physical_critical_starts_have_explicit_semantics() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let c = &model.components()[0];
    let baseline = model
        .baseline_phase_envelope(
            &[1.],
            EnvelopeOptions {
                maximum_points: 1,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(baseline.end, EnvelopeEnd::SourceCriticalParameters);
    assert_eq!(baseline.points.len(), 1);
    assert_eq!(baseline.points[0].temperature, c.critical_temperature);
    assert_eq!(baseline.points[0].pressure, c.critical_pressure);
    for start in [
        PureSaturationStart::Temperature(c.critical_temperature),
        PureSaturationStart::Pressure(c.critical_pressure),
    ] {
        let curve = model
            .pure_saturation_curve(
                0,
                PureSaturationOptions {
                    start,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(curve.end, PureSaturationEnd::CriticalPoint);
        assert_eq!(curve.points.len(), 1);
        let envelope_start = match start {
            PureSaturationStart::Temperature(t) => EnvelopeConstraint::Temperature(t),
            PureSaturationStart::Pressure(p) => EnvelopeConstraint::Pressure(p),
        };
        let envelope = model
            .phase_envelope(
                &[1.],
                EnvelopeOptions {
                    start: envelope_start,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(envelope.end, EnvelopeEnd::CriticalPoint);
        assert_eq!(envelope.points.len(), 1);
    }
}
