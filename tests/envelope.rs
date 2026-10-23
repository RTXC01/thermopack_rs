use thermopack_rs::envelope::{EnvelopeConstraint, EnvelopeEnd, EnvelopeOptions};
use thermopack_rs::{Cubic, CubicKind, SaturationKind};

#[test]
fn quantum_envelope_initialization_keeps_the_gas_liquid_branch() {
    let model = Cubic::quantum(&["H2", "NE"]).unwrap();
    let curve = model
        .phase_envelope(
            &[0.5, 0.5],
            EnvelopeOptions {
                start: EnvelopeConstraint::Pressure(1e6),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
    let first = &curve.points[0];
    assert!((first.temperature - 33.05026190854).abs() < 1e-5);
    assert!(first.vapor_volume > 5. * first.liquid_volume);
    assert!(curve.points.iter().all(|p| p.pressure < 2e6));
}

fn close_source(
    reference: &[f64],
    native: &thermopack_rs::MixtureSaturation,
    label: &str,
    index: usize,
) {
    let values: Vec<_> = [
        native.temperature,
        native.pressure,
        f64::from(native.kind == SaturationKind::Dew),
        native.liquid_volume,
        native.vapor_volume,
    ]
    .into_iter()
    .chain(native.liquid.iter().copied())
    .chain(native.vapor.iter().copied())
    .collect();
    for (j, (&a, b)) in reference.iter().zip(values).enumerate() {
        let floor = match j {
            0 => 2e-5,
            1 => 2.,
            2 => 0.,
            3 | 4 => 1e-10,
            _ => 2e-8,
        };
        assert!(
            (a - b).abs() <= floor + 2e-6 * a.abs().max(b.abs()),
            "{label} {index} {j}: {b} vs {a}"
        );
    }
}

#[test]
fn independently_traced_curves_recover_all_source_envelope_states() {
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<&str, Vec<(usize, Vec<f64>)>> = BTreeMap::new();
    for line in include_str!("fixtures/envelope.csv").lines().skip(1) {
        let f: Vec<_> = line.split(';').collect();
        groups.entry(f[0]).or_default().push((
            f[1].parse().unwrap(),
            f[2].split(',').map(|s| s.parse().unwrap()).collect(),
        ));
    }
    let mut cases = 0;
    let mut states = 0;
    let mut complete_grids = 0;
    for line in include_str!("fixtures/envelope_curves.csv").lines().skip(1) {
        let f: Vec<_> = line.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let s: Vec<f64> = c[5].split(':').map(|s| s.parse().unwrap()).collect();
        let names: Vec<_> = c[2].split(',').collect();
        let z: Vec<f64> = c[4].split(',').map(|s| s.parse().unwrap()).collect();
        let n = z.len();
        let kind = match c[1] {
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => CubicKind::PengRobinson,
        };
        let model = if c[3] == "quantum" {
            Cubic::quantum(&names).unwrap()
        } else {
            Cubic::with_options(
                kind,
                &names,
                thermopack_rs::CubicOptions {
                    volume_translation: c[3] == "translated",
                    mixing_rule: if c[3] == "unifac" {
                        thermopack_rs::MixingRule::Unifac
                    } else {
                        thermopack_rs::MixingRule::VdW
                    },
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let curve = model
            .baseline_phase_envelope(
                &z,
                EnvelopeOptions {
                    start: if s[5] > 0. {
                        EnvelopeConstraint::Temperature(s[5])
                    } else {
                        EnvelopeConstraint::Pressure(s[0])
                    },
                    initial_kind: if s[6] > 0.5 {
                        SaturationKind::Dew
                    } else {
                        SaturationKind::Bubble
                    },
                    maximum_pressure: s[1],
                    minimum_temperature: Some(s[2]),
                    maximum_step: s[3],
                    maximum_points: s[4] as usize,
                },
            )
            .unwrap();
        let expected_end = match f[2] {
            "InitialPressure" => EnvelopeEnd::InitialPressure,
            "MaximumPressure" => EnvelopeEnd::MaximumPressure,
            "MinimumTemperature" => EnvelopeEnd::MinimumTemperature,
            "PointLimit" => EnvelopeEnd::PointLimit,
            _ => EnvelopeEnd::NoConvergence,
        };
        if expected_end != EnvelopeEnd::NoConvergence {
            assert_eq!(curve.end, expected_end, "{}", f[0]);
        }
        let rows = groups.get(f[0]).unwrap();
        assert_eq!(rows.len(), f[1].parse::<usize>().unwrap());
        if f[3] == "1" {
            assert_eq!(curve.points.len(), rows.len(), "{}", f[0]);
            for ((i, a), b) in rows.iter().zip(&curve.points) {
                close_source(&a[..5 + 2 * n], b, f[0], *i);
            }
            complete_grids += 1;
        }
        for (index, a) in rows {
            let target_kind = if a[2] > 0.5 {
                SaturationKind::Dew
            } else {
                SaturationKind::Bubble
            };
            let target_sign = if target_kind == SaturationKind::Bubble {
                1.
            } else {
                -1.
            };
            let log_k = &a[12 + 3 * n..12 + 4 * n];
            let side = |p: &thermopack_rs::MixtureSaturation| {
                let sign = if p.kind == SaturationKind::Bubble {
                    1.
                } else {
                    -1.
                };
                (0..n)
                    .filter(|&i| z[i] > 0.)
                    .map(|i| sign * (p.vapor[i] / p.liquid[i]).ln() * target_sign * log_k[i])
                    .sum::<f64>()
            };
            let distance = |p: &thermopack_rs::MixtureSaturation| {
                let sign = if p.kind == SaturationKind::Bubble {
                    1.
                } else {
                    -1.
                };
                ((p.temperature / a[0]).ln() / 0.1).powi(2)
                    + ((p.pressure / a[1]).ln() / 0.2).powi(2)
                    + (0..n)
                        .filter(|&i| z[i] > 0.)
                        .map(|i| {
                            (sign * (p.vapor[i] / p.liquid[i]).ln() - target_sign * log_k[i])
                                .powi(2)
                        })
                        .sum::<f64>()
            };
            let mut seed = curve
                .points
                .iter()
                .filter(|p| side(p) > 0.)
                .min_by(|a, b| distance(a).total_cmp(&distance(b)))
                .or_else(|| {
                    curve
                        .points
                        .iter()
                        .min_by(|a, b| distance(a).total_cmp(&distance(b)))
                })
                .unwrap()
                .clone();
            if seed.kind != target_kind {
                std::mem::swap(&mut seed.liquid, &mut seed.vapor);
                std::mem::swap(&mut seed.liquid_volume, &mut seed.vapor_volume);
                seed.kind = target_kind;
            }
            let coordinate = a[11 + 3 * n] as usize;
            let constraint = if coordinate < n {
                EnvelopeConstraint::LogK {
                    component: coordinate,
                    value: log_k[coordinate],
                }
            } else if coordinate == n {
                EnvelopeConstraint::Temperature(a[0])
            } else {
                EnvelopeConstraint::Pressure(a[1])
            };
            let point = model
                .refine_envelope_state(&seed, constraint)
                .unwrap_or_else(|e| panic!("{} {index}: {e}", f[0]));
            close_source(&a[..5 + 2 * n], &point, f[0], *index);
            assert!(point.fugacity_residual < 1.01e-9);
            states += 1;
        }
        let physical = model
            .phase_envelope(
                &z,
                EnvelopeOptions {
                    start: if s[5] > 0. {
                        EnvelopeConstraint::Temperature(s[5])
                    } else {
                        EnvelopeConstraint::Pressure(s[0])
                    },
                    initial_kind: if s[6] > 0.5 {
                        SaturationKind::Dew
                    } else {
                        SaturationKind::Bubble
                    },
                    maximum_pressure: s[1],
                    minimum_temperature: Some(s[2]),
                    maximum_step: s[3],
                    maximum_points: s[4] as usize,
                },
            )
            .unwrap();
        if expected_end != EnvelopeEnd::NoConvergence {
            assert_ne!(
                physical.end,
                EnvelopeEnd::NoConvergence,
                "{} {:?}",
                f[0],
                physical.failure
            );
        }
        for point in &physical.points {
            assert!(point.liquid_volume <= point.vapor_volume);
            assert!(point.temperature >= s[2] - 1e-7 && point.pressure <= s[1] + 2.);
            let mut l = vec![0.; n];
            let mut v = vec![0.; n];
            let liquid = model.prepare(point.temperature, &point.liquid).unwrap();
            let vapor = model.prepare(point.temperature, &point.vapor).unwrap();
            assert!(
                (liquid.pressure(point.liquid_volume).unwrap().value - point.pressure).abs()
                    < 2. + 1e-7 * point.pressure
            );
            assert!(
                (vapor.pressure(point.vapor_volume).unwrap().value - point.pressure).abs()
                    < 2. + 1e-7 * point.pressure
            );
            liquid
                .ln_fugacity_tp_into(point.pressure, point.liquid_volume, &mut l)
                .unwrap();
            vapor
                .ln_fugacity_tp_into(point.pressure, point.vapor_volume, &mut v)
                .unwrap();
            for i in 0..n {
                if z[i] == 0. {
                    assert_eq!((point.liquid[i], point.vapor[i]), (0., 0.));
                } else {
                    assert!(
                        (point.liquid[i].ln() + l[i] - point.vapor[i].ln() - v[i]).abs() < 3e-8,
                        "{}",
                        f[0]
                    );
                }
            }
        }
        cases += 1;
    }
    assert_eq!((cases, states, complete_grids), (98, 6185, 71));
}

#[test]
fn invalid_envelope_inputs_and_seeds_are_rejected() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    for options in [
        EnvelopeOptions {
            maximum_points: 0,
            ..Default::default()
        },
        EnvelopeOptions {
            maximum_pressure: 0.,
            ..Default::default()
        },
        EnvelopeOptions {
            maximum_step: f64::NAN,
            ..Default::default()
        },
        EnvelopeOptions {
            minimum_temperature: Some(0.),
            ..Default::default()
        },
        EnvelopeOptions {
            start: EnvelopeConstraint::LogK {
                component: 0,
                value: 1.,
            },
            ..Default::default()
        },
    ] {
        assert!(model.phase_envelope(&[0.5, 0.5], options).is_err());
    }
    for z in [&[1.][..], &[0., 0.], &[1., -0.1]] {
        assert!(model.phase_envelope(z, EnvelopeOptions::default()).is_err());
    }
    let seed = model.bubble_pressure(180., &[0.5, 0.5]).unwrap();
    for constraint in [
        EnvelopeConstraint::Pressure(-1.),
        EnvelopeConstraint::Temperature(f64::NAN),
        EnvelopeConstraint::LogK {
            component: 2,
            value: 0.,
        },
        EnvelopeConstraint::LogK {
            component: 0,
            value: f64::INFINITY,
        },
    ] {
        assert!(model.refine_envelope_state(&seed, constraint).is_err());
    }
    let mut bad = seed;
    bad.vapor[0] = 0.;
    assert!(
        model
            .refine_envelope_state(&bad, EnvelopeConstraint::Pressure(1e5))
            .is_err()
    );
}

#[test]
fn trace_complete_mixture_curves() {
    for (kind, names, z) in [
        (CubicKind::PengRobinson, vec!["C1", "C2"], vec![0.5, 0.5]),
        (
            CubicKind::SoaveRedlichKwong,
            vec!["C1", "C2"],
            vec![0.5, 0.5],
        ),
        (
            CubicKind::PengRobinson,
            vec!["C1", "C2", "C3"],
            vec![0.2, 0.3, 0.5],
        ),
        (CubicKind::PengRobinson, vec!["CO2", "C1"], vec![0.8, 0.2]),
    ] {
        let model = Cubic::new(kind, &names).unwrap();
        let curve = model
            .phase_envelope(&z, EnvelopeOptions::default())
            .unwrap();
        eprintln!(
            "{kind:?} {names:?}: {:?} n={} failure={:?} bubble={} Tmax={} Pmax={}",
            curve.end,
            curve.points.len(),
            curve.failure,
            curve
                .points
                .iter()
                .filter(|p| p.kind == SaturationKind::Bubble)
                .count(),
            curve
                .points
                .iter()
                .map(|p| p.temperature)
                .fold(0., f64::max),
            curve.points.iter().map(|p| p.pressure).fold(0., f64::max)
        );
        assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
        assert!(
            curve
                .points
                .iter()
                .any(|p| p.kind == SaturationKind::Bubble)
        );
        assert!((curve.points.last().unwrap().pressure - 1e5).abs() < 1e-4);
        for p in &curve.points {
            assert!(p.liquid_volume <= p.vapor_volume);
            let n = z.len();
            let mut a = vec![0.; n];
            let mut b = vec![0.; n];
            let l = model.prepare(p.temperature, &p.liquid).unwrap();
            let v = model.prepare(p.temperature, &p.vapor).unwrap();
            l.ln_fugacity_tp_into(p.pressure, p.liquid_volume, &mut a)
                .unwrap();
            v.ln_fugacity_tp_into(p.pressure, p.vapor_volume, &mut b)
                .unwrap();
            for i in 0..n {
                assert!((p.liquid[i].ln() + a[i] - p.vapor[i].ln() - b[i]).abs() < 3e-8);
            }
            assert!(p.fugacity_residual < 1e-9);
        }
    }
}

#[test]
fn fixed_log_k_corrector_moves_a_warm_saturation_state() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let seed = model.bubble_pressure(180., &[0.5, 0.5]).unwrap();
    let target = (seed.vapor[0] / seed.liquid[0]).ln() + 0.01;
    let point = model
        .refine_envelope_state(
            &seed,
            EnvelopeConstraint::LogK {
                component: 0,
                value: target,
            },
        )
        .unwrap();
    assert!(((point.vapor[0] / point.liquid[0]).ln() - target).abs() < 2e-9);
    assert!((point.temperature - seed.temperature).abs() > 0.01);
}

#[test]
fn physical_pressure_temperature_and_point_limits_are_explicit() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    for (options, expected) in [
        (
            EnvelopeOptions {
                maximum_pressure: 2e6,
                ..Default::default()
            },
            EnvelopeEnd::MaximumPressure,
        ),
        (
            EnvelopeOptions {
                minimum_temperature: Some(150.),
                ..Default::default()
            },
            EnvelopeEnd::MinimumTemperature,
        ),
        (
            EnvelopeOptions {
                maximum_points: 3,
                ..Default::default()
            },
            EnvelopeEnd::PointLimit,
        ),
    ] {
        let curve = model.phase_envelope(&[0.5, 0.5], options).unwrap();
        assert_eq!(curve.end, expected);
        assert!(curve.failure.is_none());
        let last = curve.points.last().unwrap();
        if expected == EnvelopeEnd::MaximumPressure {
            assert!((last.pressure - 2e6).abs() < 1e-4);
            assert!(curve.points.len() < 100);
        }
        if expected == EnvelopeEnd::MinimumTemperature {
            assert!((last.temperature - 150.).abs() < 1e-8);
        }
        if expected == EnvelopeEnd::PointLimit {
            assert_eq!(curve.points.len(), 3);
        }
    }
    let cryogenic = Cubic::new(CubicKind::PengRobinson, &["N2", "C1"]).unwrap();
    let curve = cryogenic
        .phase_envelope(
            &[0.95, 0.05],
            EnvelopeOptions {
                minimum_temperature: Some(20.),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(curve.end, EnvelopeEnd::MinimumTemperature);
    assert!((curve.points.last().unwrap().temperature - 80.).abs() < 1e-8);
}
