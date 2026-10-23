use thermopack_rs::envelope::{EnvelopeEnd, EnvelopeOptions};
use thermopack_rs::{Cubic, CubicKind};

#[test]
fn fixed_quality_curves_preserve_feed_close_vle_and_complement_phase_labels() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0.5, 0.5, 0.];
    for beta in [0.25, 0.5, 0.75] {
        let curve = model
            .quality_envelope(&z, beta, EnvelopeOptions::default())
            .unwrap();
        assert_eq!(
            curve.end,
            EnvelopeEnd::InitialPressure,
            "{:?}",
            curve.failure
        );
        assert!(curve.points.len() > 20);
        assert!(curve.failure.is_none());
        assert!((curve.points.last().unwrap().pressure - 1e5).abs() < 0.01);
        for p in &curve.points {
            assert_eq!(p.feed, z);
            assert!(p.liquid_volume <= p.vapor_volume);
            assert!(
                (p.vapor_fraction - beta).abs() < 1e-12
                    || (p.vapor_fraction - (1. - beta)).abs() < 1e-12
            );
            assert!(p.fugacity_residual < 1e-9);
            assert!(p.material_balance_residual < 1e-9);
            assert_eq!(p.liquid[2], 0.);
            assert_eq!(p.vapor[2], 0.);
            for (i, &zi) in z.iter().enumerate() {
                assert!(
                    ((1. - p.vapor_fraction) * p.liquid[i] + p.vapor_fraction * p.vapor[i] - zi)
                        .abs()
                        < 1e-9
                );
            }
        }
    }
}

#[test]
fn pressure_limits_point_caps_and_endpoint_dispatch_are_explicit() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let curve = model
        .quality_envelope(
            &[0.5, 0.5],
            0.4,
            EnvelopeOptions {
                maximum_pressure: 2e6,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(
        curve.end,
        EnvelopeEnd::MaximumPressure,
        "{:?}",
        curve.failure
    );
    assert!((curve.points.last().unwrap().pressure - 2e6).abs() < 0.01);
    let capped = model
        .quality_envelope(
            &[0.5, 0.5],
            0.4,
            EnvelopeOptions {
                maximum_points: 3,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(capped.end, EnvelopeEnd::PointLimit);
    assert_eq!(capped.points.len(), 3);
    for beta in [0., 1.] {
        let curve = model
            .quality_envelope(
                &[0.5, 0.5],
                beta,
                EnvelopeOptions {
                    maximum_points: 3,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(curve.points[0].vapor_fraction, beta);
    }
    let curve = model
        .quality_envelope(&[1., 0.], 0.4, EnvelopeOptions::default())
        .unwrap();
    assert_eq!(curve.end, EnvelopeEnd::CriticalPoint);
    assert!(curve.points.iter().all(|p| p.vapor_fraction == 0.4));
    let rounded = model
        .baseline_quality_envelope(
            &[1., 1e-18],
            0.4,
            EnvelopeOptions {
                maximum_points: 3,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(rounded.end, EnvelopeEnd::SourceCriticalParameters);
    assert!(rounded.points.iter().all(|p| p.vapor_fraction == 0.));
    assert!(
        model
            .quality_envelope(&[0.5, 0.5], f64::NAN, EnvelopeOptions::default())
            .is_err()
    );
}

#[test]
fn physical_temperature_limit_closes_before_a_solver_bound() {
    let mut model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    model.set_temperature_bounds(80., 160.).unwrap();
    let curve = model
        .quality_envelope(&[0.5, 0.5], 0.4, EnvelopeOptions::default())
        .unwrap();
    assert_eq!(
        curve.end,
        EnvelopeEnd::MaximumTemperature,
        "{:?}",
        curve.failure
    );
    assert!((curve.points.last().unwrap().temperature - 160.).abs() < 1e-7);
    assert!(
        curve
            .points
            .iter()
            .all(|p| p.temperature <= 160. * (1. + 1e-10))
    );
}

#[test]
fn quantum_curve_endings_are_qualified_against_complete_source_returns() {
    let model = Cubic::quantum(&["H2", "D2"]).unwrap();
    for beta in [0.25, 0.5, 0.75] {
        let curve = model
            .baseline_quality_envelope(&[0.5, 0.5], beta, EnvelopeOptions::default())
            .unwrap();
        assert_eq!(
            curve.end,
            EnvelopeEnd::InitialPressure,
            "beta={beta}: {:?}, last={:?}",
            curve.failure,
            curve.points.last()
        );
        let physical = model
            .quality_envelope(&[0.5, 0.5], beta, EnvelopeOptions::default())
            .unwrap();
        assert_eq!(
            physical.end,
            EnvelopeEnd::InitialPressure,
            "beta={beta}: {:?}",
            physical.failure
        );
    }
}

#[test]
fn complete_source_grids_and_independent_warm_recovery_keep_raw_gaps_explicit() {
    use std::collections::BTreeMap;
    use thermopack_rs::envelope::EnvelopeConstraint;
    use thermopack_rs::quality::QualityState;
    use thermopack_rs::{CubicOptions, MixingRule};
    let mut groups = BTreeMap::<&str, Vec<Vec<f64>>>::new();
    for row in include_str!("fixtures/quality_envelope.csv")
        .lines()
        .skip(1)
    {
        let c: Vec<_> = row.split(';').collect();
        groups
            .entry(c[0])
            .or_default()
            .push(c[2].split(',').map(|v| v.parse().unwrap()).collect());
    }
    let compare = |a: &[f64], b: &QualityState| {
        [
            b.temperature,
            b.pressure,
            b.vapor_fraction,
            b.liquid_volume,
            b.vapor_volume,
        ]
        .into_iter()
        .chain(b.liquid.iter().copied())
        .chain(b.vapor.iter().copied())
        .enumerate()
        .all(|(j, b)| {
            let floor = match j {
                0 => 2e-5,
                1 => 2.,
                2 => 0.,
                3 | 4 => 1e-10,
                _ => 2e-8,
            };
            (a[j] - b).abs() <= floor + 2e-6 * a[j].abs().max(b.abs())
        })
    };
    let (mut cases, mut states, mut grids, mut gaps) = (0, 0, 0, 0);
    for row in include_str!("fixtures/quality_envelope_curves.csv")
        .lines()
        .skip(1)
    {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let names: Vec<_> = f[2].split(',').collect();
        let z: Vec<f64> = f[4].split(',').map(|v| v.parse().unwrap()).collect();
        let n = z.len();
        let o: Vec<f64> = f[5].split(':').map(|v| v.parse().unwrap()).collect();
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
                    volume_translation: f[3] == "translated",
                    mixing_rule: if f[3] == "unifac" {
                        MixingRule::Unifac
                    } else {
                        MixingRule::VdW
                    },
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let curve = model
            .baseline_quality_envelope(
                &z,
                o[6],
                EnvelopeOptions {
                    start: if o[5] > 0. {
                        EnvelopeConstraint::Temperature(o[5])
                    } else {
                        EnvelopeConstraint::Pressure(o[0])
                    },
                    maximum_pressure: o[1],
                    minimum_temperature: Some(o[2]),
                    maximum_step: o[3],
                    maximum_points: o[4] as usize,
                    ..Default::default()
                },
            )
            .unwrap();
        let rows = &groups[c[0]];
        assert_eq!(rows.len(), c[1].parse::<usize>().unwrap());
        if c[3] == "1" {
            assert_eq!(rows.len(), curve.points.len(), "{}", c[0]);
            for (i, (a, b)) in rows.iter().zip(&curve.points).enumerate() {
                assert!(compare(a, b), "{} node {i}", c[0]);
            }
            grids += 1;
            states += rows.len();
            cases += 1;
            continue;
        }
        for (index, a) in rows.iter().enumerate() {
            let k = &a[12 + 3 * n..12 + 4 * n];
            let coordinate = a[11 + 3 * n] as usize;
            let mut best = None;
            for p in &curve.points {
                for flip in [false, true] {
                    let beta = if flip {
                        1. - p.vapor_fraction
                    } else {
                        p.vapor_fraction
                    };
                    if (beta - a[2]).abs() > 1e-12 {
                        continue;
                    }
                    let mut guess = p.guess();
                    if flip {
                        for k in &mut guess.ratios {
                            *k = 1. / *k;
                        }
                    }
                    let distance = ((p.temperature / a[0]).ln() / 0.1).powi(2)
                        + ((p.pressure / a[1]).ln() / 0.2).powi(2)
                        + guess
                            .ratios
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| z[*i] > 0.)
                            .map(|(i, v)| (v.ln() - k[i]).powi(2))
                            .sum::<f64>();
                    if best.as_ref().is_none_or(|(d, _)| distance < *d) {
                        best = Some((distance, guess));
                    }
                }
            }
            let constraint = if coordinate < n {
                EnvelopeConstraint::LogK {
                    component: coordinate,
                    value: k[coordinate],
                }
            } else if coordinate == n {
                EnvelopeConstraint::Temperature(a[0])
            } else {
                EnvelopeConstraint::Pressure(a[1])
            };
            let p = model
                .baseline_quality_state_with_guess(&z, a[2], constraint, &best.unwrap().1)
                .unwrap();
            let p = if coordinate < n {
                model
                    .quality_state_with_guess(&z, a[2], constraint, &p.guess())
                    .unwrap_or(p)
            } else {
                p
            };
            if !compare(a, &p) {
                let known = (f[1..4] == ["PR", "H2,NE", "quantum"]
                    && o[6] == 0.25
                    && (20..28).contains(&index))
                    || (f[1..4] == ["SRK", "C1,C3", "unifac"] && o[6] == 0.5 && index == 26);
                assert!(known, "{} node {index}: {p:?}", c[0]);
                gaps += 1;
            } else {
                states += 1;
            }
            assert!(p.fugacity_residual < 1e-8);
            assert!(p.material_balance_residual < 1e-8);
        }
        cases += 1;
    }
    assert_eq!((cases, grids), (156, 148));
    assert_eq!(states + gaps, 8019);
    assert!(gaps <= 9);
}
