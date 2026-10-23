use thermopack_rs::envelope::{EnvelopeConstraint, EnvelopeExtremum, EnvelopeOptions};
use thermopack_rs::{Cubic, CubicKind};

#[test]
fn independently_source_backed_extrema_match_local_and_automatic_refinement() {
    use thermopack_rs::{CubicOptions, MixingRule, MixtureSaturation, SaturationKind};
    let (mut total, mut local, mut flipped) = (0, 0, 0);
    for row in include_str!("fixtures/envelope_extrema.csv")
        .lines()
        .skip(1)
    {
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
                    volume_translation: f[3] == "translated",
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
        let z: Vec<f64> = f[4].split(',').map(|v| v.parse().unwrap()).collect();
        let beta: f64 = f[6].parse().unwrap();
        let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
        let kind = if f[5] == "P" {
            EnvelopeExtremum::Pressure
        } else {
            EnvelopeExtremum::Temperature
        };
        let mut x: Vec<_> = z
            .iter()
            .zip(&g[2..])
            .map(|(z, k)| z / (1. - beta + beta * k))
            .collect();
        let mut y: Vec<_> = x.iter().zip(&g[2..]).map(|(x, k)| x * k).collect();
        let sx = x.iter().sum::<f64>();
        let sy = y.iter().sum::<f64>();
        x.iter_mut().for_each(|v| *v /= sx);
        y.iter_mut().for_each(|v| *v /= sy);
        let guess = MixtureSaturation {
            kind: if beta == 0. {
                SaturationKind::Bubble
            } else {
                SaturationKind::Dew
            },
            temperature: g[0],
            pressure: g[1],
            liquid: x,
            vapor: y,
            liquid_volume: 1.,
            vapor_volume: 1.,
            dpdt: 0.,
            iterations: 0,
            fugacity_residual: 0.,
        };
        let expected: Vec<f64> = c[1].split(',').map(|v| v.parse().unwrap()).collect();
        let compare = |p: &MixtureSaturation, allow_flip: bool| {
            let flip = f64::from(p.kind == SaturationKind::Dew) != expected[2];
            assert!(!flip || allow_flip);
            let (liquid, vapor, vl, vv) = if flip {
                (&p.vapor, &p.liquid, p.vapor_volume, p.liquid_volume)
            } else {
                (&p.liquid, &p.vapor, p.liquid_volume, p.vapor_volume)
            };
            let mut values = vec![p.temperature, p.pressure, expected[2], vl, vv];
            values.extend(liquid);
            values.extend(vapor);
            for (j, (&a, b)) in expected.iter().zip(values).enumerate() {
                let floor = match j {
                    0 => 2e-5,
                    1 => 2.,
                    2 => 0.,
                    3 | 4 => 1e-10,
                    _ => 2e-8,
                };
                assert!(
                    b.is_finite() && (a - b).abs() <= floor + 2e-6 * a.abs().max(b.abs()),
                    "{} {j}: {a} vs {b}",
                    c[0]
                );
            }
            assert!(p.fugacity_residual < 1e-9);
            for (xs, v) in [(&p.liquid, p.liquid_volume), (&p.vapor, p.vapor_volume)] {
                for (&zi, &xi) in z.iter().zip(xs) {
                    if zi == 0. {
                        assert_eq!(xi, 0.);
                    }
                }
                assert!(
                    (model
                        .prepare(p.temperature, xs)
                        .unwrap()
                        .pressure(v)
                        .unwrap()
                        .value
                        - p.pressure)
                        .abs()
                        < 2. + 1e-7 * p.pressure
                );
            }
            flip
        };
        match model.refine_envelope_extremum(&z, &guess, kind) {
            Ok(p) => {
                compare(&p.state, false);
                assert!(p.stationarity_residual.abs() < 1e-9);
                local += 1;
            }
            Err(_) => {
                assert_eq!(&f[1..7], &["PR", "C1,C2", "classic", "0.5,0.5", "P", "0.0"]);
                assert_eq!(f[8], "100000.0,15000000.0,80.0,180.0,1.0,0.25,1000");
            }
        }
        let o: Vec<f64> = f[8].split(',').map(|v| v.parse().unwrap()).collect();
        let curve = model
            .phase_envelope(
                &z,
                EnvelopeOptions {
                    start: if o[3] > 0. {
                        EnvelopeConstraint::Temperature(o[3])
                    } else {
                        EnvelopeConstraint::Pressure(o[0])
                    },
                    initial_kind: if o[4] == 0. {
                        SaturationKind::Bubble
                    } else {
                        SaturationKind::Dew
                    },
                    maximum_pressure: o[1],
                    minimum_temperature: Some(o[2]),
                    maximum_step: o[5],
                    maximum_points: o[6] as usize,
                },
            )
            .unwrap();
        let maxima = model.envelope_extrema(&z, &curve).unwrap();
        let p = if kind == EnvelopeExtremum::Pressure {
            maxima.pressure
        } else {
            maxima.temperature
        }
        .unwrap();
        flipped += usize::from(compare(&p.state, true));
        assert!(p.stationarity_residual.abs() < 1e-9);
        for side in [-1., 1.] {
            let constraint = if kind == EnvelopeExtremum::Pressure {
                EnvelopeConstraint::Temperature(p.state.temperature * (1. + side * 5e-5))
            } else {
                EnvelopeConstraint::Pressure(p.state.pressure * (1. + side * 5e-5))
            };
            let neighbor = model.refine_envelope_state(&p.state, constraint).unwrap();
            assert!(
                if kind == EnvelopeExtremum::Pressure {
                    neighbor.pressure < p.state.pressure
                } else {
                    neighbor.temperature < p.state.temperature
                },
                "{} neighbor {side}",
                c[0]
            );
        }
        total += 1;
    }
    assert_eq!((total, flipped), (147, 97));
    assert!((146..=147).contains(&local));
}

#[test]
fn refine_detected_pressure_and_temperature_peaks_and_check_both_sides() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.5, 0.5, 0.];
        let curve = model
            .phase_envelope(&z, EnvelopeOptions::default())
            .unwrap();
        let extrema = model.envelope_extrema(&z, &curve).unwrap();
        assert!(
            extrema.failures.is_empty(),
            "{kind:?} {:?}",
            extrema.failures
        );
        assert!(
            extrema.pressure.is_some() && extrema.temperature.is_some(),
            "{kind:?} {extrema:?}"
        );
        for root in [extrema.pressure.unwrap(), extrema.temperature.unwrap()] {
            let p = &root.state;
            assert!(root.stationarity_residual.abs() < 1e-10);
            assert!(p.fugacity_residual < 1e-9);
            assert_eq!(p.liquid[2], 0.);
            assert_eq!(p.vapor[2], 0.);
            for sign in [-1., 1.] {
                let constraint = if root.kind == EnvelopeExtremum::Pressure {
                    EnvelopeConstraint::Temperature(p.temperature + sign * 0.01)
                } else {
                    EnvelopeConstraint::Pressure(p.pressure + sign * 1000.)
                };
                let neighbor = model.refine_envelope_state(p, constraint).unwrap();
                if root.kind == EnvelopeExtremum::Pressure {
                    assert!(neighbor.pressure < p.pressure);
                } else {
                    assert!(neighbor.temperature < p.temperature);
                }
            }
        }
    }
}

#[test]
fn truncated_curves_do_not_invent_extrema_and_coincident_or_pure_guesses_fail() {
    use thermopack_rs::SaturationKind;
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let curve = model
        .phase_envelope(
            &[0.5, 0.5],
            EnvelopeOptions {
                maximum_points: 3,
                ..Default::default()
            },
        )
        .unwrap();
    let extrema = model.envelope_extrema(&[0.5, 0.5], &curve).unwrap();
    assert!(
        extrema.pressure.is_none() && extrema.temperature.is_none() && extrema.failures.is_empty()
    );
    let mut guess = curve.points[0].clone();
    guess.liquid = vec![0.5, 0.5];
    guess.vapor = vec![0.5, 0.5];
    assert!(
        model
            .refine_envelope_extremum(&[0.5, 0.5], &guess, EnvelopeExtremum::Pressure)
            .is_err()
    );
    guess.kind = SaturationKind::Bubble;
    assert!(
        model
            .refine_envelope_extremum(&[1., 0.], &guess, EnvelopeExtremum::Temperature)
            .is_err()
    );
    guess.temperature = f64::NAN;
    assert!(
        model
            .refine_envelope_extremum(&[0.5, 0.5], &guess, EnvelopeExtremum::Temperature)
            .is_err()
    );
}
