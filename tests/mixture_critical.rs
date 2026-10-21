//! Exact-source critical coordinates and Hessians, plus independent derivatives.
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule, critical::critical_point};

fn kind(name: &str) -> CubicKind {
    match name {
        "PR" => CubicKind::PengRobinson,
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => panic!("{name}"),
    }
}

#[test]
fn source_critical_points_and_tv_hessians_match() {
    let mut count = 0;
    for row in include_str!("fixtures/mixture_critical.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let names: Vec<_> = c[2].split(',').collect();
        let z: Vec<f64> = c[3].split(',').map(|x| x.parse().unwrap()).collect();
        let ws = thermopack_rs::wong_sandler::WongSandlerData::new(names.len());
        let model = Cubic::with_options(
            kind(c[1]),
            &names,
            CubicOptions {
                volume_translation: c[4] == "translated",
                wong_sandler: (c[4] == "ws").then_some(&ws),
                mixing_rule: if c[4] == "ws" {
                    MixingRule::WongSandler
                } else {
                    MixingRule::VdW
                },
                ..Default::default()
            },
        )
        .unwrap();
        let expected: Vec<f64> = f[1..4].iter().map(|x| x.parse().unwrap()).collect();
        let point = critical_point(&model, &z).unwrap();
        for (i, actual) in [point.temperature, point.volume, point.pressure]
            .into_iter()
            .enumerate()
        {
            assert!(
                (actual - expected[i]).abs() <= [1e-6, 1e-12, 1e-2][i] + 1e-6 * expected[i].abs(),
                "{} {i}: {actual} {}",
                f[0],
                expected[i]
            );
        }
        let mut h = vec![0.0; z.len() * z.len()];
        model
            .prepare(expected[0], &z)
            .unwrap()
            .log_fugacity_tv_into(expected[1], &mut vec![0.0; z.len()], &mut h)
            .unwrap();
        let source_h: Vec<f64> = f[4].split(',').map(|x| x.parse().unwrap()).collect();
        let active: Vec<_> = z
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| (x > 0.0).then_some(i))
            .collect();
        let mut k = 0;
        for &i in &active {
            for &j in &active {
                assert!(
                    (h[i * z.len() + j] - source_h[k]).abs() < 1e-9 * (1.0 + source_h[k].abs()),
                    "{} hessian {i} {j}: {} {}",
                    f[0],
                    h[i * z.len() + j],
                    source_h[k]
                );
                k += 1;
            }
        }
        let native = model
            .mixture_critical_conditions(expected[0], expected[1], &z)
            .unwrap();
        let eigen: f64 = f[5].parse().unwrap();
        let third: f64 = f[6].parse().unwrap();
        assert!((native.eigenvalue - eigen).abs() < 1e-9, "{} eigen", f[0]);
        assert!(
            (native.third_derivative - third).abs() < 1e-8,
            "{} third {} {third}",
            f[0],
            native.third_derivative
        );
        count += 1;
    }
    assert_eq!(count, 197);
}

#[test]
fn critical_conditions_close_for_binary_and_ternary_mixtures() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for (names, z) in [
            (vec!["C1", "C2"], vec![0.5, 0.5]),
            (vec!["C1", "C2", "C3"], vec![0.2, 0.3, 0.5]),
            (vec!["C1", "C2", "C3"], vec![0.0, 0.3, 0.7]),
        ] {
            let model = Cubic::new(kind, &names).unwrap();
            let point = critical_point(&model, &z).unwrap();
            let c = model
                .mixture_critical_conditions(point.temperature, point.volume, &z)
                .unwrap();
            assert!(c.eigenvalue.abs() < 3e-9, "{kind:?} {z:?} {c:?}");
            assert!(c.third_derivative.abs() < 3e-9);
            assert!(point.pressure > 0.0);
            assert!(
                model
                    .prepare(point.temperature, &z)
                    .unwrap()
                    .pressure(point.volume)
                    .unwrap()
                    .dv
                    < 0.0
            );
        }
    }
}

#[test]
fn tv_hessian_matches_mole_differences_and_tp_log_fugacity() {
    for shift in [false, true] {
        let mut model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
        if shift {
            model.enable_volume_translation("DEFAULT").unwrap();
        }
        let z = [0.2, 0.3, 0.5];
        for temperature in [180.0, 300.0, 450.0] {
            let state = model.prepare(temperature, &z).unwrap();
            for v in [state.covolume() * 1.5, 0.001] {
                let mut f = [0.0; 3];
                let mut h = [0.0; 9];
                state.log_fugacity_tv_into(v, &mut f, &mut h).unwrap();
                let p = state.pressure(v).unwrap().value;
                if p > 0.0 {
                    let mut phi = [0.0; 3];
                    state.ln_fugacity_into(v, &mut phi).unwrap();
                    for i in 0..3 {
                        assert!((f[i] - phi[i] - (z[i] * p).ln()).abs() < 1e-11);
                    }
                }
                for j in 0..3 {
                    let step = 1e-5;
                    let mut up = z;
                    let mut down = z;
                    up[j] += step;
                    down[j] -= step;
                    let mut values = Vec::new();
                    for mut amounts in [down, up] {
                        let total: f64 = amounts.iter().sum();
                        for x in &mut amounts {
                            *x /= total;
                        }
                        let trial = model.prepare(temperature, &amounts).unwrap();
                        let mut ft = [0.0; 3];
                        trial
                            .log_fugacity_tv_into(v / total, &mut ft, &mut [0.0; 9])
                            .unwrap();
                        values.push(ft);
                    }
                    for i in 0..3 {
                        let numeric = (values[1][i] - values[0][i]) / (2.0 * step);
                        assert!(
                            (numeric - h[i * 3 + j]).abs() < 1e-6 * (1.0 + numeric.abs()),
                            "{shift} {temperature} {v} {i} {j}: {numeric} {}",
                            h[i * 3 + j]
                        );
                        assert!((h[i * 3 + j] - h[j * 3 + i]).abs() < 1e-13);
                    }
                }
            }
        }
    }
    // The TV expression is regular at a pure mechanical spinodal/critical point.
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2"]).unwrap();
    let point = critical_point(&model, &[1.0]).unwrap();
    let mut f = [0.0];
    let mut h = [0.0];
    model
        .prepare(point.temperature, &[1.0])
        .unwrap()
        .log_fugacity_tv_into(point.volume, &mut f, &mut h)
        .unwrap();
    assert!(h[0].abs() < 1e-8);
    assert!(f[0].is_finite());
}

#[test]
fn critical_direction_is_a_null_mode_and_changes_stability() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0.2, 0.3, 0.5];
    let point = critical_point(&model, &z).unwrap();
    let c = model
        .mixture_critical_conditions(point.temperature, point.volume, &z)
        .unwrap();
    let mut f = [0.0; 3];
    let mut h = [0.0; 9];
    model
        .prepare(point.temperature, &z)
        .unwrap()
        .log_fugacity_tv_into(point.volume, &mut f, &mut h)
        .unwrap();
    assert!((c.eigenvector.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-14);
    for i in 0..3 {
        let value = (0..3)
            .map(|j| z[i].sqrt() * h[i * 3 + j] * z[j].sqrt() * c.eigenvector[j])
            .sum::<f64>();
        assert!((value - c.eigenvalue * c.eigenvector[i]).abs() < 1e-13);
    }
    let below = model
        .mixture_critical_conditions(point.temperature - 0.1, point.volume, &z)
        .unwrap();
    let above = model
        .mixture_critical_conditions(point.temperature + 0.1, point.volume, &z)
        .unwrap();
    assert!(below.eigenvalue < 0.0 && above.eigenvalue > 0.0);
    // Independent second difference of the physical log fugacity, along the
    // eigenvector at fixed TOTAL volume, checks the analytic third derivative.
    for temperature in [point.temperature, point.temperature + 10.0] {
        let local = model
            .mixture_critical_conditions(temperature, point.volume, &z)
            .unwrap();
        let direction: Vec<_> = z
            .iter()
            .zip(&local.eigenvector)
            .map(|(z, u)| z.sqrt() * u)
            .collect();
        let step = 2e-4;
        let mut contractions = [0.0; 3];
        for (k, s) in [-step, 0.0, step].into_iter().enumerate() {
            let mut amounts: Vec<_> = z.iter().zip(&direction).map(|(z, d)| z + s * d).collect();
            let total: f64 = amounts.iter().sum();
            for a in &mut amounts {
                *a /= total;
            }
            model
                .prepare(temperature, &amounts)
                .unwrap()
                .log_fugacity_tv_into(point.volume / total, &mut f, &mut h)
                .unwrap();
            contractions[k] = f.iter().zip(&direction).map(|(f, d)| f * d).sum();
        }
        let numeric = (contractions[2] - 2.0 * contractions[1] + contractions[0]) / (step * step);
        assert!(
            (numeric - local.third_derivative).abs() < 2e-6 * (1.0 + local.third_derivative.abs()),
            "{numeric} {}",
            local.third_derivative
        );
    }
    // A mixture critical point is not its fixed-composition pressure inflection.
    let pseudo = model.pseudo_critical(&z).unwrap();
    assert!((pseudo.temperature - point.temperature).abs() > 1.0);
}

#[test]
fn translation_permutation_pure_limit_and_explicit_seeds() {
    use thermopack_rs::mixture_critical::critical_point_with_guess;
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0.2, 0.3, 0.5];
    let point = critical_point(&model, &z).unwrap();
    let perturbed =
        critical_point_with_guess(&model, &z, point.temperature * 1.03, point.volume * 1.1)
            .unwrap();
    assert!((perturbed.temperature - point.temperature).abs() < 1e-6);
    assert!((perturbed.volume / point.volume - 1.0).abs() < 1e-7);
    let reversed = Cubic::new(CubicKind::PengRobinson, &["C3", "C2", "C1"]).unwrap();
    let reordered = critical_point(&reversed, &[0.5, 0.3, 0.2]).unwrap();
    assert!((reordered.temperature - point.temperature).abs() < 1e-6);
    assert!((reordered.pressure / point.pressure - 1.0).abs() < 1e-8);
    let mut translated = model.clone();
    translated.enable_volume_translation("DEFAULT").unwrap();
    let shifted = critical_point(&translated, &z).unwrap();
    assert!((shifted.temperature - point.temperature).abs() < 1e-6);
    assert!((shifted.pressure / point.pressure - 1.0).abs() < 1e-8);
    let correction = translated
        .prepare(point.temperature, &z)
        .unwrap()
        .volume_translation()
        .value;
    assert!(((shifted.volume + correction) / point.volume - 1.0).abs() < 1e-7);
    let pure = critical_point(&model, &[0.0, 1.0, 0.0]).unwrap();
    let critical = critical_point_with_guess(
        &model,
        &[0.0, 1.0, 0.0],
        pure.temperature * 1.01,
        pure.volume * 1.03,
    )
    .unwrap();
    assert!((pure.temperature - critical.temperature).abs() < 1e-6);
    let c = model
        .mixture_critical_conditions(critical.temperature, critical.volume, &[0.0, 1.0, 0.0])
        .unwrap();
    assert_eq!(c.eigenvector, [0.0, 1.0, 0.0]);
    for invalid in [
        vec![0.5, 0.5],
        vec![0.2, -0.1, 0.9],
        vec![0.0; 3],
        vec![f64::NAN, 0.0, 1.0],
    ] {
        assert!(critical_point(&model, &invalid).is_err());
    }
    for (t, v) in [
        (0.0, point.volume),
        (point.temperature, 0.0),
        (f64::NAN, point.volume),
        (point.temperature, 1e-9),
    ] {
        assert!(critical_point_with_guess(&model, &z, t, v).is_err());
    }
    let activity = Cubic::with_options(
        CubicKind::PengRobinson,
        &["C1", "C3"],
        CubicOptions {
            mixing_rule: MixingRule::Unifac,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(matches!(
        critical_point(&activity, &[0.5, 0.5]),
        Err(thermopack_rs::Error::Unsupported(_))
    ));
}
