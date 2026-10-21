use thermopack_rs::{Cubic, CubicKind, MixtureSaturation};

#[test]
fn all_qualified_source_bubble_and_dew_states_match() {
    let mut count = 0;
    for row in include_str!("fixtures/mixture_saturation.csv")
        .lines()
        .skip(1)
    {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let options: Vec<_> = f[5].split('@').collect();
        let kind = match f[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            _ => CubicKind::VanDerWaals,
        };
        let model = Cubic::with_options(
            kind,
            &f[1].split(',').collect::<Vec<_>>(),
            thermopack_rs::CubicOptions {
                volume_translation: options[0].starts_with("V_"),
                mixing_rule: if options.contains(&"unifac") {
                    thermopack_rs::MixingRule::Unifac
                } else {
                    thermopack_rs::MixingRule::VdW
                },
                ..Default::default()
            },
        )
        .unwrap();
        let known = f[2].parse().unwrap();
        let z: Vec<f64> = f[4].split(',').map(|s| s.parse().unwrap()).collect();
        let result = match options[0].trim_start_matches("V_") {
            "MIX_BP" => model.bubble_pressure(known, &z),
            "MIX_DP" => model.dew_pressure(known, &z),
            "MIX_BT" => model.bubble_temperature(known, &z),
            _ => model.dew_temperature(known, &z),
        }
        .unwrap_or_else(|e| panic!("{}: {e}", c[0]));
        for (i, a) in [
            result.temperature,
            result.pressure,
            result.liquid_volume,
            result.vapor_volume,
        ]
        .into_iter()
        .enumerate()
        {
            close(
                a,
                c[i + 1].parse().unwrap(),
                [1e-6, 1e-2, 1e-12, 1e-12][i],
                2e-7,
            );
        }
        for (actual, expected) in [(&result.liquid, c[5]), (&result.vapor, c[6])] {
            let expected: Vec<f64> = expected.split(',').map(|v| v.parse().unwrap()).collect();
            assert_eq!(actual.len(), expected.len());
            for (&a, b) in actual.iter().zip(expected) {
                close(a, b, 1e-9, 2e-7);
            }
        }
        check(&model, &result);
        count += 1;
    }
    assert_eq!(count, 959);
}

fn close(a: f64, b: f64, abs: f64, rel: f64) {
    assert!(
        (a - b).abs() <= abs + rel * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn check(model: &Cubic, s: &MixtureSaturation) {
    close(s.liquid.iter().sum(), 1., 1e-12, 0.);
    close(s.vapor.iter().sum(), 1., 1e-12, 0.);
    assert!(s.vapor_volume > s.liquid_volume);
    let l = model.prepare(s.temperature, &s.liquid).unwrap();
    let v = model.prepare(s.temperature, &s.vapor).unwrap();
    close(
        l.pressure(s.liquid_volume).unwrap().value,
        s.pressure,
        1e-3,
        1e-8,
    );
    close(
        v.pressure(s.vapor_volume).unwrap().value,
        s.pressure,
        1e-3,
        1e-8,
    );
    assert!(l.pressure(s.liquid_volume).unwrap().dv < 0.);
    assert!(v.pressure(s.vapor_volume).unwrap().dv < 0.);
    let mut fl = vec![0.; s.liquid.len()];
    let mut fv = fl.clone();
    l.ln_fugacity_tp_into(s.pressure, s.liquid_volume, &mut fl)
        .unwrap();
    v.ln_fugacity_tp_into(s.pressure, s.vapor_volume, &mut fv)
        .unwrap();
    let mut residual: f64 = 0.;
    for i in 0..fl.len() {
        if s.liquid[i] > 0. && s.vapor[i] > 0. {
            residual = residual.max((s.liquid[i].ln() + fl[i] - s.vapor[i].ln() - fv[i]).abs());
        } else {
            assert_eq!(s.liquid[i], 0.);
            assert_eq!(s.vapor[i], 0.);
        }
    }
    assert!(residual < 3e-10);
    close(residual, s.fugacity_residual, 1e-12, 1e-5);
}

#[test]
fn saturation_tangents_and_inverse_specifications_close_independently() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2"]).unwrap();
        let z = [0.5, 0.5];
        let dt = 0.001;
        for bubble in [true, false] {
            let calc = |t| {
                if bubble {
                    model.bubble_pressure(t, &z)
                } else {
                    model.dew_pressure(t, &z)
                }
            };
            let s = calc(180.).unwrap();
            check(&model, &s);
            let finite =
                (calc(180. + dt).unwrap().pressure - calc(180. - dt).unwrap().pressure) / (2. * dt);
            close(s.dpdt, finite, 1e-3, 1e-7);
            let inverse = if bubble {
                model.bubble_temperature(s.pressure, &z)
            } else {
                model.dew_temperature(s.pressure, &z)
            }
            .unwrap();
            check(&model, &inverse);
            close(inverse.temperature, 180., 1e-6, 1e-9);
            close(inverse.dpdt, s.dpdt, 1e-4, 1e-7);
        }
    }
}

#[test]
fn explicit_guesses_can_follow_both_retrograde_pressure_branches() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0.7, 0.2, 0.1];
    let first = model.dew_pressure(270., &z).unwrap();
    let mut upper = first.clone();
    upper.pressure = 7273058.052139043;
    upper.liquid = vec![0.4719084073363501, 0.2841578027626974, 0.24393378990095263];
    let other = model.refine_saturation_pressure(270., &upper).unwrap();
    check(&model, &first);
    check(&model, &other);
    close(first.pressure, 5464939.50895043, 1., 2e-7);
    assert!(other.pressure > 7e6 && other.dpdt < 0. && first.dpdt > 0.);
    let dt = 0.001;
    let hi = model.refine_saturation_pressure(270. + dt, &other).unwrap();
    let lo = model.refine_saturation_pressure(270. - dt, &other).unwrap();
    close(
        other.dpdt,
        (hi.pressure - lo.pressure) / (2. * dt),
        1e-2,
        1e-6,
    );
}

#[test]
fn high_pressure_temperature_seed_has_a_nontrivial_solution() {
    let model = Cubic::new(CubicKind::VanDerWaals, &["CO2", "C3"]).unwrap();
    let p = 5719761.632444506;
    let z = [0.8, 0.2];
    let s = model.dew_temperature_with_guess(p, &z, 310.).unwrap();
    check(&model, &s);
    close(s.temperature, 310., 1e-5, 1e-8);
    let fixed = model.prepare(s.temperature, &s.vapor).unwrap();
    let trial = model.prepare(s.temperature, &s.liquid).unwrap();
    let (mut ft, mut fp, mut fnn) = (vec![0.; 2], vec![0.; 2], vec![0.; 4]);
    let mut tt = ft.clone();
    let mut tp = fp.clone();
    let mut tn = fnn.clone();
    fixed
        .fugacity_derivatives_into(s.vapor_volume, &mut ft, &mut fp, &mut fnn)
        .unwrap();
    trial
        .fugacity_derivatives_into(s.liquid_volume, &mut tt, &mut tp, &mut tn)
        .unwrap();
    let slope: f64 = s
        .liquid
        .iter()
        .enumerate()
        .map(|(i, v)| v * s.temperature * (ft[i] - tt[i]))
        .sum();
    assert!(slope < 0.);
    let default = model.dew_temperature(p, &z).unwrap();
    check(&model, &default);
    close(default.temperature, s.temperature, 1e-5, 1e-8);
}

#[test]
fn absent_components_pure_limits_and_invalid_inputs_are_explicit() {
    let m = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0., 1., 0.];
    let pure = Cubic::new(CubicKind::PengRobinson, &["C2"])
        .unwrap()
        .saturation_pressure(200.)
        .unwrap();
    for s in [
        m.bubble_pressure(200., &z).unwrap(),
        m.dew_pressure(200., &z).unwrap(),
    ] {
        check(&m, &s);
        assert_eq!(s.liquid, z);
        assert_eq!(s.vapor, z);
        close(s.pressure, pure.pressure, 1e-4, 1e-8);
    }
    let mixed = m.dew_pressure(180., &[0.5, 0.5, 0.]).unwrap();
    check(&m, &mixed);
    let next = m
        .refine_saturation_temperature(mixed.pressure, &mixed)
        .unwrap();
    check(&m, &next);
    close(next.temperature, 180., 1e-6, 1e-8);
    for bad in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(m.bubble_pressure(bad, &z).is_err());
        assert!(m.dew_temperature(bad, &z).is_err());
        assert!(m.bubble_pressure_with_guess(200., &z, bad).is_err());
    }
    for bad in [
        &[0.5, 0.5][..],
        &[0.5, 0.4, 0.][..],
        &[-0.1, 0.6, 0.5][..],
        &[0., 0., 0.][..],
    ] {
        assert!(m.bubble_pressure(200., bad).is_err());
    }
    let mut invalid = mixed.clone();
    invalid.liquid[2] = 1e-6;
    assert!(m.refine_saturation_pressure(180., &invalid).is_err());
}

#[test]
fn inverse_dew_solve_rejects_the_sources_coincident_phase_false_convergence() {
    let model = Cubic::new(CubicKind::VanDerWaals, &["C1", "C2"]).unwrap();
    let z = [0.8, 0.2];
    let forward = model.dew_pressure(220., &z).unwrap();
    close(forward.pressure, 4723572.124716964, 1e-2, 2e-7);
    let inverse = model.dew_temperature(forward.pressure, &z).unwrap();
    close(inverse.temperature, 220., 1e-5, 1e-8);
    check(&model, &inverse);
    assert!((inverse.vapor_volume - inverse.liquid_volume) / inverse.vapor_volume > 0.1);
}
