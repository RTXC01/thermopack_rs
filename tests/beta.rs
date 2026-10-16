use thermopack_rs::{
    Beta, BetaCorrelation, Cubic, CubicKind, CubicOptions, Phase, TranslationKind,
    VolumeTranslation,
};

#[test]
fn cryogenic_caloric_bracketing_resolves_a_nonmonotonic_curve() {
    let model = Cubic::with_options(
        CubicKind::PengRobinson,
        &["NE"],
        CubicOptions {
            parameter_reference: "QuantumCubic",
            ..Default::default()
        },
    )
    .unwrap();
    let options = thermopack_rs::CaloricFlashOptions {
        minimum_temperature: 0.5,
        maximum_temperature: 999.,
        ..Default::default()
    };
    let state = model.prepare(8., &[1.]).unwrap();
    let v = state.volume(1e3, Phase::Liquid).unwrap();
    let h = state.enthalpy(v).unwrap();
    let s = state.entropy(v).unwrap();
    for flash in [
        model.ph_flash_with_options(1e3, &[1.], h, options),
        model.ps_flash_with_options(1e3, &[1.], s, options),
    ] {
        assert!((flash.unwrap().temperature - 8.).abs() < 1e-6);
    }
    // The fork selects the warm closing root for a target also attained at 3 K.
    let ph = model
        .ph_flash_with_options(1e6, &[1.], -7733.87431208972, options)
        .unwrap();
    assert!((ph.temperature - 14.215322191045736).abs() < 1e-6);
    let props = model.equilibrium_properties(&ph).unwrap();
    assert!((props.enthalpy + 7733.87431208972).abs() < 1e-7);
}

#[test]
fn quantum_helper_uses_combined_references_and_owned_bounds() {
    let mut model = Cubic::quantum(&["NE"]).unwrap();
    assert_eq!(model.temperature_bounds(), (2., 999.));
    assert_eq!(model.caloric_options().minimum_temperature, 2.);
    let plain = Cubic::new(CubicKind::PengRobinson, &["NE"]).unwrap();
    assert_eq!(plain.temperature_bounds(), (80., 999.));
    assert!(model.set_temperature_bounds(30., 20.).is_err());
    assert_eq!(model.temperature_bounds(), (2., 999.));
    model.set_temperature_bounds(0.5, 500.).unwrap();
    assert_eq!(model.caloric_options().maximum_temperature, 500.);
    let state = model.tp_flash(20., 1e5, &[1.]).unwrap();
    let h = model.equilibrium_properties(&state).unwrap().enthalpy;
    assert!((model.ph_flash(1e5, &[1.], h).unwrap().temperature - 20.).abs() < 1e-6);
}

#[test]
fn extreme_wilson_estimates_and_absent_components_keep_rr_endpoints() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "CO2"]).unwrap();
    assert!(model.wilson_k(0.5, 1e5).unwrap().iter().all(|k| *k == 0.));
    let beta =
        thermopack_rs::numerics::rachford_rice(&[0., 0.4, 0.6], &[1e-300, 1e300, 1e300]).unwrap();
    assert_eq!(beta, 1.);
}

#[test]
fn beta_derivatives_normalization_and_configuration() {
    for params in [
        [0., 0.],
        [1.4912, 3.2634],
        [3.0696, 12.682],
        [0.4673, 2.4634],
    ] {
        let beta = Beta::new(BetaCorrelation::Quantum, params).unwrap();
        let tc = 33.145;
        assert!((beta.evaluate(tc, tc).unwrap().value - 1.).abs() < 1e-15);
        for t in [2., 10., 80., 300.] {
            let dt = t * 1e-4;
            let value = beta.evaluate(t, tc).unwrap();
            let plus = beta.evaluate(t + dt, tc).unwrap();
            let minus = beta.evaluate(t - dt, tc).unwrap();
            assert!(
                ((plus.value - minus.value) / (2. * dt) - value.dt).abs()
                    < 1e-10 + 1e-6 * value.dt.abs()
            );
            assert!(
                ((plus.dt - minus.dt) / (2. * dt) - value.dtt).abs()
                    < 1e-11 + 1e-6 * value.dtt.abs()
            );
        }
    }
    assert!(Beta::new(BetaCorrelation::Quantum, [f64::NAN, 0.]).is_err());
    assert!(
        Beta::new(BetaCorrelation::Quantum, [1., -5.])
            .unwrap()
            .evaluate(5., 30.)
            .is_err()
    );
    let quantum = Cubic::with_options(
        CubicKind::PengRobinson,
        &["HE", "H2"],
        CubicOptions {
            parameter_reference: "QuantumCubic",
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(quantum.components()[0].critical_temperature, 5.1953);
    assert_eq!(quantum.components()[1].critical_temperature, 33.145);
    assert_eq!(quantum.beta(1).unwrap().parameters, [3.0696, 12.682]);
    assert_eq!(quantum.components()[0].molar_mass, 4.002602);
    assert_eq!(thermopack_rs::component("HE").unwrap().molar_mass, 4.003);
    assert_eq!(
        thermopack_rs::component("HE").unwrap().critical_temperature,
        5.1953
    );
}

#[test]
fn variable_covolume_pressure_caloric_and_fugacity_derivatives() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let mut model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        model.set_kij(0, 1, 0.07).unwrap();
        model.set_lij(1, 2, 0.04).unwrap();
        for i in 0..3 {
            model
                .set_beta(
                    i,
                    Beta::new(BetaCorrelation::Quantum, [3. + i as f64, 12.]).unwrap(),
                )
                .unwrap();
            model
                .set_volume_translation(
                    i,
                    VolumeTranslation::new(
                        TranslationKind::Quadratic,
                        [1e-6, 1e-9 * (i + 1) as f64, 1e-12, 0., 0., 0.],
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        for shifted in [false, true] {
            model.set_volume_translation_enabled(shifted);
            for t in [250., 400.] {
                for p in [1e5, 5e6] {
                    let z = [0.5, 0.3, 0.2];
                    let state = model.prepare(t, &z).unwrap();
                    let v = state.volume(p, Phase::Vapor).unwrap();
                    let dt = 0.001;
                    let plus = model.prepare(t + dt, &z).unwrap();
                    let minus = model.prepare(t - dt, &z).unwrap();
                    let pt = (plus.pressure(v).unwrap().value - minus.pressure(v).unwrap().value)
                        / (2. * dt);
                    let analytic = state.pressure(v).unwrap().dt;
                    assert!(
                        (pt - analytic).abs() < 1e-6 * analytic.abs(),
                        "{kind:?} shifted={shifted} T={t}"
                    );
                    let hp = plus
                        .enthalpy(plus.volume(p, Phase::Vapor).unwrap())
                        .unwrap();
                    let hm = minus
                        .enthalpy(minus.volume(p, Phase::Vapor).unwrap())
                        .unwrap();
                    let up = plus.enthalpy(v).unwrap() - plus.pressure(v).unwrap().value * v;
                    let um = minus.enthalpy(v).unwrap() - minus.pressure(v).unwrap().value * v;
                    let (cp, cv) = state.heat_capacities(v).unwrap();
                    assert!(
                        ((hp - hm) / (2. * dt) - cp).abs() < 1e-5 * cp.abs(),
                        "Cp {kind:?}"
                    );
                    assert!(
                        ((up - um) / (2. * dt) - cv).abs() < 1e-5 * cv.abs(),
                        "Cv {kind:?}"
                    );
                    let (mut ft, mut fp, mut fnn) = ([0.; 3], [0.; 3], [0.; 9]);
                    state
                        .fugacity_derivatives_into(v, &mut ft, &mut fp, &mut fnn)
                        .unwrap();
                    let (mut phip, mut phim) = ([0.; 3], [0.; 3]);
                    plus.ln_fugacity_tp_into(p, plus.volume(p, Phase::Vapor).unwrap(), &mut phip)
                        .unwrap();
                    minus
                        .ln_fugacity_tp_into(p, minus.volume(p, Phase::Vapor).unwrap(), &mut phim)
                        .unwrap();
                    for i in 0..3 {
                        assert!(
                            ((phip[i] - phim[i]) / (2. * dt) - ft[i]).abs()
                                < 1e-9 + 1e-5 * ft[i].abs(),
                            "fugacity T {kind:?} component {i}"
                        );
                    }
                    let dp = p * 1e-5;
                    state
                        .ln_fugacity_tp_into(
                            p + dp,
                            state.volume(p + dp, Phase::Vapor).unwrap(),
                            &mut phip,
                        )
                        .unwrap();
                    state
                        .ln_fugacity_tp_into(
                            p - dp,
                            state.volume(p - dp, Phase::Vapor).unwrap(),
                            &mut phim,
                        )
                        .unwrap();
                    for i in 0..3 {
                        assert!(
                            ((phip[i] - phim[i]) / (2. * dp) - fp[i]).abs()
                                < 1e-12 + 1e-5 * fp[i].abs()
                        );
                    }
                    for j in 0..3 {
                        let dn = 1e-5;
                        let zp: Vec<_> = z
                            .iter()
                            .enumerate()
                            .map(|(i, x)| (x + if i == j { dn } else { 0. }) / (1. + dn))
                            .collect();
                        let zm: Vec<_> = z
                            .iter()
                            .enumerate()
                            .map(|(i, x)| (x - if i == j { dn } else { 0. }) / (1. - dn))
                            .collect();
                        let sp = model.prepare(t, &zp).unwrap();
                        let sm = model.prepare(t, &zm).unwrap();
                        sp.ln_fugacity_tp_into(p, sp.volume(p, Phase::Vapor).unwrap(), &mut phip)
                            .unwrap();
                        sm.ln_fugacity_tp_into(p, sm.volume(p, Phase::Vapor).unwrap(), &mut phim)
                            .unwrap();
                        for i in 0..3 {
                            assert!(
                                ((phip[i] - phim[i]) / (2. * dn) - fnn[i * 3 + j]).abs()
                                    < 1e-7 + 1e-5 * fnn[i * 3 + j].abs()
                            );
                        }
                    }
                }
            }
        }
    }
}
