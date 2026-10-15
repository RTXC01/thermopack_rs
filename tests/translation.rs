use thermopack_rs::{Cubic, CubicKind, CubicOptions, Phase, R, TranslationKind, VolumeTranslation};

fn translated(kind: CubicKind, correlation: TranslationKind) -> Cubic {
    let mut model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
    for i in 0..3 {
        let scale = 1. + 0.2 * i as f64;
        let coefficients = [2e-6, 1e-8, -1e-11, 1e-14, -1e-17, 1e-20].map(|v| v * scale);
        model
            .set_volume_translation(
                i,
                VolumeTranslation::new(correlation, coefficients).unwrap(),
            )
            .unwrap();
    }
    model.set_volume_translation_enabled(true);
    model
}

#[test]
fn translated_values_and_derivatives_obey_thermodynamic_identities() {
    let z = [0.5, 0.3, 0.2];
    let t = 300.;
    let p = 5e6;
    let dt = 0.001;
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for correlation in [
            TranslationKind::Constant,
            TranslationKind::Linear,
            TranslationKind::Quadratic,
            TranslationKind::Quintic,
        ] {
            let model = translated(kind, correlation);
            let plain = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
            let state = model.prepare(t, &z).unwrap();
            let base = plain.prepare(t, &z).unwrap();
            let shift = state.volume_translation();
            let ve = base.volume(p, Phase::Vapor).unwrap();
            let v = state.volume(p, Phase::Vapor).unwrap();
            assert!((v + shift.value - ve).abs() < 1e-14);
            assert!(
                (state.enthalpy(v).unwrap()
                    - base.enthalpy(ve).unwrap()
                    - p * (t * shift.dt - shift.value))
                    .abs()
                    < 1e-7
            );
            assert!(
                (state.entropy(v).unwrap() - base.entropy(ve).unwrap() - p * shift.dt).abs() < 1e-9
            );
            let (cp, cv) = state.heat_capacities(v).unwrap();
            let plus = model.prepare(t + dt, &z).unwrap();
            let minus = model.prepare(t - dt, &z).unwrap();
            let hp = plus
                .enthalpy(plus.volume(p, Phase::Vapor).unwrap())
                .unwrap();
            let hm = minus
                .enthalpy(minus.volume(p, Phase::Vapor).unwrap())
                .unwrap();
            assert!(
                ((hp - hm) / (2. * dt) - cp).abs() < 1e-5 * cp.abs(),
                "Cp {kind:?} {correlation:?}"
            );
            let up = plus.enthalpy(v).unwrap() - plus.pressure(v).unwrap().value * v;
            let um = minus.enthalpy(v).unwrap() - minus.pressure(v).unwrap().value * v;
            assert!(
                ((up - um) / (2. * dt) - cv).abs() < 1e-5 * cv.abs(),
                "Cv {kind:?} {correlation:?}"
            );
            let at = state.pressure(v).unwrap();
            let pt =
                (plus.pressure(v).unwrap().value - minus.pressure(v).unwrap().value) / (2. * dt);
            assert!((pt - at.dt).abs() < 1e-6 * at.dt.abs());
            let (mut phi, mut phi0, mut ft, mut fp, mut fnn) = (
                vec![0.; 3],
                vec![0.; 3],
                vec![0.; 3],
                vec![0.; 3],
                vec![0.; 9],
            );
            state.ln_fugacity_tp_into(p, v, &mut phi).unwrap();
            base.ln_fugacity_tp_into(p, ve, &mut phi0).unwrap();
            state
                .fugacity_derivatives_into(v, &mut ft, &mut fp, &mut fnn)
                .unwrap();
            let mut phip = [0.; 3];
            let mut phim = [0.; 3];
            plus.ln_fugacity_tp_into(p, plus.volume(p, Phase::Vapor).unwrap(), &mut phip)
                .unwrap();
            minus
                .ln_fugacity_tp_into(p, minus.volume(p, Phase::Vapor).unwrap(), &mut phim)
                .unwrap();
            for i in 0..3 {
                let ci = model
                    .volume_translation(i)
                    .unwrap()
                    .evaluate(t)
                    .unwrap()
                    .value;
                assert!((phi[i] - phi0[i] + p * ci / (R * t)).abs() < 1e-10);
                assert!(
                    ((phip[i] - phim[i]) / (2. * dt) - ft[i]).abs() < 1e-8 + 1e-6 * ft[i].abs()
                );
            }
            let mut zero_dt = [0.; 3];
            let mut zero_dp = [0.; 3];
            let mut base_nn = [0.; 9];
            base.fugacity_derivatives_into(ve, &mut zero_dt, &mut zero_dp, &mut base_nn)
                .unwrap();
            for (a, b) in fnn.iter().zip(base_nn) {
                assert!((a - b).abs() < 1e-9);
            }
        }
    }
}

#[test]
fn translation_configuration_and_complete_cubic_constructor() {
    let mut model = Cubic::new(CubicKind::PengRobinson, &["CO2"]).unwrap();
    model
        .set_volume_translation(
            0,
            VolumeTranslation::new(TranslationKind::Constant, [1e-5, 0., 0., 0., 0., 0.]).unwrap(),
        )
        .unwrap();
    assert!(!model.volume_translation_enabled());
    let t = 250.;
    let p = 1e6;
    let original = model
        .prepare(t, &[1.])
        .unwrap()
        .volume(p, Phase::Vapor)
        .unwrap();
    model.set_volume_translation_enabled(true);
    assert!(
        (model
            .prepare(t, &[1.])
            .unwrap()
            .volume(p, Phase::Vapor)
            .unwrap()
            + 1e-5
            - original)
            .abs()
            < 1e-14
    );
    model.disable_volume_translation();
    assert_eq!(model.volume_translation(0).unwrap().coefficients[0], 1e-5);
    assert_eq!(
        model
            .prepare(t, &[1.])
            .unwrap()
            .volume(p, Phase::Vapor)
            .unwrap(),
        original
    );
    let tcpr = Cubic::with_options(
        CubicKind::PengRobinson,
        &["C1", "CO2"],
        CubicOptions {
            parameter_reference: "tcPR",
            ..Default::default()
        },
    )
    .unwrap();
    assert!(tcpr.volume_translation_enabled());
    assert_eq!(
        tcpr.alpha(0).unwrap().correlation,
        thermopack_rs::AlphaCorrelation::Twu
    );
    assert!(
        Cubic::with_options(
            CubicKind::PengRobinson,
            &["HE"],
            CubicOptions {
                parameter_reference: "QuantumCubic",
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn translated_saturation_and_caloric_quality_preserve_phase_equilibrium() {
    for kind in [CubicKind::PengRobinson, CubicKind::SoaveRedlichKwong] {
        let base = Cubic::new(kind, &["CO2"]).unwrap();
        let mut model = base.clone();
        model
            .set_volume_translation(
                0,
                VolumeTranslation::new(TranslationKind::Quadratic, [2e-6, 1e-8, 1e-11, 0., 0., 0.])
                    .unwrap(),
            )
            .unwrap();
        model.set_volume_translation_enabled(true);
        let sat = model.saturation_pressure(245.).unwrap();
        let plain = base.saturation_pressure(245.).unwrap();
        assert!((sat.pressure / plain.pressure - 1.).abs() < 1e-12);
        let state = model.prepare(245., &[1.]).unwrap();
        let h = 0.6 * state.enthalpy(sat.liquid_volume).unwrap()
            + 0.4 * state.enthalpy(sat.vapor_volume).unwrap();
        let flash = model.ph_flash(sat.pressure, &[1.], h).unwrap();
        assert!((flash.temperature - 245.).abs() < 1e-6);
        assert!((flash.vapor_fraction - 0.4).abs() < 1e-7);
        let mixture = translated(kind, TranslationKind::Quintic);
        let z = [0.5, 0.3, 0.2];
        let t = 220.;
        let p = 1e6;
        let props = mixture
            .equilibrium_properties(&mixture.tp_flash(t, p, &z).unwrap())
            .unwrap();
        assert!((mixture.ph_flash(p, &z, props.enthalpy).unwrap().temperature - t).abs() < 1e-6);
        assert!((mixture.ps_flash(p, &z, props.entropy).unwrap().temperature - t).abs() < 1e-6);
    }
}
