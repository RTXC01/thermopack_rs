use thermopack_rs::{Alpha, AlphaCorrelation as C, Cubic, CubicKind};

#[test]
fn alpha_variants_match_eos_properties_and_fugacity_derivatives() {
    let absolute = [1e-14, 1e-2, 1e-6, 1e-2, 1e-7, 1e-9, 1e-7, 1e-9, 1e-7];
    let mut count = 0;
    for row in include_str!("fixtures/alpha_properties.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let config: Vec<_> = fields[5].split('@').collect();
        let alpha_case = format!(
            "{}|{}|{}|{}|{}|ALPHA",
            fields[0],
            fields[1],
            fields[2],
            config[1],
            config[2..].join("@")
        );
        let m = model(&alpha_case);
        let z = fields[4]
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect::<Vec<f64>>();
        let n = z.len();
        let state = m.prepare(fields[2].parse().unwrap(), &z).unwrap();
        let p = fields[3].parse().unwrap();
        let v = state
            .volume(
                p,
                if config[0].ends_with('L') {
                    thermopack_rs::Phase::Liquid
                } else {
                    thermopack_rs::Phase::Vapor
                },
            )
            .unwrap();
        let pressure = state.pressure(v).unwrap();
        let residual = state.residual(v).unwrap();
        let mut actual = vec![
            v,
            pressure.value,
            pressure.dt,
            pressure.dv,
            residual.enthalpy,
            residual.entropy,
            state.enthalpy(v).unwrap(),
            state.entropy(v).unwrap(),
            state.baseline_enthalpy_derivative_tp(v).unwrap(),
        ];
        let (mut phi, mut dt, mut dp, mut dn) =
            (vec![0.; n], vec![0.; n], vec![0.; n], vec![0.; n * n]);
        state.ln_fugacity_into(v, &mut phi).unwrap();
        state
            .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
            .unwrap();
        actual.extend(phi);
        actual.extend(dt);
        actual.extend(dp);
        actual.extend(dn);
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let abs = if i < 9 {
                absolute[i]
            } else if i < 9 + 2 * n {
                1e-9
            } else if i < 9 + 3 * n {
                1e-11
            } else {
                1e-8
            };
            let rel = if i == 0 {
                2e-9
            } else if i < 9 + n {
                2e-8
            } else {
                2e-7
            };
            assert!(
                (a - b).abs() <= abs + rel * a.abs().max(b.abs()),
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 3240);
}

#[test]
fn configured_alpha_derivatives_and_caloric_flashes() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for (correlation, params) in [
            (C::Twu, vec![0.2, 0.9, 2.1]),
            (C::MathiasCopeman, vec![0.8, 0.2, -0.1]),
            (C::SchmidtWenzel, vec![0.45]),
        ] {
            let mut pure = Cubic::new(kind, &["CO2"]).unwrap();
            pure.set_alpha(0, correlation, &params).unwrap();
            let sat = pure.saturation_pressure(245.).unwrap();
            let state = pure.prepare(245., &[1.]).unwrap();
            let h = 0.6 * state.enthalpy(sat.liquid_volume).unwrap()
                + 0.4 * state.enthalpy(sat.vapor_volume).unwrap();
            let s = 0.6 * state.entropy(sat.liquid_volume).unwrap()
                + 0.4 * state.entropy(sat.vapor_volume).unwrap();
            for flash in [
                pure.ph_flash(sat.pressure, &[1.], h).unwrap(),
                pure.ps_flash(sat.pressure, &[1.], s).unwrap(),
            ] {
                assert!((flash.temperature - 245.).abs() < 1e-6);
                assert!((flash.vapor_fraction - 0.4).abs() < 1e-7);
            }
            let mut m = Cubic::new(kind, &["CO2", "C1"]).unwrap();
            for i in 0..2 {
                m.set_alpha(i, correlation, &params).unwrap();
            }
            m.set_kij(0, 1, 0.12).unwrap();
            m.set_lij(0, 1, 0.03).unwrap();
            let z = [0.8, 0.2];
            let t = 300.;
            let p = 2e6;
            let state = m.prepare(t, &z).unwrap();
            let v = state.volume(p, thermopack_rs::Phase::Vapor).unwrap();
            let at = state.pressure(v).unwrap();
            let dt = 0.001;
            let plus = m.prepare(t + dt, &z).unwrap();
            let minus = m.prepare(t - dt, &z).unwrap();
            let fd =
                (plus.pressure(v).unwrap().value - minus.pressure(v).unwrap().value) / (2. * dt);
            assert!((fd - at.dt).abs() < 1e-6 * at.dt.abs());
            let hp = plus
                .enthalpy(plus.volume(p, thermopack_rs::Phase::Vapor).unwrap())
                .unwrap();
            let hm = minus
                .enthalpy(minus.volume(p, thermopack_rs::Phase::Vapor).unwrap())
                .unwrap();
            let cp = state.heat_capacities(v).unwrap().0;
            assert!(((hp - hm) / (2. * dt) - cp).abs() < 1e-6 * cp.abs());
            let feed = m.tp_flash(t, p, &z).unwrap();
            let properties = m.equilibrium_properties(&feed).unwrap();
            for flash in [
                m.ph_flash(p, &z, properties.enthalpy).unwrap(),
                m.ps_flash(p, &z, properties.entropy).unwrap(),
            ] {
                assert!(
                    (flash.temperature - t).abs() < 1e-6,
                    "{kind:?} {correlation:?}"
                );
            }
        }
    }
}

fn model(case: &str) -> Cubic {
    let fields: Vec<_> = case.split('|').collect();
    let kind = match fields[0] {
        "PR" => CubicKind::PengRobinson,
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => panic!("model"),
    };
    let config: Vec<_> = fields[4].split('@').collect();
    let mut model = Cubic::with_alpha(
        kind,
        &fields[1].split(',').collect::<Vec<_>>(),
        C::from_name(fields[3]).unwrap(),
        config[0],
    )
    .unwrap();
    for setter in &config[1..] {
        let fields: Vec<_> = setter.split(':').collect();
        let params = fields[2]
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect::<Vec<_>>();
        model
            .set_alpha(
                fields[0].parse().unwrap(),
                C::from_name(fields[1]).unwrap(),
                &params,
            )
            .unwrap();
    }
    model
}

#[test]
fn alpha_initialization_parameters_and_derivatives_match_source() {
    let absolute = [0., 1e-14, 1e-14, 1e-14, 1e-13, 1e-13, 1e-16];
    let mut count = 0;
    for row in include_str!("fixtures/alpha.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let m = model(columns[0]);
        let mut actual = Vec::new();
        for (i, component) in m.components().iter().enumerate() {
            let alpha = m.alpha(i).unwrap();
            let value = alpha
                .evaluate(fields[2].parse().unwrap(), component.critical_temperature)
                .unwrap();
            actual.push(alpha.correlation as i32 as f64);
            actual.extend(alpha.parameters);
            actual.extend([value.value, value.dt, value.dtt]);
        }
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i % 7 == 0 {
                0.
            } else {
                absolute[i % 7] + 2e-12 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 4232);
}

#[test]
fn analytic_alpha_derivatives_and_critical_limits() {
    for name in [
        "CLASSIC",
        "TWU",
        "MC",
        "CLASSICFIT",
        "GD",
        "RK",
        "SOAVE",
        "PR",
        "PT",
        "SW",
        "VDW",
        "PR78",
    ] {
        let correlation = C::from_name(name).unwrap();
        let params = match correlation {
            C::Twu => vec![0.2, 0.9, 2.1],
            C::MathiasCopeman => vec![0.8, 0.2, -0.1],
            C::SchmidtWenzel => vec![0.45],
            _ => vec![0.7],
        };
        let alpha = Alpha::new(correlation, &params).unwrap();
        let tc = 304.2;
        assert!((alpha.evaluate(tc, tc).unwrap().value - 1.).abs() < 1e-14);
        for tr in [0.3, 0.7, 1.2, 2.] {
            let t = tr * tc;
            let dt = t * 1e-5;
            let at = alpha.evaluate(t, tc).unwrap();
            let plus = alpha.evaluate(t + dt, tc).unwrap();
            let minus = alpha.evaluate(t - dt, tc).unwrap();
            assert!(
                ((plus.value - minus.value) / (2. * dt) - at.dt).abs()
                    < 1e-9 / tc + 1e-7 * at.dt.abs(),
                "{name} {tr} dt"
            );
            assert!(
                ((plus.dt - minus.dt) / (2. * dt) - at.dtt).abs()
                    < 1e-9 / (tc * tc) + 1e-7 * at.dtt.abs(),
                "{name} {tr} dtt"
            );
        }
    }
}

#[test]
fn source_fallbacks_and_unsupported_paths_are_explicit() {
    let pressure_only =
        Cubic::with_alpha(CubicKind::PengRobinson, &["BUTANAL"], C::Twu, "tcPR").unwrap();
    let state = pressure_only.prepare(400., &[1.]).unwrap();
    let volume = state.volume(1e6, thermopack_rs::Phase::Vapor).unwrap();
    assert!(state.pressure(volume).unwrap().value.is_finite());
    assert!(state.enthalpy(volume).is_err());
    let gerg = Cubic::with_alpha(CubicKind::PengRobinson, &["C1"], C::Gerg, "DEFAULT").unwrap();
    assert_eq!(gerg.alpha(0).unwrap().correlation, C::PengRobinson);
    assert!(gerg.alpha(0).unwrap().used_fallback);
    let explicit =
        Cubic::with_alpha(CubicKind::PengRobinson, &["C1"], C::PengRobinson, "DEFAULT").unwrap();
    assert_eq!(explicit.alpha(0).unwrap().parameters[0], 0.);
    assert_ne!(gerg.alpha(0).unwrap().parameters[0], 0.);
    for eos in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let umr = Cubic::with_alpha(eos, &["C1"], C::Umr, "DEFAULT").unwrap();
        assert!(umr.prepare(300., &[1.]).is_err());
        let mut explicit = Cubic::new(eos, &["C1"]).unwrap();
        explicit.set_alpha(0, C::Gerg, &[0.2, 0.9, 2.1]).unwrap();
        assert!(explicit.prepare(300., &[1.]).is_err());
    }
    let mut model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let before = model.alpha(0).unwrap().parameters;
    assert!(model.set_alpha(0, C::Twu, &[1.]).is_err());
    assert!(model.set_alpha(0, C::PengRobinson, &[f64::NAN]).is_err());
    assert!(model.set_alpha(1, C::PengRobinson, &[0.7]).is_err());
    assert_eq!(model.alpha(0).unwrap().parameters, before);
    assert!(
        Alpha::new(C::Classic, &[0.7])
            .unwrap()
            .evaluate(-1., 300.)
            .is_err()
    );
}
