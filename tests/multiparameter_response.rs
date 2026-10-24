use thermopack_rs::{FlashPhase, Multiparameter, Phase};

#[test]
fn homogeneous_response_matches_local_fixed_pressure_difference() {
    let eos = Multiparameter::new("C3").unwrap();
    let flash = eos.tp_flash(350., 1.0e6).unwrap();
    assert_ne!(flash.phase as i32, 0);
    let r = eos.equilibrium_response(&flash).unwrap().unwrap();
    let dt = 0.02;
    let lo = eos.tp_flash(350. - dt, 1.0e6).unwrap();
    let hi = eos.tp_flash(350. + dt, 1.0e6).unwrap();
    let a = eos.equilibrium_properties(&lo).unwrap();
    let b = eos.equilibrium_properties(&hi).unwrap();
    assert!(((b.enthalpy - a.enthalpy) / (2. * dt) - r.heat_capacity).abs() < 2e-3);
    assert!(((b.volume - a.volume) / (2. * dt) - r.volume_dt).abs() < 1e-10);
    let sound = eos.equilibrium_speed_of_sound(&flash).unwrap();
    let ds = r.volume_dp + flash.temperature * r.volume_dt * r.volume_dt / r.heat_capacity;
    let phase = if flash.phase == FlashPhase::Vapor {
        Phase::Vapor
    } else {
        Phase::Liquid
    };
    let v = eos
        .prepare(flash.temperature)
        .unwrap()
        .volume(flash.pressure, phase)
        .unwrap();
    let expected = (-v * v / (eos.molar_mass() * ds)).sqrt();
    assert!((sound - expected).abs() < 1e-9 * sound);
    assert!(sound.is_finite() && sound > 0.);
    assert!(
        eos.baseline_equilibrium_speed_of_sound(&flash)
            .unwrap()
            .is_finite()
    );
    assert!(
        eos.equilibrium_response(&eos.saturated_flash(300., 0.5).unwrap())
            .unwrap()
            .is_none()
    );
    assert!(ds < 0.);
}

#[test]
fn coexistence_equilibrium_sound_is_finite_at_quality_endpoints_and_interior() {
    let eos = Multiparameter::new("R134A").unwrap();
    for beta in [0., 0.3, 1.] {
        let flash = eos.saturated_flash(300., beta).unwrap();
        let sound = eos.equilibrium_speed_of_sound(&flash).unwrap();
        assert!(
            sound.is_finite() && sound > 0.,
            "beta={beta}, sound={sound}"
        );
    }
}

#[test]
fn coexistence_sound_matches_isentropic_saturation_difference() {
    let eos = Multiparameter::new("R134A").unwrap();
    let t = 300.;
    let beta = 0.3;
    let flash = eos.saturated_flash(t, beta).unwrap();
    let sound = eos.equilibrium_speed_of_sound(&flash).unwrap();
    let sat = eos.saturation_pressure(t).unwrap();
    let prep = eos.prepare(t).unwrap();
    let sl = prep.entropy(sat.liquid_volume).unwrap();
    let sv = prep.entropy(sat.vapor_volume).unwrap();
    let target = sl + beta * (sv - sl);
    let dt = 0.01;
    let mix_v = |tt: f64| {
        let s = eos.saturation_pressure(tt).unwrap();
        let p = eos.prepare(tt).unwrap();
        let sl = p.entropy(s.liquid_volume).unwrap();
        let sv = p.entropy(s.vapor_volume).unwrap();
        let q = (target - sl) / (sv - sl);
        (s.pressure, (1. - q) * s.liquid_volume + q * s.vapor_volume)
    };
    let (pl, vl) = mix_v(t - dt);
    let (ph, vh) = mix_v(t + dt);
    let dvdp = (vh - vl) / (ph - pl);
    let v = (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume;
    let expected = (-v * v / (eos.molar_mass() * dvdp)).sqrt();
    assert!(
        (sound - expected).abs() / sound < 2e-5,
        "sound={sound}, fd={expected}"
    );
}

#[test]
fn malformed_non_pure_flash_is_rejected() {
    let eos = Multiparameter::new("C3").unwrap();
    let mut flash = eos.tp_flash(350., 1e6).unwrap();
    flash.feed[0] = 0.9;
    assert!(eos.equilibrium_response(&flash).is_err());
    assert!(eos.equilibrium_speed_of_sound(&flash).is_err());

    let mut split = eos.saturated_flash(300., 0.5).unwrap();
    split.liquid_fraction = 1.2;
    split.vapor_fraction = -0.2;
    assert!(eos.equilibrium_response(&split).is_err());
    assert!(eos.equilibrium_speed_of_sound(&split).is_err());
}

#[test]
fn response_recomputes_volumes_and_rejects_inconsistent_coexistence() {
    let eos = Multiparameter::new("C3").unwrap();
    for mut flash in [
        eos.tp_flash(350., 1e6).unwrap(),
        eos.saturated_flash(300., 0.3).unwrap(),
    ] {
        let expected = eos.equilibrium_speed_of_sound(&flash).unwrap();
        let baseline = eos.baseline_equilibrium_speed_of_sound(&flash).unwrap();
        flash.liquid_volume = Some(f64::NAN);
        flash.vapor_volume = Some(-1.);
        assert_eq!(eos.equilibrium_speed_of_sound(&flash).unwrap(), expected);
        assert_eq!(
            eos.baseline_equilibrium_speed_of_sound(&flash).unwrap(),
            baseline
        );
    }
    let split = eos.saturated_flash(300., 0.3).unwrap();
    for pressure in [f64::NAN, f64::INFINITY, 0., split.pressure * 1.01] {
        let mut bad = split.clone();
        bad.pressure = pressure;
        assert!(eos.equilibrium_response(&bad).is_err());
        assert!(eos.equilibrium_speed_of_sound(&bad).is_err());
        assert!(eos.baseline_equilibrium_speed_of_sound(&bad).is_err());
    }
    for beta in [f64::NAN, f64::INFINITY, -1e-10, 1. + 1e-10] {
        let mut bad = split.clone();
        bad.vapor_fraction = beta;
        bad.liquid_fraction = 1. - beta;
        assert!(eos.equilibrium_speed_of_sound(&bad).is_err());
        assert!(eos.baseline_equilibrium_speed_of_sound(&bad).is_err());
    }
}

#[test]
fn source_tp_sound_uses_the_component_database_mass() {
    let eos = Multiparameter::lennard_jones("AR", thermopack_rs::LennardJonesKind::Full).unwrap();
    let flash = eos.tp_flash(150., 1.0e5).unwrap();
    let physical = eos.equilibrium_speed_of_sound(&flash).unwrap();
    let baseline = eos.baseline_equilibrium_speed_of_sound(&flash).unwrap();
    let database_mass = thermopack_rs::component("AR").unwrap().molar_mass * 1e-3;
    let expected = (eos.molar_mass() / database_mass).sqrt();
    assert!((baseline / physical - expected).abs() < 1e-12);
}

#[test]
fn response_and_sound_cover_all_native_multiparameter_variants() {
    let mut models = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
        .into_iter()
        .map(|name| Multiparameter::new(name).unwrap())
        .collect::<Vec<_>>();
    for gas in ["NE", "AR", "KR", "XE"] {
        for kind in [
            thermopack_rs::LennardJonesKind::Full,
            thermopack_rs::LennardJonesKind::TruncatedShifted,
        ] {
            models.push(Multiparameter::lennard_jones(gas, kind).unwrap());
        }
    }
    models.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
    assert_eq!(models.len(), 14);
    for eos in models {
        let t = 1.4 * eos.tc();
        let v = 1. / eos.rc();
        let p = eos.prepare(t).unwrap().pressure(v).unwrap().value;
        assert!(p > 0.);
        let flash = eos.tp_flash(t, p).unwrap();
        assert!(
            eos.equilibrium_response(&flash).unwrap().is_some(),
            "{}",
            eos.record().ident
        );
        let c = eos.equilibrium_speed_of_sound(&flash).unwrap();
        assert!(c.is_finite() && c > 0., "{}", eos.record().ident);
    }
}

fn fixture_model(config: &str) -> Multiparameter {
    let p: Vec<_> = config.split(':').collect();
    match p[0] {
        "LJ_MEOS" => {
            Multiparameter::lennard_jones(p[1], thermopack_rs::LennardJonesKind::Full).unwrap()
        }
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(p[1], thermopack_rs::LennardJonesKind::TruncatedShifted)
                .unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => Multiparameter::new(p[1]).unwrap(),
    }
}

fn replay_source(fixture: &str, count: usize) {
    let mut rows = 0;
    for row in fixture.lines().skip(1) {
        let f: Vec<_> = row.split(';').collect();
        assert_eq!(f.len(), 8);
        let case: Vec<_> = f[0].split('|').collect();
        let model = fixture_model(case[1]);
        let t = case[2].parse().unwrap();
        let p = case[3].parse().unwrap();
        let mut flash = if case[4].starts_with("HOM") {
            let mut v = model.tp_flash(t, p).unwrap();
            v.phase = match case[4] {
                "HOM:L" => FlashPhase::Liquid,
                "HOM:V" => FlashPhase::Vapor,
                _ => FlashPhase::Single,
            };
            v
        } else {
            model.saturated_flash(t, case[4].parse().unwrap()).unwrap()
        };
        flash.pressure = p;
        let response = model.equilibrium_response(&flash).unwrap();
        let (cp, vt, vp) = response.map_or((0., 0., 0.), |r| {
            (r.heat_capacity, r.volume_dt, r.volume_dp)
        });
        let actual = [
            f64::from(response.is_some()),
            cp,
            vt,
            vp,
            model.equilibrium_speed_of_sound(&flash).unwrap(),
            model.baseline_equilibrium_speed_of_sound(&flash).unwrap(),
            flash.pressure,
        ];
        let floors = [0., 1e-8, 1e-12, 1e-16, 1e-7, 1e-6, 1e-6];
        for (i, value) in actual.into_iter().enumerate() {
            let expected: f64 = f[i + 1].parse().unwrap();
            let relative = [0., 2e-7, 2e-7, 2e-7, 2e-8, 2e-6, 2e-8][i];
            assert!(
                (value - expected).abs() <= floors[i] + relative * value.abs().max(expected.abs()),
                "{} field {i}: {value} != {expected}",
                f[0]
            );
        }
        rows += 1;
    }
    assert_eq!(rows, count);
}

#[test]
fn original_source_expected_values_replay_without_native_gold() {
    replay_source(
        include_str!("fixtures/multiparameter_response_source.csv"),
        854,
    );
}

#[test]
fn isentropic_curve_differences_cover_all_variants_and_phase_fractions() {
    let mut models: Vec<_> = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
        .iter()
        .map(|s| Multiparameter::new(s).unwrap())
        .collect();
    for kind in [
        thermopack_rs::LennardJonesKind::Full,
        thermopack_rs::LennardJonesKind::TruncatedShifted,
    ] {
        for c in ["NE", "AR", "KR", "XE"] {
            models.push(Multiparameter::lennard_jones(c, kind).unwrap());
        }
    }
    models.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
    for model in models {
        for tr in [0.5, 0.8, 0.95, 0.99, 0.999] {
            let t = tr * model.tc();
            let sat = model.saturation_pressure(t).unwrap();
            let prep = model.prepare(t).unwrap();
            let sl = prep.entropy(sat.liquid_volume).unwrap();
            let sv = prep.entropy(sat.vapor_volume).unwrap();
            for beta in [0.1, 0.5, 0.9] {
                let target = (1. - beta) * sl + beta * sv;
                let volume = (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume;
                let h = (3e-6 * t).min(0.001 * (model.tc() - t));
                let point = |tt: f64| {
                    let s = model.saturation_pressure(tt).unwrap();
                    let p = model.prepare(tt).unwrap();
                    let l = p.entropy(s.liquid_volume).unwrap();
                    let v = p.entropy(s.vapor_volume).unwrap();
                    let b = (target - l) / (v - l);
                    assert!((0. ..=1.).contains(&b));
                    (s.pressure, (1. - b) * s.liquid_volume + b * s.vapor_volume)
                };
                let (pl, vl) = point(t - h);
                let (pu, vu) = point(t + h);
                let derivative = (vu - vl) / (pu - pl);
                assert!(derivative < 0.);
                let fd = (-volume * volume / (model.molar_mass() * derivative)).sqrt();
                let sound = model
                    .equilibrium_speed_of_sound(&model.saturated_flash(t, beta).unwrap())
                    .unwrap();
                assert!(
                    (fd / sound - 1.).abs() < 3e-4,
                    "{} tr{tr} beta{beta}: {fd} != {sound}",
                    model.record().ident
                );
            }
        }
    }
}
