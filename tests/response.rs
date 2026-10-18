mod support;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, FlashPhase};

#[test]
fn reported_equilibrium_derivatives_match_source_helpers() {
    let absolute = [0., 1e-6, 1e-10, 1e-14];
    let mut count = 0;
    for row in include_str!("fixtures/equilibrium_response.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = support::fixture_model(&fields);
        let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
        let flash = model
            .baseline_tp_flash(fields[2].parse().unwrap(), fields[3].parse().unwrap(), &z)
            .unwrap();
        let response = model
            .baseline_equilibrium_response(&flash)
            .unwrap()
            .unwrap();
        let actual = [
            flash.phase as i32 as f64,
            response.heat_capacity,
            response.volume_dt,
            response.volume_dp,
        ];
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = absolute[i]
                + if i == 0 {
                    0.
                } else {
                    2e-7 * a.abs().max(b.abs())
                };
            assert!(
                (a - b).abs() <= tolerance,
                "{} field {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 435);
}

#[test]
fn equilibrium_response_differentiates_tp_caloric_and_volume_values() {
    let z = [0.5, 0.3, 0.2];
    let mut count = 0;
    let mut two_phase = 0;
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for translated in [false, true] {
            let mut model = Cubic::with_options(
                kind,
                &["C1", "C2", "C3"],
                CubicOptions {
                    volume_translation: translated && kind != CubicKind::VanDerWaals,
                    ..Default::default()
                },
            )
            .unwrap();
            model.set_kij(0, 1, 0.04).unwrap();
            model.set_lij(1, 2, 0.02).unwrap();
            for t in [180., 250., 350.] {
                for p in [1e5, 1e6, 8e6] {
                    let flash = model.tp_flash(t, p, &z).unwrap();
                    let property = model.equilibrium_properties(&flash).unwrap();
                    let response = model.equilibrium_response(&flash).unwrap().unwrap();
                    let at = |t, p| {
                        model
                            .equilibrium_properties(&model.tp_flash(t, p, &z).unwrap())
                            .unwrap()
                    };
                    let dt = 1e-3;
                    let dp = p * 1e-4;
                    let lo_t = at(t - dt, p);
                    let hi_t = at(t + dt, p);
                    let lo_p = at(t, p - dp);
                    let hi_p = at(t, p + dp);
                    let check = |a: f64, b: f64, absolute: f64| {
                        assert!(
                            (a - b).abs() <= absolute + 3e-5 * a.abs().max(b.abs()),
                            "{kind:?} T={t} P={p}: {a}, {b}"
                        )
                    };
                    check(
                        response.heat_capacity,
                        (hi_t.enthalpy - lo_t.enthalpy) / (2. * dt),
                        1e-6,
                    );
                    check(
                        response.heat_capacity,
                        property.heat_capacity.unwrap(),
                        1e-9,
                    );
                    check(
                        response.volume_dt,
                        (hi_t.volume - lo_t.volume) / (2. * dt),
                        1e-11,
                    );
                    check(
                        response.volume_dp,
                        (hi_p.volume - lo_p.volume) / (2. * dp),
                        1e-15,
                    );
                    check(
                        response.internal_energy_dt(p),
                        (hi_t.internal_energy - lo_t.internal_energy) / (2. * dt),
                        1e-6,
                    );
                    check(
                        response.internal_energy_dp(t, p),
                        (hi_p.internal_energy - lo_p.internal_energy) / (2. * dp),
                        1e-8,
                    );
                    check(
                        -response.volume_dt,
                        (hi_p.entropy - lo_p.entropy) / (2. * dp),
                        1e-10,
                    );
                    assert!(response.isentropic_volume_dp(t) < 0.);
                    assert!(response.heat_capacity_v(t) > 0.);
                    if flash.phase == FlashPhase::TwoPhase {
                        two_phase += 1;
                    }
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 72);
    assert!(two_phase > 10);
}

#[test]
fn pure_latent_heat_has_no_independent_tp_response_and_zero_feed_stays_absent() {
    let pure = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let sat = pure.saturation_pressure(150.).unwrap();
    let l = pure
        .prepare(150., &[1.])
        .unwrap()
        .enthalpy(sat.liquid_volume)
        .unwrap();
    let v = pure
        .prepare(150., &[1.])
        .unwrap()
        .enthalpy(sat.vapor_volume)
        .unwrap();
    let flash = pure.ph_flash(sat.pressure, &[1.], 0.5 * (l + v)).unwrap();
    assert_eq!(flash.phase, FlashPhase::TwoPhase);
    assert!(pure.equilibrium_response(&flash).unwrap().is_none());
    let model = Cubic::new(CubicKind::PengRobinson, &["N2", "C1", "C2", "C3"]).unwrap();
    let z = [0., 0.5, 0.3, 0.2];
    let flash = model.tp_flash(200., 1e6, &z).unwrap();
    let response = model.equilibrium_response(&flash).unwrap().unwrap();
    assert!(response.heat_capacity > 0. && response.volume_dp < 0.);
}
