mod support;
use thermopack_rs::{Cubic, CubicKind, FlashPhase, FlashResult};

#[test]
fn isochoric_and_tv_caloric_properties_match_the_exact_source() {
    for (fixture, expected_count, tv) in [
        (include_str!("fixtures/caloric_tv.csv"), 141, true),
        (include_str!("fixtures/isochoric_flash.csv"), 409, false),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let model = support::fixture_model(&fields);
            let value = fields[2].parse().unwrap();
            let volume = fields[3].parse().unwrap();
            let z: Vec<_> = fields[4].split(',').map(|v| v.parse().unwrap()).collect();
            let actual = if tv {
                let p = model
                    .prepare(value, &z)
                    .unwrap()
                    .caloric_tv(volume)
                    .unwrap();
                vec![p.pressure, p.internal_energy, p.entropy]
            } else {
                let config: Vec<_> = fields[5].split('@').collect();
                let mut options = model.volume_flash_options();
                if let Some(guess) = config.iter().find(|op| op.starts_with("guess:")) {
                    let values: Vec<_> = guess.split(':').collect();
                    options.temperature_guess = Some(values[1].parse().unwrap());
                    options.pressure_guess = Some(values[2].parse().unwrap());
                }
                let f = if config[0] == "V_UV" {
                    model.uv_flash_with_options(&z, value, volume, options)
                } else {
                    model.sv_flash_with_options(&z, value, volume, options)
                }
                .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
                let props = model.equilibrium_properties(&f).unwrap();
                let returned = if config[0] == "V_UV" {
                    props.internal_energy
                } else {
                    props.entropy
                };
                let scale = thermopack_rs::R
                    * if config[0] == "V_UV" {
                        f.temperature
                    } else {
                        1.
                    };
                assert!((returned - value).abs() < 1e-10 * scale);
                assert!((props.volume - volume).abs() < 1e-10 * volume);
                let two = f.phase == FlashPhase::TwoPhase;
                let mut values = vec![
                    f.temperature,
                    f.pressure,
                    f64::from(two),
                    if two { f.vapor_fraction } else { 0. },
                    if two { f.liquid_fraction } else { 1. },
                ];
                if two {
                    assert!(f.fugacity_residual <= options.scaled_tolerance);
                    for (i, &zi) in z.iter().enumerate() {
                        assert!(
                            (f.liquid_fraction * f.liquid[i] + f.vapor_fraction * f.vapor[i] - zi)
                                .abs()
                                < 1e-14
                        );
                    }
                    values.extend(f.liquid);
                    values.extend(f.vapor);
                } else {
                    values.extend(z.clone());
                    values.extend(vec![0.; z.len()]);
                }
                values
            };
            assert_eq!(actual.len(), columns.len() - 1);
            for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                let tolerance = if tv {
                    [1e-4, 1e-6, 1e-8][i] + 2e-9 * a.abs().max(b.abs())
                } else if i == 2 {
                    0.
                } else {
                    (match i {
                        0 => 2e-5,
                        1 => 0.5,
                        _ => 2e-7,
                    }) + 2e-7 * a.abs().max(b.abs())
                };
                assert!(
                    (a - b).abs() <= tolerance,
                    "{} field {i}: {a}, {b}",
                    columns[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, expected_count);
    }
}

fn kinds() -> [CubicKind; 4] {
    [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ]
}

#[test]
fn direct_tv_caloric_properties_match_tp_and_thermodynamic_derivatives() {
    for kind in kinds() {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.5, 0.3, 0.2];
        for t in [150., 300., 700.] {
            let state = model.prepare(t, &z).unwrap();
            for volume in [8e-5, 1e-4, 1e-3] {
                let props = state.caloric_tv(volume).unwrap();
                let pressure = state.pressure(volume).unwrap();
                assert!((props.pressure - pressure.value).abs() < 1e-8);
                if pressure.value > 0. {
                    assert!(
                        (props.internal_energy
                            - (state.enthalpy(volume).unwrap() - pressure.value * volume))
                            .abs()
                            < 1e-7
                    );
                    assert!((props.entropy - state.entropy(volume).unwrap()).abs() < 1e-9);
                }
                let dt = 1e-3;
                let dv = volume * 1e-5;
                let lo = model
                    .prepare(t - dt, &z)
                    .unwrap()
                    .caloric_tv(volume)
                    .unwrap();
                let hi = model
                    .prepare(t + dt, &z)
                    .unwrap()
                    .caloric_tv(volume)
                    .unwrap();
                let cv = (hi.internal_energy - lo.internal_energy) / (2. * dt);
                assert!((cv - props.heat_capacity_v).abs() < 1e-5);
                let lo = state.caloric_tv(volume - dv).unwrap();
                let hi = state.caloric_tv(volume + dv).unwrap();
                let sv = (hi.entropy - lo.entropy) / (2. * dv);
                assert!((sv - pressure.dt).abs() < 1e-4 + 3e-7 * pressure.dt.abs());
            }
        }
    }
}

fn latent(model: &Cubic, t: f64, bv: f64) -> FlashResult {
    let sat = model.saturation_pressure(t).unwrap();
    FlashResult {
        temperature: t,
        pressure: sat.pressure,
        feed: vec![1.],
        phase: FlashPhase::TwoPhase,
        liquid: vec![1.],
        vapor: vec![1.],
        vapor_fraction: bv,
        liquid_fraction: 1. - bv,
        liquid_volume: Some(sat.liquid_volume),
        vapor_volume: Some(sat.vapor_volume),
        iterations: 0,
        fugacity_residual: 0.,
    }
}

#[test]
fn uv_sv_recover_mixture_and_pure_latent_states() {
    for kind in kinds() {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.5, 0.3, 0.2];
        for (t, p) in [(200., 1e6), (350., 1e6), (200., 8e6)] {
            let original = model.tp_flash(t, p, &z).unwrap();
            let properties = model.equilibrium_properties(&original).unwrap();
            for entropy in [false, true] {
                let recovered = if entropy {
                    model.sv_flash(&z, properties.entropy, properties.volume)
                } else {
                    model.uv_flash(&z, properties.internal_energy, properties.volume)
                }
                .unwrap_or_else(|e| panic!("{kind:?} T={t} P={p} entropy={entropy}: {e}"));
                assert!(
                    (recovered.temperature - t).abs() < 1e-6,
                    "{kind:?} {t}: {}",
                    recovered.temperature
                );
                assert!(
                    (recovered.pressure - p).abs() < 0.1,
                    "{kind:?} {p}: {}",
                    recovered.pressure
                );
            }
        }
        let pure = Cubic::new(kind, &["C1"]).unwrap();
        for quality in [0.1, 0.5, 0.9] {
            let original = latent(&pure, 150., quality);
            let properties = pure.equilibrium_properties(&original).unwrap();
            for entropy in [false, true] {
                let recovered = if entropy {
                    pure.sv_flash(&[1.], properties.entropy, properties.volume)
                } else {
                    pure.uv_flash(&[1.], properties.internal_energy, properties.volume)
                }
                .unwrap_or_else(|e| panic!("{kind:?} quality={quality} entropy={entropy}: {e}"));
                assert_eq!(recovered.phase, FlashPhase::TwoPhase);
                assert!((recovered.temperature - 150.).abs() < 1e-7);
                assert!((recovered.vapor_fraction - quality).abs() < 1e-8);
            }
        }
    }
}

#[test]
fn invalid_specifications_and_excluded_states_return_errors() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let z = [0.5, 0.3, 0.2];
    let props = model
        .equilibrium_properties(&model.tp_flash(200., 1e6, &z).unwrap())
        .unwrap();
    for v in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(model.uv_flash(&z, props.internal_energy, v).is_err());
    }
    assert!(model.sv_flash(&z, f64::NAN, props.volume).is_err());
    assert!(
        model
            .uv_flash(&[1., 1., 1.], props.internal_energy, props.volume)
            .is_err()
    );
    let mut options = model.volume_flash_options();
    options.minimum_temperature = 250.;
    options.maximum_temperature = 260.;
    assert!(
        model
            .uv_flash_with_options(&z, props.internal_energy, props.volume, options)
            .is_err()
    );
    options = model.volume_flash_options();
    options.pressure_guess = Some(f64::INFINITY);
    assert!(
        model
            .sv_flash_with_options(&z, props.entropy, props.volume, options)
            .is_err()
    );
    options = model.volume_flash_options();
    options.max_iterations = 0;
    assert!(
        model
            .uv_flash_with_options(&z, props.internal_energy, props.volume, options)
            .is_err()
    );
}

#[test]
fn tv_caloric_derivatives_include_translation_and_covolume_temperature_terms() {
    for kind in kinds() {
        let mut model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        for i in 0..3 {
            model
                .set_beta(
                    i,
                    thermopack_rs::Beta::new(
                        thermopack_rs::BetaCorrelation::Quantum,
                        [3. + i as f64, 12.],
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        if kind != CubicKind::VanDerWaals {
            model.init_peneloux_volume_translation("Default").unwrap();
            for i in 0..3 {
                model
                    .set_volume_translation(
                        i,
                        thermopack_rs::VolumeTranslation::new(
                            thermopack_rs::TranslationKind::Quintic,
                            [2e-6, 1e-8, -1e-11, 1e-14, -1e-17, 1e-20],
                        )
                        .unwrap(),
                    )
                    .unwrap();
            }
        }
        let z = [0.5, 0.3, 0.2];
        for t in [150., 300., 700.] {
            let v = 1e-3;
            let dt = 1e-3;
            let dv = 1e-8;
            let state = model.prepare(t, &z).unwrap();
            let value = state.caloric_tv(v).unwrap();
            let lo = model.prepare(t - dt, &z).unwrap().caloric_tv(v).unwrap();
            let hi = model.prepare(t + dt, &z).unwrap().caloric_tv(v).unwrap();
            assert!(
                ((hi.internal_energy - lo.internal_energy) / (2. * dt) - value.heat_capacity_v)
                    .abs()
                    < 1e-5
            );
            assert!(
                ((hi.entropy - lo.entropy) / (2. * dt) - value.heat_capacity_v / t).abs() < 1e-7
            );
            let lo = state.caloric_tv(v - dv).unwrap();
            let hi = state.caloric_tv(v + dv).unwrap();
            assert!(
                ((hi.entropy - lo.entropy) / (2. * dv) - state.pressure(v).unwrap().dt).abs()
                    < 1e-4
            );
        }
    }
}
