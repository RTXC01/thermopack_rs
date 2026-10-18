use thermopack_rs::{CaloricFlashOptions, FlashPhase, Pets, R};

fn options(t: f64) -> CaloricFlashOptions {
    CaloricFlashOptions {
        temperature_guess: Some(t * 1.2),
        minimum_temperature: 1.,
        maximum_temperature: 2000.,
        ..Default::default()
    }
}

#[test]
fn pure_flashes_match_qualified_source_results() {
    for (fixture, expected) in [
        (include_str!("fixtures/pets_tp.csv"), 126),
        (include_str!("fixtures/pets_ph.csv"), 395),
        (include_str!("fixtures/pets_ps.csv"), 402),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let params: Vec<f64> = fields[1].split(',').map(|v| v.parse().unwrap()).collect();
            let model = Pets::with_parameters(params[0], params[1]).unwrap();
            let t = fields[2].parse().unwrap();
            let p = fields[3].parse().unwrap();
            let config: Vec<_> = fields[5].split('@').collect();
            let options = CaloricFlashOptions {
                temperature_guess: if config[2] == "none" {
                    None
                } else {
                    Some(config[2].parse().unwrap())
                },
                minimum_temperature: 1.,
                maximum_temperature: 2000.,
                ..Default::default()
            };
            let flash = match config[0] {
                "PETS_TP" => model.tp_flash(t, p),
                "PETS_PH" => model.ph_flash_with_options(p, config[1].parse().unwrap(), options),
                "PETS_PS" => model.ps_flash_with_options(p, config[1].parse().unwrap(), options),
                _ => panic!("unknown fixture mode"),
            }
            .unwrap();
            let state = model.equilibrium_properties(&flash).unwrap();
            let two = flash.phase == FlashPhase::TwoPhase
                && flash.vapor_fraction > 1e-8
                && flash.vapor_fraction < 1. - 1e-8;
            let actual = [
                flash.temperature,
                p,
                f64::from(two),
                if two { flash.vapor_fraction } else { 0. },
                if two {
                    flash.liquid_volume.unwrap()
                } else {
                    state.volume
                },
                if two {
                    flash.vapor_volume.unwrap()
                } else {
                    state.volume
                },
                state.enthalpy,
                state.entropy,
            ];
            let atol = [1e-5, 0.01, 0., 1e-7, 1e-9, 1e-9, 1e-3, 1e-5];
            for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                assert!(
                    (a - b).abs() <= atol[i] + 2e-7 * a.abs().max(b.abs()),
                    "{} field {i}: {a} != {b}",
                    columns[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, expected);
    }
}

#[test]
fn homogeneous_ph_ps_recover_temperature_and_close_specifications() {
    for (sigma, epsilon) in [(3.42e-10, 136.), (2.9e-10, 100.), (4.1e-10, 240.)] {
        let model = Pets::with_parameters(sigma, epsilon).unwrap();
        for tr in [0.5, 0.8, 1., 1.1, 1.5, 3.] {
            let t = tr * epsilon;
            for p in [1., 1e3, 1e5, 1e6, 5e6, 1e8] {
                let seed = model.tp_flash(t, p).unwrap();
                let state = model.equilibrium_properties(&seed).unwrap();
                for entropy in [false, true] {
                    let recovered = if entropy {
                        model.ps_flash_with_options(p, state.entropy, options(t))
                    } else {
                        model.ph_flash_with_options(p, state.enthalpy, options(t))
                    }
                    .unwrap();
                    assert!(
                        (recovered.temperature - t).abs() < 1e-7,
                        "{t}, {p}, entropy={entropy}: {}",
                        recovered.temperature
                    );
                    let result = model.equilibrium_properties(&recovered).unwrap();
                    assert!((result.enthalpy - state.enthalpy).abs() / (R * t) < 1e-9);
                    assert!((result.entropy - state.entropy).abs() / R < 1e-9);
                    assert_eq!(recovered.liquid, [1.]);
                    assert_eq!(recovered.vapor, [1.]);
                }
            }
        }
    }
}

#[test]
fn latent_quality_closes_caloric_and_material_balances() {
    let model = Pets::new().unwrap();
    for t in [80., 100., 130., 147.] {
        for beta in [0., 0.1, 0.5, 0.9, 1.] {
            let seed = model.saturated_flash(t, beta).unwrap();
            let state = model.equilibrium_properties(&seed).unwrap();
            assert!(state.heat_capacity.is_none());
            assert!(
                (state.enthalpy - state.internal_energy - seed.pressure * state.volume).abs()
                    < 1e-10
            );
            for entropy in [false, true] {
                let recovered = if entropy {
                    model.ps_flash_with_options(seed.pressure, state.entropy, options(t))
                } else {
                    model.ph_flash_with_options(seed.pressure, state.enthalpy, options(t))
                }
                .unwrap();
                assert!((recovered.temperature - t).abs() < 1e-7);
                if beta > 0. && beta < 1. {
                    assert_eq!(recovered.phase, FlashPhase::TwoPhase);
                    assert!((recovered.vapor_fraction - beta).abs() < 1e-8);
                    assert!(
                        (recovered.liquid_fraction + recovered.vapor_fraction - 1.).abs() < 1e-15
                    );
                }
            }
        }
    }
}

#[test]
fn invalid_specifications_and_off_curve_splits_are_rejected() {
    let model = Pets::new().unwrap();
    let mut flash = model.saturated_flash(100., 0.5).unwrap();
    flash.pressure *= 1.01;
    assert!(model.equilibrium_properties(&flash).is_err());
    assert!(model.saturated_flash(100., f64::NAN).is_err());
    assert!(model.saturated_flash(100., -0.1).is_err());
    assert!(model.ph_flash(1e5, f64::NAN).is_err());
    assert!(model.ps_flash(0., 1.).is_err());
    assert!(
        model
            .ph_flash_with_options(
                1e5,
                1.,
                CaloricFlashOptions {
                    max_iterations: 0,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        model
            .ps_flash_with_options(
                1e5,
                1.,
                CaloricFlashOptions {
                    temperature_guess: Some(f64::INFINITY),
                    ..Default::default()
                }
            )
            .is_err()
    );
    let seed = model.tp_flash(300., 1e5).unwrap();
    let state = model.equilibrium_properties(&seed).unwrap();
    assert!(
        model
            .ph_flash_with_options(
                1e5,
                state.enthalpy,
                CaloricFlashOptions {
                    minimum_temperature: 400.,
                    maximum_temperature: 500.,
                    ..Default::default()
                }
            )
            .is_err()
    );
    let mut malformed = seed;
    malformed.feed = vec![0.5, 0.5];
    assert!(model.equilibrium_properties(&malformed).is_err());
}

#[test]
fn caloric_setters_are_owned_and_used_by_homogeneous_and_latent_flashes() {
    let default = Pets::new().unwrap();
    let mut changed = default.clone();
    changed
        .ideal_gas_mut()
        .set_heat_capacity(8, &[3.5])
        .unwrap();
    changed.ideal_gas_mut().set_formation_enthalpy(1e5).unwrap();
    changed
        .ideal_gas_mut()
        .set_standard_entropy(300., thermopack_rs::EntropyReference::Atmosphere)
        .unwrap();
    assert_eq!(default.ideal_gas().formation_enthalpy(), 0.);
    for latent in [false, true] {
        let seed = if latent {
            changed.saturated_flash(100., 0.4)
        } else {
            changed.tp_flash(300., 1e5)
        }
        .unwrap();
        let properties = changed.equilibrium_properties(&seed).unwrap();
        for entropy in [false, true] {
            let flash = if entropy {
                changed.ps_flash(seed.pressure, properties.entropy)
            } else {
                changed.ph_flash(seed.pressure, properties.enthalpy)
            }
            .unwrap();
            assert!((flash.temperature - seed.temperature).abs() < 1e-7);
            if latent {
                assert!((flash.vapor_fraction - 0.4).abs() < 1e-8);
            }
        }
    }
}
