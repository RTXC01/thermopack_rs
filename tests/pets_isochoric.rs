use thermopack_rs::{FlashPhase, Pets, R, VolumeFlashOptions};

fn options(t: f64) -> VolumeFlashOptions {
    VolumeFlashOptions {
        temperature_guess: Some(t * 1.2),
        minimum_temperature: 1.,
        maximum_temperature: 2000.,
        minimum_pressure: 1e-12,
        maximum_pressure: 1e10,
        ..Default::default()
    }
}

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn uv_sv_match_qualified_exact_source_states() {
    for (fixture, expected) in [
        (include_str!("fixtures/pets_uv.csv"), 224),
        (include_str!("fixtures/pets_sv.csv"), 175),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let params: Vec<f64> = fields[1].split(',').map(|v| v.parse().unwrap()).collect();
            let model = Pets::with_parameters(params[0], params[1]).unwrap();
            let target = fields[2].parse().unwrap();
            let volume = fields[3].parse().unwrap();
            let config: Vec<_> = fields[5].split('@').collect();
            let options = VolumeFlashOptions {
                temperature_guess: if config[1] == "none" {
                    None
                } else {
                    Some(config[1].parse().unwrap())
                },
                pressure_guess: if config[2] == "none" {
                    None
                } else {
                    Some(config[2].parse().unwrap())
                },
                minimum_temperature: 1.,
                maximum_temperature: 2000.,
                minimum_pressure: 1e-12,
                maximum_pressure: 1e10,
                ..Default::default()
            };
            let flash = if config[0] == "PETS_UV" {
                model.uv_flash_with_options(target, volume, options)
            } else {
                model.sv_flash_with_options(target, volume, options)
            }
            .unwrap();
            let state = model.equilibrium_tv(flash.temperature, volume).unwrap();
            let two = flash.phase == FlashPhase::TwoPhase
                && flash.vapor_fraction > 1e-8
                && flash.vapor_fraction < 1. - 1e-8;
            let p = state.properties;
            let actual = [
                flash.temperature,
                flash.pressure,
                f64::from(two),
                if two { flash.vapor_fraction } else { 0. },
                p.volume,
                p.internal_energy,
                p.entropy,
                p.enthalpy,
            ];
            let atol = [2e-5, 0.5, 0., 2e-7, 1e-9, 1e-3, 1e-5, 1e-3];
            for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                assert!(
                    (a - b).abs() <= atol[i] + 2e-7 * a.abs().max(b.abs()),
                    "{} field {i}: {a} != {b}",
                    columns[0]
                );
            }
            let residual = if config[0] == "PETS_UV" {
                (p.internal_energy - target) / (R * flash.temperature)
            } else {
                (p.entropy - target) / R
            };
            assert!(
                residual.abs() < 2e-12,
                "{}: residual={residual}",
                columns[0]
            );
            count += 1;
        }
        assert_eq!(count, expected);
    }
}

#[test]
fn pure_volume_heat_capacity_includes_phase_redistribution() {
    for (sigma, eps) in [(3.42e-10, 136.), (2.9e-10, 100.), (4.1e-10, 240.)] {
        let model = Pets::with_parameters(sigma, eps).unwrap();
        let tc = model.critical_point().unwrap().temperature;
        for tr in [0.5, 0.7, 0.9, 0.99, 0.999] {
            let t = tr * tc;
            let sat = model.saturation_pressure(t).unwrap();
            for beta in [0.01, 0.1, 0.5, 0.9, 0.99] {
                let v = (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume;
                let eq = model.equilibrium_tv(t, v).unwrap();
                assert_eq!(eq.flash.phase, FlashPhase::TwoPhase);
                close(eq.flash.vapor_fraction, beta, 1e-10, 1e-10);
                assert!(eq.heat_capacity_v > 0. && eq.properties.heat_capacity.is_none());
                let dt = 1e-5 * t;
                let lo = model.equilibrium_tv(t - dt, v).unwrap();
                let hi = model.equilibrium_tv(t + dt, v).unwrap();
                let cv =
                    (hi.properties.internal_energy - lo.properties.internal_energy) / (2. * dt);
                let entropy_t = (hi.properties.entropy - lo.properties.entropy) / (2. * dt);
                close(eq.heat_capacity_v, cv, 1e-4, 2e-5);
                close(eq.heat_capacity_v / t, entropy_t, 1e-6, 2e-5);
                close(
                    eq.properties.enthalpy - eq.properties.internal_energy,
                    eq.flash.pressure * v,
                    1e-9,
                    1e-10,
                );
            }
        }
    }
}

#[test]
fn uv_sv_recover_homogeneous_and_latent_states() {
    for (sigma, eps) in [(3.42e-10, 136.), (2.9e-10, 100.), (4.1e-10, 240.)] {
        let model = Pets::with_parameters(sigma, eps).unwrap();
        let tc = model.critical_point().unwrap().temperature;
        let mut seeds = Vec::new();
        for tr in [0.5, 0.8, 1., 1.5, 3.] {
            for p in [1e3, 1e5, 1e6, 1e8] {
                let flash = model.tp_flash(eps * tr, p).unwrap();
                let props = model.equilibrium_properties(&flash).unwrap();
                seeds.push((flash.temperature, props.volume));
            }
        }
        for tr in [0.5, 0.7, 0.9, 0.99, 0.999] {
            let sat = model.saturation_pressure(tc * tr).unwrap();
            for beta in [0., 0.1, 0.5, 0.9, 1.] {
                seeds.push((
                    sat.temperature,
                    (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume,
                ));
            }
        }
        for (t, v) in seeds {
            let state = model.equilibrium_tv(t, v).unwrap();
            for entropy in [false, true] {
                let recovered = if entropy {
                    model.sv_flash_with_options(state.properties.entropy, v, options(t))
                } else {
                    model.uv_flash_with_options(state.properties.internal_energy, v, options(t))
                }
                .unwrap();
                close(recovered.temperature, t, 1e-7, 1e-10);
                close(recovered.pressure, state.flash.pressure, 1e-4, 1e-8);
                let result = model.equilibrium_tv(recovered.temperature, v).unwrap();
                assert!(
                    (result.properties.internal_energy - state.properties.internal_energy).abs()
                        / (R * t)
                        < 2e-11
                );
                assert!((result.properties.entropy - state.properties.entropy).abs() / R < 2e-11);
            }
        }
    }
}

#[test]
fn homogeneous_cv_reference_ownership_bounds_and_invalid_inputs() {
    let original = Pets::new().unwrap();
    let mut model = original.clone();
    model.ideal_gas_mut().set_heat_capacity(8, &[3.5]).unwrap();
    model.ideal_gas_mut().set_formation_enthalpy(1e5).unwrap();
    model
        .ideal_gas_mut()
        .set_standard_entropy(300., thermopack_rs::EntropyReference::Bar)
        .unwrap();
    assert_eq!(original.ideal_gas().formation_enthalpy(), 0.);
    let v = 0.001;
    let state = model.equilibrium_tv(300., v).unwrap();
    let dt = 1e-3;
    let lo = model.equilibrium_tv(300. - dt, v).unwrap();
    let hi = model.equilibrium_tv(300. + dt, v).unwrap();
    close(
        state.heat_capacity_v,
        (hi.properties.internal_energy - lo.properties.internal_energy) / (2. * dt),
        1e-5,
        1e-8,
    );
    for entropy in [false, true] {
        let flash = if entropy {
            model.sv_flash(state.properties.entropy, v)
        } else {
            model.uv_flash(state.properties.internal_energy, v)
        }
        .unwrap();
        close(flash.temperature, 300., 1e-7, 1e-10);
    }
    let latent = model.equilibrium_tv(100., v).unwrap();
    assert_eq!(latent.flash.phase, FlashPhase::TwoPhase);
    close(
        model
            .uv_flash(latent.properties.internal_energy, v)
            .unwrap()
            .temperature,
        100.,
        1e-7,
        1e-10,
    );
    close(
        model
            .sv_flash(latent.properties.entropy, v)
            .unwrap()
            .temperature,
        100.,
        1e-7,
        1e-10,
    );
    let c = original.critical_point().unwrap();
    let critical = original.equilibrium_tv(c.temperature, c.volume).unwrap();
    assert!(critical.heat_capacity_v.is_finite() && critical.heat_capacity_v > 0.);
    close(critical.flash.pressure, c.pressure, 1e-5, 1e-10);
    assert!(model.equilibrium_tv(100., 0.).is_err());
    assert!(model.equilibrium_tv(100., 1e-9).is_err());
    assert!(model.uv_flash(f64::NAN, v).is_err());
    assert!(model.sv_flash(1., f64::INFINITY).is_err());
    assert!(
        model
            .uv_flash_with_options(
                1.,
                v,
                VolumeFlashOptions {
                    max_iterations: 0,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        model
            .sv_flash_with_options(
                1.,
                v,
                VolumeFlashOptions {
                    pressure_guess: Some(f64::NAN),
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        model
            .uv_flash_with_options(
                state.properties.internal_energy,
                v,
                VolumeFlashOptions {
                    minimum_pressure: 1.,
                    maximum_pressure: 10.,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        model
            .sv_flash_with_options(
                state.properties.entropy,
                v,
                VolumeFlashOptions {
                    minimum_temperature: 400.,
                    maximum_temperature: 500.,
                    ..Default::default()
                }
            )
            .is_err()
    );
}

#[test]
fn nonfinite_caloric_values_and_capacity_cannot_converge() {
    for (correlation, coefficients) in [
        (7, &[0.; 5][..]),
        (4, &[1e308][..]),
        (7, &[30., 1., 1e308, 0., 1.][..]),
    ] {
        let mut model = Pets::new().unwrap();
        model
            .ideal_gas_mut()
            .set_heat_capacity(correlation, coefficients)
            .unwrap();
        assert!(model.equilibrium_tv(300., 1e-3).is_err());
        assert!(model.uv_flash(0., 1e-3).is_err());
        assert!(model.sv_flash(0., 1e-3).is_err());
    }
}
