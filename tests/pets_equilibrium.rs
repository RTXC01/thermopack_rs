use thermopack_rs::{Pets, Phase, R};

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn equilibrium_paths_match_exact_source_fixtures() {
    for (fixture, expected, atol) in [
        (
            include_str!("fixtures/pets_volume_baseline.csv"),
            294,
            &[1e-9, 1e-3, 1e-3, 1.][..],
        ),
        (
            include_str!("fixtures/pets_volume.csv"),
            237,
            &[1e-10, 1e-10, 1e-14][..],
        ),
        (
            include_str!("fixtures/pets_critical.csv"),
            3,
            &[1e-5, 0.1, 1e-11, 1e-7][..],
        ),
        (
            include_str!("fixtures/pets_saturation.csv"),
            18,
            &[1e-5, 0.01, 1e-10, 1e-9, 0.01][..],
        ),
        (
            include_str!("fixtures/pets_saturation_temperature.csv"),
            18,
            &[1e-5, 0.01, 1e-10, 1e-9, 0.01][..],
        ),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let params: Vec<f64> = fields[1].split(',').map(|x| x.parse().unwrap()).collect();
            let model = Pets::with_parameters(params[0], params[1]).unwrap();
            let t = fields[2].parse().unwrap();
            let p = fields[3].parse().unwrap();
            let config: Vec<_> = fields[5].split('@').collect();
            let values = match config[0] {
                "PETS_VOLUME_BASELINE" | "PETS_VOLUME" => {
                    let phase = if config[1] == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    };
                    let v = if config[0] == "PETS_VOLUME_BASELINE" {
                        model.baseline_volume_iteration(t, p, phase).unwrap().volume
                    } else {
                        model.volume(t, p, phase).unwrap()
                    };
                    let state = model.pressure_tv(t, v, 1.).unwrap();
                    if config[0] == "PETS_VOLUME_BASELINE" {
                        vec![v, state.value, state.dt, state.dv]
                    } else {
                        vec![v, -state.dt / state.dv, 1. / state.dv]
                    }
                }
                "PETS_CRITICAL" => {
                    let c = model.critical_point().unwrap();
                    vec![c.temperature, c.pressure, c.volume, c.compressibility]
                }
                "PETS_SAT" | "PETS_SAT_T" => {
                    let s = if config[0] == "PETS_SAT" {
                        model.saturation_pressure(t)
                    } else {
                        model.saturation_temperature(p)
                    }
                    .unwrap();
                    vec![
                        s.temperature,
                        s.pressure,
                        s.liquid_volume,
                        s.vapor_volume,
                        s.dpdt,
                    ]
                }
                _ => panic!("unknown fixture"),
            };
            assert_eq!(values.len(), columns.len() - 1);
            for (i, (&a, b)) in values.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                assert!(
                    (a - b).abs() <= atol[i] + 2e-7 * a.abs().max(b.abs()),
                    "{} field {i}: {a}, {b}",
                    columns[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, expected);
    }
}

#[test]
fn critical_and_saturation_satisfy_the_native_equations() {
    let model = Pets::new().unwrap();
    let critical = model.critical_point().unwrap();
    assert!((136. ..163.2).contains(&critical.temperature));
    let pressure = model
        .pressure_tv(critical.temperature, critical.volume, 1.)
        .unwrap();
    assert!((pressure.value - critical.pressure).abs() < 1e-6);
    assert!((pressure.dv * critical.volume.powi(2) / (R * critical.temperature)).abs() < 1e-10);
    let dv = critical.volume * 1e-5;
    let lo = model
        .pressure_tv(critical.temperature, critical.volume - dv, 1.)
        .unwrap();
    let hi = model
        .pressure_tv(critical.temperature, critical.volume + dv, 1.)
        .unwrap();
    let curvature = (hi.dv - lo.dv) / (2. * dv);
    assert!((curvature * critical.volume.powi(3) / (R * critical.temperature)).abs() < 1e-6);
    for tr in [0.5, 0.7, 0.9, 0.99, 0.999] {
        let t = critical.temperature * tr;
        let sat = model.saturation_pressure(t).unwrap();
        assert!(sat.liquid_volume < sat.vapor_volume && sat.dpdt > 0.);
        for v in [sat.liquid_volume, sat.vapor_volume] {
            let p = model.pressure_tv(t, v, 1.).unwrap();
            assert!((p.value - sat.pressure).abs() < 1e-5 + 1e-10 * sat.pressure);
            assert!(p.dv < 0.);
        }
        let l = model.caloric_tv(t, sat.liquid_volume, 1.).unwrap();
        let v = model.caloric_tv(t, sat.vapor_volume, 1.).unwrap();
        assert!((l.chemical_potential - v.chemical_potential).abs() < 1e-7);
        let eq = model.prepare_equilibrium(t).unwrap();
        let low = eq.volume(sat.pressure * 0.99, Phase::MinimumGibbs).unwrap();
        let high = eq.volume(sat.pressure * 1.01, Phase::MinimumGibbs).unwrap();
        assert!(low > sat.vapor_volume * (1. - 1e-9) && high < sat.liquid_volume * (1. + 1e-9));
    }
}

#[test]
fn root_response_and_clapeyron_slope_agree_with_independent_differences() {
    let model = Pets::new().unwrap();
    let critical = model.critical_point().unwrap();
    for tr in [0.5, 0.7, 0.9, 0.99, 0.999] {
        let t = tr * critical.temperature;
        let sat = model.saturation_pressure(t).unwrap();
        let dt = 1e-4;
        let derivative = (model.saturation_pressure(t + dt).unwrap().pressure
            - model.saturation_pressure(t - dt).unwrap().pressure)
            / (2. * dt);
        close(derivative, sat.dpdt, 0.01, 2e-6);
        let inverse = model.saturation_temperature(sat.pressure).unwrap();
        close(inverse.temperature, t, 1e-7, 1e-10);
        for (p, phase) in [
            (sat.pressure * 0.5, Phase::Vapor),
            (sat.pressure * 2., Phase::Liquid),
        ] {
            let v = model.volume(t, p, phase).unwrap();
            let pressure = model.pressure_tv(t, v, 1.).unwrap();
            let dp = p * 1e-5;
            let vt = (model.volume(t + dt, p, phase).unwrap()
                - model.volume(t - dt, p, phase).unwrap())
                / (2. * dt);
            let vp = (model.volume(t, p + dp, phase).unwrap()
                - model.volume(t, p - dp, phase).unwrap())
                / (2. * dp);
            close(vt, -pressure.dt / pressure.dv, 1e-9, 2e-5);
            close(vp, 1. / pressure.dv, 1e-14, 2e-5);
        }
    }
}

#[test]
fn critical_scaling_reference_independence_and_unclosed_source_iterates() {
    let mut model = Pets::new().unwrap();
    let c = model.critical_point().unwrap();
    let s = model.saturation_pressure(c.temperature * 0.7).unwrap();
    model.ideal_gas_mut().set_formation_enthalpy(1e20).unwrap();
    let unchanged = model.saturation_pressure(c.temperature * 0.7).unwrap();
    close(s.pressure, unchanged.pressure, 1e-7, 1e-12);
    close(s.dpdt, unchanged.dpdt, 1e-7, 1e-12);
    let (sigma, eps) = model.parameters();
    model.set_parameters(2. * sigma, 1.5 * eps).unwrap();
    let scaled = model.critical_point().unwrap();
    close(scaled.temperature, c.temperature * 1.5, 1e-10, 1e-12);
    close(scaled.volume, c.volume * 8., 1e-14, 1e-12);
    close(scaled.pressure, c.pressure * 1.5 / 8., 1e-7, 1e-12);
    let model = Pets::new().unwrap();
    let last = model
        .baseline_volume_iteration(68., 0.001, Phase::Vapor)
        .unwrap();
    assert!(!last.pressure_converged && last.iterations <= 41);
    assert!(last.volume < 1e-3);
    let vapor = model.volume(68., 0.001, Phase::Vapor).unwrap();
    assert!(vapor > 1e5);
    close(
        model.pressure_tv(68., vapor, 1.).unwrap().value,
        0.001,
        1e-15,
        1e-10,
    );
    let c = model.critical_point().unwrap();
    close(
        model
            .volume(c.temperature, c.pressure, Phase::MinimumGibbs)
            .unwrap(),
        c.volume,
        1e-15,
        1e-12,
    );
    assert!(model.saturation_pressure(c.temperature).is_err());
    assert!(model.saturation_temperature(c.pressure).is_err());
    assert!(model.volume(300., 0., Phase::Liquid).is_err());
}
