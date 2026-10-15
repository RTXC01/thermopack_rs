use thermopack_rs::{Cubic, CubicKind, EntropyReference, IdealGas, R, component, components};

fn configure(gas: &mut IdealGas, actions: &str) {
    for action in actions.split('/') {
        let fields: Vec<_> = action.split(':').collect();
        match fields[0] {
            "cp" => gas
                .set_heat_capacity(
                    fields[1].parse().unwrap(),
                    &fields[2]
                        .split(',')
                        .filter(|v| !v.is_empty())
                        .map(|v| v.parse().unwrap())
                        .collect::<Vec<_>>(),
                )
                .unwrap(),
            "h" => gas
                .set_formation_enthalpy(fields[1].parse().unwrap())
                .unwrap(),
            "s" => gas
                .set_standard_entropy(
                    fields[1].parse().unwrap(),
                    if fields[2].parse::<f64>().unwrap() == 101325. {
                        EntropyReference::Atmosphere
                    } else {
                        EntropyReference::Bar
                    },
                )
                .unwrap(),
            _ => panic!("unknown action"),
        }
    }
}

#[test]
fn configured_ideal_correlations_and_reference_state_parity() {
    let absolute = [1e-7, 1e-9, 1e-9, 1e-10, 1e-12, 1e-10, 1e-10];
    let mut count = 0;
    for row in include_str!("fixtures/ideal_configuration.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let t: f64 = fields[2].parse().unwrap();
        let p: f64 = fields[3].parse().unwrap();
        let mut gas = IdealGas::new(component(fields[1]).unwrap()).unwrap();
        configure(&mut gas, fields[4]);
        let properties = gas.properties(t, p).unwrap();
        let (corr, params) = gas.heat_capacity_correlation();
        let mut values = vec![
            properties.enthalpy,
            properties.entropy,
            properties.cp,
            properties.cp / t,
            -R / p,
            gas.formation_enthalpy(),
            gas.standard_entropy(),
            corr as f64,
        ];
        values.extend(params);
        assert_eq!(values.len(), columns.len() - 1);
        for (i, (a, b)) in values.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i < 7 {
                absolute[i] + 2e-8 * a.abs().max(b.abs())
            } else {
                0.0
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        let mut model = Cubic::new(CubicKind::PengRobinson, &[fields[1]]).unwrap();
        configure(model.ideal_gas_mut(0).unwrap(), fields[4]);
        let active = model.prepare(t, &[1.0]).unwrap().ideal(p).unwrap();
        assert_eq!(active.enthalpy, properties.enthalpy);
        assert_eq!(active.entropy, properties.entropy);
        assert_eq!(active.cp, properties.cp);
        count += 1;
    }
    assert_eq!(count, 720);
}

#[test]
fn ideal_configuration_is_owned_and_invalid_setters_are_atomic() {
    let mut model = Cubic::new(CubicKind::PengRobinson, &["C1", "CO2"]).unwrap();
    let other = model.clone();
    let initial = other.ideal_gas(0).unwrap().properties(300., 1e5).unwrap();
    let gas = model.ideal_gas_mut(0).unwrap();
    gas.set_heat_capacity(8, &[2.5]).unwrap();
    gas.set_formation_enthalpy(-1234.).unwrap();
    gas.set_standard_entropy(123., EntropyReference::Atmosphere)
        .unwrap();
    let at_reference = gas.properties(298.15, 101325.).unwrap();
    assert!((at_reference.enthalpy + 1234.).abs() < 1e-9);
    assert!((at_reference.entropy - 123.).abs() < 1e-9);
    let parameters = *gas.heat_capacity_correlation().1;
    for result in [
        gas.set_heat_capacity(0, &[]),
        gas.set_heat_capacity(8, &[f64::NAN]),
        gas.set_heat_capacity(8, &[1.; 22]),
        gas.set_heat_capacity(10, &[]),
        gas.set_formation_enthalpy(f64::INFINITY),
        gas.set_standard_entropy(f64::NAN, EntropyReference::Bar),
    ] {
        assert!(result.is_err());
    }
    assert_eq!(gas.heat_capacity_correlation().0, 8);
    assert_eq!(*gas.heat_capacity_correlation().1, parameters);
    assert_eq!(
        gas.properties(298.15, 101325.).unwrap().enthalpy,
        at_reference.enthalpy
    );
    assert_eq!(
        other
            .ideal_gas(0)
            .unwrap()
            .properties(300., 1e5)
            .unwrap()
            .enthalpy,
        initial.enthalpy
    );
    assert!(model.ideal_gas(2).is_err());
}

#[test]
fn configured_caloric_values_close_native_flash_targets() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let mut model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        for i in 0..3 {
            let gas = model.ideal_gas_mut(i).unwrap();
            gas.set_heat_capacity(8, &[3.0 + i as f64]).unwrap();
            gas.set_formation_enthalpy(-10000. * (i + 1) as f64)
                .unwrap();
            gas.set_standard_entropy(100. + 10. * i as f64, EntropyReference::Bar)
                .unwrap();
        }
        for t in [180., 200., 250.] {
            let z = [0.5, 0.3, 0.2];
            let props = model
                .equilibrium_properties(&model.tp_flash(t, 1e6, &z).unwrap())
                .unwrap();
            let ph = model.ph_flash(1e6, &z, props.enthalpy).unwrap();
            let ps = model.ps_flash(1e6, &z, props.entropy).unwrap();
            assert!((ph.temperature - t).abs() < 1e-6, "PH {kind:?} {t}");
            assert!((ps.temperature - t).abs() < 1e-6, "PS {kind:?} {t}");
        }
    }
}

#[test]
fn dippr_shared_exponentials_keep_low_and_high_temperature_limits() {
    let mut gas = IdealGas::new(component("C1").unwrap()).unwrap();
    gas.set_heat_capacity(7, &[30000., 20000., 700., 15000., 300.])
        .unwrap();
    assert!((gas.properties(1., 1e5).unwrap().cp - 30.).abs() < 1e-12);
    assert!((gas.properties(1e9, 1e5).unwrap().cp - 50.).abs() < 1e-10);
    for t in [1., 1400., 1e9] {
        let dt = t * 1e-5;
        let p = gas.properties(t + dt, 1e5).unwrap();
        let m = gas.properties(t - dt, 1e5).unwrap();
        let at = gas.properties(t, 1e5).unwrap();
        assert!(at.enthalpy.is_finite() && at.entropy.is_finite());
        assert!(((p.enthalpy - m.enthalpy) / (2. * dt) - at.cp).abs() < 1e-6 * at.cp);
        assert!(((p.entropy - m.entropy) / (2. * dt) - at.cp / t).abs() < 1e-6 * at.cp / t);
    }
}

#[test]
fn all_finite_default_ideal_reference_cases() {
    let absolute = [1e-7, 1e-9, 1e-9, 1e-10, 1e-12];
    let mut count = 0;
    for row in include_str!("fixtures/ideal.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let t: f64 = fields[2].parse().unwrap();
        let p: f64 = fields[3].parse().unwrap();
        let ideal = component(fields[1]).unwrap().ideal(t, p).unwrap();
        let actual = [
            ideal.enthalpy,
            ideal.entropy,
            ideal.cp,
            ideal.cp / t,
            -R / p,
        ];
        for (i, (a, reference)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = reference.parse().unwrap();
            assert!(
                (a - b).abs() <= absolute[i] + 2e-8 * a.abs().max(b.abs()),
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 2592);
}

#[test]
fn ideal_value_derivatives_for_every_supported_default_fluid() {
    let mut supported = 0;
    for c in components() {
        if c.ident == "BUTANAL" {
            // The source default is correlation 100 (invalid), and STOPs.
            assert!(c.ideal(300.0, 1e5).is_err());
            continue;
        }
        for t in [80.0, 298.15, 500.0] {
            let ideal = c.ideal(t, 1e5).unwrap();
            let dt = t * 1e-5;
            let plus = c.ideal(t + dt, 1e5).unwrap();
            let minus = c.ideal(t - dt, 1e5).unwrap();
            let dhdt = (plus.enthalpy - minus.enthalpy) / (2.0 * dt);
            let dsdt = (plus.entropy - minus.entropy) / (2.0 * dt);
            let capacity = ideal.enthalpy_dt;
            assert!(
                (dhdt - capacity).abs() <= 1e-6 * capacity.abs().max(1.0),
                "{} {t}: H derivative {dhdt}, {capacity}",
                c.ident
            );
            assert!(
                (dsdt - capacity / t).abs() <= 1e-6 * (capacity / t).abs().max(1.0),
                "{} {t}: S derivative {dsdt}, {}",
                c.ident,
                capacity / t
            );
        }
        supported += 1;
    }
    assert_eq!(supported, 108);
}

#[test]
fn hydrogen_source_values_and_actual_derivatives_remain_explicit() {
    for name in ["N-H2", "E-H2"] {
        for kind in [
            CubicKind::PengRobinson,
            CubicKind::SoaveRedlichKwong,
            CubicKind::RedlichKwong,
            CubicKind::VanDerWaals,
        ] {
            let model = Cubic::new(kind, &[name]).unwrap();
            let state = model.prepare(150.0, &[1.0]).unwrap();
            let v = state.volume(1e6, thermopack_rs::Phase::Vapor).unwrap();
            let ideal = state.ideal(1e6).unwrap();
            assert_eq!(ideal.enthalpy_dt, 1000.0 * ideal.cp);
            let cp = state.heat_capacities(v).unwrap().0;
            let cp_baseline = state.baseline_enthalpy_derivative_tp(v).unwrap();
            if kind != CubicKind::VanDerWaals {
                assert!((cp - cp_baseline - 999.0 * ideal.cp).abs() < 1e-8);
            }
            let dt = 0.001;
            let h_at = |t| {
                let s = model.prepare(t, &[1.0]).unwrap();
                s.enthalpy(s.volume(1e6, thermopack_rs::Phase::Vapor).unwrap())
                    .unwrap()
            };
            assert!(((h_at(150.0 + dt) - h_at(150.0 - dt)) / (2.0 * dt) - cp).abs() < 1e-6 * cp);
            let flash = model.tp_flash(150.0, 1e6, &[1.0]).unwrap();
            let props = model.equilibrium_properties(&flash).unwrap();
            assert!(
                (model
                    .ph_flash(1e6, &[1.0], props.enthalpy)
                    .unwrap()
                    .temperature
                    - 150.0)
                    .abs()
                    < 1e-6
            );
            assert!(
                (model
                    .ps_flash(1e6, &[1.0], props.entropy)
                    .unwrap()
                    .temperature
                    - 150.0)
                    .abs()
                    < 1e-6
            );
        }
    }
}
