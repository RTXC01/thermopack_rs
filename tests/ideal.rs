use thermopack_rs::{Cubic, CubicKind, R, component, components};

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
