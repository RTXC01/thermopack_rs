mod support;
use support::fixture_model;
use thermopack_rs::{BetaCorrelation, Phase};

#[test]
fn quantum_cubic_pure_saturation_matches_source() {
    let absolute = [1e-6, 1e-3, 1e-13, 1e-11, 1e-3];
    let mut count = 0;
    for row in include_str!("fixtures/quantum_sat.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let sat = model
            .saturation_pressure(fields[2].parse().unwrap())
            .unwrap();
        let actual = [
            sat.temperature,
            sat.pressure,
            sat.liquid_volume,
            sat.vapor_volume,
            sat.dpdt,
        ];
        for (i, (a, b)) in actual.iter().zip(&columns[2..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            assert!(
                (a - b).abs() <= absolute[i] + 2e-7 * a.abs().max(b.abs()),
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 36);
}

#[test]
fn beta_parameters_derivatives_and_constructor_data_match_source() {
    let mut count = 0;
    for row in include_str!("fixtures/beta_parameters.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let t = fields[2].parse().unwrap();
        let mut actual = Vec::new();
        for (i, component) in model.components().iter().enumerate() {
            let beta = model.beta(i).unwrap();
            let value = beta.evaluate(t, component.critical_temperature).unwrap();
            actual.push(beta.correlation as i32 as f64);
            actual.extend(if beta.correlation == BetaCorrelation::Classic {
                [0.; 2]
            } else {
                beta.parameters
            });
            actual.extend([
                value.value,
                value.dt,
                value.dtt,
                component.critical_temperature,
                component.critical_pressure,
                component.molar_mass,
            ]);
        }
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if [0, 1, 2, 6, 7, 8].contains(&(i % 9)) {
                0.
            } else {
                1e-14 + 2e-12 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 4036);
}

#[test]
fn beta_and_quantum_cubic_properties_and_fugacity_derivatives_match_source() {
    let absolute = [1e-14, 1e-2, 1e-6, 1e-2, 1e-7, 1e-9, 1e-7, 1e-9, 1e-7];
    let mut count = 0;
    for row in include_str!("fixtures/beta_properties.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let z: Vec<f64> = fields[4].split(',').map(|v| v.parse().unwrap()).collect();
        let n = z.len();
        let state = model.prepare(fields[2].parse().unwrap(), &z).unwrap();
        let phase = if fields[5].split('@').next().unwrap().ends_with('L') {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let v = state.volume(fields[3].parse().unwrap(), phase).unwrap();
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
            let ab = if i < 9 {
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
                (a - b).abs() <= ab + rel * a.abs().max(b.abs()),
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 4512);
}
