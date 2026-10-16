mod support;
use support::fixture_model;
use thermopack_rs::{Phase, TranslationKind};

#[test]
fn translation_parameters_match_source_database_and_polynomials() {
    let mut count = 0;
    for row in include_str!("fixtures/translation_parameters.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let t = fields[2].parse().unwrap();
        let mut actual = Vec::new();
        for i in 0..model.components().len() {
            let parameter = model.volume_translation(i).unwrap();
            let curve = parameter.evaluate(t).unwrap();
            actual.extend([
                if model.volume_translation_enabled() {
                    1.
                } else {
                    0.
                },
                parameter.kind as u8 as f64,
            ]);
            actual.extend(parameter.coefficients);
            actual.extend([
                curve.value,
                curve.dt,
                curve.dtt,
                if parameter.kind == TranslationKind::Constant {
                    0.
                } else {
                    1.
                },
            ]);
        }
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if [0, 1, 11].contains(&(i % 12)) {
                0.
            } else {
                1e-22 + 2e-12 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 1977);
}

#[test]
fn translated_cubic_properties_and_derivatives_match_source() {
    let absolute = [1e-14, 1e-2, 1e-6, 1e-2, 1e-7, 1e-9, 1e-7, 1e-9, 1e-7];
    let mut count = 0;
    for row in include_str!("fixtures/translation_properties.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let z = fields[4]
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect::<Vec<f64>>();
        let n = z.len();
        let state = model.prepare(fields[2].parse().unwrap(), &z).unwrap();
        let p = fields[3].parse().unwrap();
        let phase = if fields[5].split('@').next().unwrap().ends_with('L') {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let v = state.volume(p, phase).unwrap();
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
    assert_eq!(count, 3120);
}

#[test]
fn translated_flashes_match_source_and_close_caloric_specifications() {
    for (fixture, expected_count) in [
        (include_str!("fixtures/translated_tp.csv"), 1130),
        (include_str!("fixtures/translated_ph.csv"), 1293),
        (include_str!("fixtures/translated_ps.csv"), 1326),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let model = fixture_model(&fields);
            let mode = fields[5].split('@').next().unwrap();
            let target: f64 = fields[2].parse().unwrap();
            let p = fields[3].parse().unwrap();
            let z: Vec<f64> = fields[4].split(',').map(|v| v.parse().unwrap()).collect();
            let flash = match mode {
                "V_TP" => model.tp_flash(target, p, &z),
                "V_PH" => model.ph_flash(p, &z, target),
                "V_PS" => model.ps_flash(p, &z, target),
                _ => panic!("flash mode"),
            }
            .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
            let mut actual = Vec::new();
            if mode != "V_TP" {
                actual.push(flash.temperature);
            }
            let flag_index = actual.len();
            actual.extend([
                flash.phase as i32 as f64,
                flash.vapor_fraction,
                flash.liquid_fraction,
            ]);
            actual.extend(flash.liquid.iter().chain(&flash.vapor));
            assert_eq!(actual.len(), columns.len() - 2);
            for (i, (a, b)) in actual.iter().zip(&columns[2..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                let tolerance = if i == flag_index {
                    0.
                } else {
                    2e-6 + 2e-7 * a.abs().max(b.abs())
                };
                assert!(
                    (a - b).abs() <= tolerance,
                    "{} quantity {i}: {a}, {b}",
                    columns[0]
                );
            }
            if flash.phase == thermopack_rs::FlashPhase::TwoPhase {
                assert!((flash.vapor_fraction + flash.liquid_fraction - 1.).abs() < 1e-10);
                for (i, &zi) in z.iter().enumerate() {
                    let balance = flash.vapor_fraction * flash.vapor[i]
                        + flash.liquid_fraction * flash.liquid[i];
                    assert!((balance - zi).abs() < 1e-9, "{} component {i}", columns[0]);
                }
            }
            if mode != "V_TP" {
                let props = model.equilibrium_properties(&flash).unwrap();
                let (value, scale) = if mode == "V_PH" {
                    (props.enthalpy, thermopack_rs::R * flash.temperature)
                } else {
                    (props.entropy, thermopack_rs::R)
                };
                assert!(
                    (value - target).abs() <= 2e-8 * scale,
                    "{} caloric closure",
                    columns[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, expected_count);
    }
}

#[test]
fn translated_saturation_matches_source() {
    let mut count = 0;
    let absolute = [1e-6, 1e-3, 1e-13, 1e-11, 1e-3];
    for row in include_str!("fixtures/translated_sat.csv").lines().skip(1) {
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
        assert_eq!(actual.len(), columns.len() - 2);
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
    assert_eq!(count, 72);
}
