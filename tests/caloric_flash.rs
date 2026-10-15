use thermopack_rs::{Cubic, CubicKind};

#[test]
fn saturation_temperature_and_pseudo_critical() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for name in ["C1", "CO2", "C3"] {
            let model = Cubic::new(kind, &[name]).unwrap();
            let tc = model.components()[0].critical_temperature;
            let pc = model.components()[0].critical_pressure;
            let critical = model.pseudo_critical(&[1.0]).unwrap();
            assert!((critical.temperature / tc - 1.0).abs() < 1e-12);
            assert!((critical.pressure / pc - 1.0).abs() < 1e-12);
            for tr in [0.65, 0.8, 0.95] {
                let sat = model.saturation_pressure(tc * tr).unwrap();
                let inverse = model.saturation_temperature(sat.pressure).unwrap();
                assert!(
                    (inverse.temperature - tc * tr).abs() < 1e-7 * tc,
                    "{kind:?} {name} {tr}"
                );
            }
        }
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let critical = model.pseudo_critical(&[0.5, 0.3, 0.2]).unwrap();
        let state = model
            .prepare(critical.temperature, &[0.5, 0.3, 0.2])
            .unwrap();
        assert!(
            (state.pressure(critical.volume).unwrap().value / critical.pressure - 1.0).abs()
                < 1e-12
        );
        let pv = state.pressure(critical.volume).unwrap().dv;
        assert!(pv.abs() < 1e-8 * critical.pressure / critical.volume);
    }
}

#[test]
fn equilibrium_heat_capacity_and_caloric_round_trips() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.5, 0.3, 0.2];
        for t in [160.0, 200.0, 240.0, 300.0, 450.0] {
            let p = 1e6;
            let flash = model.tp_flash(t, p, &z).unwrap();
            let properties = model.equilibrium_properties(&flash).unwrap();
            let dt = 0.001;
            let plus = model
                .equilibrium_properties(&model.tp_flash(t + dt, p, &z).unwrap())
                .unwrap();
            let minus = model
                .equilibrium_properties(&model.tp_flash(t - dt, p, &z).unwrap())
                .unwrap();
            let fd = (plus.enthalpy - minus.enthalpy) / (2.0 * dt);
            let cp = properties.heat_capacity.unwrap();
            let dsdt = (plus.entropy - minus.entropy) / (2.0 * dt);
            assert!(
                (cp - fd).abs() < 1e-5 * cp.abs(),
                "{kind:?} {t}: Cp={cp}, finite difference={fd}"
            );
            assert!((cp / t - dsdt).abs() < 1e-5 * (cp / t).abs());
            let ph = model
                .ph_flash(p, &z, properties.enthalpy)
                .unwrap_or_else(|e| panic!("PH {kind:?} {t}: {e}"));
            let ps = model
                .ps_flash(p, &z, properties.entropy)
                .unwrap_or_else(|e| panic!("PS {kind:?} {t}: {e}"));
            assert!(
                (ph.temperature - t).abs() < 1e-5,
                "PH {kind:?} {t}: {}",
                ph.temperature
            );
            assert!(
                (ps.temperature - t).abs() < 1e-5,
                "PS {kind:?} {t}: {}",
                ps.temperature
            );
        }
    }
}

#[test]
fn pure_latent_heat_quality() {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2"]).unwrap();
    let sat = model.saturation_pressure(250.0).unwrap();
    let state = model.prepare(250.0, &[1.0]).unwrap();
    let hl = state.enthalpy(sat.liquid_volume).unwrap();
    let hv = state.enthalpy(sat.vapor_volume).unwrap();
    let sl = state.entropy(sat.liquid_volume).unwrap();
    let sv = state.entropy(sat.vapor_volume).unwrap();
    for beta in [0.0, 0.1, 0.5, 0.9, 1.0] {
        let ph = model
            .ph_flash(sat.pressure, &[1.0], hl * (1.0 - beta) + hv * beta)
            .unwrap();
        let ps = model
            .ps_flash(sat.pressure, &[1.0], sl * (1.0 - beta) + sv * beta)
            .unwrap();
        for flash in [ph, ps] {
            assert!((flash.temperature - 250.0).abs() < 1e-7);
            assert!((flash.vapor_fraction - beta).abs() < 1e-8);
        }
    }
}

#[test]
fn fixed_composition_helpers_match_exact_source() {
    let mut count = 0;
    for row in include_str!("fixtures/initialization.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let kind = match fields[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => panic!("unknown model"),
        };
        let model = Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>()).unwrap();
        let actual = if fields[3] == "PSEUDO" {
            let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
            let critical = model.pseudo_critical(&z).unwrap();
            vec![
                critical.temperature,
                critical.pressure,
                critical.volume,
                critical.compressibility,
            ]
        } else {
            let sat = model
                .saturation_temperature(fields[2].parse().unwrap())
                .unwrap();
            vec![
                sat.temperature,
                sat.pressure,
                sat.liquid_volume,
                sat.vapor_volume,
                sat.dpdt,
            ]
        };
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let abs = if i == 2 || i == 3 { 1e-12 } else { 1e-6 };
            assert!(
                (a - b).abs() <= abs + 2e-8 * a.abs().max(b.abs()),
                "{} quantity {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 64);
}

#[test]
fn matches_committed_caloric_reference_cases() {
    let mut cases = 0;
    for row in include_str!("fixtures/caloric_flash.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let kind = match fields[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => panic!("unknown fixture model"),
        };
        let model = Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>()).unwrap();
        let p: f64 = fields[3].parse().unwrap();
        let z: Vec<f64> = fields[4].split(',').map(|v| v.parse().unwrap()).collect();
        let specification: f64 = fields[2].parse().unwrap();
        let flash = if fields[5] == "PH" {
            model.ph_flash(p, &z, specification)
        } else {
            model.ps_flash(p, &z, specification)
        }
        .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
        let mut values = vec![
            flash.temperature,
            flash.phase as i32 as f64,
            flash.vapor_fraction,
            flash.liquid_fraction,
        ];
        values.extend(flash.liquid.iter().chain(&flash.vapor));
        assert_eq!(values.len(), columns.len() - 2);
        for (i, (actual, expected)) in values.iter().zip(&columns[2..]).enumerate() {
            let expected: f64 = expected.parse().unwrap();
            let tolerance = if i == 1 {
                0.0
            } else {
                2e-6 + 2e-7 * actual.abs().max(expected.abs())
            };
            assert!(
                (actual - expected).abs() <= tolerance,
                "{} quantity {i}: Rust {actual}, reference {expected}",
                columns[0]
            );
        }
        let props = model.equilibrium_properties(&flash).unwrap();
        let actual = if fields[5] == "PH" {
            props.enthalpy
        } else {
            props.entropy
        };
        let scale = thermopack_rs::R
            * if fields[5] == "PH" {
                flash.temperature
            } else {
                1.0
            };
        assert!(
            (actual - specification).abs() <= 2e-8 * scale,
            "{} caloric closure",
            columns[0]
        );
        cases += 1;
    }
    assert_eq!(cases, 697);
}
