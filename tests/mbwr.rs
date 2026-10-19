use thermopack_rs::{Mbwr, MbwrKind, R, mbwr_records};

fn values(model: &Mbwr, t: f64, v: f64, n: f64) -> Vec<f64> {
    let prepared = model.prepare(t).unwrap();
    let a = prepared.baseline_residual_density(n / v).unwrap();
    let f = prepared.baseline_residual_helmholtz(v, n).unwrap();
    let p = prepared.baseline_pressure_tv(v, n).unwrap();
    let d = prepared.pressure_density(n / v).unwrap();
    vec![
        a.value, a.dt, a.drho, a.dtt, a.dtrho, a.drhorho, f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv,
        f.dtn, f.dvv, f.dvn, f.dnn, p.value, p.dt, p.dv, p.dn, d.value, d.dt, d.drho, d.dtt,
        d.dtrho, d.drhorho,
    ]
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn all_original_parameters_and_kernel_outputs_match_source() {
    for (fixture, expected) in [
        (include_str!("fixtures/mbwr_parameters.csv"), 65),
        (include_str!("fixtures/mbwr_integrated.csv"), 65),
        (include_str!("fixtures/mbwr_core.csv"), 6760),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let model = Mbwr::from_record(fields[1]).unwrap();
            if fields[5] == "MBWR_INTEGRATED" {
                let actual: Vec<_> = model
                    .baseline_helmholtz_terms()
                    .iter()
                    .map(|t| t.coefficient)
                    .collect();
                assert_eq!(actual.len(), columns.len() - 1);
                for (&a, b) in actual.iter().zip(&columns[1..]) {
                    assert_eq!(a, b.parse::<f64>().unwrap());
                }
            } else if fields[5] == "MBWR_PARAMETERS" {
                let (p, gamma) = model.parameters();
                let mut actual = vec![model.kind().coefficient_count() as f64];
                actual.extend(p);
                actual.push(gamma);
                assert_eq!(actual.len(), columns.len() - 1);
                for (&a, b) in actual.iter().zip(&columns[1..]) {
                    assert_eq!(a, b.parse::<f64>().unwrap());
                }
            } else {
                let actual = values(
                    &model,
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                    fields[4].parse().unwrap(),
                );
                let atol = [
                    2e-8, 2e-9, 1e-10, 2e-10, 1e-11, 1e-13, 1e-7, 1e-8, 1e-7, 1e-7, 1e-9, 1e-8,
                    1e-8, 1e-4, 1e-7, 1e-7, 0.001, 0.001, 0.1, 0.001, 0.001, 0.001, 1e-5, 1e-5,
                    1e-7, 1e-9,
                ];
                assert_eq!(actual.len(), columns.len() - 1);
                for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                    let b: f64 = b.parse().unwrap();
                    assert!(
                        (a - b).abs() <= atol[i] + 2e-8 * a.abs().max(b.abs()),
                        "{} field {i}: {a} != {b}",
                        columns[0]
                    );
                }
            }
            count += 1;
        }
        assert_eq!(count, expected);
    }
}

#[test]
fn density_derivatives_and_pressure_obey_finite_differences_and_helmholtz_identity() {
    for record in mbwr_records() {
        let model = Mbwr::from_record(record.source_name).unwrap();
        let t = 300.;
        let rho = 1000.;
        let dt = 0.001;
        let dr = 0.01;
        let a = model.residual_density(t, rho).unwrap();
        let lo = model.residual_density(t - dt, rho).unwrap();
        let hi = model.residual_density(t + dt, rho).unwrap();
        close(a.dt, (hi.value - lo.value) / (2. * dt), 1e-7, 1e-5);
        close(a.dtt, (hi.dt - lo.dt) / (2. * dt), 1e-9, 1e-5);
        close(a.dtrho, (hi.drho - lo.drho) / (2. * dt), 1e-9, 1e-5);
        let lo = model.residual_density(t, rho - dr).unwrap();
        let hi = model.residual_density(t, rho + dr).unwrap();
        close(a.drho, (hi.value - lo.value) / (2. * dr), 1e-8, 1e-5);
        close(a.drhorho, (hi.drho - lo.drho) / (2. * dr), 1e-10, 1e-5);
        close(a.dtrho, (hi.dt - lo.dt) / (2. * dr), 1e-9, 1e-5);
        let p = model.pressure_density(t, rho).unwrap();
        close(p.value, R * t * rho * (1. + rho * a.drho), 0.001, 1e-11);
        close(
            p.dt,
            R * rho * (1. + rho * a.drho) + R * t * rho * rho * a.dtrho,
            0.001,
            1e-11,
        );
        close(
            p.drho,
            R * t * (1. + 2. * rho * a.drho + rho * rho * a.drhorho),
            0.001,
            1e-11,
        );
        let lo = model.pressure_density(t - dt, rho).unwrap();
        let hi = model.pressure_density(t + dt, rho).unwrap();
        close(p.dt, (hi.value - lo.value) / (2. * dt), 0.01, 1e-5);
        close(p.dtt, (hi.dt - lo.dt) / (2. * dt), 0.001, 1e-5);
        close(p.dtrho, (hi.drho - lo.drho) / (2. * dt), 1e-5, 1e-5);
        let lo = model.pressure_density(t, rho - dr).unwrap();
        let hi = model.pressure_density(t, rho + dr).unwrap();
        close(p.drho, (hi.value - lo.value) / (2. * dr), 0.01, 1e-5);
        close(p.drhorho, (hi.drho - lo.drho) / (2. * dr), 1e-5, 1e-5);
        close(p.dtrho, (hi.dt - lo.dt) / (2. * dr), 1e-5, 1e-5);
    }
}

#[test]
fn extensive_derivatives_euler_scaling_and_integral_switch_are_consistent() {
    for source in ["mbwrcx1", "mbwrcx3", "mbwrcx12", "mbwrcx58", "mbwrcx64"] {
        let model = Mbwr::from_record(source).unwrap();
        let t = 400.;
        let v = 1e-3;
        let n = 2.4;
        let f = model.residual_helmholtz(t, v, n).unwrap();
        close(f.value, n * f.dn + v * f.dv, 1e-9, 1e-12);
        let scale = 3.5;
        let g = model.residual_helmholtz(t, v * scale, n * scale).unwrap();
        close(g.value, scale * f.value, 1e-9, 1e-12);
        close(g.dv, f.dv, 1e-7, 1e-12);
        close(g.dnn, f.dnn / scale, 1e-9, 1e-12);
        let dt = 0.001;
        let dv = v * 1e-5;
        let dn = n * 1e-5;
        let lo = model.residual_helmholtz(t - dt, v, n).unwrap();
        let hi = model.residual_helmholtz(t + dt, v, n).unwrap();
        close(f.dt, (hi.value - lo.value) / (2. * dt), 1e-8, 1e-5);
        close(f.dtt, (hi.dt - lo.dt) / (2. * dt), 1e-10, 1e-5);
        close(f.dtv, (hi.dv - lo.dv) / (2. * dt), 1e-5, 1e-5);
        close(f.dtn, (hi.dn - lo.dn) / (2. * dt), 1e-8, 1e-5);
        let lo = model.residual_helmholtz(t, v - dv, n).unwrap();
        let hi = model.residual_helmholtz(t, v + dv, n).unwrap();
        close(f.dv, (hi.value - lo.value) / (2. * dv), 1e-4, 1e-5);
        close(f.dvv, (hi.dv - lo.dv) / (2. * dv), 0.1, 1e-5);
        close(f.dvn, (hi.dn - lo.dn) / (2. * dv), 1e-4, 1e-5);
        let lo = model.residual_helmholtz(t, v, n - dn).unwrap();
        let hi = model.residual_helmholtz(t, v, n + dn).unwrap();
        close(f.dn, (hi.value - lo.value) / (2. * dn), 1e-8, 1e-5);
        close(f.dnn, (hi.dn - lo.dn) / (2. * dn), 1e-8, 1e-5);
        let gamma = model.parameters().1;
        let rho = 1000. * (4. / gamma).sqrt();
        let lo = model.residual_density(t, rho * (1. - 1e-7)).unwrap();
        let hi = model.residual_density(t, rho * (1. + 1e-7)).unwrap();
        close(
            (hi.value - lo.value) / (rho * 2e-7),
            model.residual_density(t, rho).unwrap().drho,
            1e-7,
            1e-5,
        );
    }
}

#[test]
fn source_order_owned_parameters_zero_density_and_invalid_inputs() {
    assert_eq!(mbwr_records().len(), 65);
    assert_eq!(
        Mbwr::new(MbwrKind::Mbwr32, "methane")
            .unwrap()
            .record()
            .unwrap()
            .source_name,
        "mbwrcx1"
    );
    assert_eq!(
        Mbwr::new(MbwrKind::Bender19, "C1")
            .unwrap()
            .record()
            .unwrap()
            .source_name,
        "mbwrcx58"
    );
    let original = Mbwr::from_record("mbwrcx1").unwrap();
    let mut coefficients = original.parameters().0.to_vec();
    let copy =
        Mbwr::with_parameters(MbwrKind::Mbwr32, &coefficients, original.parameters().1).unwrap();
    coefficients[0] = 0.;
    assert_eq!(original.parameters().0, copy.parameters().0);
    assert!(copy.record().is_none());
    for record in mbwr_records() {
        let model = Mbwr::from_record(record.source_name).unwrap();
        let a = model.residual_density(300., 0.).unwrap();
        assert_eq!([a.value, a.dt, a.dtt], [0.; 3]);
        let (b, c) = model.virial_coefficients(300.).unwrap();
        assert!(b.is_finite() && c.is_finite());
        close(
            model.residual_density(300., 1e-8).unwrap().value / 1e-8,
            b,
            1e-10,
            1e-7,
        );
        let p = model.pressure_density(300., 0.).unwrap();
        assert_eq!(p.value, 0.);
        assert_eq!(p.drho, R * 300.);
    }
    let ideal = Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 1.).unwrap();
    let f = ideal.residual_density(300., 1000.).unwrap();
    assert_eq!([f.value, f.dt, f.drho, f.dtt, f.dtrho, f.drhorho], [0.; 6]);
    close(
        ideal.pressure_density(300., 1000.).unwrap().value,
        R * 300. * 1000.,
        0.,
        0.,
    );
    assert!(Mbwr::with_parameters(MbwrKind::Mbwr32, &[0.; 19], 1.).is_err());
    assert!(Mbwr::with_parameters(MbwrKind::Bender19, &[f64::NAN; 19], 1.).is_err());
    assert!(Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 0.).is_err());
    assert!(Mbwr::from_record("unknown").is_err());
    assert!(Mbwr::new(MbwrKind::Mbwr32, "unknown").is_err());
    assert!(original.prepare(0.).is_err());
    assert!(original.residual_density(300., -1.).is_err());
    assert!(original.pressure_density(300., f64::INFINITY).is_err());
    assert!(original.residual_helmholtz(300., 1., 0.).is_err());
}

#[test]
fn stable_integral_derivatives_remain_consistent_in_dilute_cryogenic_extrapolations() {
    for source in ["mbwrcx1", "mbwrcx2", "mbwrcx7", "mbwrcx10", "mbwrcx11"] {
        let model = Mbwr::from_record(source).unwrap();
        for t in [0.5, 1., 5.] {
            for rho in [1e-5, 0.01, 1.] {
                let a = model.residual_density(t, rho).unwrap();
                let dt = t * 1e-5;
                let dr = rho * 1e-5;
                let lo = model.residual_density(t - dt, rho).unwrap();
                let hi = model.residual_density(t + dt, rho).unwrap();
                close(a.dt, (hi.value - lo.value) / (2. * dt), 1e-6, 2e-7);
                close(a.dtt, (hi.dt - lo.dt) / (2. * dt), 1e-5, 2e-7);
                let lo = model.residual_density(t, rho - dr).unwrap();
                let hi = model.residual_density(t, rho + dr).unwrap();
                close(a.drho, (hi.value - lo.value) / (2. * dr), 1e-5, 2e-7);
                let p = model.pressure_density(t, rho).unwrap();
                close(p.value, R * t * rho * (1. + rho * a.drho), 1e-7, 2e-11);
            }
        }
    }
}
