use thermopack_rs::{DensityHelmholtz, Pets, PureHelmholtz};

fn density_values(a: DensityHelmholtz) -> [f64; 6] {
    [a.value, a.dt, a.drho, a.dtt, a.dtrho, a.drhorho]
}
fn extensive_values(f: PureHelmholtz) -> [f64; 10] {
    [
        f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn,
    ]
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn pets_matches_the_exact_source() {
    let absolute = [1e-22, 1e-23, 1e-24]
        .into_iter()
        .chain([1e-11, 1e-11, 1e-16, 1e-10, 1e-15, 1e-20].repeat(3))
        .chain(
            [
                1e-10, 1e-10, 1e-8, 1e-10, 1e-9, 1e-7, 1e-9, 1e-3, 1e-7, 1e-9,
            ]
            .repeat(2),
        )
        .chain([
            1e-3, 1e-4, 1., 1e-3, 1e-5, 1e-7, 1e-5, 1e-5, 1e-6, 1e-5, 1e-6, 1e-15, 1e-11, 1e-9,
            1e-6,
        ])
        .collect::<Vec<_>>();
    assert_eq!(absolute.len(), 56);
    for (fixture, count, core) in [
        (include_str!("fixtures/pets_core.csv"), 588, true),
        (include_str!("fixtures/pets_potential.csv"), 16, false),
    ] {
        let mut actual_count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let params: Vec<f64> = fields[1].split(',').map(|v| v.parse().unwrap()).collect();
            let model = Pets::with_parameters(params[0], params[1]).unwrap();
            let actual = if core {
                let t = fields[2].parse().unwrap();
                let v = fields[3].parse().unwrap();
                let n = fields[4].parse().unwrap();
                let prepared = model.prepare(t).unwrap();
                let diameter = prepared.diameter();
                let a = prepared.contributions(n / v).unwrap();
                let f = prepared.residual_helmholtz(v, n).unwrap();
                let c = model.lj_truncation_correction(t, v, n).unwrap();
                let p = prepared.caloric_tv(v, n).unwrap();
                let mut values = vec![diameter.value, diameter.dt, diameter.dtt];
                for a in [a.total, a.hard_sphere, a.dispersion] {
                    values.extend(density_values(a));
                }
                for f in [f, c] {
                    values.extend(extensive_values(f));
                }
                values.extend([
                    p.pressure.value,
                    p.pressure.dt,
                    p.pressure.dv,
                    p.pressure.dn,
                    p.internal_energy,
                    p.entropy,
                    p.helmholtz_energy,
                    p.enthalpy,
                    p.heat_capacity_v,
                    p.chemical_potential,
                ]);
                let (b, c) = model.baseline_virial_coefficients(t).unwrap();
                values.extend([
                    prepared.baseline_speed_of_sound_tv(v, n).unwrap(),
                    b,
                    c,
                    p.ln_fugacity_coefficient.unwrap_or(0.),
                    p.heat_capacity_p.unwrap_or(f64::NAN),
                ]);
                values
            } else {
                let r: Vec<f64> = fields[2].split(',').map(|v| v.parse().unwrap()).collect();
                let mut values = vec![0.; r.len()];
                let mut force = vec![0.; r.len()];
                model
                    .potential_into(&r, &mut values, Some(&mut force))
                    .unwrap();
                values.extend(force);
                values
            };
            assert_eq!(actual.len(), columns.len() - 1);
            for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                let atol = if core {
                    absolute[i]
                } else if i < actual.len() / 2 {
                    1e-10
                } else {
                    1e-3
                };
                assert!(
                    (a - b).abs() <= atol + 5e-10 * a.abs().max(b.abs()),
                    "{} field {i}: {a}, {b}",
                    columns[0]
                );
            }
            actual_count += 1;
        }
        assert_eq!(actual_count, count);
    }
}

#[test]
fn extensive_derivatives_obey_finite_differences_and_euler_identities() {
    let model = Pets::new().unwrap();
    for t in [80., 150., 300., 1000.] {
        for rho in [100., 10000., 25000.] {
            for n in [0.3, 1., 3.5] {
                let v = n / rho;
                let f = model.residual_helmholtz(t, v, n).unwrap();
                let dt = t * 1e-5;
                let dv = v * 1e-5;
                let dn = n * 1e-5;
                let lo = model.residual_helmholtz(t - dt, v, n).unwrap();
                let hi = model.residual_helmholtz(t + dt, v, n).unwrap();
                for (a, b) in [
                    ((hi.value - lo.value) / (2. * dt), f.dt),
                    ((hi.dt - lo.dt) / (2. * dt), f.dtt),
                    ((hi.dv - lo.dv) / (2. * dt), f.dtv),
                    ((hi.dn - lo.dn) / (2. * dt), f.dtn),
                ] {
                    close(a, b, 1e-9, 3e-6);
                }
                let lo = model.residual_helmholtz(t, v - dv, n).unwrap();
                let hi = model.residual_helmholtz(t, v + dv, n).unwrap();
                for (a, b) in [
                    ((hi.value - lo.value) / (2. * dv), f.dv),
                    ((hi.dv - lo.dv) / (2. * dv), f.dvv),
                    ((hi.dt - lo.dt) / (2. * dv), f.dtv),
                    ((hi.dn - lo.dn) / (2. * dv), f.dvn),
                ] {
                    close(a, b, 1e-7, 3e-6);
                }
                let lo = model.residual_helmholtz(t, v, n - dn).unwrap();
                let hi = model.residual_helmholtz(t, v, n + dn).unwrap();
                for (a, b) in [
                    ((hi.value - lo.value) / (2. * dn), f.dn),
                    ((hi.dv - lo.dv) / (2. * dn), f.dvn),
                    ((hi.dt - lo.dt) / (2. * dn), f.dtn),
                    ((hi.dn - lo.dn) / (2. * dn), f.dnn),
                ] {
                    close(a, b, 1e-8, 3e-6);
                }
                close(v * f.dv + n * f.dn, f.value, 1e-12, 2e-12);
                close(v * f.dtv + n * f.dtn, f.dt, 1e-12, 2e-12);
                close(v * f.dvv, -n * f.dvn, 1e-7, 2e-12);
                close(v * f.dvn, -n * f.dnn, 1e-10, 2e-12);
            }
        }
    }
}

#[test]
fn caloric_and_mechanical_properties_obey_thermodynamic_identities() {
    let model = Pets::new().unwrap();
    for t in [80., 150., 300.] {
        for rho in [100., 10000., 25000.] {
            let n = 3.5;
            let v = n / rho;
            let p = model.caloric_tv(t, v, n).unwrap();
            let dt = t * 1e-5;
            let dv = v * 1e-5;
            let dn = n * 1e-5;
            let lo = model.caloric_tv(t - dt, v, n).unwrap();
            let hi = model.caloric_tv(t + dt, v, n).unwrap();
            close(
                (hi.internal_energy - lo.internal_energy) / (2. * dt),
                p.heat_capacity_v,
                1e-7,
                2e-6,
            );
            close(
                (hi.helmholtz_energy - lo.helmholtz_energy) / (2. * dt),
                -p.entropy,
                1e-7,
                2e-6,
            );
            close(
                (hi.pressure.value - lo.pressure.value) / (2. * dt),
                p.pressure.dt,
                1e-5,
                2e-6,
            );
            let lo = model.caloric_tv(t, v - dv, n).unwrap();
            let hi = model.caloric_tv(t, v + dv, n).unwrap();
            close(
                (hi.helmholtz_energy - lo.helmholtz_energy) / (2. * dv),
                -p.pressure.value,
                1e-4,
                2e-6,
            );
            close(
                (hi.entropy - lo.entropy) / (2. * dv),
                p.pressure.dt,
                1e-5,
                2e-6,
            );
            close(
                (hi.internal_energy - lo.internal_energy) / (2. * dv),
                t * p.pressure.dt - p.pressure.value,
                1e-4,
                2e-6,
            );
            let lo = model.caloric_tv(t, v, n - dn).unwrap();
            let hi = model.caloric_tv(t, v, n + dn).unwrap();
            close(
                (hi.helmholtz_energy - lo.helmholtz_energy) / (2. * dn),
                p.chemical_potential,
                1e-6,
                2e-6,
            );
            close(
                p.enthalpy - t * p.entropy,
                n * p.chemical_potential,
                1e-8,
                2e-12,
            );
            close(v * p.pressure.dv, -n * p.pressure.dn, 1e-5, 2e-12);
            let scaled = model.caloric_tv(t, v * 7.3, n * 7.3).unwrap();
            close(scaled.internal_energy, p.internal_energy * 7.3, 1e-8, 2e-12);
            close(scaled.entropy, p.entropy * 7.3, 1e-8, 2e-12);
            close(scaled.chemical_potential, p.chemical_potential, 1e-8, 2e-12);
        }
    }
}

#[test]
fn parameters_potential_cutoff_and_domains_are_explicit() {
    let mut model = Pets::new().unwrap();
    let other = model.clone();
    let before = model.parameters();
    assert!(model.set_parameters(0., 100.).is_err());
    assert_eq!(model.parameters(), before);
    model.set_parameters(4e-10, 250.).unwrap();
    assert_eq!(other.parameters(), before);
    let (sigma, eps) = model.parameters();
    let r = [sigma, 3. * sigma, sigma];
    let mut values = [99.; 3];
    let mut forces = [99.; 3];
    model
        .potential_into(&r, &mut values, Some(&mut forces))
        .unwrap();
    assert_eq!(values[1..], [0., 0.]);
    assert_eq!(forces[1..], [0., 0.]);
    close(forces[0], 24. * eps / sigma, 1e-4, 1e-14);
    assert_eq!(model.potential(2.5 * sigma).unwrap(), (0., 0.));
    let r = 1.1 * sigma;
    let dr = sigma * 1e-5;
    let force = model.potential(r).unwrap().1;
    close(
        -(model.potential(r + dr).unwrap().0 - model.potential(r - dr).unwrap().0) / (2. * dr),
        force,
        1.,
        2e-7,
    );
    assert!(model.prepare(0.).is_err());
    assert!(model.prepare(f64::NAN).is_err());
    assert!(model.residual_density(300., -1.).is_err());
    assert!(model.residual_density(300., 1e10).is_err());
    assert!(model.residual_helmholtz(300., 0., 1.).is_err());
    assert!(model.residual_helmholtz(300., 1., 0.).is_err());
    let vacuum = model.residual_density(300., 0.).unwrap();
    assert_eq!(vacuum.value, 0.);
    assert_eq!(vacuum.dt, 0.);
    assert!(vacuum.drho.is_finite());
}

#[test]
fn sound_matches_independent_constant_entropy_perturbations() {
    let model = Pets::new().unwrap();
    let mass = thermopack_rs::component("AR").unwrap().molar_mass * 1e-3;
    for rho in [100., 10000., 25000.] {
        let (t, n) = (300., 3.5);
        let v = n / rho;
        let dv = v * 1e-5;
        let entropy = model.caloric_tv(t, v, n).unwrap().entropy;
        let pressure = |volume| {
            let mut temperature = t;
            for _ in 0..5 {
                let p = model.caloric_tv(temperature, volume, n).unwrap();
                temperature -= (p.entropy - entropy) * temperature / p.heat_capacity_v;
            }
            model.pressure_tv(temperature, volume, n).unwrap().value
        };
        let dpdv = (pressure(v + dv) - pressure(v - dv)) / (2. * dv);
        let fd = (-v * v / (n * mass) * dpdv).sqrt();
        let sos = model.speed_of_sound_tv(t, v, n).unwrap();
        close(sos, fd, 1e-5, 1e-6);
        close(
            sos,
            model.baseline_speed_of_sound_tv(t, v, n).unwrap(),
            1e-12,
            1e-12,
        );
        close(
            sos,
            model.speed_of_sound_tv(t, v * 7.3, n * 7.3).unwrap(),
            1e-12,
            1e-12,
        );
        let state = model.caloric_tv(t, v, n).unwrap();
        let enthalpy_at = |temperature| {
            let mut volume = v;
            for _ in 0..6 {
                let p = model.pressure_tv(temperature, volume, n).unwrap();
                volume -= (p.value - state.pressure.value) / p.dv;
            }
            model.caloric_tv(temperature, volume, n).unwrap().enthalpy
        };
        let dt = 0.003;
        close(
            (enthalpy_at(t + dt) - enthalpy_at(t - dt)) / (2. * dt),
            state.heat_capacity_p.unwrap(),
            1e-5,
            1e-6,
        );
    }
}

#[test]
fn analytic_virial_limits_and_owned_caloric_references() {
    let mut model = Pets::new().unwrap();
    for t in [2.72, 40.8, 136., 300., 1360.] {
        let (b, c) = model.virial_coefficients(t).unwrap();
        let near = model.residual_density(t, 0.01).unwrap();
        close(near.value / 0.01, b, 1e-10, 1e-5);
        close(near.drhorho, c, 1e-13, 1e-5);
        // Finite-volume C loses digits by cancellation; bound its rounding
        // amplification using the magnitude of the two terms being subtracted.
        let (source_b, source_c) = model.baseline_virial_coefficients(t).unwrap();
        close(b, source_b, 1e-11, 1e-9);
        close(c, source_c, 16. * f64::EPSILON * b.abs() * 1e7, 1e-9);
    }
    let before = model.caloric_tv(300., 0.003, 3.).unwrap();
    let reference = model.ideal_gas().formation_enthalpy();
    model
        .ideal_gas_mut()
        .set_formation_enthalpy(reference + 3000.)
        .unwrap();
    let after = model.caloric_tv(300., 0.003, 3.).unwrap();
    close(
        after.internal_energy - before.internal_energy,
        9000.,
        1e-8,
        1e-12,
    );
    close(
        after.chemical_potential - before.chemical_potential,
        3000.,
        1e-8,
        1e-12,
    );
    close(after.pressure.value, before.pressure.value, 1e-9, 1e-12);
}
