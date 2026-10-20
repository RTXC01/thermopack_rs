use thermopack_rs::{LennardJonesKind, Multiparameter, Phase};

fn model(config: &str) -> Multiparameter {
    let c: Vec<_> = config.split(':').collect();
    match c[0] {
        "LJ_MEOS" => Multiparameter::lennard_jones(c[1], LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(c[1], LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => Multiparameter::new(c[1]).unwrap(),
    }
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn models() -> Vec<Multiparameter> {
    include_str!("fixtures/meos_metadata.csv")
        .lines()
        .skip(1)
        .map(|r| model(r.split('|').nth(1).unwrap()))
        .collect()
}

#[test]
fn all_constructor_metadata_matches_source() {
    let mut count = 0;
    for row in include_str!("fixtures/meos_metadata.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let config = c[0].split('|').nth(1).unwrap();
        let eos = model(config);
        let r = eos.record();
        let (mt, mp) = eos.lennard_jones_limits().unwrap_or((0., 0.));
        let (sig, eps) = eos
            .lennard_jones_parameters()
            .map_or((0., 0.), |p| (p.sigma, p.epsilon_over_k));
        let sizes = if config.starts_with("NIST") {
            [0.; 3]
        } else {
            [
                r.pol_terms.len() as f64,
                (r.pol_terms.len() + r.exp_terms.len()) as f64,
                (r.pol_terms.len() + r.exp_terms.len() + r.expexp_terms.len()) as f64,
            ]
        };
        let actual = [
            r.tc,
            r.pc,
            r.rc,
            r.acf,
            eos.baseline_molar_mass(),
            r.t_triple,
            r.p_triple,
            r.rholiq_triple,
            r.rhovap_triple,
            r.rgas_fit,
            eos.gas_constant(),
            mt,
            mp,
            sig,
            eps,
            sizes[0],
            sizes[1],
            sizes[2],
        ];
        for (i, a) in actual.into_iter().enumerate() {
            let atol = if i == 13 { 1e-22 } else { 1e-10 };
            close(a, c[i + 1].parse().unwrap(), atol, 2e-10);
        }
        count += 1;
    }
    assert_eq!(count, 14);
}

#[test]
fn tv_derivatives_caloric_identities_and_scalar_apis_close_independently() {
    for eos in models() {
        let t = 1.4 * eos.tc();
        let prep = eos.prepare(t).unwrap();
        for delta in [0.15, 1.8] {
            let v = 1. / (eos.rc() * delta);
            let s = prep.state(v).unwrap();
            let dt = t * 1e-5;
            let dv = v * 1e-5;
            let hi = eos.prepare(t + dt).unwrap().state(v).unwrap();
            let lo = eos.prepare(t - dt).unwrap().state(v).unwrap();
            close(
                s.entropy,
                -(hi.helmholtz_energy - lo.helmholtz_energy) / (2. * dt),
                1e-7,
                2e-8,
            );
            close(
                s.heat_capacity_v,
                (hi.internal_energy - lo.internal_energy) / (2. * dt),
                1e-7,
                2e-8,
            );
            close(
                s.pressure.dt,
                (hi.pressure.value - lo.pressure.value) / (2. * dt),
                1e-6,
                2e-8,
            );
            let hi = prep.state(v + dv).unwrap();
            let lo = prep.state(v - dv).unwrap();
            close(
                s.pressure.value,
                -(hi.helmholtz_energy - lo.helmholtz_energy) / (2. * dv),
                1e-5,
                2e-8,
            );
            close(
                s.pressure.dv,
                (hi.pressure.value - lo.pressure.value) / (2. * dv),
                1e-3,
                2e-8,
            );
            close(
                s.enthalpy,
                s.internal_energy + s.pressure.value * v,
                1e-10,
                1e-12,
            );
            close(s.gibbs_energy, s.enthalpy - t * s.entropy, 1e-10, 1e-12);
            close(
                s.heat_capacity_p,
                s.heat_capacity_v - t * s.pressure.dt.powi(2) / s.pressure.dv,
                1e-10,
                1e-12,
            );
            let w = s.speed_of_sound.unwrap();
            close(
                w * w,
                -v * v * s.pressure.dv * s.heat_capacity_p / (eos.molar_mass() * s.heat_capacity_v),
                1e-9,
                1e-12,
            );
            close(prep.pressure(v).unwrap().value, s.pressure.value, 0., 1e-14);
            close(prep.pressure(v).unwrap().dt, s.pressure.dt, 0., 1e-14);
            close(prep.pressure(v).unwrap().dv, s.pressure.dv, 0., 1e-14);
            close(prep.enthalpy(v).unwrap(), s.enthalpy, 1e-10, 1e-13);
            close(prep.entropy(v).unwrap(), s.entropy, 1e-10, 1e-13);
            close(
                prep.internal_energy(v).unwrap(),
                s.internal_energy,
                1e-10,
                1e-13,
            );
            close(
                prep.heat_capacity_v(v).unwrap(),
                s.heat_capacity_v,
                1e-10,
                1e-13,
            );
            close(
                prep.heat_capacity_p(v).unwrap(),
                s.heat_capacity_p,
                1e-10,
                1e-13,
            );
            close(prep.speed_of_sound(v).unwrap(), w, 1e-10, 1e-13);
            close(
                prep.baseline_speed_of_sound(v).unwrap(),
                s.baseline_speed_of_sound.unwrap(),
                1e-12,
                1e-13,
            );
            close(
                prep.ln_fugacity(v).unwrap(),
                s.ln_fugacity.unwrap(),
                1e-12,
                1e-13,
            );
        }
    }
}

#[test]
fn extensive_t_v_n_hessians_close_independent_differences_and_scaling() {
    for eos in models() {
        let t = 1.4 * eos.tc();
        let n = 2.3;
        let v = n / (0.3 * eos.rc());
        let prep = eos.prepare(t).unwrap();
        for ideal in [false, true] {
            let calc = |t, v, n| {
                let p = eos.prepare(t).unwrap();
                if ideal {
                    p.ideal_helmholtz(v, n)
                } else {
                    p.residual_helmholtz(v, n)
                }
                .unwrap()
            };
            let q = calc(t, v, n);
            for (axis, dx) in [t * 1e-5, v * 1e-5, n * 1e-5].into_iter().enumerate() {
                let mut a = [t, v, n];
                let mut b = a;
                a[axis] += dx;
                b[axis] -= dx;
                let a = calc(a[0], a[1], a[2]);
                let b = calc(b[0], b[1], b[2]);
                close(
                    [q.dt, q.dv, q.dn][axis],
                    (a.value - b.value) / (2. * dx),
                    1e-7,
                    1e-7,
                );
                let column = match axis {
                    0 => [q.dtt, q.dtv, q.dtn],
                    1 => [q.dtv, q.dvv, q.dvn],
                    _ => [q.dtn, q.dvn, q.dnn],
                };
                for (expected, (a, b)) in column
                    .into_iter()
                    .zip([a.dt, a.dv, a.dn].into_iter().zip([b.dt, b.dv, b.dn]))
                {
                    close(expected, (a - b) / (2. * dx), 1e-6, 1e-7);
                }
            }
            let scaled = if ideal {
                prep.ideal_helmholtz(3. * v, 3. * n)
            } else {
                prep.residual_helmholtz(3. * v, 3. * n)
            }
            .unwrap();
            close(scaled.value, 3. * q.value, 1e-10, 1e-13);
            close(scaled.dn, q.dn, 1e-10, 1e-13);
        }
    }
}

#[test]
fn lj_reduced_reference_units_mass_convention_and_invalid_inputs_are_explicit() {
    let raw = Multiparameter::reduced_lennard_jones_truncated_shifted();
    let v = 0.8 / 0.001;
    let s = raw.prepare(0.8).unwrap().ideal_jet(v).unwrap();
    close(s.d_tau, -1., 1e-14, 1e-14);
    close(s.d_tau - s.alpha, 0., 1e-14, 1e-14);
    assert_eq!(raw.gas_constant(), 1.);
    assert_eq!(raw.molar_mass(), 1.);
    for eos in models() {
        if eos.lennard_jones_parameters().is_some() {
            let w = eos
                .prepare(eos.tc() * 2.)
                .unwrap()
                .state(1. / eos.rc())
                .unwrap();
            close(
                w.speed_of_sound.unwrap() / w.baseline_speed_of_sound.unwrap(),
                1000_f64.sqrt(),
                1e-12,
                1e-13,
            );
            let m = eos.record().molar_mass;
            close(m * 1000., eos.baseline_molar_mass(), 1e-12, 1e-13);
            for phase in [Phase::Liquid, Phase::Vapor] {
                assert_eq!(
                    eos.prepare(eos.tc() * 1.1)
                        .unwrap()
                        .sat_delta_estimate(phase)
                        .unwrap(),
                    1.
                );
            }
        }
        for bad in [0., -1., f64::NAN, f64::INFINITY] {
            assert!(eos.prepare(bad).is_err());
        }
        let p = eos.prepare(eos.tc()).unwrap();
        for bad in [0., -1., f64::NAN, f64::INFINITY] {
            assert!(p.state(bad).is_err());
            assert!(p.ideal_jet(bad).is_err());
            assert!(p.residual_jet(bad).is_err());
            assert!(p.residual_helmholtz(1., bad).is_err());
            assert!(p.ideal_helmholtz(bad, 1.).is_err());
        }
    }
    assert!(Multiparameter::lennard_jones("C1", LennardJonesKind::Full).is_err());
}
