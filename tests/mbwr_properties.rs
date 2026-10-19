use thermopack_rs::{Mbwr, MbwrKind, MbwrResidualTp, R, mbwr_records};

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn point(model: &Mbwr, t: f64, p: f64, n: f64, seed: f64) -> MbwrResidualTp {
    let mut v = seed;
    for _ in 0..30 {
        let state = model.pressure_tv(t, v, n).unwrap();
        let error = state.value - p;
        if error.abs() <= 1e-13 * p {
            return model.residual_tp(t, v, n).unwrap();
        }
        v = (v - error / state.dv).clamp(v * 0.5, v * 2.);
    }
    panic!("local pressure root did not converge");
}

#[test]
fn physical_tp_derivatives_satisfy_finite_differences_and_thermodynamic_identities() {
    for record in mbwr_records() {
        let model = Mbwr::from_record(record.source_name).unwrap();
        let t = 800.;
        let v = 1e-3;
        let n = 1.;
        let state = model.residual_tp(t, v, n).unwrap();
        let p = state.pressure;
        let dt = 0.001;
        let dp = p * 1e-5;
        let dn = 1e-5;
        close(
            state.enthalpy.value - t * state.entropy.value,
            state.gibbs_energy,
            1e-8,
            1e-10,
        );
        close(
            state.gibbs_energy,
            n * R * t * state.ln_fugacity_coefficient.value,
            1e-8,
            1e-12,
        );
        close(
            state.entropy.derivatives.unwrap().dt,
            state.enthalpy.derivatives.unwrap().dt / t,
            1e-14,
            1e-12,
        );
        let lo = point(&model, t - dt, p, n, v);
        let hi = point(&model, t + dt, p, n, v);
        close(
            (hi.gibbs_energy - lo.gibbs_energy) / (2. * dt),
            -state.entropy.value,
            1e-6,
            1e-5,
        );
        for (a, l, h) in [
            (state.enthalpy, lo.enthalpy, hi.enthalpy),
            (state.entropy, lo.entropy, hi.entropy),
            (
                state.compressibility,
                lo.compressibility,
                hi.compressibility,
            ),
            (
                state.ln_fugacity_coefficient,
                lo.ln_fugacity_coefficient,
                hi.ln_fugacity_coefficient,
            ),
        ] {
            close(
                a.derivatives.unwrap().dt,
                (h.value - l.value) / (2. * dt),
                1e-6,
                1e-5,
            );
        }
        let lo = point(&model, t, p - dp, n, v);
        let hi = point(&model, t, p + dp, n, v);
        close(
            (hi.gibbs_energy - lo.gibbs_energy) / (2. * dp),
            v - n * R * t / p,
            1e-9,
            1e-5,
        );
        for (a, l, h) in [
            (state.enthalpy, lo.enthalpy, hi.enthalpy),
            (state.entropy, lo.entropy, hi.entropy),
            (
                state.compressibility,
                lo.compressibility,
                hi.compressibility,
            ),
            (
                state.ln_fugacity_coefficient,
                lo.ln_fugacity_coefficient,
                hi.ln_fugacity_coefficient,
            ),
        ] {
            close(
                a.derivatives.unwrap().dp,
                (h.value - l.value) / (2. * dp),
                1e-9,
                1e-5,
            );
        }
        let lo = point(&model, t, p, n - dn, v * (n - dn) / n);
        let hi = point(&model, t, p, n + dn, v * (n + dn) / n);
        for (a, l, h) in [
            (state.enthalpy, lo.enthalpy, hi.enthalpy),
            (state.entropy, lo.entropy, hi.entropy),
            (
                state.compressibility,
                lo.compressibility,
                hi.compressibility,
            ),
            (
                state.ln_fugacity_coefficient,
                lo.ln_fugacity_coefficient,
                hi.ln_fugacity_coefficient,
            ),
        ] {
            close(
                a.derivatives.unwrap().dn,
                (h.value - l.value) / (2. * dn),
                1e-6,
                1e-5,
            );
        }
    }
}

#[test]
fn residual_reference_scaling_and_source_entropy_signs_are_explicit() {
    let ideal = Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 1.).unwrap();
    let t = 300.;
    let v = 1e-3;
    let n = 2.;
    let physical = ideal.residual_tp(t, v, n).unwrap();
    for q in [
        physical.enthalpy,
        physical.entropy,
        physical.ln_fugacity_coefficient,
    ] {
        close(q.value, 0., 1e-12, 0.);
        let d = q.derivatives.unwrap();
        close(d.dt, 0., 1e-12, 0.);
        close(d.dp, 0., 1e-12, 0.);
        close(d.dn, 0., 1e-12, 0.);
    }
    let baseline = ideal
        .baseline_residual_tp(t, physical.pressure, v, n)
        .unwrap();
    let sd = baseline.entropy.derivatives.unwrap();
    close(sd.dt, -2. * n * R / t, 1e-12, 1e-12);
    close(sd.dp, 2. * v / t, 1e-12, 1e-12);
    close(baseline.enthalpy.derivatives.unwrap().dt, 0., 1e-12, 0.);
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let a = model.residual_tp(t, v, n).unwrap();
    let b = model.residual_tp(t, v * 3.5, n * 3.5).unwrap();
    close(a.pressure, b.pressure, 1e-7, 1e-12);
    close(b.enthalpy.value, 3.5 * a.enthalpy.value, 1e-8, 1e-12);
    close(b.entropy.value, 3.5 * a.entropy.value, 1e-8, 1e-12);
    close(b.gibbs_energy, 3.5 * a.gibbs_energy, 1e-8, 1e-12);
    close(
        b.ln_fugacity_coefficient.value,
        a.ln_fugacity_coefficient.value,
        1e-12,
        1e-12,
    );
    assert_eq!(a.compressibility.derivatives.unwrap().dn, 0.);
    assert_eq!(a.ln_fugacity_coefficient.derivatives.unwrap().dn, 0.);
}

#[test]
fn singular_responses_preserve_values_and_invalid_inputs_are_rejected() {
    let mut coefficients = [0.; 19];
    coefficients[0] = -R / 2.;
    let model = Mbwr::with_parameters(MbwrKind::Bender19, &coefficients, 1.).unwrap();
    let q = model.residual_tp(300., 1e-3, 1.).unwrap();
    assert!(
        q.gibbs_energy.is_finite() && q.enthalpy.value.is_finite() && q.entropy.value.is_finite()
    );
    assert!(
        q.enthalpy.derivatives.is_none()
            && q.entropy.derivatives.is_none()
            && q.compressibility.derivatives.is_none()
    );
    assert!(q.ln_fugacity_coefficient.derivatives.is_some());
    assert!(model.baseline_residual_tp(300., 0., 1e-3, 1.).is_err());
    assert!(model.residual_tp(300., 0., 1.).is_err());
    assert!(model.residual_tp(0., 1., 1.).is_err());
    assert!(model.baseline_residual_tp(300., f64::NAN, 1., 1.).is_err());
    assert!(model.residual_tp(300., 1., 0.).is_err());
    assert!(model.residual_tp(300., 1e-4, 1.).is_err());
}
